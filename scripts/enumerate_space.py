"""Measure the funnel across the whole combination space.

The library decomposed out of PROTAC-8K spans ~1.65e8 formal combinations.
Enumerating that to find out how many survive filtering would take days and is
unnecessary, because the two stages have very different shapes:

* **The attachment rule is separable.** Whether a bond may be formed depends
  only on the elements of the two atoms being joined, so the number of
  survivors is a sum over element pairs and can be counted *exactly* in
  milliseconds, with no molecule built. The honest headline number is free.

* **Structure alerts and properties are not separable.** They depend on the
  assembled molecule, so these are estimated from a uniform random sample,
  reported with a Wilson interval. A sample of 20,000 pins the pass rate to
  well under a percentage point, which is finer than the question deserves.

Usage::

    python scripts/enumerate_space.py [--sample 20000] [--seed 0]
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from rdkit import RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")

from lasprotacs.assemble import assemble  # noqa: E402
from lasprotacs.rules import attachment_atoms, bond_ok, new_alerts  # noqa: E402
from lasprotacs.score import developability, properties  # noqa: E402

SYMBOL = {6: "C", 7: "N", 8: "O", 16: "S", 15: "P", 9: "F"}


def exact_attachment_count(warheads, linkers, e3s):
    """Exact number of combinations passing the attachment rule.

    ``attachment_ok`` factorises: the warhead side and the E3 side are
    independent, and each depends only on a pair of element numbers. So group
    the components by the element at their attachment point and count pairs.
    """
    wh = Counter()
    for c in warheads:
        a = attachment_atoms(c["smiles"]).get(1)
        if a:
            wh[a] += 1
    e3 = Counter()
    for c in e3s:
        a = attachment_atoms(c["smiles"]).get(2)
        if a:
            e3[a] += 1
    lk = Counter()
    for c in linkers:
        m = attachment_atoms(c["smiles"])
        if 1 in m and 2 in m:
            lk[(m[1], m[2])] += 1

    total = 0
    per_linker_end = Counter()
    for (l1, l2), nl in lk.items():
        left = sum(n for e, n in wh.items() if bond_ok(e, l1))
        right = sum(n for e, n in e3.items() if bond_ok(e, l2))
        total += nl * left * right
        per_linker_end[(SYMBOL.get(l1, l1), SYMBOL.get(l2, l2))] += nl
    return total, wh, lk, e3, per_linker_end


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    lib = json.loads((ROOT / "data" / "components.json").read_text())
    W, L, E = lib["warheads"], lib["linkers"], lib["e3s"]
    raw = len(W) * len(L) * len(E)
    print(f"library   : {len(W)} warheads x {len(L)} linkers x {len(E)} E3 ligands")
    print(f"raw space : {raw:,}")

    passed, wh, lk, e3, ends = exact_attachment_count(W, L, E)
    print(f"\nattachment rule (exact, no enumeration): {passed:,}  "
          f"({100 * passed / raw:.1f}%)")
    print("  warhead exit atoms :",
          {SYMBOL.get(k, k): v for k, v in wh.most_common()})
    print("  E3 exit atoms      :",
          {SYMBOL.get(k, k): v for k, v in e3.most_common()})
    print("  linker ends (warhead-side, E3-side):",
          {f"{a}/{b}": n for (a, b), n in ends.most_common(6)})

    rng = random.Random(args.seed)
    n = built = clean = 0
    devs = []
    print(f"\nsampling {args.sample:,} combinations for the post-assembly stages")
    while n < args.sample:
        w, l, e = rng.choice(W), rng.choice(L), rng.choice(E)
        wa, la, ea = (attachment_atoms(w["smiles"]), attachment_atoms(l["smiles"]),
                      attachment_atoms(e["smiles"]))
        if 1 not in wa or 1 not in la or 2 not in la or 2 not in ea:
            continue
        if not (bond_ok(wa[1], la[1]) and bond_ok(ea[2], la[2])):
            continue  # already counted exactly above
        n += 1
        smi = assemble(w["smiles"], l["smiles"], e["smiles"])
        if smi is None:
            continue
        built += 1
        props = properties(smi)
        if props is None:
            continue
        if not new_alerts(smi, [w["smiles"], l["smiles"], e["smiles"]]):
            clean += 1
            devs.append(developability(props))

    lo, hi = wilson(clean, n)
    print(f"  assembled cleanly : {built}/{n}")
    print(f"  no new alerts     : {clean}/{n} = {100 * clean / n:.1f}% "
          f"[95% CI {100 * lo:.1f}-{100 * hi:.1f}]")
    print(f"  => of the {passed:,} that pass attachment, "
          f"~{int(passed * clean / n):,} survive the alert filter")
    if devs:
        devs.sort()
        print(f"  developability    : median {devs[len(devs) // 2]:.2f}, "
              f"top decile >= {devs[int(0.9 * len(devs))]:.2f}")

    out = {
        "raw": raw,
        "attachment": passed,
        "attachment_exact": True,
        "sample": n,
        "sample_clean": clean,
        "clean_rate": clean / n,
        "clean_ci": [lo, hi],
        "clean_estimate": int(passed * clean / n),
    }
    (ROOT / "data" / "funnel.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {ROOT / 'data' / 'funnel.json'}")


if __name__ == "__main__":
    main()
