# Equivariant G-complex simplification

Load HAP in a fresh GAP session and use `ReduceGComplex(L)` to return a
reduced complex. The input uses the HAP G-complex interface: `dimension`,
`boundary`, `stabilizer`, `action`, `elts`, and `group`.

The algorithm repeatedly replaces balanced stars that pass the implemented
disk-collapse criterion. `HomeomorphicToDisk` implements that criterion;
it is not a general algorithm for recognizing arbitrary topological disks.
The replacement does not construct a contracting homotopy (`homotopy = fail`).

Progress messages are off by default. Enable them with:

```gap
SetInfoLevel(InfoWunden, 1);
```

Declarations are in `decsWunden.gd`; `initWunden.gi` loads the seven modules
in dependency order. 

Run the small regression suites in a fresh GAP session:

```gap
LoadPackage("hap");
ReadPackage("hap", "tst/Wunden/core.g");
ReadPackage("hap", "tst/Wunden/integration.g");
```

Optional, slower homology regressions:

```gap
ReadPackage("hap", "tst/Wunden/sl3zs.g");
ReadPackage("hap", "tst/Wunden/bianchi.g");
```

These test scripts assign global example variables; use a disposable session.
The examples are regression checks, not a proof for all input complexes.

## Bruhat--Tits buildings and restriction

```gap
LoadPackage("hap");
B := BruhatTitsBuildingSL(3, 2);;
# K must be a finite-index subgroup of B!.group.
S := RestrictNonFreeResolution(B, K);;
```

`BruhatTitsBuildingSL(n,p)` requires an integer `n >= 2` and prime `p`.
It supplies the chamber cell orbits, signed boundaries, and infinite cell
stabilizers for SL(n,Z[1/p]). Stabilizer membership is tested by preservation
of the stored lattice bases (`BTBuildingLatticeBases`). No contracting
homotopy is supplied. Routines requiring finite stabilizers do not apply.

`RestrictNonFreeResolution(R,K)` requires finite index, terminating subgroup
membership tests, and the same group-element representation for K and R.
The input must have a finite nonnegative recorded length. Infinite-index
inputs may cause enumeration not to terminate; subgroup inclusion is only
checked automatically when the ambient group is known finite.

The output includes `toOriginal`, `fromOriginal`, `restrictionCells`, and
`cosetRepresentatives`. It transports an available contracting homotopy and
preserves `fail` when none is available. Restriction may extend R!.elts and
populate caches. `RestrictionCellStabilizerData` supports exact membership
in the restricted cell stabilizers.

Run the regression suites in a disposable session:

```gap
ReadPackage("hap", "tst/Wunden/test_building_general.g");
ReadPackage("hap", "tst/Wunden/test_restrict_nonfree.g");
```

The building suite reads `building_sl3_reference.g`, a test-only copy of the
original rank-three constructor, for independent comparison. That constructor
is not loaded as part of the package implementation. The tests assign global
names and are currently run manually rather than through HAP's test runner.
Do not also read the original standalone implementation scripts in the same
session: their global assignments conflict with the package declarations.
