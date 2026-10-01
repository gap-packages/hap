This file describes changes in the HAP package.

## Unreleased

- Add a first implementation of congruence subgroups of SL(4,Z), and faster
  Γ₀(N) coset methods for SL(3,Z) (#161)
- Update `AmbientPosition` and `AmbientTree` for SL(3,Z)
- Add `FiniteProjectiveSpace3` and an alternative construction
  `FiniteProjectivePlane_alt` of finite projective planes
- Update the method for congruence subgroups of SL(2,Z) of composite level
  (#160)
- Require the congruence package

## 1.79 (2026-09-01)

- Rewrite `IsAspherical` to use the polymaking interface; it now catches and
  reports errors (#158)
- Add `DualComplex`, `TensorWithIntegersMod2Torsion`, `HomToRationals`,
  `ResolutionFiniteSubgroup_NonFree` and an `EulerCharacteristic` method for
  regular CW complexes
- Add `VoronoiGenerators`, `VoronoiWord` and related functions
- Fix `GroupHomology(f,n,p)` when `f` does not map a Sylow p-subgroup of its
  source into the chosen Sylow p-subgroup of its target
- Rename `CuspidalCohomologyHomomorphism` to `InteriorCohomologyHomomorphism`

## 1.78 (2026-07-12)

- Janitorial changes

## 1.77 (2026-07-12)

- Fix the version requirement for the congruence package (#155)
- Add `GComplexToFiniteRegularCWRegion` and `RegularCWCellClosure`

## 1.76 (2026-07-08)

- Add `ContractibleGcomplex` data for SL(4,Z)
- List the needed system packages (graphviz, ImageMagick, Singular, polymake) in
  `PackageInfo.g` (#154)

## 1.75 (2026-03-24)

- Add a new implementation of congruence subgroups of SL(2,Z) and SL(3,Z),
  including Γ₀ subgroups of SL(3,Z), with `CongruenceSubgroupGamma0`,
  `CongruenceSubgroupGamma1`, `PrincipalCongruenceSubgroup`, `AmbientTree`,
  `AmbientTransversal`, `IndexInAmbientGroup` and `LevelOfCongruenceSubgroup`
  (#146, #149, #150, #151, #152)
- Add `FiniteProjectiveLine`, `FiniteProjectivePlane` and
  `GComplexToRegularCWComplex`
- Require the polymaking package, which does not need polymake installed, and
  fix warnings when it is not loaded (#153)
- Fix brackets in some precomputed barcode files (#145)

## 1.74 (2026-01-10)

- Adapt to the changed return type of Singular's `hilb` in newer Singular
  versions (#142, #143)
- Remove the `POLYMAKE_PATH` variable (#140)
- Require the SmallGrp package
- Improve transversals of congruence subgroups
- Add `SparseIdentityMat`, `SparseMatAddToEntry` and `SparseMatConcatenation`

## 1.73 (2025-12-22)

- Fix `CuspidalCohomologyHomomorphism`

## 1.72 (2025-12-21)

- Janitorial changes

## 1.71 (2025-12-21)

- Require GAP >= 4.12 (#137)
- Remove the legacy polymake interface; use the polymaking package and its
  `POLYMAKE_COMMAND` to locate polymake (#136, #138)
- Fix `BarComplexOfMonoid` and the visualisation in `OrbitPolytope`
- Add `BianchiGcomplex`, `WeakCommutativityCommutatorGroup` and
  `SymmetricCommutativityCommutatorGroup`

## 1.70 (2025-07-19)

- Extend the Bianchi groups chapter of the tutorial

## 1.69 (2025-07-17)

- Improve functions for Bianchi groups

## 1.68 (2025-07-09)

- Janitorial changes

## 1.67 (2025-07-09)

- Speed up `GroupHomology`
- Extend `SimplicialMap` from inclusions to arbitrary simplicial maps
- Improve Swan's algorithm for Bianchi groups and the Hecke operator functions;
  add `ResolutionAbelianBianchiSubgroup`
- Fix arithmetic with quadratic numbers
- Add `ModPCohomologyPresentationBounds`, `ResolutionFiniteCcGroup` and
  `ResolutionInfiniteCcGroup`

## 1.66 (2024-10-24)

- Do not read files that no longer exist (#124)

## 1.65 (2024-07-29)

- Update precomputed data for Bianchi groups

## 1.64 (2024-07-28)

- Add functions for Bianchi groups, unimodular pairs and quadratic numbers, such
  as `BianchiPolyhedron`, `CoverOfUnimodularPairs`, `DisplayUnimodularPairs` and
  `QuadraticNumber`
- Add `ParallelPersistentBettiNumbers`

## 1.63 (2024-03-20)

- Extend persistent Betti number functions
- Add `ClassifyingSpaceFiniteGroup`, `RegularCWComplexReordered`,
  `CupProductOfRegularCWComplexModP` and suspensions
- Speed up the Chevalley-Eilenberg complex

## 1.62 (2024-02-01)

- Add `ChevalleyEilenbergComplexOfModule`, `FiltrationTerms` and
  `ExpandedComplex`
- Define NC versions of the `PreImages...` functions if missing (#121)

## 1.61 (2024-01-02)

- Remove `InitialObject` and `TerminalObject` (#118)
- Do not call `ObjectifyWithAttributes` on a string (#117)

## 1.60 (2023-10-15)

- Janitorial changes

## 1.59 (2023-10-15)

- Add functions for homomorphisms on cohomology, such as
  `PrimePartDerivedFunctorHomomorphism` and `DirectProductOfGroupHomomorphisms`
- Add a new Bockstein implementation
- Add `PersistentBettiNumbersViaContractions`

## 1.58 (2023-08-06)

- Add `DirectProductOfSimplicialComplexes`
- Improve Bockstein computations for spaces

## 1.57 (2023-07-25)

- Add `Suspension`
- Add Bockstein homomorphisms for spaces

## 1.56 (2023-05-24)

- Add `RegularCWAssociahedron` and functions for planar binary trees
- Add `LowDimensionalCupProduct`, `CupProductMatrix`,
  `SignatureOfSymmetricMatrix` and `DiagonalChainMap`
- Add `PoincareBipyramidCWComplex` and `DisplayVectorField`

## 1.55 (2023-04-17)

- Add `ManifoldType`, `PoincareDodecahedronCWComplex`,
  `PoincareOctahedronCWComplex` and `PoincarePrismCWComplex`

## 1.54 (2023-03-19)

- Add `BarycentricallySimplifiedComplex`, `NonManifoldVertices`, `RemoveStar`,
  `ThreeManifoldWithBoundary` and `PoincareCubeCWComplexNS`
- Make `compile.sh` easier to use (#109)

## 1.53 (2023-02-27)

- Add `IsClosedManifold` and `PoincareCubeCWComplex`

## 1.52 (2023-02-11)

- Janitorial changes

## 1.51 (2023-02-11)

- Do not call `MutableCopyMat` on vectors (#107)

## 1.50 (2023-02-02)

- Add `PSubgroupSimplicialComplex`, `PSubgroupGChainComplex` and
  `HomologicalGroupDecomposition` for equivariant Quillen complexes
- Add `TensorWithModPModule`
- Fix a bad entity in the tutorial (#102)

## 1.49 (2023-01-07)

- Add `PrimePartDerivedTwistedFunctor`

## 1.48 (2023-01-02)

- Add a diagonal for regular CW complexes and improve cup products
- Add `RegularCWCube`, `RegularCWSimplex`, `RegularCWPolygon`,
  `RegularCWPermutahedron`, `DirectProductOfRegularCWComplexesLazy`,
  `ComposeCWMaps`, `ChainComplexWithChainHomotopy`, `QuotientChainMap` and
  `HomToModPModule`
- Change `CohomologicalData`
- Avoid hard coded assumptions about paths (#101)
- Fix a bad entity in the tutorial (#94)

## 1.47 (2022-08-14)

- Add `CohomologicalData`, `HeckeOperator`, `ResolutionSL2ZConjugated`,
  `MinimizeRingRelations` and `TransferCochainMap`

## 1.46 (2022-07-25)

- Janitorial changes

## 1.45 (2022-07-24)

- Change `ResolutionAbelianGroup` and `ResolutionSpaceGroup`; add
  `IsPeriodicSpaceGroup`
- Fix ImageMagick error "pixels are not authentic"
- Fix `ViewPureCubicalKnot`

## 1.44 (2022-07-06)

- Add `ResolutionSpaceGroup` and `CrystallographicComplex`

## 1.43 (2022-06-29)

- Janitorial changes

## 1.42 (2022-06-28)

- Make `GroupHomology` and `GroupCohomology` handle almost crystallographic pcp
  groups

## 1.41 (2022-06-02)

- Janitorial changes

## 1.40 (2022-06-01)

- Add `LinkingFormHomotopyInvariant` and `LinkingFormHomeomorphismInvariant`,
  replacing `LinkingFormInvariant`
- Add `RandomArc2Presentation`, `ViewArc2Presentation`, `KinkArc2Presentation`
  and `NumberOfCrossingsInArc2Presentation`
- Add `VertexLink`, `VertexStar`, `Tube`, `SequentialRegularCWComplexComplement`
  and `ChainComplexHomeomorphismEquivalenceOfRegularCWComplex`

## 1.39 (2022-04-20)

- Add `DijkgraafWittenInvariant`, `LinkingForm`, `ThreeManifoldViaDehnSurgery`,
  `FundamentalGroupWithPathReps` and `RegularCWComplexWithRemovedCell`

## 1.38 (2022-03-09)

- Add `IdentifyKnot` and `ReadLinkImageAsGaussCode`

## 1.37 (2022-02-11)

- Fix broken links in the documentation (#58)

## 1.36 (2022-02-10)

- Janitorial changes

## 1.35 (2022-02-09)

- Add functions for non-free resolutions
- Make `StarGraph` an operation (#59)

## 1.34 (2021-07-20)

- Add cohomology rings of simplicial complexes (`CohomologyRing`); improve cup
  products and their documentation
- Rename `BaryCentricSubdivision` to `BarycentricSubdivision`
- Fix compatibility with polymake 4

## 1.33 (2021-06-30)

- Add functions for simplicial manifolds and intersection forms, such as
  `ClosedSurface`, `ConnectedSum`, `WedgeSum`, `Sphere`,
  `ComplexProjectiveSpace` and `SimplicialK3Surface`

## 1.32 (2021-06-16)

- Fix the tutorial address

## 1.31 (2021-06-16)

- Add `WeakCommutativityGroup`, `SymmetricCommutativityGroup` and
  `Nil3TensorSquare`
- Add `ArcDiagramToTubularSurface`, `LiftColouredSurface`,
  `ViewColouredArcDiagram` and `GModuleAsGOuterGroup`
- Fix `GroupHomology` and extend the short manual
- Fix the declaration of `IsHapQuotientElementRep` (#53)

## 1.30 (2021-04-07)

- Add functions for Eilenberg-MacLane spaces:
  `EilenbergMacLaneSimplicialFreeAbelianGroup`,
  `HomologySimplicialFreeAbelianGroup` and
  `CohomologySimplicialFreeAbelianGroup`
- Add `BarComplexOfMonoid`, `E1HomologyPage`, `E1CohomologyPage` and
  `PathObjectForChainComplex`
- Upgrade the polymake functions; suggest the Polymaking package
- Update `ResolutionAbelianGroup`

## 1.29 (2021-01-07)

- Add `SpunLinkComplement`; rename `SpunAboutInitialHyperplane` to
  `SpunAboutHyperplane`

## 1.28 (2021-01-05)

- Add functionality for Hecke operators (`HeckeOperatorWeight2`)
- Add functions for regular CW complexes, such as `ClosureCWCell`,
  `IntersectionCWSubcomplex`, `PathComponentsCWSubcomplex` and
  `RegularCWComplexComplement`
- Add `ChainComplexToSparseChainComplex` and `SparseChainComplexToChainComplex`

## 1.27 (2020-05-20)

- Fix a bug with congruence subgroups and update the manual
- Fix the package URL (#34)

## 1.26 (2020-05-04)

- Add functions for congruence subgroups and quadratic number fields, such as
  `QuadraticNumberField`, `RingOfQuadraticIntegers`,
  `ResolutionSL2QuadraticIntegers`, `IndexInSL2Z` and `HeckeOperator`
- Fix `HomogeneousPolynomials`

## 1.25 (2020-01-25)

- Do not call `ObjectifyWithAttributes` on a string (#27)

## 1.24 (2019-12-11)

## 1.23 (2019-11-14)

## 1.22 (2019-11-14)

## 1.21 (2019-11-14)

## 1.20 (2019-07-17)

## 1.18 (2018-11-25)

## 1.17 (2018-11-12)

## 1.16 (2018-11-05)

## 1.15 (2018-09-19)

## 1.13 (2018-09-19)

## 1.12.7 (2018-09-18)

## 1.12.6 (2018-04-16)

## 1.12.5 (2017-11-21)

## 1.12.4 (2017-11-21)

## 1.12.2 (2017-09-20)

## 1.12.1 (2017-09-07)

## 1.12.0 (2017-08-27)

## 1.11.15 (2017-02-21)

## 1.11.14 (2017-02-04)

## 1.11.12 (2015-11-11)

## 1.11.11 (2015-11-11)

## 1.11.13 (2015-11-03)

## 1.11.7 (2015-10-28)

## 1.11.6 (2015-10-28)

## 1.11.5 (2015-10-28)

## 1.11.4 (2015-10-26)

## 1.11.3 (2015-10-23)

## 1.11.1 (2015-07-29)

## 1.11 (2015-05-18)

## 1.10.15 (2013-12-07)

## 1.10.14.3 (2013-11-21)

## 1.10.14.2 (2013-11-21)

## 1.10.14.1 (2013-11-04)

## 1.10.14 (2013-10-22)

## 1.10.13 (2013-08-12)

## 1.10.12 (2013-07-08)

## 1.10.11 (2013-07-01)

## 1.10.10.2 (2013-03-05)

## 1.10.10.1 (2013-03-05)

## 1.10.9.6 (2013-03-05)

## 1.10.9.5 (2013-01-24)

## 1.10.9.4 (2013-01-24)

## 1.10.9.3 (2013-01-24)

## 1.10.9.2 (2013-01-24)

## 1.10.9.1 (2013-01-02)

## 1.10.9 (2013-01-02)

## 1.10.8 (2012-06-19)

## 1.10.7 (2012-06-19)

## 1.10.6.1 (2012-06-12)

## 1.10.6 (2012-06-12)

## 1.10.5 (2012-06-12)

## 1.10.4 (2012-06-03)

## 1.10.3 (2012-06-03)

## 1.10.2 (2012-05-29)

## 1.10.1 (2012-05-15)

## 1.10.0 (2012-03-16)

## 1.9.5 (2012-03-16)

## 1.9.3 (2010-12-13)

## 1.9.2 (2010-07-21)

## 1.9.4 (2010-04-26)

## 1.9.1 (2009-12-16)

## 1.9 (2009-10-10)

## 1.8.9 (2008-12-15)

## 1.8.8 (2008-07-28)

## 1.8.7 (2008-07-13)

## 1.8.6 (2008-02-01)

## 1.8.5 (2008-01-20)

## 1.8.4 (2007-12-04)

## 1.8.3 (2007-10-07)

## 1.8.2 (2007-09-06)

## 1.8 (2007-08-16)

## 1.7.4 (2007-04-10)

## 1.7.3 (2007-03-09)

## 1.7 (2006-09-01)

## 1.6 (2006-08-15)

## 1.5 (2006-06-11)

## 1.4 (2006-05-03)

## 1.3 (2006-03-24)

## 1.2 (2006-03-07)

## 1.1 (2006-01-31)
