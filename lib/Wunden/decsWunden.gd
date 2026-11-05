# Declarations for equivariant G-complex simplification.
DeclareInfoClass("InfoWunden");
DeclareGlobalFunction("kStar");
DeclareGlobalFunction("ForgetG");
DeclareGlobalFunction("PreIsGBalanced");
DeclareGlobalFunction("IsGBalanced");
DeclareGlobalFunction("Stabstar");
DeclareGlobalFunction("EntireBoundary");
DeclareGlobalFunction("Cobound");
DeclareGlobalFunction("ShouldRemove");
DeclareGlobalFunction("JoinCells");
DeclareGlobalFunction("SimplifyLattice");
DeclareGlobalFunction("HomeomorphicToDisk");
DeclareGlobalFunction("Replace2");
DeclareGlobalFunction("ReduceGComplex");

# Bruhat--Tits buildings and finite-index restriction.
DeclareAttribute("BTBuildingLatticeBases", IsGroup);
DeclareAttribute("RestrictionCellStabilizerData", IsGroup);
DeclareGlobalFunction("BruhatTitsBuildingSL");
DeclareGlobalFunction("RestrictNonFreeResolution");
