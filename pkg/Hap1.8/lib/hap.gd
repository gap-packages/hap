#(C) Graham Ellis, 2005-2006.

## FREE G MODULES ###################################################
DeclareGlobalFunction("Negate");
DeclareGlobalFunction("NegateWord");
DeclareGlobalFunction("AlgebraicReduction");
DeclareGlobalFunction("AddFreeWords");
DeclareGlobalFunction("AddFreeWordsModP");
DeclareGlobalFunction("PrintZGword");
DeclareGlobalFunction("TietzeReduction");
DeclareGlobalFunction("MultiplyWord");
DeclareGlobalFunction("WordModP");

## FpG MODULES#######################################################
DeclareGlobalFunction("DesuspensionFpGModule");
DeclareGlobalFunction("RadicalOfFpGModule");
DeclareGlobalFunction("GeneratorsOfFpGModule");
DeclareGlobalFunction("ComplementaryBasis");
DeclareGlobalFunction("FpGModule");
DeclareGlobalFunction("DirectSumOfFpGModules");
DeclareGlobalFunction("ResolutionFpGModule");
DeclareGlobalFunction("IsFpGModuleHomomorphismData");
DeclareGlobalFunction("FpGModuleDualBasis");
DeclareGlobalFunction("MultipleOfFpGModule");
DeclareGlobalFunction("IntersectionOfFpGModules");
DeclareGlobalFunction("SumOfFpGModules");
DeclareGlobalFunction("ProjectedFpGModule");
DeclareGlobalFunction("RandomHomomorphismOfFpGModules");
DeclareGlobalFunction("VectorsToFpGModuleWords");
DeclareGlobalFunction("FpGModuleHomomorphismNC");
DeclareGlobalFunction("FpGModuleHomomorphism");
DeclareGlobalFunction("ImageOfFpGModuleHomomorphism");
DeclareGlobalFunction("CompositionOfFpGModuleHomomorphisms");


## NONABELIAN TENSOR ################################################
DeclareGlobalFunction("NonabelianTensorSquare");
DeclareGlobalFunction("NonabelianTensorProduct");
DeclareGlobalFunction("ThirdHomotopyGroupOfSuspensionB");
DeclareGlobalFunction("NonabelianExteriorProduct");
DeclareGlobalFunction("RelativeSchurMultiplier");
DeclareGlobalFunction("EpiCentre");
DeclareGlobalFunction("UpperEpicentralSeries");
DeclareGlobalFunction("BaerInvariant");
DeclareGlobalFunction("TensorCentre");
DeclareGlobalFunction("ThirdHomotopyGroupOfSuspensionB_alt");
DeclareGlobalFunction("NonabelianTensorSquare_inf");


## RESOLUTIONS ######################################################
DeclareGlobalFunction("ResolutionFiniteGroup");
DeclareGlobalFunction("ResolutionSmallFpGroup");
DeclareGlobalFunction("PresentationOfResolution");
DeclareGlobalFunction("ResolutionFiniteSubgroup");
DeclareGlobalFunction("ResolutionSubgroup");
DeclareGlobalFunction("ResolutionAsphericalPresentation");
DeclareGlobalFunction("ResolutionAbelianGroup");
DeclareGlobalFunction("ResolutionAlmostCrystalGroup");
DeclareGlobalFunction("ResolutionAlmostCrystalQuotient");
DeclareGlobalFunction("CayleyGraphDisplay");

## RESOLUTIONS MOD P ################################################
DeclareGlobalFunction("ResolutionPrimePowerGroup");
DeclareGlobalFunction("RankPrimeHomology");
DeclareGlobalFunction("RankHomologyPGroup");
DeclareGlobalFunction("PoincareSeries");
DeclareGlobalFunction("PoincareSeries_alt");
DeclareGlobalFunction("PoincareSeriesApproximation");
DeclareGlobalFunction("PoincareSeriesPrimePart");
DeclareGlobalFunction("ExpansionOfRationalFunction");
DeclareGlobalFunction("EfficientNormalSubgroups");


## FUNCTORS #########################################################
DeclareGlobalFunction("EvaluateProperty");
DeclareGlobalFunction("EvenSubgroup");
DeclareGlobalFunction("EquivariantChainMap");
DeclareGlobalFunction("TensorWithIntegers");
DeclareGlobalFunction("TensorWithIntegersModP");
DeclareGlobalFunction("PrimePartDerivedFunctor");
DeclareGlobalFunction("ReduceGenerators");
DeclareGlobalFunction("HomToIntegers");
DeclareGlobalFunction("HomToIntegralModule");
DeclareGlobalFunction("PermToMatrixGroup");
DeclareGlobalFunction("AbelianInvariantsToTorsionCoefficients");
DeclareGlobalFunction("TorsionGeneratorsAbelianGroup");
DeclareGlobalFunction("TensorWithRationals");
DeclareGlobalFunction("BigStepLCS");
DeclareGlobalFunction("HomToIntegersModP");
DeclareGlobalFunction("Coclass");
DeclareGlobalFunction("BoundaryMatrix");
DeclareGlobalFunction("Prank");
DeclareGlobalFunction("PrankAlt");
DeclareOperation("Compose",[IsGroupHomomorphism,IsGroupHomomorphism]);


## PERTURBATIONS ####################################################
DeclareGlobalFunction("TwistedTensorProduct");
DeclareGlobalFunction("ResolutionFiniteExtension");
DeclareGlobalFunction("ResolutionNormalSeries");
DeclareGlobalFunction("ResolutionDirectProduct");
DeclareGlobalFunction("ResolutionExtension");
DeclareGlobalFunction("ResolutionFiniteDirectProduct");
DeclareGlobalFunction("ResolutionSubnormalSeries");



## ARTIN COXETER ####################################################
DeclareGlobalFunction("CoxeterDiagramMatrix");
DeclareGlobalFunction("CoxeterDiagramVertices");
DeclareGlobalFunction("CoxeterDiagramFpArtinGroup");
DeclareGlobalFunction("CoxeterDiagramFpCoxeterGroup");
DeclareGlobalFunction("CoxeterSubDiagram");
DeclareGlobalFunction("CoxeterDiagramComponents");
DeclareGlobalFunction("CoxeterDiagramDegree");
DeclareGlobalFunction("CoxeterDiagramIsSpherical");
DeclareGlobalFunction("ResolutionArtinGroup");
DeclareGlobalFunction("CoxeterDiagramDisplay");

## HOMOLOGY #########################################################
DeclareOperation("Homology",[IsObject,IsObject]);
#DeclareGlobalFunction("Homology");
DeclareGlobalFunction("IntegralHomology");
DeclareGlobalFunction("ModularHomology");
DeclareGlobalFunction("GroupHomology");
DeclareGlobalFunction("IntegralCohomology");
DeclareOperation("Cohomology",[IsObject,IsObject]);
#DeclareGlobalFunction("Cohomology");
DeclareGlobalFunction("Syzygy");
DeclareGlobalFunction("CR_IntegralCycleToClass");
DeclareGlobalFunction("CocycleCondition");
DeclareGlobalFunction("StandardCocycle");
DeclareGlobalFunction("IsSuperperfect");
DeclareGlobalFunction("ModularCohomology");
DeclareOperation("SolutionsMatDestructive",
			[IsOrdinaryMatrix,IsOrdinaryMatrix]);
DeclareGlobalFunction("HomologyPb");
DeclareGlobalFunction("HomologyPrimePart");
DeclareGlobalFunction("CohomologyPrimePart");
DeclareGlobalFunction("GroupCohomology");
DeclareGlobalFunction("IntegralHomologyOfChainComplex");
DeclareGlobalFunction("IntegralCohomologyOfCochainComplex");

## RINGS ############################################################
DeclareGlobalFunction("CR_IntegralCohomology");
DeclareGlobalFunction("CR_ChainMapFromCocycle");
DeclareGlobalFunction("CR_CocyclesAndCoboundaries");
DeclareGlobalFunction("CR_IntegralClassToCocycle");
DeclareGlobalFunction("CR_IntegralCocycleToClass");
DeclareGlobalFunction("IntegralCupProduct");
DeclareGlobalFunction("IntegralRingGenerators");

## ModPRings ########################################################
DeclareGlobalFunction("ModPCohomologyRing");
DeclareGlobalFunction("ModPCohomologyRing_part_1");
DeclareGlobalFunction("ModPCohomologyRing_part_2");

DeclareGlobalFunction("ModPRingGeneratorsAlt");
DeclareGlobalFunction("ModPRingGenerators");
DeclareGlobalFunction("ModPCohomologyGenerators");


## CURVATURE ########################################################
DeclareGlobalFunction("IsAspherical");
DeclareGlobalFunction("PolytopalGenerators");
DeclareGlobalFunction("VectorStabilizer");
DeclareGlobalFunction("PolytopalComplex");
DeclareGlobalFunction("OrbitPolytope");

## POLYCYLIC ########################################################
DeclareGlobalFunction("ResolutionAbelianPcpGroup");
DeclareGlobalFunction("ResolutionNilpotentGroup");

## OBJECTIFICATION ##################################################
DeclareOperation("Target",[IsObject]);
DeclareOperation("Map",[IsObject]);
DeclareOperation("BoundaryMap",[IsObject]);
DeclareOperation("GroupOfResolution",[IsObject]);

## GRAPHS OF GROUPS #################################################
DeclareGlobalFunction("GraphOfGroupsDisplay");
DeclareGlobalFunction("ResolutionGraphOfGroups");
DeclareGlobalFunction("GraphOfGroupsTest");

## LIE ALGEBRAS #####################################################
DeclareGlobalFunction("ChevalleyEilenbergComplex");
DeclareGlobalFunction("IsLieAlgebraHomomorphism");
DeclareGlobalFunction("LieAlgebraHomology");
DeclareGlobalFunction("LowerCentralSeriesLieAlgebra");
DeclareGlobalFunction("LeibnizComplex");
DeclareGlobalFunction("LeibnizAlgebraHomology");

## MANUAL ###########################################################
DeclareGlobalFunction("MakeHAPManual");

## MEAT AXE #########################################################
DeclareGlobalFunction("GeneratorsOfMtxModule");
DeclareGlobalFunction("DesuspensionMtxModule");
DeclareGlobalFunction("FpG_to_MtxModule");

## STREAMS ##########################################################
DeclareGlobalFunction("ChildProcess");
DeclareGlobalFunction("ChildClose");
DeclareGlobalFunction("ChildFunction");
DeclareGlobalFunction("ChildCommand");
DeclareGlobalFunction("ChildRead");
DeclareGlobalFunction("ChildReadEval");
DeclareGlobalFunction("NextAvailableChild");
DeclareGlobalFunction("ParallelList");
DeclareGlobalFunction("ChildPut");
DeclareGlobalFunction("ChildGet");
DeclareGlobalFunction("IsAvailableChild");



## PSEUDOLISTS ######################################################
DeclareGlobalFunction("ListToPseudoList");
DeclareCategory("IsPseudoList",IsObject);
DeclareRepresentation(  "IsPseudoListRep",
                        IsComponentObjectRep,
                        ["elts",
                         "pos" ]);

## TOPOLOGICAL DATA ANALYSIS ########################################
DeclareGlobalFunction("MatrixToChainComplex");
DeclareGlobalFunction("ReadImageAsMatrix");
DeclareGlobalFunction("ContractMatrix");
DeclareGlobalFunction("BettiNumbersOfMatrix");
DeclareGlobalFunction("ThickenedMatrix");
DeclareGlobalFunction("ComplementedMatrix");
DeclareGlobalFunction("WriteImageFile");
DeclareGlobalFunction("ViewMatrix");
DeclareGlobalFunction("DeleteAndContract");
DeclareGlobalFunction("BoundaryOfMatrix");
DeclareGlobalFunction("SingularityMatrix");

## TOPOLOGICAL SPACES ###############################################
DeclareGlobalFunction("MatrixToTopologicalSpace");
DeclareGlobalFunction("ReadImageAsTopologicalSpace");
DeclareGlobalFunction("WriteTopologicalSpaceAsImage");
DeclareGlobalFunction("ViewTopologicalSpace");
DeclareGlobalFunction("ConcatenatedTopologicalSpace");
DeclareGlobalFunction("PathComponent");
DeclareGlobalFunction("ZerothBettiNumber");
DeclareOperation("BettiNumbers",[IsObject,IsInt]);
DeclareGlobalFunction("ThickenedTopologicalSpace");
DeclareGlobalFunction("ComplementTopologicalSpace");
DeclareGlobalFunction("BoundaryTopologicalSpace");
DeclareGlobalFunction("SetT");
DeclareGlobalFunction("EvalT");
DeclareGlobalFunction("ContractTopologicalSpace");
DeclareGlobalFunction("SingularChainComplex");
DeclareGlobalFunction("BoundarySingularities");
