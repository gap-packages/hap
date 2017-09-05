#(C) Graham Ellis, 2005-2006.

DeclareGlobalFunction("EvaluateProperty");
ReadPackage("HAP","lib/Objectifications/types.gd");
ReadPackage("HAP","lib/PureComplexes/complexTypes.gd");
ReadPackage("HAP","lib/GOuterGroups/goutergroup.gd");
ReadPackage("HAP","lib/SimplicialGroups/simpTypes.gd");
ReadPackage("HAP","lib/SimplicialGroups/hapbar.gd");
ReadPackage("HAP","lib/RegularCWComplexes/cwTypes.gd");
ReadPackage("HAP","lib/RegularCWComplexes/cocycle.gd");
ReadPackage("HAP","lib/Sparse/sparse.gd");
ReadPackage("HAP","lib/ArithmeticGroups/arithTypes.gd");
ReadPackage("HAP","lib/TorsionSubcomplexes/torsionsubcomplex.gd");
ReadPackage("HAP","lib/RahmSanchez/DavisComplex.gd");
ReadPackage("HAP","lib/GOuterGroups/functorialGouter.gd");



## FREE G MODULES ###################################################
DeclareGlobalFunction("Negate");
DeclareGlobalFunction("NegateWord");
DeclareGlobalFunction("AlgebraicReduction");
DeclareGlobalFunction("AlgebraicReduction_alt");
DeclareGlobalFunction("AddFreeWords");
DeclareGlobalFunction("AppendFreeWord");
DeclareGlobalFunction("AddFreeWordsModP");
DeclareGlobalFunction("PrintZGword");
DeclareGlobalFunction("TietzeReduction");
DeclareGlobalFunction("MultiplyWord");
DeclareGlobalFunction("WordModP");
DeclareGlobalFunction("OppositeGroup");
DeclareGlobalFunction("QuotientGroup");
DeclareGlobalFunction("ResolutionBoundaryOfWord");

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
DeclareGlobalFunction("GroupAlgebraAsFpGModule");
DeclareGlobalFunction("MaximalSubmodulesOfFpGModule");
DeclareGlobalFunction("MaximalSubmoduleOfFpGModule");
DeclareGlobalFunction("RadicalSeriesOfFpGModule");
DeclareOperation("RadicalSeries",[IsHapFPGModule]);
DeclareOperation("Radical", [IsHapFPGModule]);
DeclareGlobalFunction("CompositionSeriesOfFpGModule");
DeclareGlobalFunction("Classify");
DeclareGlobalFunction("RefineClassification");
DeclareGlobalFunction("FpGModuleSection");


## NONABELIAN TENSOR ################################################
DeclareOperation("HAP_EquivalenceClasses",[IsList,IsFunction]);
DeclareGlobalFunction("NonabelianTensorSquare");
DeclareGlobalFunction("NonabelianTensorSquareAsCatOneGroup");
DeclareGlobalFunction("NonabelianTensorSquareAsCrossedModule");
DeclareGlobalFunction("NonabelianSymmetricSquare");
DeclareGlobalFunction("NonabelianSymmetricSquare_inf");
DeclareGlobalFunction("SymmetricCentre");
DeclareGlobalFunction("NonabelianTensorProduct");
DeclareGlobalFunction("NonabelianTensorProduct_Inf");
DeclareGlobalFunction("NonabelianTensorProduct_alt");
DeclareGlobalFunction("ThirdHomotopyGroupOfSuspensionB");
DeclareGlobalFunction("NonabelianSymmetricKernel");
DeclareGlobalFunction("NonabelianExteriorProduct");
DeclareGlobalFunction("RelativeSchurMultiplier");
DeclareGlobalFunction("EpiCentre");
DeclareGlobalFunction("UpperEpicentralSeries");
DeclareGlobalFunction("BaerInvariant");
DeclareGlobalFunction("TensorCentre");
DeclareGlobalFunction("ThirdHomotopyGroupOfSuspensionB_alt");
DeclareGlobalFunction("NonabelianSymmetricKernel_alt");
DeclareGlobalFunction("NonabelianTensorSquare_inf");
DeclareGlobalFunction("BogomolovMultiplier");
DeclareGlobalFunction("BogomolovMultiplier_viaTensorSquare");
DeclareGlobalFunction("Bogomology");
DeclareGlobalFunction("AreIsoclinic");
DeclareGlobalFunction("PartialIsoclinismClasses");
DeclareGlobalFunction("IsoclinismClasses");
DeclareGlobalFunction("StrongGeneratorsOfDerivedSubgroup");
DeclareGlobalFunction("StrongGeneratorsOfDerivedSubgroup_alt");



## RESOLUTIONS ######################################################
DeclareGlobalFunction("ResolutionArithmeticGroup");
DeclareGlobalFunction("ResolutionGenericGroup");
DeclareOperation("Resolution",[IsGroup,IsInt]);
DeclareGlobalFunction("ResolutionFiniteGroup");
DeclareGlobalFunction("ResolutionSmallFpGroup");
DeclareGlobalFunction("ResolutionSmallGroup");
DeclareGlobalFunction("PresentationOfResolution");
DeclareGlobalFunction("PresentationOfResolution_alt");
DeclareGlobalFunction("ResolutionToResolutionOfFpGroup");

DeclareGlobalFunction("ResolutionFiniteSubgroup");
DeclareGlobalFunction("ResolutionSubgroup");
DeclareGlobalFunction("ResolutionAsphericalPresentation");
DeclareGlobalFunction("ResolutionAbelianGroup");
DeclareGlobalFunction("ResolutionAlmostCrystalGroup");
DeclareGlobalFunction("RelativeCentralQuotientSpaceGroup");
DeclareGlobalFunction("ResolutionAlmostCrystalQuotient");
DeclareGlobalFunction("CayleyGraphOfGroupDisplay");
DeclareGlobalFunction("CayleyGraphOfGroup");
DeclareGlobalFunction("TietzeReducedResolution");
DeclareGlobalFunction("HAPTietzeReduction_OneLevel");
DeclareGlobalFunction("HAPTietzeReduction_OneStep");
DeclareGlobalFunction("HAPTietzeReduction_Inf");
#DeclareGlobalFunction("TietzeReducedResolution_alt");
DeclareGlobalFunction("RecalculateIncidenceNumbers");


## RESOLUTIONS MOD P ################################################
DeclareGlobalFunction("ResolutionPrimePowerGroup");
DeclareGlobalFunction("ResolutionPrimePowerGroupSparse");
DeclareGlobalFunction("RankPrimeHomology");
DeclareGlobalFunction("RankHomologyPGroup");
DeclareGlobalFunction("NumberGeneratorsOfGroupHomology");
DeclareGlobalFunction("PoincareSeries");
DeclareGlobalFunction("PoincareSeries_alt");
DeclareGlobalFunction("PoincareSeriesApproximation");
DeclareGlobalFunction("PoincareSeriesPrimePart");
DeclareGlobalFunction("ExpansionOfRationalFunction");
DeclareGlobalFunction("EfficientNormalSubgroups");
DeclareGlobalFunction("RadicalSeriesOfResolution");


## FUNCTORS #########################################################
#DeclareGlobalFunction("EvaluateProperty");
DeclareGlobalFunction("EvenSubgroup");
DeclareGlobalFunction("EquivariantChainMap");
DeclareGlobalFunction("ModularEquivariantChainMap");
DeclareGlobalFunction("TensorWithIntegers");
DeclareGlobalFunction("TensorWithIntegersSparse");
DeclareGlobalFunction("FilteredTensorWithIntegers");
DeclareGlobalFunction("FilteredTensorWithIntegersModP");
DeclareGlobalFunction("TensorWithIntegersModP");# partial doc
DeclareGlobalFunction("TensorWithIntegersModPSparse");
DeclareGlobalFunction("TensorWithTwistedIntegers");
DeclareGlobalFunction("TensorWithTwistedIntegersModP");
DeclareGlobalFunction("PrimePartDerivedFunctor");
DeclareGlobalFunction("PrimePartDerivedFunctorViaSubgroupChain");
DeclareGlobalFunction("ReduceGenerators");
DeclareGlobalFunction("ReduceGenerators_alt");
DeclareGlobalFunction("HomToIntegers"); #partial doc
DeclareGlobalFunction("HomToInt_ChainComplex");#<HomToIntegers
DeclareGlobalFunction("HomToInt_ChainMap");#<HomToIntegers
DeclareGlobalFunction("HomToIntegralModule");
DeclareGlobalFunction("TensorWithIntegralModule");
DeclareGlobalFunction("PermToMatrixGroup");
DeclareGlobalFunction("AbelianInvariantsToTorsionCoefficients");
DeclareGlobalFunction("TorsionGeneratorsAbelianGroup");
DeclareGlobalFunction("TensorWithRationals");
DeclareGlobalFunction("BigStepLCS");
DeclareGlobalFunction("HomToIntegersModP");
DeclareGlobalFunction("CoClass");
DeclareGlobalFunction("BoundaryMatrix");
DeclareGlobalFunction("Prank");
DeclareGlobalFunction("PrankAlt");
DeclareGlobalFunction("PCentre");
DeclareGlobalFunction("PUpperCentralSeries");
DeclareOperation("Compose",[IsGroupHomomorphism,IsGroupHomomorphism]);
DeclareGlobalFunction("CanonicalRightCountableCosetElement");
DeclareProperty("IsHAPRationalMatrixGroup",IsMatrixGroup);
DeclareProperty("IsHAPRationalSpecialLinearGroup",IsMatrixGroup);
DeclareGlobalFunction("SL2Z");
DeclareGlobalFunction("KernelWG");
DeclareGlobalFunction("ScatterPlot");#doc
DeclareGlobalFunction("PushoutOfFpGroups");#<coproduct
DeclareOperation("Pushout",[IsGroupHomomorphism,IsGroupHomomorphism]);
DeclareGlobalFunction("Fp2PcpAbelianGroupHomomorphism");
DeclareGlobalFunction("IsIsomorphismOfAbelianFpGroups");
DeclareGlobalFunction("SignedPermutationGroup");



## PERTURBATIONS ####################################################
DeclareGlobalFunction("TwistedTensorProduct");
DeclareGlobalFunction("ResolutionFiniteExtension");
DeclareGlobalFunction("ResolutionNormalSeries");
DeclareGlobalFunction("ResolutionDirectProduct");
DeclareGlobalFunction("ResolutionExtension");
DeclareGlobalFunction("ResolutionFiniteDirectProduct");
DeclareGlobalFunction("ResolutionSubnormalSeries");
DeclareGlobalFunction("FreeGResolution");
DeclareGlobalFunction("ContractibleGcomplex");
DeclareGlobalFunction("ContractibleSL2ZComplex");
DeclareGlobalFunction("QuotientOfContractibleGcomplex");
DeclareGlobalFunction("ExtendScalars");
DeclareGlobalFunction("InduceScalars");
DeclareGlobalFunction("TwistedResolution");
DeclareGlobalFunction("CoxeterComplex");
DeclareGlobalFunction("CoxeterComplex_alt");
DeclareGlobalFunction("ResolutionCoxeterGroup");
DeclareGlobalFunction("CyclesOfFilteredChainComplex");
DeclareGlobalFunction("BoundariesOfFilteredChainComplex");
#DeclareGlobalFunction("ResolutionGTree");
DeclareGlobalFunction("ResolutionSL2Z");
DeclareGlobalFunction("ResolutionSL2Z_alt");


## ARTIN COXETER ####################################################
DeclareGlobalFunction("CoxeterDiagramMatrix");
DeclareGlobalFunction("CoxeterDiagramVertices");
DeclareGlobalFunction("CoxeterDiagramFpArtinGroup");
DeclareGlobalFunction("CoxeterDiagramFpCoxeterGroup");
DeclareGlobalFunction("CoxeterDiagramMatCoxeterGroup");
DeclareGlobalFunction("CoxeterSubDiagram");
#DeclareOperation("CoxeterMatrix",[IsList]);
DeclareAttribute("CoxeterMatrix",IsList);
DeclareGlobalFunction("CoxeterDiagramComponents");
DeclareGlobalFunction("CoxeterDiagramDegree");
DeclareGlobalFunction("CoxeterDiagramIsSpherical");
DeclareGlobalFunction("ResolutionArtinGroup");
DeclareGlobalFunction("CoxeterDiagramDisplay");
DeclareGlobalFunction("NoncrossingPartitionsLatticeDisplay");
DeclareGlobalFunction("CoxeterWythoffComplex");

## HOMOLOGY #########################################################
DeclareOperation("Homology",[IsHapChain,IsInt]);
DeclareOperation("Homology",[IsHapGComplex,IsInt]);
DeclareOperation("PathComponents",[IsMagma]);
#DeclareOperation("Homology",[IsObject,IsInt]);
DeclareOperation("PersistentHomology",[IsList,IsInt, IsInt]);
DeclareGlobalFunction("BarCode");#doc
DeclareGlobalFunction("BarCodeDisplay");#doc
DeclareGlobalFunction("BarCodeCompactDisplay");#doc
DeclareGlobalFunction("PersistentHomologyOfPureCubicalComplex");
DeclareGlobalFunction("PersistentHomologyOfPureCubicalComplex_Alt");
DeclareGlobalFunction("ZZPersistentHomologyOfPureCubicalComplex");
DeclareGlobalFunction("PersistentHomologyOfQuotientGroupSeries");
DeclareGlobalFunction("PersistentCohomologyOfQuotientGroupSeries");
DeclareGlobalFunction("PersistentHomologyOfFilteredPureCubicalComplex");
DeclareGlobalFunction("PersistentHomologyOfFilteredPureCubicalComplex_alt");
DeclareGlobalFunction("PersistentHomologyOfFilteredChainComplex");
DeclareGlobalFunction("NormalSeriesToQuotientHomomorphisms");
DeclareGlobalFunction("LinearHomomorphismsPersistenceMat");
DeclareGlobalFunction("LinearHomomorphismsZZPersistenceMat");
DeclareGlobalFunction("PersistentHomologyOfQuotientGroupSeries_Int");
DeclareGlobalFunction("PersistentHomologyOfSubGroupSeries");
DeclareGlobalFunction("TruncatedGComplex");
DeclareGlobalFunction("UniversalBarCode");
DeclareGlobalFunction("UniversalBarCodeEval");
DeclareGlobalFunction("HomologyVectorSpace");
DeclareGlobalFunction("IntegralHomology");#<Homology
DeclareGlobalFunction("ModularHomology");#<Homology
DeclareGlobalFunction("GroupHomology");
DeclareGlobalFunction("RelativeGroupHomology");
DeclareGlobalFunction("RipsHomology");
DeclareGlobalFunction("IntegralCohomology");#<Cohomology
#DeclareOperation("Cohomology",[IsObject,IsObject]);
DeclareOperation("Cohomology",[IsHapCochain,IsInt]);#doc
DeclareOperation("Cohomology",[IsHapGCocomplex,IsInt]);
DeclareOperation("CupProduct",[IsHapRegularCWComplex]);
DeclareGlobalFunction("HAP_CupProductOfPresentation");
#DeclareGlobalFunction("Cohomology");
DeclareGlobalFunction("Syzygy");
DeclareGlobalFunction("CR_IntegralCycleToClass");
DeclareGlobalFunction("CocycleCondition");
DeclareGlobalFunction("StandardCocycle");
DeclareGlobalFunction("IsSuperperfect");
DeclareGlobalFunction("ModularCohomology");#<Homology
DeclareOperation("SolutionsMatDestructive",
			[IsOrdinaryMatrix,IsOrdinaryMatrix]);
DeclareGlobalFunction("HomologyPb");
DeclareGlobalFunction("HomologyPbs");
DeclareGlobalFunction("HomologyPrimePart");
DeclareGlobalFunction("CohomologyPrimePart");
DeclareGlobalFunction("GroupCohomology");
DeclareGlobalFunction("IntegralHomologyOfChainComplex");#<Homology
DeclareGlobalFunction("IntegralCohomologyOfCochainComplex");#<Cohomology
DeclareGlobalFunction("LefschetzNumberOfChainMap");
DeclareOperation("LefschetzNumber",[IsObject]);
DeclareGlobalFunction("HomologyOfPureCubicalComplex");#<Homology

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


## CURVATURE & POLYTOPES ############################################
DeclareGlobalFunction("IsAspherical");
DeclareAttribute("StarGraph",IsFpGroup);
DeclareGlobalFunction("PolytopalGenerators");
DeclareGlobalFunction("VectorStabilizer");
DeclareGlobalFunction("PolytopalComplex");#doc
DeclareGlobalFunction("OrbitPolytope");#doc
DeclareGlobalFunction("RegularCWOrbitPolytope");#doc
DeclareGlobalFunction("HAPRegularCWPolytope");#<RegularCWPolytope
DeclareOperation("RegularCWPolytope",[IsList]);#doc
DeclareGlobalFunction("PolymakeFaceLattice");
DeclareGlobalFunction("HAP_TzPair");
DeclareGlobalFunction("HAP_AddGenerator");
DeclareGlobalFunction("SmoothedFpGroup");


## POLYCYLIC ########################################################
DeclareGlobalFunction("ResolutionAbelianPcpGroup");
DeclareGlobalFunction("ResolutionNilpotentGroup");

## OBJECTIFICATION ##################################################
#DeclareOperation("Target",[IsObject]);
DeclareAttribute("Target",IsObject);
DeclareOperation("Map",[IsObject]);
DeclareOperation("BoundaryMap",[IsHapRegularCWComplex]);
DeclareOperation("GroupOfResolution",[IsObject]);

## GRAPHS OF GROUPS #################################################
DeclareGlobalFunction("GraphOfGroupsDisplay");
DeclareGlobalFunction("ResolutionGraphOfGroups");
DeclareGlobalFunction("GraphOfGroupsTest");
DeclareGlobalFunction("GraphOfResolutionsTest");
DeclareGlobalFunction("GraphOfResolutionsToGroups");
DeclareGlobalFunction("GraphOfResolutions");
DeclareGlobalFunction("GraphOfResolutionsDisplay");
DeclareGlobalFunction("TreeOfResolutionsToContractibleGcomplex");
DeclareGlobalFunction("TreeOfGroupsToContractibleGcomplex");
DeclareGlobalFunction("ConjugatedResolution");


## LIE ALGEBRAS #####################################################
DeclareGlobalFunction("ChevalleyEilenbergComplex");
DeclareGlobalFunction("IsLieAlgebraHomomorphism");
DeclareGlobalFunction("LieAlgebraHomology");
DeclareGlobalFunction("LowerCentralSeriesLieAlgebra");
DeclareGlobalFunction("LeibnizComplex");
DeclareGlobalFunction("LeibnizAlgebraHomology");
DeclareGlobalFunction("LieTensorSquare");
DeclareGlobalFunction("LieCoveringHomomorphism");
DeclareGlobalFunction("LeibnizQuasiCoveringHomomorphism");
DeclareGlobalFunction("LieExteriorSquare");
DeclareGlobalFunction("LieEpiCentre");
DeclareGlobalFunction("LieTensorCentre");

## MANUAL ###########################################################
DeclareGlobalFunction("MakeHAPManual");

## MEAT AXE #########################################################
DeclareGlobalFunction("GeneratorsOfMtxModule");
DeclareGlobalFunction("DesuspensionMtxModule");
DeclareGlobalFunction("FpG_to_MtxModule");

## STREAMS ##########################################################
DeclareGlobalFunction("ChildProcess");
DeclareGlobalFunction("ChildClose");
DeclareGlobalFunction("ChildRestart");
DeclareGlobalFunction("ChildFunction");
DeclareGlobalFunction("ChildCommand");
DeclareGlobalFunction("ChildRead");
DeclareGlobalFunction("ChildReadEval");
DeclareGlobalFunction("NextAvailableChild");
DeclareGlobalFunction("ParallelList");
DeclareGlobalFunction("ChildPut");
DeclareGlobalFunction("ChildGet");
DeclareGlobalFunction("IsAvailableChild");
DeclareOperation("HAPPrintTo",[IsString,IsObject]);
DeclareOperation("HAPRead",[IsString]);



## PSEUDOLISTS ######################################################
DeclareGlobalFunction("ListToPseudoList");
DeclareCategory("IsPseudoList",IsObject);
DeclareRepresentation(  "IsPseudoListRep",
                        IsComponentObjectRep,
                        ["elts",
                         "pos" ]);
DeclareProperty("IsPseudoListWithFunction",IsPseudoList);
DeclareGlobalFunction("LazyList");


## PURE COMPLEXES  ###############################################
DeclareOperation("Dimensions",[IsHapPureCubicalComplex]);
DeclareGlobalFunction("PureComplex");
DeclareGlobalFunction("IsPureComplex");
DeclareGlobalFunction("PureComplexToSimplicialComplex");#<Nerve
DeclareGlobalFunction("UnitBall");
DeclareGlobalFunction("UnitCubicalBall");
DeclareGlobalFunction("UnitPermutahedralBall");
DeclareGlobalFunction("ThickenedPureComplex");#<PureComplexThickened
DeclareGlobalFunction("ComplementOfPureComplex");#<PureComplexComplement
DeclareGlobalFunction("PureComplexUnion");#doc
DeclareGlobalFunction("PureComplexIntersection");#doc
DeclareGlobalFunction("PureComplexDifference");#doc
DeclareOperation("PureComplexSubcomplex",[IsHapPureCubicalComplex,IsList]);
DeclareGlobalFunction("HAP_PureComplexSubcomplex");
DeclareOperation("PureComplexMeet",[IsHapPureCubicalComplex,IsHapPureCubicalComplex]);
DeclareGlobalFunction("BoundaryOfPureComplex");#<PureComplexBoundary

## POLYTOPAL COMPLEXES  ###############################################
DeclareOperation("SimplicialComplex",[IsList]);#doc
DeclareGlobalFunction("CubicalToPermutahedralArray");
DeclareGlobalFunction("PermutahedralToCubicalArray");
DeclareGlobalFunction("ContractibleSubcomplexOfPureCubicalComplex");
DeclareGlobalFunction("ContractibleSubcomplexOfSimplicialComplex");
DeclareOperation("ContractibleSubcomplex",[IsHapPureCubicalComplex]);
DeclareGlobalFunction("AcyclicSubcomplexOfPureCubicalComplex");
DeclareGlobalFunction("IntegerSimplicialComplex");
DeclareGlobalFunction("HomotopyEquivalentMaximalPureCubicalSubcomplex");
DeclareGlobalFunction("HomotopyEquivalentMaximalPureSubcomplex");
DeclareGlobalFunction("HomotopyEquivalentMinimalPureCubicalSubcomplex");
DeclareGlobalFunction("HomotopyEquivalentMinimalPureSubcomplex");
DeclareAttribute("EulerCharacteristic",IsHapPureCubicalComplex);
DeclareAttribute("EulerCharacteristic",IsHapCubicalComplex);
DeclareAttribute("EulerCharacteristic",IsHapSimplicialComplex);
DeclareOperation("ContractedComplex",[IsObject]);#doc
DeclareGlobalFunction("ReadImageAsPureCubicalComplex");#doc
DeclareGlobalFunction("ReadImageAsFilteredPureCubicalComplex");#doc
DeclareGlobalFunction("ReadLinkImageAsPureCubicalComplex");
DeclareGlobalFunction("ReadMatrixAsPureCubicalComplex");
DeclareGlobalFunction("ReadImageSequenceAsPureCubicalComplex");
DeclareGlobalFunction("WritePureCubicalComplexAsImage");
DeclareGlobalFunction("ViewPureCubicalComplex");#<Display
DeclareGlobalFunction("ViewPureComplex");#<Display
DeclareGlobalFunction("View3dPureComplex");#<Display
DeclareGlobalFunction("PureCubicalComplex");#doc
DeclareGlobalFunction("CubicalComplex");#doc
DeclareGlobalFunction("PermutahedralComplexToRegularCWComplex");#<RegularCWComplex
DeclareGlobalFunction("PurePermutahedralComplex");#doc
DeclareGlobalFunction("PureCubicalComplexUnion");#<PureComplexUnion
DeclareGlobalFunction("PureCubicalComplexDifference");#<PureComplexDifference
DeclareGlobalFunction("PureCubicalComplexIntersection");#<PureComplexIntersection
DeclareGlobalFunction("PureCubicalComplexToCubicalComplex");
DeclareGlobalFunction("FilteredPureCubicalComplexToCubicalComplex");
DeclareGlobalFunction("ConcentricallyFilteredPureCubicalComplex");#<ConcentricFiltration
DeclareGlobalFunction("ContractedFilteredPureCubicalComplex");#<ContractedComplex
DeclareOperation("ChainComplex",[IsObject]);#doc
DeclareOperation("ChainMap",[IsHapRegularCWMap]);#doc
DeclareOperation("CochainComplex",[IsObject]);
DeclareOperation("CoboundaryMatrix",[IsHapCochainComplex,IsInt]);
DeclareOperation("ChainComplex",[IsObject, IsBool]);
DeclareOperation("ChainComplexOfPair",[IsObject,IsObject]);
DeclareOperation("SparseChainComplexOfPair",[IsObject,IsObject]);
DeclareGlobalFunction("ChainComplexOfCubicalComplex");#<ChainComplex
DeclareGlobalFunction("SparseChainComplexOfCubicalComplex");
DeclareGlobalFunction("SparseFilteredChainComplexOfFilteredCubicalComplex");
DeclareGlobalFunction("ChainComplexOfCubicalPair");
DeclareGlobalFunction("SparseChainComplexOfCubicalPair");
DeclareGlobalFunction("HAP_PureCubicalPairToCWMap");
DeclareGlobalFunction("HAP_SimplicialPairToCWMap");
DeclareOperation("ExcisedPair",[IsHapPureCubicalComplex,IsHapPureCubicalComplex]);#doc
DeclareGlobalFunction("ExcisedPureCubicalPair");#<ExcisedPair
DeclareGlobalFunction("ExcisedPureCubicalPair_dim_2");#<ExcisedPair
DeclareGlobalFunction("ChainComplexOfSimplicialPair");
DeclareGlobalFunction("ChainMapOfCubicalPairs");#<ChainMap
DeclareGlobalFunction("SparseChainMapOfCubicalPairs");
DeclareGlobalFunction("ChainComplexOfSimplicialComplex");
DeclareGlobalFunction("SparseChainComplexOfSimplicialComplex");
DeclareGlobalFunction("SparseFilteredChainComplexOfFilteredSimplicialComplex");
DeclareGlobalFunction("ChainMapOfSimplicialMap");#<ChainMap
DeclareGlobalFunction("SkeletonOfSimplicialComplex");
DeclareGlobalFunction("CechComplexOfPureCubicalComplex");
DeclareOperation("Nerve",[IsHapPureCubicalComplex]);#doc
DeclareOperation("Nerve",[IsHapPureCubicalComplex,IsInt]);#doc
DeclareOperation("Nerve",[IsHapPurePermutahedralComplex]);#doc
DeclareOperation("Nerve",[IsHapPurePermutahedralComplex,IsInt]);#doc
DeclareGlobalFunction("QuillenComplex");#doc
DeclareGlobalFunction("GChainComplex");
DeclareGlobalFunction("SimplicialMap");
DeclareGlobalFunction("SimplicialMapNC");
DeclareGlobalFunction("MaximalSimplicesToSimplicialComplex");#<SimplicialComplex
DeclareGlobalFunction("SimplicesToSimplicialComplex");#<SimplicialComplex
DeclareGlobalFunction("MaximalSimplicesOfSimplicialComplex");
DeclareGlobalFunction("ContractSimplicialComplex");#<ContractedSimplicialComplex
DeclareGlobalFunction("PathComponentOfPureCubicalComplex");
DeclareGlobalFunction("PathComponentOfPureComplex");
DeclareOperation("PathComponent",[IsHapPureCubicalComplex,IsInt]);
DeclareOperation("BettiNumber",[IsHapRegularCWComplex,IsInt]);#doc
DeclareOperation("PersistentBettiNumbers",[IsHapFilteredRegularCWComplex,IsInt]);#doc
DeclareOperation("Bettinumbers",[IsObject,IsInt]);
DeclareGlobalFunction("BettinumbersOfPureCubicalComplex_dim_2");#<BettiNumber
DeclareGlobalFunction("ThickenedPureCubicalComplex");#<PureComplexThickened
DeclareGlobalFunction("ThickenedHEPureCubicalComplex");
DeclareGlobalFunction("ThickenedPureCubicalComplex_dim2");#<PureComplexThickened
DeclareGlobalFunction("BoundaryOfPureCubicalComplex");
DeclareGlobalFunction("ContractPureCubicalComplex");#<ContractedComplex
DeclareGlobalFunction("ContractPureComplex");#<ContractedComplex
DeclareGlobalFunction("ZigZagContractedPureCubicalComplex");#<ZigZagContractedComplex
DeclareGlobalFunction("ZigZagContractedPureComplex");#<ZigZagContractedComplex
DeclareGlobalFunction("ZigZagContractedFilteredPureCubicalComplex");#<ZigZagContractedComplex
DeclareOperation("ZigZagContractedComplex",[IsHapPureCubicalComplex]);#doc
DeclareGlobalFunction("CropPureCubicalComplex");
DeclareGlobalFunction("CropPureComplex");
DeclareGlobalFunction("ComplementOfPureCubicalComplex");#<PureComplexComplement
DeclareGlobalFunction("SingularitiesOfPureCubicalComplex");
DeclareGlobalFunction("ChainInclusionOfPureCubicalPair");
DeclareGlobalFunction("DirectProductOfPureCubicalComplexes");#<DirectProduct
DeclareGlobalFunction("PureCubicalComplexToTextFile");
DeclareGlobalFunction("CoreducedChainComplex");
#DeclareGlobalFunction("ReducedChainComplex");
DeclareGlobalFunction("2CoreducedChainComplex");
DeclareGlobalFunction("SuspendedChainComplex");
DeclareGlobalFunction("ReducedSuspendedChainComplex");
DeclareGlobalFunction("RipsChainComplex");
DeclareGlobalFunction("VectorsToSymmetricMatrix");
DeclareGlobalFunction("SymmetricMatrixToIncidenceMatrix");
DeclareGlobalFunction("DensityMat");
DeclareGlobalFunction("IncidenceMatrixToGraph");
DeclareGlobalFunction("SymmetricMatrixToGraph");#doc
DeclareGlobalFunction("SymmetricMatrixToFilteredGraph");#doc
DeclareGlobalFunction("PermGroupToFilteredGraph");
DeclareGlobalFunction("GraphOfSimplicialComplex");#<Graph
DeclareOperation("Graph",[IsHapSimplicialComplex]);#doc
DeclareGlobalFunction("PathComponentsOfGraph");
DeclareGlobalFunction("PathComponentsOfSimplicialComplex");
DeclareGlobalVariable("PathComponentOfSimplicialComplex");
DeclareGlobalFunction("PathComponentsOfSimplicialComplex_alt");
DeclareGlobalFunction("ContractGraph");#<ContractedComplex
DeclareGlobalFunction("SimplicialNerveOfGraph");#<CliqueComplex
DeclareGlobalFunction("SimplicialNerveOfTwoComplex");#<CliqueComplex
DeclareOperation("CliqueComplex",[IsHapSimplicialComplex, IsInt]);#doc
DeclareGlobalFunction("SimplicialNerveOfFilteredGraph");#<CliqueComplex
DeclareGlobalFunction("GraphDisplay");#<Display
DeclareGlobalVariable("ViewGraph");#<Display
DeclareGlobalFunction("SymmetricMatDisplay");
DeclareGlobalFunction("SkeletonOfCubicalComplex");
DeclareGlobalFunction("MorseFiltration");
DeclareGlobalFunction("ContractCubicalComplex_dim2");#<ContractedComplex
DeclareGlobalFunction("ContractCubicalComplex_dim3");#<ContractedComplex
DeclareGlobalFunction("ContractCubicalComplex");#<ContractedComplex
DeclareGlobalFunction("DVFReducedCubicalComplex");#<ContractedComplex
DeclareGlobalFunction("BoundingPureCubicalComplex");
DeclareGlobalFunction("BoundingPureComplex");
DeclareGlobalFunction("SuspensionOfPureCubicalComplex");
DeclareGlobalFunction("ThickeningFiltration");#doc
DeclareOperation("ConcentricFiltration",[IsHapPureCubicalComplex,IsInt]);#doc
DeclareGlobalFunction("Dendrogram");
DeclareOperation("FiltrationTerm",[IsHapPureCubicalComplex,IsInt]);
DeclareGlobalFunction("FiltrationTermOfPureCubicalComplex");
DeclareGlobalFunction("FiltrationTermOfRegularCWComplex");
DeclareGlobalFunction("FiltrationTermOfGraph");
DeclareGlobalFunction("DisplayDendrogram");#doc
DeclareGlobalFunction("DisplayDendrogramMat");#doc
DeclareGlobalFunction("DendrogramMat");#doc
DeclareGlobalFunction("ComplementOfFilteredPureCubicalComplex");
DeclareGlobalFunction("DendrogramToPersistenceMat");
DeclareGlobalFunction("BarCodeOfSymmetricMatrix");
DeclareGlobalFunction("BarCodeOfFilteredPureCubicalComplex");
DeclareGlobalFunction("HenonOrbit");#doc
DeclareGlobalFunction("RandomCubeOfPureCubicalComplex");
DeclareGlobalFunction("RandomCellOfPureComplex");
DeclareOperation("PureComplexRandomCell",[IsHapPureCubicalComplex]);
DeclareGlobalFunction("RandomSimplicialGraph");
DeclareGlobalFunction("RandomSimplicialTwoComplex");
DeclareGlobalFunction("FirstHomologySimplicialTwoComplex");
DeclareGlobalFunction("Mapper");
DeclareGlobalFunction("Mapper_alt");
DeclareGlobalFunction("VectorsToOneSkeleton");


########################## ARRAYS ################################
DeclareGlobalFunction("ArrayValue");
DeclareGlobalFunction("ArrayValueFunctions");
DeclareGlobalFunction("ArrayValueKD");
DeclareGlobalFunction("ArrayToPureCubicalComplex");
DeclareGlobalFunction("FrameArray");
DeclareGlobalFunction("FramedPureCubicalComplex");
DeclareGlobalFunction("UnframeArray");
DeclareGlobalFunction("PermuteArray");
DeclareGlobalFunction("ArraySum");
DeclareGlobalFunction("ArrayDimension");
DeclareGlobalFunction("ArrayDimensions");
DeclareGlobalFunction("ContractArray");
DeclareGlobalFunction("ContractPermArray");
DeclareGlobalFunction("HAP_Binlisttoint");
DeclareGlobalFunction("HAP_PermBinlisttoint");
DeclareGlobalFunction("ContractMatrix");
DeclareGlobalFunction("ContractPermMatrix");
DeclareGlobalFunction("ContractibleSubArray");
DeclareGlobalFunction("HomotopyEquivalentLargerSubArray");
DeclareGlobalFunction("HomotopyEquivalentLargerSubPermArray");
DeclareGlobalFunction("HomotopyEquivalentLargerSubArray3D");
DeclareGlobalFunction("HomotopyEquivalentLargerSubPermArray3D");
DeclareGlobalFunction("HomotopyEquivalentSmallerSubArray");
DeclareGlobalFunction("HomotopyEquivalentSmallerSubPermArray");
DeclareGlobalFunction("HomotopyEquivalentSmallerSubArray3D");
DeclareGlobalFunction("HomotopyEquivalentSmallerSubPermArray3D");
DeclareGlobalFunction("ContractibleSubMatrix");
DeclareGlobalFunction("HomotopyEquivalentLargerSubMatrix");
DeclareGlobalFunction("HomotopyEquivalentLargerSubPermMatrix");
DeclareGlobalFunction("HomotopyEquivalentSmallerSubMatrix");
DeclareGlobalFunction("HomotopyEquivalentSmallerSubPermMatrix");
DeclareGlobalFunction("ArrayAssign");
DeclareGlobalFunction("UnboundedArrayAssign");
DeclareGlobalFunction("ArrayAssignFunctions");
DeclareGlobalFunction("ArrayIterate");
DeclareGlobalFunction("ArrayIterateBreak");
DeclareGlobalFunction("IsContractibleCube_higherdims");
DeclareGlobalFunction("Array");
DeclareGlobalFunction("BinaryArrayToTextFile");


## CAT ONE GROUPS ###################################################
DeclareGlobalFunction("XmodToHAP");
DeclareGlobalFunction("AutomorphismGroupAsCatOneGroup");
DeclareGlobalFunction("GModuleAsCatOneGroup");
DeclareOperation("HomotopyGroup",[IsHapCatOneGroup,IsInt]);
DeclareOperation("HomotopyGroup",[IsHapSimplicialGroup,IsInt]);
DeclareOperation("HomotopyModule",[IsHapCatOneGroup,IsInt]);
DeclareOperation("MooreComplex",[IsObject]);
DeclareGlobalFunction("HasTrivialPostnikovInvariant");
DeclareGlobalFunction("IdentityAmongRelatorsDisplay");
DeclareGlobalFunction("IdentityAmongRelators");
DeclareGlobalFunction("NormalSubgroupAsCatOneGroup");
DeclareGlobalFunction("QuasiIsomorph");
#DeclareGlobalFunction("IsQuasiMorphismCatOneGroups");
DeclareGlobalFunction("QuotientQuasiIsomorph");
DeclareGlobalFunction("SubQuasiIsomorph");
DeclareGlobalFunction("SylowSubgroupOfCatOneGroup");
DeclareGlobalFunction("CrossedInvariant");
AutomorphismGroupAsCrossedModule:=CrossedModuleByAutomorphismGroup;
MakeReadOnlyGlobal("AutomorphismGroupAsCrossedModule");

## COMMUTATIVE DIAGRAMS #############################################
DeclareGlobalFunction("HomomorphismChainToCommutativeDiagram");
DeclareGlobalFunction("NerveOfCommutativeDiagram");
DeclareGlobalFunction("GroupHomologyOfCommutativeDiagram");
DeclareGlobalFunction("PersistentHomologyOfCommutativeDiagramOfPGroups");
DeclareGlobalFunction("NormalSeriesToQuotientDiagram");

## REGULAR CW-COMPLEXES ################################################
DeclareGlobalFunction("ReadImageAsWeightFunction");
DeclareGlobalFunction("EulerIntegral");
DeclareOperation("RegularCWComplex",[IsList]);#doc
DeclareOperation("FilteredRegularCWComplex",[IsHapFilteredSimplicialComplex]);
DeclareOperation("RegularCWMap",[IsHapPureCubicalComplex,IsHapPureCubicalComplex]);
DeclareGlobalFunction("HAPRegularCWComplex");#<RegularCWComplex
DeclareGlobalFunction("OrientRegularCWComplex");
DeclareGlobalFunction("SimplifiedRegularCWComplex");#<SimplifiedComplex
DeclareOperation("SimplifiedComplex",[IsHapRegularCWComplex]);#doc
DeclareGlobalFunction("ContractedRegularCWComplex");#<ContractedComplex
DeclareGlobalFunction("DeformationRetract");
DeclareGlobalFunction("SimplicialComplexToRegularCWComplex");#<RegularCWComplex
DeclareGlobalFunction("GraphOfRegularCWComplex");#Not yet implemented
DeclareGlobalFunction("HomotopyGraph");
DeclareGlobalFunction("CubicalComplexToRegularCWComplex");#<RegularCWComplex
DeclareGlobalFunction("HAPContractRegularCWComplex");
DeclareGlobalFunction("HAPCocontractRegularCWComplex");
DeclareGlobalFunction("HAPRemoveCellFromRegularCWComplex");
DeclareGlobalFunction("CriticalCellsOfRegularCWComplex");#<CriticalCells
DeclareOperation("CriticalCells",[IsHapRegularCWComplex]);#doc
DeclareGlobalFunction("CocriticalCellsOfRegularCWComplex");
DeclareGlobalFunction("HAPContractRegularCWComplex_Alt");
DeclareGlobalFunction("ChainComplexOfRegularCWComplex");
DeclareGlobalFunction("ChainComplexOfRegularCWComplexWithVectorField");
DeclareGlobalFunction("FundamentalGroupOfRegularCWComplex");
DeclareGlobalFunction("FundamentalGroupOfRegularCWMap");
DeclareGlobalFunction("FundamentalGroupSimplicialTwoComplex");
DeclareOperation("FundamentalGroup",[IsHapRegularCWComplex]);
DeclareOperation("FundamentalGroup",[IsHapRegularCWComplex,IsInt]);
DeclareOperation("FundamentalGroupOfQuotient",[IsHapEquivariantCWComplex]);
DeclareOperation("ChainComplexOfQuotient",[IsHapEquivariantCWComplex]);
DeclareGlobalFunction("RestrictedEquivariantCWComplex");
DeclareGlobalFunction("ResolutionAffineCrystGroup");
DeclareGlobalFunction("ResolutionToEquivariantCWComplex");
DeclareGlobalFunction("EquivariantCWComplexToResolution");
DeclareGlobalFunction("EquivariantEuclideanSpace");
DeclareGlobalFunction("EquivariantOrbitPolytope");
DeclareGlobalFunction("EquivariantTwoComplex");
DeclareGlobalFunction("HAPRemoveVectorField");
DeclareGlobalFunction("IsPureRegularCWComplex");
DeclareGlobalFunction("BoundaryOfPureRegularCWComplex");
DeclareGlobalFunction("BoundaryPairOfPureRegularCWComplex");
DeclareGlobalFunction("HAP_Sequence2Boundaries");
DeclareOperation("PiZero",[IsHapRegularCWComplex]);#doc
DeclareOperation("PiZero",[IsHapSimplicialComplex]);#doc
DeclareOperation("PiZero",[IsHapGraph]);#doc
DeclareGlobalFunction("PiZeroOfRegularCWComplex");#doc
DeclareGlobalFunction("TruncatedRegularCWComplex");
DeclareGlobalFunction("HomotopyTruncation");
DeclareGlobalFunction("VerticesOfRegularCWCell");
DeclareGlobalFunction("BoundaryOfRegularCWCell");
DeclareGlobalFunction("FilteredCubicalComplexToFilteredRegularCWComplex");
DeclareGlobalFunction("StructuralCopyOfFilteredRegularCWComplex");
DeclareGlobalFunction("SparseChainComplexOfFilteredRegularCWComplex");
DeclareGlobalFunction("HAPContractFilteredRegularCWComplex");
DeclareGlobalFunction("ContractedFilteredRegularCWComplex");
DeclareGlobalFunction("DirectProductOfRegularCWComplexes");
DeclareGlobalFunction("DiagonalApproximation");
DeclareGlobalFunction("CWMap2ChainMap");
DeclareGlobalFunction("NonRegularCWBoundary");

## KNOTS ############################################################
DeclareGlobalFunction("PureCubicalKnot");#doc
DeclareGlobalFunction("GaussCodeOfPureCubicalKnot");
DeclareOperation("DisplayArcPresentation",[IsHapPureCubicalComplex]);#doc
DeclareGlobalFunction("WirtingerGroup");
DeclareGlobalFunction("WirtingerGroup_gc");
DeclareGlobalFunction("NumberOfPrimeKnots");
DeclareGlobalFunction("ReflectedCubicalKnot");
DeclareGlobalFunction("ArcPresentation");
DeclareGlobalFunction("PurePermutahedralKnot");#doc
DeclareGlobalFunction("ViewPureCubicalKnot");
DeclareGlobalFunction("KnotGroup");
DeclareGlobalFunction("KnotSum");
DeclareOperation("KnotReflection",[IsHapPureCubicalComplex]);
DeclareGlobalFunction("AlexanderMatrix");
DeclareGlobalFunction("AlexanderPolynomial");
DeclareGlobalFunction("ReadPDBfileAsPureCubicalComplex");#doc
DeclareGlobalFunction("ReadPDBfileAsPurePermutahedralComplex");#doc
DeclareGlobalFunction("ReadCSVfileAsPureCubicalKnot");
DeclareGlobalFunction("DisplayCSVknotFile");
DeclareGlobalFunction("DisplayPDBfile");
DeclareGlobalFunction("ProjectionOfPureCubicalComplex");

## METRICS #########################################################
DeclareOperation("CayleyMetric",[IsPerm,IsPerm,IsInt]);
DeclareOperation("CayleyMetric",[IsPerm,IsPerm]);
DeclareOperation("KendallMetric",[IsPerm,IsPerm,IsInt]);
DeclareOperation("KendallMetric",[IsPerm,IsPerm]);
DeclareOperation("HammingMetric",[IsPerm,IsPerm,IsInt]);
DeclareOperation("HammingMetric",[IsPerm,IsPerm]);
DeclareOperation("EuclideanSquaredMetric",[IsVector,IsVector]);
DeclareOperation("ManhattanMetric",[IsVector,IsVector]);
DeclareOperation("EuclideanApproximatedMetric",[IsVector,IsVector]);
DeclareGlobalFunction("IsMetricMatrix");

## SPARSE ###########################################################
DeclareGlobalFunction("SparseMattoMat");
DeclareGlobalFunction("SparseMat");
DeclareGlobalFunction("SparseRowMult");
DeclareGlobalFunction("SparseRowInterchange");
DeclareGlobalFunction("SparseRowAdd");
DeclareGlobalFunction("SparseSemiEchelon");
DeclareGlobalFunction("SparseRowReduce");
DeclareGlobalFunction("SparseChainComplexOfRegularCWComplex");
DeclareGlobalFunction("SparseChainComplexOfRegularCWComplexWithVectorField");
DeclareGlobalFunction("SparseBoundaryMatrix");
DeclareOperation("SparseChainComplex",[IsHapRegularCWComplex]);
DeclareGlobalFunction("NullspaceSparseMatDestructive");
DeclareGlobalFunction("TransposeOfSparseMat");
DeclareGlobalFunction("ReverseSparseMat");
DeclareGlobalFunction("PersistentHomologyOfFilteredSparseChainComplex");
DeclareGlobalFunction("FilteredChainComplexToFilteredSparseChainComplex");



## OTHER ############################################################
ReadPackage("HAP","lib/CategoryTheory/categories.gd");
ReadPackage("HAP","lib/CategoryTheory/commutativeDiagrams.gd");
ReadPackage("HAP","lib/Knots/knots.gd");
ReadPackage("HAP","lib/HapPrime/derivation.gd");
ReadPackage("HAP","lib/HapPrime/singular.gd");
ReadPackage("HAP","lib/HapPrime/gradedalgebra.gd");
ReadPackage("HAP","lib/HapPrime/polynomials.gd");
ReadPackage("HAP","lib/HapPrime/ringhomomorphism.gd");
ReadPackage("HAP","lib/HapPrime/rings.gd");
ReadPackage("HAP","lib/HapPrime/happrime.gd");
ReadPackage("HAP","lib/ArithmeticGroups/crystTypes.gd");
ReadPackage("HAP","lib/CohomologyOperations/cohoOps.gd");
ReadPackage("HAP","lib/Quandles/quandles.gd");



