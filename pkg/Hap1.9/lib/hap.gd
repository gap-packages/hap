#(C) Graham Ellis, 2005-2006.

DeclareGlobalFunction("EvaluateProperty");
ReadPackage("HAP","lib/Objectifications/types.gd");
ReadPackage("HAP","lib/PolyComplexes/complexTypes.gd");
ReadPackage("HAP","lib/GOuterGroups/goutergroup.gd");
ReadPackage("HAP","lib/SimplicialGroups/simpTypes.gd");
ReadPackage("HAP","lib/RegularCWSpaces/cwTypes.gd");



## FREE G MODULES ###################################################
DeclareGlobalFunction("Negate");
DeclareGlobalFunction("NegateWord");
DeclareGlobalFunction("AlgebraicReduction");
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
DeclareGlobalFunction("CompositionSeriesOfFpGModule");
DeclareGlobalFunction("Classify");
DeclareGlobalFunction("RefineClassification");


## NONABELIAN TENSOR ################################################
DeclareGlobalFunction("NonabelianTensorSquare");
DeclareGlobalFunction("NonabelianSymmetricSquare");
DeclareGlobalFunction("NonabelianSymmetricSquare_inf");
DeclareGlobalFunction("SymmetricCentre");
DeclareGlobalFunction("NonabelianTensorProduct");
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


## RESOLUTIONS ######################################################
DeclareGlobalFunction("ResolutionGenericGroup");
DeclareGlobalFunction("ResolutionFiniteGroup");
DeclareGlobalFunction("ResolutionSmallFpGroup");
DeclareGlobalFunction("PresentationOfResolution");
DeclareGlobalFunction("PresentationOfResolution_alt");
DeclareGlobalFunction("ResolutionToResolutionOfFpGroup");
DeclareGlobalFunction("ResolutionFiniteSubgroup");
DeclareGlobalFunction("ResolutionSubgroup");
DeclareGlobalFunction("ResolutionAsphericalPresentation");
DeclareGlobalFunction("ResolutionAbelianGroup");
DeclareGlobalFunction("ResolutionAlmostCrystalGroup");
DeclareGlobalFunction("ResolutionAlmostCrystalQuotient");
DeclareGlobalFunction("CayleyGraphDisplay");
DeclareGlobalFunction("TietzeReducedResolution");
DeclareGlobalFunction("RecalculateIncidenceNumbers");

## RESOLUTIONS MOD P ################################################
DeclareGlobalFunction("ResolutionPrimePowerGroup");
DeclareGlobalFunction("RankPrimeHomology");
DeclareGlobalFunction("RankHomologyPGroup");
DeclareGlobalFunction("NumberGeneratorsOfGroupHomology");
DeclareGlobalFunction("PoincareSeries");
DeclareGlobalFunction("PoincareSeries_alt");
DeclareGlobalFunction("PoincareSeriesApproximation");
DeclareGlobalFunction("PoincareSeriesPrimePart");
DeclareGlobalFunction("ExpansionOfRationalFunction");
DeclareGlobalFunction("EfficientNormalSubgroups");


## FUNCTORS #########################################################
#DeclareGlobalFunction("EvaluateProperty");
DeclareGlobalFunction("EvenSubgroup");
DeclareGlobalFunction("EquivariantChainMap");
DeclareGlobalFunction("ModularEquivariantChainMap");
DeclareGlobalFunction("TensorWithIntegers");
DeclareGlobalFunction("FilteredTensorWithIntegers");
DeclareGlobalFunction("FilteredTensorWithIntegersModP");
DeclareGlobalFunction("TensorWithIntegersModP");
DeclareGlobalFunction("TensorWithTwistedIntegers");
DeclareGlobalFunction("TensorWithTwistedIntegersModP");
DeclareGlobalFunction("PrimePartDerivedFunctor");
DeclareGlobalFunction("ReduceGenerators");
DeclareGlobalFunction("HomToIntegers");
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
DeclareGlobalFunction("QuotientOfContractibleGcomplex");
DeclareGlobalFunction("ExtendScalars");
DeclareGlobalFunction("InduceScalars");
DeclareGlobalFunction("TwistedResolution");
DeclareGlobalFunction("CoxeterComplex");
DeclareGlobalFunction("ResolutionCoxeterGroup");
DeclareGlobalFunction("CyclesOfFilteredChainComplex");
DeclareGlobalFunction("BoundariesOfFilteredChainComplex");


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
DeclareGlobalFunction("NoncrossingPartitionsLatticeDisplay");
DeclareGlobalFunction("CoxeterWythoffComplex");

## HOMOLOGY #########################################################
DeclareOperation("Homology",[IsHapChain,IsInt]);
DeclareOperation("Homology",[IsHapGComplex,IsInt]);
#DeclareOperation("Homology",[IsObject,IsInt]);
DeclareOperation("PersistentHomology",[IsList,IsInt, IsInt]);
DeclareGlobalFunction("BarCode");
DeclareGlobalFunction("BarCodeDisplay");
DeclareGlobalFunction("PersistentHomologyOfPureCubicalComplex");
DeclareGlobalFunction("PersistentHomologyOfPureCubicalComplex_Alt");
DeclareGlobalFunction("ZZPersistentHomologyOfPureCubicalComplex");
DeclareGlobalFunction("PersistentHomologyOfQuotientGroupSeries");
DeclareGlobalFunction("PersistentCohomologyOfQuotientGroupSeries");
DeclareGlobalFunction("NormalSeriesToQuotientHomomorphisms");
DeclareGlobalFunction("LinearHomomorphismsPersistenceMat");
DeclareGlobalFunction("LinearHomomorphismsZZPersistenceMat");
DeclareGlobalFunction("PersistentHomologyOfQuotientGroupSeries_Int");
DeclareGlobalFunction("PersistentHomologyOfSubGroupSeries");
DeclareGlobalFunction("PersistentHomologyOfFilteredChainComplex");
DeclareGlobalFunction("TruncatedGComplex");
DeclareGlobalFunction("UniversalBarCode");
DeclareGlobalFunction("UniversalBarCodeEval");
DeclareGlobalFunction("HomologyVectorSpace");
DeclareGlobalFunction("IntegralHomology");
DeclareGlobalFunction("ModularHomology");
DeclareGlobalFunction("GroupHomology");
DeclareGlobalFunction("RipsHomology");
DeclareGlobalFunction("IntegralCohomology");
#DeclareOperation("Cohomology",[IsObject,IsObject]);
DeclareOperation("Cohomology",[IsHapCochain,IsInt]);
DeclareOperation("Cohomology",[IsHapGCocomplex,IsInt]);
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
DeclareGlobalFunction("LefschetzNumberOfChainMap");
DeclareOperation("LefschetzNumber",[IsObject]);

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
DeclareGlobalFunction("PolytopalGenerators");
DeclareGlobalFunction("VectorStabilizer");
DeclareGlobalFunction("PolytopalComplex");
DeclareGlobalFunction("OrbitPolytope");


## POLYCYLIC ########################################################
DeclareGlobalFunction("ResolutionAbelianPcpGroup");
DeclareGlobalFunction("ResolutionNilpotentGroup");

## OBJECTIFICATION ##################################################
#DeclareOperation("Target",[IsObject]);
DeclareAttribute("Target",IsObject);
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

#ReadPackage("HAP","lib/Objectifications/types.gd");


## POLYTOPAL COMPLEXES  ###############################################
DeclareGlobalFunction("ContractibleSubcomplexOfPureCubicalComplex");
DeclareGlobalFunction("ContractibleSubcomplexOfSimplicialComplex");
DeclareGlobalFunction("AcyclicSubcomplexOfPureCubicalComplex");
DeclareGlobalFunction("IntegerSimplicialComplex");
DeclareGlobalFunction("HomotopyEquivalentMaximalPureCubicalSubcomplex");
DeclareGlobalFunction("HomotopyEquivalentMinimalPureCubicalSubcomplex");
#DeclareOperation("EulerCharacteristic",[IsObject]);
#DeclareAttribute("EulerCharacteristic",IsObject);
#ReadPackage("HAP","lib/PolyComplexes/complexTypes.gd");
DeclareAttribute("EulerCharacteristic",IsHapPureCubicalComplex);
DeclareAttribute("EulerCharacteristic",IsHapCubicalComplex);
DeclareAttribute("EulerCharacteristic",IsHapSimplicialComplex);
DeclareOperation("ContractedComplex",[IsObject]);
DeclareGlobalFunction("ReadImageAsPureCubicalComplex");
DeclareGlobalFunction("ReadMatrixAsPureCubicalComplex");
DeclareGlobalFunction("ReadImageSequenceAsPureCubicalComplex");
DeclareGlobalFunction("WritePureCubicalComplexAsImage");
DeclareGlobalFunction("ViewPureCubicalComplex");
DeclareGlobalFunction("PureCubicalComplex");
DeclareGlobalFunction("PureCubicalComplexUnion");
DeclareGlobalFunction("PureCubicalComplexDifference");
DeclareGlobalFunction("PureCubicalComplexIntersection");
DeclareGlobalFunction("PureCubicalComplexToCubicalComplex");
DeclareOperation("ChainComplex",[IsObject]);
DeclareOperation("ChainComplexOfPair",[IsObject,IsObject]);
DeclareGlobalFunction("ChainComplexOfCubicalComplex");
DeclareGlobalFunction("ChainComplexOfCubicalPair");
DeclareGlobalFunction("ExcisedPureCubicalPair");
DeclareGlobalFunction("ExcisedPureCubicalPair_dim_2");
DeclareGlobalFunction("ChainComplexOfSimplicialPair");
DeclareGlobalFunction("ChainMapOfCubicalPairs");
DeclareGlobalFunction("ChainComplexOfSimplicialComplex");
DeclareGlobalFunction("ChainMapOfSimplicialMap");
DeclareGlobalFunction("SkeletonOfSimplicialComplex");
DeclareGlobalFunction("CechComplexOfPureCubicalComplex");
DeclareGlobalFunction("PureComplexToSimplicialComplex");
DeclareGlobalFunction("QuillenComplex");
DeclareGlobalFunction("SimplicialMap");
DeclareGlobalFunction("SimplicialMapNC");
DeclareGlobalFunction("MaximalSimplicesToSimplicialComplex");
DeclareGlobalFunction("SimplicesToSimplicialComplex");
DeclareGlobalFunction("MaximalSimplicesOfSimplicialComplex");
DeclareGlobalFunction("ContractSimplicialComplex");
DeclareGlobalFunction("PathComponentOfPureCubicalComplex");
DeclareOperation("Bettinumbers",[IsObject,IsInt]);
DeclareGlobalFunction("BettinumbersOfPureCubicalComplex_dim_2");
DeclareGlobalFunction("ThickenedPureCubicalComplex");
DeclareGlobalFunction("ThickenedPureCubicalComplex_dim2");
DeclareGlobalFunction("BoundaryOfPureCubicalComplex");
DeclareGlobalFunction("ContractPureCubicalComplex");
DeclareGlobalFunction("ZigZagContractedPureCubicalComplex");
DeclareGlobalFunction("CropPureCubicalComplex");
#DeclareGlobalFunction("SingularChainComplex");
DeclareGlobalFunction("ComplementOfPureCubicalComplex");
DeclareGlobalFunction("SingularitiesOfPureCubicalComplex");
DeclareGlobalFunction("ChainInclusionOfPureCubicalPair");
DeclareGlobalFunction("DirectProductOfPureCubicalComplexes");
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
DeclareGlobalFunction("SymmetricMatrixToGraph");
DeclareGlobalFunction("GraphOfSimplicialComplex");
DeclareGlobalFunction("PathComponentsOfGraph");
DeclareGlobalFunction("PathComponentsOfSimplicialComplex");
DeclareGlobalFunction("PathComponentsOfSimplicialComplex_alt");
DeclareGlobalFunction("ContractGraph");
DeclareGlobalFunction("SimplicialNerveOfGraph");
DeclareGlobalFunction("GraphDisplay");
DeclareGlobalFunction("SkeletonOfCubicalComplex");
DeclareGlobalFunction("MorseFiltration");
DeclareGlobalFunction("ContractCubicalComplex_dim2");
DeclareGlobalFunction("ContractCubicalComplex_dim3");
DeclareGlobalFunction("ContractCubicalComplex");
DeclareGlobalFunction("DVFReducedCubicalComplex");
DeclareGlobalFunction("BoundingPureCubicalComplex");



########################## ARRAYS ################################
DeclareGlobalFunction("ArrayValue");
DeclareGlobalFunction("ArrayValueFunctions");
DeclareGlobalFunction("ArrayValueKD");
DeclareGlobalFunction("ArrayToPureCubicalComplex");
DeclareGlobalFunction("FrameArray");
DeclareGlobalFunction("UnframeArray");
DeclareGlobalFunction("PermuteArray");
DeclareGlobalFunction("ArraySum");
DeclareGlobalFunction("ArrayDimension");
DeclareGlobalFunction("ArrayDimensions");
DeclareGlobalFunction("ContractArray");
DeclareGlobalFunction("HAP_Binlisttoint");
DeclareGlobalFunction("ContractMatrix");
DeclareGlobalFunction("ContractibleSubArray");
DeclareGlobalFunction("HomotopyEquivalentLargerSubArray");
DeclareGlobalFunction("HomotopyEquivalentLargerSubArray3D");
DeclareGlobalFunction("HomotopyEquivalentSmallerSubArray");
DeclareGlobalFunction("HomotopyEquivalentSmallerSubArray3D");
DeclareGlobalFunction("ContractibleSubMatrix");
DeclareGlobalFunction("HomotopyEquivalentLargerSubMatrix");
DeclareGlobalFunction("HomotopyEquivalentSmallerSubMatrix");
DeclareGlobalFunction("ArrayAssign");
DeclareGlobalFunction("UnboundedArrayAssign");
DeclareGlobalFunction("ArrayAssignFunctions");
DeclareGlobalFunction("ArrayIterate");
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
DeclareGlobalFunction("QuotientQuasiIsomorph");
DeclareGlobalFunction("SubQuasiIsomorph");
DeclareGlobalFunction("SylowSubgroupOfCatOneGroup");

## COMMUTATIVE DIAGRAMS #############################################
DeclareGlobalFunction("HomomorphismChainToCommutativeDiagram");
DeclareGlobalFunction("NerveOfCommutativeDiagram");
DeclareGlobalFunction("GroupHomologyOfCommutativeDiagram");
DeclareGlobalFunction("PersistentHomologyOfCommutativeDiagramOfPGroups");
DeclareGlobalFunction("NormalSeriesToQuotientDiagram");

## SIMPLICIAL GROUPS ################################################
DeclareGlobalFunction("NerveOfCatOneGroup");
DeclareGlobalFunction("BarResolutionBoundary");
DeclareGlobalFunction("BarResolutionEquivalence");
DeclareGlobalFunction("BarComplexBoundary");
DeclareGlobalFunction("BarComplexEquivalence");
DeclareGlobalFunction("ReduceAlg");
DeclareGlobalFunction("AddWord");
DeclareGlobalFunction("ChainComplexOfSimplicialGroup");
DeclareGlobalFunction("EilenbergMacLaneSimplicialGroup");
DeclareGlobalFunction("EilenbergMacLaneSimplicialGroup_alt");
DeclareGlobalFunction("TensorProductOfChainComplexes");


## REGULAR CW-SPACES ################################################
DeclareGlobalFunction("SimplicialComplexToRegularCWSpace");
DeclareGlobalFunction("GraphOfRegularCWSpace");#Not yet implemented
DeclareGlobalFunction("CubicalComplexToRegularCWSpace");#Not yet implemented
DeclareGlobalFunction("ContractRegularCWSpace");
DeclareGlobalFunction("RemoveCellFromRegularCWSpace");
DeclareGlobalFunction("CriticalCellsOfRegularCWSpace");
DeclareGlobalFunction("ContractRegularCWSpace_Alt");

## OTHER ############################################################
#ReadPackage("HAP","lib/Objectifications/types.gi");
#ReadPackage("HAP","lib/TopologicalSpaces/topTypes.gd");
#ReadPackage("HAP","lib/PolyComplexes/complexTypes.gd");
#ReadPackage("HAP","lib/GOuterGroups/goutergroup.gd");
ReadPackage("HAP","lib/CategoryTheory/categories.gd");
ReadPackage("HAP","lib/CategoryTheory/commutativeDiagrams.gd");
ReadPackage("HAP","lib/HapPrime/derivation.gd");
ReadPackage("HAP","lib/HapPrime/singular.gd");
ReadPackage("HAP","lib/HapPrime/gradedalgebra.gd");
ReadPackage("HAP","lib/HapPrime/polynomials.gd");
ReadPackage("HAP","lib/HapPrime/ringhomomorphism.gd");
ReadPackage("HAP","lib/HapPrime/rings.gd");
ReadPackage("HAP","lib/HapPrime/happrime.gd");

