#####################################################################
##
##  read.g                  HAP library               Graham Ellis 
##
#####################################################################

HAPconstant:=2;	
SetInfoLevel(InfoWarning,0); #We shouldn't really do this!
ReadPackage("HAP","boolean");
HAP_ROOT:=DirectoriesPackageLibrary("HAP");
HAP_ROOT:=Filename(HAP_ROOT,".");
HAP_ROOT:=HAP_ROOT{[1..Length(HAP_ROOT)-1]};
MakeReadOnlyGlobal("HAP_ROOT");
#POLYMAKE_PATH:=Concatenation(HAP_ROOT,"Polymake/polymakeLegacy ");
#POLYMAKE_PATH:=Filename( DirectoriesSystemPrograms( ), "polymake" );
#MakeReadOnlyGlobal("POLYMAKE_PATH");
#MakeReadWriteGlobal("POLYMAKE_COMMAND");
#POLYMAKE_COMMAND:=POLYMAKE_PATH;
#MakeReadOnlyGlobal(POLYMAKE_COMMAND);


#################################
#################################
ReadPackageHap:=function(file)
local str, i, ln;

if COMPILED=false then
ReadPackage("HAP", file);
else

str:=SplitString(file,'/');
ln:=Length(str);

for i in [2..ln] do
str[i]:=Concatenation("/",str[i]);
od;

str:=Concatenation(str{[1..ln-1]},["/Compiled"],[str[ln]]);
str:=Concatenation(str);
str:=Concatenation(HAP_ROOT{[1..Length(HAP_ROOT)-4]},str);
str:=str{[1..Length(str)-2]};
str:=Concatenation(str,"la.so");
if IsExistingFile(str) then 
LoadDynamicModule(str); 
else
ReadPackage("HAP", file);
fi;
fi;
end;
#################################
#################################

ReadPackageHap( "lib/TitlePage/copyright.gap");

################# NQ COMMANDS #######################################
if  IsPackageMarkedForLoading("nq","1.1") then 
HAP_NqEpimorphismNilpotentQuotient:=NqEpimorphismNilpotentQuotient;
else
if IsPackageMarkedForLoading("nql","1.0") then 
HAP_NqEpimorphismNilpotentQuotient:=NqEpimorphismNilpotentQuotientLpGroup;
else
HAP_NqEpimorphismNilpotentQuotient:=EpimorphismNilpotentQuotient;
fi;
fi;
################# NQ COMMANDS DONE ###############################

################# SIMPHOM COMMANDS ##################################
if not IsPackageMarkedForLoading("homology","0.0") then 
SMInvariantFactors:=function(M); return fail; end;
InfoHomology:=function(M); return fail; end;
else SetInfoLevel(InfoHomology,0);
ReadPackageHap( "lib/Homology/probHomology.gi");
ReadPackageHap( "lib/Homology/sparseprobHomology.gi");
fi;
################ SIMPHOM COMMANDS DONR ##############################

################# EDIM COMMANDS #######################################
if not IsPackageMarkedForLoading("edim","1.2.2") then 
ElementaryDivisorsPPartRk:=function(G); return fail; end;
fi;
################# EDIM COMMANDS DONE ###############################

################# GAPDOC COMMANDS #######################################
#if not IsPackageMarkedForLoading("gapdoc","0.0") then 
#MakeGAPDocDoc:=function(G); return fail; end;
#fi;
################# GAPDOC COMMANDS DONE ###############################

################# CONGRUENE COMMANDS #######################################
if not IsPackageMarkedForLoading("congruence","0.0") then
CongruenceSubgroupGamma0:=function(m); return fail; end;
fi;
################# CONGRUENCE COMMANDS DONE ###############################

################# POLYMAKING COMMANDS ####################################
if IsPackageMarkedForLoading("polymaking","0.0") then
ReadPackageHap( "lib/Polymake/convexCWspace.gi");
ReadPackageHap( "lib/Polymake/fix.gi");
fi;
################# POLYMAKING COMMANDS DONE #################################

################# HAPCRYST COMMANDS ####################################
if IsPackageMarkedForLoading("HAPcryst","0.0") then
ReadPackageHap( "lib/RegularCWComplexes/equivariantCW.gi");
fi;
################# POLYMAKING COMMANDS DONE #################################

################# XMOD COMMANDS ############################################
if IsPackageMarkedForLoading("xmod","0.0") then
ReadPackageHap( "lib/SimplicialGroups/identity.gi");
fi;
################# XMOD COMMANDS DONE #######################################

ReadPackageHap( "lib/TitlePage/makeHapMan.gi");

################# OBJECTIFICATIONS ###############################
ReadPackageHap( "lib/Objectifications/basicMethods.gi");
################# OBJECTIFICATIONS DONE ##########################



################# ACLIB & NQ COMMANDS #####################################
if IsPackageMarkedForLoading("aclib","1.1") then
ReadPackageHap( "lib/Resolutions/resACgroup.gi");
ReadPackageHap( "lib/Resolutions/resACquotient.gi");
else
IsAlmostCrystallographic:=function(G); return fail; end;
fi;
if IsPackageMarkedForLoading("nq","1.1") then
ReadPackageHap( "lib/NonabelianTensor/epiNilGrp.gi");
ReadPackageHap( "lib/NonabelianTensor/multNilGrp.gi");
ReadPackageHap( "lib/NonabelianTensor/tensorSquareInf.gi");
ReadPackageHap( "lib/NonabelianTensor/nil2TensorSquare.gi");
ReadPackageHap( "lib/NonabelianTensor/symmetricSquareInf.gi");
ReadPackageHap( "lib/NonabelianTensor/tensor.gi");
fi;

################# ACLIB DONE #########################################

##################### RENAME GAP FUNCTIONS ##########################
AbsInt_HAP:=AbsInt;
MakeReadOnlyGlobal("AbsInt_HAP");
SignInt_HAP:=SignInt;
MakeReadOnlyGlobal("SignInt_HAP");

##################### FREE G MODULES ################################
ReadPackageHap( "lib/FreeGmodules/wordOperations.gi");
ReadPackageHap( "lib/FreeGmodules/tietze.gi");

##################### FPG MODULES ##################################
ReadPackageHap( "lib/FpGmodules/fpgbasics.gi");
ReadPackageHap( "lib/FpGmodules/resfpgmod.gi");
ReadPackageHap( "lib/FpGmodules/homs.gi");

##################### NONABELIAN TENSOR #############################
ReadPackageHap( "lib/NonabelianTensor/tensorSquare.gi");
ReadPackageHap( "lib/NonabelianTensor/tensorPair.gi");
ReadPackageHap( "lib/NonabelianTensor/tensorPairInf.gi");
ReadPackageHap( "lib/NonabelianTensor/tensorPair.alt");
ReadPackageHap( "lib/NonabelianTensor/exteriorProduct.gi");
ReadPackageHap( "lib/NonabelianTensor/SBG.gi");
ReadPackageHap( "lib/NonabelianTensor/symmetricSquare.gi");
#ReadPackageHap( "lib/NonabelianTensor/symmetricSquareInf.gi");
ReadPackageHap( "lib/NonabelianTensor/bogomolov.gi");
ReadPackageHap( "lib/NonabelianTensor/equivalenceclasses.gi");
ReadPackageHap( "lib/NonabelianTensor/weak.gi");


##################### RESOLUTIONS ###################################
ReadPackageHap( "lib/Resolutions/resAspherical.gi");
ReadPackageHap( "lib/Resolutions/resAbGroup.gi");
ReadPackageHap( "lib/Resolutions/resFiniteGroup.gi");
ReadPackageHap( "lib/Resolutions/barComplexMonoid.gi");
ReadPackageHap( "lib/Resolutions/resSmallFpGroup.gi");
ReadPackageHap( "lib/Resolutions/presentation.gi");
ReadPackageHap( "lib/Resolutions/resSubgroup.gi");
ReadPackageHap( "lib/Resolutions/resInfSubgroup.gi");
ReadPackageHap( "lib/Resolutions/resGeneric.gi");
ReadPackageHap( "lib/Resolutions/tietzered.gi");
ReadPackageHap( "lib/Resolutions/coreducedRes.gi");
ReadPackageHap( "lib/Resolutions/pseudoLists.gi");
ReadPackageHap( "lib/Resolutions/resSL2Z.gi");
ReadPackageHap( "lib/Resolutions/resBianchi.gi");
ReadPackageHap( "lib/Resolutions/gens.gi");

##################### RESOLUTIONS MOD P #############################
ReadPackageHap( "lib/ResolutionsModP/resPrimeGroup.gi");
ReadPackageHap( "lib/ResolutionsModP/resPrimeGroupSparse.gi");
ReadPackageHap( "lib/ResolutionsModP/ranksPrimeGroup.gi");
ReadPackageHap( "lib/ResolutionsModP/poincare.gi");
#ReadPackageHap( "lib/ResolutionsModP/primepart.gi");
ReadPackageHap( "lib/ResolutionsModP/radical.gi");

##################### FUNCTORS ######################################
ReadPackageHap( "lib/Functors/permMatrix.gi");
ReadPackageHap( "lib/Functors/homToZmodule.gi");
ReadPackageHap( "lib/Functors/tensorWithZ.gi");
ReadPackageHap( "lib/Functors/tensorWithZmodule.gi");
ReadPackageHap( "lib/Functors/tensorWithTwistedZ.gi");
ReadPackageHap( "lib/Functors/tensorWithTwistedZmodP.gi");
ReadPackageHap( "lib/Functors/tensorWithZmodP.gi");
ReadPackageHap( "lib/Functors/various.gi");
ReadPackageHap( "lib/Functors/equiChainMap.gi");
ReadPackageHap( "lib/Functors/modularEquiChainMap.gi");
ReadPackageHap( "lib/Functors/primePartDerived.gi");
ReadPackageHap( "lib/Functors/primePartDerivedvsgc.gi");
ReadPackageHap( "lib/Functors/homToZ.gi");
ReadPackageHap( "lib/Functors/tensorWithRationals.gi");
ReadPackageHap( "lib/Functors/homToZmodP.gi");
ReadPackageHap( "lib/Functors/homtint.gi");
ReadPackageHap( "lib/Functors/transfer.gi");
ReadPackageHap( "lib/Functors/bianchiHomogeneousPolys.gi");
ReadPackageHap( "lib/Functors/homogeneousPolys.gi");
ReadPackageHap( "lib/Functors/simplify.gi");

##################### CONGRUENCE ######################################
ReadPackageHap( "lib/Congruence/cong.gi");
ReadPackageHap( "lib/Congruence/hecke.gi");
ReadPackageHap( "lib/Congruence/sl2subgroups.gi");
ReadPackageHap( "lib/Congruence/sl2word.gi");
ReadPackageHap( "lib/Congruence/cuspidal.gi");
ReadPackageHap( "lib/Congruence/quadraticIntegers.gi");
ReadPackageHap( "lib/Congruence/sl2o.gi");
ReadPackageHap( "lib/Congruence/volume.gi");





##################### HOMOLOGY ######################################
ReadPackageHap( "lib/Homology/integralHomology.gi");
#ReadPackageHap( "lib/Homology/integralHomology.may2017working");
ReadPackageHap( "lib/Homology/lefschetz.gi");
ReadPackageHap( "lib/Homology/modularHomology.gi");
ReadPackageHap( "lib/Homology/modularHomologyVectSpace.gi");
ReadPackageHap( "lib/Homology/homology.gi");
ReadPackageHap( "lib/Homology/cat1homology.gi");
ReadPackageHap( "lib/Homology/groupHomology.gi");
ReadPackageHap( "lib/Homology/integralCohomology.gi");
ReadPackageHap( "lib/Homology/cohomology.gi");
ReadPackageHap( "lib/Homology/syzygy.gi");
ReadPackageHap( "lib/Homology/cycles.gi");
ReadPackageHap( "lib/Homology/cocycleCondition.gi");
ReadPackageHap( "lib/Homology/isSuperperfect.gi");
ReadPackageHap( "lib/Homology/modularCohomology.gi");
ReadPackageHap( "lib/Homology/solutionsMat.gi");
ReadPackageHap( "lib/Homology/groupCohomology.gi");
ReadPackageHap( "lib/Homology/integralHomologyObj.gi");
ReadPackageHap( "lib/Homology/integralCohomologyObj.gi");
ReadPackageHap( "lib/Homology/persistent.gi");

##################### PERTURBATIONS #################################
ReadPackageHap( "lib/Perturbations/resExtension.gi");
ReadPackageHap( "lib/Perturbations/resDirectProd.gi");
ReadPackageHap( "lib/Perturbations/twistedTensorProduct.gi");
ReadPackageHap( "lib/Perturbations/resFiniteExt.gi");
ReadPackageHap( "lib/Perturbations/resNormalSer.gi");
ReadPackageHap( "lib/Perturbations/resFiniteDirectProd.gi");
ReadPackageHap( "lib/Perturbations/resSubNormSeries.gi");
ReadPackageHap( "lib/Perturbations/freeRes.gi");
ReadPackageHap( "lib/Perturbations/dutour.gi");
ReadPackageHap( "lib/Perturbations/contractibleSL2Zcomplex.gi");
ReadPackageHap( "lib/Perturbations/contractibleSL2ZcomplexALT.gi");
ReadPackageHap( "lib/Perturbations/filteredChainComplex.gi");
ReadPackageHap( "lib/Perturbations/resDirectProdLazy.gi"); 
ReadPackageHap( "lib/Perturbations/conrad.gi");



#################### ARTIN COXETER ##################################
ReadPackageHap( "lib/ArtinCoxeter/diagrams.gi");
ReadPackageHap( "lib/ArtinCoxeter/resArtin.gi");
ReadPackageHap( "lib/ArtinCoxeter/coxeterWythoff.gi");
ReadPackageHap( "lib/ArtinCoxeter/noncrossing.gi");


#################### COHOMOLOGY RINGS ###############################
ReadPackageHap( "lib/Rings/intCoh.gi");
ReadPackageHap( "lib/Rings/cocycleChainMap.gi");
ReadPackageHap( "lib/Rings/cupProduct.gi");
ReadPackageHap( "lib/Rings/integralGens.gi");

################### POLYMAKE #######################################
ReadPackageHap( "lib/Polymake/aspherical.gi");
ReadPackageHap( "lib/Polymake/polyGens.gi");
ReadPackageHap( "lib/Polymake/stabilizer.gi");
ReadPackageHap( "lib/Polymake/polyFaces.gi");
ReadPackageHap( "lib/Polymake/orbitPoly.gi");
#ReadPackageHap( "lib/Polymake/fix.gi");
ReadPackageHap( "lib/Polymake/TZ.gi");

################### POLYCYLIC ######################################
ReadPackageHap( "lib/Polycyclic/resAbPcpGroup.gi");
ReadPackageHap( "lib/Polycyclic/resNilpotentPcpGrp.gi");



################### MOD P RINGS ####################################
ReadPackageHap( "lib/ModPRings/record.gi");
ReadPackageHap( "lib/ModPRings/recordPart1.gi");
ReadPackageHap( "lib/ModPRings/recordPartII.gi");


################### GRAPHS OF GROUPS ###############################
ReadPackageHap( "lib/GraphsOfGroups/graphs.gi");
ReadPackageHap( "lib/GraphsOfGroups/resGraph.gi");
ReadPackageHap( "lib/GraphsOfGroups/graphOfResolutions.gi");


################### TEST ###########################################
ReadPackageHap( "tst/test.g");

################### STREAMS ########################################
ReadPackageHap("lib/Streams/streams.gi");
ReadPackageHap("lib/Streams/HAPexport.gi");
ReadPackageHap("lib/Streams/HAPimport.gi");

################### RESOLUTIONS (CONTD) ############################
ReadPackageHap("lib/Resolutions/cayley.gi");

################### Lie Algebras ###################################
ReadPackageHap("lib/LieAlgebras/chevalleyEilenberg.gi");
ReadPackageHap("lib/LieAlgebras/isLieHom.gi");
ReadPackageHap("lib/LieAlgebras/groupToLie.gi");
ReadPackageHap("lib/LieAlgebras/leibniz.gi");
ReadPackageHap("lib/LieAlgebras/LieTensorSquare.gi");
ReadPackageHap("lib/LieAlgebras/LieCover.gi");
ReadPackageHap("lib/LieAlgebras/LeibnizQuasiCover.gi");
ReadPackageHap("lib/LieAlgebras/LieExteriorSquare.gi");

################### MEAT AXE #######################################
ReadPackageHap("lib/FpGmodules/meataxe.gi");

################## PURE COMPLEXES #############################
ReadPackageHap("lib/PureComplexes/pureComplexes.gi");


################## POLYTOPAL COMPLEXES #############################
ReadPackageHap("lib/PolyComplexes/arrayOps.gi");
ReadPackageHap("lib/PolyComplexes/pureCubicalComplexes.gi");
ReadPackageHap("lib/PolyComplexes/sparseCubicalComplexes.gi");
ReadPackageHap("lib/PolyComplexes/chainComplexes.gi");
ReadPackageHap("lib/PolyComplexes/twoDimensional.gi");
ReadPackageHap("lib/PolyComplexes/twoDimensionalPerm.gi");
ReadPackageHap("lib/PolyComplexes/threeDimensional.gi");
ReadPackageHap("lib/PolyComplexes/threeDimensionalPerm.gi");
ReadPackageHap("lib/PolyComplexes/dvf.gi");
ReadPackageHap("lib/PolyComplexes/rips.gi");

ReadPackageHap("lib/PolyComplexes/simplicialComplexes.gi");
ReadPackageHap("lib/PolyComplexes/groupComplexes.gi");
ReadPackageHap("lib/PolyComplexes/cluster.gi");
ReadPackageHap("lib/PolyComplexes/clique.gi");
ReadPackageHap("lib/PolyComplexes/metrics.gi");
ReadPackageHap("lib/PolyComplexes/graphviz.gi");
ReadPackageHap("lib/PolyComplexes/hap2chomp.gi");
ReadPackageHap("lib/PolyComplexes/filteredCubical.gi");
ReadPackageHap("lib/PolyComplexes/purePermutahedralComplexes.gi");


################## CATEGORY THEORY #################################
ReadPackageHap("lib/CategoryTheory/categories.gi");
ReadPackageHap("lib/CategoryTheory/commutativeDiagrams.gi");

################## CAT ONE GROUPS ##################################
ReadPackageHap("lib/CatGroups/CatConstructions.gi");
ReadPackageHap("lib/CatGroups/CatBasic.gi");
ReadPackageHap("lib/CatGroups/identities.gi");
ReadPackageHap("lib/CatGroups/algIdentities.gi");
ReadPackageHap("lib/CatGroups/CrossedInvariant.gi");
ReadPackageHap("lib/CatGroups/eilenberg.gi");
ReadPackageHap("lib/CatGroups/coeilenberg.gi");

################## G-OUTER GROUPS ##################################
ReadPackageHap("lib/GOuterGroups/goutergroup.gi");
ReadPackageHap("lib/GOuterGroups/homtogouter.gi");
ReadPackageHap("lib/GOuterGroups/functorialGouter.gi");
ReadPackageHap("lib/GOuterGroups/hadamard.gi");

################## SIMPLICIAL GROUPS ###############################
#ReadPackageHap("lib/SimplicialGroups/eilenbergMacLane.gi");
#ReadPackageHap("lib/SimplicialGroups/nerveOfCatOneGroup.gi");
#ReadPackageHap("lib/SimplicialGroups/mooreComplex.gi");
#ReadPackageHap("lib/SimplicialGroups/barResolutionEquivalence.gi");
#ReadPackageHap("lib/SimplicialGroups/barComplexEquivalence.gi");
#ReadPackageHap("lib/SimplicialGroups/chainComplexOfSimplicialGroup.gi");
#ReadPackageHap("lib/SimplicialGroups/tensor2chains.gi");
#ReadPackageHap("lib/SimplicialGroups/homotopyLowerCenterSeries.gi");
#ReadPackageHap("lib/SimplicialGroups/crossedModule.gi");
#ReadPackageHap("lib/SimplicialGroups/quasiIsomorph.gi");
#ReadPackageHap("lib/SimplicialGroups/homology.gi");
#ReadPackageHap("lib/SimplicialGroups/persistentHomology.gi");
#ReadPackageHap("lib/SimplicialGroups/dataCatOneGroups.gi");
#ReadPackageHap("lib/SimplicialGroups/catOneGroupsByGroup.gi");
#ReadPackageHap("lib/SimplicialGroups/quasiCatOneGroup.gi");

ReadPackageHap("lib/SimplicialGroups/dataCatOneGroups.data");
ReadPackageHap("lib/SimplicialGroups/dataTwoTypes.data");
ReadPackageHap("lib/SimplicialGroups/eilenbergMacLane.gi");
ReadPackageHap("lib/SimplicialGroups/nerveOfCatOneGroup.gi");
ReadPackageHap("lib/SimplicialGroups/mooreComplex.gi");
ReadPackageHap("lib/SimplicialGroups/barResolutionEquivalence.gi");
ReadPackageHap("lib/SimplicialGroups/barComplexEquivalence.gi");
ReadPackageHap("lib/SimplicialGroups/chainComplexOfSimplicialGroup.gi");
ReadPackageHap("lib/SimplicialGroups/tensor2chains.gi");
ReadPackageHap("lib/SimplicialGroups/homotopyLowerCenterSeries.gi");
ReadPackageHap("lib/SimplicialGroups/quasiIsomorph.gi");
ReadPackageHap("lib/SimplicialGroups/homology.gi");
ReadPackageHap("lib/SimplicialGroups/persistentHomology.gi");
ReadPackageHap("lib/SimplicialGroups/crossedModule.gi");
ReadPackageHap("lib/SimplicialGroups/catOneGroup.gi");
ReadPackageHap("lib/SimplicialGroups/twoTypes.gi");




################## REGULAR CW_COMPLEXES ###############################
ReadPackageHap("lib/RegularCWComplexes/basicRegular.gi");
#ReadPackageHap("lib/RegularCWComplexes/contractAlt.gi");
ReadPackageHap("lib/RegularCWComplexes/fundamental.gi");
ReadPackageHap("lib/RegularCWComplexes/cocontract.gi");
ReadPackageHap("lib/RegularCWComplexes/piZero.gi");
ReadPackageHap("lib/RegularCWComplexes/filteredCW.gi");
ReadPackageHap("lib/RegularCWComplexes/directproduct.gi");
ReadPackageHap("lib/RegularCWComplexes/cocycle.gi");
ReadPackageHap("lib/RegularCWComplexes/universalCover.gi");
ReadPackageHap("lib/RegularCWComplexes/spin.gi");
ReadPackageHap("lib/RegularCWComplexes/spunknotcomp.gi");
ReadPackageHap("lib/RegularCWComplexes/grannyknot.gi");

if IsPackageMarkedForLoading("congruence","0.0") then
############### ARITHMETIC GROUPS#########################
 ReadPackageHap("lib/ArithmeticGroups/arithVarious.gi");
 ReadPackageHap("lib/ArithmeticGroups/arithmeticOps.gi");
 ReadPackageHap("lib/ArithmeticGroups/sl2zngens.gi");
 ReadPackageHap("lib/ArithmeticGroups/cplGTree.gi");
 ReadPackageHap("lib/ArithmeticGroups/resGTree.gi");
 ReadPackageHap("lib/ArithmeticGroups/sl2zres.gi");
ReadPackageHap("lib/ArithmeticGroups/sl2zresalt.gi");
# ReadPackageHap("lib/ArithmeticGroups/resDirectProd.gi");
ReadPackageHap("lib/ArithmeticGroups/barycentric.gi");
######################################################
fi;

################## KNOTS ##########################################
ReadPackageHap("/lib/Knots/knotdata.gi");
ReadPackageHap("/lib/Knots/cubicalKnot.gi");
ReadPackageHap("/lib/Knots/csvknot.gi");
ReadPackageHap("/lib/Knots/surface.gi");

################## SPARSE ##########################################
ReadPackageHap("lib/Sparse/sparse.gi");


################# CRYSTALLOGRAPHIC GROUPS #####################
if IsPackageMarkedForLoading("HAPcryst","0.0") then
ReadPackageHap("lib/ArithmeticGroups/crystGbasis.gi");
#ReadPackageHap("lib/ArithmeticGroups/crystGbasisnew.gi");
ReadPackageHap("lib/ArithmeticGroups/crystVarious.gi");
ReadPackageHap("lib/ArithmeticGroups/crystGcomplex.gi");
#ReadPackageHap("lib/ArithmeticGroups/crystGcomplexnew.gi");
ReadPackageHap("lib/ArithmeticGroups/freeZGRes.gi");
fi;

ReadPackageHap("lib/Operations/hapOps.gi");

################## HAP PRIME ##################################
if IsPackageMarkedForLoading("singular","06.07.23") then
ReadPackageHap("lib/HapPrime/singular.gi");
ReadPackageHap("lib/HapPrime/rings.gi");
ReadPackageHap("lib/HapPrime/ringhomomorphism.gi");
ReadPackageHap("lib/HapPrime/gradedalgebra.gi");
ReadPackageHap("lib/HapPrime/polynomials.gi");
ReadPackageHap("lib/HapPrime/happrime.gi");
ReadPackageHap("lib/HapPrime/derivation.gi");
#ReadPackageHap("lib/HapPrime/poincare.gd");
#ReadPackageHap("lib/HapPrime/poincare.gi");
fi;


ReadPackageHap("lib/Homology/BarCodes/barcode.gi");
ReadPackageHap("lib/TorsionSubcomplexes/torsioninit.gi");
ReadPackageHap("lib/RahmSanchez/DavisComplex.gi");

################## Cohomology Operations #####################
ReadPackageHap("lib/CohomologyOperations/cohomology_homomorphism.gi");
ReadPackageHap("lib/CohomologyOperations/connecting_homomorphism.gi");
ReadPackageHap("lib/CohomologyOperations/steenrod.gi");
ReadPackageHap("lib/CohomologyOperations/mycupi.gi");
ReadPackageHap("lib/CohomologyOperations/toplevelsquares.gi");
ReadPackageHap("lib/CohomologyOperations/detection.gi");
ReadPackageHap("lib/CohomologyOperations/cohodata.gi");
ReadPackageHap("lib/CohomologyOperations/stiefel.gi");

################## HAP COCYCLIC ##############################
if not IsPackageMarkedForLoading("HAPcocyclic","0") then
ReadPackageHap("lib/HapCocyclic/init.g");
ReadPackageHap("lib/HapCocyclic/read.g");
fi;

################## QUANDLES ##################################
ReadPackageHap("lib/Quandles/quandles.gi");
ReadPackageHap("lib/Quandles/planarDiagramData.gi");
ReadPackageHap("lib/Quandles/quandleKnots.gi");
ReadPackageHap("lib/Quandles/quandleOrbits.gi");
ReadPackageHap("lib/Quandles/isoreps.gi");

################## MANIFOLDS #################################
ReadPackageHap("lib/Manifolds/manifolds.gi");
ReadPackageHap("lib/Manifolds/simplicialcup.gi");

################## KELVIN ##################################
ReadPackageHap("lib/Kelvin/init.g");



SetInfoLevel(InfoWarning,1); #This is GAP's default level

##AFTER THOUGHT######################################################
ReadPackageHap("lib/hap_afterthought.gd");



