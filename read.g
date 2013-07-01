#####################################################################
##
##  read.g                  HAP library               Graham Ellis 
##
#####################################################################

HAPconstant:=2;	
SetInfoLevel(InfoWarning,0); #We shouldn't really do this!

#From version 4.5 of GAP we'll use the new function IsPackageMarkedForLoading
if not CompareVersionNumbers(VERSION,"4.5") then
        IsPackageMarkedForLoading:=function(ver,num) local b;
        b:=LoadPackage(ver,num,false);
        if b=true then return b; else return false; fi;
        end;
fi;

ReadPackage("HAP","boolean");
ReadPackage("HAP", "lib/TitlePage/copyright.gap");
HAP_ROOT:=DirectoriesPackageLibrary("HAP");
HAP_ROOT:=Filename(HAP_ROOT,".");
HAP_ROOT:=HAP_ROOT{[1..Length(HAP_ROOT)-1]};
MakeReadOnlyGlobal("HAP_ROOT");

################### POLYCYLIC COMMANDS ##############################
## Most HAP functions should work on pcp groups if the polycyclic package 
## is installed. Otherwise we need to give a meaning to certain commands
## defined in the polycyclic package.
#if not IsPackageMarkedForLoading("polycyclic","1.1") then 
#DeclareOperation("NaturalHomomorphism",[IsGroup,IsGroup]);
#IsPcpGroup:=function(G);return false;end;
#Collector:=function(x);return fail; end;
#PcpGroupByCollector:=function(x);return fail; end;
#Igs:=function(x);return fail; end;
#GenExpList:=function(x);return fail; end;
#HeisenbergPcpGroup:=function(x);return fail; end;
#Pcp:=function(G,D); return fail; end;
#IsAlmostCrystallographic:=function(G); return fail; end;
#GeneratorsOfPcp:=function(G); return fail; end;
#IsomorphismPcpGroup:=function(G);return fail;end;
#AbelianPcpGroup:=function(G);return fail;end;
#fi;
################## POLYCYCLIC COMMANDS DONE #########################

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
ReadPackage("HAP", "lib/Homology/probHomology.gi");
ReadPackage("HAP", "lib/Homology/sparseprobHomology.gi");
fi;
################ SIMPHOM COMMANDS DONR ##############################

################# EDIM COMMANDS #######################################
if not IsPackageMarkedForLoading("edim","1.2.2") then 
ElementaryDivisorsPPartRk:=function(G); return fail; end;
fi;
################# EDIM COMMANDS DONE ###############################

################# GAPDOC COMMANDS #######################################
if not IsPackageMarkedForLoading("gapdoc","0.0") then 
MakeGAPDocDoc:=function(G); return fail; end;
fi;
################# GAPDOC COMMANDS DONE ###############################

################# CONGRUENE COMMANDS #######################################
if not IsPackageMarkedForLoading("congruence","0.0") then
CongruenceSubgroupGamma0:=function(m); return fail; end;
fi;
################# CONGRUENCE COMMANDS DONE ###############################

################# POLYMAKING COMMANDS ####################################
if IsPackageMarkedForLoading("polymaking","0.0") then
ReadPackage("HAP", "lib/Polymake/convexCWspace.gi");
fi;
################# POLYMAKING COMMANDS DONE #################################

################# HAPCRYST COMMANDS ####################################
if IsPackageMarkedForLoading("HAPcryst","0.0") then
ReadPackage("HAP", "lib/RegularCWComplexes/equivariantCW.gi");
fi;
################# POLYMAKING COMMANDS DONE #################################



ReadPackage("HAP", "lib/TitlePage/makeHapMan.gi");

################# OBJECTIFICATIONS ###############################
ReadPackage("HAP", "lib/Objectifications/basicMethods.gi");
################# OBJECTIFICATIONS DONE ##########################

if COMPILED=true then
ReadPackage("HAP","lib/compiledVersion.gap");
fi;

if COMPILED=false then

################# ACLIB & NQ COMMANDS #####################################
if IsPackageMarkedForLoading("aclib","1.1") then
ReadPackage("HAP", "lib/Resolutions/resACgroup.gi");
ReadPackage("HAP", "lib/Resolutions/resACquotient.gi");
else
IsAlmostCrystallographic:=function(G); return fail; end;
fi;
if IsPackageMarkedForLoading("nq","1.1") then
ReadPackage("HAP", "lib/NonabelianTensor/epiNilGrp.gi");
ReadPackage("HAP", "lib/NonabelianTensor/multNilGrp.gi");
ReadPackage("HAP", "lib/NonabelianTensor/tensorSquareInf.gi");
ReadPackage("HAP", "lib/NonabelianTensor/symmetricSquareInf.gi");
fi;

################# ACLIB DONE #########################################

##################### RENAME GAP FUNCTIONS ##########################
AbsInt_HAP:=AbsInt;
MakeReadOnlyGlobal("AbsInt_HAP");
SignInt_HAP:=SignInt;
MakeReadOnlyGlobal("SignInt_HAP");

##################### FREE G MODULES ################################
ReadPackage("HAP", "lib/FreeGmodules/wordOperations.gi");
ReadPackage("HAP", "lib/FreeGmodules/tietze.gi");

##################### FPG MODULES ##################################
ReadPackage("HAP", "lib/FpGmodules/fpgbasics.gi");
ReadPackage("HAP", "lib/FpGmodules/resfpgmod.gi");
ReadPackage("HAP", "lib/FpGmodules/homs.gi");

##################### NONABELIAN TENSOR #############################
ReadPackage("HAP", "lib/NonabelianTensor/tensorSquare.gi");
ReadPackage("HAP", "lib/NonabelianTensor/tensorPair.gi");
ReadPackage("HAP", "lib/NonabelianTensor/exteriorProduct.gi");
ReadPackage("HAP", "lib/NonabelianTensor/SBG.gi");
ReadPackage("HAP", "lib/NonabelianTensor/symmetricSquare.gi");
#ReadPackage("HAP", "lib/NonabelianTensor/symmetricSquareInf.gi");
ReadPackage("HAP", "lib/NonabelianTensor/bogomolov.gi");

##################### RESOLUTIONS ###################################
ReadPackage("HAP", "lib/Resolutions/resAspherical.gi");
ReadPackage("HAP", "lib/Resolutions/resAbGroup.gi");
ReadPackage("HAP", "lib/Resolutions/resFiniteGroup.gi");
ReadPackage("HAP", "lib/Resolutions/resSmallFpGroup.gi");
ReadPackage("HAP", "lib/Resolutions/presentation.gi");
ReadPackage("HAP", "lib/Resolutions/resSubgroup.gi");
ReadPackage("HAP", "lib/Resolutions/resInfSubgroup.gi");
ReadPackage("HAP", "lib/Resolutions/resGeneric.gi");
ReadPackage("HAP", "lib/Resolutions/tietzered.gi");
ReadPackage("HAP", "lib/Resolutions/coreducedRes.gi");
ReadPackage("HAP", "lib/Resolutions/pseudoLists.gi");
ReadPackage("HAP", "lib/Resolutions/resSL2Z.gi");

##################### RESOLUTIONS MOD P #############################
ReadPackage("HAP", "lib/ResolutionsModP/resPrimeGroup.gi");
ReadPackage("HAP", "lib/ResolutionsModP/ranksPrimeGroup.gi");
ReadPackage("HAP", "lib/ResolutionsModP/poincare.gi");
#ReadPackage("HAP", "lib/ResolutionsModP/primepart.gi");

##################### FUNCTORS ######################################
ReadPackage("HAP", "lib/Functors/permMatrix.gi");
ReadPackage("HAP", "lib/Functors/homToZmodule.gi");
ReadPackage("HAP", "lib/Functors/tensorWithZ.gi");
ReadPackage("HAP", "lib/Functors/tensorWithZmodule.gi");
ReadPackage("HAP", "lib/Functors/tensorWithTwistedZ.gi");
ReadPackage("HAP", "lib/Functors/tensorWithTwistedZmodP.gi");
ReadPackage("HAP", "lib/Functors/tensorWithZmodP.gi");
ReadPackage("HAP", "lib/Functors/various.gi");
ReadPackage("HAP", "lib/Functors/equiChainMap.gi");
ReadPackage("HAP", "lib/Functors/modularEquiChainMap.gi");
ReadPackage("HAP", "lib/Functors/primePartDerived.gi");
ReadPackage("HAP", "lib/Functors/homToZ.gi");
ReadPackage("HAP", "lib/Functors/tensorWithRationals.gi");
ReadPackage("HAP", "lib/Functors/homToZmodP.gi");


##################### HOMOLOGY ######################################
ReadPackage("HAP", "lib/Homology/integralHomology.gi");
ReadPackage("HAP", "lib/Homology/lefschetz.gi");
ReadPackage("HAP", "lib/Homology/modularHomology.gi");
ReadPackage("HAP", "lib/Homology/modularHomologyVectSpace.gi");
ReadPackage("HAP", "lib/Homology/homology.gi");
ReadPackage("HAP", "lib/Homology/cat1homology.gi");
ReadPackage("HAP", "lib/Homology/groupHomology.gi");
ReadPackage("HAP", "lib/Homology/integralCohomology.gi");
ReadPackage("HAP", "lib/Homology/cohomology.gi");
ReadPackage("HAP", "lib/Homology/syzygy.gi");
ReadPackage("HAP", "lib/Homology/cycles.gi");
ReadPackage("HAP", "lib/Homology/cocycleCondition.gi");
ReadPackage("HAP", "lib/Homology/isSuperperfect.gi");
ReadPackage("HAP", "lib/Homology/modularCohomology.gi");
ReadPackage("HAP", "lib/Homology/solutionsMat.gi");
ReadPackage("HAP", "lib/Homology/groupCohomology.gi");
ReadPackage("HAP", "lib/Homology/integralHomologyObj.gi");
ReadPackage("HAP", "lib/Homology/integralCohomologyObj.gi");
ReadPackage("HAP", "lib/Homology/persistent.gi");

##################### PERTURBATIONS #################################
ReadPackage("HAP", "lib/Perturbations/resExtension.gi");
ReadPackage("HAP", "lib/Perturbations/resDirectProd.gi");
ReadPackage("HAP", "lib/Perturbations/twistedTensorProduct.gi");
ReadPackage("HAP", "lib/Perturbations/resFiniteExt.gi");
ReadPackage("HAP", "lib/Perturbations/resNormalSer.gi");
ReadPackage("HAP", "lib/Perturbations/resFiniteDirectProd.gi");
ReadPackage("HAP", "lib/Perturbations/resSubNormSeries.gi");
ReadPackage("HAP", "lib/Perturbations/freeRes.gi");
ReadPackage("HAP", "lib/Perturbations/dutour.gi");
ReadPackage("HAP", "lib/Perturbations/contractibleSL2Zcomplex.gi");
ReadPackage("HAP", "lib/Perturbations/filteredChainComplex.gi");



#################### ARTIN COXETER ##################################
ReadPackage("HAP", "lib/ArtinCoxeter/diagrams.gi");
ReadPackage("HAP", "lib/ArtinCoxeter/resArtin.gi");
ReadPackage("HAP", "lib/ArtinCoxeter/coxeterWythoff.gi");
ReadPackage("HAP", "lib/ArtinCoxeter/noncrossing.gi");


#################### COHOMOLOGY RINGS ###############################
ReadPackage("HAP", "lib/Rings/intCoh.gi");
ReadPackage("HAP", "lib/Rings/cocycleChainMap.gi");
ReadPackage("HAP", "lib/Rings/cupProduct.gi");
ReadPackage("HAP", "lib/Rings/integralGens.gi");

################### POLYMAKE #######################################
ReadPackage("HAP", "lib/Polymake/aspherical.gi");
ReadPackage("HAP", "lib/Polymake/polyGens.gi");
ReadPackage("HAP", "lib/Polymake/stabilizer.gi");
ReadPackage("HAP", "lib/Polymake/polyFaces.gi");
ReadPackage("HAP", "lib/Polymake/orbitPoly.gi");

################### POLYCYLIC ######################################
ReadPackage("HAP", "lib/Polycyclic/resAbPcpGroup.gi");
ReadPackage("HAP", "lib/Polycyclic/resNilpotentPcpGrp.gi");



################### MOD P RINGS ####################################
ReadPackage("HAP", "lib/ModPRings/record.gi");
ReadPackage("HAP", "lib/ModPRings/recordPart1.gi");
ReadPackage("HAP", "lib/ModPRings/recordPartII.gi");


################### GRAPHS OF GROUPS ###############################
ReadPackage("HAP", "lib/GraphsOfGroups/graphs.gi");
ReadPackage("HAP", "lib/GraphsOfGroups/resGraph.gi");
ReadPackage("HAP", "lib/GraphsOfGroups/graphOfResolutions.gi");

fi;

################### TEST ###########################################
ReadPackage("HAP", "test/test.gap");

################### STREAMS ########################################
ReadPackage("HAP","lib/Streams/streams.gi");
ReadPackage("HAP","lib/Streams/HAPexport.gi");
ReadPackage("HAP","lib/Streams/HAPimport.gi");

################### RESOLUTIONS (CONTD) ############################
ReadPackage("HAP","lib/Resolutions/cayley.gi");

################### Lie Algebras ###################################
ReadPackage("HAP","lib/LieAlgebras/chevalleyEilenberg.gi");
ReadPackage("HAP","lib/LieAlgebras/isLieHom.gi");
ReadPackage("HAP","lib/LieAlgebras/groupToLie.gi");
ReadPackage("HAP","lib/LieAlgebras/leibniz.gi");
ReadPackage("HAP","lib/LieAlgebras/LieTensorSquare.gi");
ReadPackage("HAP","lib/LieAlgebras/LieCover.gi");
ReadPackage("HAP","lib/LieAlgebras/LeibnizQuasiCover.gi");
ReadPackage("HAP","lib/LieAlgebras/LieExteriorSquare.gi");

if COMPILED=false then
################### MEAT AXE #######################################
ReadPackage("HAP","lib/FpGmodules/meataxe.gi");

################## POLYTOPAL COMPLEXES #############################
ReadPackage("HAP","lib/PolyComplexes/arrayOps.gi");
ReadPackage("HAP","lib/PolyComplexes/pureCubicalComplexes.gi");
ReadPackage("HAP","lib/PolyComplexes/sparseCubicalComplexes.gi");
ReadPackage("HAP","lib/PolyComplexes/chainComplexes.gi");
ReadPackage("HAP","lib/PolyComplexes/twoDimensional.gi");
ReadPackage("HAP","lib/PolyComplexes/threeDimensional.gi");
ReadPackage("HAP","lib/PolyComplexes/dvf.gi");
ReadPackage("HAP","lib/PolyComplexes/rips.gi");

ReadPackage("HAP","lib/PolyComplexes/simplicialComplexes.gi");
ReadPackage("HAP","lib/PolyComplexes/groupComplexes.gi");
ReadPackage("HAP","lib/PolyComplexes/cluster.gi");
ReadPackage("HAP","lib/PolyComplexes/metrics.gi");
ReadPackage("HAP","lib/PolyComplexes/graphviz.gi");
ReadPackage("HAP","lib/PolyComplexes/hap2chomp.gi");
ReadPackage("HAP","lib/PolyComplexes/filteredCubical.gi");
fi;


################## CATEGORY THEORY #################################
ReadPackage("HAP","lib/CategoryTheory/categories.gi");
ReadPackage("HAP","lib/CategoryTheory/commutativeDiagrams.gi");

################## CAT ONE GROUPS ##################################
ReadPackage("HAP","lib/CatGroups/CatConstructions.gi");
ReadPackage("HAP","lib/CatGroups/CatBasic.gi");
ReadPackage("HAP","lib/CatGroups/identities.gi");
ReadPackage("HAP","lib/CatGroups/algIdentities.gi");

################## G-OUTER GROUPS ##################################
ReadPackage("HAP","lib/GOuterGroups/goutergroup.gi");
ReadPackage("HAP","lib/GOuterGroups/homtogouter.gi");

if COMPILED=false then
################## SIMPLICIAL GROUPS ###############################
#ReadPackage("HAP","lib/SimplicialGroups/Kpinmap.gi");
ReadPackage("HAP","lib/SimplicialGroups/eilenbergMacLane.gi");
#ReadPackage("HAP","lib/SimplicialGroups/simplicialmap.gi");
ReadPackage("HAP","lib/SimplicialGroups/nerveOfCatOneGroup.gi");
ReadPackage("HAP","lib/SimplicialGroups/mooreComplex.gi");
ReadPackage("HAP","lib/SimplicialGroups/barresolution.gi");
ReadPackage("HAP","lib/SimplicialGroups/barcomplex.gi");
ReadPackage("HAP","lib/SimplicialGroups/chainComplexOfSimplicialGroup.gi");
ReadPackage("HAP","lib/SimplicialGroups/Kpin.gi");
ReadPackage("HAP","lib/SimplicialGroups/tensor2chains.gi");
#ReadPackage("HAP","lib/SimplicialGroups/lowerCentralSeriesOfCatOneGroup.gi");
ReadPackage("HAP","lib/SimplicialGroups/homotopyLowerCenterSeries.gi");
ReadPackage("HAP","lib/SimplicialGroups/quasiIsomorph.gi");


################## REGULAR CW_COMPLEXES ###############################
ReadPackage("HAP","lib/RegularCWComplexes/basicRegular.gi");
#ReadPackage("HAP","lib/RegularCWComplexes/contractAlt.gi");
ReadPackage("HAP","lib/RegularCWComplexes/fundamental.gi");
ReadPackage("HAP","lib/RegularCWComplexes/cocontract.gi");
fi;

if IsPackageMarkedForLoading("congruence","0.0") then
############### ARITHMETIC GROUPS#########################
 ReadPackage("HAP","lib/ArithmeticGroups/arithVarious.gi");
 ReadPackage("HAP","lib/ArithmeticGroups/arithmeticOps.gi");
 ReadPackage("HAP","lib/ArithmeticGroups/sl2zngens.gi");
 ReadPackage("HAP","lib/ArithmeticGroups/cplGTree.gi");
 ReadPackage("HAP","lib/ArithmeticGroups/resGTree.gi");
 ReadPackage("HAP","lib/ArithmeticGroups/sl2zres.gi");
ReadPackage("HAP","lib/ArithmeticGroups/sl2zresalt.gi");
# ReadPackage("HAP","lib/ArithmeticGroups/resDirectProd.gi");
######################################################
fi;

################## KNOTS ##########################################
ReadPackage("HAP","/lib/Knots/knotdata.gi");
ReadPackage("HAP","/lib/Knots/cubicalKnot.gi");

################## SPARSE ##########################################
ReadPackage("HAP","lib/Sparse/sparse.gi");


################# CRYSTALLOGRAPHIC GROUPS #####################
if IsPackageMarkedForLoading("HAPcryst","0.0") then
ReadPackage("HAP","lib/ArithmeticGroups/crystGbasis.gi");
ReadPackage("HAP","lib/ArithmeticGroups/crystGbasisnew.gi");
ReadPackage("HAP","lib/ArithmeticGroups/crystVarious.gi");
ReadPackage("HAP","lib/ArithmeticGroups/crystGcomplex.gi");
ReadPackage("HAP","lib/ArithmeticGroups/crystGcomplexnew.gi");
fi;



################## HAP PRIME ##################################
if IsPackageMarkedForLoading("singular","06.07.23") then
ReadPackage("HAP","lib/HapPrime/singular.gi");
ReadPackage("HAP","lib/HapPrime/rings.gi");
ReadPackage("HAP","lib/HapPrime/ringhomomorphism.gi");
ReadPackage("HAP","lib/HapPrime/gradedalgebra.gi");
ReadPackage("HAP","lib/HapPrime/polynomials.gi");
ReadPackage("HAP","lib/HapPrime/happrime.gi");
ReadPackage("HAP","lib/HapPrime/derivation.gi");
fi;


ReadPackage("HAP", "lib/Homology/BarCodes/barcode.gi");
ReadPackage("HAP","lib/TorsionSubcomplexes/TorsionSubcomplexes.gi");

SetInfoLevel(InfoWarning,1); #This is GAP's default level
