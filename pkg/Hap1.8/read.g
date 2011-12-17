#####################################################################
##
##  read.g                  HAP library               Graham Ellis 
##
#####################################################################
#ReadPackage("HAP","lib/hap.gd");


HAPconstant:=2;	
SetInfoLevel(InfoWarning,0); #We shouldn't really do this!

ReadPackage("HAP","boolean");
#ReadPackage("HAP", "lib/TitlePage/title.gap");
ReadPackage("HAP", "lib/TitlePage/copyright.gap");
#ReadPackage("HAP", "lib/TitlePage/makeHapMan.gi");


################### POLYCYLIC COMMANDS ##############################
## Most functions should work on pcp groups if the polycyclic package 
## is installed. Otherwise we need to give a meaning to certain commands
## defined in the polycyclic package.
Bool:=LoadPackage("polycyclic","0.0",false);
if Bool=fail then
NaturalHomomorphism:=NaturalHomomorphismByNormalSubgroupNC;
IsPcpGroup:=function(G);return false;end;
Collector:=function(x);return fail; end;
PcpGroupByCollector:=function(x);return fail; end;
Igs:=function(x);return fail; end;
GenExpList:=function(x);return fail; end;
HeisenbergPcpGroup:=function(x);return fail; end;
Pcp:=function(G,D); return fail; end;
IsAlmostCrystallographic:=function(G); return fail; end;
GeneratorsOfPcp:=function(G); return fail; end;
IsomorphismPcpGroup:=function(G);return fail;end;
fi;
################## POLYCICLIC COMMANDS DONE #########################

################# NQ COMMANDS #######################################
Bool:=LoadPackage("nq","0.0",false);
if Bool=fail then
NqEpimorphismNilpotentQuotient:=function(G); return fail; end;
fi;
################# NQ COMMANDS DONE ###############################

################# SIMPHOM COMMANDS ##################################
Bool:=LoadPackage("homology","0.0",false);
if Bool=fail then
SMInvariantFactors:=function(M); return fail; end;
else SetInfoLevel(InfoHomology,0);
ReadPackage("HAP", "lib/Homology/probHomology.gi");
fi;
################ SIMPHOM COMMANDS DONR ##############################

################# EDIM COMMANDS #######################################
Bool:=LoadPackage("edim","0.0",false);
if Bool=fail then
ElementaryDivisorsPPartRk:=function(G); return fail; end;
fi;
################# EDIM COMMANDS DONE ###############################

################# GAPDOC COMMANDS #######################################
Bool:=LoadPackage("gapdoc","0.0",false);
if Bool=fail then
MakeGAPDocDoc:=function(G); return fail; end;
fi;
################# GAPDOC COMMANDS DONE ###############################

ReadPackage("HAP", "lib/TitlePage/makeHapMan.gi");

################# OBJECTIFICATIONS ###############################
#ReadPackage("HAP", "lib/Objectifications/types.gi");
ReadPackage("HAP", "lib/Objectifications/basicMethods.gi");
################# OBJECTIFICATIONS DONE ##########################

if COMPILED=true then
ReadPackage("HAP","lib/compiledVersion.gap");
fi;


if COMPILED=false then
##################### FREE G MODULES ################################
ReadPackage("HAP", "lib/FreeGmodules/wordOperations.gi");
ReadPackage("HAP", "lib/FreeGmodules/tietze.gi");

##################### FP G MODULES ##################################
ReadPackage("HAP", "lib/FpGmodules/fpgbasics.gi");
ReadPackage("HAP", "lib/FpGmodules/resfpgmod.gi");
ReadPackage("HAP", "lib/FpGmodules/homs.gi");

##################### NONABELIAN TENSOR #############################
ReadPackage("HAP", "lib/NonabelianTensor/tensorSquare.gi");
ReadPackage("HAP", "lib/NonabelianTensor/tensorPair.gi");
ReadPackage("HAP", "lib/NonabelianTensor/exteriorProduct.gi");
ReadPackage("HAP", "lib/NonabelianTensor/SBG.gi");

if LoadPackage("nq","0.0",false)=true then
ReadPackage("HAP", "lib/NonabelianTensor/epiNilGrp.gi");
ReadPackage("HAP", "lib/NonabelianTensor/multNilGrp.gi");
ReadPackage("HAP", "lib/NonabelianTensor/tensorSquareInf.gi");
fi;

##################### RESOLUTIONS ###################################
ReadPackage("HAP", "lib/Resolutions/resAspherical.gi");
ReadPackage("HAP", "lib/Resolutions/resAbGroup.gi");
ReadPackage("HAP", "lib/Resolutions/resFiniteGroup.gi");
ReadPackage("HAP", "lib/Resolutions/resSmallFpGroup.gi");
ReadPackage("HAP", "lib/Resolutions/presentation.gi");
ReadPackage("HAP", "lib/Resolutions/resSubgroup.gi");
ReadPackage("HAP", "lib/Resolutions/resInfSubgroup.gi");
ReadPackage("HAP", "lib/Resolutions/pseudoLists.gi");

if LoadPackage("aclib","0.0",false)=true then
ReadPackage("HAP", "lib/Resolutions/resACgroup.gi");
ReadPackage("HAP", "lib/Resolutions/resACquotient.gi");
fi;


##################### RESOLUTIONS MOD P #############################
ReadPackage("HAP", "lib/ResolutionsModP/resPrimeGroup.gi");
ReadPackage("HAP", "lib/ResolutionsModP/ranksPrimeGroup.gi");
ReadPackage("HAP", "lib/ResolutionsModP/poincare.gi");


##################### FUNCTORS ######################################
ReadPackage("HAP", "lib/Functors/permMatrix.gi");
ReadPackage("HAP", "lib/Functors/homToZmodule.gi");
ReadPackage("HAP", "lib/Functors/tensorWithZ.gi");
ReadPackage("HAP", "lib/Functors/tensorWithZmodP.gi");
ReadPackage("HAP", "lib/Functors/various.gi");
ReadPackage("HAP", "lib/Functors/equiChainMap.gi");
ReadPackage("HAP", "lib/Functors/primePartDerived.gi");
ReadPackage("HAP", "lib/Functors/homToZ.gi");
ReadPackage("HAP", "lib/Functors/tensorWithRationals.gi");
ReadPackage("HAP", "lib/Functors/homToZmodP.gi");


##################### HOMOLOGY ######################################
ReadPackage("HAP", "lib/Homology/integralHomology.gi");
ReadPackage("HAP", "lib/Homology/modularHomology.gi");
ReadPackage("HAP", "lib/Homology/homology.gi");
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




##################### PERTURBATIONS #################################
ReadPackage("HAP", "lib/Perturbations/resExtension.gi");
ReadPackage("HAP", "lib/Perturbations/resDirectProd.gi");
ReadPackage("HAP", "lib/Perturbations/twistedTensorProduct.gi");
ReadPackage("HAP", "lib/Perturbations/resFiniteExt.gi");
ReadPackage("HAP", "lib/Perturbations/resNormalSer.gi");
ReadPackage("HAP", "lib/Perturbations/resFiniteDirectProd.gi");
ReadPackage("HAP", "lib/Perturbations/resSubNormSeries.gi");



#################### ARTIN COXETER ##################################
ReadPackage("HAP", "lib/ArtinCoxeter/diagrams.gi");
ReadPackage("HAP", "lib/ArtinCoxeter/resArtin.gi");

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

fi;

################### GRAPHS OF GROUPS ###############################
ReadPackage("HAP", "lib/GraphsOfGroups/graphs.gi");
ReadPackage("HAP", "lib/GraphsOfGroups/resGraph.gi");

################### TEST ###########################################
ReadPackage("HAP", "test/test.gap");

################### STREAMS ########################################
ReadPackage("HAP","lib/Streams/streams.gi");


################### RESOLUTIONS (CONTD) ############################
ReadPackage("HAP","lib/Resolutions/cayley.gi");

################### Lie Algebras ###################################
ReadPackage("HAP","lib/LieAlgebras/chevalleyEilenberg.gi");
ReadPackage("HAP","lib/LieAlgebras/isLieHom.gi");
ReadPackage("HAP","lib/LieAlgebras/groupToLie.gi");
ReadPackage("HAP","lib/LieAlgebras/leibniz.gi");

if COMPILED=false then
################### MEAT AXE #######################################
ReadPackage("HAP","lib/FpGmodules/meataxe.gi");
fi;

################## TDA #############################################
ReadPackage("HAP","lib/TDA/tda.gi");
ReadPackage("HAP","lib/TDA/cont.gi");
