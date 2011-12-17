#####################################################################
##
##  read.g                  HAP library               Graham Ellis 
##
#####################################################################

COMPILED:=false;   

if COMPILED=true then
ReadPackage("HAP","lib/compiledVersion.gap");
fi;

#ReadPackage("HAP", "lib/TitlePage/title.gap");
ReadPackage("HAP", "lib/TitlePage/copyright.gap");

##################### FREE G MODULES ################################
if COMPILED=false then
ReadPackage("HAP", "lib/FreeGmodules/wordOperations.gi");
ReadPackage("HAP", "lib/FreeGmodules/tietze.gi");
fi;

##################### NONABELIAN TENSOR #############################
ReadPackage("HAP", "lib/NonabelianTensor/tensorSquare.gi");
ReadPackage("HAP", "lib/NonabelianTensor/exteriorProduct.gi");

##################### RESOLUTIONS ###################################
ReadPackage("HAP", "lib/Resolutions/resAspherical.gi");
ReadPackage("HAP", "lib/Resolutions/resAbGroup.gi");
if COMPILED=false then
ReadPackage("HAP", "lib/Resolutions/resFiniteGroup.gi");
ReadPackage("HAP", "lib/Resolutions/resSmallFpGroup.gi");
ReadPackage("HAP", "lib/Resolutions/presentation.gi");
ReadPackage("HAP", "lib/Resolutions/resSubgroup.gi");
ReadPackage("HAP", "lib/Resolutions/resInfSubgroup.gi");
fi;

##################### RESOLUTIONS MOD P #############################
if COMPILED=false then
ReadPackage("HAP", "lib/ResolutionsModP/resPrimeGroup.gi");
fi;

##################### FUNCTORS ######################################
ReadPackage("HAP", "lib/Functors/permMatrix.gi");
ReadPackage("HAP", "lib/Functors/homToZmodule.gi");
if COMPILED=false then
ReadPackage("HAP", "lib/Functors/tensorWithZ.gi");
ReadPackage("HAP", "lib/Functors/tensorWithZmodP.gi");
ReadPackage("HAP", "lib/Functors/various.gi");
ReadPackage("HAP", "lib/Functors/equiChainMap.gi");
ReadPackage("HAP", "lib/Functors/primePartDerived.gi");
ReadPackage("HAP", "lib/Functors/homToZ.gi");
fi;

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

##################### PERTURBATIONS #################################
if COMPILED=false then
ReadPackage("HAP", "lib/Perturbations/resDirectProd.gi");
ReadPackage("HAP", "lib/Perturbations/twistedTensorProduct.gi");
ReadPackage("HAP", "lib/Perturbations/resFiniteExt.gi");
ReadPackage("HAP", "lib/Perturbations/resNormalSer.gi");
fi;

#################### ARTIN COXETER ##################################
ReadPackage("HAP", "lib/ArtinCoxeter/diagrams.gi");
ReadPackage("HAP", "lib/ArtinCoxeter/resArtin.gi");

#################### COHOMOLOGY RINGS ###############################
if COMPILED=false then
ReadPackage("HAP", "lib/Rings/intCoh.gi");
ReadPackage("HAP", "lib/Rings/cocycleChainMap.gi");
ReadPackage("HAP", "lib/Rings/cupProduct.gi");
ReadPackage("HAP", "lib/Rings/integralGens.gi");
fi;

################### POLYMAKE #######################################
ReadPackage("HAP", "lib/Polymake/aspherical.gi");
ReadPackage("HAP", "lib/Polymake/polyGens.gi");
ReadPackage("HAP", "lib/Polymake/stabilizer.gi");
ReadPackage("HAP", "lib/Polymake/polyFaces.gi");
