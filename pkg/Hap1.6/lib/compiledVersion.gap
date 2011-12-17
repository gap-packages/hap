#(C) Graham Ellis, 2005-2006

########## FREE G MODULES ###########################################

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/FreeGmodules/Compiled/tietze.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/FreeGmodules/Compiled/wordOperations.so");


##################### NONABELIAN TENSOR #############################
LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/NonabelianTensor/Compiled/tensorSquare.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/NonabelianTensor/Compiled/exteriorProduct.so");

if LoadPackage("nq")=true then

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/NonabelianTensor/Compiled/epiNilGrp.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/NonabelianTensor/Compiled/multNilGrp.so");

fi;


########## RESOLUTIONS ##############################################

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Resolutions/Compiled/resFiniteGroup.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Resolutions/Compiled/presentation.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Resolutions/Compiled/resInfSubgroup.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Resolutions/Compiled/resSubgroup.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Resolutions/Compiled/resSmallFpGroup.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Resolutions/Compiled/resAbGroup.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Resolutions/Compiled/resAspherical.so");

if LoadPackage("aclib")=true then

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Resolutions/Compiled/resACgroup.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Resolutions/Compiled/resACquotient.so");

fi;


######### RESOLUTIONS MOD P #########################################

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/ResolutionsModP/Compiled/resPrimeGroup.so");


######### PERTURBATIONS #############################################

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Perturbations/Compiled/resFiniteExt.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Perturbations/Compiled/resNormalSer.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Perturbations/Compiled/twistedTensorProduct.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Perturbations/Compiled/resDirectProd.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Perturbations/Compiled/resFiniteDirectProd.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Perturbations/Compiled/resExtension.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Perturbations/Compiled/resSubNormSeries.so");



#################### ARTIN COXETER ##################################
LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/ArtinCoxeter/Compiled/diagrams.so");


LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/ArtinCoxeter/Compiled/resArtin.so");


######### FUNCTORS ##################################################

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Functors/Compiled/equiChainMap.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Functors/Compiled/homToZ.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Functors/Compiled/primePartDerived.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Functors/Compiled/tensorWithZmodP.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Functors/Compiled/tensorWithZ.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Functors/Compiled/various.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Functors/Compiled/permMatrix.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Functors/Compiled/homToZmodule.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Functors/Compiled/tensorWithRationals.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Functors/Compiled/homToZmodP.so");


#####################################################################


##################### HOMOLOGY ######################################
LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Homology/Compiled/cocycleCondition.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Homology/Compiled/integralHomology.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Homology/Compiled/modularHomology.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Homology/Compiled/homology.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Homology/Compiled/groupHomology.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Homology/Compiled/integralCohomology.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Homology/Compiled/cohomology.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Homology/Compiled/syzygy.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Homology/Compiled/cycles.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Homology/Compiled/isSuperperfect.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Homology/Compiled/modularCohomology.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Homology/Compiled/solutionsMat.so");

######### COHOMOLOGY RINGS ##########################################

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Rings/Compiled/cocycleChainMap.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Rings/Compiled/cupProduct.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Rings/Compiled/intCoh.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Rings/Compiled/integralGens.so");


################### POLYMAKE #######################################
LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Polymake/Compiled/aspherical.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Polymake/Compiled/polyGens.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Polymake/Compiled/stabilizer.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Polymake/Compiled/polyFaces.so");

################### POLYCYLIC ######################################
LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Polycyclic/Compiled/resAbPcpGroup.so");

LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/Polycyclic/Compiled/resNilpotentPcpGrp.so");

################## MOD P RINGS #####################################
LoadDynamicModule(
"/home/graham/pkg/Hap1.6/lib/ModPRings/Compiled/record.so");
