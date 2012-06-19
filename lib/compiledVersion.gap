#(C) Graham Ellis, 2005-2006

tmppath:=HAP_ROOT{[1..Length(HAP_ROOT)-4]};

########## COMPILED GAP FUNCTIONS ###################################
LoadDynamicModule(
Concatenation(tmppath,"lib/CompiledGAP/Compiled/absint.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/CompiledGAP/Compiled/signint.so"));

########## FREE G MODULES ###########################################

LoadDynamicModule(
Concatenation(tmppath,"lib/FreeGmodules/Compiled/tietze.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/FreeGmodules/Compiled/wordOperations.so"));

######### FPG MODULES ###############################################
LoadDynamicModule(
Concatenation(tmppath,"lib/FpGmodules/Compiled/fpgbasics.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/FpGmodules/Compiled/resfpgmod.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/FpGmodules/Compiled/homs.so"));

######### MEATAXE ##################################################
LoadDynamicModule(
Concatenation(tmppath,"lib/FpGmodules/Compiled/meataxe.so"));



##################### NONABELIAN TENSOR #############################
LoadDynamicModule(
Concatenation(tmppath,"lib/NonabelianTensor/Compiled/tensorSquare.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/NonabelianTensor/Compiled/symmetricSquare.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/NonabelianTensor/Compiled/tensorPair.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/NonabelianTensor/Compiled/exteriorProduct.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/NonabelianTensor/Compiled/SBG.so"));

if LoadPackage("nq")=true then

LoadDynamicModule(
Concatenation(tmppath,"lib/NonabelianTensor/Compiled/epiNilGrp.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/NonabelianTensor/Compiled/multNilGrp.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/NonabelianTensor/Compiled/tensorSquareInf.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/NonabelianTensor/Compiled/symmetricSquareInf.so"));
fi;


########## RESOLUTIONS ##############################################

LoadDynamicModule(
Concatenation(tmppath,"lib/Resolutions/Compiled/resFiniteGroup.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Resolutions/Compiled/resGeneric.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Resolutions/Compiled/presentation.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Resolutions/Compiled/resInfSubgroup.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Resolutions/Compiled/resSubgroup.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Resolutions/Compiled/resSmallFpGroup.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Resolutions/Compiled/resAbGroup.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Resolutions/Compiled/resAspherical.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Resolutions/Compiled/pseudoLists.so"));

if LoadPackage("aclib")=true then

LoadDynamicModule(
Concatenation(tmppath,"lib/Resolutions/Compiled/resACgroup.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Resolutions/Compiled/resACquotient.so"));

fi;


########## RESOLUTIONS MOD P #########################################

LoadDynamicModule(
Concatenation(tmppath,"lib/ResolutionsModP/Compiled/resPrimeGroup.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/ResolutionsModP/Compiled/poincare.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/ResolutionsModP/Compiled/ranksPrimeGroup.so"));


######### PERTURBATIONS #############################################

LoadDynamicModule(
Concatenation(tmppath,"lib/Perturbations/Compiled/resFiniteExt.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Perturbations/Compiled/resNormalSer.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Perturbations/Compiled/twistedTensorProduct.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Perturbations/Compiled/resDirectProd.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Perturbations/Compiled/resFiniteDirectProd.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Perturbations/Compiled/resExtension.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Perturbations/Compiled/resSubNormSeries.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Perturbations/Compiled/freeRes.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Perturbations/Compiled/dutour.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Perturbations/Compiled/contractibleSL2Zcomplex.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Perturbations/Compiled/filteredChainComplex.so"));





#################### ARTIN COXETER ##################################
LoadDynamicModule(
Concatenation(tmppath,"lib/ArtinCoxeter/Compiled/diagrams.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/ArtinCoxeter/Compiled/resArtin.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/ArtinCoxeter/Compiled/coxeterWythoff.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/ArtinCoxeter/Compiled/noncrossing.so"));



######### FUNCTORS ##################################################

LoadDynamicModule(
Concatenation(tmppath,"lib/Functors/Compiled/equiChainMap.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Functors/Compiled/modularEquiChainMap.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Functors/Compiled/homToZ.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Functors/Compiled/primePartDerived.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Functors/Compiled/tensorWithZmodP.so"));

#LoadDynamicModule( Concatenation(tmppath,"lib/Functors/Compiled/tensorWithZModule.so"));
ReadPackage("HAP", "lib/Functors/Compiled/tensorWithZModule.gi");

LoadDynamicModule(
Concatenation(tmppath,"lib/Functors/Compiled/tensorWithTwistedZmodP.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Functors/Compiled/tensorWithTwistedZ.so"));


LoadDynamicModule(
Concatenation(tmppath,"lib/Functors/Compiled/tensorWithZ.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Functors/Compiled/various.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Functors/Compiled/permMatrix.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Functors/Compiled/homToZmodule.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Functors/Compiled/tensorWithRationals.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Functors/Compiled/homToZmodP.so"));


#####################################################################


##################### HOMOLOGY ######################################
LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/cocycleCondition.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/integralHomology.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/lefschetz.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/integralHomologyObj.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/modularHomology.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/modularHomologyVectSpace.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/homology.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/groupHomology.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/groupCohomology.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/integralCohomology.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/integralCohomologyObj.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/cohomology.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/syzygy.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/cycles.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/isSuperperfect.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/modularCohomology.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/solutionsMat.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Homology/Compiled/persistent.so"));


######### COHOMOLOGY RINGS ##########################################

LoadDynamicModule(
Concatenation(tmppath,
"lib/Rings/Compiled/cocycleChainMap.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Rings/Compiled/cupProduct.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Rings/Compiled/intCoh.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Rings/Compiled/integralGens.so"));


################### POLYMAKE #######################################
LoadDynamicModule(
Concatenation(tmppath,"lib/Polymake/Compiled/aspherical.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Polymake/Compiled/polyGens.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Polymake/Compiled/stabilizer.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Polymake/Compiled/polyFaces.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Polymake/Compiled/orbitPoly.so"));


################### POLYCYLIC ######################################
LoadDynamicModule(
Concatenation(tmppath,"lib/Polycyclic/Compiled/resAbPcpGroup.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/Polycyclic/Compiled/resNilpotentPcpGrp.so"));


################## MOD P RINGS #####################################
LoadDynamicModule(
Concatenation(tmppath,"lib/ModPRings/Compiled/record.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/ModPRings/Compiled/recordPart1.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/ModPRings/Compiled/recordPartII.so"));

################### GRAPHS OF GROUPS ###############################
LoadDynamicModule(
Concatenation(tmppath, "lib/GraphsOfGroups/Compiled/graphs.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/GraphsOfGroups/Compiled/resGraph.so"));

LoadDynamicModule(
Concatenation(tmppath, "lib/GraphsOfGroups/Compiled/graphOfResolutions.so"));

################### SIMPLICIAL GROUPS ##############################
LoadDynamicModule(
Concatenation(tmppath, "lib/SimplicialGroups/Compiled/nerveCat1Group.so"));

LoadDynamicModule(
Concatenation(tmppath, "lib/SimplicialGroups/Compiled/mooreComplex.so"));

LoadDynamicModule(
Concatenation(tmppath, "lib/SimplicialGroups/Compiled/barresolution.so"));

LoadDynamicModule(
Concatenation(tmppath, "lib/SimplicialGroups/Compiled/chaincomplexofsimplicialgroup.so"));

LoadDynamicModule(
Concatenation(tmppath, "lib/SimplicialGroups/Compiled/Kpin.so"));

LoadDynamicModule(
Concatenation(tmppath, "lib/SimplicialGroups/Compiled/eilen_mac.so"));

LoadDynamicModule(
Concatenation(tmppath, "lib/SimplicialGroups/Compiled/tensor2chains.so"));

LoadDynamicModule(
Concatenation(tmppath, "lib/SimplicialGroups/Compiled/barcomplex.so"));

################## POLYTOPAL COMPLEXES ################################
LoadDynamicModule(
Concatenation(tmppath,"lib/PolyComplexes/Compiled/arrayOps.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/PolyComplexes/Compiled/pureCubicalComplexes.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/PolyComplexes/Compiled/chainComplexes.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/PolyComplexes/Compiled/twoDimensional.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/PolyComplexes/Compiled/threeDimensional.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/PolyComplexes/Compiled/dvf.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/PolyComplexes/Compiled/rips.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/PolyComplexes/Compiled/simplicialComplexes.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/PolyComplexes/Compiled/groupComplexes.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/PolyComplexes/Compiled/cluster.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/PolyComplexes/Compiled/hap2chomp.so"));


#####################REGULAR CW SPACES ##############################
LoadDynamicModule(
Concatenation(tmppath,"lib/RegularCWSpaces/Compiled/basicRegular.so"));

LoadDynamicModule(
Concatenation(tmppath,"lib/RegularCWSpaces/Compiled/contractAlt.so"));

