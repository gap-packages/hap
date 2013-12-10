
DeclareOperation("Homology",[IsHapCrossedModule,IsInt]);
DeclareOperation("Homology",[IsHapCatOneGroup,IsInt]);
DeclareOperation("Homology",[IsHapSimplicialGroup,IsInt]);

DeclareOperation("PersistentHomology",[IsHapCrossedModule,IsInt]);
DeclareOperation("PersistentHomology",[IsHapCatOneGroup,IsInt]);

############ Crossed Module ####################################
DeclareOperation("HomotopyGroup",[IsHapCrossedModule,IsInt]);
DeclareAttribute("Size",IsHapCrossedModule);
DeclareAttribute("Order",IsHapCrossedModule);
DeclareGlobalFunction("CatOneGroupByCrossedModule");
DeclareGlobalFunction("CrossedModuleByCatOneGroup");
DeclareGlobalFunction("CrossedModuleByAutomorphismGroup");
DeclareGlobalFunction("CrossedModuleByNormalSubgroup");
DeclareGlobalFunction("HomotopyLowerCentralSeriesOfCrossedModule");
DeclareGlobalFunction("PersistentHomologyOfCrossedModule");
DeclareGlobalFunction("HomotopyLowerCentralSeries");
DeclareGlobalFunction("HomotopyCrossedModule");
DeclareGlobalFunction("NumberSmallCrossedModules");
DeclareGlobalFunction("SmallCrossedModule");
DeclareGlobalFunction("IdCrossedModule");
DeclareGlobalFunction("IsomorphismCrossedModules");

################## Cat One Group ###############################
DeclareGlobalFunction("NerveOfCatOneGroup");
##DeclareGlobalFunction("XmodToHAP");
##DeclareGlobalFunction("SubQuasiIsomorph");
##DeclareGlobalFunction("QuotientQuasiIsomorph");
##DeclareGlobalFunction("QuasiIsomorph");
DeclareGlobalFunction("IsomorphismCatOneGroups");
DeclareGlobalFunction("CatOneGroupsByGroup");
DeclareGlobalFunction("HomotopyCatOneGroup");
DeclareGlobalFunction("IdCatOneGroup");
DeclareGlobalFunction("SmallCatOneGroup");
DeclareGlobalFunction("NumberSmallCatOneGroups");

#################### Two Types ##########################
DeclareGlobalFunction("NumberSmallQuasiCatOneGroups");
DeclareGlobalFunction("SmallQuasiCatOneGroup");
DeclareGlobalFunction("IdQuasiCatOneGroup");
DeclareGlobalFunction("NumberSmallQuasiCrossedModules");
DeclareGlobalFunction("SmallQuasiCrossedModule");
DeclareGlobalFunction("IdQuasiCrossedModule");



DeclareGlobalFunction("BarResolutionEquivalence");
DeclareGlobalFunction("BarComplexEquivalence");
DeclareGlobalFunction("ChainComplexOfSimplicialGroup");
##DeclareGlobalFunction("SylowSubgroupOfCatOneGroup");
DeclareGlobalFunction("TensorProductOfChainComplexes");
DeclareGlobalFunction("EilenbergMacLaneSimplicialGroup");






