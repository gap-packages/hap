#############################################################################
##
##  init.g                  HAP library                  Graham Ellis 
##
############################################################################

if not IsBound(HapGlobalDeclarationsAreAlreadyLoaded) then
ReadPackage("HAP","lib/hap.gd");

#ReadPackage("HAP","lib/Objectifications/types.gi");
#ReadPackage("HAP","lib/TopologicalSpaces/topTypes.gd");

HapGlobalDeclarationsAreAlreadyLoaded:=true;
fi;


