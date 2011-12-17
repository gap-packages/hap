#############################################################################
##
##  init.g                  HAP library                  Graham Ellis 
##
############################################################################

if not IsBound(HapGlobalDeclarationsAreAlreadyLoaded) then
ReadPackage("HAP","lib/hap.gd");
HapGlobalDeclarationsAreAlreadyLoaded:=true;
MakeReadOnlyGlobal("HapGlobalDeclarationsAreAlreadyLoaded");
fi;

ReadPackage("HAP","/lib/externalSoftware.gap");


