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

##I introduce the NC versions of PreImages...
if not IsBound( PreImagesNC ) then
    BindGlobal( "PreImagesNC", PreImages );
fi;
if not IsBound( PreImagesElmNC ) then
    BindGlobal( "PreImagesElmNC", PreImagesElm );
fi;
if not IsBound( PreImagesSetNC ) then
    BindGlobal( "PreImagesSetNC", PreImagesSet );
fi;
if not IsBound( PreImagesRepresentativeNC ) then
    BindGlobal( "PreImagesRepresentativeNC", PreImagesRepresentative );
fi;

ReadPackage("HAP","/lib/externalSoftware.gap");


