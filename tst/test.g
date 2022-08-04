#################################################
#################################################
InstallGlobalFunction(HapFile,
function(str)
local exfile;

exfile:=DirectoriesPackageLibrary("HAP","tst/examples");;
exfile:=Filename(exfile,str);;

return exfile;
end);
#################################################
#################################################

#################################################
#################################################
InstallGlobalFunction(HapExample,
function(str)
local exfile, input, linefile, x;

exfile:=DirectoriesPackageLibrary("HAP","tst/examples");;
exfile:=Filename(exfile, str);;

input:=InputTextFile(exfile);
x:=ReadLine(input);
while not x=fail do
Print("gap> ");
Print(x);
x:=ReadLine(input);
od;
Print("\n\n");
Read(exfile);
end);
#################################################
#################################################






#################################################
#################################################
InstallGlobalFunction(TestHapBook,
function()
ReadPackage("HAP","tst/testallextra.g");
end);
#################################################
#################################################

#################################################
#################################################
InstallGlobalFunction(TestHap,
function()
ReadPackage("HAP","tst/testall.g");
end);
#################################################
#################################################

#################################################
#################################################
InstallGlobalFunction(TestHapQuick,
function()
ReadPackage("HAP","tst/testquick.g");
end);
#################################################
#################################################


