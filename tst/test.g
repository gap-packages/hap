#################################################
#################################################
InstallGlobalFunction(HapFile,
function(str)
local exfile, input, linefile, x;

exfile:=DirectoriesPackageLibrary( "HAP","")[1];;
exfile:=Filename(exfile,"");;
exfile:=Concatenation(exfile,"tst/examples/");;
exfile:=Concatenation(exfile,str);

return exfile;
end);
#################################################
#################################################

#################################################
#################################################
InstallGlobalFunction(HapExample,
function(str)
local exfile, input, linefile, x;

exfile:=DirectoriesPackageLibrary( "HAP","")[1];;
exfile:=Filename(exfile,"");;
exfile:=Concatenation(exfile,"tst/examples/");;
exfile:=Concatenation(exfile,str);

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
local dir;

dir:=DirectoriesPackageLibrary( "HAP","")[1];;
dir:=Filename(dir,"");;
dir:=Concatenation(dir,"tst/testallextra.g");;
Read(dir);

end);
#################################################
#################################################

#################################################
#################################################
InstallGlobalFunction(TestHap,
function()
local dir;

dir:=DirectoriesPackageLibrary( "HAP","")[1];;
dir:=Filename(dir,"");;
dir:=Concatenation(dir,"tst/testall.g");;
Read(dir);

end);
#################################################
#################################################


