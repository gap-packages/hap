
#################################################
#################################################
TestHapBook:=function()
local dir;

dir:=DirectoriesPackageLibrary( "HAP","")[1];
dir:=Filename(dir,"");
dir:=Concatenation(dir,"test/files/");


if Test(Concatenation(dir,"1.1.1.tst")) then Print(true," 1.1.1\n");
else Print("error in 1.1.1\n"); fi;

if Test(Concatenation(dir,"1.2.1.tst")) then Print(true," 1.2.1\n");
else Print("error in 1.2.1\n"); fi;

if Test(Concatenation(dir,"1.2.2.tst")) then Print(true," 1.2.2\n");
else Print("error in 1.2.2\n"); fi;

if Test(Concatenation(dir,"1.2.4.tst")) then Print(true," 1.2.4\n");
else Print("error in 1.2.4\n"); fi;

if Test(Concatenation(dir,"1.2.5.tst")) then Print(true," 1.2.5\n");
else Print("error in 1.2.5\n"); fi;

if Test(Concatenation(dir,"1.3.1.tst")) then Print(true," 1.3.1\n");
else Print("error in 1.3.1\n"); fi;

if Test(Concatenation(dir,"1.3.2.tst")) then Print(true," 1.3.2\n");
else Print("error in 1.3.2\n"); fi;

if Test(Concatenation(dir,"1.3.3.tst")) then Print(true," 1.3.3\n");
else Print("error in 1.3.3\n"); fi;

if Test(Concatenation(dir,"1.3.4.tst")) then Print(true," 1.3.4\n");
else Print("error in 1.3.4\n"); fi;

if Test(Concatenation(dir,"1.4.1.tst")) then Print(true," 1.4.1\n");
else Print("error in 1.4.1\n"); fi;

if Test(Concatenation(dir,"1.4.3.tst")) then Print(true," 1.4.3\n");
else Print("error in 1.4.3\n"); fi;

if Test(Concatenation(dir,"1.5.3.tst")) then Print(true," 1.5.3\n");
else Print("error in 1.5.3\n"); fi;

if Test(Concatenation(dir,"1.6.1.tst")) then Print(true," 1.6.1\n");
else Print("error in 1.6.1\n"); fi;

if Test(Concatenation(dir,"1.6.2.tst")) then Print(true," 1.6.2\n");
else Print("error in 1.6.2\n"); fi;

if Test(Concatenation(dir,"1.6.4.tst")) then Print(true," 1.6.4\n");
else Print("error in 1.6.4\n"); fi;

if Test(Concatenation(dir,"1.6.5.tst")) then Print(true," 1.6.5\n");
else Print("error in 1.6.5\n"); fi;

if Test(Concatenation(dir,"1.6.7.tst")) then Print(true," 1.6.7\n");
else Print("error in 1.6.7\n"); fi;

if Test(Concatenation(dir,"1.7.2.tst")) then Print(true," 1.7.2\n");
else Print("error in 1.7.2\n"); fi;

if Test(Concatenation(dir,"1.7.3.tst")) then Print(true," 1.7.3\n");
else Print("error in 1.7.3\n"); fi;

if Test(Concatenation(dir,"1.8.1.tst")) then Print(true," 1.8.1\n");
else Print("error in 1.8.1\n"); fi;

if Test(Concatenation(dir,"1.8.2.tst")) then Print(true," 1.8.2\n");
else Print("error in 1.8.2\n"); fi;

if Test(Concatenation(dir,"1.9.1.tst")) then Print(true," 1.9.1\n");
else Print("error in 1.9.1\n"); fi;

if Test(Concatenation(dir,"2.1.4.tst")) then Print(true," 2.1.4\n");
else Print("error in 2.1.4\n"); fi;

if Test(Concatenation(dir,"2.1.5.tst")) then Print(true," 2.1.5\n");
else Print("error in 2.1.5\n"); fi;

if Test(Concatenation(dir,"2.1.6.tst")) then Print(true," 2.1.6\n");
else Print("error in 2.1.6\n"); fi;


end;
#################################################
#################################################

