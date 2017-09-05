
#################################################
#################################################
TestHapBook:=function()
local dir;

dir:=DirectoriesPackageLibrary( "HAP","")[1];;
dir:=Filename(dir,"");;
dir:=Concatenation(dir,"test/files/");;


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

if Test(Concatenation(dir,"2.2.1.tst")) then Print(true," 2.2.1\n");
else Print("error in 2.2.1\n"); fi;

if Test(Concatenation(dir,"2.4.2.tst")) then Print(true," 2.4.2\n");
else Print("error in 2.4.2\n"); fi;

if Test(Concatenation(dir,"2.4.3.tst")) then Print(true," 2.4.3\n");
else Print("error in 2.4.3\n"); fi;

if Test(Concatenation(dir,"2.4.4.tst")) then Print(true," 2.4.4\n");
else Print("error in 2.4.4\n"); fi;

if Test(Concatenation(dir,"2.4.5.tst")) then Print(true," 2.4.5\n");
else Print("error in 2.4.5\n"); fi;

if Test(Concatenation(dir,"2.4.6.tst")) then Print(true," 2.4.6\n");
else Print("error in 2.4.6\n"); fi;

if Test(Concatenation(dir,"2.4.8.tst")) then Print(true," 2.4.8\n");
else Print("error in 2.4.8\n"); fi;

if Test(Concatenation(dir,"2.4.9.tst")) then Print(true," 2.4.9\n");
else Print("error in 2.4.9\n"); fi;

if Test(Concatenation(dir,"2.5.1.tst")) then Print(true," 2.5.1\n");
else Print("error in 2.5.1\n"); fi;

if Test(Concatenation(dir,"2.6.1.tst")) then Print(true," 2.6.1\n");
else Print("error in 2.6.1\n"); fi;

if Test(Concatenation(dir,"2.6.2.tst")) then Print(true," 2.6.2\n");
else Print("error in 2.6.2\n"); fi;

if Test(Concatenation(dir,"2.6.3.tst")) then Print(true," 2.6.3\n");
else Print("error in 2.6.3\n"); fi;

if Test(Concatenation(dir,"2.7.1.tst")) then Print(true," 2.7.1\n");
else Print("error in 2.7.1\n"); fi;

if Test(Concatenation(dir,"2.7.1a.tst")) then Print(true," 2.7.1a\n");
else Print("error in 2.7.1a\n"); fi;

if Test(Concatenation(dir,"2.7.2.tst")) then Print(true," 2.7.2\n");
else Print("error in 2.7.2\n"); fi;

if Test(Concatenation(dir,"2.7.3.tst")) then Print(true," 2.7.3\n");
else Print("error in 2.7.3\n"); fi;

if Test(Concatenation(dir,"2.7.4.tst")) then Print(true," 2.7.4\n");
else Print("error in 2.7.4\n"); fi;

if Test(Concatenation(dir,"2.7.5.tst")) then Print(true," 2.7.5\n");
else Print("error in 2.7.5\n"); fi;

if Test(Concatenation(dir,"3.1.2.tst")) then Print(true," 3.1.2\n");
else Print("error in 3.1.2\n"); fi;

if Test(Concatenation(dir,"3.1.3.tst")) then Print(true," 3.1.3\n");
else Print("error in 3.1.3\n"); fi;

if Test(Concatenation(dir,"3.1.4.tst")) then Print(true," 3.1.4\n");
else Print("error in 3.1.4\n"); fi;

if Test(Concatenation(dir,"3.2.1.tst")) then Print(true," 3.2.1\n");
else Print("error in 3.2.1\n"); fi;

if Test(Concatenation(dir,"3.2.2.tst")) then Print(true," 3.2.2\n");
else Print("error in 3.2.2\n"); fi;

if Test(Concatenation(dir,"3.2.3.tst")) then Print(true," 3.2.3\n");
else Print("error in 3.2.3\n"); fi;

if Test(Concatenation(dir,"3.3.1.tst")) then Print(true," 3.3.1\n");
else Print("error in 3.3.1\n"); fi;

if Test(Concatenation(dir,"3.3.2.tst")) then Print(true," 3.3.2\n");
else Print("error in 3.3.2\n"); fi;

if Test(Concatenation(dir,"3.3.3.tst")) then Print(true," 3.3.3\n");
else Print("error in 3.3.3\n"); fi;

if Test(Concatenation(dir,"3.3.4.tst")) then Print(true," 3.3.4\n");
else Print("error in 3.3.4\n"); fi;

if Test(Concatenation(dir,"3.3.5.tst")) then Print(true," 3.3.5\n");
else Print("error in 3.3.5\n"); fi;

if Test(Concatenation(dir,"3.3.6.tst")) then Print(true," 3.3.6\n");
else Print("error in 3.3.6\n"); fi;

if Test(Concatenation(dir,"3.3.7.tst")) then Print(true," 3.3.7\n");
else Print("error in 3.3.7\n"); fi;

if Test(Concatenation(dir,"3.3.8.tst")) then Print(true," 3.3.8\n");
else Print("error in 3.3.8\n"); fi;

if Test(Concatenation(dir,"3.3.9.tst")) then Print(true," 3.3.9\n");
else Print("error in 3.3.9\n"); fi;

if Test(Concatenation(dir,"3.3.10.tst")) then Print(true," 3.3.10\n");
else Print("error in 3.3.10\n"); fi;

if Test(Concatenation(dir,"3.4.2.tst")) then Print(true," 3.4.2\n");
else Print("error in 3.4.2\n"); fi;

if Test(Concatenation(dir,"3.5.1.tst")) then Print(true," 3.5.1\n");
else Print("error in 3.5.1\n"); fi;

if Test(Concatenation(dir,"3.5.2.tst")) then Print(true," 3.5.2\n");
else Print("error in 3.5.2\n"); fi;

if Test(Concatenation(dir,"3.5.3.tst")) then Print(true," 3.5.3\n");
else Print("error in 3.5.3\n"); fi;

if Test(Concatenation(dir,"3.5.4.tst")) then Print(true," 3.5.4\n");
else Print("error in 3.5.4\n"); fi;

if Test(Concatenation(dir,"3.5.5.tst")) then Print(true," 3.5.5\n");
else Print("error in 3.5.5\n"); fi;

if Test(Concatenation(dir,"3.6.1.tst")) then Print(true," 3.6.1\n");
else Print("error in 3.6.1\n"); fi;

if Test(Concatenation(dir,"3.6.2.tst")) then Print(true," 3.6.2\n");
else Print("error in 3.6.2\n"); fi;

if Test(Concatenation(dir,"3.7.1.tst")) then Print(true," 3.7.1\n");
else Print("error in 3.7.1\n"); fi;

if Test(Concatenation(dir,"3.7.4.tst")) then Print(true," 3.7.4\n");
else Print("error in 3.7.4\n"); fi;

if Test(Concatenation(dir,"3.7.5.tst")) then Print(true," 3.7.5\n");
else Print("error in 3.7.5\n"); fi;

if Test(Concatenation(dir,"3.8.1.tst")) then Print(true," 3.8.1\n");
else Print("error in 3.8.1\n"); fi;

if Test(Concatenation(dir,"3.9.1.tst")) then Print(true," 3.9.1\n");
else Print("error in 3.9.1\n"); fi;

if Test(Concatenation(dir,"3.9.3.tst")) then Print(true," 3.9.3\n");
else Print("error in 3.9.3\n"); fi;

end;
#################################################
#################################################

