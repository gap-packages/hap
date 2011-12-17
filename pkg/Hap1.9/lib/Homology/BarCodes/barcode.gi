HAPBARCODE:=0;

#############################################################
InstallGlobalFunction(UniversalBarCode,
function(str,n,d)
local B, i, j, M, file,PolRing,x;

if not str in
["DerivedSeries", "LowerCentralSeries", "UpperCentralSeries", "PCentralSeries",
"PUpperCentralSeries"] then 
Print("The first argument must be a string equal to one of \"DerivedSeries\", \"LowerCentralSeries\", \"UpperCentralSeries\", \"PCentralseries\", \"PUpperCentralSeries\" \n");
return fail;
fi;

 
PolRing:=HapConstantPolRing;
x:=GeneratorsOfAlgebra(PolRing)[2];

file:=Concatenation("lib/Homology/BarCodes/",str,"_",String(n),".gi");
ReadPackage("HAP",file);

M:=HAPBARCODE[d];

B:=[];
for i in [1..Length(M)] do
B[i]:=[];
for j in [1..Length(M[i])] do
B[i][j]:=ValuePol(M[i][j][1],x)/ValuePol(M[i][j][2],x);
od;
od;

return B;
end);
#############################################################

#############################################################
InstallGlobalFunction(UniversalBarCodeEval,
function(BB,k)
local i,j,B;

B:=StructuralCopy(BB);
for i in [1..Length(B)] do
for j in [1..Length(B[i])] do
if not IsInt(B[i][j]) then
B[i][j]:=ExpansionOfRationalFunction(B[i][j],k)[k];
fi;
od;
od;

return B;
end);
#############################################################
