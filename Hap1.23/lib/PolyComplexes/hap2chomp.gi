#(C) Graham Ellis 2009

#####################################################################
#####################################################################
InstallGlobalFunction(BinaryArrayToTextFile,
function(file,AA)
local
	i,j,k,l,m,A, AppendTo, PrintTo;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

if IsHapPureCubicalComplex(AA) then
A:=AA!.binaryArray;
else
A:=AA;
fi;

PrintTo(file,"");

######
if ArrayDimension(A)=1 then 
for i in [1..ArrayDimensions(A)[1]] do 
if A[i]=1 then
AppendTo(file,
[i,i+1],"\n");
fi;
od;
return true;
fi;
######

######
if ArrayDimension(A)=2 then
for i in [1..ArrayDimensions(A)[2]] do
for j in [1..ArrayDimensions(A)[1]] do
if A[i][j]=1 then
AppendTo(file,
Filtered(String([i,i+1]),a->not a=' '),
"x",
Filtered(String([j,j+1]),a->not a=' '),
"\n");
fi;
od;
od;
return true;
fi;
######

######
if ArrayDimension(A)=3 then
for i in [1..ArrayDimensions(A)[3]] do
for j in [1..ArrayDimensions(A)[2]] do
for k in [1..ArrayDimensions(A)[1]] do
if A[i][j][k]=1 then
AppendTo(file,
Filtered(String([i,i+1]),a->not a=' '),
"x",
Filtered(String([j,j+1]),a->not a=' '),
"x",
Filtered(String([k,k+1]),a->not a=' '),
"\n");
fi;
od;
od;
od;
return true;
fi;
######

######
if ArrayDimension(A)=4 then
for i in [1..ArrayDimensions(A)[4]] do
for j in [1..ArrayDimensions(A)[3]] do
for k in [1..ArrayDimensions(A)[2]] do
for l in [1..ArrayDimensions(A)[1]] do
if A[i][j][k][l]=1 then
AppendTo(file,
Filtered(String([i,i+1]),a->not a=' '),
"x",
Filtered(String([j,j+1]),a->not a=' '),
"x",
Filtered(String([k,k+1]),a->not a=' '),
"x",
Filtered(String([l,l+1]),a->not a=' '),
"\n");
fi;
od;
od;
od;
od;
return true;
fi;
######

######
if ArrayDimension(A)=5 then
for i in [1..ArrayDimensions(A)[5]] do
for j in [1..ArrayDimensions(A)[4]] do
for k in [1..ArrayDimensions(A)[3]] do
for l in [1..ArrayDimensions(A)[2]] do
for m in [1..ArrayDimensions(A)[1]] do
if A[i][j][k][l][m]=1 then
AppendTo(file,
Filtered(String([i,i+1]),a->not a=' '),
"x",
Filtered(String([j,j+1]),a->not a=' '),
"x",
Filtered(String([k,k+1]),a->not a=' '),
"x",
Filtered(String([l,l+1]),a->not a=' '),
"x",
Filtered(String([m,m+1]),a->not a=' '),
"\n");
fi;
od;
od;
od;
od;
od;
return true;
fi;
######




end);
#####################################################################
#####################################################################




