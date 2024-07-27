#(C) 2009 Graham Ellis

#####################################################################
#####################################################################
InstallGlobalFunction(ArrayToPureCubicalComplex,
function(AA,threshold)

# Inputs an array of integers where 0 represents black 
# and 753 represents white. It returns a pure cubical complex
# where all integers below the threshold are converted to 1, and
# higher integers are converted to 0, in the binary list f the cubical
# complex. 

local
        int2binary,A;

##############################
int2binary:=function(A)
local i;
if ArrayDimension(A)=1 then 
for i in [1..Length(A)] do
if A[i] <threshold then A[i]:=1; else
A[i]:=0;fi;
od;return A;
fi;
return List(A,x->int2binary(StructuralCopy(x)));
end;
##############################

A:=int2binary(AA);

return Objectify(HapPureCubicalComplex,
           rec(
           binaryArray:=A,
           properties:=[
           ["dimension",ArrayDimension(A)],
           ["arraySize",ArrayDimensions(A)]]
           ));
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(PureCubicalComplex,
function(AA)
# Inputs a binary array and returns a pure cubical complex.
local
        A;

A:=StructuralCopy(AA);

return Objectify(HapPureCubicalComplex,
           rec(
           binaryArray:=A,
           properties:=[
           ["dimension",ArrayDimension(A)],
           ["arraySize",ArrayDimensions(A)]]
           ));
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(FilteredPureCubicalComplex,
function(AA,BB)
# Inputs a binary array AA and a filtration array BB and returns a pure cubical complex.
local
        A, B;

A:=StructuralCopy(AA);
B:=StructuralCopy(BB);

return Objectify(HapFilteredPureCubicalComplex,
           rec(
           binaryArray:=A,
           filtration:=B,
           properties:=[
           ["dimension",ArrayDimension(A)],
           ["arraySize",ArrayDimensions(A)]]
           ));
end);
#####################################################################
#####################################################################



#####################################################################
#####################################################################
InstallOtherMethod(Dimension,
"Dimension of pure cubical complex",
[IsHapPureCubicalComplex],
function(f) return EvaluateProperty(f,"dimension");
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Dimension,
"Dimension of pure cubical complex",
[IsHapCubicalComplex],
function(f) return EvaluateProperty(f,"dimension");
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallOtherMethod(Size,
"Volume of a pure cubical complex",
[IsHapPureCubicalComplex],
function(f) return Sum(Flat(f!.binaryArray));
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Size,
"Volume of a pure cubical complex",
[IsHapCubicalComplex],
function(f) return Sum(Flat(f!.binaryArray));
end);
#####################################################################
#####################################################################




#####################################################################
#####################################################################
InstallMethod(EulerCharacteristic,
"Euler characteristic  of a pure  cubical complex",
[IsHapPureCubicalComplex],
function(M) local  D;
D:=ChainComplexOfCubicalComplex(
PureCubicalComplexToCubicalComplex(M),true);
return
Sum(List([0..Dimension(M)],i->((-1)^i)*D!.dimension(i)));
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallOtherMethod(EulerCharacteristic,
"Euler characteristic  of a cubical complex",
[IsHapCubicalComplex],
function(M) local  C;
C:=ChainComplexOfCubicalComplex(M,true);
return
Sum(List([0..Dimension(M)],i->((-1)^i)*C!.dimension(i)));
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(EulerCharacteristic,
"Euler characteristic  of a chain complex",
[IsHapChainComplex],
function(C);

return
Sum(List([0..Length(C)],i->((-1)^i)*C!.dimension(i)));
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Length,
"Length of a chain complex",
[IsHapChainComplex],
function(C);

return
EvaluateProperty(C,"length");
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Length,
"Length of a cochain complex",
[IsHapCochainComplex],
function(C);

return
EvaluateProperty(C,"length");
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallOtherMethod(Length,
"Length of a sparse chain complex",
[IsHapSparseChainComplex],
function(C);

return
EvaluateProperty(C,"length");
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallOtherMethod(EulerCharacteristic,
"Euler characteristic  of a sparse chain complex",
[IsHapSparseChainComplex],
function(C);

return
Sum(List([0..Length(C)],i->((-1)^i)*C!.dimension(i)));
end);
#####################################################################
#####################################################################





#################################################################
#################################################################
InstallGlobalFunction(ReadMatrixAsPureCubicalComplex,
function(AA)
local PCC,B, fn, Dims,x, Dim, i, Iter,A;

A:=SSortedList(AA);

Dims:=[];
Dim:=Length(A[1]);
for i in [1..Dim] do
Add(Dims,3+Maximum(List(A,x->x[i])));
od;

B:=List([1..1+Dims[Length(Dims)]],a->0);
for i in Reversed([1..Dim-1]) do
B:=List([1..Dims[i]],a->StructuralCopy(B));
od;

##############
if Dim=1 then
for x in A do
B[1+x[1]]:=1;
od;
fi;
##############
if Dim=2 then
for x in A do
B[1+x[1]][1+x[2]]:=1;
od;
fi;
##############
if Dim=3 then
for x in A do
B[1+x[1]][1+x[2]][1+x[3]]:=1;
od;
fi;
##############
if Dim=4 then
for x in A do
B[1+x[1]][1+x[2]][1+x[3]][1+x[4]]:=1;
od;
fi;
##############
if Dim=5 then
for x in A do
B[1+x[1]][1+x[2]][1+x[3]][1+x[4]][1+x[5]]:=1;
od;
fi;
##############
if Dim=6 then
for x in A do
B[1+x[1]][1+x[2]][1+x[3]][1+x[4]][1+x[5]][1+x[6]]:=1;
od;
fi;
##############
if Dim=7 then
for x in A do
B[1+x[1]][1+x[2]][1+x[3]][1+x[4]][1+x[5]][1+x[6]][1+x[7]]:=1;
od;
fi;
##############
if Dim=8 then
for x in A do
B[1+x[1]][1+x[2]][1+x[3]][1+x[4]][1+x[5]][1+x[6]][1+x[7]][1+x[8]]:=1;
od;
fi;
##############
if Dim>8 then
Print("This function is not yet implemented for matrices with more than 8 columns.\n");
return fail;
fi;

return PureCubicalComplex(B);;
end);
#####################################################################
#####################################################################



HAPAAA:=0;
#################################################################
#################################################################
InstallGlobalFunction(ReadImageAsPureCubicalComplex,
function(file,threshold)
# Inputs a string "file" that points either to a single image file, or 
# to a list of suitable image files, and returns either a 2-dimensional 
# or 3-dimensional cubical complex.
local f,i,x,prog,B,A,AA,pth,file1, file2, tmpdir;

pth:=HAP_ROOT;
prog:=Concatenation(pth,"PolyComplexes/prog");


##################################
if IsString(file) then

tmpdir := DirectoryTemporary();;
file1:=Filename( tmpdir , "im.txt" );
file2:=Filename( tmpdir , "im.g" );

i:=Concatenation("convert -colorspace RGB -depth 8  ",file," ",file1);
Exec(i);
i:=Concatenation("perl ",prog," ",file1," > ",file2);
Exec(i);

Exec(Concatenation("rm ",file1));
Read(file2);

B:=StructuralCopy(HAPAAA);
HAPAAA:=0;

A:=[];
for i in [1..B[1]] do
A[i]:=List([1..B[2]],a->0);
od;

for x in B{[3..Length(B)-1]} do
A[x[1]][x[2]]:=x[3];
od;
if not IsInt(threshold) then return A; fi;
return ArrayToPureCubicalComplex(A,threshold);
fi;
#################################  

##################################
if IsList(file) then
AA:=[];
for f in file do
tmpdir := DirectoryTemporary();;
file1:=Filename( tmpdir , "im.txt" );
file2:=Filename( tmpdir , "im.g" );
#i:=Concatenation("convert ",f," /tmp/im.txt");
i:=Concatenation("convert ",f," ",file1);
Exec(i);
#i:=Concatenation("perl ",prog," /tmp/im.txt >/tmp/im.g");
i:=Concatenation("perl ",prog," ",file1," > ",file2 );
Exec(i);

#Read("/tmp/im.g");
#Exec("rm /tmp/im.g");
#Exec("rm /tmp/im.txt");
Read(file2);
B:=StructuralCopy(HAPAAA);
A:=[];
for i in [1..B[1]] do
A[i]:=List([1..B[2]],a->0);
od;
for x in B{[3..Length(B)-1]} do
A[x[1]][x[2]]:=x[3];
od;
Add(AA,StructuralCopy(A));
HAPAAA:=0;
od;

return ArrayToPureCubicalComplex(AA,threshold);
fi;
#################################

end);
#################################################################
#################################################################


#################################################################
#################################################################
InstallGlobalFunction(ReadImageSequenceAsPureCubicalComplex,
function(directory,threshold)
local L,file,linecommand,instr,f;

file:="/tmp/imagesfile.txt";
L:=[];

linecommand:=Concatenation("ls ",directory," > ",file);
Exec(linecommand);

instr:=InputTextFile(file);

f:=ReadLine(instr);
while not f=fail do
Add(L,Concatenation(directory,"/",f{[1..Length(f)-1]}));
f:=ReadLine(instr);
od;

return 
ReadImageAsPureCubicalComplex(L,threshold);
end);
#################################################################
#################################################################


#################################################################
#################################################################
InstallGlobalFunction(PureCubicalComplexToCubicalComplex,
function(M)
# Converts a pure cubical complex to a cubical complex
local A,x,i,dim,dim1,dims,ball,b,dimsM,ArrayValueDim, ArrayValueDim1,
ArrayAssignDim,ArrayIt,dimsSet,Fun;

if IsHapFilteredPureCubicalComplex(M) then return
FilteredPureCubicalComplexToCubicalComplex(M); fi;
  
dim:=Dimension(M);
dim1:=dim-1;
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim-1);
ArrayAssignDim:=ArrayAssignFunctions(dim);
ArrayIt:=ArrayIterate(dim);
dimsM:=EvaluateProperty(M,"arraySize");
dims:=List(EvaluateProperty(M,"arraySize"),n->2*n+1);
A:=List([1..dims[1]],a->0);
for i in [2..Dimension(M)] do
A:=List([1..dims[i]],a->StructuralCopy(A));
od;

ball:=Cartesian(List([1..dim],i->[-1,0,1]));
dimsSet:=List([1..dim],x->[1..dimsM[x]]);

#########################
Fun:=function(i);
if ArrayValueDim(M!.binaryArray,i)=1 then
for b in ball do
x:=2*i+b;
ArrayAssignDim(A,x,1);
od;
fi;
end;
#########################
ArrayIt(dimsSet,Fun);

return Objectify(HapCubicalComplex,
           rec(
           binaryArray:=A,
           properties:=[
           ["dimension",dim],
           ["arraySize",dims]]
           ));

end);
#################################################################
#################################################################


#################################################################
#################################################################
InstallGlobalFunction(ChainComplexOfCubicalComplex,
function(arg)
local M,faces,dims,dim,dim1,x,y,z,numevens,Dimsn,Boundary,BinLst,LstBin,
ArrayValueDim,ArrayValueDim1,ArrayAssignDim,ArrayIt,dimsSet,Fun,LF;
# Inputs a cubical complex and returns the associated cellular 
# chain complex.

M:=arg[1];

if not IsHapCubicalComplex(M) then
Print("This function must be applied to a cubical complex.\n");
return fail;
fi;


faces:=List([0..Dimension(M)+10],i->0);
dims:=EvaluateProperty(M,"arraySize");
dimsSet:=List(dims,d->[1..d]);
dim:=Dimension(M);
dim1:=dim-1;
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim1);
ArrayAssignDim:=ArrayAssignFunctions(dim);
ArrayIt:=ArrayIterate(dim);


BinLst:=M!.binaryArray*0; #BinLst will be an array where the value 
                         #of an entry is the position of the cell 
                         #corresponding to that entry in the list 
                         #of all cells of similar dimension.
LstBin:=List([0..Dimension(M)],i->[]);	         
                         #LstBin is a list of lists. The i-th entry of the
                         #j-th list is the coordinates in M!.binaryArray
                         #of the i-th cell in dimension j-1.


##############################
numevens:=function(x);
return Length(Filtered(x,i->IsEvenInt(i)));
end;
##############################

if Length(arg)=1 then

#######################
Fun:=function(x) local y,z;
if ArrayValueDim(M!.binaryArray,x)=1 then
y:=numevens(x)+1;
faces[y]:=faces[y]+1;
ArrayAssignDim(BinLst,x,faces[y]);
LstBin[y][faces[y]]:=x;
fi;
end;
#######################
ArrayIt(dimsSet,Fun);

else

#######################
Fun:=function(x) local y;
if ArrayValueDim(M!.binaryArray,x)=1 then
y:=numevens(x)+1;
faces[y]:=faces[y]+1;
fi;
end;
#######################
ArrayIt(dimsSet,Fun);

fi;

#######################################
Dimsn:=function(n);
return faces[n+1];
end;
#######################################

LF:=[];
for x in [1..Length(faces)] do
LF[x]:=0*[1..faces[x]];
od;

if EvaluateProperty(M,"nonregular")=true then
#######################################
Boundary:=function(n,j) #THIS ONLY WORKS CORRECTLY MOD 2 AT PRESENT
local x,poscells,negcells,nn,a,b,cnt;
poscells:=[];
negcells:=[];

cnt:=0;
nn:=Reversed(LstBin[n+1][j]);
for x in [1..Length(nn)] do
  if IsEvenInt(nn[x]) then
    cnt:=cnt+1;
    a:=StructuralCopy(nn);
    a[x]:=a[x]+1;
    b:=StructuralCopy(nn);
    b[x]:=b[x]-1;
    if IsOddInt(cnt) then
        Append(poscells,M!.rewrite([a])); 
        Append(negcells,M!.rewrite([b]));
    else
        Append(poscells,M!.rewrite([b]));
        Append(negcells,M!.rewrite([a]));
    fi;
  fi;
od;

Apply(poscells,x->ArrayValueDim(BinLst,Reversed(x)));
Apply(negcells,x->ArrayValueDim(BinLst,Reversed(x)));

#nn:=List([1..faces[n]],i->0);
nn:=0*LF[n];

for x in poscells do
  nn[x]:=nn[x]-1;
od;
for x in negcells do
  nn[x]:=nn[x]+1;
od;

return nn;

end;
#######################################
else
#######################################
Boundary:=function(n,j)
local x,poscells,negcells,nn,a,b,cnt;
poscells:=[];
negcells:=[];

cnt:=0;
nn:=LstBin[n+1][j];
for x in [1..Length(nn)] do
  if IsEvenInt(nn[x]) then
    cnt:=cnt+1;
    a:=StructuralCopy(nn);
    a[x]:=a[x]+1;
    b:=StructuralCopy(nn);
    b[x]:=b[x]-1;
    if IsOddInt(cnt) then
        Add(poscells,a); 
        Add(negcells,b);
    else
        Add(poscells,b);
        Add(negcells,a);
    fi;
  fi;
od;

Apply(poscells,x->ArrayValueDim(BinLst,x));
Apply(negcells,x->ArrayValueDim(BinLst,x));
nn:=List([1..faces[n]],i->0);
#nn:=(LF[n]);
Apply(nn,x->0);

for x in poscells do
  nn[x]:=-1;
od;
for x in negcells do
  nn[x]:=1;
od;

return nn;
end;
#######################################
fi;

return
Objectify(HapChainComplex,
           rec(
           dimension:=Dimsn,
	   boundary:=Boundary,
	   coordinateToPosition:=BinLst,
	   positionToCoordinate:=LstBin,
           properties:=[
	   ["length",EvaluateProperty(M,"dimension")],
	   ["type","chainComplex"],
           ["characteristic",0]]
           ));

end);
#################################################################
#################################################################


#####################################################################
#####################################################################
InstallMethod(ChainComplex,
"Cellular chain complex of a cubical complex",
[IsHapCubicalComplex],
function(M) 
return ChainComplexOfCubicalComplex(M);;
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallMethod(ChainComplexOfPair,
"Cellular chain complex of a pair of cubical complexes",
[IsHapCubicalComplex,IsHapCubicalComplex],
function(M,S)
return ChainComplexOfCubicalPair(M,S);;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(ChainComplexOfPair,
"Cellular chain complex of a pair of simplicial complexes",
[IsHapSimplicialComplex,IsHapSimplicialComplex],
function(M,S)
return ChainComplexOfSimplicialPair(M,S);;
end);
#####################################################################
#####################################################################



#####################################################################
#####################################################################
InstallOtherMethod(ChainComplexOfPair,
"Cellular chain complex of a pair of pure cubical complexes",
[IsHapPureCubicalComplex,IsHapPureCubicalComplex],
function(MM,SS)
local M,S,P;
P:=ExcisedPureCubicalPair(MM,SS);
M:=P[1];S:=P[2];
#M:=MM;S:=SS;
return ChainComplexOfCubicalPair(
PureCubicalComplexToCubicalComplex(M),
PureCubicalComplexToCubicalComplex(S));;
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallOtherMethod(ChainComplex,
"Cellular chain complex of a pure cubical complex",
[IsHapPureCubicalComplex],
function(M)
return ChainComplexOfCubicalComplex(PureCubicalComplexToCubicalComplex(M));;
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(ChainComplexOfCubicalPair,
function(M,S)
local faces,dim,dim1,dims,index,x,y,z,numevens,Dimsn,Boundary,BinLst,LstBin,
ArrayValueDim,ArrayValueDim1,Generator2Coordinates,Coordinates2Generator,NN;
# Inputs a pair of cubical complexes and returns the  cellular
# chain complex of " M/S ".

faces:=List([0..Dimension(M)],i->0);
dims:=EvaluateProperty(M,"arraySize");
index:=Cartesian(List([1..Dimension(M)],x->[1..dims[x]]));
##
##
##For n=2,3,4 I should write quicker code that does not construct index
##
##
##
dim:=Dimension(M);
dim1:=dim-1;
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim1);

BinLst:=M!.binaryArray*0; #BinLst will be an array where the value
                         #of an entry is the position of the cell
                         #corresponding to that entry in the list
                         #of all cells of similar dimension.
LstBin:=List([0..Dimension(M)],i->[]);
                         #LstBin is a list of lists. The i-th entry of the
                         #j-th list is the coordinates in M!.binaryArray
                         #of the i-th cell in dimension j-1.


##############################
numevens:=function(x);
return Length(Filtered(x,i->IsEvenInt(i)));
end;
##############################

for x in index do
if 
	ArrayValueDim(M!.binaryArray,x)=1 
	and ArrayValueDim(S!.binaryArray,x)=0
then
y:=numevens(x)+1;
faces[y]:=faces[y]+1;
z:=ArrayValueDim1(BinLst,x{[2..Length(x)]});
z[x[1]]:=faces[y];
LstBin[y][faces[y]]:=x;
fi;
od;

#######################################
Dimsn:=function(n);
return faces[n+1];
end;
#######################################

NN:=List([1..dim],i-> List([1..faces[i]],a->0));
#######################################
Boundary:=function(n,j)
local x,poscells,negcells,nn,a,b,cnt;
poscells:=[];
negcells:=[];

cnt:=0;
nn:=LstBin[n+1][j];
for x in [1..Length(nn)] do
if IsEvenInt(nn[x]) then
cnt:=cnt+1;
a:=StructuralCopy(nn);
a[x]:=a[x]+1;
b:=StructuralCopy(nn);
b[x]:=b[x]-1;
if IsOddInt(cnt) then
if ArrayValueDim(S!.binaryArray,a)=0 then
Add(poscells,a);
fi;
if ArrayValueDim(S!.binaryArray,b)=0 then
Add(negcells,b);
fi;
else
if ArrayValueDim(S!.binaryArray,b)=0 then
Add(poscells,b);
fi;
if ArrayValueDim(S!.binaryArray,a)=0 then
Add(negcells,a);
fi;
fi;
fi;
od;

Apply(poscells,x->ArrayValueDim(BinLst,x));
Apply(negcells,x->ArrayValueDim(BinLst,x));

#nn:=List([1..faces[n]],i->0);
nn:=StructuralCopy(NN[n]);
for x in poscells do
nn[x]:=-1;
od;
for x in negcells do
nn[x]:=1;
od;

return nn;
end;
#######################################

#######################################
Generator2Coordinates:=function(n,j);
return LstBin[n+1][j];
end;
#######################################

#######################################
Coordinates2Generator:=function(x);
return ArrayValueDim(BinLst,x);
end;
#######################################

return
Objectify(HapChainComplex,
           rec(
           dimension:=Dimsn,
           boundary:=Boundary,
	   generator2Coordinates:=Generator2Coordinates,
	   coordinates2Generator:=Coordinates2Generator, 
           properties:=[
           ["length",EvaluateProperty(M,"dimension")],
           ["type","chainComplex"],
           ["characteristic",0]]
           ));

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(ChainMapOfCubicalPairs,
function(arg)
			#X is a subspace of Y.
			#S is a contractible subspace of X (and T). 
			#T is a contractible subspace of Y.
			#We return C(X/S)-->C(Y/T).
local
	X,S,Y,T,C,D,map;

###########INPUT####################
if Length(arg)=4 then
X:=arg[1];
S:=arg[2];
Y:=arg[3];
T:=arg[4];
C:=ChainComplexOfPair(X,S);
D:=ChainComplexOfPair(Y,T);
fi;

if Length(arg)=2 then
C:=arg[1];
D:=arg[2];
fi;
####################################

###############################################
map:=function(v,n)
local x,w,i;

w:=List([1..D!.dimension(n)],a->0);

for i in [1..Length(v)] do
if not v[i]=0 then
x:=C!.generator2Coordinates(n,i);
x:=D!.coordinates2Generator(x);
if not x=0 then w[x]:=v[i]; fi;
fi;
od;

return w;
end;
###############################################

return Objectify(HapChainMap,
        rec(
           source:=C,
           target:=D,
           mapping:=map,
           properties:=[ ["type","chainMap"],
           ["characteristic", 0]
           ]));

end);
#####################################################################
##################################################################### 


######################################################################
######################################################################
InstallGlobalFunction(ContractibleSubcomplexOfPureCubicalComplex,
function(T)
local A,P,mx,Pmx;

if not IsHapPureCubicalComplex(T) then
Print("This function can only be applied to pure cubical complexes.\n");
return fail;
fi;

if Size(T)=0 then return T; fi;

P:=PathComponentOfPureCubicalComplex(T,0);
P:=List([1..P],x->Size(PathComponentOfPureCubicalComplex(T,x)));
mx:=Position(P,Maximum(P));
Pmx:=PathComponentOfPureCubicalComplex(T,mx);
A:=ContractibleSubArray(Pmx!.binaryArray);


return Objectify(HapPureCubicalComplex,
	rec(
	binaryArray:=A,
	properties:=StructuralCopy(T!.properties)
	));

end);
######################################################################
######################################################################


######################################################################
######################################################################
InstallGlobalFunction(AcyclicSubcomplexOfPureCubicalComplex,
function(T)
local A,P,x;

if not IsHapPureCubicalComplex(T) then
Print("This function can only be applied to pure cubical complexes.\n");
return fail;
fi;

P:=PathComponentOfPureCubicalComplex(T,0);
A:=0*T!.binaryArray;

for x in T!.pathReps do
ArrayAssign(A,x,1);
od;

A:=PureCubicalComplex(A);;
A:=HomotopyEquivalentMaximalPureCubicalSubcomplex(T,A);

return A;

end);
######################################################################
######################################################################


#################################################################
#################################################################
InstallGlobalFunction(WritePureCubicalComplexAsImage,
function(T,file,ext)
local
        AppendTo, PrintTo, A,i,j,rows,cols,colour,filetxt;
# This writes a 2-dimensional cubical complex to an image file.

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

if not Dimension(T)=2 then
Print("There is no method for viewing a pure cubical complex of dimension ",
Dimension(T),".\n"); return fail; fi;

A:=T!.binaryArray;
rows:=Length(A);;
cols:=Length(A[1]);

filetxt:=Concatenation(file,".txt");

PrintTo(filetxt,"# ImageMagick pixel enumeration: ",
Length(A),",",Length(A[1]),",255,RGB\n");

for i in [1..rows] do
for j in [1..cols] do
if A[i][j]=0 then colour:="(255,255,255)"; else colour:="(0,0,100)";fi;
AppendTo(filetxt,i-1,",",j-1,": ",colour,"\n");
od;
od;

i:=Concatenation("convert ",filetxt," ",file,".",ext);
Exec(i);
i:=Concatenation("rm ",filetxt);
Exec(i);

end);
#################################################################
#################################################################


#################################################################
#################################################################
InstallGlobalFunction(ViewPureCubicalComplex,
function(arg)
local i,A,viewer,T,file,file1,file2,tmpdir;
# Function for viewing 2-dimensional pure cubical complexes.

T:=PureCubicalComplex(1*arg[1]!.binaryArray);

if not (Dimension(T)=2 or Dimension(T)=3) then
Print("There is no method for viewing a pure cubical complex of dimension ",
Dimension(T),".\n"); return fail; fi;

if Dimension(T)=3 then

View3dPureComplex(T);

else

T!.binaryArray:=FrameArray(T!.binaryArray);
T!.binaryArray:=FrameArray(T!.binaryArray);

if Length(arg)>1 then viewer:=arg[2];
else viewer:=DISPLAY_PATH;fi;
tmpdir := DirectoryTemporary();;
file:=Filename( tmpdir , "HAPtmpImage.png" );
file1:=file{[1..Length(file)-4]};
file2:=file{[Length(file)-2..Length(file)]};
#WritePureCubicalComplexAsImage(T,"/tmp/HAPtmpImage","png");
WritePureCubicalComplexAsImage(T,file1,file2);
#Exec(Concatenation(viewer," ","/tmp/HAPtmpImage.png"));
Exec(Concatenation(viewer," ",file));
Sleep(2);
#Exec(Concatenation("rm  ","/tmp/HAPtmpImage.png"));

fi;
end);
#################################################################
#################################################################


#################################################################
#################################################################
InstallGlobalFunction(CechComplexOfPureCubicalComplex,
function(arg)
local
	M,DIM,Vertices,
	Simplices,
	SimplicesLst,
	NrSimplices,
	EnumeratedSimplex,
 	S,
	VertexCoordinates,
        dim,dim1,dims, 
	CartProd,
	Ball,
	Star,
	EdgeStar,
	ArrayValueDim,
	ArrayValueDim1,
	v, x, y, U, V, i;

M:=arg[1];
if Length(arg)>1 then DIM:=1+arg[2]; else DIM:=10000; fi; ##SLOPPY!!
if not IsHapPureCubicalComplex(M) then
Print("This function can only be applied to a pure cubical complex.\n");
return fail; fi;

dim:=Dimension(M);
dim1:=dim-1;
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim1);
dims:=EvaluateProperty(M,"arraySize");
S:=StructuralCopy(M!.binaryArray);
Vertices:=0;
VertexCoordinates:=[];
Ball:=Cartesian(List([1..dim],i->[-1,1]));
CartProd:=Cartesian(List([1..dim],a->[1..dims[a]]));

for x in CartProd do
if ArrayValueDim(M!.binaryArray,x)=1 then Vertices:=Vertices+1;
y:=ArrayValueDim1(S,x{[2..dim]});
y[x[1]]:=Vertices;
VertexCoordinates[Vertices]:=x;
fi;
od;

Vertices:=[1..Vertices];
SimplicesLst:=List([1..1000],i->[]);  #VERY SLOPPY!!!
SimplicesLst[1]:=List(Vertices,i->[i]);

#################################################################
NrSimplices:=function(n);
return Length(SimplicesLst[n+1]);
end;
#################################################################

#################################################################
Simplices:=function(n,i);
return SimplicesLst[n+1][i];
end;
#################################################################

#################################################################
Star:=function(x)
# For an odd  array point x in the corresponding cubical complex T 
# this function returns the list of all integers in Vertices
# corresponding to those pure cells in T that touch the vertex
# T!.binaryArray[x].
local y,z,verts;

verts:=[];

for y in Ball do
z:=(x+y)/2; #z is even
if (not 0 in z) and (not -1 in dims - z) then
     if ArrayValueDim(M!.binaryArray,z)=1 then
     Add(verts, ArrayValueDim(S,z)); fi;
fi;
od;

return verts;
end;
#################################################################

for v in Cartesian(List([1..dim],a->List([0..dims[a]],j->2*j+1))) do
U:=Star(v);
for i in [2..Minimum(DIM,Length(U))] do
V:=Combinations(U,i);

for y in V do
Add(SimplicesLst[i],y);
od;
od;
od;

for i in [1..Length(SimplicesLst)] do
SimplicesLst[i]:=SSortedList(SimplicesLst[i]);
od;

#############################################
EnumeratedSimplex:=function(v);
return
Position(SimplicesLst[Length(v)],v);
end;
#############################################

return
Objectify(HapSimplicialComplex,
           rec(
           vertices:=Vertices,
           simplices:=Simplices,
           simplicesLst:=SimplicesLst,
	   nrSimplices:=NrSimplices,
	   enumeratedSimplex:=EnumeratedSimplex,
           properties:=[
           ["dimension",PositionProperty(SimplicesLst,IsEmpty)-2]
	   ]
           ));

end);
#################################################################
#################################################################



#################################################################
#################################################################
InstallGlobalFunction(HomologyOfPureCubicalComplex,
function(arg)
local M,S,T,C,H,dim,b,P,i,n,p;  #Only use this function on path-connected spaces

M:=arg[1];

######################
if Length(arg)>1 then
if arg[2]=0 then
return [1..PathComponentOfPureCubicalComplex(M,0)]*0;
fi;
fi;
######################

dim:=Dimension(M);
T:=rec();
T.binaryArray:=StructuralCopy(M!.binaryArray);
T.properties:=StructuralCopy(M!.properties);
T:=Objectify(HapPureCubicalComplex,T);
#T:=ZigZagContractedPureCubicalComplex(T);
#T:=CropPureCubicalComplex(T);

#S:=ContractibleSubcomplexOfPureCubicalComplex(T);
S:=AcyclicSubcomplexOfPureCubicalComplex(T);
C:=ChainComplexOfPair(T,S);

#####################
if Length(arg)=1 then
H:=[];
H[1]:=List([1..PathComponentOfPureCubicalComplex(T,0)],a->0);
for i in [1..dim-1] do
H[i+1]:=Homology(C,i);
od;
H[dim+1]:=[];
return H;
fi;
#####################

#####################
if Length(arg)=2 then 
n:=arg[2];
return Homology(C,n);
fi;
#####################

#####################
if Length(arg)=3 then
n:=arg[2];
p:=arg[3];
return Homology(TensorWithIntegersModP(C,p),n);
fi;
#####################

end);
#################################################################
#################################################################


#################################################################
#################################################################
InstallOtherMethod(Homology,
"Homology of a pure cubical complex",
[IsHapPureCubicalComplex],
function(MM)
local M,i,k,b,P,H,tmp,dim;

M:=ZigZagContractedPureCubicalComplex(MM);
M:=CropPureCubicalComplex(M);
############################
if Dimension(M)=2 then
H:=List(BettinumbersOfPureCubicalComplex_dim_2(M),b->[1..b]*0);
return H;
fi;
############################

dim:=Dimension(M);
H:=List([0..dim],x->[]);
b:=PathComponentOfPureCubicalComplex(M,0);
if b=0 then return H; fi;

for i in [1..b] do
P:=PathComponentOfPureCubicalComplex(M,i);
tmp:=HomologyOfPureCubicalComplex(P);
for k in [0..dim] do
Append(H[k+1],tmp[k+1]);
od;
od;

return H;
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallOtherMethod(Homology,
"Homology of a pure permutahedral complex",
[IsHapPurePermutahedralComplex, IsInt],
function(MM,n)
local M;

if Dimension(MM)<=3 then M:=ZigZagContractedComplex(MM);
else M:=MM; fi;

M:=RegularCWComplex(M);
return Homology(M,n);
end);
#################################################################
#################################################################




#################################################################
#################################################################
InstallOtherMethod(Bettinumbers,
"Ranks of homologies of a pure cubical complex",
[IsHapPureCubicalComplex],
function(M)
local H,T;

############################
if Dimension(M)=2 then
return BettinumbersOfPureCubicalComplex_dim_2(M);
fi;
############################

H:=Homology(M);
H:=List(H,x->Length(Filtered(x,a->a=0)));
return H;
end);
#################################################################
#################################################################


#################################################################
#################################################################
InstallOtherMethod(Bettinumbers,
"Rank of n-th homology of a pure cubical complex",
[IsHapPureCubicalComplex,IsInt],
function(M,n)
local H;
if n=0 then return PathComponentOfPureCubicalComplex(M,0); fi;
if Dimension(M)=2 then return Bettinumbers(M)[n+1]; fi;
H:=Homology(M,n);
H:=Length(Filtered(H,a->a=0));
return H;
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallOtherMethod(Bettinumbers,
"Rank of n-th homology of a chain complex",
[IsHapChainComplex,IsInt],
function(M,n)
local H, p;
if EvaluateProperty(M,"characteristic")<=0  then
return Homology(TensorWithRationals(M),n);
fi;
p:=EvaluateProperty(M,"characteristic");
return Homology(TensorWithIntegersModP(M,p),n);
end);
#################################################################
#################################################################





#################################################################
#################################################################
InstallOtherMethod(Homology,
"Homology of a pure cubical complex",
[IsHapPureCubicalComplex,IsInt,IsInt],
function(MM,n,p) 
local M,i,k,b,P,H,tmp,dim;


############################
if Dimension(MM)=2 then
if p=0 then
return [1..BettinumbersOfPureCubicalComplex_dim_2(MM)[n+1]]*0;
fi;
return BettinumbersOfPureCubicalComplex_dim_2(MM)[n+1]; 
fi;
############################
M:=ZigZagContractedPureCubicalComplex(MM);
M:=CropPureCubicalComplex(M);
dim:=Dimension(M);
b:=PathComponentOfPureCubicalComplex(M,0);
if b=0 then return []; fi;

if p=0 then
H:=[];
for i in [1..b] do
P:=PathComponentOfPureCubicalComplex(M,i);
tmp:=HomologyOfPureCubicalComplex(P,n,p);
Append(H,tmp);
od;
return H;
fi;

H:=0;
for i in [1..b] do
P:=PathComponentOfPureCubicalComplex(M,i);
tmp:=HomologyOfPureCubicalComplex(P,n,p);
H:=H+tmp;
od;
return H;




end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallOtherMethod(Homology,
"Homology of a pure cubical complex",
[IsHapPureCubicalComplex,IsInt],
function(M,n)
return Homology(M,n,0);
end);
#################################################################
#################################################################

#################################################################
################################################################
InstallGlobalFunction(SingularitiesOfPureCubicalComplex,
function(MM,radius,tolerance)
local
        M,B,
        CART, dim,dim1,dims,
        IsSingular,
	Circle, Ball, dimsBall,
	ArrayValueDim,
	ArrayValueDim1,
        x,z;

#############################################
if not IsHapPureCubicalComplex(MM) then
Print("This function must be applied to a pure cubical complex.\n");
return fail;
fi;
#############################################

M:=BoundaryOfPureCubicalComplex(MM);
#M:=ContractPureCubicalComplex(M);
B:=M!.binaryArray*0;

dim:=Dimension(M);
dim1:=dim-1;
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim1);
dims:=EvaluateProperty(M,"arraySize");
CART:=Cartesian(List([1..dim],a->[1..dims[a]]));

Circle:=Cartesian(List([1..dim],i->[-radius..radius]));
Ball:=Cartesian(List([1..dim],i->[-radius+1..radius-1]));
Circle:=Filtered(Circle,x->not x in Ball);

Ball:=[1..2*radius+1];
for x in [2..dim] do
Ball:=List(Ball,x->StructuralCopy(Ball));
od;
Ball:=Ball*0;

for x in Circle do
z:=ArrayValueDim1(Ball,x{[2..Length(x)]}+radius+1);
z[x[1]+radius+1]:=1;
od;


Ball:=rec(
	  binaryArray:=Ball,
	  properties:=[
	  ["dimension",dim],
	  ["arraySize", List([1..dim],a->2*radius+1)]]);
Ball:=Objectify(HapPureCubicalComplex,Ball);


dimsBall:=EvaluateProperty(Ball,"arraySize");


########################
IsSingular:=function(x)
local w,y,yy,z,B,ball,L,w1,y1;
if ArrayValueDim(M!.binaryArray,x)=0 then return false;fi;
B:=StructuralCopy(Ball!.binaryArray);
ball:=ArrayToPureCubicalComplex(B,1);
ball!.binaryArray:=B;

for y in Circle do
z:=x+y;yy:=y+radius+1;
if (not  Minimum(z)<1) and (not  Minimum(dims - z)<0) then
  if ArrayValueDim(M!.binaryArray,z)=1 then
  w:=ArrayValueDim1(ball!.binaryArray,yy{[2..Length(yy)]});
  w[yy[1]]:=0;
  fi;
fi;

od;


yy:=PathComponentOfPureCubicalComplex(ball,0);  
L:=List([1..yy],a->PathComponentOfPureCubicalComplex(ball,a));
L:=List(L,a->Flat(a!.binaryArray));
L:=List(L, a->Filtered(a,i->i=1));
L:=List(L,a->Length(a));
w:=Maximum(L);
w1:=Position(L,w);
L[w1]:=0;
y:=Maximum(L);

if 100*AbsInt(y-w)/Length(Circle)>tolerance then return true; else return false; fi;

end;
########################

for x in CART do
if IsSingular(x) then
z:=ArrayValueDim1(B,x{[2..dim]});
z[x[1]]:=1;
fi;
od;

return Objectify(HapPureCubicalComplex,
           rec(
           binaryArray:=B,
           properties:=[
           ["dimension",Dimension(M)],
           ["arraySize",dims]]
           ));

end);
#################################################################
#################################################################


#################################################################
#################################################################
InstallGlobalFunction(ChainInclusionOfPureCubicalPair,
function(M,N)
local
	map,C,D,k,U;

C:=ChainComplex(M);
D:=ChainComplex(N);

U:=[];
for k in [0..EvaluateProperty(D,"length")] do
U[k+1]:=List([1..D!.dimension(k)],i->0);
od;

###########################
map:=function(v,k)
local u,i;

u:=U[k+1]*0;

for i in [1..C!.dimension(k)] do
if not v[i]=0 then
u[ ArrayValue(D!.coordinateToPosition,C!.positionToCoordinate[k+1][i]) ]:=v[i];
fi;
od;

return u;
end;
###########################


return Objectify(HapChainMap,
        rec(
           source:=C,
           target:=D,
           mapping:=map,
           properties:=[ ["type","chainMap"],
           ["characteristic", 0],
           ]));

end);
#################################################################
#################################################################


#################################################################
#################################################################
InstallGlobalFunction(DirectProductOfPureCubicalComplexes,
function(M,N)
local
	D,
	dimM,dimN,dim,
	dimsM,dimsN,dims,
	CART,
ArrayValueDimM,
ArrayValueDimN,
ArrayValueDim1,
	x,w;

dimM:=Dimension(M);
dimN:=Dimension(N);
dim:=dimM+dimN;
ArrayValueDimM:=ArrayValueFunctions(dimM);
ArrayValueDimN:=ArrayValueFunctions(dimN);
ArrayValueDim1:=ArrayValueFunctions(dim-1);

dimsM:=EvaluateProperty(M,"arraySize");
dimsN:=EvaluateProperty(N,"arraySize");
dims:=Concatenation(dimsM,dimsN);
CART:=Cartesian(List([1..dim],i->[1..dims[i]]));
D:=[1..dims[1]]*0;
for x in [2..dim] do
D:=List([1..dims[x]],i->StructuralCopy(D));
od;


for x in CART do
if ArrayValueDimM(M!.binaryArray,x{[1..dimM]})=1 and ArrayValueDimN(N!.binaryArray,x{[dimM+1..dimM+dimN]})=1
then
w:=ArrayValueDim1(D,x{[2..dim]});
w[x[1]]:=1;
fi;
od;

D:=PureCubicalComplex(D);


return D;
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(PureCubicalComplexToTextFile,
function(file,M)
local
	dim,dims,CART,x,i, AppendTo, PrintTo;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

dim:=Dimension(M);
dims:=EvaluateProperty(M,"arraySize");
CART:=Cartesian(List(dims,d->[1..d]));

for x in CART do
if ArrayValue(M!.binaryArray,x)=1 then 
for i in [1..dim] do
AppendTo(file,String(x[i])," ");
od;
AppendTo(file,"\n");
fi;
od;

end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(ExcisedPureCubicalPair,
function(M,S)
local B,intS,N, M1, S1;

if Dimension(M)=2 then return ExcisedPureCubicalPair_dim_2(M,S); fi;

#B:=BoundaryOfPureCubicalComplex(S);
#intS:=PureCubicalComplexDifference(S,B);
#N:=PureCubicalComplexDifference(M,intS);

B:=PureComplexDifference(M,S);
B:=ThickenedPureComplex(B);
B:=PureCubicalComplexDifference(S,B);
S1:=PureComplexDifference(S,B);
M1:=PureComplexDifference(M,B);

#return [N,B];
return [M1,S1];
end);
#################################################################
#################################################################



###################################################
###################################################
InstallGlobalFunction(MorseFiltration,
function(arg)
local gradient, M,N,d,bool,dim, dims, G, B, t, J;

#################################
gradient:=function(arg)
local A,N,J,bool,B,dim,dims,cart, x, ArrayAss;

A:=arg[1];
N:=arg[2];
J:=arg[3];
if Length(arg)>3 then bool:=arg[4]; else bool:=true; fi;

B:=StructuralCopy(A);
dim:=ArrayDimension(A);
dims:=ArrayDimensions(A);

cart:=List(dims,x->[1..x]);
cart:=Cartesian(cart);

ArrayAss:=ArrayAssignFunctions(dim);

if bool then
    for x in cart do
    if x[N]<=J then
    ArrayAss(B,x,0);
    fi;
    od;
else
    for x in cart do
    if x[N]>=J then
    ArrayAss(B,x,0);
    fi;
    od;
fi;


return B;
end;
#################################


M:=arg[1];
N:=arg[2];
d:=arg[3];
if Length(arg)>3 then bool:=arg[4]; else bool:=true; fi;

dim:=ArrayDimension(M!.binaryArray);
dims:=ArrayDimensions(M!.binaryArray);
G:=[];
t:=Int(dims[N]/d);
J:=0;

while J<dims[N] do
B:=StructuralCopy(M!.binaryArray);
B:=gradient(B,N,J,bool);
B:=PureCubicalComplex(B);
Add(G,B);
J:=J+t;
od;

G:=Filtered(G,x->Size(x)>0);
if bool then
return Reversed(G);
else
return G;
fi;

end);
##########################################
##########################################

##########################################
##########################################
InstallGlobalFunction(SkeletonOfCubicalComplex,
function(MM,K)
local
	M, A,dimsSet,dimsM,Fun,numevens,ArrayIt,ArrayAssignDim,
        ArrayValueDim,dim;

if IsHapPureCubicalComplex(MM) then
M:=PureCubicalComplexToCubicalComplex(MM);
else
M:=MM;
fi;

if not IsHapCubicalComplex(M) then
Print("This functions must be applied to a cubical, or pure cubical, complex.\n");
return fail;
fi;

dim:=Dimension(M);
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayAssignDim:=ArrayAssignFunctions(dim);
ArrayIt:=ArrayIterate(dim);
dimsM:=EvaluateProperty(M,"arraySize");
dimsSet:=List([1..dim],x->[1..dimsM[x]]);


##############################
numevens:=function(x);
return Length(Filtered(x,i->IsEvenInt(i)));
end;
##############################

#######################
Fun:=function(x) local y,z;
if ArrayValueDim(M!.binaryArray,x)=1 then
y:=numevens(x);
if y<=K then
 ArrayAssignDim(A,x,1);
fi;
fi;
end;
#######################

A:=StructuralCopy(M!.binaryArray);
A:=0*A;
ArrayIt(dimsSet,Fun);

return Objectify(HapCubicalComplex,
rec(
   binaryArray:=A,
   properties:=StructuralCopy(M!.properties)));

end);
##########################################
##########################################

####################################################
####################################################
InstallGlobalFunction(SuspensionOfPureCubicalComplex,
function(M)
local
	A,B,C;

A:=M!.binaryArray;
B:=0*A;
B:=1+B;
C:=[B,A,B];
Apply(C,a->StructuralCopy(a));

return PureCubicalComplex(C);
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(ReadLinkImageAsPureCubicalComplex,
function(arg)
local file,thk,M, arcs, S, T, K, tsld, C, sC, sCnew;

file:=arg[1];
if Length(arg)>1 then thk:=arg[2]; else thk:=5; fi;
 
M:=ReadImageAsPureCubicalComplex(file,1/2);
tsld:=(Maximum(Flat(M))+Minimum(Flat(M)))/2;
M:=ArrayToPureCubicalComplex(M,tsld);
M:=PureCubicalComplex(FrameArray(M!.binaryArray));
M:=PureCubicalComplex(FrameArray(M!.binaryArray));

if Length(arg)=1 then
sC:=0;;
C:=PureComplexComplement(M);; sCnew:=Size(C); thk:=1;

while sC<sCnew do
thk:=thk+1;
C:=ThickenedPureComplex(C);
sC:=sCnew; sCnew:=Size(C);
od;

fi;;

arcs:=PathComponentOfPureCubicalComplex(M,0);
S:=SingularitiesOfPureCubicalComplex(M,thk,50);
while PathComponentOfPureCubicalComplex(S,0)>2*arcs or
Length(Homology(S,1))>0 do
S:=ThickenedPureCubicalComplex(S);
od;

#Sanity check
if not PathComponentOfPureCubicalComplex(S,0)=2*arcs then return fail; fi;

T:=ThickenedPureCubicalComplex(S);
while PathComponentOfPureCubicalComplex(T,0)>arcs or
Length(Homology(T,1))>0 do
T:=ThickenedPureCubicalComplex(T);
od;

#Sanity check
if not PathComponentOfPureCubicalComplex(T,0)=arcs then return fail; fi;

K:=PureCubicalComplex(
   FrameArray([T!.binaryArray,S!.binaryArray,S!.binaryArray,M!.binaryArray])
   );

K:=PureCubicalComplex(FrameArray(K!.binaryArray));


return K;
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(ReadLinkImageAsGaussCode,
function(arg)
local file,thk,M, arcs, edges, vertices, S, T, K, tsld, C, sC, sCnew, 
EdgeArc, MT, Boundaries, TBoundaries, i,j,P,Q, R, Gauss, a, e, v, OverUnder,
edge, V, pairs, cnt, ii, A, L, mx, mn, Orientation, Heads, Signs;

file:=arg[1];
if Length(arg)>1 then thk:=arg[2]; else thk:=5; fi;

M:=ReadImageAsPureCubicalComplex(file,1/2);
tsld:=(Maximum(Flat(M))+Minimum(Flat(M)))/2;
M:=ArrayToPureCubicalComplex(M,tsld);
A:=M!.binaryArray;
L:=Filtered([1..Length(A)],i->not IsZero(A[i]));
mn:=Maximum(L[1]-1,1);
mx:=Minimum(L[Length(L)]+1, Length(A));
A:=A{[mn..mx]};
A:=1*TransposedMat(A);
L:=Filtered([1..Length(A)],i->not IsZero(A[i]));
mn:=Maximum(L[1]-1,1);
mx:=Minimum(L[Length(L)]+1, Length(A));
A:=A{[mn..mx]};
A:=1*TransposedMat(A);

M!.binaryArray:=A;
M:=PureCubicalComplex(FrameArray(M!.binaryArray));
M:=PureCubicalComplex(FrameArray(M!.binaryArray));

if Length(arg)=1 then
sC:=0;;
C:=PureComplexComplement(M);; sCnew:=Size(C); thk:=1;

while sC<sCnew do
thk:=thk+1;
C:=ThickenedPureComplex(C);
sC:=sCnew; sCnew:=Size(C);
od;

fi;;

arcs:=PathComponentOfPureCubicalComplex(M,0);
S:=SingularitiesOfPureCubicalComplex(M,thk,50);
while PathComponentOfPureCubicalComplex(S,0)>2*arcs or
Length(Homology(S,1))>0 do
S:=ThickenedPureCubicalComplex(S);
od;

#Sanity check
if not PathComponentOfPureCubicalComplex(S,0)=2*arcs then return fail; fi;

T:=ThickenedPureCubicalComplex(S);
while PathComponentOfPureCubicalComplex(T,0)>arcs or
Length(Homology(T,1))>0 do
T:=ThickenedPureCubicalComplex(T);
od;
T:=ThickenedPureCubicalComplex(T);
T:=ThickenedPureCubicalComplex(T);

#Sanity check
if not PathComponentOfPureCubicalComplex(T,0)=arcs then return fail; fi;

vertices:=PathComponentOfPureCubicalComplex(T,0);
MT:=PureComplexDifference(M,T);
MT:=PureComplexThickened(MT);
edges:=PathComponentOfPureCubicalComplex(MT,0);
if not 2*vertices = edges then return fail; fi;

EdgeArc:=NullMat(edges,arcs);
Boundaries:=NullMat(edges,vertices);
for e in [1..edges] do
P:=PathComponentOfPureCubicalComplex(MT,e);

for a in [1..arcs] do
if Size(PureComplexIntersection(PathComponentOfPureCubicalComplex(M,a),P))>0
then EdgeArc[e][a]:=1; fi;
od;

for v in [1..vertices] do
Q:=PathComponentOfPureCubicalComplex(T,v);
R:=PureComplexIntersection(P,Q);
if Size(R) >0  then 
Boundaries[e][v]:=1;
fi;
od;
od;

TBoundaries:=TransposedMat(Boundaries);
##############################
OverUnder:=function(c)
local quad, quad1, quad2, pairs, i;
quad:=Filtered([1..edges], e->TBoundaries[c][e]=1);
quad1:=List(quad,e->Position(EdgeArc[e],1));
quad2:=Collected(quad1);
pairs:=[[],[]];
for i in [1..Length(quad)] do
if [quad1[i],2] in quad2 then Add(pairs[1],quad[i]);
else Add(pairs[2],quad[i]); fi;
od;
if Length(pairs[2])=1 then 
   if Sum(Boundaries[pairs[1][1]])=1 then pairs[1][1]:=1*pairs[2][1]; fi; 
   if Sum(Boundaries[pairs[1][2]])=1 then pairs[1][2]:=1*pairs[2][1]; fi;
fi;
return pairs;
end;
##############################

Heads:=List([1..vertices],i->[]);

Gauss:=[1];
edge:=OverUnder(1);
edge:=edge[1][1];
Add(Heads[1],edge);
i:=0;
while i<2*vertices do
   i:=i+1;
   V:=Filtered([1..vertices],v->Boundaries[edge][v]=1);
   V:=Filtered(V,v->not v=AbsInt(Gauss[i]));
   V:=V[1];
   pairs:=OverUnder(V);
   if edge in pairs[1] then 
      Add(Gauss,V);
      if Length(pairs[2])=1 then Add(Gauss,-V); i:=i+1; fi;
      if pairs[1][1]=edge then edge:=pairs[1][2]; Heads[V][1]:=edge; else edge:=pairs[1][1]; Heads[V][1]:=edge; fi;
   else
      Add(Gauss,-V);
      if pairs[2][1]=edge then edge:=pairs[2][2]; Heads[V][2]:=edge; else edge:=pairs[2][1]; Heads[V][2]:=edge;fi;
   fi;
od;

Gauss:=Gauss{[2..2*vertices+1]}; #A silly fudge!!!

#########################################
Orientation:=function(c)
local maxX, maxY, minX, minY, P, OU, OU2, L, A, i, j, ornt,orient, EDG;
P:=PathComponentOfPureCubicalComplex(T,c);
A:=P!.binaryArray;
L:=Filtered([1..Length(A)],i->not IsZero(A[i]));
minY:=Minimum(L);
maxY:=Maximum(L);
A:=TransposedMat(A);
L:=Filtered([1..Length(A)],i->not IsZero(A[i]));
minX:=Minimum(L);
maxX:=Maximum(L);
OU:=SSortedList(Flat(OverUnder(c)));
OU2:=List(OU,i->PathComponentOfPureCubicalComplex(MT,i));
ornt:=[];
for i in [minX..maxX] do
for j in [1..Length(OU)] do
if OU2[j]!.binaryArray[minY][i]>0 then Add(ornt,OU[j]); fi;
od;
od;
for i in [minY..maxY] do
for j in [1..Length(OU)] do
if OU2[j]!.binaryArray[i][maxX]>0 then Add(ornt,OU[j]); fi;
od;
od;
for i in Reversed([minX..maxX]) do
for j in [1..Length(OU)] do
if OU2[j]!.binaryArray[maxY][i]>0 then Add(ornt,OU[j]); fi;
od;
od;
for i in Reversed([minY..maxY]) do
for j in [1..Length(OU)] do
if OU2[j]!.binaryArray[i][minX]>0 then Add(ornt,OU[j]); fi;
od;
od;

orient:=[];
for i in ornt do
if not i in orient then Add(orient,i); fi;
od;

Add(orient,orient[1]);
return orient;

end;
#########################################
L:=List([1..vertices],Orientation);

Signs:=[];
for i in [1..vertices] do
if Length(Heads[i])<2 then Signs[i]:=1;
else 
for j in [1..4] do
if [ L[i][j] , L[i][j+1] ] = Heads[i] then Signs[i]:=1; break; fi;
if [ L[i][j+1] , L[i][j] ] = Heads[i] then Signs[i]:=-1; break; fi;
od;
fi;
od;

return [[Gauss],Signs];
end);
###################################################
###################################################


##########################################
##########################################
#HenonOrbit:=function(x,A,B,K,M,L)
InstallGlobalFunction(HenonOrbit,
function(arg)
local x,A,B,K,M,L,henon, Henon, Roundoff,P,N,DIMS;

x:=arg[1];
A:=arg[2];
B:=arg[3];
if Length(arg)>3 then K:=arg[4];
else K:=10^6;fi;
if Length(arg)>4 then M:=arg[5];
else M:=1800;fi;
if Length(arg)>5 then L:=arg[6];
else L:=20;fi;
N:=10^L;

#############################
henon:=function(x)
local z;
z:=[];
z[1]:=x[2]+1-A*x[1]^2;
z[2]:=B*x[1];
return z;
end;
#############################


#############################
Roundoff:=function(x)
local r1,r2;
r1:=N*x[1];
r1:=Int(r1);
r1:=r1/N;

r2:=N*x[2];
r2:=Int(r2);
r2:=r2/N;

return [r1,r2];
end;
#############################

#############################
DIMS:=function(xx,k)
local Xmax, Ymax, Xmin, Ymin, x,i;

x:=Roundoff(xx);
Xmax:=Int(M*(x[1]+2));
Xmin:=Int(M*(x[1]-2));
Ymax:=Int(M*(x[2]+2));
Ymin:=Int(M*(x[2]-2));
for i in [1..k] do
x:=Roundoff(henon(x));
Xmax:=Maximum(Xmax,Int(M*(x[1]+2)));
Ymax:=Maximum(Ymax,Int(M*(x[2]+2)));
Xmin:=Minimum(Xmin,Int(M*(x[1]-2)));
Ymin:=Minimum(Ymin,Int(M*(x[2]-2)));
od;
return [Xmin,Xmax,Ymin,Ymax];
end;
#############################

if Length(arg)<=6 then
#############################
Henon:=function(x,k)
local z,i,A, D;

D:=DIMS(x,k);
A:=NullMat(D[2]-D[1]+10,D[4]-D[3]+10);

x:=Roundoff(x);

for i in [1..2000] do
x:=Roundoff(henon(x));
od;

for i in [1..k] do
x:=Roundoff(henon(x));

A[Int(M*(x[1]+2))-D[1]+1][Int(D[4]-M*(x[2]+2))+1]:=1;
od;

return A;
end;
#############################
else
#############################
Henon:=function(x,k)
local z,i,A, D, DD;

DD:=arg[7];
D:=[];
D[2]:=Int(M*(DD[2]+2));
D[4]:=Int(M*(DD[4]+2));
D[1]:=Int(M*(DD[1]-2));
D[3]:=Int(M*(DD[3]-2));

A:=NullMat(D[2]-D[1]+10,D[4]-D[3]+10);

x:=Roundoff(x);

for i in [1..2000] do
x:=Roundoff(henon(x));
od;

for i in [1..k] do
x:=Roundoff(henon(x));
if x[1]>DD[1] and x[1]<DD[2] and x[2]>DD[3] and x[2]<DD[4] then

A[Int(M*(x[1]+2))-D[1]+1][Int(D[4]-M*(x[2]+2))+1]:=1;

fi;
od;

return A;
end;
#############################
fi;
P:=Henon(x,K);
P:=PureCubicalComplex(P);
P:=CropPureCubicalComplex(P);
P:=FrameArray(P!.binaryArray);
P:=PureCubicalComplex(P);

return P;

end);
##############################################
##############################################

###################################################
###################################################
InstallGlobalFunction(RandomCubeOfPureCubicalComplex,
function(M)
local
        B,
        cart, CART, dim,dims,
        ArrayValueDim,
        ArrayAssignDim,
        x,z,cnt, ran;

#############################################
if not IsHapPureCubicalComplex(M) then
Print("This function must be applied to a pure cubical complex.\n");
return fail;
fi;
#############################################

dim:=Dimension(M);
ran:=Random([1..Size(M)]);
cnt:=0;

ArrayValueDim:=ArrayValueFunctions(dim);
ArrayAssignDim:=ArrayAssignFunctions(dim);
dims:=EvaluateProperty(M,"arraySize");
B:=0*M!.binaryArray;
CART:=Cartesian(List([1..dim],a->[1..dims[a]]));
cart:=Cartesian(List([1..dim],a->[-1,0,1]));

for x in CART  do
if ArrayValueDim(M!.binaryArray,x)=1 then cnt:=cnt+1; fi;
if cnt=ran then ArrayAssignDim(B,x,1); break;fi;
od;

return PureCubicalComplex(B);

end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(FramedPureCubicalComplex,
function(M);

return PureCubicalComplex(FrameArray(M!.binaryArray));

end);
#################################################################
#################################################################


