#(C) 2009 Graham Ellis



#################################################################
#################################################################
InstallGlobalFunction(SparseChainComplexOfCubicalComplex,
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
nn:=[];

for x in poscells do
  Add(nn, [x,-1]);
od;
for x in negcells do
  Add(nn, [x,1]);
od;

return nn;
end;
#######################################

return
Objectify(HapSparseChainComplex,
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
InstallGlobalFunction(SparseChainComplexOfCubicalPair,
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

nn:=[];

for x in poscells do
  Add(nn, [x,-1]);
od;
for x in negcells do
  Add(nn, [x,1]);
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
Objectify(HapSparseChainComplex,
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
InstallGlobalFunction(SparseChainMapOfCubicalPairs,
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
C:=SparseChainComplexOfPair(X,S);
D:=SparseChainComplexOfPair(Y,T);
fi;

if Length(arg)=2 then
C:=arg[1];
D:=arg[2];
fi;
####################################

###############################################
map:=function(V,n)
local x,w,v,W;

W:=[];

for v in V do
x:=C!.generator2Coordinates(n,v[1]);
x:=D!.coordinates2Generator(x);
if not v[2]=0 then Add(W,[x,v[2]]);   fi;
od;

return W;
end;
###############################################

return Objectify(HapSparseChainMap,
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




#####################################################################
#####################################################################
InstallMethod(SparseChainComplexOfPair,
"Cellular chain complex of a pair of cubical complexes",
[IsHapCubicalComplex,IsHapCubicalComplex],
function(M,S)
return SparseChainComplexOfCubicalPair(M,S);;
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallOtherMethod(SparseChainComplexOfPair,
"Cellular chain complex of a pair of pure cubical complexes",
[IsHapPureCubicalComplex,IsHapPureCubicalComplex],
function(MM,SS)
local M,S,P;
P:=ExcisedPureCubicalPair(MM,SS);
M:=P[1];S:=P[2];
#M:=MM;S:=SS;
return SparseChainComplexOfCubicalPair(
PureCubicalComplexToCubicalComplex(M),
PureCubicalComplexToCubicalComplex(S));;
end);
#####################################################################
#####################################################################

#################################################################
#################################################################
InstallGlobalFunction(FilteredChainComplexToFilteredSparseChainComplex,
function(C)
local Boundary;

########################
Boundary:=function(n,k)
local v,i,bnd;
v:=C!.boundary(n,k);
bnd:=[];

for i in [1..Length(v)] do
if not IsZero(v[i]) then
Add(bnd, [i,v[i]]);
fi;
od;

return bnd;
end;
########################

return
Objectify(HapFilteredSparseChainComplex,
           rec(
           dimension:=C!.dimension,
           boundary:=Boundary,
           filteredDimension:=C!.filteredDimension,
           properties:=[
           ["length",EvaluateProperty(C,"length")],
           ["type","chainComplex"],
           ["characteristic",0],
           ["filtration_length",EvaluateProperty(C,"filtration_length")]]
           ));

end);
#####################################################################
#####################################################################



#################################################################
#################################################################
InstallGlobalFunction(SparseFilteredChainComplexOfFilteredCubicalComplex,
function(arg)
local M,faces,dims,dim,dim1,x,y,z,n,r,numevens,Dimsn,Boundary,BinLst,LstBin,
ArrayValueDim,ArrayValueDim1,ArrayAssignDim,ArrayIt,dimsSet,Fun,LF,
flen, filtmat,  permutation, invpermutation, cnt, FilteredDimension;
# Inputs a cubical complex and returns the associated cellular
# chain complex.

M:=arg[1];

if not IsHapFilteredCubicalComplex(M) then
Print("This function must be applied to a filtered cubical complex.\n");
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

flen:=Maximum(Flat(M!.filtration));
filtmat:=[];
for n in [0..dim] do
filtmat[n+1]:=List([1..flen],x->[]);
od;

#######################
Fun:=function(x) local y,z,k;
if ArrayValueDim(M!.binaryArray,x)=1 then
y:=numevens(x)+1;
faces[y]:=faces[y]+1;
ArrayAssignDim(BinLst,x,faces[y]);
LstBin[y][faces[y]]:=x;
k:=ArrayValueDim(M!.filtration,x);
Add(filtmat[y][k],faces[y]);
fi;
end;
#######################
ArrayIt(dimsSet,Fun);

permutation:=[];  #old --> new
invpermutation:=[];#new --> old
for n in [1..dim+1] do
cnt:=0;
permutation[n]:=[];
invpermutation[n]:=[];
for r in [1..flen] do
for x in filtmat[n][r] do
cnt:=cnt+1;
permutation[n][x]:=cnt; 
invpermutation[n][cnt]:=x;
od;
od;
od;

#######################################
Dimsn:=function(n);
return faces[n+1];
end;
#######################################

LF:=[];
for x in [1..Length(faces)] do
LF[x]:=0*[1..faces[x]];
od;

#######################################
Boundary:=function(n,jj)
local x,poscells,negcells,nn,a,b,cnt,j;
poscells:=[];
negcells:=[];
j:=invpermutation[n+1][jj];
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
nn:=[];

for x in poscells do
  Add(nn, [x,-1]);
od;
for x in negcells do
  Add(nn, [x,1]);
od;
nn:=List(nn,x->[permutation[n][x[1]],x[2]]);
return nn;
end;
#######################################

#######################################
FilteredDimension:=function(r,n);

return 
Sum(List([1..r],i->Length(filtmat[n+1][i]) ));

end;
#######################################


return
Objectify(HapFilteredSparseChainComplex,
           rec(
           dimension:=Dimsn,
           boundary:=Boundary,
           filteredDimension:=FilteredDimension,
           coordinateToPosition:=BinLst,
           positionToCoordinate:=LstBin,
           properties:=[
           ["length",EvaluateProperty(M,"dimension")],
           ["type","chainComplex"],
           ["characteristic",0],
           ["filtration_length",flen]]
           ));

end);
#################################################################
#################################################################


###############################################################
###############################################################
InstallGlobalFunction(PersistentHomologyOfFilteredSparseChainComplex,
function(C,N)
local BN,BN1, BNlist, BN1list, flen, fdimsN, fdimsN1,
      PHmat, L, n,i,j, r,s;

flen:=EvaluateProperty(C,"filtration_length");
fdimsN:=List([1..flen],i->C!.filteredDimension(i,N));
fdimsN1:=List([1..flen],i->C!.filteredDimension(i,N+1));

if N>0 then
BN:=SparseBoundaryMatrix(C,N);

#Print("Before: ",Sum(List(Filtered(BN!.mat,r->IsBound(r)),Length)),"\n");

for r in  fdimsN do
if IsBound(BN!.heads) then Unbind(BN!.heads); fi;
BN!.rows:=r;
SparseSemiEchelon(BN);
od;

#Print("After: ",Sum(List(Filtered(BN!.mat,r->IsBound(r)),Length)),"\n");

BNlist:=[];
for n in [1..BN!.rows] do
if IsBound(BN!.mat[n]) then
BNlist[n]:=Maximum(List(BN!.mat[n],x->x[1]));
else BNlist[n]:=0;
fi;
od;
Unbind(BN);

#else
#BNlist:=fdimsN;

fi;

BN1:=SparseBoundaryMatrix(C,N+1);
#Print("Before: ",Sum(List(Filtered(BN1!.mat,r->IsBound(r)),Length)),"\n");

ReverseSparseMat(BN1);
for r in  fdimsN1 do
if IsBound(BN1!.heads) then Unbind(BN1!.heads); fi;
BN1!.rows:=r;
SparseSemiEchelon(BN1);
od;
ReverseSparseMat(BN1);

#Print("After: ",Sum(List(Filtered(BN1!.mat,r->IsBound(r)),Length)),"\n");

BN1list:=[];
for n in [1..BN1!.rows] do
if IsBound(BN1!.mat[n]) then
BN1list[n]:=Maximum(List(BN1!.mat[n],x->x[1]));
else BN1list[n]:=0;
fi;
od;
Unbind(BN1);

PHmat:=NullMat(flen,flen);

for i in [1..flen] do
for j in [i..flen] do
r:=fdimsN[i]; s:=fdimsN1[j];
if N>0 then
L:=Filtered(BNlist{[1..r]}, x->x=0);
PHmat[i][j]:=  Length(L);
else
PHmat[i][j]:= r;
fi;
L:=Filtered(BN1list{[1..s]}, x-> x<=r and (not x=0));
PHmat[i][j]:=PHmat[i][j]-Length(L);
od;od;

return PHmat;

end);
###############################################################
###############################################################





