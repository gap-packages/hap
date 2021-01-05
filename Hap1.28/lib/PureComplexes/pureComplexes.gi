#(C) 2009 Graham Ellis

#############################
#############################
InstallGlobalFunction(RandomCellOfPureComplex,
function(M)
local dims, A, B,x,U,v,  i;

if not IsPureComplex(M) then
Print("Function must be applied to a pure cubical or pure permutahedral complex.\n");
return fail;
fi;

if ArraySum(M!.binaryArray)=0 then return fail; fi;

dims:=EvaluateProperty(M,"arraySize");
dims:=Reversed(dims);
A:=M!.binaryArray;
B:=0*A;

while true do
U:=A;
v:=[];
for i in [1..Length(dims)] do
x:=Random([1..dims[i]]);
Add(v,x);
U:=U[x];
od;
if U=1 then break;  fi;
od;

U:=B;
for i in v{[1..Length(v)-1]} do
U:=U[i];
od;
U[v[Length(v)]]:=1;

return PureComplex(M,B);
end);
#############################
#############################

#############################
#############################
InstallGlobalFunction(IsPureComplex,
function(M);

if IsHapPureCubicalComplex(M) or IsHapPurePermutahedralComplex(M) then
return true;
fi;

return false;

end);
#############################
#############################


#############################
#############################
InstallGlobalFunction(PureComplex,
function(M,A)
local       record;

record:=   rec(
           binaryArray:=A*1,
           properties:=[
           ["dimension",ArrayDimension(A)],
           ["arraySize",ArrayDimensions(A)]]
           );

if IsHapPureCubicalComplex(M) then
return Objectify(HapPureCubicalComplex,record); fi;

if IsHapPurePermutahedralComplex(M) then
return Objectify(HapPurePermutahedralComplex,record); fi;

return fail;

end);
#############################
#############################


#############################
#############################
InstallGlobalFunction(UnitCubicalBall,
function(dim) local Ball;
Ball:=Cartesian(List([1..dim],i->[-1,0,1]));
RemoveSet(Ball,List([1..dim],i->0));
return Ball;
end);
#############################
#############################

#############################
#############################
InstallGlobalFunction(UnitPermutahedralBall,
function(dim)
local  n,i,B,U,A;

if dim=2 then return [[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0]]; fi;
if dim=3 then return
[[0,0,-1],[1,0,-1],[0,1,-1],[1,1,-1],[0,-1,0],[1,-1,0], [-1,0,0],[1,0,0],[-1,1,0],[0,1,0],[-1,-1,1],[0,-1,1], [-1,0,1],[0,0,1]];
fi;
if dim=4 then return [[0,0,0,1],[0,0,1,1],[0,1,1,1],[0,1,0,1],[1,0,0,0], [1,-1,0,0],[1,0,-1,0], [1,-1,-1,0],[1,-1,-1,-1],[1,0,-1,-1],[1,-1,0,-1], [0,-1,-1,-1],[1,0,0,-1],[0,0,-1,-1], [0,-1,0,-1],[0,-1,-1,0],[0,0,0,-1], [0,0,-1,0],[0,-1,0,0],[0,1,0,0],[0,1,1,0],[-1,1,0,0], [-1,1,1,0],[-1,1,0,1], [0,0,1,0],[-1,0,1,0],[-1,0,0,0],[-1,0,0,1],[-1,0,1,1],[-1,1,1,1]]; fi;
#dim>4 do the following

n:=dim+1;

A:=List([1..n],i->List([1..n],j->1));
for i in [1..n] do
A[i][i]:=-n+1;
od;

U:=Filtered(Combinations(A),x->not Length(x) in [0,n]);;
U:=List(U,x->Sum(x));

U:=SolutionsMatDestructive(A,U);
U:=List(U,x->x{[1..n-1]});
return U;
end);
#############################
#############################

#############################
#############################
InstallGlobalFunction(UnitBall,
function(M);
if IsHapPureCubicalComplex(M) then return UnitCubicalBall(Dimension(M));fi;
if IsHapPurePermutahedralComplex(M) then return 
SSortedList(UnitPermutahedralBall(Dimension(M))); fi;
return UnitPermutahedralBall(Dimension(M)); 	#This is so Fintan's 
						#package still works!
end);
#############################
#############################

#####################################################################
#####################################################################
InstallMethod(Nerve,
"Nerve of lattice complex",
[IsHapPureCubicalComplex],
function(M);
return PureComplexToSimplicialComplex(M);
end);

InstallMethod(Nerve,
"Nerve of lattice complex",
[IsHapPureCubicalComplex,IsInt],
function(M,n);
return PureComplexToSimplicialComplex(M,n);
end); 

InstallMethod(Nerve,
"Nerve of lattice complex",
[IsHapPurePermutahedralComplex],
function(M);
return PureComplexToSimplicialComplex(M);
end); 

InstallMethod(Nerve,
"Nerve of lattice complex",
[IsHapPurePermutahedralComplex,IsInt],
function(M,n);
return PureComplexToSimplicialComplex(M,n);
end);


#####################################################################
#####################################################################




##################################################
##################################################
InstallGlobalFunction(PureComplexToSimplicialComplex,
function(arg)
local M, DIM, AO,A,dim,dims,
      ArrayValueDim,
      #CartProd,
      dimSet, ArrayIt, FN,
      Vertices, VertexCoordinates,ArrayValueDim1,
      Ball, Balls,
      SimplicesLst, Simplices, NrSimplices, EnumeratedSimplex,
       b, i, j, t, t1, t2, v, x, y;

M:=arg[1];
if Length(arg)=2 then DIM:=arg[2];
else DIM:=Dimension(M); fi;

#################################
if not IsPureComplex(M) then
Print("This function must be applied to a pure cubical or pure permutahedral complex.\n");
return fail; fi;
#################################

AO:=FrameArray(M!.binaryArray);
A:=StructuralCopy(AO);
dim:=ArrayDimension(A);
dims:=ArrayDimensions(A);
Vertices:=0;
VertexCoordinates:=[];
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim-1);
#CartProd:=Cartesian(List([1..dim],a->[2..dims[a]-1]));


Ball:=UnitBall(M);


#############################
Balls:=[];
Balls[1]:=Ball;
for t in [2..DIM] do
  Balls[t]:=Cartesian(Balls[t-1],Ball);
  if t>2 then
    Balls[t]:=List(Balls[t],x->Concatenation(x[1],[x[2]]));
  fi;
  Balls[t]:=Filtered(Balls[t],x->x[t-1]>x[t]);
  for i in [1..t-1] do
    Balls[t]:=Filtered(Balls[t],x->x[i]-x[t] in Ball);
  od;
od;
#############################

#for x in CartProd do
FN:=function(x)
local y;
  if ArrayValueDim(AO,x)=1 then Vertices:=Vertices+1;
    y:=ArrayValueDim1(A,x{[2..dim]});
    y[x[1]]:=Vertices;
    VertexCoordinates[Vertices]:=x;
  fi;
end;
#od;

dimSet:=List([1..dim],x->[2..dims[x]-1]);
ArrayIt:=ArrayIterate(dim);
ArrayIt(dimSet,FN);


Vertices:=[1..Vertices];
SimplicesLst:=List([1..1000],i->[]);  #VERY SLOPPY!!!

if DIM>=0 then
  SimplicesLst[1]:=List(Vertices,i->[i]);
fi;

if DIM>=1 then
for v in Vertices do
  x:=VertexCoordinates[v];
  for b in Ball do
    t:= ArrayValueDim(A,b+x);
    if t>v then Add(SimplicesLst[2],[v,t]); fi;
  od;
od;
fi;


if DIM>=2 then
for j in [2..DIM] do

for v in Vertices do
  x:=VertexCoordinates[v];
  for b in Balls[j] do
    t:=List([1..j],i->ArrayValueDim(A,b[i]+x));
    if not 0 in t then
       Add(SimplicesLst[j+1],SortedList(Concatenation([v],t)));
    fi;
  od;
od;
SimplicesLst[j+1]:=SSortedList(SimplicesLst[j+1]);
od;
fi;

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

#############################################
EnumeratedSimplex:=function(v);
return Position(SimplicesLst[Length(v)],v);
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
##################################################
##################################################


#################################################################
#################################################################
InstallGlobalFunction(ThickenedPureComplex,
function(M)
local
        B,
        cart, CART, dim,dim1,dims,
        Thicken,
        ArrayValueDim,
        ArrayValueDim1,
        dimSet, ArrayIt,
        x,z, record;

#################################
if not IsPureComplex(M) then
Print("This function must be applied to a pure cubical or pure permutahedral complex.\n");
return fail; fi;
#################################

dim:=Dimension(M);

if dim=2 and IsHapPureCubicalComplex(M) then 
return ThickenedPureCubicalComplex_dim2(M); fi;

dim1:=dim-1;
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim1);
dims:=EvaluateProperty(M,"arraySize");
B:=StructuralCopy(M!.binaryArray);
#cart:=Cartesian(List([1..dim],a->[-1,0,1]));
cart:=UnitBall(M);

########################
Thicken:=function(y)
local  x, z, w;

if ArrayValueDim(M!.binaryArray,y)=0 then return false;fi;

for x in cart do
z:=x+y;
if (not 0 in z) and (not -1 in dims - z) then
  w:=ArrayValueDim1(B,z{[2..Length(z)]});
w[z[1]]:=1;
fi;
od;

end;
########################


dimSet:=List([1..dim],x->[1..dims[x]]);
ArrayIt:=ArrayIterate(dim);
ArrayIt(dimSet,Thicken);

return PureComplex(M,B);

end);
#################################################################
#################################################################

PureComplexThickened:=ThickenedPureComplex;
MakeReadOnlyGlobal("PureComplexThickened");

#################################################################
#################################################################
InstallGlobalFunction(ThickenedPureCubicalComplex,
function(M)
return ThickenedPureComplex(M);
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(ComplementOfPureComplex,
function(M)
local
        B,
        CART, dim,dim1,dims,
        Opp,
        ArrayValueDim,
        ArrayValueDim1,
        dimSet,ArrayIt,
        x,z;

#################################
if not IsPureComplex(M) then
Print("This function must be applied to a pure cubical or pure permutahedral complex.\n");
 return fail; fi;
#################################

dim:=Dimension(M);
dim1:=dim-1;
dims:=EvaluateProperty(M,"arraySize");
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim1);
B:=StructuralCopy(M!.binaryArray);

########################
Opp:=function(y)
local z;
z:=ArrayValueDim1(B,y{[2..Length(y)]});
if ArrayValueDim(M!.binaryArray,y)=0 then
z[y[1]]:=1;
else
z[y[1]]:=0;
fi;
end;
########################

dimSet:=List([1..dim],x->[1..dims[x]]);
ArrayIt:=ArrayIterate(dim);
ArrayIt(dimSet,Opp);

return PureComplex(M,B);

end);
#################################################################
#################################################################

PureComplexComplement:=ComplementOfPureComplex;
MakeReadOnlyGlobal("PureComplexComplement");

PureComplexBoundary:=BoundaryOfPureComplex;
MakeReadOnlyGlobal("PureComplexBoundary");


#################################################################
#################################################################
InstallGlobalFunction(ComplementOfPureCubicalComplex,
function(M)
return ComplementOfPureComplex(M);
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(PureComplexUnion,
function(M,N)
local
        D,
        dim,dims,
        fn,
        ArrayValueDim,
        ArrayValueDim1,
        ArrayIt, dimSet,
        x,w,d;

###################################
if not
(IsHapPureCubicalComplex(M)
and
IsHapPureCubicalComplex(N))
or
(IsHapPurePermutahedralComplex(M)
and
IsHapPurePermutahedralComplex(N))
then
Print("This function must be applied to a pair of pure cubical complexes or a pair of pure permutahedral complexes.\n");
return fail;
fi;
if not
EvaluateProperty(M,"arraySize")=
EvaluateProperty(N,"arraySize")
then
Print("The pure complexes have different array sizes.\n");
return fail;
fi;
###################################

D:=PureComplex(M,M!.binaryArray);
dim:=Dimension(D);
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim-1);
dims:=EvaluateProperty(D,"arraySize");


#for x in CART do
####################
fn:=function(x);
if ArrayValueDim(N!.binaryArray,x)=1 then
w:=ArrayValueDim1(D!.binaryArray,x{[2..dim]});
w[x[1]]:=1;
fi;
end;
####################
#od;

dimSet:=List([1..dim],x->[1..dims[x]]);
ArrayIt:=ArrayIterate(dim);
ArrayIt(dimSet,fn);

return D;
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(PureCubicalComplexUnion,
function(M,N)
return PureComplexUnion(M,N);
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(PureComplexIntersection,
function(M,N)
local
        D,
        dim,dims,
        fn, dimSet,
        ArrayValueDim,
        ArrayValueDim1,
        ArrayIt,
        x,w,d;

###################################
if not
((IsHapPureCubicalComplex(M)
and
IsHapPureCubicalComplex(N))
or
(IsHapPurePermutahedralComplex(M)
and
IsHapPurePermutahedralComplex(N))
)
then
Print("This function must be applied to a pair of pure cubical complexes or a pair of pure permutahedral complexes.\n");
return fail;
fi;
if not
EvaluateProperty(M,"arraySize")=
EvaluateProperty(N,"arraySize")
then
Print("The pure complexes have different array sizes.\n");
return fail;
fi;
###################################

D:=PureComplex(M,M!.binaryArray*0);
dim:=Dimension(D);
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim-1);
dims:=EvaluateProperty(D,"arraySize");


####################
fn:=function(x);
if ArrayValueDim(N!.binaryArray,x)=1 and ArrayValueDim(M!.binaryArray,x)=1 then
w:=ArrayValueDim1(D!.binaryArray,x{[2..dim]});
w[x[1]]:=1;
fi;
end;
####################

dimSet:=List([1..dim],x->[1..dims[x]]);
ArrayIt:=ArrayIterate(dim);
ArrayIt(dimSet,fn);

return D;
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(PureCubicalComplexIntersection,
function(M,N)
return PureComplexIntersection(M,N);
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(PureComplexDifference,
function(M,N)
local
        D,
        dim,dims,
        CART,
        ArrayValueDim,
        ArrayValueDim1,
        dimSet,ArrayIt, Opp,
        x,w,d;

###################################
if not
(
(IsHapPureCubicalComplex(M)
and
IsHapPureCubicalComplex(N))
or
(IsHapPurePermutahedralComplex(M)
and
IsHapPurePermutahedralComplex(N))
)
then
Print("This function must be applied to a pair of pure cubical complexes or a pair of pure permutahedral complexes.\n");
return fail;
fi;
if not
EvaluateProperty(M,"arraySize")=
EvaluateProperty(N,"arraySize")
then
Print("The pure complexes have different array sizes.\n");
return fail;
fi;
###################################

D:=PureComplex(M,M!.binaryArray);
dim:=Dimension(D);
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim-1);
dims:=EvaluateProperty(D,"arraySize");

#####################
Opp:=function(x)
local w;
if ArrayValueDim(N!.binaryArray,x)=1 then
w:=ArrayValueDim1(D!.binaryArray,x{[2..dim]});
w[x[1]]:=0;
fi;
end;
#####################

dimSet:=List([1..dim],x->[1..dims[x]]);
ArrayIt:=ArrayIterate(dim);
ArrayIt(dimSet,Opp);

return D;
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(PureCubicalComplexDifference,
function(M,N)
return PureComplexDifference(M,N);
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(BoundaryOfPureComplex,
function(M)
local
        B,
        cart,dim,dims,
        InBoundary,
        ArrayValueDim,
        ArrayAssignDim,ArrayIt,
        Fun,
        dimsSet,
        x,z;

#############################################
if IsHapRegularCWComplex(M) then
if not IsPureRegularCWComplex(M) then
Print("This function must be applied to a pure cubical, permutahedral or regular CW complex.\n");
return fail;
else
return
BoundaryOfPureRegularCWComplex(M);
fi;
fi;
#############################################


#############################################
if not IsPureComplex(M) then
Print("This function must be applied to a pure cubical, permutahedral or regular CW complex.\n");
return fail;
fi;
#############################################

dim:=Dimension(M);
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayAssignDim:=ArrayAssignFunctions(dim);
ArrayIt:=ArrayIterate(dim);
dims:=EvaluateProperty(M,"arraySize");
B:=M!.binaryArray*0;

dimsSet:=List([1..dim],a->[1..dims[a]]);
#cart:=Cartesian(List([1..dim],a->[-1,0,1]));
cart:=UnitBall(M);

########################
InBoundary:=function(y)
local  x, z;

if ArrayValueDim(M!.binaryArray,y)=0 then return false;fi;

for x in cart do
z:=x+y;
if (not 0 in z) and (not -1 in dims - z) then
   if ArrayValueDim(M!.binaryArray,z)=0 then return true; fi;
fi;
od;

return false;
end;
########################


####################
Fun:=function(x);
if InBoundary(x) then
ArrayAssignDim(B,x,1);
fi;
end;
####################

ArrayIt(dimsSet,Fun);

return PureComplex(M,B);

end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(BoundaryOfPureCubicalComplex,
function(M)
return BoundaryOfPureComplex(M);
end);
#################################################################
#################################################################


######################################################################
######################################################################
InstallGlobalFunction(HomotopyEquivalentMinimalPureSubcomplex,
function(T,S)
local A;

if 
((not IsHapPureCubicalComplex(T))
and
(not IsHapPureCubicalComplex(S)))
and 
((not IsHapPurePermutahedralComplex(T))
and
(not IsHapPurePermutahedralComplex(S)))
then
Print("This function can only be applied to pure cubical or permutahedral complexes.\n");
return fail;
fi;

#####################
if IsHapPureCubicalComplex(T) then
A:=HomotopyEquivalentSmallerSubArray(T!.binaryArray,S!.binaryArray);
return PureCubicalComplex(A);
fi;
#####################
#####################
if IsHapPurePermutahedralComplex(T) then
A:=HomotopyEquivalentSmallerSubPermArray(T!.binaryArray,S!.binaryArray);
return PurePermutahedralComplex(A);
fi;
#####################

end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(HomotopyEquivalentMinimalPureCubicalSubcomplex,
function(T,S);
return HomotopyEquivalentMinimalPureSubcomplex(T,S);
end);
######################################################################
######################################################################


######################################################################
######################################################################
InstallGlobalFunction(ContractPureComplex,
function(T);


#############################################
if not IsPureComplex(T) then
Print("This function must be applied to a pure cubical complex or pure permutahedral complex.\n");
return fail;
fi;
#############################################

if EvaluateProperty(T,"contracted")=true then return T; fi;

if IsHapPureCubicalComplex(T) then
T!.binaryArray:=ContractArray(T!.binaryArray);
Add(T!.properties,["contracted",true]);
return T;
fi;

if IsHapPurePermutahedralComplex(T) then
T!.binaryArray:=ContractPermArray(T!.binaryArray);
Add(T!.properties,["contracted",true]);
return T;
fi;

end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(ContractPureCubicalComplex,
function(T);
ContractPureComplex(T);
end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(HomotopyEquivalentMaximalPureSubcomplex,
function(T,S)
local A;

#############################################
if not IsPureComplex(T) then
Print("This function must be applied to a pure cubical complex or pure permutahedral complex.\n");
return fail;
fi;
#############################################


if IsHapPureCubicalComplex(T) then
A:=HomotopyEquivalentLargerSubArray(T!.binaryArray,S!.binaryArray);

return PureCubicalComplex(A);
fi;

if IsHapPurePermutahedralComplex(T) then
A:=HomotopyEquivalentLargerSubPermArray(T!.binaryArray,S!.binaryArray);
return PurePermutahedralComplex(A);
fi;


end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(HomotopyEquivalentMaximalPureCubicalSubcomplex,
function(T,S);
return HomotopyEquivalentMaximalPureSubcomplex(T,S);
end);
######################################################################
######################################################################

##################################################
##################################################
InstallGlobalFunction(CropPureComplex,
function(M)
local A,B,x,dim,dims,dimsSet,firsts,lasts,Fun,
      ArrayValueDim,ArrayIt,ArrayAssignDim,
      d, NewDimsSet;

#############################################
if not IsPureComplex(M) then
Print("This function must be applied to a pure cubical complex or pure permutahedral complex.\n");
return fail;
fi;
#############################################


A:=M!.binaryArray;
dim:=ArrayDimension(A);
dims:=ArrayDimensions(A);
dimsSet:=List(dims,d->[1..d]);
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayIt:=ArrayIterate(dim);
ArrayAssignDim:=ArrayAssignFunctions(dim);


firsts:=List([1..dim],i->infinity);
lasts:=List([1..dim],i->0);

################
Fun:=function(x);
if ArrayValueDim(A,x)=1 then
firsts:=List([1..dim],i->Minimum(x[i],firsts[i]));
lasts:=List([1..dim],i->Maximum(x[i],lasts[i]));
fi;
end;
################
ArrayIt(dimsSet,Fun);


NewDimsSet:=List([1..dim],n->[1..lasts[n]-firsts[n]+1]);
B:=0;
for d in [1..dim] do
B:=List(NewDimsSet[d],i->StructuralCopy(B));
od;

firsts:=firsts-List([1..Length(firsts)],i->1);

################
Fun:=function(x);
ArrayAssignDim(B,x,ArrayValueDim(A,x+firsts));
end;
################
ArrayIt(NewDimsSet,Fun);

if IsHapPureCubicalComplex(M) then
return PureCubicalComplex(B);
fi;

if IsHapPurePermutahedralComplex(M) then
return PurePermutahedralComplex(B);
fi;

end);
##################################################
##################################################

######################################################################
######################################################################
InstallGlobalFunction(CropPureCubicalComplex,
function(T);
return CropPureComplex(T);
end);
######################################################################
######################################################################

##################################################
##################################################
InstallGlobalFunction(BoundingPureComplex,
function(M)
local A,B,x,dim,dims,dimsSet,firsts,lasts,Fun,
      ArrayValueDim,ArrayIt,ArrayAssignDim;

#############################################
if not IsPureComplex(M) then
Print("This function must be applied to a pure cubical complex or pure permutahedral complex.\n");
return fail;
fi;
#############################################


A:=M!.binaryArray;
B:=A*0;
dim:=ArrayDimension(A);
dims:=ArrayDimensions(A);
dimsSet:=List(dims,d->[1..d]);
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayIt:=ArrayIterate(dim);
ArrayAssignDim:=ArrayAssignFunctions(dim);


firsts:=List([1..dim],i->infinity);
lasts:=List([1..dim],i->0);

################
Fun:=function(x);
if ArrayValueDim(A,x)=1 then
firsts:=List([1..dim],i->Minimum(x[i],firsts[i]));
lasts:=List([1..dim],i->Maximum(x[i],lasts[i]));
fi;
end;
################
ArrayIt(dimsSet,Fun);


################
Fun:=function(x);
if not false in List([1..dim],d-> firsts[d] <= x[d])
and
not false in List([1..dim],d-> lasts[d] >= x[d])
then
ArrayAssignDim(B,x,1);
fi;
end;
################
ArrayIt(dimsSet,Fun);

if IsHapPureCubicalComplex(M) then
return PureCubicalComplex(B);
fi;

if IsHapPurePermutahedralComplex(M) then
return PurePermutahedralComplex(B);
fi;

end);
##################################################
##################################################

######################################################################
######################################################################
InstallGlobalFunction(BoundingPureCubicalComplex,
function(T);
return BoundingPureComplex(T);
end);
######################################################################
######################################################################

#####################################################################
#####################################################################
InstallMethod(ContractedComplex,
"Contracted pure complex",
[IsObject],
function(M) local  C,A;

#############################################
if IsHapRegularCWComplex(M) then 
return ContractedRegularCWComplex(M);
fi;
#############################################


#############################################
if not IsPureComplex(M) then
Print("This function must be applied to a pure cubical complex or pure permutahedral complex or regular CW-complex.\n");
return fail;
fi;
#############################################


A:=StructuralCopy(M!.binaryArray);;
C:=rec();
C.properties:=StructuralCopy(M!.properties);

if IsHapPureCubicalComplex(M) then
C:=Objectify(HapPureCubicalComplex,C);
fi;

if IsHapPurePermutahedralComplex(M) then
C:=Objectify(HapPurePermutahedralComplex,C);
fi;


if not EvaluateProperty(M,"contracted")=true then
if IsHapPureCubicalComplex(M) then
A:=ContractArray(A);
fi;
if IsHapPurePermutahedralComplex(M) then
A:=ContractPermArray(A);
fi;
Add(C!.properties,["contracted",true]);
fi;

C!.binaryArray:=A;
return C;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallMethod(CochainComplex,
"Cochain complex for pure cubical complexes",
[IsHapPureCubicalComplex],
function(M) local  C;
C:=ChainComplex(M);
return HomToIntegers(C);
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallMethod(CochainComplex,
"Cochain complex for pure permutahedral complexes",
[IsHapPurePermutahedralComplex],
function(M) local  C;
C:=ChainComplex(M);
return HomToIntegers(C);
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallMethod(CochainComplex,
"Cochain complex for simplicial complexes",
[IsHapSimplicialComplex],
function(M) local  C;
C:=ChainComplex(M);
return HomToIntegers(C);
end);
#####################################################################
#####################################################################




##########################################
##########################################
InstallGlobalFunction(ZigZagContractedPureComplex,
function(arg)
local MM,A,B,M,N,i,d,dim;

MM:=arg[1];

#############################################
if not IsPureComplex(MM) then
Print("This function must be applied to a pure cubical complex or pure permutahedral complex.\n");
return fail;
fi;
#############################################


dim:=Dimension(MM);
M:=ContractedComplex(MM);
M:=CropPureComplex(M);

A:=M!.binaryArray;

   #########################
if IsHapPureCubicalComplex(MM) then
   for i in [2..Length(A)] do
      if A[i]=A[i-1] then A[i-1]:=0; fi;
   od;
   A:=Filtered(A,x-> not x=0);

   for d in [2..dim] do
   A:=PermuteArray(A,(1,d));
   for i in [2..Length(A)] do
      if A[i]=A[i-1] then A[i-1]:=0; fi;
   od;
   A:=Filtered(A,x-> not x=0);
   od;

M:=PureCubicalComplex(A);
fi;
   #########################

if IsHapPurePermutahedralComplex(MM) then
M:=PurePermutahedralComplex(A);
fi;

B:=BoundingPureComplex(M);
N:=HomotopyEquivalentMaximalPureSubcomplex(B,M);
N:=ContractedComplex(N);

while Size(N) < Size(M) do

   M:=CropPureComplex(N);
   B:=BoundingPureComplex(M);
   N:=HomotopyEquivalentMaximalPureSubcomplex(B,M);
   A:=N!.binaryArray;

if IsHapPureCubicalComplex(MM) then
   #########################
   for i in [2..Length(A)] do
      if A[i]=A[i-1] then A[i-1]:=0; fi;
   od;
   A:=Filtered(A,x-> not x=0);

   for d in [2..dim] do
   A:=PermuteArray(A,(1,d));
   for i in [2..Length(A)] do
      if A[i]=A[i-1] then A[i-1]:=0; fi;
   od;
   A:=Filtered(A,x-> not x=0);
   od;
   #########################

   N:=PureCubicalComplex(A);
fi;

if IsHapPurePermutahedralComplex(MM) then
   N:=PurePermutahedralComplex(A);
fi;

   ContractPureComplex(N);
if Length(arg)>1 then return CropPureComplex(N); fi;
od;

return M;
end);
##########################################
##########################################


######################################################################
######################################################################
InstallGlobalFunction(ZigZagContractedPureCubicalComplex,
function(T);
return ZigZagContractedPureComplex(T);
end);
######################################################################
######################################################################


############################################################
############################################################
InstallGlobalFunction(View3dPureComplex,
function(M)
local a1,a2,a3,A,  B, BB, squares, T, i, j, k, s, t,  VtoS, tmpdir, file,
    AppendTo, PrintTo;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

B:=M!.binaryArray;
A:=[];

for i in [1..Length(B)] do
for j in [1..Length(B[1])] do
for k in [1..Length(B[1][1])] do
if B[i][j][k]>0 then Add(A,[i,j,k]); fi;
od;od;od;

##############
VtoS:=function(V);
return Concatenation("(" , String(V[1]) , "," , String(V[2]) , "," , String(V[3]) , ")");
end;
##############

if IsHapPurePermutahedralComplex(M) then
################################
squares:=[];
squares[1]:=[ [1,0,2], [2,0,1], [2,-1,0], [1,-2,0], [0,-2,1], [0,-1,2] ];
squares[2]:=[ [0,1,2], [1,0,2], [0,-1,2], [-1,0,2] ];
squares[3]:=[ [-1,-2,0], [0,-2,-1], [1,-2,0], [0,-2,1] ];
squares[4]:=[ [-1,2,0], [0,2,-1], [0,1,-2], [-1,0,-2], [-2,0,-1], [-2,1,0] ];
squares[5]:=[ [0,-1,-2], [-1,0,-2], [-2,0,-1], [-2,-1,0], [-1,-2,0], [0,-2,-1]];
squares[6]:=[ [-2,-1,0], [-2,0,-1], [-2,1,0], [-2,0,1] ];
squares[7]:=[ [0,-1,2], [-1,0,2], [-2,0,1], [-2,-1,0], [-1,-2,0], [0,-2,1] ];
squares[8]:=[ [0,1,2], [0,2,1], [-1,2,0], [-2,1,0], [-2,0,1], [-1,0,2] ];
squares[9]:=[ [0,-1,-2], [-1,0,-2], [0,1,-2], [1,0,-2] ];
squares[10]:=[ [1,2,0], [0,2,1], [-1,2,0], [0,2,-1] ];
squares[11]:=[ [2,-1,0], [2,0,-1], [1,0,-2], [0,-1,-2], [0,-2,-1], [1,-2,0] ];
squares[12]:=[ [1,2,0], [2,1,0], [2,0,-1], [1,0,-2], [0,1,-2], [0,2,-1] ];
squares[13]:=[ [2,1,0], [2,0,1], [2,-1,0], [2,0,-1] ];
squares[14]:=[ [0,1,2], [1,0,2], [2,0,1], [2,1,0], [1,2,0], [0,2,1] ];
a1:=[1,1,-1];
a2:=[1,1,1];
a3:=[2,0,0];
T:=2*[a1,a2,a3];
######################################
fi;
if IsHapPureCubicalComplex(M) then
######################################
squares:=[];
squares[1]:=[ [0,0,0], [1,0,0], [1,1,0], [0,1,0] ];
squares[2]:=[ [0,0,1], [1,0,1], [1,1,1], [0,1,1] ];
squares[3]:=[ [0,0,0], [1,0,0], [1,0,1], [0,0,1] ];
squares[4]:=[ [0,1,0], [1,1,0], [1,1,1], [0,1,1] ];
squares[5]:=[ [0,0,0], [0,1,0], [0,1,1], [0,0,1] ];
squares[6]:=[ [1,0,0], [1,1,0], [1,1,1], [1,0,1] ];

T:=[[1,0,0],[0,1,0],[0,0,1]];;
######################################
fi;


tmpdir := DirectoryTemporary();;
file:=Filename( tmpdir , "tmp.asy" );

PrintTo(file, "import three;\n\n");
AppendTo(file, "size(500);\n\n");
AppendTo(file, "defaultpen(0.2);\n\n");

for i in [1..Length(A)] do
BB:=A[i][1]*T[1]+A[i][2]*T[2]+A[i][3]*T[3];

for j in [1..Length(squares)] do
s:=  BB+squares[j];
AppendTo(file,"path3[] g=");
AppendTo(file,VtoS(s[1]));
AppendTo(file,"--" );
AppendTo(file, VtoS(s[2]));
AppendTo(file,"--");
AppendTo(file,VtoS(s[3]));
AppendTo(file,"--");
AppendTo(file,VtoS(s[4]));
if Length(s)>4 then
AppendTo(file,"--");
AppendTo(file,VtoS(s[5]));
AppendTo(file,"--");
AppendTo(file,VtoS(s[6]));
fi;
AppendTo(file,"--cycle;\n");

AppendTo(file, "draw(surface(g),green+opacity(0.2));\n");
AppendTo(file, "draw(g,black);\n");

od;
od;

Exec( Concatenation( ASY_PATH, "-V ", file) );

#RemoveFile(file);
#file:=Filename(tmpdir,"");
#RemoveFile(file);

end);
#############################################################

#############################################################
#############################################################
InstallGlobalFunction(ViewPureComplex,
function(M);

if IsHapPureCubicalComplex(M) then
if Dimension(M)=3 then View3dPureComplex(M); fi;
if Dimension(M)=2 then ViewPureCubicalComplex(M); fi;
fi;

if IsHapPurePermutahedralComplex(M) then
if Dimension(M)=3 then View3dPureComplex(M); fi;
#if Dimension(M)=2 then ViewPurePermutahedralComplex(M); fi;
fi;

end);
#############################################################
#############################################################


#################################################################
#################################################################
InstallGlobalFunction(PathComponentOfPureComplex,
function(M,N)
local
        PathCompBinList,dim,dims, dimsSet,
        ArrayValueDim,ArrayValueDim1, ArrayAssignDim,
        ArrayIt, ArrayItBreak, revdimsSet,Fun,
        w,P,x,z,i,n;
n:=N+1;
dims:=EvaluateProperty(M,"arraySize");
revdimsSet:=List(dims,d->Reversed([2..d+1]));
dim:=Dimension(M);

ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim-1);
ArrayAssignDim:=ArrayAssignFunctions(dim);
ArrayIt:=ArrayIterate(dim);
ArrayItBreak:=ArrayIterateBreak(dim);

#############################################
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  #
PathCompBinList:=function()
local B,ColourNeighbours,ColourComponent,cart,CART,
NEWLYCOLOURED,GetStart,start,colour,ONE;

ONE:=List([1..Dimension(M)],i->1);
B:=StructuralCopy(FrameArray(M!.binaryArray));
#cart:=Cartesian(List([1..Dimension(M)],i->[-1,0,1]));
cart:=UnitBall(M);
RemoveSet(cart,List([1..Dimension(M)],i->0));

M!.pathReps:=[];

################################
ColourNeighbours:=function(x,j)
local w,y,z,bool;
bool:=false;
if ArrayValueDim(B,x)=j then
        for y in cart do
        z:=x+y;

        if ArrayValueDim(B,z)=1 then
        ArrayAssignDim(B,z,j);
        bool:=true;
        Add(NEWLYCOLOURED,z);
        fi;
        od;
fi;

return bool;
end;
################################

################################
ColourComponent:=function(j)
local bool,x,CPNC;
bool:=true;

while bool do
bool:=false;
CPNC:=ShallowCopy(NEWLYCOLOURED);
for x in CPNC do
if ColourNeighbours(x,j) then bool:=true; fi;;
Unbind(NEWLYCOLOURED[Position(NEWLYCOLOURED,x)]);
od;
NEWLYCOLOURED:=SSortedList(NEWLYCOLOURED);
od;

return bool;
end;
################################

################################
GetStart:=function()
local Fun,start,x;

start:=fail;

Fun:=function(x);
if ArrayValueDim(B,x)=1 then start:=x; return true; else return false;  fi;
end;

x:=ArrayItBreak(revdimsSet,Fun);

if not x=fail then revdimsSet[1]:=Reversed([2..x[1]]); fi;

if not start=fail then Add(M!.pathReps,x - ONE); fi;

return start;
end;
################################

colour:=1;
start:=GetStart();
while not start=fail do
colour:=colour+1;
ArrayAssignDim(B,start,colour);
NEWLYCOLOURED:=[start];
ColourComponent(colour);
start:=GetStart();


od;

M!.pathCompBinList:=UnframeArray(B);
M!.zeroBetti:=colour-1;
end;
##  ##  ##  ##  ##  ##  ##  ##  ##  ##  ##  #
#############################################

if not "pathCompBinList" in NamesOfComponents(M) then
PathCompBinList();
fi;

######################################
if N=0 then return M!.zeroBetti;fi;###
######################################

Fun:=function(z); if z=n then return 1;else return 0;fi;end;
P:=Array(M!.pathCompBinList,Fun);

if IsHapPureCubicalComplex(M) then
return Objectify(HapPureCubicalComplex,
           rec(
           binaryArray:=P,
           properties:=[
           ["dimension",Dimension(M)],
           ["arraySize",dims]]
           ));
fi;
if IsHapPurePermutahedralComplex(M) then
return Objectify(HapPurePermutahedralComplex,
           rec(
           binaryArray:=P,
           properties:=[
           ["dimension",Dimension(M)],
           ["arraySize",dims]]
           ));
fi;


end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(PathComponentOfPureCubicalComplex,
function(M,N);

return PathComponentOfPureComplex(M,N);

end);
#################################################################
#################################################################

