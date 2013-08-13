#(C) 2009 Graham Ellis

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
           binaryArray:=A,
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
UnitPermutahedralBall(Dimension(M)); fi;
return UnitPermutahedralBall(Dimension(M)); 	#This is so Fintan's 
						#package still works!
end);
#############################
#############################



##################################################
##################################################
InstallGlobalFunction(PureComplexToSimplicialComplex,
function(arg)
local M, DIM, AO,A,dim,dims,
      ArrayValueDim,
      CartProd,
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
CartProd:=Cartesian(List([1..dim],a->[2..dims[a]-1]));


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

for x in CartProd do
  if ArrayValueDim(AO,x)=1 then Vertices:=Vertices+1;
    y:=ArrayValueDim1(A,x{[2..dim]});
    y[x[1]]:=Vertices;
    VertexCoordinates[Vertices]:=x;
  fi;
od;

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
or
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

return Objectify(HapPureCubicalComplex,
        rec(
        binaryArray:=A,
        properties:=StructuralCopy(T!.properties)
        ));
fi;
#####################
#####################
if IsHapPurePermutahedralComplex(T) then
A:=HomotopyEquivalentSmallerSubPermArray(T!.binaryArray,S!.binaryArray);

return Objectify(HapPurePermutahedralComplex,
        rec(
        binaryArray:=A,
        properties:=StructuralCopy(T!.properties)
        ));
fi;
#####################

end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(HomotopyEquivalentMinimalPureCubicalSubcomplex,
function(T,S);
return HomotopyEquivalentMinimalPureCubicalSubcomplex(T,S);
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

return Objectify(HapPureCubicalComplex,
        rec(
        binaryArray:=A,
        properties:=StructuralCopy(T!.properties)
        ));
fi;

if IsHapPurePermutahedralComplex(T) then
A:=HomotopyEquivalentLargerSubPermArray(T!.binaryArray,S!.binaryArray);

return Objectify(HapPurePermutahedralComplex,
        rec(
        binaryArray:=A,
        properties:=StructuralCopy(T!.properties)
        ));
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

