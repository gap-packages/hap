
#####################################################################
#####################################################################
InstallGlobalFunction(PurePermutahedralComplex,
function(AA)
# Inputs a binary array and returns a pure permutahedral complex.
local
        A;

A:=StructuralCopy(AA);

return Objectify(HapPurePermutahedralComplex,
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
InstallOtherMethod(Dimension,
"Dimension of pure permutahedral complex",
[IsHapPurePermutahedralComplex],
function(f) return EvaluateProperty(f,"dimension");
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Size,
"Volume of a pure permutahedral complex",
[IsHapPurePermutahedralComplex],
function(f) return Sum(Flat(f!.binaryArray));
end);
#####################################################################
#####################################################################

##################################################
##################################################
InstallGlobalFunction(PermutahedralComplexToRegularCWComplex, 
function(M)
local AO,A,MO, dim,dims, DIM,
      ArrayValueDim, FN,
      #CartProd, 
      dimSet, ArrayIt, 
      Vertices, MVertices, VertexCoordinates,ArrayValueDim1,
      IsMSimplex,
      Ball, Balls,
      SimplicesLst,  EnumeratedSimplex,
      bnd, pos, orien, TMP, tmp, cnt,
       b, i, j, t, t1, t2, d, dd, v, x, y, Y;

#################################
if not IsHapPurePermutahedralComplex(M) then
Print("This function must be applied to a pure permutahedral complex.\n");
return fail; fi;
#################################

MO:=FrameArray(M!.binaryArray)*1;
MO:=FrameArray(MO);
AO:=FrameArray(M!.binaryArray)*1;
AO:=FrameArray(AO);
AO:=PurePermutahedralComplex(AO);
AO:=ThickenedPureComplex(AO);
AO:=AO!.binaryArray;
A:=AO*1;
dim:=ArrayDimension(A);
DIM:=dim;
dims:=ArrayDimensions(A);
Vertices:=0;
VertexCoordinates:=[];
MVertices:=[];
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

#######################
IsMSimplex:=function(x);
if true in List(x, i->MVertices[i])
then return true;
else return false;
fi;
end;
########################

#for x in CartProd do
FN:=function(x)
local y;
  if ArrayValueDim(AO,x)=1 then Vertices:=Vertices+1;
    y:=ArrayValueDim1(A,x{[2..dim]});
    y[x[1]]:=Vertices;
    VertexCoordinates[Vertices]:=x;
    if ArrayValueDim(MO,x)=1 then
    MVertices[Vertices]:=true;
    else
    MVertices[Vertices]:=false;
    fi;
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
  SimplicesLst[1]:=Filtered(SimplicesLst[1],IsMSimplex);
  SimplicesLst[1]:=SSortedList(SimplicesLst[1]);
fi;

if DIM>=1 then
for v in Vertices do
  x:=VertexCoordinates[v];
  for b in Ball do
    t:= ArrayValueDim(A,b+x);
    if t>v then Add(SimplicesLst[2],[v,t]); fi;
  od;
od;
  SimplicesLst[2]:=Filtered(SimplicesLst[2],IsMSimplex);
SimplicesLst[2]:=SSortedList(SimplicesLst[2]);
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
SimplicesLst[j+1]:=Filtered(SimplicesLst[j+1],IsMSimplex);
SimplicesLst[j+1]:=SSortedList(SimplicesLst[j+1]);

od;
fi;
#############################################
EnumeratedSimplex:=function(v)
local pos;
pos:=Position(TMP[v[1]],v);
if pos=fail then return fail; else
return tmp[v[1]][ pos ];
fi;
end;
#############################################

bnd:=[];
orien:=[];

bnd[1]:=List(SimplicesLst[dim+1],x->[1,0]);
orien[1]:=List(SimplicesLst[dim+1],x->[1]);

for d in [2..dim+1] do
bnd[d]:=List(SimplicesLst[dim+2-d],i->[]);
orien[d]:=List(SimplicesLst[dim+2-d],i->[]);
TMP:=List(Vertices,i->[]);;
tmp:=List(Vertices,i->[]);;;
cnt:=1;
for x in SimplicesLst[dim+2-d] do
Add(TMP[x[1]],x);
Add(tmp[x[1]],cnt);
cnt:=cnt+1;
od;

for i in [1..Length(SimplicesLst[dim+3-d])] do
x:=SimplicesLst[dim+3-d][i];
for j in [1..Length(x)] do
b:=Filtered(x,a->not a=x[j]);
pos:=EnumeratedSimplex(b);
if not pos=fail then
Add(bnd[d][pos],i);
  if IsOddInt(j) then Add(orien[d][pos],1);
  else Add(orien[d][pos],-1);
  fi;
fi;
od;
od;

Apply(bnd[d],x->Concatenation([Length(x)],x));
od;

Add(bnd,[]);
Add(orien,[]);

Y:=RegularCWComplex(bnd,orien);
OrientRegularCWComplex(Y);
HAP_Sequence2Boundaries(Y);

return Y;

end);
##################################################

#####################################################################
#####################################################################
InstallMethod(ChainComplex,
"Cellular chain complex of a pure permutahedral complex",
[IsHapPurePermutahedralComplex],
function(M)
local Y;
Y:=PermutahedralComplexToRegularCWComplex(M);
return ChainComplex(Y);;
end);
#####################################################################
#####################################################################

#####################################################
#####################################################
InstallGlobalFunction(CubicalToPermutahedralArray,
function(D)
local sq2, sq3, sq6, A, B, v, w, P, LN;

sq2:=Float("1.414213562":PrecisionFloat:=1000000);
sq3:=Float("1.732050808":PrecisionFloat:=1000000);
sq6:=sq2*sq3;

if Length(D)=0 then return D; fi;
LN:=Length(D[1]);

if LN=2 then
A:=[[Rat(1/sq2), -Rat(1/sq2),0], [Rat(1/sq6),Rat(1/sq6),-Rat(2/sq6)]];
B:=[[1,-1,0],[0,-1,1],[0,0,1]]^-1;
fi;

if LN=3 then
A:=[[Rat(1/sq2), -Rat(1/sq2),0,0],
    [1/2, 1/2,-1/2, -1/2],
    [0,0,Rat(1/sq2), -Rat(1/sq2)]];
B:=[[1,-1,0,0],[0,-1,1,0],[0,0,1,-1],[0,0,0,1]]^-1;
fi;


P:=[];

for v in D do
w:=v*A;
w:=w*B;
Add(P, w{[1..LN]});
od;

return P;
end);
####################################################
####################################################


#####################################################
#####################################################
InstallGlobalFunction(PermutahedralToCubicalArray,
function(D)
local sq2, sq3, sq6, A, B, BAS, v, w, P, LN;

sq2:=Float("1.414213562":PrecisionFloat:=1000000);
sq3:=Float("1.732050808":PrecisionFloat:=1000000);
sq6:=sq2*sq3;

if Length(D)=0 then return D; fi;
LN:=Length(D[1]);

if LN=2 then
A:=[[Rat(1/sq2), -Rat(1/sq2),0],
    [Rat(1/sq6),Rat(1/sq6),-Rat(2/sq6)],
    [0,0,1]]^-1;

BAS:=[[1,-1,0],[0,-1,1]];
fi;

if LN=3 then
A:=[[Rat(1/sq2), -Rat(1/sq2),0,0],
    [1/2, 1/2,-1/2, -1/2],
    [0,0,Rat(1/sq2), -Rat(1/sq2)],
    [0,0,0,1]]^-1;
BAS:=[[1,-1,0,0], [1,1,-1,-1], [0,0,1,-1]];

fi;

P:=[];

for v in D do
w:=v*BAS;
w:=w*A;
w:=w{[1..LN]};
Add(P,w);
od;

return P;
end);
####################################################
####################################################

