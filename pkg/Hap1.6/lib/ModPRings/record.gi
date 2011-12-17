#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(ModPCohomologyRing,
function(arg)
local
	G,R,lngth,A,Degree,
	GeneratorLifts,
	Lift,
	WordToVectorList,
	prime,
	VectorListToWord,
	GactZG,
	GactZGlist,
	BoundaryMatrix,
	BoundaryMatrices,
	EchelonMatrices,
	MT,
	eltsG,
	pp,
	one,
	InverseFlat,
	IntPairList,
	IntToPair,
	PairToInt,
	ComposeGens,
	ComposeGensAugmented,
	SCT,
	dim,
	SolutionMatBoundaryMatrices,
	Echelonize,
	M,n,i,j,k;

if Length(arg)=1 then R:=arg[1]; 
G:=R!.group; fi;

if Length(arg)=2 then
G:=arg[1];
R:=ResolutionPrimePowerGroup(G,arg[2]); fi;

pp:=Order(G);

if not 
	(
	IsPrimePowerInt(pp)
	and
	EvaluateProperty(R,"isMinimal")=true 
	and
	EvaluateProperty(R,"characteristic")>0
	)
	then 
Print("The resolution must be minimal, for a prime-power group, and over a finite field. \n");
return fail;
fi;

prime:=Factors(pp)[1];
one:=Identity(GF(prime));

eltsG:=Elements(G);
MT:=MultiplicationTable(eltsG);

IntPairList:=[];
for n in [1..Length(R)] do
for i in [1..R!.dimension(n)] do
Append(IntPairList,[[n,i]]);
od;
od;
dim:=Length(IntPairList);

#####################################################################
IntToPair:=function(k);
return IntPairList[k]; 
end;
#####################################################################

#####################################################################
PairToInt:=function(x);
return Position(IntPairList,x);
end;
#####################################################################

#####################################################################
BoundaryMatrix:=function(k)     #Returns the matrix of d_k:R_k->R_k-1
local M, b, i, j,v;
                                #M is actually the transpose of the matrix!

M:=[];

for i in [1..R!.dimension(k)] do
v:=WordToVectorList(R!.boundary(k,i),k-1);
for j in [1..pp] do
M[j + (i-1)*pp]:=Flat(GactZGlist(j,v));
od;
od;

return M*one;
end;
#####################################################################


#####################################################################
WordToVectorList:=function(w,k) #w is a word in R_k.
local v,x,a;                    #v is a list of vectors mod p.
v:=List([1..R!.dimension(k)],x->List([1..pp],y->0) );

for x in w do
a:=AbsoluteValue(x[1]);
v[a][x[2]]:=v[a][x[2]] + SignInt(x[1]);
od;

return v mod prime ;
end;
#####################################################################


#####################################################################
VectorListToWord:=function(v)
local w, i, x;

w:=[];
for x in [1..Length(v)] do
for i in [1..Length(v[x])] do
if not v[x][i]=0 then
Append(w, MultiplyWord(v[x][i],[ [x,i]   ]));
fi;
od;
od;

return w;
end;
#####################################################################

#####################################################################
GactZG:=function(g,v)
local u,h;
u:=[];
for h in [1..Length(v)] do
u[MT[g][h]]:=v[h];
od;
return u;
end;
#####################################################################

#####################################################################
GactZGlist:=function(g,w)
local v,gw;
return List(w,v->GactZG(g,v));
end;
#####################################################################

#####################################################################
InverseFlat:=function(v)
local w,x,cnt,i;

w:=[];
cnt:=0;
while cnt< Length(v) do
x:=[];
for i in [cnt+1..cnt+pp] do
Append(x,[v[i]]);
cnt:=cnt+1;
od;
Append(w,[x]);
od;
return w;
end;
#####################################################################



BoundaryMatrices:=List([1..Length(R)],
i->BoundaryMatrix(i)*one);

for M in BoundaryMatrices do
ConvertToMatrixRepNC(M);
od;

EchelonMatrices:=[];

#####################################################################
Echelonize:=function()
local i,  Mt, T;

for i in [1..Length(R)] do
#We want to solve XM=W, so work on (Mt)(Xt)=(Wt)
# where 
Mt:=TransposedMat(BoundaryMatrices[i]);
ConvertToMatrixRepNC(Mt);
T:=SemiEchelonMatTransformation(Mt);
EchelonMatrices[i]:=[T.coeffs,T.heads,T.vectors];
od;

end;
#####################################################################
Echelonize();

#####################################################################
SolutionMatBoundaryMatrices:=function(m,w)
local h,u,v,i,cnt,pos,col,row,diff;

ConvertToVectorRep(w);
v:=EchelonMatrices[m][1]*w;
h:=StructuralCopy(Reversed(EchelonMatrices[m][2]));
u:=List([1..Length(h)],i->0*one);

while not IsZero(h) do

col:=PositionProperty(h,x->not x=0);
row:=h[col];
h[col]:=0;
diff:=EchelonMatrices[m][3][row]*u;
pos:=Length(u)+1-col;
u[pos]:=v[row]-diff;

od;

return u;
end;
#####################################################################

GeneratorLifts:=[1..Length(R)];
GeneratorLifts:=List(GeneratorLifts,n->[1..Dimension(R)(n)]);

	#The cohomology generators are represented by pairs [n,i] where n
	#lies between 1 and Length(R), and i lies between 1 and 
	#Dimension(R)(n).
	
	#GeneratorLifts[n,i] is a list [F[0],F[1],...,F[Length(R)-n]
	#where F[m] represents a ZG-homomorphism R_{n+m}-->R_{m}. Here F[m] is
	#actually a list [w_1,...,w_d] of words in R_m, with w_j the
	#image of the j-th generator of R_{n+m}. The ZG-homomorphism is 
	#induced by the standard cocycle representing the generator [n,i].

#####################################################################
Lift:=function(n,i) 	#This function calculates GeneratorLifts[n,i].
local
	Fm,j,w,m,s,x;

GeneratorLifts[n][i]:=[];

Fm:=[];  #m=0
for j in [1..Dimension(R)(n)] do
if j=i then Fm[j]:=[[1,1]];
else Fm[j]:=[];fi;
od;
GeneratorLifts[n][i][1]:=Fm;

for m in [1..Minimum(Length(R)-n,Int(Length(R)/2),n)] do
Fm:=[];

for j in [1..Dimension(R)(n+m)] do
w:=Flat(WordToVectorList([],m-1));
for x in R!.boundary(n+m,j) do
s:=MultiplyWord( SignInt(x[1]),
GeneratorLifts[n][i][m][AbsoluteValue(x[1])]);
s:=WordToVectorList(s,m-1);
s:=Flat(GactZGlist(x[2],s));
w:=w+s;
od;
#w:=SolutionMat(BoundaryMatrices[m],w*one);
w:=SolutionMatBoundaryMatrices(m,w*one);
Apply(w,i->IntFFE(i));
w:=VectorListToWord(InverseFlat(w));
#w:=R!.partialHomotopy(m-1,w);
Fm[j]:=w;
od;
GeneratorLifts[n][i][m+1]:=Fm;

od;
end;
#####################################################################

for n in [1..Length(R)] do
for i in [1..R!.dimension(n)] do
Lift(n,i);
od;
od;

#####################################################################
ComposeGens:=function(x,y)
local xx,yy,prdct,prod,j,z;


if y>x then prod:=ComposeGens(y,x);
if IsOddInt(IntToPair(x)[1]*IntToPair(y)[1]) then  #graded commutativity!!
for j in [1..Length(prod)/2] do
prod[2*j-1]:=-prod[2*j-1];
od;
fi;
return prod;
fi;

xx:=IntToPair(x);
yy:=IntToPair(y);

if xx[1]+yy[1]>Length(R) then return []; fi;



prdct:=List([1..R!.dimension(xx[1]+yy[1])],i->0);

for j in [1..R!.dimension(xx[1]+yy[1])] do
for z in GeneratorLifts[xx[1]][xx[2]][yy[1]+1][j] do
if AbsoluteValue(z[1])=yy[2] then 
prdct[j]:=prdct[j]+SignInt(z[1]);
fi;
od;
od;

prdct:= prdct mod prime;
prod:=[];

j:=PairToInt([xx[1]+yy[1],1])-1;
for z in [1..Length(prdct)] do
if not prdct[z]=0 then
Append(prod,[prdct[z],j+z]);
fi;
od;

return prod; 
end;
#####################################################################

#####################################################################
ComposeGensAugmented:=function(i,j) #ComposeGens and all other above
				    #functions do not incorporate 
local prod,k,yy;		    #H^0(G,Z_p). ComposeGensAugmented
				    #does so, and thus involves a 
				    #clumsy dimension shift.
if i>1 and j>1 then
prod:=ComposeGens(i-1,j-1);
for k in [1..Length(prod)/2] do
prod[2*k]:=prod[2*k]+1;
od;
return prod;
fi;

if i=1 then return [1,j]; fi;
if j=1 then return [1,i]; fi;
end;
#####################################################################

SCT:=EmptySCTable(dim+1,Zero(GF(prime)));

for i in [1..dim+1] do
for j in [1..dim+1] do
SetEntrySCTable(SCT,i,j,ComposeGensAugmented(i,j));
od;
od;

A:=AlgebraByStructureConstants(GF(prime),SCT);


#####################################################################
Degree:=function(x)
local i;

i:=Position(GeneratorsOfAlgebra(A),x);

if i=1 then return 0; fi;

return IntToPair(i-1)[1];
end;
#####################################################################

A!.degree:=Degree;

return A;
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(ModPRingGenerators,
function(A)
local S, gens, gensA, x, dim, dimgensA;

S:=GeneratorsOfAlgebra(A);
dim:=Dimension(A);
gens:=[S[1]];
gensA:=Subalgebra(A,gens);
dimgensA:= Dimension(gensA);

for x in S do

if not x in gensA then
Append(gens,[x]);
gensA:=Subalgebra(A,gens);
dimgensA:=Dimension(gensA);
fi;

if dimgensA=dim then return gens; fi;

od;


end);
#####################################################################
#####################################################################
