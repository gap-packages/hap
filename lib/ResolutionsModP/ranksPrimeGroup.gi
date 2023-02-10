#(C) Graham Ellis, 2005-2006


#####################################################################
#####################################################################
InstallGlobalFunction(RankPrimeHomology,
function(arg)
local
	G,n,
	eltsG,
	gensG,
	Dimension,
	DimList,
	Boundary,
	PseudoBoundary,
	PseudoBoundaryAsVec,
	WordToVectorList,
	VectorListToWord,
	prime, pp,
	BoundaryMatrix,
	MT,
	GactMat,
	ZGbasisOfKernel,
	one,
	InverseFlat,
	ComplementaryBasis,
	zero,
	pcgens,
	g,h,i,x,xx,xxx,tmp;


G:=arg[1];
if arg[2]=-1 then n:=1000;
else
n:=arg[2];
fi;

tmp:=SSortedList(Factors(Order(G)));
if Length(tmp)>1 then 
Print("This function can only be applied to small prime-power groups. \n");
return fail;
fi;
prime:=tmp[1];
pp:=Order(G);
one:=Identity(GaloisField(prime));
zero:=0*one;

gensG:=ReduceGenerators(GeneratorsOfGroup(G),G);
eltsG:=Elements(G);
MT:=MultiplicationTable(eltsG);

pcgens:=Pcgs(G);
pcgens:=ReduceGenerators(pcgens,G);
pcgens:=List(pcgens,x->Position(eltsG,x));

DimList:=[];
PseudoBoundary:=[];
PseudoBoundaryAsVec:=[];
for i in [1..n] do
PseudoBoundary[i]:=[];
PseudoBoundaryAsVec[i]:=[];
od;

#####################################################################
Dimension:=function(i);
if i<0 then return 0; fi;
if i=0 then return 1; fi;
return DimList[i];
end;
#####################################################################

#####################################################################
Boundary:=function(i,j);
if i<=0 then return []; fi;
if j>0 then
return PseudoBoundary[i][j]; 
else return
 NegateWord(PseudoBoundary[i][-j]);
 fi;
end;
#####################################################################

#####################################################################
WordToVectorList:=function(w,k)	#w is a word in R_k. 
local v,x,a;			#v is a list of vectors mod p.
#v:=List([1..Dimension(k)],x->List([1..pp],y->0) );
v:=ListWithIdenticalEntries(Dimension(k),ListWithIdenticalEntries(pp,0) );


for x in w do
a:=AbsoluteValue(x[1]);
v[a][x[2]]:=v[a][x[2]] + SignInt(x[1]);
od;

return v mod prime;
end;
#####################################################################

for x in gensG do
Add(PseudoBoundary[1], [[-1,1],[1,Position(eltsG,x)]]  );
Add(PseudoBoundaryAsVec[1], Flat(WordToVectorList
( [[-1,1],[1,Position(eltsG,x)]] ,0))*one);
od;
DimList[1]:=Length(PseudoBoundary[1]);

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
GactMat:=function(g,tB)
local k,q,h,C;

C:=[];

k:=0;
for q in [0..(-1+Length(tB)/pp)] do

for h in [1..pp] do
C[k+MT[g][h]]:=tB[k+h];
od;
k:=k+pp;
od;

ConvertToMatrixRepNC(C);

return C;
end;
#####################################################################



#####################################################################
BoundaryMatrix:=function(k)	#Returns the matrix of d_k:R_k->R_k-1
local B,M,r, b, i, g,j,gB;		
				#M is actually the transpose of the matrix!
B:=TransposedMat(PseudoBoundaryAsVec[k]);

M:=[];

for g in [1..pp] do
gB:=TransposedMat(GactMat(g,B));

for i in [0..Dimension(k)-1] do
M[i*pp+g]:=gB[i+1];
od;
od;


return M;
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
Add(x,v[i]);
cnt:=cnt+1;
od;
Add(w,x);
od;
return w;
end;
#####################################################################

#####################################################################
ComplementaryBasis:=function(arg)
local B, NS, BC, heads,ln, i, v,zeroheads;


B:=arg[1];
if Length(arg)>1 then
NS:=arg[2];
fi;

ConvertToMatrixRepNC(B,prime); 
heads:=SemiEchelonMat(B).heads;
ln:=Length(B[1]);
BC:=[];

zeroheads:=Filtered([1..ln],i->heads[i]=0);

if Length(arg)=1 then
	for i in zeroheads do
	v:=ListWithIdenticalEntries(ln,zero);
	v[i]:=one;
	ConvertToVectorRep(v,prime);
	Add(BC,v);
	od;
else
	for i in zeroheads do
	v:=NS.vectors[NS.heads[i]];
	ConvertToVectorRep(v,prime);
	Add(BC,v);
	od;
fi;



return BC;
end;
#####################################################################

#####################################################################
ZGbasisOfKernel:=function(k)		#The workhorse!
local  	tB,i, v, g, h, b,bb, ln, B, B1, B2,NS, 
	Bcomp, BndMat;

BndMat:=BoundaryMatrix(k);
PseudoBoundaryAsVec[k]:=0;
PseudoBoundary[k]:=0;

ConvertToMatrixRepNC(BndMat,prime); 
NS:=SemiEchelonMat(NullspaceMat(BndMat));

tB:=TransposedMat(NS.vectors);
Bcomp:=ComplementaryBasis(NS.vectors);


for g in pcgens do     	 
Append(Bcomp,SemiEchelonMat(TransposedMat(tB-GactMat(g,tB))).vectors);
od;							



B1:=ComplementaryBasis(Bcomp,NS);
NS:=0;
#B1:=List(B1,v->List(v,x->IntFFE(x)));
Apply(B1,v->List(v,x->IntFFE(x)));
return B1;
end;
#####################################################################

for i in [2..n] do


for x in ZGbasisOfKernel(i-1) do
Add(PseudoBoundary[i], VectorListToWord(InverseFlat(x))   );
Add(PseudoBoundaryAsVec[i], x*one   );
od;
DimList[i]:=Length(PseudoBoundary[i]);


if i>15 and arg[2]=-1 then
x:=PoincareSeries(List([0..i],j->Dimension(j)),i+1);
xx:=PoincareSeries(List([0..i-1],j->Dimension(j)),i);
xxx:=PoincareSeries(List([0..i-2],j->Dimension(j)),i-1);
if x=xx and xx=xxx and (not x=fail) then return x; fi;
fi;

od;

#####################################################################
#####################################################################

return Dimension ;
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(EfficientNormalSubgroups,
function(arg)
local G, dim, EffSubgroups, N, x, R,S,T,D;

###########################################################
G:=arg[1];
if Length(arg)>1 then dim:=arg[2]; else dim:=4; fi;
###########################################################


###########################################################
if not IsPrimePowerInt(Order(G)) then
Print("This function applies to prime-power groups only.");
return fail;
fi;
###########################################################


N:=NormalSubgroups(G);
N:=Filtered(N,x->(Order(x)>1 and Order(x)<Order(G)));

EffSubgroups:=[];
R:=RankPrimeHomology(G,dim);;
R:=List([0..dim],i->R(i));

for x in N do
S:=ResolutionPrimePowerGroup(x,dim);
T:=ResolutionPrimePowerGroup(G/x,dim);
D:=ResolutionDirectProduct(S,T);
D:=List([0..dim],i->D!.dimension(i));

if D=R then Add(EffSubgroups,x); fi;
od;

return EffSubgroups;
end);
#####################################################################
#####################################################################
