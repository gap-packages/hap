#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(ResolutionPrimePowerGroup,
function(arg)
local
	G,n,
	eltsG,
	gensG,
	Dimension,
	Boundary,
	Homotopy,
	PseudoBoundary,
	WordToVectorList,
	VectorListToWord,
	prime, pp,
	BoundaryMatrix,
	MT,
	GactZG,
	GactZGlist,
	ZGbasisOfKernel,
	one,
	InverseFlat,
	ComplementaryBasis,
	zero,
	pcgens,
	g,h,i,x,tmp;

G:=arg[1];
n:=arg[2];
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
if Length(arg)>2 then
pcgens:=ReduceGenerators(pcgens,G);
fi;
pcgens:=List(pcgens,x->Position(eltsG,x));

PseudoBoundary:=[];
for i in [1..n] do
PseudoBoundary[i]:=[];
od;

for x in gensG do
Append(PseudoBoundary[1], [  [[-1,1],[1,Position(eltsG,x)]]  ]);
od;

#####################################################################
Dimension:=function(i);
if i<0 then return 0; fi;
if i=0 then return 1; fi;
return Length(PseudoBoundary[i]);
end;
#####################################################################

#####################################################################
Boundary:=function(i,j);
if i<0 then return []; fi;
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
v:=List([1..Dimension(k)],x->List([1..pp],y->0) );

for x in w do
a:=AbsoluteValue(x[1]);
v[a][x[2]]:=v[a][x[2]] + SignInt(x[1]);
od;

return v mod prime;
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
BoundaryMatrix:=function(k)	#Returns the matrix of d_k:R_k->R_k-1
local M, b, i, j,v;		
				#M is actually the transpose of the matrix!

M:=[];

for i in [1..Dimension(k)] do
v:=WordToVectorList(Boundary(k,i),k-1);
for j in [1..pp] do
M[j + (i-1)*pp]:=Flat(GactZGlist(j,v));
od;
od;

return M*one;
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

#####################################################################
ComplementaryBasis:=function(arg)
local B, NS, BC, heads,ln, i, v;

B:=arg[1];
if Length(arg)>1 then
NS:=arg[2];
fi;

ConvertToMatrixRepNC(B,prime); 
heads:=SemiEchelonMat(B).heads;
ln:=Length(B[1]);
BC:=[];

if Length(arg)=1 then
	for i in [1..ln] do
	if heads[i]=0 then
	v:=List([1..ln], x->zero);
	v[i]:=one;
	Append(BC,[v]);
	fi;
	od;
else
 	for i in [1..ln] do
        if heads[i]=0 then
	v:=NS.vectors[NS.heads[i]];
	Append(BC,[v]);
	fi;
	od;
fi;

return BC;
end;
#####################################################################

#####################################################################
ZGbasisOfKernel:=function(k)		#The workhorse!
local  	i, v, g, h, b, ln, B, B1, B2,NS, 
	Bcomp, Bfrattini, BfrattComp, BasInts, IF,BndMat;

IF:=InverseFlat;
BndMat:=BoundaryMatrix(k);
ConvertToMatrixRepNC(BndMat,prime); 
NS:=SemiEchelonMat(NullspaceMat(BndMat));
B:=NS.vectors;
Bcomp:=ComplementaryBasis(B);

Bfrattini:=[];
for b in B do
for g in pcgens do     #I should check the maths here!
#for g in [2..pp] do
Append(Bfrattini,[b-  Flat(GactZGlist(g,IF(b)))]);
od;
od;

if  Length(arg)=1 and prime>2 then 
ln:=Length(B[1]);
for b in B do
for g in [2..pp] do
if Order(eltsG[g])=prime then
v:=List([1..ln],x->0);
for i in [1..prime] do
v:=v+ Flat(GactZGlist(Position(eltsG,eltsG[g]^i),IF(b)));
od;
Append(Bfrattini,[v]);
fi;
od;
od;
fi;

BfrattComp:=Concatenation(Bfrattini,Bcomp);

B1:=ComplementaryBasis(BfrattComp,NS);
#B2:=[];
#for b in B1 do
#i:=PositionProperty(b,x->(not IsZero(x)));
#Append(B2, [  B[NS.heads[i]] ]);
#od;

BasInts:=[];
for v in B1 do
Append(BasInts,[List(v,x->IntFFE(x))]);
od;

return List(BasInts,v->InverseFlat(v));
end;
#####################################################################

for i in [2..n] do
for x in ZGbasisOfKernel(i-1) do
#x:=TietzeReduction(PseudoBoundary[i],x);
Append(PseudoBoundary[i], [VectorListToWord(x)]   );
od;
od;


#####################################################################
Homotopy:=function(k,w)		#assume w is in kernel d_n
local v;	

v:=Flat(WordToVectorList(w,k));
v:=SolutionMat((BoundaryMatrix(k+1)),v*one);
Apply(v,i->IntFFE(i));

if not v=fail then v:=VectorListToWord(InverseFlat(v)); fi;
return v;
end;
#####################################################################

return Objectify(HapResolution,
	        rec(
		dimension:=Dimension,
		boundary:=Boundary,
		homotopy:=fail,
		partialHomotopy:=Homotopy,
		elts:=eltsG,
		group:=G,
		properties:=
			[["length",n],
			 ["reduced",true],
			 ["type","resolution"],
			 ["characteristic",prime],
			 ["isMinimal",true]]));
end);
#####################################################################
#####################################################################
