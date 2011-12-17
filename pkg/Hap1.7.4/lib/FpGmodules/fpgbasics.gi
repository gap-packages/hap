#(C) Graham Ellis, 2005-2007


#####################################################################
#####################################################################
InstallGlobalFunction(FpGCyclesModule,
function(arg)
local
	R,G,kk,Mdule, Action,
	PseudoBoundaryAsVec,
	WordToVectorList,
	prime, pp,
	BoundaryMatrix,
	MT,
	GactMat,
	one,
	zero,
	n;
R:=arg[1];
G:=R!.group;
kk:=arg[2];
prime:=PrimePGroup(G);
pp:=Order(G);
one:=Identity(GaloisField(prime));
zero:=0*one;
MT:=MultiplicationTable(R!.elts);

PseudoBoundaryAsVec:=[];

#####################################################################
WordToVectorList:=function(w,k)	#w is a word in R_k. 
local v,x,a;			#v is a list of vectors mod p.
v:=List([1..R!.dimension(k)],i->List([1..pp],j->0) );

for x in w do
a:=AbsoluteValue(x[1]);
v[a][x[2]]:=v[a][x[2]] + SignInt(x[1]);
od;

return v mod prime;
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
Action:=function(g,B);
return TransposedMat(GactMat(
Position(R!.elts,g),
TransposedMat(B)));
end;
#####################################################################

PseudoBoundaryAsVec:=[];
for  n in [1..R!.dimension(kk)] do
PseudoBoundaryAsVec[n]:=Flat(WordToVectorList(R!.boundary(kk,n),kk-1));
od;

#####################################################################
BoundaryMatrix:=function()	#Returns the matrix of d_k:R_k->R_k-1
local B,M,r, b, i, g,j,gB;		
				#M is actually the transpose of the matrix!
B:=TransposedMat(PseudoBoundaryAsVec);

M:=[];

for g in [1..pp] do
gB:=TransposedMat(GactMat(g,B));

for i in [0..R!.dimension(kk)-1] do
M[i*pp+g]:=gB[i+1];
od;
od;

M:=M*one;
ConvertToMatrixRep(M);
return (M);
end;
#####################################################################

Mdule:=SemiEchelonMat(NullspaceMat(BoundaryMatrix())).vectors;

#####################################################################


return Objectify(HapFPGModule,
	        rec(
		group:=G,
		matrix:=Mdule,
		action:=Action,
		dimension:=Length(Mdule),
		ambientDimension:=R!.dimension(kk)*pp
	));
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(RadicalOfFpGModule,
function(M)
local G, B, prime, g,radB, gens;

if M!.matrix=[] then return M; fi;

G:=M!.group;
prime:=PrimePGroup(G);
#B:=MutableCopyMat(M!.matrix);
B:=M!.matrix;
ConvertToMatrixRepNC(B,prime);

radB:=[];
gens:=GeneratorsOfGroup(G);
gens:=ReduceGenerators(gens,G);

for g in gens do
Append(radB,SemiEchelonMat(B-M!.action(g,B)).vectors);
od;

if Length(radB)>0 then
radB:=SemiEchelonMat(radB).vectors;
fi;

return Objectify(HapFPGModule,
	rec(
	group:=M!.group,
	matrix:=radB,
	action:=M!.action,
	dimension:=Length(radB),
	ambientDimension:=M!.ambientDimension));

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(GeneratorsOfFpGModule,
function(M)
local
	G,prime,pp,one,zero,rad,gens;

if "generators" in NamesOfComponents(M) then
return M!.generators;
fi;

G:=M!.group;
prime:=PrimePGroup(G);
pp:=Order(G);
one:=Identity(GaloisField(prime));
zero:=0*one;

rad:=RadicalOfFpGModule(M);
rad!.matrix:=Concatenation(rad!.matrix,ComplementaryBasis([M!.matrix,M]));

gens:=ComplementaryBasis([
rad!.matrix, M],
SemiEchelonMat(M!.matrix));

M!.generators:=gens;

return gens;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(ComplementaryBasis,
function(arg)
local 
	M,G,prime,pp,one,zero,B, NS, BC, heads,ln, i, v,zeroheads;

M:=arg[1][2];
G:=M!.group;
prime:=PrimePGroup(G);
pp:=Order(G);
one:=Identity(GaloisField(prime));
zero:=0*one;


B:=StructuralCopy(arg[1][1]);
if Length(arg)>1 then
NS:=StructuralCopy(arg[2]);
fi;

ConvertToMatrixRepNC(B,prime);
B:=MutableCopyMat(B);
heads:=SemiEchelonMatDestructive(B).heads;
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
end);
#####################################################################
#####################################################################
										#####################################################################
#####################################################################
InstallGlobalFunction(FpGModule,
function(A,G)
local
	M,GA,one,ambdim,pp,prime,MT,GactMat,Action,Elts,g;

pp:=Order(G);
prime:=PrimePGroup(G);
ambdim:=Length(A[1]);
Elts:=Elements(G);
MT:=MultiplicationTable(Elements(G));
one:=Identity(GF(prime));

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
Action:=function(g,B);
return TransposedMat(GactMat(
Position(Elts,g),
TransposedMat(B)));
end;
#####################################################################

GA:=[];

for g in G do
A:=MutableCopyMat(A);
Append(A,Action(g,A));
A:=SemiEchelonMat(A).vectors;
od;

M:=rec(
	group:=G,
	matrix:=A,
	dimension:=Length(A),
	ambientDimension:=ambdim,
	action:=Action
	);

return Objectify(HapFPGModule, M);
end);
#####################################################################
#####################################################################
										
										#####################################################################
#####################################################################
InstallGlobalFunction(DirectSumOfFpGModules,
function(arg)
local  M,N,A,B,G, zero, prime, v,w;

if Length(arg)=1 and IsList(arg) then

if Length(arg[1])=2 
then return DirectSumOfFpGModules(arg[1][1],arg[1][2]);
fi;

M:=arg[1][1];
N:=arg[1]{[2..Length(arg[1])]};
return 
DirectSumOfFpGModules(M,DirectSumOfFpGModules(N));
fi;

M:=arg[1];
N:=arg[2];
G:=M!.group;
if not G=N!.group then
Print("Modules must be over the same group ring. \n");
return fail; fi;

prime:=PrimePGroup(G);
zero:=Zero(GF(prime));

A:=[];
B:=[];

for w in M!.matrix do
v:=MutableCopyMat([w])[1];
while Length(v)<M!.ambientDimension + N!.ambientDimension do
Add(v,zero);
od;
Add(A,v);
od;

for w in N!.matrix do
v:=Reversed(MutableCopyMat([w])[1]);
while Length(v)<M!.ambientDimension + N!.ambientDimension do
Add(v,zero);
od;
v:=Reversed(v);
Add(B,v);
od;

Append(A,B);
A:=SemiEchelonMat(A).vectors;

return Objectify(HapFPGModule,

rec(
	group:=G,
	matrix:=A,
	dimension:=Length(A),
	ambientDimension:=M!.ambientDimension + N!.ambientDimension,
	action:=M!.action)
	);
	


end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(IsFpGModuleHomomorphismData,
function(M,N,F)
local
	R, Boundary, i, gens, vector, w,x, elts;

elts:=Elements(N!.group);
R:=ResolutionFpGModule(M,2);

if not elts=R!.elts then return fail; fi;

Boundary:=[];
for i in [1..R!.dimension(2)] do
Add(Boundary,R!.boundary(2,i));
od;

#gens:=GeneratorsOfFpGModule(N);

for w in Boundary do
vector:=N!.matrix[1]*Zero(GF(PrimePGroup(N!.group)));

for x in w do
vector:=vector+SignInt(x[1])*N!.action(elts[x[2]],[F[AbsInt(x[1])]])[1];
od;

if not IsZero(vector) then return false; fi;
od;

return true;
end);
#####################################################################
#####################################################################
