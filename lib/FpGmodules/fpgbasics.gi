#(C) Graham Ellis, 2005-2007


#####################################################################
#####################################################################
InstallGlobalFunction(DesuspensionFpGModule,
function(arg)
local
	R,G,kk,KK,Mdule, Action,
	PseudoBoundaryAsVec,
	WordToVectorList,
	prime, pp,
	BoundaryMatrix,
	MT,
	GactMat,
	one,
	zero,
	n;

if IsHapFPGModule(arg[1]) then
R:=ResolutionFpGModule(arg[1],arg[2]);
else
R:=arg[1];
fi;
G:=R!.group;
kk:=arg[2];
KK:=Maximum(kk,1);
prime:=EvaluateProperty(R,"characteristic");
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
for  n in [1..R!.dimension(KK)] do
PseudoBoundaryAsVec[n]:=Flat(WordToVectorList(R!.boundary(KK,n),KK-1));
od;

#####################################################################
BoundaryMatrix:=function()	#Returns the matrix of d_k:R_k->R_k-1
local B,M,r, b, i, g,j,gB;		
				#M is actually the transpose of the matrix!
B:=TransposedMat(PseudoBoundaryAsVec);

M:=[];

for g in [1..pp] do
gB:=TransposedMat(GactMat(g,B));

for i in [0..R!.dimension(KK)-1] do
M[i*pp+g]:=gB[i+1];
od;
od;

M:=M*one;
ConvertToMatrixRep(M);
return (M);
end;
#####################################################################

if kk=0 then
Mdule:=SemiEchelonMat(BoundaryMatrix()).vectors;
else
Mdule:=SemiEchelonMat(NullspaceMat(BoundaryMatrix())).vectors;
fi;

#####################################################################


return Objectify(HapFPGModule,
	        rec(
		group:=G,
		matrix:=Mdule,
		action:=Action,
		dimension:=Length(Mdule),
		ambientDimension:=R!.dimension(kk)*pp,
		characteristic:=prime
	));
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(RadicalOfFpGModule,
function(M)
local  G, B, prime, g,radB, gens;

if M!.matrix=[] then return M; fi;

G:=M!.group;
prime:=M!.characteristic;

B:=M!.matrix;
ConvertToMatrixRepNC(B,prime);

radB:=[];
gens:=GeneratorsOfGroup(SylowSubgroup(G,prime));
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
	ambientDimension:=M!.ambientDimension,
	characteristic:=prime));

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(MultipleOfFpGModule,
function(wrd,M)
local G, B, BB,prime, g, gens;

if M!.matrix=[] then return M; fi;

G:=M!.group;
prime:=M!.characteristic;
B:=MutableCopyMat(M!.matrix);
BB:=B*0;
ConvertToMatrixRepNC(B,prime);

for g in wrd do
BB:=BB+ M!.action(g,B);
od;

if Length(BB)>0 then
BB:=SemiEchelonMat(BB).vectors;
fi;

return BB;

end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(GeneratorsOfFpGModule,
function(M)
local
	G,prime,pp,one,zero,rad,gens,  v, NS, rrb,rrb1, Reducedgens;

if "generators" in NamesOfComponents(M) then
return M!.generators;
fi;

G:=M!.group;

prime:=M!.characteristic; 
pp:=Order(G);
one:=Identity(GaloisField(prime));
zero:=0*one;

rad:=RadicalOfFpGModule(M);
rad!.matrix:=Concatenation(rad!.matrix,ComplementaryBasis([M!.matrix,M]));

gens:=ComplementaryBasis([ rad!.matrix, M],
SemiEchelonMat(M!.matrix));

NS:=Length(M!.matrix);

if not IsPrimePowerInt(Order(G)) then
        Reducedgens:=[];
        rrb:=0;
        rrb1:=0;
         for v in gens do
                if rrb1=NS then break; fi;
                rrb1:=FpGModule
                (Concatenation(Reducedgens,[v]),G,prime)!.dimension;
                if rrb1>rrb then
                rrb:=rrb1;
                Reducedgens:=Concatenation(Reducedgens,[v]);
                fi;
        od;

gens:=StructuralCopy(Reducedgens);
                for v in gens do
	        Reducedgens:=Filtered(Reducedgens,x->not x=v);
	        rrb:=FpGModule
	        (Reducedgens,G,prime)!.dimension;
	        if rrb<NS then Add(Reducedgens,v); fi;
	        od;

gens:=Reducedgens;
fi;

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

pp:=Order(G);

prime:=M!.characteristic; 

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
function(arg)
local
	A,G,M,GA,one,ambdim,pp,prime,MT,GactMat,Action,Elts,g;

A:=arg[1];
G:=arg[2];
pp:=Order(G);
if Length(arg)>2 then prime:=arg[3];
else 
if IsPrimePowerInt(Order(G)) then prime:=PrimePGroup(G); 
else Print("A prime must be entered for non-prime-power groups.\n"); 
return fail; fi;
fi;
if Length(A)>0 then
ambdim:=Length(A[1]);
else ambdim:=0;fi;
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

if Length(A)>0 then
for g in G do
A:=MutableCopyMat(A);
Append(A,Action(g,A));
A:=SemiEchelonMat(A).vectors;
od;
fi;

M:=rec(
	group:=G,
	matrix:=A,
	dimension:=Length(A),
	ambientDimension:=ambdim,
	action:=Action,
	characteristic:=prime
	);

return Objectify(HapFPGModule, M);
end);
#####################################################################
#####################################################################
										
										#####################################################################
#####################################################################
InstallGlobalFunction(FpGModuleHomomorphismNC,
function(M,N,A);

return Objectify(HapFPGModuleHomomorphism,

rec(	source:=M,
	target:=N,
	matrix:=A )
);

end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(FpGModuleHomomorphism,
function(M,N,A);

if not IsFpGModuleHomomorphismData(M,N,A) then
return fail;
else
return FpGModuleHomomorphismNC(M,N,A); fi;
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
prime:=M!.characteristic;

if not (G=N!.group and prime=N!.characteristic) then
Print("Modules must be over the same group ring. \n");
return fail; fi;


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
	action:=M!.action,
	characteristic:=prime)
	);
	
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(IsFpGModuleHomomorphismData,
function(M,N,F)
local
	R, Boundary, i,  vector, w,x, elts;

elts:=Elements(N!.group);
R:=ResolutionFpGModule(M,2);

if not (elts=R!.elts and M!.characteristic = N!.characteristic) then 
Print("The two modules are over different group rings.\n");
return fail; fi;

Boundary:=[];
for i in [1..R!.dimension(2)] do
Add(Boundary,R!.boundary(2,i));
od;

for w in Boundary do
vector:=N!.matrix[1]*Zero(GF(N!.characteristic));

for x in w do
vector:=vector+SignInt(x[1])*N!.action(elts[x[2]],[F[AbsInt(x[1])]])[1];
od;

if not IsZero(vector) then return false; fi;
od;

return true;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(IntersectionOfFpGModules,
function(M,N)
local V;

if not (M!.group = N!.group and M!.characteristic=N!.characteristic)
or not M!.ambientDimension = N!.ambientDimension
then 
Print("Modules are not compatible. \n"); return fail;
fi;

V:=SumIntersectionMat(M!.matrix,N!.matrix)[2];


return Objectify(HapFPGModule,

rec(
        group:=M!.group,
        matrix:=V,
        dimension:=Length(V),
        ambientDimension:=M!.ambientDimension,
        action:=M!.action,
	characteristic:=M!.characteristic)
        );

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(SumOfFpGModules,
function(M,N)
local V;

if not (M!.group = N!.group and M!.characteristic=N!.characteristic)
or not M!.ambientDimension = N!.ambientDimension
then
Print("Modules are not compatible. \n"); return fail;
fi;

V:=SumIntersectionMat(M!.matrix,N!.matrix)[1];


return Objectify(HapFPGModule,

rec(
        group:=M!.group,
        matrix:=V,
        dimension:=Length(V),
        ambientDimension:=M!.ambientDimension,
        action:=M!.action,
	characteristic:=M!.characteristic)
    );

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(VectorsToFpGModuleWords,
function(M,v)
local				#v is a list of vectors
	G,gens, Ggens, 
	VecToWord,g,w;

G:=M!.group;
gens:=GeneratorsOfFpGModule(M);

Ggens:=[];
for g in G do
Append(Ggens, M!.action(g,gens));
od;
ConvertToMatrixRep(Ggens);

########################################
VecToWord:=function(w)
local wrd,i,g ;

wrd:=[];

for i in [0..Length(gens)-1] do
for g in [1..Order(G)] do
if not IsZero(w[i*Order(G)+g]) then
Append(wrd,
MultiplyWord(IntFFE(w[i*Order(G)+g]),[[i+1,g]])
);fi;
od;od;

return wrd;
end;
########################################

w:=SolutionsMatDestructive(MutableCopyMat(Ggens),MutableCopyMat(v));
Apply(w,y->VecToWord(y));
return w;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(\+,
"for FpG homomorphisms",
[IsHapFPGModuleHomomorphism,IsHapFPGModuleHomomorphism],
function(x,y)
if not (Source(x)=Source(y) and Target(x)=Target(y)) then return fail; fi;
return FpGModuleHomomorphismNC(x!.source,x!.target,x!.matrix + x!.matrix);
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(SumOp,
"for FpG homomorphisms",
[IsHapFPGModuleHomomorphism,IsHapFPGModuleHomomorphism],
function(x,y)
if not (Source(x)=Source(y) and Target(x)=Target(y)) then return fail; fi;
return FpGModuleHomomorphismNC(x!.source,x!.target,x!.matrix + x!.matrix);
end);
#####################################################################
#####################################################################



#####################################################################
#####################################################################
InstallOtherMethod(Source,
"for FpG homomorphisms",
[IsHapFPGModuleHomomorphism],
function(x)
return x!.source;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Target,
"for FpG homomorphisms",
[IsHapFPGModuleHomomorphism],
function(x)
return x!.target;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(ImageOfFpGModuleHomomorphism,
function(x)
return FpGModule(x!.matrix,Target(x)!.group);
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Rank,
"for FpG homomorphisms",
[IsHapFPGModuleHomomorphism],
function(x)
return Dimension(FpGModule(x!.matrix,Target(x)!.group));
end);
#####################################################################
#####################################################################


