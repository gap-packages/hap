#(C) 2007 Graham Ellis

#####################################################################
#####################################################################
InstallGlobalFunction(FpGModuleDualBasis,
		function(M)
		local 	
		R, G, Boundary, SortedBoundary,
		one, zero, gens, EltsG, VecToMat,
		MatRep,b,i, BoundToPmat, B,V;

R:=ResolutionFpGModule(M,2);
Boundary:=List( [1..R!.dimension(2)], i-> R!.boundary(2,i));


G:=M!.group;
EltsG:=Elements(G);
gens:=GeneratorsOfFpGModule(M);
one:=Identity(GF(M!.characteristic));
zero:=0*one;

#################################################
VecToMat:=function(v)
local m,i;


m:=[];
for i in [1..Length(v)/Order(G)] do
m[i]:=v{[(i-1)*Order(G)+1..i*Order(G)]};
od;

return m;
end;
#################################################
#################################################
MatRep:=function(g)
local A;
A:=IdentityMat(Order(G));
return (M!.action(g,A)*one);
end;
#################################################

SortedBoundary:=[];
for b in Boundary do
Append(SortedBoundary,[
List([1..R!.dimension(1)], i->List(Filtered(b,x->x[1]=i),y->y[2]))
]);
od;

Unbind(Boundary);

########################################################
BoundToPmat:=function(b)
local PMat, GMat,w,g;

PMat:=[];

for w in [1..R!.dimension(1)] do
GMat:=NullMat(Order(G),Order(G),GF(M!.characteristic));

for g in b[w] do
GMat:=GMat+MatRep(EltsG[g]);
od;

Append(PMat,GMat);
od;

ConvertToMatrixRep(PMat);

return PMat;
end;
########################################################

b:=SortedBoundary[1];
B:=NullspaceMat(BoundToPmat(b));

for b in SortedBoundary{[2..R!.dimension(2)]} do

B:=SumIntersectionMat(B,NullspaceMat(BoundToPmat(b)))[2];

od;

B:=List(B,b->VecToMat(b));

return 

rec(
freeModule:=FpGModule(IdentityMat(Order(G),GF(M!.characteristic)),G),

basis:=B);
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(ProjectedFpGModule,
function(M,k)
local
	SE,Hds,A,s,t,G;

G:=M!.group;
s:=1+(k-1)*Order(G);
t:=k*Order(G);


ConvertToMatrixRep(M!.matrix);
SE:=SemiEchelonMat(M!.matrix);
Hds:=Filtered(SE.heads,i->(i>s));


A:=List(M!.matrix{Hds}, v -> v{[s..t]});
ConvertToMatrixRep(A);
if not Length(A) =0 then
A:=SemiEchelonMat(A).vectors;
fi;

return Objectify(HapFPGModule,
                rec(
                group:=G,
                matrix:=A,
                action:=M!.action,
                dimension:=Length(A),
                ambientDimension:=Order(G),
		characteristic:=M!.characteristic
	        ));

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(RandomHomomorphismOfFpGModules,
function(M,N)
local 
	G,DualBasisM, Map, ProjN, GFP, VecToMat, PreImVec,
	PreimageVector, SemiEchN, Mapp,
	ExpandVector,zero, RandomVector,n,l,dd,j,B,A;

Map:=[];
G:=M!.group;
GFP:=GF(M!.characteristic);
zero:=Zero(GFP);
if not G=N!.group then
Print("The modules must be over a common group.\n");
return fail;
fi;

DualBasisM:=FpGModuleDualBasis(M).basis;
dd:=Length(DualBasisM[1]);
Apply(DualBasisM,m->Flat(m));
l:=Length(DualBasisM[1]);
ConvertToMatrixRep(DualBasisM);

ConvertToMatrixRep(N!.matrix);
SemiEchN:=SemiEchelonMat(N!.matrix);

#################################################
VecToMat:=function(v)
local m,i;


m:=[];
for i in [1..Length(v)/Order(G)] do
m[i]:=v{[(i-1)*Order(G)+1..i*Order(G)]};
od;

return m;
end;
#################################################

#############################################################
RandomVector:=function(A)
local v,i;

if Length(A)=0 then
return ListWithIdenticalEntries(l,zero);
fi;

v:=Random(GFP)*A[1];

for i in [2..Length(A)] do
v:=v+Random(GFP)*A[i];
od;

return v;
end;
#############################################################


#############################################################
ExpandVector:=function(v,k)
local s,t;
s:=(k-1)*Order(G); t:=l-k*Order(G);
return Concatenation(
ListWithIdenticalEntries(s,zero),
v,
ListWithIdenticalEntries(t ,zero));
end;
#############################################################

#############################################################
PreimageVector:=function(v,k)
local pv, vv,B, Hds,frst,i;
	#Inputs a "suitable" vector v and returns a vector
	#pv which: (1) lies in the module N, and (2) is such that the vectors
	#pv and v agree in their first k entries.

	#This function is stupidly clumsy. Easy modifications will
	#speed it up.
	
B:=SemiEchN.vectors;
Hds:=SemiEchN.heads;

vv:=Concatenation(v,ListWithIdenticalEntries(N!.ambientDimension-Length(v),zero));
pv:=ListWithIdenticalEntries(N!.ambientDimension,zero);

frst:=PositionProperty(vv,x-> not IsZero(x));

while (frst <= k) and not (frst = fail)  do
vv:=vv+B[Hds[frst]];
pv:=pv+B[Hds[frst]];

frst:=PositionProperty(vv,x-> not IsZero(x));
od;

return pv;

end;
#############################################################

for n in [1..N!.ambientDimension/Order(G)] do
ProjN:=ProjectedFpGModule(N,n)!.matrix;
B:=[];

Mapp:=List(Map,m->TransposedMat(m));
Mapp:=TransposedMat(Mapp);
Mapp:=List(Mapp,m->Flat(m));
Mapp:=List(Mapp,v->PreimageVector(v,(n-1)*Order(G)));

for j in [1..dd] do
Append(B,List(ProjN,v->ExpandVector(v,j)));
od;

A:=VecToMat(RandomVector(SumIntersectionMat(B,DualBasisM)[2]));
A:=A+Mapp;


Append(Map,[TransposedMat(A)]);
Unbind(A);


od;

Map:=List(Map,m->TransposedMat(m));
Map:=TransposedMat(Map);
Map:=List(Map,m->Flat(m));

return FpGModuleHomomorphismNC(M,M,Map);

end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(CompositionOfFpGModuleHomomorphisms,
function(g,f)			#homomorphisms f:L->M, g:M->N
local 
	fg, Imagef;

if not Source(g)=Target(f) then
Print("Homomorphisms can't be composed.\n");
return fail; fi;

Imagef:=VectorsToFpGModuleWords(f!.target,f!.matrix);

return Imagef; 

#OBVIOUSLY TO BE FINISHED

end);
#####################################################################
#####################################################################
