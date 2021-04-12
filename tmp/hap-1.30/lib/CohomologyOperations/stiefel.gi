
###################################################
###################################################
InstallGlobalFunction(PolytopalRepresentationComplex,
function(arg)
local phi, v, G, imG, gens, P,
      Dimension, Boundary, EltsG, StabilizerSubgroup,StabAction,ln,
      BoundaryRec,StabilizerSubgroupRec, StabActionRec, iso, bnd,n,k  ;

# Inputs (G,v) or (phi,v). In the first case G must be a matrix group or a permutation
# group which we convert to a permutation matrix group. In the second case phi:G-->Q
# must be a homomorphism to a matrix or permutation group Q.

v:=arg[2];

if IsGroup(arg[1]) then
G:=arg[1];
gens:=GeneratorsOfGroup(G);
phi:=GroupHomomorphismByImages(G,G,gens,gens);
else
phi:=arg[1];
fi;

G:=Source(phi);
imG:=Image(phi);
gens:=GeneratorsOfGroup(imG);
gens:=List(gens,x->x*TransposedMat(x));
gens:=List(gens,x->x=One(imG));
if false in gens then 
Print("The representation must be orthogonal.\n");
return fail;
fi;

P:=PolytopalComplex(imG,v);

if IsPermGroup(imG)then
iso:=PermToMatrixGroup(imG);
else 
gens:=GeneratorsOfGroup(imG);
iso:=GroupHomomorphismByImages(imG,imG,gens,gens);
fi; 

Dimension:=P!.dimension;
ln:=Length(P);
EltsG:=Elements(G);

################################
BoundaryRec:=[];
for n in [1..ln] do
BoundaryRec[n]:=[];
for k in [1..Dimension(n)] do
bnd:=1*P!.boundary(n,k);
Apply(bnd,x->[x[1], Position(EltsG, 
PreImagesRepresentative(phi, PreImagesRepresentative(iso,P!.elts[x[2]])))  ]);
Add(BoundaryRec[n],bnd);
od;
od;
Boundary:=function(n,k);
if n>ln or n<1 then return []; fi;
if k>0 then return BoundaryRec[n][k];fi;
if k<0 then return NegateWord(BoundaryRec[n][-k]);fi;
end;
################################

################################
StabilizerSubgroupRec:=[];
for n in [0..ln] do
StabilizerSubgroupRec[n+1]:=[];
for k in [1..Dimension(n)] do
gens:=GeneratorsOfGroup(P!.stabilizer(n,k));
gens:=List(gens,x-> 
PreImages(phi, PreImagesRepresentative(iso,x)));
gens:=List(gens,x->List(x,a->a));
gens:=Concatenation(gens);
if Length(gens)=0 then StabilizerSubgroupRec[n+1][k]:=Group(One(G));
else
StabilizerSubgroupRec[n+1][k]:=Group(gens);
fi;
od;
od;
StabilizerSubgroup:=function(n,k);
return StabilizerSubgroupRec[n+1][AbsInt(k)];
end;
################################

################################
StabAction:=function(n,k,h)
local hh;
hh:=Position(P!.elts,Image(iso,Image(phi,EltsG[h])));
return P!.action(n,k,hh);
end;
################################


return Objectify(HapNonFreeResolution,
           rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=fail,
            elts:=EltsG,
            group:=G,
            #standardWord:=StandardWord,
            stabilizer:=StabilizerSubgroup,
            #basis:=StabilizerBasis,
            action:=StabAction,
            #hasse:=Hasse,
            properties:=
             [["type","resolution"],
              ["length",ln],
              ["characteristic", 0] ]));

end);
###################################################
###################################################


###################################################
###################################################
InstallGlobalFunction(HAP_StiefelWhitney,
function(arg)
local G,v,n,P,R,T,homid,ThomR,CRhomCT,CThomCR,CR, CT,trans, iso,pos, CTmatCR,
      TRmapping,stief,one,zero,d,u,uu,A,B,B1,B2,i,j;

# Inputs (G,v) or (G,v,n) or (phi,v) or (phi,v,n). In the first two cases G must be a matrix group 
# or a permutation group which we convert to a permutation matrix group. In the second case phi:G-->Q
# must be a homomorphism to a matrix or permutation group Q.

if IsGroup(arg[1]) then G:=arg[1];
else G:=Source(arg[1]);
fi;

v:=arg[2];

P:=PolytopalRepresentationComplex(arg[1],v);

if IsBound(arg[3]) then n:=arg[3]; else 
n:=0; while P!.dimension(n)>0 do n:=n+1; od; 
fi;

T:=FreeGResolution(P,n+1,2);
R:=ResolutionPrimePowerGroup(G,n+1);
homid:=GroupHomomorphismByImages(G,G,  
                            GeneratorsOfGroup(G),GeneratorsOfGroup(G));
ThomR:=EquivariantChainMap(T,R,homid);
CRhomCT:=HomToIntegersModP(ThomR,2);
CR:=Source(CRhomCT);
CT:=Target(CRhomCT);

####
##Now compute cochain map CThomCR:CT--->CR
CTmatCR:=[];
for i in [0..n+1] do
u:=List([1..R!.dimension(i)],x->0);
A:=[];
for j in [1..R!.dimension(i)] do
u:=0*u;
u[j]:=One(GF(2));
Add(A,CRhomCT!.mapping(u,i));
od;
B1:=SemiEchelonMatTransformation(A);

uu:=List([1..Length(B1.heads)],x->0);
B2:=[];
for j in [1..Length(B1.heads)] do
uu:=0*uu;
if B1.heads[j]=0 then
uu[j]:=One(GF(2));
Add(B2,1*uu);
fi;
od;

B:=SolutionsMatDestructive(One(GF(2))*Concatenation(A,B2),
                                      One(GF(2))*IdentityMat(Length(A[1])));
B:=List(B,x->x{[1..Length(A)]});


Add(CTmatCR,B);
od;
# The cochain map CT--->CR sends vector v to v*CTmatCR 

#######
TRmapping:=function(v,n);
return v*CTmatCR[n+1];
end;
#######
CThomCR:=rec(source:=CT, target:=CR,properties:=CRhomCT!.properties, mapping:=TRmapping);
CThomCR:=Objectify(HapCochainMap,CThomCR);
##DONE
####


####
####
#We now construct a function stief(v,k):CR^k--CR^(n+d), v|-->w
#where d is the dimension of the polytope P
d:=0;
while P!.dimension(d+1)>0 do d:=d+1; od;
zero:=Zero(GF(2));
one:=One(GF(2));

####
stief:=function(v,k)
local w;
w:=zero*[1..T!.filteredDimension(d-1,k+d)];
Append(w,one*v);
return w*CTmatCR[1+k+d];
end;
####


return [stief,d,R,P] ;


end);
###################################################
###################################################


###################################################
###################################################
InstallGlobalFunction(FundamentalMultiplesOfStiefelWhitneyClasses,
function(arg)
local bool,G,v,A,P,N,S,R,stief,d,L,i,k,j,fund,Bas,Bas1,Bas2,swc,swc1,PIRep,u,a,b,w;

# Inputs either (G,v,A) where G is either a matrix/permutation group or a group representation.

G:=arg[1];
v:=arg[2];
A:=arg[3];
if Length(arg)=4 then bool:=arg[4]; else bool:=false; fi;

S:=HAP_StiefelWhitney(G,v);
P:=S[4];
N:=0; while P!.dimension(N+1)>0 do N:=N+1; od; 
stief:=S[1];
d:=S[2];
R:=S[3];

Bas:=Basis(A);

L:=[];
for k in [0..N] do
for j in [1..R!.dimension(k)] do
w:=Zero(A);
if k+d<=N then
u:=Zero(GF(2))*[1..R!.dimension(k)];
u[j]:=One(GF(2));
u:=stief(u,k);

a:=Sum(List([0..d+k-1],s->R!.dimension(s)));
b:=a+R!.dimension(d+k);

for i in [a+1..b] do
w:=w+u[i-a]*Bas[i];
od;

fi;

Add(L,w); 
od;
od;

fund:=L[1];
#if IsZero(fund) then
#Print("Method is not successful on this example.\n"); return fail;
#fi;

swc:=[];
for i in [0..N] do
Add(swc,Sq(A,i,fund));
od;

if not bool then
return swc;
fi;


swc1:=[];
for i in [0..N] do

Bas1:=Filtered(Basis(A),x->A!.degree(x)=i);
Bas2:=List(Bas1,i->i*fund);
L:=LeftModuleGeneralMappingByImages(Subspace(A,Bas1),Subspace(A,Bas2),Bas1,Bas2);

#####
PIRep:=function(L,x);
if x=fail then return x; fi;
if x=Zero(A) then return Elements(Kernel(L)); fi;
return Elements(PreImages(L,x));
end;
#####

Add(swc1,PIRep(L,swc[i+1]));

od;
return swc1;

end);
###################################################
###################################################

###################################################
###################################################
InstallGlobalFunction(CohomologyHomomorphismOfRepresentation,
function(G,v,A)
local N,S,R,stiefel,d,L,i,k,j,Bas,w,u,a,b;

# The input G can be a matrix/permutation group or a group representation.
# A is the Mod2 Steenrod algebra up to some degree.
# v is a vector on which the group (representation) acts.

N:=Maximum(List(Basis(A),x->A!.degree(x)));
S:=HAP_StiefelWhitney(G,v,N);
stiefel:=S[1];
d:=S[2];
R:=S[3];

Bas:=Basis(A);

L:=[];
for k in [0..N] do
for j in [1..R!.dimension(k)] do
w:=Zero(A);
if k+d<=N then
u:=Zero(GF(2))*[1..R!.dimension(k)];
u[j]:=One(GF(2));
u:=stiefel(u,k);

a:=Sum(List([0..d+k-1],s->R!.dimension(s)));
b:=a+R!.dimension(d+k);

for i in [a+1..b] do
w:=w+u[i-a]*Bas[i];
od;

fi;

Add(L,w);
od;
od;

L:=LeftModuleGeneralMappingByImages(A,A,Bas,L);
A!.StiefelWhitneyHomomorphism:=L;


return L;

end);
###################################################
###################################################

