
###############################################
###############################################
InstallGlobalFunction(HAP_PHI,
function(G,NN)
local N, phiRec, phiRecT, Diag,Diagrec, C2, C2G, GG, RG, RGG, RC2, RC2G, 
phi_0,phi,Tau,tauhom,tau,taurec,MultGG, FirstEmbedding, SecondEmbedding, Fproj, Sproj,FirstProjection,
SecondProjection, n, MT;
N:=NN+1;
## G is a group
## N is a positive integer
## This function will return a triple [phi,RC2G,RGG] where RC2G is a
## free resolution for C_2xG, RGG is a free resolution for GxG, and  
## phi(n,a,i) is a function that inputs an integer n in [0..N],
## inputs an element  a=(a1,a2) in C_2xG, and inputs an integer i
## corresponding to the i-th generator of degree n in the resolution
## RC2G. If i is negative it represents minus the i-th generator.
## The function phi(n,a,i) outputs a list of lists 
##        [[m_1,g_1], [m_2,g_2], ..., [m_k,g_k]] 
## representing, in the usual fashion, an element of degree n in RGG.
## The function phi represents a chain map phi:RC2G --> RGG.

C2:=Group((1,2));
RG:=ResolutionPrimePowerGroup(G,N);
RC2:=ResolutionPrimePowerGroup(C2,N);
RGG:=ResolutionFiniteDirectProduct(RG,RG);
RC2G:=ResolutionFiniteDirectProduct(RC2,RG);
#RGG:=ResolutionDirectProductLazy(RG,RG);
#RC2G:=ResolutionDirectProductLazy(RC2,RG);

C2G:=RC2G!.group;
FirstProjection:=Projection(C2G,1);
SecondProjection:=Projection(C2G,2);
GG:=RGG!.group;
FirstEmbedding:=Embedding(GG,1);
SecondEmbedding:=Embedding(GG,2);
Fproj:=Projection(GG,1);
Sproj:=Projection(GG,2);
tauhom:=GroupHomomorphismByFunction(GG,GG,x->
ImageElm(SecondEmbedding,ImageElm(Fproj,x))*ImageElm(FirstEmbedding,ImageElm(Sproj,x))  );

################################
tau:=function(k);
return Position(RGG!.elts,ImageElm(tauhom,RGG!.elts[k]));
end;
################################
taurec:=List([1..Order(GG)],i->tau(i));

phiRec:=List([1..N+1],i->[]);
phiRecT:=List([1..N+1],i->[]);

################################
Diag:=function(x);
# G-->GxG, x--> (x,x)
return Position(RGG!.elts,ImageElm(FirstEmbedding,x)*ImageElm(SecondEmbedding,x));
end;
################################
Diagrec:=List([1..Order(G)],i->Diag(Elements(G)[i]));
Diag:=function(x); return Diagrec[Position(Elements(G),x)];end;

################################
Tau:=function(n,kk)
local  x,y,z,k;
k:=AbsInt(kk);
x:=RGG!.Int2Vector(n,k); #x=[p,q,r,s]
y:=RGG!.Vector2Int(x[2],x[1],x[4],x[3]);
return y;  #working mod 2 so this is ok
end;
################################

MT:=MultiplicationTable(GG);
################################
MultGG:=function(i,j); return MT[i][j]; end;
################################

################################
phi_0:=function(a,ii)
local w, a2,i;
## We assume that R!.dimension(0)=1. This is correct since
## R:=ResolutionFiniteGroup(G,N). For other resolutions we'll need
## to generalize this function.
i:=AbsInt(ii);
## So i=1 since R!.dimension(0)=1.
## I think the appropriate embedding G --> GG is a2 --> (a2,a2).
a2:=ImageElm(SecondProjection,a);
w:=[  [i, Diag(a2)  ]  ];  
return w; #ok since we work mod 2
end;
################################

################################
phi:=function(n,a,ii)
local w, a1,a2,i;

###############
if n=0 then return phi_0(a,ii); fi;
###############

i:=AbsInt(ii);
a1:=ImageElm(FirstProjection,a);
a2:=ImageElm(SecondProjection,a);

###############Don't calculate phi(n,1,i) twice! 
if
(Order(a1)=1 and (not IsBound(phiRec[n+1][i])))
or
(Order(a1)=2 and (not IsBound(phiRecT[n+1][i])))
then

w:=1*RC2G!.boundary(n,i);       
w:=List(w,x->  phi(n-1,RC2G!.elts[x[2]],x[1]));
w:=1*Concatenation(w);
w:=AlgebraicReduction(w,2); #######
w:=List(w,x->  RGG!.homotopy(n-1,x));
w:=1*Concatenation(w);
w:=AlgebraicReduction(w,2); #######
if Order(a1)=2 then 
w:=List(w,x-> [Tau(n,x[1]),taurec[x[2]]]);
fi;

if Order(a1)=1 then phiRec[n+1][i]:=w;
else phiRecT[n+1][i]:=w; fi;

fi;
############

if Order(a1)=1 then w:=1*phiRec[n+1][i];
else w:=1*phiRecT[n+1][i]; fi;

Apply(w,x->[x[1],MultGG(Diag(a2),x[2])]); 

return w;  #ok since we work mod 2   
end;
################################


return [phi, RC2G, RGG, RG];
end);
###############################################
###############################################



###############################################
###############################################
InstallGlobalFunction(Mod2SteenrodAlgebra,
function(arg)
local G,N,RG, RGG, RC2G, BasA, Bas, bas, cupgen, cup, sqr, phi, L, A, i,n,one;
# This function inputs a finite 2-group G and an integer N.
# It returns the steenrod algebra up to degree N, including the
# cup-i product version of the Steenrod square.

if Length(arg)=2 then
G:=arg[1];
N:=arg[2];
fi;

L:=HAP_PHI(G,N);
phi:=L[1];
RC2G:=L[2];
RGG:=L[3];
one:=Identity(RC2G!.group);
RG:=L[4]; 
A:=ModPCohomologyRing(RG);
BasA:=CanonicalBasis(A);
A!.maxdeg:=Maximum(List(BasA,x->A!.degree(x)));

bas:=[];
Bas:=[];
for n in [0..N] do
bas[n+1]:=Filtered([1..Length(BasA)],i->A!.degree(BasA[i])=n);
Bas[n+1]:=Filtered(BasA,b->A!.degree(b)=n);
od;

###############################################
cupgen:=function(n,i,j,k,m)
local v,w,a,x;
# This returns the m-th coordinate of the cup-i product e_j u e_k
# where e_j and e_k are generators of H^n(G).

v:=RC2G!.Vector2Int(i,2*n-i,1,m);
v:=phi(2*n,one,v); 
w:=RGG!.Vector2Int(n,n,j,k); 
a:=0;
for x in v do 
if AbsInt(x[1])=w then a:=a+1 mod 2; fi;
od;

return a*Bas[2*n-i+1][m];
end;
###############################################


###############################################
cup:=function(i,w)
local n, v, c,j,k,m;
v:=Zero(A);
n:=A!.degree(w);
if 2*n<= i then return v; fi;
c:=Coefficients(BasA,w);
c:=c{bas[n+1]};
c:=Filtered([1..Length(c)],i->not IsZero(c[i])); 

for j in c do
for k in c do
for m in [1..RG!.dimension(2*n-i)] do
v:=v+cupgen(n,i,j,k,m);
od;od;od;

return v;
end;
###############################################

###############################################
sqr:=function(i,w)
local n;
if i=0 then return w; fi;
n:=A!.degree(w); 
if n<i then return Zero(A); fi;
return cup(n-i,w);
end;
###############################################

A!.squares:=List([0..N],i->(w->sqr(i,w)));
A!.res:=RG;
return A;
end);
###############################################
###############################################





