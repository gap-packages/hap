###############################################################
###############################################################
InstallGlobalFunction(HAP_coho_isoms,
function(R,S,A,Hn,n)
local G,a,b,x,bhomc,B,C,psi,delta,GhomG,SmapR,BasA,Basn,
      gensHnA, HnAgrp, gensHnAgrp, hn, natn, Kn, Cn,fun,fun2,
      genshn,imgenshn,hn2HnA,xx,HnA2hn,imgensHnA,RmapSdual,SmapRdual,matRS,
      matSR, w,i,j,jj, HnAgrp2hn, hn2HnAgrp,p;

#R is a minimal resolution of GF(2) for the 2-group G.
#S is any resolution of GF(2) (or even Z) for the 2-group G.
#A is the algebra constructed from ModPCohomologyRing(R).
#Hn is H^n(Hom(S,K)) where K is any GOuter group representation of
#the trivial G-module GF(2)
#We'll return a pair of isomorphisms [hn2HnA, HnA2hn] where hn is the
#cohomology *group* H^n(S,K) and HnA is a sub magma of A. 

p:=EvaluateProperty(R,"characteristic");

#########################
if not IsPrimeInt(p) then 
Print("The first input resolution must be a minimal resolution of prime characteristic.\n");
return fail;
fi;
#########################

G:=R!.group;
if not G=S!.group then return fail; fi;
GhomG:=GroupHomomorphismByFunction(G,G,x->x);
SmapR:=EquivariantChainMap(S,R,GhomG);

matSR:=NullMat(S!.dimension(n),R!.dimension(n));
for i in [1..S!.dimension(n)] do
w:=SmapR!.mapping([[i,1]],n);
w:=List(w,x->x[1]);
for j in w do
jj:=AbsInt(j);
matSR[i][jj]:=matSR[i][jj]+SignInt(j) mod p;
od;
od;

matSR:=TransposedMat(matSR) mod p;

###################################
SmapRdual:=function(w,n)
local v,x,j,a;
v:=[];

for x in w do
for j in [1..S!.dimension(n)] do
#if not matSR[AbsInt(x[1])][j]=0 then Add(v,[j,1]); fi;
a:=matSR[AbsInt(x[1])][j];
if not a=0 then 
v:=AddFreeWords(v,MultiplyWord(a,[[j,1]]),p); fi;
od;
od;

return AlgebraicReduction(v,p);
end;
###################################


BasA:=CanonicalBasis(A);
Basn:=Filtered([1..Length(BasA)],i->A!.degree(BasA[i])=n);
gensHnA:=Filtered(BasA,i->A!.degree(i)=n);

hn:=Hn!.ActedGroup;
natn:=Hn!.nat; natn:=natn!.Mapping;
Kn:=Source(natn);
Cn:=Kn!.ParentAttr;

genshn:=Pcgs(hn);

xx:=GeneratorsOfGroup(Source(Embedding(Cn,1)))[1];
###################
fun:=function(w)
local v,y;
v:=One(Cn);
for y in w do
v:=v*Image(Embedding(Cn,AbsInt(y[1])),xx)^SignInt(y[1]);
od;
return v;
end;
###################

imgensHnA:=List([1..R!.dimension(n)],i->SmapRdual([[i,1]],n));
imgensHnA:=List(imgensHnA,x->fun(x));
imgensHnA:=List(imgensHnA,x->Image(natn,x));

##################
HnA2hn:=function(w)
local c, v, i;
c:=Coefficients(BasA,w);
c:=c{Basn};
v:=One(hn);
for i in [1..Length(c)] do
#if not IsZero(c[i]) then v:=v*imgensHnA[i]; fi;
v:=v*imgensHnA[i]^IntFFE(c[i]);
od;
return v;
end;
##################

HnAgrp:=AbelianGroup(List(gensHnA,i->p));
gensHnAgrp:=Pcgs(HnAgrp);
HnAgrp2hn:=GroupHomomorphismByImages(HnAgrp,hn,gensHnAgrp,
List([1..Length(gensHnA)],i-> HnA2hn(gensHnA[i])) );

hn2HnAgrp:=GroupHomomorphismByImages(hn,HnAgrp, genshn,
List(genshn,x->PreImagesRepresentative(HnAgrp2hn,x)));

###################
hn2HnA:=function(w)   
local v, z, i;
v:=Image(hn2HnAgrp,w);
v:=ExponentsOfPcElement(gensHnAgrp,v);
z:=Zero(A);
for i in [1..Length(v)] do
z:=z+v[i]*gensHnA[i];
od;
return z;
end;
###################


return [hn2HnA, HnA2hn];
end);
###############################################################
###############################################################

###############################################################
###############################################################
InstallGlobalFunction(ModPSteenrodAlgebra,
function(arg)
local G,N,R,A,S,p,
      x,a,b,B,C,psi,bhomc,delta,i, Sq0,Bok,AhomH,HhomA,K,maxdeg;

G:=arg[1];
p:=PrimePGroup(G);
N:=arg[2];
if Length(arg)>2 then R:=arg[3]; else
R:=ResolutionPrimePowerGroup(G,N+1);fi;
R!.properties[PositionProperty(R!.properties,x->x[1]="length")][2]:=N;
A:=ModPCohomologyRing(R);
if Length(arg)>3 then S:=arg[4]; 
else
S:=ResolutionGenericGroup(G,N+1); fi;

x:=(1,2,3,4);;
x:=[1..p^2];
x:=x{[2..p^2]};
x[p^2]:=1;
x:=PermList(x);
a:=Group(x^p);;
b:=Group(x);;
bhomc:=NaturalHomomorphismByNormalSubgroup(b,a);
B:=TrivialGModuleAsGOuterGroup(G,b);
C:=TrivialGModuleAsGOuterGroup(G,Image(bhomc));
psi:=GOuterGroupHomomorphism();
psi!.Source:=B;
psi!.Target:=C;
psi!.Mapping:=bhomc;

delta:=[];
for i in [1..N-1] do
delta[i]:=ConnectingCohomologyHomomorphism(psi,i,S);;
od;

AhomH:=[];
for i in [1..N-1] do
K:=Source(delta[i]);
AhomH[i]:=HAP_coho_isoms(R,S,A,K,i)[2];
od;

HhomA:=[];
for i in [1..N-1] do
K:=Target(delta[i]);
HhomA[i]:=HAP_coho_isoms(R,S,A,K,i+1)[1];
od;

########################
Bok:=function(w)
local n, v,del,iso;
n:=A!.degree(w);
if n=0 then return Zero(A); fi;
iso:=AhomH[n];
v:=iso(w);
del:=delta[n];
del:=del!.Mapping;
v:=Image(del,v);
iso:=HhomA[n];
v:=iso(v);
return v;
end;
########################

########################
Sq0:=function(w);
return w;
end;
########################
A!.squares:=[Sq0];
if p=2 then A!.squares[2]:=Bok; fi;
A!.bockstein:=Bok;
A!.maxdeg:=Maximum(List(CanonicalBasis(A),x->A!.degree(x)));
return A;
end);
###############################################################
###############################################################

###############################################################
###############################################################
InstallMethod(Sq,
"steenrod squares for Mod 2 cohomology rings",
[IsAlgebra,IsInt,IsObject],
function(A,n,w)
local W, M, v, x, i, MAX, sqq, a,b;

####################################################
## This function makes use of the Cartan relations. At present it
## does not make any use of the Adems relations.
####################################################

#### Are Steenrod squares defined at all?###########
if not IsBound(A!.squares) then
#Print("Steenrod squares are not defined for this algebra.\n");
return fail;
fi;
####################################################

#### n=0 ###########################################
if n=0 then return w; fi;
####################################################

M:=HAP_MultiplicativeGenerators(A);
W:=M[2](w);

#### Sq^n=0 if n> degrees of all homogeneous parts##
if Length(W)=0 then MAX:=0; else
MAX:=Maximum(List(Flat(W),A!.degree)); fi;
if n>MAX then return Zero(A); fi;
####################################################

#### Sq^n(w) not defined if maxdeg<MAX+n #############
if A!.maxdeg<MAX+n then
#Print("Steenrod square image has too high a degree.\n");
return fail;
fi;
###################################################### 

#### additivity: apply to the homogeneous parts ####
if Length(W)>1 then
    v:=Zero(A);
    for x in W do
        v:=v+Sq(A,n,Product(x));
    od;
    return v;
fi;
####################################################

#### So now W is homogeneous #######################
#### We remove outer brackets of W  ################
W:=W[1];

##### if Degree(W)=n then Sq^n(W)=W^2 ##############
if n=Sum(List(W,A!.degree)) then W:=Product(W);
return W^2; fi;
####################################################

### Length(W)>1 : so W is a product of generators ##
if Length(W)>1 then
v:=Zero(A);
for i in [0..n] do
a:=Sq(A,i,W[1]);
if a=fail then return fail; fi;
b:=Sq(A,n-i,Product(W{[2..Length(W)]}));
if b=fail then return fail; fi;
v:=v+ a*b;
#v:=Sq(A,i,W[1])*Sq(A,n-i,Product(W{[2..Length(W)]}));
od;
return v;
fi;
####################################################

#### Now W is a list of just one ring generator ####
if not IsBound(A!.squares[n+1]) then return fail; fi;
sqq:=A!.squares[n+1];

return sqq(w);
end);
##################################################################
##################################################################

###############################################################
###############################################################
InstallMethod(Bockstein,
"Bockstein for Mod p cohomology rings",
[IsAlgebra,IsObject],
function(A,w)
local W, M, v, x, MAX;

#### Is the Bockstein defined at all?##############
if not IsBound(A!.bockstein) then
Print("The Bockstein operation is not defined for this algebra.\n");
return fail;
fi;
####################################################

M:=HAP_MultiplicativeGenerators(A);
W:=M[2](w);

####################################################
if Length(W)=0 then MAX:=0; else
MAX:=Maximum(List(Flat(W),A!.degree)); fi;
if A!.maxdeg<MAX+1 then
#Print("Bockstein image has too high a degree.\n");
return fail;
fi;
####################################################


#### additivity: apply to the homogeneous parts ####
if Length(W)>1 then
    v:=Zero(A);
    for x in W do
        v:=v+Bockstein(A,Product(x));
    od;
    return v;
fi;
####################################################

#### So now W is homogeneous #######################
#### We remove outer brackets of W and multiply  ###
W:=W[1];
W:=Product(W);
return A!.bockstein(W); 
####################################################

end);
####################################################;
####################################################
