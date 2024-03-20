###############################################################
###############################################################
InstallGlobalFunction(HAP_coho_isoms,
function(R,S,A,Hn,n)
local G,a,b,x,bhomc,B,C,psi,delta,GhomG,SmapR,BasA,Basn,
      gensHnA, HnAgrp, gensHnAgrp, hn, natn, Kn, Cn,fun,fun2,
      genshn,imgenshn,hn2HnA,xx,HnA2hn,imgensHnA,RmapSdual,SmapRdual,matRS,
      matSR, w,i,j,jj, HnAgrp2hn, hn2HnAgrp,p;

#R is a minimal resolution of GF(p) for the p-group G.
#S is any resolution over Z for the p-group G.
#A is the algebra constructed from ModPCohomologyRing(R).
#Hn is H^n(Hom(S,K)) where K is any GOuter group representation of
#the trivial G-module GF(p)
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
a:=SignInt(x[1])*matSR[AbsInt(x[1])][j];
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

#for x in GeneratorsOfGroup(hn) do
#if not x = HnA2hn(hn2HnA(x)) then Print("Ooops!\n"); fi;
#od;

return [hn2HnA, HnA2hn];
end);
###############################################################
###############################################################

###############################################################
###############################################################
InstallGlobalFunction(ModPSteenrodAlgebra,
function(arg)
local G,N,R,A,S,p,mx,
      x,a,b,B,C,psi,bhomc,delta,i, Sq0,Bok,AhomH,HhomA,K,maxdeg;

G:=arg[1];
p:=PrimePGroup(G);
N:=arg[2];
if Length(arg)>2 then R:=arg[3]; else
R:=ResolutionPrimePowerGroup(G,N+1);fi;
R!.properties[PositionProperty(R!.properties,x->x[1]="length")][2]:=N;
A:=ModPCohomologyRing(R);
mx:=ModPRingGenerators(A);;
mx:=List(mx,A!.degree);;
mx:=Maximum(mx)+1;
mx:=Minimum(mx,N);
if Length(arg)>3 then S:=arg[4]; 
else
S:=ResolutionGenericGroup(G,mx+2); fi;

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
for i in [1..Minimum(mx,N-1)] do
delta[i]:=ConnectingCohomologyHomomorphism(psi,i,S);;
od;

AhomH:=[];
for i in [1..Minimum(mx,N-1)] do
K:=Source(delta[i]);
AhomH[i]:=HAP_coho_isoms(R,S,A,K,i)[2];
od;

HhomA:=[];
for i in [1..Minimum(mx,N-1)] do
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
A!.AhomH:=AhomH;
A!.HhomA:=HhomA;
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
local W, WW, M, v, x, i, MAX, sqq, a,b,V;

####################################################
## This function makes use of the Cartan relations. At present it
## does not make any use of the Adems relations.
####################################################

if not Characteristic(A)=2 then return fail; fi;

#### Is the Bockstein defined?##############
if not IsBound(A!.squares) and IsBound(A!.chainComplex) then
A!.maxdeg:=Length(A!.chainComplex);
A!.complete:=true;
A!.squares:=[];
A!.squares[1]:=function(x); return x; end;
A!.bockstein:=HAP_bockstein(A);
A!.squares[2]:= A!.bockstein;
fi;
####################################################

#### Are Steenrod squares defined at all?###########
if not IsBound(A!.squares) then
return fail;
fi;
####################################################

#### n=0 ###########################################
if n=0 then return w; fi;
####################################################

M:=HAP_MultiplicativeGenerators(A);
W:=M[3](w);

#### Sq^n=0 if n> degrees of all homogeneous parts##
if Length(W)=0 then MAX:=0; else
WW:=List(W,  x->List(x,b->A!.degree(b)) );
WW:=List(WW, x->Sum(x) );
MAX:=Maximum(WW); fi;
if n>MAX then return Zero(A); fi;
####################################################

#### Sq^n(w) not defined if maxdeg<MAX+n #############
if A!.maxdeg<MAX+n then
#Print("Steenrod square image has too high a degree.\n");
if IsBound(A!.complete) then
  if A!.complete=true then return Zero(A); fi;
fi;
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
V:=W[1];

##### if Degree(W)=n then Sq^n(W)=W^2 ##############
if n=Sum(List(V,A!.degree)) then V:=Product(V);
return V^2; fi;
####################################################

### Length(V)>1 : so V is a product of generators ##
if Length(V)>1 then
v:=Zero(A);
for i in [0..n] do
a:=Sq(A,i,V[1]);
if a=fail then return fail; fi;
b:=Sq(A,n-i,Product(V{[2..Length(V)]}));
if b=fail then return fail; fi;
v:=v+ a*b;
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
local W, WW, V, M, i, v, x, MAX, a, b,c, gens, gensbas;;

#### Is the Bockstein defined?##############
if not IsBound(A!.bockstein) and IsBound(A!.chainComplex) then
A!.bockstein:=HAP_bockstein(A);
fi;
####################################################


#### Is the Bockstein defined at all?##############
if not IsBound(A!.bockstein) then
Print("The Bockstein operation is not defined for this algebra.\n");
return fail;
fi;
####################################################
#
#### Bockstein for CW complexes ####################
if not IsBound(A!.maxdeg) then
return A!.bockstein(w);
fi;
####################################################

###### If w=0 then return 0 ########################
if IsZero(w) then return w; fi;
####################################################

M:=HAP_MultiplicativeGenerators(A);
W:=M[3](w);
gens:=ModPRingGenerators(A);
gensbas:=Basis(Submodule(A,gens),gens);
if not IsBound(A!.bocksteinrec) then
A!.bocksteinrec:=[]; fi;

####################################################
if Length(W)=0 then MAX:=0; else
WW:=List(W,  x->List(x,b->A!.degree(b)) );
WW:=List(WW, x->Sum(x) );
MAX:=Maximum(WW); fi;
if A!.maxdeg<MAX+1 then
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

###### Now Length(W)=1 #############################
V:=W[1];
if Length(V)=1 then 
c:=Coefficients(gensbas,V[1]);
c:=Sum(c);
i:=Position(gens,c^-1*V[1]);
if not IsBound(A!.bocksteinrec[i]) then
A!.bocksteinrec[i]:=A!.bockstein(c^-1*V[1]); fi;
return c*A!.bocksteinrec[i]; fi;

a:=V[1];
b:=Product(V{[2..Length(V)]});
return A!.bockstein(a)*b + (-1)^(A!.degree(a))*a*Bockstein(A,b);

####################################################

end);
####################################################;
####################################################

########################################################
########################################################
InstallOtherMethod(Bockstein,
"Bockstein for chain complexes",
[IsHapChainComplex,IsInt,IsInt],
function(C,n,prime)
return HAP_chain_bockstein(C,n,prime);
end);
########################################################
########################################################


########################################################
########################################################
InstallGlobalFunction(BocksteinHomology,
function(A,n)
local Bas, gensn, gensnm1, B, Z ;

Bas:=Basis(A);;
gensn:=Filtered(Bas,x->A!.degree(x)=n);
gensnm1:=Filtered(Bas,x->A!.degree(x)=n-1);
B:=Submodule(A,List(gensnm1, x->Bockstein(A,x)));;
B:=Dimension(B);
Z:=Submodule(A,List(gensn, x->Bockstein(A,x)));;
Z:=Dimension(Submodule(A,gensn)) - Dimension(Z);

return Z-B;
end);
########################################################
########################################################

########################################################
########################################################
InstallGlobalFunction(PrintAlgebraWordAsPolynomial,
function(arg)
local  A,w,M, B, e, c,d,x,p,i,j,str;

A:=arg[1];
w:=arg[2];

M:=HAP_MultiplicativeGenerators(A);
e:=M[2](w);

B:=List(M[1],x->Product(x));
B:=Basis(A,B);

str:=[];
c:=Coefficients(B,w);
d:=Filtered([1..Length(c)], i-> not IsZero(c[i]));
for i in d do
     if i<Length(c) then
        if not IsOne(c[i]) then #Print(c[i],"*"); 
                                 Append(str,"*"); fi;
     else
        if not IsOne(c[i]) then #Print(c[i]); 
                                  Append(str, String(c[i])); fi;
     fi;
     for j in [1..Length(M[1][i])] do
        if j < Length(M[1][i]) then #Print(M[1][i][j],"*"); 
        Append(str,String(M[1][i][j])); Append(str,"*");
        else #Print(M[1][i][j]); 
              Append(str,String(M[1][i][j]));
        fi;
     od;
     if not i=d[Length(d)] then #Print(" + "); 
                                 Add(str,'+'); fi;
od;
if Length(arg)=3 then 
return str;
fi;
end);
########################################################
########################################################

