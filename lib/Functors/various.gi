#(C) Graham ellis, 2005-2006

#####################################################################
InstallGlobalFunction(HAP_PrintTo,
function(arg)
local file, out,i;

file:=arg[1];
out:=OutputTextFile(file,false);
for i in [2..Length(arg)] do
WriteAll(out,String(arg[i]));
od;
CloseStream(out);

end);
#####################################################################

#####################################################################
InstallGlobalFunction(HAP_AppendTo,
function(arg)
local file, out,i;

file:=arg[1];
out:=OutputTextFile(file,true);
for i in [2..Length(arg)] do
WriteAll(out,String(arg[i]));
od;
CloseStream(out);
end);
#####################################################################


#####################################################################
InstallGlobalFunction(EvaluateProperty,
function(X,name)
local
	x;

if not "properties" in NamesOfComponents(X) then return fail; fi;

x:=First(X!.properties,v->v[1]=name);
if x=fail then return fail;
else return x[2]; fi;

end);
#####################################################################

#####################################################################
InstallGlobalFunction(EvenSubgroup,
function(G)
local
	x,y,gens;

gens:=[];
for x in GeneratorsOfGroup(G) do
for y in GeneratorsOfGroup(G) do
Append(gens,[x*y]);
od;
od;

return Group(gens);
end);
#####################################################################

#####################################################################
InstallGlobalFunction(ReduceGenerators,
function(gens,G)    		#I should probably use the
                                #Frattini subgroup!
local x,newgens;

for x in gens do
newgens:=DifferenceLists(gens, [x]);
#if Order(Group(Concatenation(newgens,[Identity(G)])))=Order(G) then
if Group(Concatenation(newgens,[Identity(G)]))=G then
return ReduceGenerators(newgens,G); fi;
od;
return SSortedList(gens);
end);
#####################################################################

#####################################################################
InstallGlobalFunction(ReduceGenerators_alt,
function(gens,G)             
                              
local x,newgens;

newgens:=[gens[1]];

for x in gens do
if not x in Group(newgens) then Add(newgens,x); fi;
if Order(Group(newgens))=Order(G) then break; fi;
od;

return newgens;
end);
#####################################################################


#####################################################################
InstallGlobalFunction(AbelianInvariantsToTorsionCoefficients,
function(L)
local
	primes, invariants,max,coeffs,i,x;

primes:=SSortedList(Factors(Product(L)));
invariants:=[];
for i in [1..Length(primes)] do
invariants[i]:=Filtered(L,x->(primes[i] in Factors(x)));
od;

max:=Maximum(List(invariants,x->Length(x)));

for x in invariants do
while Length(x)<max do
Append(x,[1]);
od;
od;

coeffs:=[];
for i in [1..max] do
coeffs[i]:= Product(List([1..Length(primes)],j->invariants[j][i]));
od;

return coeffs;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(TorsionGeneratorsAbelianGroup,
function(G)
local
        L,L1,gens,x,y;

L:=IndependentGeneratorsOfAbelianGroup(G);
L1:=SSortedList(StructuralCopy(L));
gens:=[];

while Length(L1)>0 do
x:=L1[1];
RemoveSet(L1,x);
	for y in L1 do
	if Gcd(Order(x),Order(y))=1 then x:=x*y; 
	RemoveSet(L1,y);
	fi;
	od;
Append(gens,[x]);
od;

return gens;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(BigStepLCS,
function(G,n)
local LCS,BSLCS,i;

LCS:=LowerCentralSeries(G);;
BSLCS:=[LCS[1]];

for i in [2..Length(LCS)] do
if Order(LCS[i])=1 or
Order(BSLCS[Length(BSLCS)])/Order(LCS[i])>n then
Append(BSLCS,[LCS[i]]);
fi;
od;

return BSLCS;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(CoClass,
function(G)
local
	facs;

facs:=SSortedList(Factors(Order(G)));
if Length(facs)>1 then 
Print("G should be a prime-power group. \n");
return fail; 
fi;

return LogInt(Order(G),facs[1])-NilpotencyClassOfGroup(G);
end);
#####################################################################

#####################################################################
InstallGlobalFunction(BoundaryMatrix,
function(C,n)
local
        M,i,j;

M:=[];
for i in [1..C!.dimension(n)] do
M[i]:=C!.boundary(n,i);
od;

return TransposedMat(M);

end);
####################################################################

#####################################################################
InstallGlobalFunction(PrankAlt,
function(G)
local H,S,x, p;

p:=SSortedList(Factors(Order(G)));

if not Length(p)=1 then
Print("G must be a finite p-group. \n");
return fail;
fi;

H:=[];
for x in G do
if Order(x)=p then Add(H,x); fi;
od;
H:=Group(H);

p:=p[1];
S:=LatticeSubgroups(H);
S:=ConjugacyClassesSubgroups(S);
S:=List(S,x->ClassElementLattice(x,1));
S:=Filtered(S,x->IsElementaryAbelian(x));
S:=List(S,x->Order(x));

return Log(Maximum(S),p);

end);
####################################################################

#####################################################################
InstallGlobalFunction(Prank,
function(G)
local
	pPowers,
	xP,
	x,y,
	p,
	AbSubGrps,
	AbSubGrps1,
	AbSubGrps2,
	S, Sx,
	A,
	N,
	N1,
	reps,
	toggle;

###########################################################
p:=SSortedList(Factors(Order(G)));
if not Length(p)=1 then
Print("G must be a finite p-group. \n");
return fail;
fi;
###########################################################

p:=PrimePGroup(G);

###########################################################
if IsAbelian(G) then return
Length(AbelianInvariants(G)); fi;
###########################################################

pPowers:=[];
xP:=[];



AbSubGrps:=[];
AbSubGrps1:=[GeneratorsOfGroup(Center(G))];
AbSubGrps2:=[];
reps:=[Center(G)];
xP:=Concatenation(List(reps,x->Elements(x)));
xP:=SSortedList(xP);

for x in G do
if Order(x)=p and not (x in xP)  then
Add(pPowers,x);
Append(xP,Elements(Group(x)));
fi;
od;


toggle:=true;
while toggle do
toggle:=false;

	for S in AbSubGrps1 do
	for x in pPowers do
	if not x in Group(S) then
	Sx:=Concatenation(S,[x]);
	A:=Group(Sx);
	if IsAbelian(A) and not (A in reps) then
	Add(AbSubGrps2,Sx); AddSet(reps,A); toggle:=true; fi;
	fi;
	od;
	od;
Append(AbSubGrps,AbSubGrps1);
AbSubGrps1:=AbSubGrps2;
AbSubGrps2:=[];

od;


Apply(AbSubGrps, x->Prank(Group(x)));

return Maximum(AbSubGrps);
end);
#####################################################################

####################################################################
#####################################################################
InstallMethod(Compose,
"for group homomorphisms xoy=x(y)",
[IsGroupHomomorphism,IsGroupHomomorphism],
function(x,y)
#if not (Source(x)=Range(y)) then return fail; fi;
return GroupHomomorphismByFunction(Source(y),Range(x),
a->ImagesRepresentative(x,ImagesRepresentative(y,a)));
end);
#####################################################################

####################################################################
#####################################################################
InstallOtherMethod(Compose,
"for chain homomorphisms FoG=F(G)",
[IsHapChainMap,IsHapChainMap],
function(F,G)
local S,T,map,char;
S:=Source(F);
T:=Target(G);
if not S=T then return fail; fi;
char:=EvaluateProperty(F,"characteristic");

map:=function(n,k);
return F!.mapping(n,G!.mapping(n,k));
end;

return Objectify(HapChainMap,
       rec(
           source:=Source(G),
           target:=Target(F),
           mapping:=map,
           properties:=[["characteristic", char],["type","chainMap"]]));

end);
#####################################################################


#####################################################################
InstallOtherMethod(Compose,
"for FpG-module homomorphisms xoy=x(y)",
[IsHapFPGModuleHomomorphism,IsHapFPGModuleHomomorphism],
function(x,y)
return CompositionOfFpGModuleHomomorphisms(x,y);
end);
#####################################################################

#####################################################################
InstallOtherMethod(Size,
"for free ZG-resolutions",
[IsHapResolution],
function(R) local n,i,s,L;
L:=[];

for n in [1..Length(R)] do
s:=0;
for i in [1..R!.dimension(n)] do
s:=s+Length(R!.boundary(n,i));
od;
Add(L,s);
od;
return L;
end);
#####################################################################


####################################################################
#####################################################################
InstallGlobalFunction(PCentre,
function(arg)
local
	G,prime, gens, gensp, C, hom;

G:=arg[1];

if not IsPGroup(G) and Length(arg)=1 then
Print("Error: input must be a p-group or a group and a prime.\n");
return fail;
fi;

if Length(arg)=2 then
prime:=arg[2];
else
prime:=PrimePGroup(G);
fi;

if Order(G)=1 then return G; fi;

C:=Centre(G);
gens:=GeneratorsOfGroup(C);
gensp:=List(gens,x->x^prime);
hom:=GroupHomomorphismByImages(C,C,gens,gensp);

return Kernel(hom);
end);
####################################################################
#####################################################################



####################################################################
#####################################################################
InstallGlobalFunction(PUpperCentralSeries,
function(arg)
local
        G,Q,S,U,P,prime, hom, sz;

G:=arg[1];

if not IsPGroup(G) and Length(arg)=1 then
Print("Error: input must be a p-group or a group and a prime.\n");
return fail;
fi;

if Length(arg)=2 then
prime:=arg[2];
else
prime:=PrimePGroup(G);
fi;


U:=[Group(Identity(G))];
P:=U[Length(U)];
sz:=0;
while Order(P)>sz do
sz:=Order(P);
hom:=NaturalHomomorphismByNormalSubgroup(G,P);
Q:=Image(hom);
P:=PreImagesSet(hom,PCentre(Q,prime));
if Order(P)>sz then
Add(U,P);
fi;
od;

return Reversed(U);
end);
####################################################################
#####################################################################

####################################################################
#####################################################################
InstallGlobalFunction(CanonicalRightCountableCosetElement,
function(U,g)
local x, h;

if IsFinite(U) or IsBound(U!.gname)
then return CanonicalRightCosetElement(U,g); fi;

if g in U then return One(U); fi;

if not IsBound(U!.ccelts) then U!.ccelts:=[]; fi;

h:=g^-1;
for x in U!.ccelts do
if x*h in U then return x; fi;
od;

Add(U!.ccelts,g);
return g;

end);
####################################################################
#####################################################################

####################################################################
#####################################################################
InstallGlobalFunction( SL2Z,
    function(g)
    local p,S, F, T, G,P,gens;

S:=[[0,-1],[1,0]];
T:=[[1,0],[1,1]];

if IsInt(g) then
############################
if (not IsPrime(g)) and (not g=1) then
Print("SL2Z(p) is not yet implemented for non primes p.\n");
return fail; fi;

P:=[[1,0],[0,g]];
G:=Group([P*S*P^-1,P*T*P^-1]);
SetName(G,Concatenation("SL(2,Z)^",String(P))  );
G!.mat:=P;
SetIsHAPRationalMatrixGroup(G,true);
SetIsHAPRationalSpecialLinearGroup(G,true);
    
return G;
############################
fi;

############################
gens:=[S,T];
F:=DenominatorRat(g); 
F:=Factors(F);
if Length(F)>Length(SSortedList(F)) then
Print("SL(Z[1/m]) is implemented only for m a square free integer.\n");
return fail;
fi;

for p in F do
P:=[[1,0],[0,p]];
Add(gens,P*S*P^-1);
Add(gens,P*T*P^-1);
od;
G:=Group(gens);
SetName(G,Concatenation("SL(2,Z[1/",String(DenominatorRat(g)),"])")  );
#G!.primes:=F[1];
G!.primes:=DenominatorRat(g);
SetIsHAPRationalMatrixGroup(G,true);
SetIsHAPRationalSpecialLinearGroup(G,true);

return G;
############################

        end);
####################################################################
#####################################################################


######################
InstallMethod( \in,
               "for SL(2,Z)_p",
              [ IsMatrix,  IsHAPRationalSpecialLinearGroup ],
function ( g, G )
local P,n,d,facs,   m,p,H,K;

if IsBound(G!.mat) then
###############
P:=G!.mat;
return 
P^-1*g*P in SL(2,Integers);
###############
fi;

if IsBound(G!.primes) then
###############
facs:=SSortedList(Factors(G!.primes));
AddSet(facs,1);
for n in Flat(g) do
d:=SSortedList(Factors(DenominatorRat(n)));
if not IsSubset(facs,d) then # and n>1 then
return false; fi;
od;

if  Determinant(g)=1 then return true; fi;
return false; 
###############
fi;


########
if IsBound(G!.coprimes) then
m:=G!.coprimes[1];
p:=G!.coprimes[2];
P:=[[1,0],[0,p]];
return P^-1*g*P in SL2Z(1/m);
fi;
###############
if IsBound(G!.levels) then
m:=G!.levels[1];
p:=G!.levels[2];
H:=SL2Z(1/m);
K:=ConjugateSL2ZGroup(H,[[1,0],[0,p]]);
return (g in H and g in K);
fi;
##############################

end );
######################


if IsPackageMarkedForLoading("congruence","0.0") then
######################
InstallMethod( \in,
               "for CongruenceSubgroupGamma0(p)",
              [ IsMatrix,  IsCongruenceSubgroupGamma0 ],
function ( g, G )
local p;
p:=G!.LevelOfCongruenceSubgroup;
if g in SL(2,Integers) then
if Determinant(g)=1 then
if IsInt(g[2][1]/p) then return true;fi;
fi;
fi;
############################
return false;
end );
fi;


############################################
InstallGlobalFunction(KernelWG,
function(phi)
local K;

K:=Kernel(phi);
if Length(GeneratorsOfGroup(K))=0 then
K:=Group(One(K));
fi;
return K;
end);
############################################

############################################
InstallGlobalFunction(ScatterPlot,
function(arg)
local LL, L, tmpdir, file, colour, i, x, xmax, xmin, ymax, ymin, xscale, yscale,
xminnew, yminnew, xmaxnew, ymaxnew, AppendTo, PrintTo;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

LL:=arg[1];

tmpdir := DirectoryTemporary();;
file:=Filename( tmpdir , "tmp.asy" );


xmax:=Maximum(List(LL,x->x[1]));
xmin:=Minimum(List(LL,x->x[1]));
ymax:=Maximum(List(LL,x->x[2]));
ymin:=Minimum(List(LL,x->x[2]));

if Length(arg)=2 then
xmax:=Maximum(xmax,ymax);
ymax:=xmax;
xmin:=Minimum(xmin,ymin);
ymin:=xmin;
fi;

xscale:=200/(xmax-xmin);
yscale:=200/(ymax-ymin);

L:=List(LL,x->[x[1]-xmin,x[2]-ymin]);
L:=List(L,x->[x[1]*xscale,x[2]*yscale]);
for i in [1..Length(LL)] do
if Length(LL[i])=3 then Add(L[i],LL[i][3]); fi;
od;

xminnew:=0;
yminnew:=0;
xmaxnew:=200;
ymaxnew:=200;

PrintTo(file,"import math; size(400,400);");

AppendTo(file,"draw((", xminnew, ",", yminnew,")--(",xmaxnew,",",yminnew,"),blue+linewidth(0.5));");

AppendTo(file,"draw((", xminnew, ",", yminnew,")--(",xminnew,",",ymaxnew,"),blue+linewidth(0.5));");

AppendTo(file,"label(\"$" , xmin , "$\" , (" , xminnew , ",", yminnew-20, "),blue);");

AppendTo(file,"label(\"$" , xmax , "$\" , (" , xmaxnew , ",", yminnew-20, "),blue);");

AppendTo(file,"label(\"$" , ymin , "$\" , (" , xminnew-20 , ",", yminnew, "),blue);");

AppendTo(file,"label(\"$" , ymax , "$\" , (" , xminnew-20 , ",", ymaxnew, "),blue);");


for x in L do
if Length(x)=3 then  colour:=x[3]; else colour:="blue"; fi;
AppendTo(file, "dot((", x[1], ",", x[2], "),",colour, "+linewidth(4));");
od;



Exec( Concatenation( ASY_PATH,"-V ", file) );

RemoveFile(file);
file:=Filename(tmpdir,"");
RemoveFile(file);

end);
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(PushoutOfFpGroups,
function(f,g)
local P, rels, FhomP, FFhomP, GhomP, GGhomP, FGhomP, F, FF, G, GG, FG,
gensF, gensFF, gensG, gensGG, gensP, x, gg;

#if not Source(f)=Source(g) then return fail; fi;

F:=Target(f);
FF:=FreeGroupOfFpGroup(F);
G:=Target(g);
GG:=FreeGroupOfFpGroup(G);
FG:=Source(f);

gensFF:=GeneratorsOfGroup(FF);
gensF:=GeneratorsOfGroup(F);
gensGG:=GeneratorsOfGroup(GG);
gensG:=GeneratorsOfGroup(G);
P:=FreeGroup(Length(gensFF)+Length(gensGG));
gensP:=GeneratorsOfGroup(P);

FFhomP:=GroupHomomorphismByImagesNC(FF,P,gensFF, gensP{[1..Length(gensFF)]});
FhomP:=GroupHomomorphismByImagesNC(F,P,gensF, gensP{[1..Length(gensFF)]});
GGhomP:=GroupHomomorphismByImagesNC(GG,P,gensGG, gensP{[1+Length(gensFF)..Length(gensP)]});
GhomP:=GroupHomomorphismByImagesNC(G,P,gensG, gensP{[1+Length(gensFF)..Length(gensP)]});

gg:=GroupHomomorphismByImagesNC(FG,G,GeneratorsOfGroup(FG),
List(GeneratorsOfGroup(Source(g)),x->Image(g,x)) );

rels:=[];
Append(rels, List(RelatorsOfFpGroup(F),x->Image(FFhomP,x)) );
Append(rels, List(RelatorsOfFpGroup(G),x->Image(GGhomP,x)) );

for x in GeneratorsOfGroup(FG) do
Add(rels,
Image(FhomP, Image(f,x))^-1 * Image(GhomP, Image(gg,x))
);
od;


return P/rels;

end);
##########################################################################
##########################################################################

##########################################################################
##########################################################################
InstallGlobalFunction(Fp2PcpAbelianGroupHomomorphism,
function(F)
local P, S, T, pS, pT, ShompS, pShomS, pThomT, ThompT, gensS,  genspS;

S:=Source(F);
T:=Target(F);

ShompS:=HAP_NqEpimorphismNilpotentQuotient(S,1);
pS:=Range(ShompS);
genspS:=GeneratorsOfGroup(pS);
gensS:=List(genspS, x->PreImagesRepresentative(ShompS,x));
pShomS:=GroupHomomorphismByImagesNC(pS,S,genspS,gensS);

ThompT:=HAP_NqEpimorphismNilpotentQuotient(T,1);
pT:=Range(ThompT);

P:=GroupHomomorphismByFunction(pS,pT,
x->Image(ThompT,Image(F, Image(pShomS,x) ) ) );

return P;
end);
##########################################################################
##########################################################################

##########################################################################
##########################################################################
InstallGlobalFunction(IsIsomorphismOfAbelianFpGroups,
function(F)
local h;

## I've adjusted this so that it will also work for Pc and Pcp groups.
if IsFpGroup(Source(F)) and IsFpGroup(Target(F)) then
h:=Fp2PcpAbelianGroupHomomorphism(F);
else h:=F;
fi;

if not IsTrivial(Kernel(h)) then return false; fi;
if not IsTrivial(Range(h)/Image(h)) then return false; fi;

return true;
end);
##########################################################################
##########################################################################

##########################################################################
##########################################################################
InstallGlobalFunction(SignedPermutationGroup,
function(n)
local A, i, gens;

gens:=[];;

for i in [1..n-1] do
Add(gens,  PermutationMat((i,i+1),n));
od;

for i in [1..n] do
A:=1*IdentityMat(n);
A[i]:=-A[i];
Add(gens,A);
od;

return Group(gens);
end);
##########################################################################
##########################################################################

##########################################################################
##########################################################################
InstallMethod(IsPeriodic,
"Tests to see if a finite group has periodic cohomology",
[IsGroup and IsFinite],
function(G)
local p,S;
if Order(G)=1 then return true; fi;
for p in SSortedList(Factors(Order(G))) do
S:=SylowSubgroup(G,p);
if not (IsCyclic(S) or IsQuaternionGroup(S))  then return false; fi;
od;
return true;
end);
##########################################################################
##########################################################################

##########################################################################
##########################################################################
InstallMethod(CohomologicalPeriod,
"Returns the period of a finite group with periodic cohomology",
[IsGroup and IsFinite],
function(G)
local p,S,L,N,C;
if not IsPeriodic(G) then
Print("The group does not have periodic cohomology.\n");
return fail; fi;

L:=[];
for p in SSortedList(Factors(Order(G))) do
S:=SylowSubgroup(G,p);
if p=2 then
if IsCyclic(S) then Add(L,2); else Add(L,4); fi;
else
N:=Normalizer(G,S);
C:=Centralizer(G,S);
Add(L,2*Order(N)/Order(C));
fi;
od;
return Lcm(L);
end);
##########################################################################
##########################################################################

###############################################
###############################################
InstallGlobalFunction(SignatureOfSymmetricMatrix,
function(A)
local cp,f, zer,pos,neg,i,j;

if not IsMatrix(A) then
Print("The argument must be a matrix.\n");
return fail;
fi;

if not A=TransposedMat(A) then
Print("The argument must be a symmetric matrix.\n");
return fail;
fi;

for i in [1..Length(A)] do
for j in [1..Length(A[1])] do
if not IsRat(A[i][j]) then
Print("The argument must be a symmetric matrix of rational numbers.\n");
return fail;
fi;
od;od;

cp:=CharacteristicPolynomial(A);
f:=CoefficientsOfUnivariatePolynomial(cp);
f:=List(f,SignInt);

zer:=PositionProperty(f,x->x<>0);
if zer=fail then zer:=Length(f)-1;
else zer:=zer-1; fi;

f:=Filtered(f,a->not IsZero(a));
f:=List([1..Length(f)-1],i->f[i]*f[i+1]);

pos:= Length(Filtered(f,IsNegInt));
neg:= Degree(cp)-pos-zer;
return rec( zero_eigenvalues:=zer, positive_eigenvalues:=pos, negative_eigenvalues:=neg, determinant:=Determinant(A) );
end);
###############################################
###############################################
