#############################################
#############################################
CommutingProbability:=function(G);
return Length(ConjugacyClasses(G))/Order(G);
end;
#############################################
#############################################

#########################################################
#########################################################
InstallGlobalFunction(BogomolovMultiplier,
function(arg)
local G, A, RA, RG, ChG, RS, CC, Pairs, gensA, H2, gensH2, 
      Q, GhomQ, RQ, ChQ, OrdImH2G, OrdH2Q, HH2, AA, PairsQ,
      g, h, x, M, n, V, P, p, bool, cnt, toggle;

G:=arg[1];
if Length(arg)=1 then toggle:="standard"; else
toggle:= arg[2]; fi;

################### Very easy cases  ###############
###################                  ###############
if IsAbelian(G) then G!.BogomolovMultiplier:=[]; return []; fi;
if CommutingProbability(G)>1/4 then G!.BogomolovMultiplier:=[]; return []; fi;
################### Very easy cases  ###############
################### finished         ###############


if not toggle in ["standard", "homology", "tensor"] then
Print("The available algorithms are \"standard\", \"homology\" and \"tensor\".\n");
return fail; fi;

if IsBound(G!.BogomolovMultiplier) then
return G!.BogomolovMultiplier; fi;


if not IsFinite(G) then
Print("At present this function is only implemented for finite groups.\n");
return fail;
fi; 

################### Non p-groups ##################
###################              ##################
if not IsPGroup(G) then
bool:=true;
P:=SSortedList(Factors(Order(G)));

for p in P do
if Length(BogomolovMultiplier(SylowSubgroup(G,p),toggle))>0 then
bool:=false; break; fi;
od;

if bool then  G!.BogomolovMultiplier:=[]; return []; fi;

##for bool=false we should implement a Cartan-Eilenberg double
##coset formula to significantly speed up the computation.

fi;
################### Non p-groups ##################
################### finished     ##################

if toggle="homology" then return Bogomology(G,2); fi;

if toggle="homology" then return BogomolovMultiplier_viaTensorSquare(G); fi;

RA:=ResolutionAbelianGroup_alt([0,0],3);
A:=RA!.group;
gensA:=GeneratorsOfGroup(A);

CC:=ConjugacyClasses(G);;

#if Length(CC)/Size(G) <150/3125 then return Bogomology(G,2); fi;

CC:=List(CC,x->Representative(x));
Pairs:=[];

for g in CC do
if not IsOne(g) then
M:=MinimalGeneratingSet(Centralizer(G,g));
for h in M do
Add(Pairs,[g,h]);
od;
fi;
od;


GhomQ:=NaturalHomomorphismByNormalSubgroup(G,Centre(G));
Q:=Range(GhomQ);

if IsPGroup(Q) then
#RQ:=ResolutionGenericGroup(Q,3);
#RQ:=ResolutionSubnormalSeries(CompositionSeries(Q),3);
RQ:=ResolutionNilpotentGroup(Q,3);
  P:=Homology(TensorWithIntegers(RQ),2);
else
  P:=GroupHomology(Q,2);
fi;

AA:=Intersection(Center(G),DerivedSubgroup(G));
OrdH2Q:=Product(P);
OrdImH2G:=OrdH2Q/Order(AA);
if OrdImH2G=0 then  G!.BogomolovMultiplier:=[]; return []; fi;

if IsPGroup(G) then
################### An easier case ###############
################### for p-groups   ###############
ChQ:=TensorWithIntegers(RQ);
H2:=TransposedMat(BoundaryMatrix(ChQ,3));
H2:=List(H2,x->x);
PairsQ:=List(Pairs,x-> [Image(GhomQ,x[1]),Image(GhomQ,x[2])]);
PairsQ:=SSortedList(PairsQ);
PairsQ:=Filtered(PairsQ,x-> not ( IsOne(x[1])  or IsOne(x[2]) ));
PairsQ:=Filtered(PairsQ,x-> not x[1]=x[2]);



V:=[1..ChQ!.dimension(2)]*0;

cnt:=0;
for x in PairsQ do
cnt:=cnt+1;
V:=V*0;
h:=GroupHomomorphismByImagesNC(A,Q,gensA,x);
h:=EquivariantChainMap(RA,RQ,h);
h:=h!.mapping([[1,1]],2);
h:=List(h,w->w[1]);
for n in h do
V[AbsInt(n)]:=V[AbsInt(n)]+SignInt(n);
od;
Add(H2,V);

if (cnt mod 10)=0 then
HH2:=SmithNormalFormIntegerMat(H2);
HH2:=List([1..Length(H2[1])],i->HH2[i][i]);
HH2:=Filtered(HH2,i->not i=1);
HH2:=Filtered(HH2,i->not i=0);
if 
OrdImH2G/( OrdH2Q/Product( HH2 )) =1 then  
G!.BogomolovMultiplier:=[]; return []; fi;
fi;

od;

H2:=SmithNormalFormIntegerMat(H2);
H2:=List([1..Length(H2[1])],i->H2[i][i]);
H2:=Filtered(H2,i->not i=1);
H2:=Filtered(H2,i->not i=0);

x:=OrdImH2G/(OrdH2Q/Product(H2));
if x =1 then G!.BogomolovMultiplier:=[]; return []; fi;
if IsPrime(x) then G!.BogomolovMultiplier:=[x]; return [x]; fi;
################### Easier case     ###############
################### finished       ###############
fi;


RG:=ResolutionGenericGroup(G,3);
ChG:=TensorWithIntegers(RG);
H2:=TransposedMat(BoundaryMatrix(ChG,3));
H2:=List(H2,x->x);

V:=[1..ChG!.dimension(2)]*0;

cnt:=0;
for x in Pairs do
cnt:=cnt+1;
V:=V*0;
h:=GroupHomomorphismByImagesNC(A,G,gensA,x);
h:=EquivariantChainMap(RA,RG,h);
h:=h!.mapping([[1,1]],2);
h:=List(h,w->w[1]);
for n in h do
V[AbsInt(n)]:=V[AbsInt(n)]+SignInt(n);
od;
Add(H2,V);

if (cnt mod 50)=0 and Filtered(SmithNormalFormIntegerMat(H2), a->Sum(a)>0)=[] 
   then  G!.BogomolovMultiplier:=[]; return [];
fi;

od;

H2:=SmithNormalFormIntegerMat(H2);
H2:=List([1..Length(H2[1])],i->H2[i][i]);
H2:=Filtered(H2,i->not i=1);
H2:=Filtered(H2,i->not i=0);

G!.BogomolovMultiplier:=H2;
return H2;
end);
###########################################################
###########################################################

#########################################################
#########################################################
InstallGlobalFunction(Bogomology,
function(G,N)
local A, RA, RG, ChG, Chh, AbSubs, HN, gensHN,
      g, h, hh, k, x, M, n, V, P, p, B, bool;

if not IsFinite(G) then
Print("At present this function is only implemented for finite groups.\n");
return fail;
fi;

if IsAbelian(G) then return []; fi;

#RG:=ResolutionGenericGroup(G,N+1);
RG:=ResolutionNilpotentGroup(G,N+1);
ChG:=TensorWithIntegers(RG);
HN:=TransposedMat(BoundaryMatrix(ChG,N+1));
HN:=List(HN,x->x);

##Probably a clumsy method for finding all maximal abelian subgroups
##
AbSubs:=ConjugacyClassesSubgroups(G);
AbSubs:=Filtered(AbSubs,x->IsAbelian(Representative(x)));
AbSubs:=List(AbSubs,a->Representative(a));
AbSubs:=List(AbSubs,b->ConjugateSubgroups(G,b));
AbSubs:=Concatenation(AbSubs);
AbSubs:=Filtered(AbSubs,x->Order(x)>1);
A:=[];
for g in AbSubs do
bool:=true;
for x in AbSubs do
if not g=x and IsSubgroup(x,g) then bool:=false; break;fi;
od;
if bool then Add(A,g); fi;
od;
AbSubs:=A;
##
##Clumsy method finished


for A in AbSubs do
h:=GroupHomomorphismByImagesNC(A,G,GeneratorsOfGroup(A),GeneratorsOfGroup(A));
RA:=ResolutionGenericGroup(A,N);
h:=EquivariantChainMap(RA,RG,h);
Chh:=TensorWithIntegers(h);
V:=[1..RA!.dimension(N)]*0;
B:=TransposedMat(BoundaryMatrix(TensorWithIntegers(RA),N));
B:=NullspaceIntMat(B);

for V in B do
hh:=Chh!.mapping(V,N);
Add(HN,hh);
od;

od;

HN:=SmithNormalFormIntegerMat(HN);
HN:=List([1..Length(HN[1])],i->HN[i][i]);
HN:=Filtered(HN,i->not i=1);
HN:=Filtered(HN,i->not i=0);
return HN;

end);
###########################################################
###########################################################

#########################################################
#########################################################
InstallGlobalFunction(AreIsoclinic,
function(arg)
local G,H,bool, GhomGZ, HhomHZ, GZ, gensGZ, isotest, Aut, Inn, OutReps,
      HZ, DG, DH, Homs, HOMS, gensDG, gensDG1, gensDH, JG,
      CG, CH, L, f, ff, g, h, i, gg,hh, x, xx, iso, iso1, gensG, TotalDegree;
#This function can return true or false.
#The function does NOT test isomorphism. So sometimes an 
#isomorphism test first will speed things up. 
#We assume that input groups have derived subgroup and central quotient
#of size <1024 and not equal to 512.

G:=arg[1];
H:=arg[2];
if Length(arg)>2 then bool:=arg[3]; else bool:=true; fi;

######################Quick positive tests########################
if G=H then return true; fi;
if IsAbelian(G) and IsAbelian(H) then return true; fi;
######################Quick positive tests done###################

######################Quick negative tests########################
DG:=DerivedSubgroup(G);
DH:=DerivedSubgroup(H);
if not IdGroup(DG)=IdGroup(DH) then return false; fi;
######################Quick negative tests done###################

GhomGZ:=NaturalHomomorphismByNormalSubgroup(G,Center(G));
GZ:=Image(GhomGZ);

HhomHZ:=NaturalHomomorphismByNormalSubgroup(H,Center(H));
HZ:=Image(HhomHZ);


######################More quick negative tests##########
if not IdGroup(GZ)=IdGroup(HZ) then return false; fi;

if bool then

CG:=Collected(List(ConjugacyClasses(G), Size));
CG:=List(CG,x->[x[1],x[2]/Order(G)]);
CH:=Collected(List(ConjugacyClasses(G), Size));
CH:=List(CH,x->[x[1],x[2]/Order(H)]);
if not CG=CH then return false; fi;

if not CommutingProbability(G)=CommutingProbability(H) then return false; fi;

TotalDegree:= G -> Sum(CharacterDegrees(G), Product)/Order(G);

if not TotalDegree(G)=TotalDegree(H) then return false; fi;

fi;
######################More quick negative tests done#####

if not bool then return false; fi;  #So in this case some isoclinic pairs
                                    #may be deemed non-isoclinic
                                    #This line actually speed up 
                                    #IsoclinismClasses()


######################Another quick positive test########
gensDG1:=StrongGeneratorsOfDerivedSubgroup(G); #this is time consuming!
gensDG:=List(gensDG1[1],x->x[1]);
JG:=gensDG1[3];
gensG:=GeneratorsOfGroup(G);

######################################
######################################
isotest:=function(iso)
local ff, gensDH, iso1,T, imG, a, b, x;

imG:=List(gensG,x->Image(iso,Image(GhomGZ,x)));
ff:=GroupHomomorphismByImagesNC(G,HZ,gensG,imG);

T:=gensDG1[2];

gensDH:=List(gensDG1[1],x->
[PreImagesRepresentative(HhomHZ,Image(ff,x[2][1])),
PreImagesRepresentative(HhomHZ,Image(ff,x[2][2]))]
);

Apply(gensDH,x->x[1]^-1*x[2]^-1*x[1]*x[2]);

#####Case when Bogomolov Multiplier is trivial#####
if T=true then
if SSortedList(gensDH)=[One(DH)] then return  true; else return false; fi;
fi;
#####Special case done#############################

iso1:=GroupHomomorphismByImagesNC(T,DH,gensDG,gensDH);

for x in JG do
if not Image(iso1,x)=One(DH) then return false; fi;
od;

return true;
end;
######################################
######################################

iso:=IsomorphismGroups(GZ,HZ);

if isotest(iso) then return true; fi;
######################Another quick positive test done#####


if not bool then return false; fi;  #So in this case some isoclinic pairs
                                    #may be deemed non-isoclinic


if not IsBound(H!.OutReps) then
######################################
Aut:=AutomorphismGroup(HZ);
Inn:=InnerAutomorphismsAutomorphismGroup(Aut);
OutReps:=RightCosets(Aut,Inn);
OutReps:=List(OutReps,x->Representative(x));
H!.OutReps:=OutReps;
######################################
fi;

OutReps:=H!.OutReps;

for f in OutReps do
iso1:=GroupHomomorphismByFunction(GZ,HZ,x->Image(f,Image(iso,x)));
if isotest(iso1) then return true; fi;
od;

return false;
end);
###########################################################
###########################################################

###########################################################
###########################################################
InstallGlobalFunction(PartialIsoclinismClasses,
function(L)
local test;

test:=function(G,H); return AreIsoclinic(G,H,false); end;

return HAP_EquivalenceClasses(L,test);

end);
###########################################################
###########################################################

###########################################################
###########################################################
InstallGlobalFunction(IsoclinismClasses,
function(L)
local inv, C, D, DD, G,H,P,TotalDegree,LL,A,M;
#We assume that each group has derived subgroup and central quotient
#of order <1024 and not equal to 512.

TotalDegree:= G -> Sum(CharacterDegrees(G), Product)/Order(G);

########
inv:=function(G)
local CG;
CG:=Collected(List(ConjugacyClasses(G), Size));
CG:=List(CG,x->[x[1],x[2]/Order(G)]);
return [IdGroup(DerivedSubgroup(G)),
        IdGroup(G/Center(G)),
        CG,
        NilpotencyClassOfGroup(G),
        CommutingProbability(G),
        TotalDegree(G)];
end;
########

LL:=Classify(L,inv);
A:=[];
        
for M in LL do
#####################################
DD:=PartialIsoclinismClasses(M);;
D:=List(DD,d->d[1]);

C:=[];

while Length(D)>0 do
G:=D[1];
P:=[];
for H in D do
if AreIsoclinic(G,H,true) then Add(P,H); fi;
od;
D:=Filtered(D,x-> not x in P);
P:=Filtered(DD,x->Length(Filtered(x,i->i in P)) >0);
P:=Concatenation(P);
Add(C,P);
od;

Append(A,C);
#####################################
od;

return A;
end);
###########################################################
###########################################################


#########################################################
#########################################################
InstallGlobalFunction(StrongGeneratorsOfDerivedSubgroup,
function(G)
local DG, sgensDG, elts, T, TS, TStmp, delta, pairing,
      GhomQ, Q, g, h, p, i, j, bool, JG, DELTA, gens, imgens;

if not IsBound(G!.StrongGeneratorsOfDerivedSubgroup) then

if Order(G)/Order(Center(G)) > 81 then
p:=StrongGeneratorsOfDerivedSubgroup_alt(G);
if not p=fail then return p; fi;
fi;

################################################
GhomQ:=NaturalHomomorphismByNormalSubgroup(G,Center(G));
Q:=Image(GhomQ);
if IsPolycyclicGroup(Q) and not IsNilpotentGroup(Q) then
T:=NonabelianTensorSquare_inf(Q);   ##Changed 19 August 2017
else
T:=NonabelianTensorSquare(Q);
fi;
delta:=T.homomorphism;
pairing:=T.pairing;
TS:=Source(delta);
TStmp:=Group(One(TS));

sgensDG:=[];
bool:=false;
for g in Q do
for h in Q do
p:=pairing(g,h);
if not p in TStmp then
Add(sgensDG,[p,[g,h]]);
TStmp:=Group(List(sgensDG,x->x[1]));
if Order(TStmp)=Order(TS) then bool:=true; break; fi;
fi;
od;
if bool then break; fi;
od;

TStmp:=List(sgensDG,x->x[1]);
TStmp:=ReduceGenerators_alt(TStmp,Group(TStmp));

sgensDG:=Filtered(sgensDG,x->x[1] in TStmp);

sgensDG:=List(sgensDG,x->[x[1],[PreImagesRepresentative(GhomQ,x[2][1]),
PreImagesRepresentative(GhomQ,x[2][2])]]);
gens:=List(sgensDG,x->x[2]);
imgens:=List(gens,x->Comm(x[1],x[2]));
DELTA:=GroupHomomorphismByImagesNC(TS,G,TStmp,imgens);
JG:=Kernel(DELTA);
JG:=ReduceGenerators_alt(GeneratorsOfGroup(JG),JG);

G!.StrongGeneratorsOfDerivedSubgroup:= [sgensDG,Source(delta),JG];
##########################
fi;

return G!.StrongGeneratorsOfDerivedSubgroup;
end);
########################################################
########################################################

########################################################
########################################################
InstallGlobalFunction(StrongGeneratorsOfDerivedSubgroup_alt,
function(G)
local F, g, h, CC, Pairs, M, T, i1, i2, JG;

if IsBound(G!.StrongGeneratorsOfDerivedSubgroup) then
return G!.StrongGeneratorsOfDerivedSubgroup; fi;

if Length(BogomolovMultiplier(G))>0 then 
Print("Using nonabelian tensor square on small group ", IdGroup(G)," This can be slow. \n");
return fail; fi;

F:=FreeProduct(G,G);
i1:=Embedding(F,1);
i2:=Embedding(F,2);

CC:=ConjugacyClasses(G);;
CC:=List(CC,x->Representative(x));
Pairs:=[];

for g in CC do
if not IsOne(g) then
M:=MinimalGeneratingSet(Centralizer(G,g));
for h in M do
Add(Pairs,[Comm(Image(i1,g),Image(i2,g)) ,[g,h]]);
od;
fi;
od;

T:=true;
JG:=true;

G!.StrongGeneratorsOfDerivedSubgroup:=[Pairs,T,JG];
return [Pairs,T,JG];
end);
########################################################
########################################################

##########################################################
##########################################################
InstallGlobalFunction(BogomolovMultiplier_viaTensorSquare,
function(G)
local GhomQ, Q, x,y, T, TS, TStmp, pairing, CC, Pairs, JG,
      p, g, h, M, K, gens, gens1, hom, bool;

GhomQ:=NaturalHomomorphismByNormalSubgroup(G,Center(G));
Q:=Image(GhomQ);

T:=NonabelianTensorSquare(Q);
TS:=Source(T.homomorphism);
pairing:=T.pairing;

CC:=ConjugacyClasses(G);;
CC:=List(CC,x->Representative(x));
Pairs:=[];

for g in CC do
if not IsOne(g) then
M:=MinimalGeneratingSet(Centralizer(G,g));
for h in M do
Add(Pairs,[g,h]);
od;
fi;
od;

Pairs:=List(Pairs,x-> [Image(GhomQ,x[1]),Image(GhomQ,x[2])]);
Pairs:=SSortedList(Pairs);
Pairs:=Filtered(Pairs,x-> not ( IsOne(x[1])  or IsOne(x[2]) ));

K:=List(Pairs, x->pairing(x[1],x[2]));
for g in GeneratorsOfGroup(Q) do
Add(K,pairing(g,g));
for h in GeneratorsOfGroup(Q) do
Add(K,pairing(g,h)*pairing(h,g));
od;
od;
K:=Group(K);
K:=NormalClosure(Parent(TS),K);

gens:=[];
TStmp:=Group(One(Q));

bool:=false;
for g in Q do
for h in Q do
p:=pairing(g,h);
if not p in TStmp then
Add(gens,[p,[g,h]]); TStmp:=Group(List(gens,x->x[1]));
if Order(TStmp)=Order(TS) then bool:=true; break; fi;
fi;
od;
if bool then break; fi;
od;

gens1:=List(gens,x->x[1]);;
gens1:=ReduceGenerators_alt(gens1,TStmp);
gens:=Filtered(gens,x->x[1] in gens1);
gens1:=List(gens,x->x[1]);
gens:=List(gens,x->Comm(PreImagesRepresentative(GhomQ,x[2][1]),
PreImagesRepresentative(GhomQ,x[2][2])));

hom:=GroupHomomorphismByImagesNC(TS,G,gens1,gens);
JG:=Kernel(hom);

return AbelianInvariants(JG/K);
end);
############################################################
############################################################

