#####################################################################
#####################################################################
InstallGlobalFunction(SimplifiedQuandlePresentation,
function(P)
local smplfy, bool, Q,subst,rename;

############################
rename:=function(s,n)
local w,x,i;
w:=1*s;
for i in [1,2] do
if IsInt(w[i]) then
if w[i]>n then w[i]:=w[i]-1; fi;
else
w[i]:=rename(w[i],n);
fi;
od;
return w;
end;
############################
############################
subst:=function(s,r)
local i,w;
# r=[i,w] with i an integer and w a list or integer.
if IsInt(s[1]) then
if s[1]=r[1] then s[1]:=r[2]; fi;
else s[1]:=subst(s[1],r);
fi;
if IsInt(s[2]) then
if s[2]=r[1] then s[2]:=r[2]; fi;
else s[2]:=subst(s[2],r);
fi;

#s:=rename(s,r[1]);
return s;
end;
############################

############################
smplfy:=function(P)
local r,rels,i;
bool:=false;
for i in [1..Length(P!.relators)] do
r:=P!.relators[i];

if IsInt(r[2]) then P!.relators[i]:=[r[2],r[1]]; r:=P!.relators[i];fi;
if IsInt(r[1]) and not r[1] in Flat(r[2]) then  bool:=true; break; fi;
od;
if bool then

rels:=Filtered(P!.relators,s->not s=r);

Apply(rels,s->subst(s,r));
Apply(rels,s->rename(s,r[1]));


P!.relators:=rels;
P!.generators:=[1..Length(P!.generators)-1];
fi;
return P;
end;
############################

bool:=true;
Q:=Objectify(HapQuandlePresentation,rec( relators:=1*P!.relators,
generators:=1*P!.generators));
while bool do
Q:=smplfy(Q);
od;
return Q;
end);
#####################################################################
#####################################################################



#####################################################
#####################################################
InstallGlobalFunction(TupleOrbitReps_perm,
function(A,Q,n,tester)
local orbs, neworbs, i, x, q, y, D, S, B;

orbs:=Orbits(A,Q);
orbs:=List(orbs,x->[[x[1]],Stabilizer(A,x[1])]);

for i in [2..n] do
neworbs:=[];
for x in orbs do
S:=x[2];
D:=[];

for q in Q do
y:=Concatenation(x[1],[q]);
Add(D,y);
od;

B:=Orbits(A,D,OnTuples);
B:=List(B,x->[x[1],Stabilizer(S,x[1],OnTuples)]);
Append(neworbs,B);

od;
orbs:=neworbs;
orbs:=Filtered(orbs,x->tester(x[1]));
od;

orbs:=List(orbs,x->x[1]);
return orbs;
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallGlobalFunction(TupleOrbitReps,
function(A,Q,n,tester)
local gens,N,permgens,elts,Aperm,a,p,orbs;

if IsPermGroup(A) then return TupleOrbitReps_perm(A,Q,n,tester); fi;


elts:=Elements(Q);
gens:=GeneratorsOfGroup(A );
N:=Size(Q);
permgens:=[];
for a in gens do
p:=List([1..N],i->Position(elts,elts[i]^a));
p:=PermList(p);
Add(permgens,p);
od;

Aperm:=Group(permgens);
permgens:=SmallGeneratingSet(Aperm);
Aperm:=Group(permgens);

orbs:=TupleOrbitReps_perm(Aperm,[1..N],n,tester);
return [orbs,Aperm];
end);
#####################################################
#####################################################

#####################################################################
#####################################################################
InstallGlobalFunction(NumberOfHomomorphisms_connected,
#function(genRelK,Q)
function(arg)
local genRelK,Q,P,multTab,tester,A,T,nbGen,rels,RELS,r,n, wrd;

genRelK:=arg[1];
Q:=arg[2];
multTab:=MultiplicationTable(Q);
P:=SimplifiedQuandlePresentation(genRelK);
nbGen:=Length(P!.generators);
rels:=P!.relators;
RELS:=[];
for n in [1..nbGen] do
RELS[n]:=[];
od;
for r in rels do
Add(RELS[Maximum(Flat(r))], r);
od;

###########################
wrd:=function(mapping,rr);
if not IsList(rr) then return mapping[rr]; fi;
return multTab[wrd(mapping,rr[1])][wrd(mapping,rr[2])];
end;
###########################

###########################
tester:=function(mapping)   #no longer used!
local R,re,x,y,z;
R:=RELS[Length(mapping)];
for re in R do
Print(re[1],"   ",re[2]," ");
Print(wrd(mapping,re[1]),"   ",wrd(mapping,re[2])," ");
x:=re[1][1];
y:=re[1][2];
z:=re[2];
if not mapping[z]=multTab[mapping[x]][mapping[y]] then return false; fi;
od;
return true;
end;
###########################

###########################
tester:=function(mapping)
local R,re;
R:=RELS[Length(mapping)];
for re in R do
if not wrd(mapping,re[1])=wrd(mapping,re[2]) then return false; fi;
od;
return true;
end;
###########################


#A:=AutomorphismGroupQuandle(Q);
A:=RightMultiplicationGroupOfQuandle(Q);
T:=TupleOrbitReps(A,Q,nbGen,tester);
if Length(arg)=3 then
return List(T[1],x->Magma(List(x,j->Elements(Q)[j])));
fi;
T:=List(T[1],x->Size(Orbit(T[2],x,OnTuples)));
return T;

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallMethod(HomomorphismsImages,
"for an fp quandle and finite connected quandle",
[IsHapQuandlePresentation,IsMagma],
function(K,Q)
return NumberOfHomomorphisms_connected(K,Q,true);
end);


#####################################################################
#####################################################################
InstallOtherMethod(NumberOfHomomorphisms,
"for an fp group and a finite group",
[IsGroup,IsGroup],
function(KK,Q)
local GG,F, K, gens, tester,A,T,elts,r,i,w,ww;

#Q is a finite group and K is an fp group
elts:=Elements(Q);
K:=Image(IsomorphismFpGroup(KK));
F:=FreeGroupOfFpGroup(K);
gens:=GeneratorsOfGroup(F);
GG:=[];
for i in [1..Length(gens)] do
GG[i]:=[];
for r in RelatorsOfFpGroup(K) do
w:=ExtRepOfObj(r);
ww:=w{2*[1..Length(w)/2]-1};
if Maximum(ww)<=i then
Add(GG[i],w); fi;
od;
od;

###########################
tester:=function(mapping)
local w, g;
for w in GG[Length(mapping)] do
g:=One(Q);
for i in [1..Length(w)/2] do
g:=g*elts[mapping[w[2*i-1]]]^w[2*i];  
od;
if not Order(g)=1 then return false; fi;
od;
return true;
end;
###########################


A:=AutomorphismGroup(Q);
T:=TupleOrbitReps(A,Q,Length(gens),tester);
T:=List(T[1],x->Size(Orbit(T[2],x,OnTuples)));
return Sum(T);

end);
#####################################################################

