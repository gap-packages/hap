####################################################
####################################################
InstallGlobalFunction(PSubgroupSimplicialComplex,
function(arg)
local
        G,p,filt,P,SubsCl, Subs, fn, cl, MaxSubs, bool, MaxSimps, K,  k,
        s,t,x,y,m,mm,tmp;

G:=arg[1];
p:=arg[2];
if Length(arg)=3 then filt:=arg[3];
else
filt:=function(G); return true; end;
fi;
######################################
if not IsPrimeInt(p) then
Print("Second variable is not a prime.\n");
return fail;
fi;
######################################

P:=SylowSubgroup(G,p);                    #CHANGED Jan 2023
SubsCl:=ConjugacyClassesSubgroups(LatticeByCyclicExtension(P,
      filt,true));;

SubsCl:=List(SubsCl,c->ClassElementLattice(c,1));
SubsCl:=Filtered(SubsCl,cl->Order(cl)>1);
SubsCl:=List(SubsCl,g->g^G);

Subs:=[];
for cl in SubsCl do
for x in [1..Size(cl)] do
Add(Subs,ClassElementLattice(cl,x));
od;
od;

Unbind(SubsCl);

###################
fn:=function(A,B);
return Order(A)>=Order(B);
end;
###################

Sort(Subs,fn);

MaxSimps:=[];
for s in Subs do
bool:=true;
for t in MaxSimps do
if IsSubgroup(t,s) then bool:=false; break; fi;
od;
if bool then Add(MaxSimps,s); fi;
od;

Unbind(Subs);
MaxSimps:=List(MaxSimps,s->[s]);

bool:=true;

while bool do
bool:=false;

for x in [1..Length(MaxSimps)] do
if IsBound(MaxSimps[x]) then

m:=MaxSimps[x];
if Order(m[Length(m)])>p then
for t in MaximalSubgroups(m[Length(m)]) do
mm:=Concatenation(m,[t]);
Add(MaxSimps,mm); if Order(t)>p then bool:=true; fi;
Unbind(MaxSimps[x]);
od;
fi;

fi;
od;
od;
MaxSimps:=Filtered(MaxSimps,i->IsBound(i));

K:= MaximalSimplicesToSimplicialComplex(MaxSimps);

for k in [1..Dimension(K)+1] do
Apply(K!.simplicesLst[k],x->SortedList(x,fn));
od;

return K;

end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(QuillenComplex,
function(G,p);
return PSubgroupSimplicialComplex(G,p,IsElementaryAbelian);
end);
####################################################
####################################################

######################################################
######################################################
InstallGlobalFunction(GChainComplex,
function(K,G)
local Ksimps,R, orbits, stabilizers, stabfn, Dim,  boundfn,
elts, inv, gg, i,j, k, x, y, m,Action ,ontuples, A, B;

elts:=Elements(G);
inv:=List(elts,x->Position(elts,x^-1));
Ksimps:=[];

for k in [1..1+Dimension(K)] do
Ksimps[k]:=List(K!.simplicesLst[k],x->SSortedList(x));
od;

#############################
Action:=function(a,b,c) return 1; end;
#############################

#############################
ontuples:=function(x,g)
local g1;
g1:=g;
return SSortedList(OnTuples(x,g1));
end;
#############################


orbits:=[];
for k in [1..1+Dimension(K)] do
orbits[k]:=OrbitsDomain(G,Ksimps[k], ontuples);
od;

stabilizers:=[];
for k in [1..1+Dimension(K)] do
stabilizers[k]:=[];
for i in [1..Length(orbits[k])] do
stabilizers[k][i]:=Stabilizer(G,orbits[k][i][1],ontuples);
od;od;

######################
Dim:=function(k);
if k<0 or k>Dimension(K) then return 0; fi;
return Length(orbits[k+1]);
end;
######################

######################
stabfn:=function(k,i);
return stabilizers[k+1][i];
end;
######################

######################
boundfn:=function(n,i)
local V,Vhat, ii, j, bnd,g,ob;

if n<=0 then return []; fi;

V:=orbits[n+1][i][1];

bnd:=[];

for j in [1..Length(V)] do
Vhat:=List(V,v->v);
RemoveSet(Vhat,V[j]);
Vhat:=SSortedList(Vhat);
ob:=fail;
for ii in [1..Length(orbits[n])] do
if Vhat in orbits[n][ii] then ob:=ii; break; fi;
od;
gg:=fail;
for g in [1..Length(elts)] do
if ontuples(orbits[n][ob][1],elts[g])=Vhat then 
gg:=g; break; fi; od;
if IsOddInt(j) then
Add(bnd,[ob,inv[gg]]);
else
Add(bnd,[-ob,inv[gg]]);
fi;
od;

return bnd;

end;
######################

R:=Objectify(HapGChainComplex,
            #HapNonFreeResolution,
            rec(
            dimension:=Dim,
            boundary:=boundfn,
            homotopy:=fail,
            elts:=elts,
            group:=G,
            stabilizer:=stabfn,
            action:=Action,
            properties:=
            [["length",1000],
             ["characteristic",0],
             ["type","chaincomplex"]]));

return R;
end);
######################################################
######################################################

######################################################
######################################################
InstallGlobalFunction(PSubgroupGChainComplex,
function(arg)
local G,p, filt,P, OS,  K, act, dm, stab, STAB, orbs, Dim, R, n,k, tmp,tmpp;

G:=arg[1];
p:=arg[2];
if Length(arg)=3 then
filt:=arg[3];
else
filt:=function(G); return true; end;
fi;

P:=SylowSubgroup(G,p);
K:=PSubgroupSimplicialComplex(P,p,filt);
dm:=Dimension(K);
K:=K!.simplicesLst;

orbs:=[];

#################################
act:=function(x,g);
return List(x,y->y^g);
end;
#################################

for n in [0..dm] do
orbs[n+1]:=[];
while Length(K[n+1])>0 do
Add(orbs[n+1],K[n+1][1]);
tmp:=Orbit(G,K[n+1][1],act);
K[n+1]:=Filtered(K[n+1],x-> not x in tmp);
od;
od;

#########################
Dim:=function(n);
if n<0 or n>dm then return 0; fi;
return Length(orbs[n+1]);
end;
#########################

STAB:=[];
for n in [0..dm] do
STAB[n+1]:=[];
for k in [1..Dim(n)] do
OS:=OrbitStabilizer(G,orbs[n+1][k],act);
STAB[n+1][k]:=OS.stabilizer;
#STAB[n+1][k]:=Stabilizer(G,orbs[n+1][k],act);
od;
od;

#########################
stab:=function(n,k);
return STAB[n+1][k];
end;
#########################

R:=Objectify(HapGChainComplex,
            rec(
            dimension:=Dim,
            boundary:=fail,
            homotopy:=fail,
            elts:=fail,
            group:=G,
            stabilizer:=stab,
            action:=fail,
            properties:=
            [["length",1000],
             ["characteristic",0],
             ["type","chaincomplex"]]));

return R;

end);
######################################################
######################################################

######################################################
######################################################
InstallGlobalFunction(HomologicalGroupDecomposition,
function(G,p)
local C, D, dm, n, k, j, iso;

C:=PSubgroupGChainComplex(G,p,IsElementaryAbelian);

dm:=0;n:=1;
while C!.dimension(n)>0 do
dm:=dm+1; n:=n+1;
od;

D:=[[],[]];
for n in [0..dm] do
if IsOddInt(n) then
for k in [1..C!.dimension(n)] do
Add(D[2],C!.stabilizer(n,k));
od;
else
for k in [1..C!.dimension(n)] do
Add(D[1],C!.stabilizer(n,k));
od;
fi;
od;

for k in [1..Length(D[1])] do
for j in [1..Length(D[2])] do
iso:=IsomorphismGroups(D[1][k],D[2][j]);
if not iso=fail then D[1][k]:=0; D[2][j]:=0; D[2]:=Filtered(D[2],a->not a=0);
break; fi;
od;
od;
D[1]:=Filtered(D[1],a->not a=0);
return D;
end);
######################################################
######################################################

