
###################################################################
###################################################################
InstallGlobalFunction(HAP_GenericSL2OSubgroup,
function()
local type, G;

    type := NewType( FamilyObj([[[1,0],[0,1]]]),
                     IsGroup and
                     IsAttributeStoringRep and
                     IsFinitelyGeneratedGroup and
                     IsMatrixGroup and
                     IsHapSL2OSubgroup);

G:=rec(
    membership:= fail,
    tree:= fail,
    generators:= fail,
    level:= fail,
    cosetRep:= fail,
    cosetPos:= fail,
    ugrp:= fail,
    name:="Congruence subgroup");

ObjectifyWithAttributes(G, type,
DimensionOfMatrixGroup, 2,
IsIntegerMatrixGroup, true,
IsFinite, false,
IsHapSL2OSubgroup, true,
IsHapSL2Subgroup, true);

return G;
end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallGlobalFunction(HAP_CongruenceSubgroupGamma0Ideal,
function(I)
local G, membership, g,x,y,a;

G:=HAP_GenericSL2OSubgroup();

###################################################
membership:=function(g);
#if not Determinant(g)=1 then return false; fi;
if not g[2][1] in I   then return false; fi;
return true; 
end;
###################################################

G!.membership:=membership;
G!.level:=I;
if Norm(I)=1 then
  G!.GeneratorsOfMagmaWithInverses:=GeneratorsOfGroup(SL2QuadraticIntegers(AssociatedRing(I)!.bianchiInteger,true));
fi;
G!.name:="CongruenceSubgroupGamma0";
G!.ugrp:=Group(IdentityMat(2));
return G;

end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallGlobalFunction(HAP_PrincipalCongruenceSubgroupIdeal,
function(I)
local G, membership, g,x,y,a;

G:=HAP_GenericSL2OSubgroup();

###################################################
membership:=function(g);
#if not Determinant(g)=1 then return false; fi;
if not g[2][1] in I   then return false; fi;
if not g[1][2] in I   then return false; fi;
if not (g[1][1] -1) in I   then return false; fi;
if not (g[2][2] -1) in I   then return false; fi;
return true;
end;
###################################################

G!.membership:=membership;
G!.level:=I;
G!.name:="PrincipalCongruenceSubgroup";
G!.ugrp:=Group(IdentityMat(2));
return G;

end);
###################################################################
###################################################################


###################################################################
###################################################################
InstallGlobalFunction(HAP_SL2OSubgroupTree_slow,
function(G)
local d, gens, one, tree,InGmodU,Ugrp,v,p,g,s,n,q,vv,P,R, sublst,
      leaves,nodes,generators,Perturb, InNodesModG, csts, elements,
      generators2,pos,i,j,PermRep,elementsNumbers;;
if G!.tree=fail then
d:=AssociatedRing(G!.level)!.bianchiInteger;
gens:=Generators(ResolutionSL2QuadraticIntegers(d,2),true);
gens:=SSortedList(gens);


R:=ResolutionSL2QuadraticIntegers(d,2,true);;
P:=PresentationOfResolution(R);
gens:= R!.elts{P.gens};

one:=IdentityMat(2);

Ugrp:=G!.ugrp;
Ugrp:=Elements(Ugrp);

tree:=[];
leaves:=NewDictionary(one,true,SL2QuadraticIntegers(d));
elementsNumbers:=NewDictionary(one,true,SL2QuadraticIntegers(d));

if Size(Ugrp)>1 then
###########################################
InGmodU:=function(g)
local x;
for x in Ugrp do
if G!.membership(x*g)  
then return true; fi;;
od;
return false;
end;
###########################################
else
     InGmodU:=G!.membership;
fi;

###########################################
InNodesModG:=function(g)
local x,y,gg,B1,B2;
gg:=g^-1;

for x in nodes do
y:=x*gg;
if InGmodU(y) 
then return LookupDictionary(elementsNumbers,x); fi;
od;

return false;
end;
###########################################


if Size(Ugrp)>0 then
###########################################
Perturb:=function(g)
local u;
for u in Ugrp do
if G!.membership(u*g) then return u*g; fi;
od;
return fail;
end;
###########################################
else
    Perturb:=function(g); return g; end;
fi;

nodes:=[one]; 
elements:=[one];
AddDictionary(elementsNumbers,one,1);
PermRep:=List(gens,i->[]);

for s in Reversed([1..Length(gens)]) do
 if not gens[s] in G then
  AddDictionary(leaves,gens[s],2);
  AddDictionary(elementsNumbers,gens[s],2);
 tree[2]:=[1,s];
 AddSet(nodes,gens[s]);
 Add(elements,gens[s]);
 PermRep[s][1]:=2;
 break;
 fi;
od;

generators:=[];


############Tree Construction########################
while Size(leaves)>0 do
vv:=leaves!.entries[1];
v:=vv[1];
    for s in [1..Length(gens)] do
        g:=v*gens[s]; 
        q:=InNodesModG(g);
        p:=LookupDictionary(leaves,v);
        if q=false then 
            AddDictionary(leaves,g,1+Size(tree));
            AddDictionary(elementsNumbers,g,1+Size(tree));
            AddSet(nodes,g);
            Add(tree,[p, s]);
            PermRep[s][p]:=Size(tree);
            Add(elements,g);
        else 
            Add(generators,[g,elements[q],v,s]);
            PermRep[s][p]:=q;
        fi;
    od;
RemoveDictionary(leaves,v);
od;
#####################################################

G!.tree:=tree;
G!.generators:=generators;
generators:=List(generators,x->Perturb(x[2]*x[1]^-1));
generators:=List(generators,x->Minimum(x,x^-1));
generators2:=SSortedList(generators);
generators2:=Filtered(generators2,x->not x=one);

G!.GeneratorsOfMagmaWithInverses:=generators2;
sublst:=Filtered([1..Length(generators)],i->generators[i] in generators2);
G!.generators:=G!.generators{sublst};
G!.elements:=elements;
G!.IndexInSL2O:=Length(tree);
G!.gens:=gens;
G!.presentation:=P;

for s in [1..Length(gens)] do
if not IsBound(PermRep[s][1]) then
PermRep[s][1]:=Filtered([1..Length(PermRep[s])],i->not i in PermRep[s])[1];
fi;
od;
G!.PermRep:=PermRep;
fi;
end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallGlobalFunction(HAP_TransversalCongruenceSubgroupsIdeal,
function(G,H)
local R, RR, lenR, x, poscan, gensH, I, t, one, P,Psorted,zero;

I:=H!.level;
one:=One(I!.AssociatedRing mod I);

if not (H!.name="CongruenceSubgroupGamma0" and IsPrime(I)) then
   HAP_SL2SubgroupTree(H);
   R:=H!.elements;
   lenR:=Length(R);
   ##########################################
   poscan:=function(x)
   local i, xinv;
   xinv:=x^-1;
   for i in [1..lenR] do
   if R[i]*xinv in H then return i; fi;
   od;
   end;
   ##########################################
else
   RR:=[]; P:=[];
   for t in RightTransversal(I) do
      Add(RR,[[1,0],[t,1]]);
   od;
   Add(RR,[[0,1],[-1,0]]);
   lenR:=Length(RR);
   P:=List(RR, x -> x[2]*one);
   P:=List(P, x->[x[1]![1],x[2]![1]]);
   Psorted:=SSortedList(P);
   R:=List(Psorted, p->RR[Position(P,p)]);
   zero:=Position(R,RR[lenR]); 
   ##########################################
   poscan:=function(x)
   local p;
   p:=x[2]*one; 
   if IsZero(p[2]) then return zero; fi;
   p:=p*(p[2]^-1); return Position(Psorted,[p[1]![1],p[2]![1]]);
   end;
   ##########################################
fi;

return Objectify( NewType( FamilyObj( G ),
                    IsHapRightTransversalSL2ZSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
          rec( group := G,
               subgroup := H,
               cosets:=R,
               poscan:=poscan ));

end);
###################################################################
###################################################################



###################################################################
###################################################################
InstallGlobalFunction(HAP_TransversalCongruenceSubgroupsIdeal_alt,
function(G,H)
local gensG, gensH, N, L, GN, HN, R, R2, x, S,  one, iso, epi,epi2, poscan;

gensH:=GeneratorsOfGroup(H);
for x in gensH do
if not x in G then
Print("The second argument is not a subgroup of the first.\n");
return fail;fi;
od;
gensG:=GeneratorsOfGroup(G);
N:=H!.level;
S:= AssociatedRing(N) mod N;

one:=One(S);
GN:=Group(gensG*one);
iso:=IsomorphismPermGroup(GN);
epi:=GroupHomomorphismByImagesNC(G,Image(iso),gensG,List(gensG*one,y->Image(iso,y)));
epi2:=GroupHomomorphismByFunction(G,Image(iso),y->Image(iso,y*one)  );

HN:=Group(List(gensH*one,y->Image(iso,y)));
R:=RightTransversal(Image(iso,GN),HN);
R2:=List(R,x->PreImagesRepresentative(epi,x));
L:=Length(R2);

##########################################
poscan:=function(x);
return PositionCanonical(R,ImagesRepresentative(epi2,x));
end;
##########################################

##########################################
#poscan:=function(x)
#local i,xx;
#xx:=x^-1;
#for i in [1..L]  do
#if R2[i]*xx in H then return i; fi;
#od;
#end;
##########################################

return Objectify( NewType( FamilyObj( G ),
                    IsHapRightTransversalSL2ZSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
          rec( group := G,
               subgroup := H,
               cosets:=R2,
               poscan:=poscan ));

end);
###################################################################
###################################################################



###################################################################
###################################################################
InstallGlobalFunction(HAP_TransversalGamma0SubgroupsIdeal,
function(G,H)
local  N, L, R, T, F, poscan, t, one, c;

N:=H!.level;
T:=RightTransversal(N);
R:=[];
for t in T do
Add(R, [[1 mod N,0 mod N],[t mod N,1 mod N]]);
od;

Add(R, [[0 mod N,1 mod N],[-1 mod N,0 mod N]]);

R:=SSortedList(R);

L:=Length(R);
F:=AssociatedRing(N) mod N;;
one:=One(F);

##########################################
poscan:=function(x)
local i,xx, a, b, c, d, A, B, C, D, pos;
xx:=x^-1;
#for i in [1..L]  do
#if R[i]*xx in H then return i; fi;
#od;
if not IsZero(x[2][2]*one) then 
   c:=(x[2][1]*one)*(x[2][2]*one)^-1; c:=c![1];
   #return Position(R,[[1 mod N,0],[c,1 mod N]]);  
   return Position(T,c);
fi; 
   #return Position(R,[[0 mod N,1 mod N],[-1 mod N,0 mod N]]); 
   return L;
end;
##########################################

return Objectify( NewType( FamilyObj( G ),
                    IsHapRightTransversalSL2ZSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
          rec( group := G,
               subgroup := H,
               cosets:=R,
               poscan:=poscan ));

end);
###################################################################
###################################################################

#####################################################
#####################################################
InstallOtherMethod(AbelianInvariants,
"for HAPSL2OSubgroups",
[IsHapSL2OSubgroup],
1000000, #Hmm!
function(H)
local P,G, CosetTable,r;

HAP_SL2OSubgroupTree_slow(H);
P:=H!.presentation;
G:=P.freeGroup/P.relators;

CosetTable:=[];
for r in H!.PermRep do
Add(CosetTable,r);
Add(CosetTable, ListPerm(PermList(r)^-1)  );
od;

return AbelianInvariantsSubgroupFpGroupRrs(G,CosetTable);

end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallMethod(AsFpGroup,
"for HAPSL2OSubgroups",
[IsHapSL2OSubgroup],
1000000, #Hmm!
function(H)
local I,d,R,P,G,gensG, gensH, HG, tree, loops,
vertex2word,x,iso,iso1,iso2, CosetTable,r;

HAP_SL2OSubgroupTree_slow(H);
P:=H!.presentation;
G:=P.freeGroup/P.relators;

gensG:=GeneratorsOfGroup(G);
HAP_SL2OSubgroupTree_slow(H);
tree:=H!.tree;
loops:=List(H!.generators,x->[Position(H!.elements,x[2]),Position(H!.elements,x[3]),x[4]]);

###################################
vertex2word:=function(v)
local word, x;
word:=One(G);
while IsBound(tree[v]) do
word:=gensG[tree[v][2]]*word;
v:=tree[v][1];
od;
return word;
end;
####################################

gensH:=[];
for x in loops do
Add(gensH, vertex2word(x[2])*gensG[x[3]]*(vertex2word(x[1]))^-1);
od;

#HG:=PresentationSubgroupRrs(G,CosetTable);
HG:=PresentationSubgroup(G,Group(gensH));
#HG!.TzOptions.printLevel:=0;
#TzEliminate(HG,Length(HG!.generators)-500);
HG:= FpGroupPresentation(HG);
return HG;
end);
#####################################################
#####################################################

