
###################################################################
###################################################################
InstallGlobalFunction(HAP_GenericSL2ZSubgroup,
function()
local type, G;

    type := NewType( FamilyObj([[[1,0],[0,1]]]),
                     IsGroup and
                     IsAttributeStoringRep and
                     IsFinitelyGeneratedGroup and
                     IsMatrixGroup and
                     IsHapSL2ZSubgroup);

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
IsHapSL2ZSubgroup, true,
IsHapSL2Subgroup, true);

return G;
end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallMethod(\in,
"membership test for Hap_SL2Subgroups",
[IsMatrix, IsHapSL2Subgroup and IsGroup],
1000000,  #There must be a better way to ensure this method is used!
function(g,G)
return  G!.membership(g);
end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallMethod(GeneratorsOfGroup,
"Generating set for Hap_SL2ZSubgroups",
[IsHapSL2Subgroup and IsGroup],
1000000,  #There must be a better way to ensure this method is used!
function(G)
HAP_SL2SubgroupTree(G);
return  G!.GeneratorsOfMagmaWithInverses;
end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallGlobalFunction(HAP_SL2SubgroupTree,
function(G);
    
if G!.tree=fail then
 if IsRing(G!.level) then
      HAP_SL2OSubgroupTree_slow(G);
 else  
   if G!.cosetRep=fail then 
      HAP_SL2ZSubgroupTree_slow(G);
   else
      HAP_SL2ZSubgroupTree_fast(G);
   fi;
 fi;
fi;
end);
###################################################################
###################################################################


###################################################################
###################################################################
InstallGlobalFunction(HAP_SL2ZSubgroupTree_slow,
function(G)
local tree,InGmodU,Ugrp,S,T,U,U1,U2,v,p,g,s,n,q,vv,gens,
      leaves,nodes,generators,Perturb, InLowDegreeNodesModG, csts;

S:=[[0,-1],[1,0]];;
T:=[[1,1],[0,1]];
U:=S*T;
U1:=S*U; 
U2:=S*U^2;
Ugrp:=G!.ugrp;
Ugrp:=Elements(Ugrp);

tree:=[];
leaves:=NewDictionary(S,true,SL(2,Integers));

if Size(Ugrp)>1 then
###########################################
InGmodU:=function(g)
local x;
for x in Ugrp do
if G!.membership(x*g) then return true; fi;;
od;
return false;
end;
###########################################
else
     InGmodU:=G!.membership;
fi;

###########################################
InLowDegreeNodesModG:=function(g)
local x,gg,B1,B2;
gg:=g^-1;

for x in nodes do
if InGmodU(x*gg) then return x; fi;
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

v:=1;;
nodes:=[];
if not S in G then v:=v+1; AddDictionary(leaves,S,v); tree[v]:=[1,0]; 
Add(nodes,S);fi;
if not U1 in G then v:=v+1; AddDictionary(leaves,U1,v); tree[v]:=[1,1]; Add(nodes,U1); fi;
if not U2 in G then v:=v+1; AddDictionary(leaves,U2,v); tree[v]:=[1,2];
Add(nodes,U2);fi;
Add(nodes,S^0);
nodes:=SSortedList(nodes);

generators:=[];
gens:=[];

############Tree Construction########################
while Size(leaves)>0 do
vv:=leaves!.entries[1];   
v:=vv[1];
    for s in [1,2] do
        if s=1 then g:=v*U1; else g:=v*U2; fi;
        q:=InLowDegreeNodesModG(g);
        if q=false then AddDictionary(leaves,g,1+Size(tree)); 
         AddSet(nodes,g);
            p:=LookupDictionary(leaves,v);
            Add(tree,[p, s]);
            else Add(generators,[v,g,q]);
        fi;
    od;
RemoveDictionary(leaves,v);
od;
#####################################################

G!.tree:=tree; 
generators:=List(generators,x->Perturb(x[3]*x[2]^-1));
generators:=List(generators,x->Minimum(x,x^-1));
generators:=SSortedList(generators);
generators:=Filtered(generators,x->not x=IdentityMat(2));
G!.GeneratorsOfMagmaWithInverses:=generators;
end);
###################################################################
###################################################################


###################################################################
###################################################################
InstallGlobalFunction(HAP_SL2ZSubgroupTree_fast,
function(G)
local tree,InGmodU,Ugrp,S,T,U,U1,U2,v,p,g,s,n,q,vv,gens,Lift,
      leaves,generators,Perturb, InLowDegreeNodesModG, csts;

S:=[[0,-1],[1,0]];;
T:=[[1,1],[0,1]];
U:=S*T;
U1:=S*U;
U2:=S*U^2;
Ugrp:=G!.ugrp;
Ugrp:=Elements(Ugrp);


tree:=[ ];

leaves:=NewDictionary(S,true,SL(2,Integers));

###########################################
csts:=[];
InLowDegreeNodesModG:=function(g);
if not IsBound(csts[G!.cosetPos(g)]) then return false; fi;
return G!.cosetRep(g);
end;
###########################################

if Size(Ugrp)>1 then
###########################################
Perturb:=function(g)
local u;
for u in Ugrp do
if G!.membership(u*g) then return u*g; fi;
od;
return [fail,g];
end;
###########################################
else 
    Perturb:=function(g); return g; end;
fi;


AddDictionary(leaves,S,2);
AddDictionary(leaves,U1,3);
AddDictionary(leaves,U2,4);

generators:=[];
gens:=[];

#####################################################
tree[G!.cosetPos(S)]:=[G!.cosetPos(S^0),0];
tree[G!.cosetPos(U1)]:=[G!.cosetPos(S^0),1];
tree[G!.cosetPos(U2)]:=[G!.cosetPos(S^0),2];

csts[G!.cosetPos(S^0)]:=1;
csts[G!.cosetPos(S)]:=1;
csts[G!.cosetPos(U1)]:=1;
csts[G!.cosetPos(U2)]:=1;
#####################################################
while Size(leaves)>0 do
vv:=leaves!.entries[1];
v:=vv[1];
    for s in [1,2] do
        if s=1 then g:=v*U1; 
        else  g:=v*U2; fi;
        q:=InLowDegreeNodesModG(g);
        if q=false then
             AddDictionary(leaves,g,1+Size(tree));
             #p:=LookupDictionary(leaves,v);
             csts[G!.cosetPos(g)]:=1;
             p:=G!.cosetPos(v);
            #Add(tree,[p, s]);
            tree[G!.cosetPos(g)]:=[p,s];
            else Add(generators,[v,g]);
                 AddSet(gens,[v,g]);
        fi;
    od;
RemoveDictionary(leaves,v);
od;
#####################################################

#####################################################
#####################################################
Lift:=function(g)
local r,root,x,ROOT;
if not IsBound(tree[G!.cosetPos(g)]) then return S^0; fi;
r:=S^0;
root:=G!.cosetPos(r);
x:=tree[G!.cosetPos(g)];

while not x[1]=root do
if x[2]=1 then r:=U1*r;fi;
if x[2]=2 then r:=U2*r; fi;
#if x[2]=0 then r:=S*r; fi;
x:=tree[x[1]];
od;
if x[2]=1 then r:=U1*r; fi;
if x[2]=2 then r:=U2*r; fi;
if x[2]=0 then r:=S*r;fi;
return r;
end;
#####################################################
#####################################################

G!.tree:=tree;
generators:=List(generators,x-> [x[1],x[2],Lift(x[2])]);

generators:=List(generators,x->Perturb(x[2]*x[3]^-1));
generators:=List(generators,x->Minimum(x,x^-1));
generators:=SSortedList(generators);
generators:=Filtered(generators,x->not x=IdentityMat(2));
G!.GeneratorsOfMagmaWithInverses:=generators;
end);
###################################################################
###################################################################


###################################################################
###################################################################
InstallGlobalFunction(HAP_PrincipalCongruenceSubgroup,
function(n)
local G,sl, sln, S, T, U, Ugrp, R,RR, membership, CosetRep, CosetPos, 
      Uelts, one,g;

if IsRing(n) then return
   HAP_PrincipalCongruenceSubgroupIdeal(n);
fi;

sl:=SL(2,Integers);
G:=HAP_GenericSL2ZSubgroup();

###################################################
membership:=function(g);
if not g in sl then return false; fi;
if not g[1][1] mod n = 1  then return false; fi;
if not g[2][2] mod n = 1  then return false; fi;
if not g[1][2] mod n = 0  then return false; fi;
if not g[2][1] mod n = 0  then return false; fi;
return true;
end;
###################################################

G!.membership:=membership;
G!.level:=n;
G!.name:="PrincipalCongruenceSubgroup";
G!.index:=n^3*Product(List(SSortedList(Factors(n)), p->1-1/p^2));
sln:=SL(2,Integers mod n);
S:=[[0,-1],[1,0]];;
T:=[[1,1],[0,1]];

##############################################
if n<=2 then
G!.ugrp:=Group((S*T)^0);
G!.name:="CongruenceSubgroupGamma0";
return G;
fi;
##############################################
one:=One(sln);
U:=S*T*one;
Ugrp:=Group(U);
Uelts:=Elements(Ugrp);
RR:=Enumerator(sln);;
R:=[];;
for g in RR do
Add(R,Minimum(List(Uelts,u->u*g)));
od;
R:=SSortedList(R);

###################################################
CosetPos:=function(g)
local gg;
gg:=g*one;
gg:=Minimum(List(Uelts,u->u*gg));
return Position(R,gg);
end;
###################################################
###################################################
CosetRep:=function(g)
local gg;
gg:=g*one;
gg:=Minimum(List(Uelts,u->u*gg));

return gg;

end;
###################################################

G!.cosetPos:=CosetPos;
G!.cosetRep:=CosetRep;
G!.ugrp:=Group(S*T);
return G;

end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallGlobalFunction(HAP_SL2TreeDisplay,
function(G)
local A, T, i, j,L;
HAP_SL2SubgroupTree(G);
T:=G!.tree;
A:=0*IdentityMat(Length(T));;
for i in [1..Length(T)] do
if IsBound(T[i]) then
A[i][T[i][1]]:=1;
A[T[i][1]][i]:=1;
fi;
od;
L:=Filtered([1..Length(A)],i-> not IsZero(A[i]));
A:=A{L};
A:=TransposedMat(A)*1;
A:=A{L};

A:=IncidenceMatrixToGraph(A);
GraphDisplay(A);
end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallGlobalFunction(HAP_CongruenceSubgroupGamma0,
function(n)
local G,sl,S,T,membership,g,x,y,a;

if IsRing(n) then
    return HAP_CongruenceSubgroupGamma0Ideal(n);
fi;

sl:=SL(2,Integers);
G:=HAP_GenericSL2ZSubgroup();

###################################################
membership:=function(g);
if not g in sl then return false; fi;
if not g[2][1] mod n = 0  then return false; 
else return true; fi;
end;
###################################################

G!.membership:=membership;
G!.level:=n;

S:=[[0,-1],[1,0]];;
T:=[[1,1],[0,1]];

G!.ugrp:=Group((S*T)^0);
G!.name:="CongruenceSubgroupGamma0";
G!.index:=n*Product(List(SSortedList(Factors(n)), p->1+1/p));
return G;

end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallGlobalFunction(HAP_TransversalCongruenceSubgroups,
function(G,H)
local gensG, gensH, N, GN, iso, HN, R, R2, x, sln, one, epi, epi2, poscan;

#AT PRESENT THIS APPROACH IS SLOW AND SO RANKED VERY LOW AS A METHOD
Print("HAP should not be using this method.\n");
gensH:=GeneratorsOfGroup(H);
for x in gensH do
if not x in G then
Print("The second argument is not a subgroup of the first.\n");
return fail;fi;
od;
gensG:=GeneratorsOfGroup(G);
N:=H!.level;

sln:=SL(2,Integers mod N);
one:=One(sln);
GN:=Group(gensG*one);
iso:=IsomorphismPermGroup(GN);
epi:=GroupHomomorphismByImagesNC(G,Image(iso),gensG,List(gensG*one,x->Image(iso,x)));
epi2:=GroupHomomorphismByFunction(G,GN,x->Image(iso,x*one));
HN:=Group(List(gensH*one,x->Image(iso,x)));

R:=RightTransversal(Image(iso),HN);
R2:=List(R,x->PreImagesRepresentative(epi,x));

##########################################
poscan:=function(x);
return PositionCanonical(R,ImagesRepresentative(epi2,x));
end;
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

############################################################
############################################################
InstallOtherMethod(RightTransversal,
"right transversal for finite index subgroups of SL(2,Integers)",
[IsMatrixGroup,IsHapSL2ZSubgroup],
0,  #There must be a better way to ensure this method is not used!
function(H,HH)
return HAP_TransversalCongruenceSubgroups(H,HH);

end);
############################################################
############################################################

############################################################
############################################################
InstallOtherMethod(RightTransversal,
"right transversal for finite index subgroups of SL2QuadraticIntegers(d))",
[IsHapSL2OSubgroup,IsHapSL2OSubgroup],
1000000, #Must be a better way than this to ensure this method
function(H,HH)
local N;

if IsPrime(HH!.level) and IsRing(HH!.level) then 
    N:=HH!.level;
    if AssociatedRing(N)!.bianchiInteger<0 then
        return HAP_TransversalGamma0SubgroupsIdeal(H,HH); 
    fi;
fi;

return HAP_TransversalCongruenceSubgroupsIdeal(H,HH);

end);
############################################################
############################################################


############################################################
############################################################
InstallMethod(IndexInSL2Z,
"Index of HAP_congruence subgroup in SL(2,Integers)",
[IsHapSL2ZSubgroup],
function(H)

return H!.index;

end);
############################################################
############################################################

############################################################
############################################################
InstallMethod(IndexInSL2O,
"Index of HAP_congruence subgroup in SL(2,Integers)",
[IsHapSL2OSubgroup],
function(H)

HAP_SL2OSubgroupTree_slow(H);
H!.index:=Length(H!.tree);

return H!.index;

end);
############################################################
############################################################


