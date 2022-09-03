
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
local tree,InGmodU,Ugrp,S,T,U,v,p,g,s,n,q,vv,gens,
      nodes, leaves, ambientGenerators, InLowDegreeNodesModG, 
      one, genTriples, vertex2word, triple2word, csts;

####################
S:=[[0,-1],[1,0]];;
T:=[[1,1],[0,1]];
U:=S*T;
one:=IdentityMat(2);
####################

ambientGenerators:=[S,S*U];
Ugrp:=G!.ugrp;
Ugrp:=Elements(Ugrp);
tree:=[1 ];
genTriples:=[];
cnt:=1;
leaves:=NewDictionary(S,true,SL(2,Integers));
nodes:=[one];
AddDictionary(leaves,one,1);

if Size(Ugrp)>1 then
###########################################
InGmodU:=function(g)
local u;
for u in Ugrp do
if G!.membership(u*g) then return true; fi;;
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




############Tree Construction########################
while Size(leaves)>0 do
vv:=leaves!.entries[1];   
v:=vv[1];
    for s in [1..Length(ambientGenerators)] do
        g:=v*ambientGenerators[s]; 
        q:=InLowDegreeNodesModG(g);
        if q=false then 
         Add(nodes,g);
         AddDictionary(leaves,g,Length(nodes));
            p:=LookupDictionary(leaves,v);
            Add(tree,[p, s]);
            else Add(genTriples,[v,g,q]);
        fi;
    od;
RemoveDictionary(leaves,v);
od;
#####################################################

#####################################################
vertex2word:=function(v)
local word;
word:=one;
while not v=1 do
word:=ambientGenerators[tree[v][2]]*word;
v:=tree[v][1];
od;
return word;
end;
#####################################################

#####################################################
triple2word:=function(x)
local u,uu,g,q,c,ans;
#Print(Position(nodes,x[3]),"  ",x[2],"\n");
ans:=[];
for u in Reversed(Ugrp) do
c:=x[2]*u*vertex2word(  Position(nodes,x[3])  )^-1;
if c in G then Add(ans,c); fi;
od;

return ans;
return fail;  #This should never happen
end;
#####################################################

G!.tree:=tree;
Unbind(tree[1]);
genTriples:=List(genTriples,x->triple2word(x));
genTriples:=Concatenation(genTriples);
genTriples:=List(genTriples,x->Minimum(x,x^-1));
genTriples:=SSortedList(genTriples);
G!.GeneratorsOfMagmaWithInverses:=genTriples;

end);
###################################################################
###################################################################


###################################################################
###################################################################
InstallGlobalFunction(HAP_SL2ZSubgroupTree_fast,
function(G)
local ambientGenerators, tree,InGmodU,Ugrp,S,T,U,v,p,g,s,n,q,vv,
      leaves,genTriples,generators, InLowDegreeNodesModG, csts, cnt,
      vertex2word, one, triple2word,i,j,u,c,a,b;

####################
S:=[[0,-1],[1,0]];;
T:=[[1,1],[0,1]];
U:=S*T;
one:=IdentityMat(2);
####################

ambientGenerators:=[S*U,S*U^2];
Ugrp:=G!.ugrp;
Ugrp:=Elements(Ugrp);

tree:=[ ];
genTriples:=[];
cnt:=1;
leaves:=NewDictionary(S,true,SL(2,Integers));
csts:=[];
csts[G!.cosetPos(one)]:=1;
AddDictionary(leaves,one,G!.cosetPos(one));

###########################################
InLowDegreeNodesModG:=function(g)
local pos;
pos:=G!.cosetPos(g);
if not IsBound(csts[pos]) then return false; else;
return pos; fi;
end;
###########################################


#########The work horse##############################
while Size(leaves)>0 do
vv:=1*leaves!.entries[1];
v:=vv[1];
p:=G!.cosetPos(v);
    for s in [1..Length(ambientGenerators)] do
        g:=ambientGenerators[s]*v;
        q:=InLowDegreeNodesModG(g);
        if q=false then
             q:=G!.cosetPos(g);
             AddDictionary(leaves,g,q);
             csts[q]:=1;
             tree[q]:=[p,s];
            else Add(genTriples,[p,s,v,g]);  #Usp=Uq mod G
        fi;
    od;
RemoveDictionary(leaves,v);
od;
#####################################################

#####################################################
vertex2word:=function(v)
local word;
word:=one;
while IsBound(tree[v]) do
word:=word*ambientGenerators[tree[v][2]];
v:=tree[v][1];
od;
return word;
end;
#####################################################

#####################################################
triple2word:=function(x)
local u,uu,g,q,c;
for u in Ugrp do
#Print(vertex2word(G!.cosetPos(x[3])) = x[3], " ");
c:=x[4]^-1*u*vertex2word(G!.cosetPos(x[4]));
if c in G then return c; fi;
od;
return fail;  #This should never happen
end;
#####################################################


genTriples:=List(genTriples,x->triple2word(x));
genTriples:=List(genTriples,x->Minimum(x,x^-1));
genTriples:=SSortedList(genTriples);
G!.tree:=tree;
G!.GeneratorsOfMagmaWithInverses:=genTriples;

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
G!.name:="PrincipalCongruenceSubgroup";
return G;
fi;
##############################################
one:=One(sln);
U:=S*T*one;
Ugrp:=Group(U); #########################################
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
G!.ugrp:=Group((S*T));
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

if false then
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
fi;

###################################################################
###################################################################
InstallGlobalFunction(HAP_TransversalCongruenceSubgroups, 
function(G,H)
local tree,InH,S,T,U,v,p,g,s,n,q,vv,gens,
      nodes, nodesinv, leaves, ambientGenerators, InLowDegreeNodesModH,
      one, poscan, nind ;

####################
S:=[[0,-1],[1,0]];;
T:=[[1,1],[0,1]];
U:=S*T;
one:=IdentityMat(2);
####################

ambientGenerators:=[S,S*U];
tree:=[1 ];
cnt:=1;
leaves:=NewDictionary(S,true,SL(2,Integers));
nodes:=[one];
AddDictionary(leaves,one,1);

InH:=H!.membership;

###########################################
InLowDegreeNodesModH:=function(g)
local x,gg,B1,B2;
gg:=g^-1;

for x in nodes do
if InH(x*gg) then return x; fi;
od;

return false;
end;
###########################################




############Tree Construction########################
while Size(leaves)>0 do
vv:=leaves!.entries[1];
v:=vv[1];
    for s in [1..Length(ambientGenerators)] do
        g:=v*ambientGenerators[s];
        q:=InLowDegreeNodesModH(g);
        if q=false then
         Add(nodes,g);
         AddDictionary(leaves,g,Length(nodes));
            p:=LookupDictionary(leaves,v);
            Add(tree,[p, s]);
        fi;
    od;
RemoveDictionary(leaves,v);
od;
#####################################################

nodes:=Filtered(nodes,g-> g in G);
nodesinv:=List(nodes,g->g^-1);
nind:=[1..Length(nodes)];

####################################################
poscan:=function(x)
local i;

for i in nind do
if InH(  x*nodesinv[i]  ) then return i; fi;
od;
return fail;
end;
####################################################

return Objectify( NewType( FamilyObj( G ),
                    IsHapRightTransversalSL2ZSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
          rec( group := G,
               subgroup := H,
               cosets:=nodes,
               poscan:=poscan ));
end);
###################################################################
###################################################################


############################################################
############################################################
InstallOtherMethod(RightTransversal,
"right transversal for finite index subgroups of SL(2,Integers)",
[IsHapSL2ZSubgroup,IsHapSL2ZSubgroup],
0,  #There must be a better way to ensure this method is not used!
#[IsHapSL2ZSubgroup,IsMatrixGroup],
#1000000,
function(H,HH)
Print("Using   HAP_TransversalCongruenceSubgroups\n");
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

if H!.tree=true then
return HAP_TransversalCongruenceSubgroupsIdeal(H,HH);
else
return HAP_TransversalCongruenceSubgroupsIdeal_alt(H,HH);
fi;

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


