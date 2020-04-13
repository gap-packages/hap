
###################################################################
###################################################################
InstallGlobalFunction(HAP_GenericSL2Subgroup,
function()
local type, G;

    type := NewType( FamilyObj([[[1,0],[0,1]]]),
                     IsGroup and
                     IsAttributeStoringRep and
                     IsFinitelyGeneratedGroup and
                     IsMatrixGroup and
                     IsHapSL2Subgroup);

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
"Generating set for Hap_SL2Subgroups",
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
   if G!.cosetRep=fail then 
      HAP_SL2SubgroupTree_slow(G);
   else
      HAP_SL2SubgroupTree_fast(G);
   fi;
fi;
end);
###################################################################
###################################################################


###################################################################
###################################################################
InstallGlobalFunction(HAP_SL2SubgroupTree_slow,
function(G)
local tree,InGmodU,Ugrp,S,T,U,U1,U2,v,p,g,s,n,q,vv,gens,
      leaves,generators,Perturb, InLeavesAndGensModG, csts;

#This function is mainly of theoretical value and not nearly as 
#fast as the fast version. Furthermore, at present it works only
#for the principal congruence subgroups.

S:=[[0,-1],[1,0]];;
T:=[[1,1],[0,1]];
U:=S*T;
U1:=S*U;
U2:=S*U^2;

Ugrp:=G!.ugrp;

Ugrp:=Elements(Ugrp);

#if U^2 in G or U^3 in G then return fail; fi;
#For the moment we only handle groups that act freely on the vertices of
#the cubic tree.

tree:=[ , [1,0],[1,1],[1,2]];
#Despite what is written in all the LEFTIST comments, 
#tree is a set of RIGHT cosets reps.
#tree is a list with i-th entry equal to an integer pair [j,s] denoting
#that there is a directed edge in the tree from vertex i to vertex j and 
#where the edge is G[i]=U^s*S*G[j] ---> G[j]. Here G[j] is the i-th 
#element of the group G. We let G[1] be the identity 
#matrix. Vertex 1 is the root of the tree and so tree[1] is unbound.
#Initially just vertices 2,3 and 4 point to a vertex (namely to vertex 1).

leaves:=NewDictionary(S,true,SL(2,Integers));
#leaves is a list of matrices corresponding to just those tree vertices that 
#currently have degree = 1 or 2. Looking up a matrix in the dictionary 
#yields an integer n where n is the tree vertex corresponding to the matrix. 

###########################################
InGmodU:=function(g)
local x;
for x in Ugrp do
if G!.membership(x*g) then return true; fi;;
od;
return false;
end;
###########################################

###########################################
InLeavesAndGensModG:=function(g)
local x,gg,B1,B2;
gg:=g^-1;

for x in leaves!.entries do
if InGmodU(x[1]*gg) then return x[1]; fi;
od;
for x in gens do
if InGmodU(x[1]*gg) then RemoveSet(gens,x); return x[1]; fi;
if InGmodU(x[3]*gg) then RemoveSet(gens,x); return x[3]; fi;
od;


return false;
end;
###########################################


###########################################
Perturb:=function(g)
local u;
for u in Ugrp do
if G!.membership(u*g) then return u*g; fi;
od;
return fail;
end;
###########################################


AddDictionary(leaves,S,2);
AddDictionary(leaves,U1,3);
AddDictionary(leaves,U2,4);

generators:=[];
gens:=[];

############Tree Construction########################
while Size(leaves)>0 do
vv:=leaves!.entries[1];   
v:=vv[1];
    for s in [1,2] do
        if s=1 then g:=v*U1; else g:=v*U2; fi;
        q:=InLeavesAndGensModG(g);
        if q=false then AddDictionary(leaves,g,1+Size(tree)); 
            p:=LookupDictionary(leaves,v);
            Add(tree,[p, s]);
            else Add(generators,[v,g,q]);
                 AddSet(gens,[v,g,q]);
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
InstallGlobalFunction(HAP_SL2SubgroupTree_fast,
function(G)
local tree,InGmodU,Ugrp,S,T,U,U1,U2,v,p,g,s,n,q,vv,gens,Lift,
      leaves,generators,Perturb, InLeavesAndGensModG, csts;

S:=[[0,-1],[1,0]];;
T:=[[1,1],[0,1]];
U:=S*T;
U1:=S*U;
U2:=S*U^2;
Ugrp:=G!.ugrp;
Ugrp:=Elements(Ugrp);

#if U^2 in G or U^3 in G then return fail; fi;
#For the moment we only handle groups that act freely on the vertices of
#the cubic tree.

tree:=[ ];
#Despite what is written in the comments, tree is a set of RIGHT cosets reps.
#tree is a list with i-th entry equal to an integer pair [j,s] denoting
#that there is a directed edge in the tree from vertex i to vertex j and
#where the edge is G[i]=U^s*S*G[j] ---> G[j]. Here G[j] is the i-th
#element of the group G. We let G[1] be the identity
#matrix. Vertex 1 is the root of the tree and so tree[1] is unbound.
#Initially just vertices 2,3 and 4 point to a vertex (namely to vertex 1).

leaves:=NewDictionary(S,true,SL(2,Integers));
#leaves is a list of matrices corresponding to just those tree vertices that
#currently have degree = 1 or 2. Looking up a matrix in the dictionary
#yields an integer n where n is the tree vertex corresponding to the matrix.

###########################################
csts:=[];
InLeavesAndGensModG:=function(g);
if not IsBound(csts[G!.cosetPos(g)]) then return false; fi;
return G!.cosetRep(g);
end;
###########################################


###########################################
Perturb:=function(g)
local u;
for u in Ugrp do
if G!.membership(u*g) then return u*g; fi;
od;
return [fail,g];
end;
###########################################


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
        q:=InLeavesAndGensModG(g);
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

sl:=SL(2,Integers);

G:=HAP_GenericSL2Subgroup();

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
sln:=SL(2,Integers mod n);
S:=[[0,-1],[1,0]];;
T:=[[1,1],[0,1]];
one:=One(sln);
U:=S*T*one;
Ugrp:=Group(U);
Uelts:=Elements(Ugrp);
#R:=RightTransversal(sln,Ugrp);
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
local G,sl, sln, S, T, U, Ugrp, R,RR, UT,membership, CosetRep, CosetPos,
      Uelts, one,g,x,y,a;

sl:=SL(2,Integers);

G:=HAP_GenericSL2Subgroup();

###################################################
membership:=function(g);
if not g in sl then return false; fi;
if not g[2][1] mod n = 0  then return false; 
else return true; fi;
end;
###################################################

G!.membership:=membership;
G!.level:=n;

sln:=SL(2,Integers mod n);
S:=[[0,-1],[1,0]];;
T:=[[1,1],[0,1]];
one:=One(sln);
U:=S*T*one;
#UT:=Filtered(Elements(sln),IsUpperTriangularMat);
UT:=[];
for a in [1..n] do
if Gcd(a,n)=1 then Add(UT,[[a,0],[0,a^-1]]*one); fi;
od;
Add(UT,[[1,1],[0,1]]*one);
UT:=Group(UT);
Uelts:=Elements(UT);
Ugrp:=Elements(Group(U^0));

RR:=RightTransversal(sln,UT);
#RR:=Enumerator(sln);;
R:=[];;
for g in RR do
x:=Minimum(List(Uelts,u->u*g));
Add(R,SSortedList(Ugrp*x));
od;
R:=SSortedList(R);

###################################################
CosetPos:=function(g)
local gg,y,a;
gg:=g*one;
gg:=Minimum(List(Uelts,u->u*gg));

a:=Position(R,SSortedList(Ugrp*gg));
if not a=fail then return a; fi;

return fail;
end;
###################################################
###################################################
CosetRep:=function(g)
local gg;
gg:=g*one;
gg:=Minimum(List(Uelts,u->u*gg));
gg:=Uelts*gg;
return gg[1];

end;
###################################################

G!.cosetPos:=CosetPos;
G!.cosetRep:=CosetRep;
G!.ugrp:=Group((S*T)^0);
G!.name:="CongruenceSubgroupGamma0";
return G;

end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallGlobalFunction(HAP_TransversalCongruenceSubgroups,
function(G,H)
local gensG, gensH, N, GN, HN, R, x, sln, one, epi;

gensH:=GeneratorsOfGroup(H);
for x in gensH do
if not x in G then
Print("The second argument is not a subgroup of the first.\n");
return fail;fi;
od;
gensG:=GeneratorsOfGroup(G);
#GG:=Group(gensG);;
N:=H!.level;

sln:=SL(2,Integers mod N);
one:=One(sln);
GN:=Group(gensG*one);
epi:=GroupHomomorphismByImagesNC(G,GN,gensG,gensG*one);
HN:=Group(gensH*one);
R:=RightTransversal(GN,HN);
R:=List(R,x->PreImagesRepresentative(epi,x));

return R;
end);
###################################################################
###################################################################

############################################################
############################################################
InstallOtherMethod(RightTransversal,
"right transversal for finite index subgroups of SL(2,Integers)",
[IsHapSL2Subgroup,IsHapSL2Subgroup],
function(H,HH)

return HAP_TransversalCongruenceSubgroups(H,HH);

end);
############################################################
############################################################

