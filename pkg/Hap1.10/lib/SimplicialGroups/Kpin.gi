#C Graham Ellis

###############################################
###############################################
InstallGlobalFunction(EilenbergMacLaneSimplicialGroup_alt,
function(G,n,dim)
local
        KPi1, KPi2;

###############################################
KPi1:=function(G,dim)
#
#Inputs a group G and returns dim+1 terms of the simplicial
#group K(G,1) with G in each dimension.
#
local hom;

hom:=GroupHomomorphismByImages(G,G,GeneratorsOfGroup(G),GeneratorsOfGroup(G));


return Objectify(HapSimplicialGroup, rec(
           groupsList:=function(i) return G; end,
           boundariesList:=function(n,i) return hom; end,
           degeneraciesList:=function(n,i) return hom; end,           
           properties:=[["length",dim]]
          ));
end;
###############################################


###############################################
KPi2:=function(G,dim)
#
#Inputs an abelian group G and returns dim+1 terms of a simplicial
#abelian group K(G,2) [the first term being indexed by 0.
#
local m, n, i, j, L, GroupsList, BoundariesList, DegeneraciesList, Properties,
Embed, Project, gens, imgens, S, T, hom, a, x;

################################
if not IsAbelian(G) then 
Print("This function must be applied to an abelian group.\n");
return fail;
fi;
################################

GroupsList:=[Group(Identity(G))];
BoundariesList:=[];
DegeneraciesList:=[];
Properties:=[["length",dim]];

################################
for n in [1..dim] do
Add(GroupsList, DirectProduct(List([1..n],i->G))  );
od;
################################

################################
Embed:=function(n,i,g);
return Image(Embedding(GroupsList[n+1],i),g);
end;

Project:=function(n,i,g);
return Image(Projection(GroupsList[n+1],i),g);
end;
################################

################################
S:=GroupsList[2];
T:=GroupsList[1];
gens:=GeneratorsOfGroup(S);
imgens:=List([1..Length(gens)],a->Identity(T));
hom:=GroupHomomorphismByImages(S,T,gens,imgens);
L:=[hom,hom];
BoundariesList:=[L];

for n in [2..dim] do
L:=[];
S:=GroupsList[n+1];
T:=GroupsList[n];
gens:=GeneratorsOfGroup(S);

for i in [0..n] do
imgens:=[];

for x in gens do
a:=Identity(T);

for j in [1..i-1] do
a:=a*Embed(n-1,j,Project(n,j,x));
od;

if i>0 and i<n then
a:=a*Embed(n-1,i,Project(n,i,x))*Embed(n-1,i,Project(n,i+1,x));
fi;

for j in [i+2..n] do
a:=a*Embed(n-1,j-1,Project(n,j,x));
od;


Add(imgens,a);
od;


hom:=GroupHomomorphismByImages(S,T,gens,imgens);
Add(L,hom);
od;

Add(BoundariesList,L);
od;
################################


return Objectify(HapSimplicialGroup, rec(
           groupsList:=function(i) return GroupsList[i+1]; end,
           boundariesList:=function(n,i) return BoundariesList[n][i+1]; end,
           degeneraciesList:=function(n,i) return DegeneraciesList[n+1][i]; end,
           properties:=Properties
          ));
end;
###############################################

if not n in [1,2] then
Print("This function is currently only implemented for K(G,1) and K(G,2). \n");
return fail;
fi;

if n=1 then return KPi1(G,dim); fi;
if n=2 then return KPi2(G,dim); fi;
end);
###############################################
###############################################
