#Read("Desktop/PREVIOUS/TopBook/DATA/simulatedData.txt");;
#DataMatrix:=Concatenation(DataMatrix,DataMatrix+[1,0,-1]);
#DataMatrix:=Concatenation(DataMatrix,DataMatrix+[0,-2,1]);
#S:=VectorsToSymmetricMatrix(DataMatrix);;
#G:=SymmetricMatrixToGraph(S,10);;

#########################################################
#########################################################
InstallGlobalFunction(RipsChainComplex,
function(G,N)
local  MaximalTree, Maximal2Complex, ChooseCriticalEdge,
       VERTICES, V, EDGES, EDGESET, vv,v,i,j,k,b,
       CRITICALVERTICES,CRITICALEDGES,
       EDGEBOUNDARIES, FACEBOUNDARIES, F, C,
       PseudoBoundary, Dimension, Boundary;

#Critical vertices and edges coloured 1 or -1. All initially critical.
#Redundant vertices and edges coloured 2

N:=N+1;

#This is not yet implemented for N>2. For the moment we use the following
#in the case N>2.

#######################
if N>2 then
F:=IncidenceMatrixToGraph(G!.incidenceMatrix);
ContractGraph(F);
C:=SimplicialNerveOfGraph(F,N);
C:=SimplicialComplexToRegularCWComplex(C);
C:=SparseChainComplex(C);
return C;
fi;
#######################


EDGES:=StructuralCopy(G!.incidenceMatrix);
V:=Length(EDGES);
VERTICES:=List([1..V],i->1);
EDGESET:=List([1..V],i->Filtered([1..V],j->EDGES[i][j]=1));

################################################
MaximalTree:=function(v)
local i,j,colour,leaves,newleaves;
#Find a maximal tree containing the free vertex v.
#Colour the edges of the tree 2, and also colour the "initial" vertices  2.

colour:=2;
VERTICES[v]:=-1;
leaves:=[v];

while Size(leaves)>0 do
newleaves:=[];
for i in leaves do
#for j in [1..V] do
for j in EDGESET[i] do
if EDGES[i][j]=1 and  VERTICES[j]=1 then 
EDGES[i][j]:=colour; EDGES[j][i]:=colour; VERTICES[j]:=colour;
Add(newleaves,j);
fi;
od;
od;
leaves:=newleaves;
od;

end;
################################################

################################################
Maximal2Complex:=function()
local i,j,k,toggle;

toggle:=true;

while toggle do
toggle:=false;
for i in [1..V] do
#for j in [i+1..V] do
for j in EDGESET[i] do

if EDGES[i][j]=1 then
for k in [1..V] do
if EDGES[i][k]+EDGES[j][k]=4 then
EDGES[i][j]:=2; EDGES[j][i]:=2; toggle:=true; break; fi;
od;
fi;

od;od;
EDGESET:=List([1..V],i->Filtered([1..V],j->EDGES[i][j]=1));
od;

end;
################################################

################################################
ChooseCriticalEdge:=function()
local i,j,k;

for i in [1..V] do
#for j in [i+1..V] do
for j in EDGESET[i] do
if EDGES[i][j]=1 then
for k in [1..V] do
if EDGES[i][k]=1 and not EDGES[j][k]=0 then
EDGES[i][k]:=-1;EDGES[k][i]:=-1;EDGES[i][j]:=2; return true; fi;

if EDGES[j][k]=1 and not EDGES[i][k]=0 then
EDGES[j][k]:=-1;EDGES[k][j]:=-1;EDGES[i][j]:=2; return true; fi;
od;

fi;
od;od;

return false;

end;
################################################

###################
Dimension:=function(n);
if n=0 then return Length(CRITICALVERTICES); fi;
if n<=N and n>0 then
return Length(PseudoBoundary[n]); fi;
return 0;
end;
###################

###################
Boundary:=function(n,i);
if n=0 then return []; fi;
return PseudoBoundary[n][i];
end;
###################

EDGEBOUNDARIES:=[];
FACEBOUNDARIES:=[];
PseudoBoundary:=[[],[]];;
 
###################
C:=Objectify(HapSparseChainComplex,
                rec(
                dimension:=Dimension,
                boundary:=Boundary,
                properties:=
                [["length",N],
                ["type", "chainComplex"],
                ["characteristic",0 ] ]));
###################

###################
CRITICALVERTICES:=[1..V];
if N=0 then return C; fi;
###################

###################
v:=Position(VERTICES,1);
while IsInt(v) do
MaximalTree(v);
v:=Position(VERTICES,1);
od;
CRITICALVERTICES:=Filtered([1..V],v->AbsInt(VERTICES[v])=1);
if N=1 then return C; fi;
###################

Maximal2Complex();
while ChooseCriticalEdge() do
Maximal2Complex();
od;

EDGES:=List(EDGES,row->List(row,x->AbsInt(x)));
CRITICALEDGES:=[];;
for i in [1..V] do
for j in [i+1..V] do
if EDGES[i][j]=1 then Add(CRITICALEDGES,[i,j]);fi;
od;od;

EDGEBOUNDARIES:=List(CRITICALEDGES,e->Intersection(e,CRITICALEDGES));
EDGEBOUNDARIES:=List(EDGEBOUNDARIES,e->List(e,i->Position(CRITICALVERTICES,i)));
FACEBOUNDARIES:=[];
F:=Length(CRITICALEDGES);

for i in [1..F] do
for j in [i+1..F] do
if Length(Intersection(CRITICALEDGES[i],CRITICALEDGES[j]))>0 then 
for k in [j+1..F] do
if Length(Intersection(CRITICALEDGES[k],CRITICALEDGES[i]))>0
and Length(Intersection(CRITICALEDGES[k],CRITICALEDGES[j]))>0 then
Add(FACEBOUNDARIES,[i,j,k]);
fi;
od;
fi;
od;od;

v:=List([1..Length(CRITICALVERTICES)],i->0);

for b in EDGEBOUNDARIES do
vv:=[];;
for i in b do
Add(vv,[AbsInt(i),SignInt(i)]);
od;
Add(PseudoBoundary[1],vv);
od;

v:=List([1..Length(CRITICALEDGES)],i->0);

for b in FACEBOUNDARIES do
vv:=[];;
for i in b do
Add(vv, [AbsInt(i), SignInt(i)]);
od;
Add(PseudoBoundary[2],vv);
od;

if N=2 then return C; fi;

end);
#########################################################
#########################################################


