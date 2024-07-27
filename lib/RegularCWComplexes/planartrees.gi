
#HAP_PlanarTreeJoin(S,T) inputs two rooted planar trees and outputs J= S V T.
#############################################
#############################################
InstallGlobalFunction(HAP_PlanarTreeJoin,
function(S,T)
local TT, SS, lS;

SS:=S+1;
SS:=Concatenation([0],SS);
lS:=Length(SS);
TT:=T+lS;
TT[1]:=1;
return Concatenation(SS,TT);

end);
#############################################
#############################################

#############################################
#############################################
InstallGlobalFunction(HAP_PlanarTreeLeaves,
function(T);
return Difference([1..Length(T)],T);
end);
#############################################
#############################################

#HAP_PlanarTreeGraft(T,S,l) grafts the root of tree S to leaf i of T
#############################################
#############################################
InstallGlobalFunction(HAP_PlanarTreeGraft,
function(T,S,l)
local L,r,s,start, middle, finish;
L:=HAP_PlanarTreeLeaves(T);
r:=L[l]*1;
s:=T[r]*1;
start:=T{[1..r-1]};
middle:=r-1+1*S;
middle[1]:=s;
finish:=Length(S)-1+1*T{[r+1..Length(T)]};
if Length(finish)>0 then
finish[1]:=s;
fi;
return Concatenation(start,middle,finish);
end);
#############################################
#############################################

#HAP_PlanarBinaryTrees(n) returns a list of the rooted planar binary trees with n
#leaves.
#############################################
#############################################
InstallGlobalFunction(HAP_PlanarBinaryTrees,
function(n)
local PBT, i,j,k,S,T;
PBT:=[];
PBT[1]:=[ [0] ];;
PBT[2]:=[   [0,1,1]   ];;
PBT[3]:=[   [0,1,2,2,1],  [0,1,1,3,3]    ];;
if n<=2 then return PBT[n]; fi;
for k in [3..n] do
PBT[k]:=[];
for i in [1..k-1] do
for S in PBT[i] do
for T in PBT[k-i] do
Add(PBT[k],HAP_PlanarTreeJoin(S,T));
od;od;od;od;
return PBT[n];
end);
#############################################
#############################################

#HAP_PlanarTreeRemoveEdge(T,i) removes the i-th edge of the tree T,
#thus reducing by 1 the number of interior vertices and leaving 
#unchanged the number of leaves.
#############################################
############################################# 
InstallGlobalFunction(HAP_PlanarTreeRemoveEdge,
function(T,i)
local A,B,fn;

if i<2 then return T; fi;
#############
fn:=function(k)
local x;
if k<i then x:= T[k];
else x:= T[k+1];
fi;
if x<i then return x; fi;
if x=i then return T[i]; fi;
return x-1;
end;
#############

return List([1..Length(T)-1],fn);

end);
#############################################
#############################################

#HAP_Associahedron(n) returns a list whose i-th term is the list
#of i-1-cells of the associahedron of dimension n-2.
#############################################
#############################################
InstallGlobalFunction(HAP_AssociahedronCells,
function(n)
local Cells, k, T,nodes, v;

Cells:=[];
Cells[1]:=HAP_PlanarBinaryTrees(n);

for k in [1..n-2] do
Cells[k+1]:=[];
for T in Cells[k] do
nodes:=SSortedList(T);
nodes:=Filtered(nodes,i->i>1);
for v in nodes do
Add(Cells[k+1], HAP_PlanarTreeRemoveEdge(T,v));
od;
od;
Cells[k+1]:=SSortedList(Cells[k+1]);
od;
return Cells;
end);
#############################################
#############################################

#############################################
#############################################
InstallGlobalFunction(HAP_AssociahedronBoundaries,
function(n)
local Cells, Boundaries,nodes,cb,k,T,i,v;
Cells:=HAP_AssociahedronCells(n);
Boundaries:=List(Cells, x-> List([1..Size(x)],i->[]));
Boundaries[1]:=List(Cells[1],x->[1,0]);

for k in [1..Length(Cells)-1] do

for i in [1..Length(Cells[k])] do
T:=Cells[k][i];
nodes:=SSortedList(T);
nodes:=Filtered(nodes,j->j>1);
for v in nodes do
cb:=HAP_PlanarTreeRemoveEdge(T,v);
Add( Boundaries[k+1][Position(Cells[k+1],cb)], i);
od;
od;
Apply(Boundaries[k+1],x->Concatenation([Length(x)],SortedList(x)));

od;

Add(Boundaries,[]);

return [Boundaries,Cells];
end);
#############################################
#############################################

#############################################
#############################################
InstallGlobalFunction(RegularCWAssociahedron,
function(arg)
local n,bounds,Y,k,B;
n:=arg[1];
bounds:=HAP_AssociahedronBoundaries(n);
Y:=RegularCWComplex(bounds[1]);
if Length(arg)=2 then return Y; fi;
Y!.trees:=bounds[2];
Y!.directed:=[];
for k in [1..Y!.nrCells(1)] do
B:=Y!.boundaries[2][k]{[2,3]};
if Y!.trees[1][B[1]] < Y!.trees[1][B[2]] then Add(Y!.directed,1*B);
else Add(Y!.directed,Reversed(1*B)); fi;
od;
CriticalCells(Y);;
return Y;
end);
#############################################
#############################################


############################################
############################################
InstallGlobalFunction(HAP_DisplayPlanarTree,
function(arg)
local Arg,T,i,n,tmpDir,treedot,treepng, cnt;

tmpDir:=DirectoryTemporary();
treedot:=Filename(tmpDir,"tree.dot");
treepng:=Filename(tmpDir,"tree.png");

Arg:=arg;
if Length(Arg)=1 then
if not IsInt(Arg[1][1]) then Arg:=Arg[1]; fi;
fi;
PrintTo(treedot,"digraph BST {\n" );
cnt:=0;
for n in [1..Length(Arg)] do

T:=Arg[n];
AppendTo(treedot,    "node [fontname=\"Arial\",fontsize=34];\n");
for i in [1..Length(T)] do
AppendTo(treedot, i+cnt," -> ", T[i]+cnt, " ;\n"   );
od;
cnt:=1+cnt+Length(T);
od;

AppendTo(treedot, "}\n"   );
Exec(Concatenation("dot -Tpng ", treedot, " > ", treepng ));
Exec(Concatenation(DISPLAY_PATH, " ", treepng));
RemoveFile(treedot);
RemoveFile(treepng);
Print(RemoveDir(Filename(tmpDir,"")));

end);
############################################
############################################

############################################
############################################
InstallGlobalFunction(HAP_PlanarTreeRemovableEdge,
function(T)
local IsRemovable, pos, k;

###################
IsRemovable:=function(k)
local pos;
if k=0 or k=1 then return false; fi;
pos:=Position(T,k);
if pos=fail then return false; fi;
if Position(T,T[k],k)=fail then return false; fi;
return true;
end;
###################

for k in SSortedList(T) do
if IsRemovable(k) then return k; fi;
od;

return fail;
end);
############################################
############################################


#############################################
#############################################
InstallGlobalFunction(HAP_DisplayVectorField,
function(Y,n)
local tmpDir, treedot,treepng, L,M,i,b, j;

L:=Y!.inverseVectorField[n+1];
M:=Y!.boundaries[n+2];
M:=List(M,x->x{[2..Length(x)]});

tmpDir:=DirectoryTemporary();
treedot:=Filename(tmpDir,"tree.dot");
treepng:=Filename(tmpDir,"tree.png");

PrintTo(treedot,"digraph BST {\n" );

AppendTo(treedot,    "node [fontname=\"Arial\"];\n");


for i in [1..Length(L)] do
if IsBound(L[i]) then
#AppendTo(treedot, i,".",Y!.height[n+1][i]," -> U", L[i], "\n");
AppendTo(treedot, i," -> U", L[i], "\n");
b:=Difference( M[L[i]], [i]);

for j in b do
#AppendTo(treedot, "U",L[i]," -> ", j,".",Y!.height[n+1][j], "\n");
AppendTo(treedot, "U",L[i]," -> ", j, "\n");
od;
fi;
od;

AppendTo(treedot, "}\n"   );
Exec(Concatenation("dot -Tpng ", treedot, " > ", treepng ));
Exec(Concatenation(DISPLAY_PATH, " ", treepng));
end);
#############################################
#############################################

#############################################
#############################################
RegularCWAssociahedronWithDiscreteVectorField:=function(n)
local Y, B, T,i,b, pos, V, k,U, Height, h, ii;
Y:=RegularCWComplex(HAP_AssociahedronBoundaries(n));
B:=HAP_AssociahedronCells(n);
V:=[];
V[1]:=[];


for i in [1..Length(B[1])] do
pos:=HAP_PlanarTreeRemovableEdge(B[1][i]);
if IsInt(pos) then
T:=HAP_PlanarTreeRemoveEdge(B[1][i],pos);
pos:=Position(B[2],T);
V[1][i]:=pos;
fi;
od;

for k in [2..n-2] do
V[k]:=[];
for i in Difference([1..Length(B[k])],SSortedList(V[k-1])) do
pos:=HAP_PlanarTreeRemovableEdge(B[k][i]);
T:=HAP_PlanarTreeRemoveEdge(B[k][i],pos);
pos:=Position(B[k+1],T);
V[k][i]:=pos;
od;
od;

U:=List(V,x->[]);
for k in [1..Length(U)] do
for i in [1..Length(V[k])] do
if IsBound(V[k][i]) then
U[k][V[k][i]]:=i;
fi;
od;
od;

Y!.inverseVectorField:=V;
Y!.vectorField:=U;
#Y!.criticalCells:=[[0,1]];

Height:=List([0..Dimension(Y)],x->[]);

Y!.height:=Height;
for i in [1..Y!.nrCells(0)] do
h:=0;
ii:=i;
while IsBound(Y!.inverseVectorField[1][ii]) do
b:=Y!.boundaries[2][Y!.inverseVectorField[1][ii]];
b:=b{[2,3]};
b:=Difference(b,[ii]);
ii:=b[1];
h:=h+1;
od;
Y!.height[1][i]:=h;
od;

for k in [1..Dimension(Y)] do
for i in [1..Y!.nrCells(k)] do
b:=Y!.boundaries[k+1][i];
b:=b{[2..b[1]]};
b:=List(b,j->Y!.height[k][j]);
b:=Maximum(b);
Y!.height[k+1][i]:=b;
od;
od;
return Y;
end;
#############################################
#############################################

################################################
################################################
InstallGlobalFunction(HAP_CriticalCellsDirected,
function(Y)
local arrows,k,e,a,aa,aaa,n,d,i,pos1,pos2,nonCriticals, ncells, n1cells,
vectorField, inverseVectorField, criticalCells, comb;
#Construct a contracting discrete vector field on
#a directed regular CW-complex and return the unique critical cell.

if not Y!.criticalCells=fail then return Y!.criticalCells; fi;
if not IsBound(Y!.directed) then return fail; fi;

arrows:=List([1..Y!.nrCells(0)],i->[]); #arrows[k] is a list of edge numbers
                                        #of edges with source vertex k.
for e in [1..Length(Y!.directed)] do
Add(arrows[Y!.directed[e][1]] , e);
od;
Apply(arrows,SSortedList);
criticalCells:=[];
vectorField:=List([1..Dimension(Y)],i->[]);
inverseVectorField:=List([0..Dimension(Y)-1],i->[]);
nonCriticals:=List([0..Dimension(Y)],i->[]);

#Add those arrows to the discrete vector field that
#have a source 0-cell and a target 1-cell.
for k in [1..Length(arrows)] do
a:=arrows[k];
if Length(a)>0 then
a:=a[1];  #Number of first edge with source vertex k.
d:=Y!.directed[a];  #The ordered pair of boundary vertices of edge a
inverseVectorField[1][d[1]]:=a;
vectorField[1][a]:=d[1];
Add(nonCriticals[1],d[1]);
Add(nonCriticals[2],a);
fi;
od;

#if false then
#Now add those arrows to the discrete vector field that
#have a source k-cell and a target k+1-cell for k>0.
#TO BE COMPLETED
for n in [1..Dimension(Y)-1] do
ncells:=[];
n1cells:=[];
for k in [1..Y!.nrCells(n)] do
a:=ClosureCWCell(Y,n,k)[2];;
a:=a[2];
Add(ncells,SSortedList(a));
od;
for k in [1..Y!.nrCells(n+1)] do
a:=ClosureCWCell(Y,n+1,k)[2];;
a:=a[2];
Add(n1cells,SSortedList(a));
od;
for k in [1..Y!.nrCells(0)] do
a:=arrows[k];
if Length(a)>=n+1 then
comb:=Combinations(a{[2..Length(a)]},n);
for aa in comb do
aaa:=Concatenation([a[1]],aa);;
pos1:=PositionProperty(ncells,x->IsSubset(x,aa));
pos2:=PositionProperty(n1cells,x->IsSubset(x,aaa));
inverseVectorField[n+1][pos1]:=pos2;
vectorField[n+1][pos2]:=pos1;
Add(nonCriticals[n+1],pos1);
Add(nonCriticals[n+2],pos2);
od;
fi;
od;od;
#END OF TO BE COMPLETED
#fi;

criticalCells:=List([0..Dimension(Y)],n->Difference([1..Y!.nrCells(n)],nonCriticals[n+1]));
criticalCells:=List([0..Dimension(Y)],n->List(criticalCells[n+1],x->[n,x]));
criticalCells:=Concatenation(criticalCells);
#pos1:=PositionProperty(arrows,x->Length(x)=0);
#criticalCells:=[[0,pos1]];

Y!.vectorField:=vectorField;
Y!.inverseVectorField:=inverseVectorField;
Y!.criticalCells:=criticalCells;

return criticalCells;;
end);
################################################
################################################

