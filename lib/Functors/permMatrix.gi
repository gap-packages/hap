#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(PermToMatrixGroup,
function(G,n)
local
	gensG, gensA, A, GhomA, M, g, row, i,j;

gensG:=GeneratorsOfGroup(G);
gensA:=[];

for g in gensG do
M:=[];
for i in [1..n] do
row :=List([1..n],x->0);
row[i^g]:=1;
Append(M,[row]);
od;
Append(gensA,[M]);
od;

A:=Group(gensA);
GhomA:=GroupHomomorphismByImagesNC(G,A,gensG,gensA);

return GhomA;


end);
#####################################################################
