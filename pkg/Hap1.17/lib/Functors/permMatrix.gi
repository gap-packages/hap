#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(PermToMatrixGroup,
function(arg)
local
	G, n, gensG, gensA, A, GhomA, M, g, row, i,j;

G:=arg[1];
if Length(arg)=2 then n:=arg[2];
else n:=Maximum(MovedPoints(G));
fi;
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
