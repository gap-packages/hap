#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(VectorStabilizer,
function(G,v)
local	S,g;


if IsPermGroup(G) then
return Stabilizer(G,v,Permuted);
fi;

S:=[];
if IsMatrixGroup(G) then
for g in G do
if g*v = v then Append(S,[g]); fi;
od;
return Group(S);
fi;

Print("G must be a matrix group or a permutation group.\n");
return fail;

end);
#####################################################################
