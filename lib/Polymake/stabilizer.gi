#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(VectorStabilizer,
function(G,v)
local	S,toggle,g,i;

S:=[];

if IsPermGroup(G) then

return Stabilizer(G,v,Permuted);

#for g in G do
#toggle:=true;
#for i in [1..Length(v)] do
#if not v[i] = v[i^g] then toggle:=false; fi;
#od;
#if toggle then Append(S,[g]); fi;
#od;
#return Group(S);
fi;

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
