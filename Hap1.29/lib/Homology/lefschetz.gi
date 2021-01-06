
##########################################################
InstallGlobalFunction(LefschetzNumberOfChainMap,
function(F)
local
	C,map,dim,lefnum,zro,v,n,i;

if not IsHapChainMap(F) then
Print("The function can only be applied to a chain map.\n");
return fail;
fi;

if not Source(F)=Target(F) then
Print("The chain map must have a common source and target.\n");
return fail;
fi;

C:=Source(F);
map:=F!.mapping;
dim:=EvaluateProperty(C,"length");
lefnum:=0;

for n in [0..dim] do
zro:=[1..C!.dimension(n)]*0;
for i in [1..C!.dimension(n)] do
v:=StructuralCopy(zro); v[i]:=1;
lefnum:=lefnum + (-1)^n*map(v,n)[i];
od;
od;

return lefnum;

end);
##########################################################

##########################################################
InstallMethod(LefschetzNumber,"Lefschetz numbers of chain map",
[IsHapChain],0,
function(F)
return LefschetzNumberOfChainMap(F);
end);
##########################################################



