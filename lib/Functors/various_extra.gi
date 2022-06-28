if IsPackageMarkedForLoading("congruence","0.0") then

######################
InstallMethod( \in,
               "for CongruenceSubgroupGamma0(p)",
              [ IsMatrix,  IsCongruenceSubgroupGamma0 ],
function ( g, G )
local p;
p:=G!.LevelOfCongruenceSubgroup;
if g in SL(2,Integers) then
if Determinant(g)=1 then
if IsInt(g[2][1]/p) then return true;fi;
fi;
fi;
end);
############################

fi;
