
#####################################################################
InstallGlobalFunction(ModPCohomologyRing,
function(arg)
local A;

if Length(arg)=1 then
A:=ModPCohomologyRing_part_1(arg[1]);fi;
if Length(arg)=2 then
A:=ModPCohomologyRing_part_1(arg[1],arg[2]);fi;
if Length(arg)=3 then
A:=ModPCohomologyRing_part_1(arg[1],arg[2],arg[3]);fi;

if arg[Length(arg)]="high" or arg[Length(arg)]="low" then
A:=ModPCohomologyRing_part_2(A); fi;

return A;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(ModPCohomologyGenerators,
function(arg)
local A;

if IsHapResolution(arg[1]) then
A:=ModPCohomologyRing_part_1(arg[1],"low");fi;
if IsGroup(arg[1]) then
A:=ModPCohomologyRing_part_1(arg[1],arg[2],"low");fi;

return [ModPRingGenerators(A),A!.degree];
end);
#####################################################################

