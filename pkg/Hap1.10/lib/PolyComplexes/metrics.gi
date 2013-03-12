##########################################
InstallGlobalFunction(KendallMetric,
function(g,h,N)
local S, i,j;
S:=0;

for i in [1..N] do
for j in [i+1..N] do
if (i^g-j^g)*(i^h-j^h)<0 then S:=S+1; fi;
od;od;

return S;
end);
##########################################

##########################################
InstallGlobalFunction(CayleyMetric,
function(g,h,N)
local S,T,R,i,a;

a:=h*g^-1;
S:=Sum(Flat(CycleStructurePerm(a)));
T:=0;;
for i in [1..N] do
if i=i^a then T:=T+1; fi;
od;

return
N-S-T;

end);
##########################################

##########################################
InstallGlobalFunction(HammingMetric,
function(g,h,N)
local T,R,i,a;

a:=h*g^-1;
T:=0;;
for i in [1..N] do
if i=i^a then T:=T+1; fi;
od;

return
N-T;

end);
##########################################


