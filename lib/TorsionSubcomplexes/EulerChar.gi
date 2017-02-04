InstallGlobalFunction(EquivariantEulerCharacteristic,
function(C)

local N,i,j,t;

N:=0;
while C!.dimension(N)>0 do
    N:=N+1;
od;
N:=N-1;
t:=0;
for i in [0..N] do
    for j in [1..C!.dimension(i)] do
        t:=t+(-1)^i*(1/Size(C!.stabilizer(i,j)));
    od;
od;

return t;
end);

