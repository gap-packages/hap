InstallGlobalFunction(EquivariantEulerCharacteristic,
function(arg)

local R,N,i,j,t,n;

if Length(arg)=1 then R:=arg[1];

N:=0;
while R!.dimension(N)>0 do
    N:=N+1;
od;
N:=N-1;
t:=0;
for i in [0..N] do
    for j in [1..R!.dimension(i)] do
        t:=t+(-1)^i*(1/Size(R!.stabilizer(i,j)));
    od;
od;

return t;

else 
R:=arg[1]; n:=arg[2];t:=0;
for i in [n..n] do
    for j in [1..R!.dimension(i)] do
        t:=t+(-1)^i*(1/Size(R!.stabilizer(i,j)));
    od;
od;
	return t;
fi;
end);

