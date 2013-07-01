#############
InstallGlobalFunction(ConjugateSL2ZGroup,
function(H,P)
    local m,gens,i,j, G,gensG,NameH,p;
p:=P[2][2];
if H=SL(2,Integers) then return SL2Z(p);fi;
NameH:=Name(H);
i:=Position(NameH,'/');
j:=Position(NameH,']');
m:=Int(NameH{[i+1..j-1]});

gens:=GeneratorsOfGroup(SL2Z(1/m));
gensG:=List(gens,x->P*x*(P^-1));
G:=Group(gensG);
SetName(G,Concatenation("SL(2,Z[",String(1/m),"])^",String(P))  );
G!.coprimes:=[m,p];
SetIsHAPRationalMatrixGroup(G,true);
SetIsHAPRationalSpecialLinearGroup(G,true);
   
return G;
end);
####
InstallGlobalFunction(CongruenceSubgroup,
function(m,p)
local H,K,G;
if m=1 then return CongruenceSubgroupGamma0(p);fi;
H:=SL2Z(1/m);
K:=ConjugateSL2ZGroup(H,[[1,0],[0,p]]);
G:=Intersection(H,K);
SetName(G,Concatenation("CongruenceSubgroup of ",Name(H)," level ",String(p)));
G!.levels:=[m,p];
SetIsHAPRationalMatrixGroup(G,true);
SetIsHAPRationalSpecialLinearGroup(G,true);

return G;
end);
##########################################


