

###################################################################
###################################################################
InstallGlobalFunction(HAP_chain_bockstein,
function(C,n,prime)
local G,A,i,x,a,b,c,h,s,bhomc,B,CM,psi,R,dim,
      D;

if not EvaluateProperty(C,"characteristic")=0 then
Print("The chain complex must be of characteristic 0.\n");
return fail;
fi;

if n=0 then 
D:=HomToIntegersModP(C,prime);
A:=Cohomology(D,0);
B:=Cohomology(D,1);
A:=AbelianGroup(List([1..A],i->prime));
B:=AbelianGroup(List([1..B],i->prime));
return GroupHomomorphismByFunction(A,B,x->One(B));
fi;

G:=CyclicGroup(1);
b:=CyclicGroup(prime^2);;
x:=GeneratorsOfGroup(b)[1];
c:=Group(x^prime);;

B:=TrivialGModuleAsGOuterGroup(G,b);
CM:=TrivialGModuleAsGOuterGroup(G,c);
bhomc:=GroupHomomorphismByImages(b,c,[x],[x^prime]);
psi:=GOuterGroupHomomorphism();
psi!.Source:=B;
psi!.Target:=CM;
psi!.Mapping:=bhomc;

dim:=Length(C);
R:=HAP_ChainComplexToEquivariantChainComplex(C);
R!.properties[1]:=["length",dim+1];

return ConnectingCohomologyHomomorphism(psi,n,R)!.Mapping;

end);
###################################################################
###################################################################





