#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(IntegralCupProduct,
function(arg)
local 
	R, u, v, p, q, P, Q, N,
	uCocycle,
	vCocycle,
	uvCocycle,
	uChainMap,
	DimensionR,
	i, w, x, sw;

	####################BEGIN TO READ THE INPUT##################
R:=arg[1];
DimensionR:=R!.dimension;
u:=arg[2];
v:=arg[3];
p:=arg[4];
q:=arg[5];

if Length(arg)>5 then P:=arg[6];
else
P:=CR_CocyclesAndCoboundaries(R,p,true); 
fi;

if Length(arg)>6 then Q:=arg[7];
else
Q:=CR_CocyclesAndCoboundaries(R,q,true);
fi;

if Length(arg)>7 then N:=arg[8];
else
N:=CR_CocyclesAndCoboundaries(R,p+q,true);
fi;
	#####################FINISHED REAQDING THE INPUT#############

uCocycle:=P.classToCocycle(u);
vCocycle:=Q.classToCocycle(v);
uChainMap:=CR_ChainMapFromCocycle(R,uCocycle,p,q);

uvCocycle:=[];
for i in [1..DimensionR(p+q)] do
w:=uChainMap([[i,1]]);
sw:=0;
	for x in w do
	sw:=sw+ SignInt(x[1])*vCocycle[AbsoluteValue(x[1])];
	od;
uvCocycle[i]:=sw;
od;

return N.cocycleToClass(uvCocycle);
end);
#####################################################################
#####################################################################
