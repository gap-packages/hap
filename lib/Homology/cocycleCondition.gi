#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(CocycleCondition,
function(R,n)
local
	M, i, v, BoundaryToVector;
M:=[];

#####################################################################
BoundaryToVector:=function(m) 	#Take R!.boundary(n+1,m)
local x, v, w;

w:=List([1..R!.dimension(n)],x->0);
v:= ShallowCopy(R!.boundary(n+1,m));
v:=List(v,x->x[1]);

for x in v do
w[AbsoluteValue(x)]:=w[AbsoluteValue(x)] + SignInt(x);
od;

return w;
end;
#####################################################################

for i in [1..R!.dimension(n+1)] do
v:=BoundaryToVector(i);
if Sum(List(v,x->AbsoluteValue(x)))>0 then
Append(M,[v]); fi;
od;

return M;
end);
#####################################################################
