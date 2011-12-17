#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
#####################################################################
InstallGlobalFunction(TensorWithIntegers,
function(X)
local 	TensorWithZ_Obj,
	TensorWithZ_Arr;

	
if EvaluateProperty(X,"characteristic")>0 then
Print("ERROR: You should use the function TensorWithIntegersModP(). \n");
return fail; fi;

#####################################################################
#####################################################################
TensorWithZ_Obj:=function(R)
local  
	BoundaryC,
	LengthC,
	M;

if "tensorWithIntRec" in NamesOfComponents(R) then
return R!.tensorWithIntRec; fi;


LengthC:=EvaluateProperty(R,"length");
M:=[1..LengthC];				

#####################################################################
#BoundaryC:=function(n,k)
#local
#	row, Mt, i, j, x, sum;

#if n <0 then return false; fi;
#if n=0 then return [0]; fi;

#if M[n]=n then 
#   Mt:=[];

#   for i in [1..R!.dimension(n-1)] do
#   row:=[];
#        for j in [1..R!.dimension(n)] do
#        sum:=0;
#                for x in R!.boundary(n,j) do
#                if AbsoluteValue(x[1])=i then
#                sum := sum + SignInt(x[1]);
#                fi;
#                od;
#        row[j]:=sum;
#        od;
#   Mt[i]:=row;
#   od;

#   M[n]:=TransposedMat(Mt);
#fi;

#return M[n][k];
#end;
#####################################################################

#####################################################################
BoundaryC:=function(n,k)
local returnvec, bound, x, i;

if n <0 then return false; fi;
if n=0 then return [0]; fi;
returnvec:=0*[1..R!.dimension(n-1)];
# 0*[1..n] is faster than List([1..n],i->0)
# in (seemingly) any case.
# For large n, NullMat(1,n)[1] is faster than 0*[1..n].

bound:=R!.boundary(n,k);
for x in [1..Size(bound)]
do
i:=AbsInt(bound[x][1]);
returnvec[i]:=returnvec[i]+SignInt(bound[x][1]);
od;

return returnvec;
end;
#####################################################################

R!.tensorWithIntRec:= Objectify(HapChainComplex,
		rec(
		dimension:=R!.dimension,
		boundary:=BoundaryC,
		properties:=
		[["length",LengthC],
		["connected",true],
		["type", "chainComplex"],
		["characteristic", 
		EvaluateProperty(R,"characteristic")] ]));
return R!.tensorWithIntRec;
end;
#####################################################################
#####################################################################





#####################################################################
#####################################################################
TensorWithZ_Arr:=function(F)
local
		R,S,RhomS,        	#R->S is an equivariant chain
		CmapR, SmapD, 		#map. C->D is the chain map
		C,D,ChomD,		#got by killing the actions.
		DimensionS,
		DimensionC,
		x;
R:=F!.source;
S:=F!.target;
DimensionS:=S!.dimension;
RhomS:=F!.mapping;
C:=TensorWithZ_Obj(R);
D:=TensorWithZ_Obj(S);
DimensionC:=C!.dimension;

#####################################################################
CmapR:=function(v,n)
local  i,j,w,x;

w:=[];

for i in [1..DimensionC(n)] do
if not v[i]=0 then
	x:=[SignInt(v[i])*i,1];
	for j in [1..AbsoluteValue(v[i])] do
	     Append(w,[x]);
	od;
fi;
od;

return w;
end;
#####################################################################

#####################################################################
SmapD:=function(w,n)
local i,x,v;

v:=[];
for i in [1..DimensionS(n)] do
v[i]:=0;
od;

for x in w do
v[AbsoluteValue(x[1])]:=v[AbsoluteValue(x[1])]+SignInt(x[1]);
od;

return v;
end;
#####################################################################

#####################################################################
ChomD:=function(v,n);
return
SmapD(RhomS(CmapR(v,n),n),n);
end;
#####################################################################


return Objectify(HapChainMap,
	rec(
	   source:=C,
	   target:=D,
	   mapping:=ChomD,
	   properties:=[ ["type","chainMap"],
	   ["characteristic", Maximum( 
	   EvaluateProperty(X!.source,"characteristic"),
	   EvaluateProperty(X!.target,"characteristic"))]
	   ]));
end;
#####################################################################
#####################################################################

if EvaluateProperty(X,"type") = "resolution" then
return TensorWithZ_Obj(X); fi;

if EvaluateProperty(X,"type") = "equivariantChainMap" then
return TensorWithZ_Arr(X); fi;

Print("ERROR: Input should be a resolution or equivariant map between resolutions. \n");

end);
#####################################################################
#####################################################################
#####################################################################



