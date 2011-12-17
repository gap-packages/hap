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
C:=TensorWithZ_Obj(R);
D:=TensorWithZ_Obj(S);
DimensionC:=C!.dimension;

#####################################################################
RhomS:=function(w,n)
local x,y,v;
v:=[];

for x in w do
y:=F!.mapping([[x[2],1]],n);
Apply(y,t->[x[1]*SignInt(t[1]),AbsInt(t[1])]);
Append(v,y);
od;

return v;
end;
#####################################################################


#####################################################################
CmapR:=function(v,n)
local  i,j,w,x;

w:=[];

for i in [1..DimensionC(n)] do
if not v[i]=0 then
	Add(w,[v[i],i]);
fi;
od;

return w;
end;
#####################################################################

#####################################################################
SmapD:=function(w,n)
local i,x,v;

v:=List([1..DimensionS(n)],x->0);

for x in w do
v[x[2]]:=v[x[2]]+x[1];
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



