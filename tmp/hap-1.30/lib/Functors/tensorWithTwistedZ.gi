#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
#####################################################################
InstallGlobalFunction(TensorWithTwistedIntegers,
function(X,rho)
local 	TensorWithZ_Obj,
	TensorWithZ_Arr;

	
if EvaluateProperty(X,"characteristic")>0 then
Print("ERROR: You should use the function TensorWithIntegersModP(). \n");
return fail; fi;

#####################################################################
#####################################################################
TensorWithZ_Obj:=function(R,rho)
local  
	BoundaryC,
	LengthC,
	M,Bool;

Bool:=false;
if "tensorWithTwistedIntRec" in NamesOfComponents(R) then
if R!.tensorWithTwistedIntRec[2]=rho then Bool:=true;fi;
fi;
if Bool then return R!.tensorWithTwistedIntRec[1]; fi;


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
returnvec[i]:=returnvec[i]+SignInt(bound[x][1])*rho(R!.elts[bound[x][2]]);
od;

return returnvec;
end;
#####################################################################

R!.tensorWithTwistedIntRec:=[];
R!.tensorWithTwistedIntRec[1]:= Objectify(HapChainComplex,
		rec(
		dimension:=R!.dimension,
		boundary:=BoundaryC,
	        twist:=rho,
		properties:=
		[["length",LengthC],
		["connected",true],
		["type", "chainComplex"],
		["characteristic", 
		EvaluateProperty(R,"characteristic")] ]));
R!.tensorWithTwistedIntRec[2]:=rho;
return R!.tensorWithTwistedIntRec[1];
end;
#####################################################################
#####################################################################





#####################################################################
#####################################################################
TensorWithZ_Arr:=function(F,rho)
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
C:=TensorWithZ_Obj(R,rho);
D:=TensorWithZ_Obj(S,rho);
DimensionC:=C!.dimension;

#####################################################################
RhomS:=function(w,n)
local x,y,v;
v:=[];

for x in w do
y:=F!.mapping([[x[2],1]],n);
Apply(y,t->[x[1]*SignInt(t[1]),AbsInt(t[1]),t[2]]);
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

v:=[];
for i in [1..DimensionS(n)] do
v[i]:=0;
od;

for x in w do
v[x[2]]:=v[x[2]]+x[1]*rho(S!.elts[x[3]]);


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
	   twist:=rho,
	   properties:=[ ["type","chainMap"],
	   ["characteristic", Maximum( 
	   EvaluateProperty(X!.source,"characteristic"),
	   EvaluateProperty(X!.target,"characteristic"))]
	   ]));
end;
#####################################################################
#####################################################################
if IsHapEquivariantChainComplex(X) then
return TensorWithZ_Obj(X,rho); fi;  #Added 12/02/2019

if EvaluateProperty(X,"type") = "resolution" then
return TensorWithZ_Obj(X,rho); fi;

if EvaluateProperty(X,"type") = "equivariantChainMap" then
return TensorWithZ_Arr(X,rho); fi;

Print("ERROR: Input should be a resolution or equivariant map between resolutions. \n");

end);
#####################################################################
#####################################################################
#####################################################################



