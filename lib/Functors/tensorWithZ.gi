#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
#####################################################################
InstallGlobalFunction(TensorWithIntegers,
function(X)
local 	TensorWithZ_Obj,
        TensorNFWithZ_Obj,
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

if IsBound(R!.action) then
return TensorNFWithZ_Obj(R); fi;

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
TensorNFWithZ_Obj:=function(R)
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

bound:=R!.boundary(n,k);
for x in [1..Size(bound)]
do
i:=AbsInt(bound[x][1]);
returnvec[i]:=returnvec[i]+SignInt(bound[x][1])*R!.action(n-1,i,bound[x][2]);
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
ChomD:=function(v,n)
local ans;
 return SmapD(RhomS(CmapR(v,n),n),n);
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

if IsHapEquivariantChainComplex(X) then
return TensorWithZ_Obj(X); fi;

if EvaluateProperty(X,"type") = "resolution" 
   or
   EvaluateProperty(X,"type") = "chaincomplex" then
return TensorWithZ_Obj(X); fi;

if EvaluateProperty(X,"type") = "equivariantChainMap" then
return TensorWithZ_Arr(X); fi;

Print("ERROR: Input should be a resolution or equivariant map between resolutions. \n");

end);
#####################################################################
#####################################################################
#####################################################################



#####################################################################
#####################################################################
InstallGlobalFunction(FilteredTensorWithIntegers,
function(R)
local
	C, ln;

if not IsBound(R!.filteredDimension) 
   or 
   not IsInt(EvaluateProperty(R,"filtration_length"))
then Print("This function must be applied to a filtered resolution.\n");
return fail;
fi;

ln:=EvaluateProperty(R,"filtration_length");

C:=TensorWithIntegers(R);

return Objectify(HapFilteredChainComplex,
rec(
                dimension:=C!.dimension,
                filteredDimension:=R!.filteredDimension,
                boundary:=C!.boundary,
                pseudoBoundary:=R!.pseudoBoundary,
                properties:=
                [["length",EvaluateProperty(R,"length")],
                ["connected",true],
                ["type", "chainComplex"],
                ["filtration_length",ln],
                ["initial_inclusion",EvaluateProperty(R,"initial_inclusion")],
                ["characteristic",
                EvaluateProperty(R,"characteristic")] ]));

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(FilteredTensorWithIntegersModP,
function(R,prime)
local
        C, ln;

if not IsBound(R!.filteredDimension)
   or
   not IsInt(EvaluateProperty(R,"filtration_length"))
then Print("This function must be applied to a filtered resolution.\n");
return fail;
fi;


ln:=EvaluateProperty(R,"filtration_length");

C:=TensorWithIntegersModP(R,prime);

return Objectify(HapFilteredChainComplex,
rec(
                dimension:=C!.dimension,
                filteredDimension:=R!.filteredDimension,
                boundary:=C!.boundary,
                pseudoBoundary:=R!.pseudoBoundary,
                properties:=
                [["length",EvaluateProperty(R,"length")],
                ["connected",true],
                ["type", "chainComplex"],
                ["filtration_length",ln],
                ["initial_inclusion",EvaluateProperty(R,"initial_inclusion")],
                ["characteristic",prime]
                 ]));

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(TensorWithIntegersSparse,
function(R)
local
        BoundaryRec,
        BoundaryC,
        LengthC;



LengthC:=EvaluateProperty(R,"length");
BoundaryRec:=List([1..LengthC],n->[] );


#####################################################################
BoundaryC:=function(n,k)
local  bound, x, i;

if n <0 then return false; fi;
if n=0 then return []; fi;

#if not IsBound(BoundaryRec[n][k]) then 
if true then 
bound:=StructuralCopy(R!.boundary(n,k));
Apply(bound,x->[x[1],1]);
bound:=AlgebraicReduction(bound);
Apply(bound,x->[AbsInt(x[1]),SignInt(x[1])]);
bound:=Collected(bound);
Apply(bound,x->[x[1][1],x[1][2]*x[2]]);
BoundaryRec[n][k]:=bound;
fi;

return BoundaryRec[n][k];
end;
#####################################################################


return  Objectify(HapSparseChainComplex,
                rec(
                dimension:=R!.dimension,
                boundary:=BoundaryC,
                properties:=
                [["length",LengthC],
                ["connected",true],
                ["type", "chainComplex"],
                ["characteristic",
                EvaluateProperty(R,"characteristic")] ]));

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(SparseChainComplexToChainComplex,
function(C)
local Boundary,zero, char;

char:=EvaluateProperty(C,"characteristic");
if char<=0 then zero:=0; else zero:=Zero(GF(char)); fi;

########################################
Boundary:=function(n,k)
local B, x;
B:=[1..C!.dimension(n-1)]*zero;
for x in C!.boundary(n,k) do
B[x[1]]:=B[x[1]]+x[2];
od;
return B;
end;
########################################
return  Objectify(HapChainComplex,
                rec(
                dimension:=C!.dimension,
                boundary:=Boundary,
                properties:=
                [["length",Length(C)],
                ["connected",EvaluateProperty(C,"connected")],
                ["type", "chainComplex"],
                ["characteristic", char] ]));

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(ChainComplexToSparseChainComplex,
function(C)
local Boundary,zero, char,bound;

char:=EvaluateProperty(C,"characteristic");
if char<=0 then zero:=0; else zero:=Zero(GF(char)); fi;

bound:=List([1..Length(C)],i->[]);
 
########################################
Boundary:=function(n,k)
local B,v,i;

if not IsBound(bound[n][k]) then 
B:=[];
v:=C!.boundary(n,k); 
for i in [1..Length(v)] do
if not v[i]=0 then Add(B,[i,v[i]]); fi;
od;
bound[n][k]:=B;
fi;

return bound[n][k];
end;
########################################
return  Objectify(HapSparseChainComplex,
                rec(
                dimension:=C!.dimension,
                boundary:=Boundary,
                properties:=
                [["length",Length(C)],
                ["connected",EvaluateProperty(C,"connected")],
                ["type", "chainComplex"],
                ["characteristic", char] ]));

end);
#####################################################################
#####################################################################

