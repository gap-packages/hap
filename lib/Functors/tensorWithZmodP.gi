#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
#####################################################################
InstallGlobalFunction(TensorWithIntegersModP,
function(X,prime)
local 		
		Tensor_Obj,
		Tensor_Arr,
		TensorChainComplex,
	        TensorChainMap;

#####################################################################
#####################################################################
Tensor_Obj:= function(R,prime)
local  
	BoundaryC,
	LengthC,
	M,
	One,
	Charact;

if "tensorWithIntModPRec" in NamesOfComponents(R) then
return R!.tensorWithIntModPRec; fi;

One:=Elements(GaloisField(prime))[2];
LengthC:=EvaluateProperty(R,"length");
M:=[1..LengthC];				

#####################################################################
BoundaryC:=function(n,k)
local
	row, Mt, i, j, x, sum;

if n <0 then return false; fi;
if n=0 then return [0]; fi;

if M[n]=n then 
   Mt:=[];

   for i in [1..R!.dimension(n-1)] do
   row:=[];
        for j in [1..R!.dimension(n)] do
        sum:=0;
                for x in R!.boundary(n,j) do
                if AbsoluteValue(x[1])=i then
                sum := sum + SignInt(x[1]);
                fi;
                od;
        row[j]:=sum*One;
        od;
   Mt[i]:=row;
   od;

   M[n]:=TransposedMat(Mt);
fi;

return M[n][k];
end;
#####################################################################

if EvaluateProperty(R,"characteristic")=0
or EvaluateProperty(R,"characteristic")=prime
then Charact:=prime;
else
Print("ERROR: You probably entered the wrong prime. \n");
return fail; fi;

R!.tensorWithIntModPRec:= Objectify(HapChainComplex,
		rec(
		dimension:=R!.dimension,
		boundary:=BoundaryC,
		properties:=
		[["length",LengthC],
		["connected",true],
		["type", "chainComplex"],
		["characteristic", Charact]
		 ]));
return R!.tensorWithIntModPRec;
end;
#####################################################################
#####################################################################

#####################################################################
#####################################################################
Tensor_Arr:=function(F,prime)
local

                R,S,RhomS,              #R->S is an equivariant chain
                CmapR, SmapD,           #map. C->D is the chain map
                C,D,ChomD,              #got by killing the actions.
                DimensionS,
                DimensionC,
		FieldToInt, one,
                x;

R:=F!.source;
S:=F!.target;
DimensionS:=S!.dimension;
RhomS:=F!.mapping;
C:=Tensor_Obj(R,prime);
D:=Tensor_Obj(S,prime);
DimensionC:=C!.dimension;

one:=Elements(GaloisField(prime))[2];

#####################################################################
FieldToInt:=function(x)
local
	i;
for i in [0..prime] do
if i*one=x then return i; fi;
od;

end;
#####################################################################

#####################################################################
CmapR:=function(v,n)
local  i,j,w,x;

w:=[];

for i in [1..DimensionC(n)] do
if not FieldToInt(v[i])=0 then
        x:=[i,1];
        for j in [1..FieldToInt(v[i])] do
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

return one*v;
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
           properties:=
	   [ ["type","chainMap"],
             ["characteristic", prime]    ]));
end;
#####################################################################
#####################################################################

#####################################################################
#####################################################################
TensorChainComplex:=function(C,prime)
local D, pos, oldboundary, newboundary, newproperties, one;

if not EvaluateProperty(C,"characteristic") in [0,prime] then
Print("Only characteristic 0 or p chain complexes can be tensored with a prime p field.\n");
return fail;
fi;

D:=rec();
D!.dimension:=StructuralCopy(C!.dimension);
D!.properties:=List(C!.properties,a->StructuralCopy(a));
pos:=PositionProperty(C!.properties,a->"characteristic" in a);
D!.properties[pos][2]:=prime;

one:=One(GF(prime));
oldboundary:=StructuralCopy(C!.boundary);
#############################
newboundary:=function(n,i);
return oldboundary(n,i)*one;
end;
#############################
D!.boundary:=newboundary;

return 
Objectify(HapChainComplex, D);
end;
####################################################################
####################################################################

####################################################################
####################################################################
TensorChainMap:=function(M,prime)
local N, Mapping,oldmapping,one,pos;

N:=rec();
N.source:=TensorWithIntegersModP(M!.source,prime);
N.target:=TensorWithIntegersModP(M!.target,prime);
N!.properties:=List(M!.properties,a->StructuralCopy(a));
pos:=PositionProperty(M!.properties,a->"characteristic" in a);
N!.properties[pos][2]:=prime;

one:=One(GF(prime));
oldmapping:=StructuralCopy(M!.mapping);
####################################
Mapping:=function(v,n);
return oldmapping(v,n)*one;
end;
####################################
N!.mapping:=Mapping;

return Objectify(HapChainMap,N);
end;
####################################################################
####################################################################



if IsHapResolution(X)  then 
return Tensor_Obj(X,prime); fi;

if IsHapEquivariantChainMap(X) then 
return Tensor_Arr(X,prime); fi;

if IsHapChainComplex(X) then
return TensorChainComplex(X,prime); fi;

if IsHapChainMap(X) then
return TensorChainMap(X,prime); fi;


return fail;
end);
#####################################################################
#####################################################################
#####################################################################


