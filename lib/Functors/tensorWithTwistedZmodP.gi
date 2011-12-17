#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
#####################################################################
InstallGlobalFunction(TensorWithTwistedIntegersModP,
function(X,prime,rho)
local 		
		Tensor_Obj,
		Tensor_Arr;

#####################################################################
#####################################################################
Tensor_Obj:= function(R,prime,rho)
local  
	BoundaryC,
	LengthC,
	M,
	One,
	Charact,
	Bool;

Bool:=false;
if "tensorWithTwistedIntModPRec" in NamesOfComponents(R) then
if R!.tensorWithTwistedIntModPRec[2]=rho 
and R!.tensorWithTwistedIntModPRec[3]=prime then Bool:=true; fi;
fi;
if Bool then return R!.tensorWithTwistedIntModPRec[1]; fi;

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
                sum := sum + SignInt(x[1])*rho(R!.elts[x[2]]);
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

R!.tensorWithTwistedIntModPRec:=[];
R!.tensorWithTwistedIntModPRec[1]:=
 Objectify(HapChainComplex,
		rec(
		dimension:=R!.dimension,
		boundary:=BoundaryC,
		twist:=rho,
		properties:=
		[["length",LengthC],
		["connected",true],
		["type", "chainComplex"],
		["characteristic", Charact]
		 ]));
R!.tensorWithTwistedIntModPRec[2]:=rho;
R!.tensorWithTwistedIntModPRec[3]:=prime;
return R!.tensorWithTwistedIntModPRec[1];
end;
#####################################################################
#####################################################################

#####################################################################
#####################################################################
Tensor_Arr:=function(F,prime,rho)
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
C:=Tensor_Obj(R,prime,rho);
D:=Tensor_Obj(S,prime,rho);
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
v[AbsoluteValue(x[1])]:=v[AbsoluteValue(x[1])]+SignInt(x[1])*rho(S!.elts[x[2]]);
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
	   twist:=rho,
           properties:=
	   [ ["type","chainMap"],
             ["characteristic", prime]    ]));
end;
#####################################################################
#####################################################################


if EvaluateProperty(X,"type")="resolution" then
return Tensor_Obj(X,prime,rho); fi;

if EvaluateProperty(X,"type")="equivariantChainMap" then
return Tensor_Arr(X,prime,rho); fi;


return fail;
end);
#####################################################################
#####################################################################
#####################################################################



