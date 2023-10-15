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
                TensorCochainComplex,
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
if EvaluateProperty(R!.tensorWithIntModPRec,"characteristic")=prime then
return R!.tensorWithIntModPRec; 
fi;
fi;

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

#####################################################################
BoundaryC:=function(n,k)
local returnvec, bound, x, i;

if n <0 then return false; fi;
if n=0 then return [0]*One; fi;
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

return returnvec*One;
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

##############
if not IsPrimeInt(prime) then
D!.boundary:=StructuralCopy(C!.boundary);
return
Objectify(HapChainComplex, D);
fi;
##############

one:=One(GF(prime));
oldboundary:=StructuralCopy(C!.boundary);
#############################
newboundary:=function(n,i);
return oldboundary(n,i)*one;
end;
#############################
D!.boundary:=newboundary;

if not IsHapFilteredChainComplex(C) then
return 
Objectify(HapChainComplex, D);
fi;

D!.filteredDimension:=C!.filteredDimension;
return
Objectify(HapFilteredChainComplex, D);

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


####################################################################
####################################################################
TensorCochainComplex:=function(C,prime)
local D, pos, oldboundary, newboundary, newproperties, one;

if not EvaluateProperty(C,"characteristic") in [0,prime] then
Print("Only characteristic 0 or p cochain complexes can be tensored with a prime p field.\n");
return fail;
fi;

D:=rec();
D!.dimension:=StructuralCopy(C!.dimension);
D!.properties:=List(C!.properties,a->StructuralCopy(a));
pos:=PositionProperty(C!.properties,a->"characteristic" in a);
D!.properties[pos][2]:=prime;

##############
if not IsPrimeInt(prime) then
D!.boundary:=StructuralCopy(C!.boundary);
return
Objectify(HapCochainComplex, D);
fi;
##############

one:=One(GF(prime));
oldboundary:=StructuralCopy(C!.boundary);
#############################
newboundary:=function(n,i);
return oldboundary(n,i)*one;
end;
#############################
D!.boundary:=newboundary;

return
Objectify(HapCochainComplex, D);

end;
####################################################################
####################################################################

if IsHapResolution(X)  then 
return Tensor_Obj(X,prime); fi;

if IsHapEquivariantChainMap(X) then 
return Tensor_Arr(X,prime); fi;

if IsHapChainComplex(X) then
return TensorChainComplex(X,prime); fi;

if IsHapCochainComplex(X) then
return TensorCochainComplex(X,prime); fi;

if IsHapChainMap(X) then
return TensorChainMap(X,prime); fi;

return fail;
end);
#####################################################################
#####################################################################
#####################################################################



#####################################################################
#####################################################################
InstallGlobalFunction(TensorWithIntegersModPSparse,
function(R,prime)
local TensorResolution, TensorChainComplex, S;

#########################################################
#########################################################
TensorResolution:=function(R,prime)
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
local  bound, x, i;

if n <0 then return false; fi;
if n=0 then return []; fi;

bound:=StructuralCopy(R!.boundary(n,k));
Apply(bound,x->[x[1],One]);
bound:=AlgebraicReduction(bound);
Apply(bound,x->[AbsInt(x[1]),SignInt(x[1])]);
bound:=Collected(bound);
Apply(bound,x->[x[1][1],x[1][2]*x[2]]);

return bound;
end;
#####################################################################



if EvaluateProperty(R,"characteristic")=0
or EvaluateProperty(R,"characteristic")=prime
then Charact:=prime;
else
Print("ERROR: You probably entered the wrong prime. \n");
return fail; fi;

return       Objectify(HapSparseChainComplex,
                rec(
                dimension:=R!.dimension,
                boundary:=BoundaryC,
                properties:=
                [["length",LengthC],
                ["connected",true],
                ["type", "chainComplex"],
                ["characteristic", Charact]
                 ]));
end;
###############################################
###############################################

###############################################
###############################################
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

##############
if not IsPrimeInt(prime) then
D!.boundary:=StructuralCopy(C!.boundary);
return
Objectify(HapSparseChainComplex, D);
fi;
##############

one:=One(GF(prime));
oldboundary:=StructuralCopy(C!.boundary);
#############################
newboundary:=function(n,i)
local L;
L:=List(oldboundary(n,i), x->[x[1],x[2]*one]);
L:=Filtered(L,x->not IsZero(x[2]));
return L;
end;
#############################
D!.boundary:=newboundary;


if not IsHapFilteredSparseChainComplex(C) then
return
Objectify(HapSparseChainComplex, D);
fi;

D!.filteredDimension:=C!.filteredDimension;
return
Objectify(HapFilteredSparseChainComplex, D);

end;
#################################################
#################################################



if IsHapResolution(R) then return TensorResolution(R,prime); fi;
if IsHapSparseChainComplex(R) then return TensorChainComplex(R,prime); fi;
if IsHapFilteredSparseChainComplex(R) then S:=TensorChainComplex(R,prime); 
S!.filteredDimension:=StructuralCopy(R!.filteredDimension); return S; fi;

end);
#####################################################################
#####################################################################


