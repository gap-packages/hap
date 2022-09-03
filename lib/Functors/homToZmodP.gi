#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
#####################################################################
InstallGlobalFunction(HomToIntegersModP,
function(X,prime)
local 	HomToZ_Obj,
	HomToZ_Arr;
	

#####################################################################
#####################################################################
HomToZ_Obj:=function(R)
local  
	BoundaryC,
	LengthC,
	M,
        myOne,
	Charact;

#One:=Elements(GaloisField(prime))[2];
myOne:=One(GF(prime));
LengthC:=EvaluateProperty(R,"length");
M:=[1..LengthC];				

#####################################################################
BoundaryC:=function(N,k)
local
	n,row, Mt, i, j, x, sum;

n:=N+1;
if n <0 then return false; fi;

if n=0 then
#return List([1..R!.dimesnion(1)],x->0); fi;		
return [1..R!.dimesnion(1)]*0; fi;

if M[n]=n then 
   Mt:=[];

   if R!.dimension(n)>0 then

   for i in [1..R!.dimension(n)] do
   row:=[];
        for j in [1..R!.dimension(n-1)] do
        sum:=0;
                for x in R!.boundary(n,i) do
                if AbsoluteValue(x[1])=j then
                sum := sum + SignInt(x[1]);
                fi;
                od;
        row[j]:=sum;
        od;
   Mt[i]:=row;
   od;

   M[n]:=TransposedMat(Mt);

   else

   row:=[];
   for j in [1..R!.dimension(n-1)] do
   row[j]:=0;
   od;
   for i in [1..R!.dimension(n-1)] do
   Append(Mt,[row]);
   od;
   M[n]:=Mt;
   fi;

fi;

return M[n][k]*myOne;
end;
#####################################################################

if EvaluateProperty(R,"characteristic")=0
or EvaluateProperty(R,"characteristic")=prime
then Charact:=prime;
else
Print("ERROR: You probably entered the wrong prime. \n");
return fail; fi;


return Objectify(HapCochainComplex,
		rec(
		dimension:=R!.dimension,
		boundary:=BoundaryC,
		properties:=
		[["length",LengthC],
		["connected",true],
		["type", "cochainComplex"],
		["characteristic", Charact] 
		]));
end;
#####################################################################
#####################################################################



#####################################################################
#####################################################################
HomToZ_Arr:=function(F)
local
                R,S,RhomS,              #R->S is an equivariant chain
                                        #map. C<-D is the chain map
                C,D,DhomC,              #got by Homing.
                DimensionC,
                DimensionD,             #Throughout the program we
                x;                      #identify R and S with their
R:=F!.source;                           #duals Hom(R,Z).
S:=F!.target;
RhomS:=F!.mapping;
C:=HomToZ_Obj(R);
D:=HomToZ_Obj(S);
DimensionC:=C!.dimension;
DimensionD:=D!.dimension;

#####################################################################
DhomC:=function(v,n)
local
                u, i,j,temp,x;
#u:=List([1..DimensionC(n)],x->0);
u:=[1..DimensionC(n)]*0;

for i in [1..DimensionD(n)] do
        for j in [1..DimensionC(n)] do
        temp:=0;
        for x in List(RhomS([[j,1]],n),y->y[1]) do
        if x=i then temp:=temp+1;fi;
        if x=-i then temp:=temp-1; fi;
        od;
        u[j]:=u[j]+temp*v[i];
        od;
od;
return u;
end;
#####################################################################



return Objectify(HapCochainMap,
        rec(
           source:=D,
           target:=C,
           mapping:=DhomC,
           properties:=[ ["type","cochainMap"],
           ["characteristic", prime
           ]
           ]));
end;
#####################################################################
#####################################################################





if EvaluateProperty(X,"type") = "resolution" then
return HomToZ_Obj(X); fi;

if EvaluateProperty(X,"type") = "equivariantChainMap" then
return HomToZ_Arr(X); fi;

if EvaluateProperty(X,"type") = "chainComplex" then
return HAP_HomToIntModP_ChainComplex(X,prime); fi;

if EvaluateProperty(X,"type") = "chainMap" then
return HAP_HomToIntModP_ChainMap(X,prime); fi;

if EvaluateProperty(X,"type") = "cochainComplex" then
return HAP_HomToIntModP_CochainComplex(X,prime); fi;

if EvaluateProperty(X,"type") = "cochainMap" then
return HAP_HomToIntModP_CochainMap(X,prime); fi;

return fail;

Print("ERROR: Input should be a resolution or equivariant map between resolutions. \n");

end);
#####################################################################
#####################################################################
#####################################################################

