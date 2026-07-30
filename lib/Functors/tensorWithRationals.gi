#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(TensorWithRationals,
function(R)
local
        C,D,map,BoundaryC,
        LengthC,
        M,
        One,
        Charact ;

##########################################################
if IsHapChainComplex(R) then
return TensorWithIntegersModP(R,-1/2);
fi;
##########################################################

##########################################################
if IsHapCochainComplex(R) then
return TensorWithIntegersModP(R,-1/2);
fi;
##########################################################

##########################################################
if IsHapCochainMap(R) then
C:=TensorWithRationals(Source(R));
D:=TensorWithRationals(Target(R));
map:=R!.mapping;
return
Objectify(HapCochainMap,
            rec(
            source:=C,
            target:=D,
            mapping:=map,
            properties:=
            [["type","cochainMap"],
             ["characteristic", -1/2]]));
fi;
##########################################################

##########################################################
if IsHapChainMap(R) then
C:=TensorWithRationals(Source(R));
D:=TensorWithRationals(Target(R));
map:=R!.mapping;
return
Objectify(HapChainMap,
            rec(
            source:=C,
            target:=D,
            mapping:=map,
            properties:=
            [["type","chainMap"],
             ["characteristic", -1/2]]));
fi;
##########################################################


##########################################################
if IsHapNonFreeResolution(R) then
return TensorNonFreeResolutionWithRationals(R);
fi;
##########################################################


One:=1;
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

return Objectify(HapChainComplex,
	    rec(
            dimension:=R!.dimension,
            boundary:=BoundaryC,
            properties:=
            [["length",LengthC],
             ["connected",true],
             ["type", "chainComplex"],
             ["characteristic", -1/2]
				     ]));
end);
#####################################################################
#####################################################################
																 
##################################################
##################################################
InstallGlobalFunction(TensorNonFreeResolutionWithRationals,
function(K)
local C,pos;

C:=TensorWithIntegersMod2Torsion(K);
pos:=PositionProperty(C!.properties,x->x[1]="characteristic");
C!.properties[pos][2]:=-1/2;
return C;

end);
###################################################
###################################################

##################################################
##################################################
InstallGlobalFunction(HomToRationals,
function(K)
local C,pos;

C:=TensorWithIntegers(K);
C:=HomToIntegers(C);
pos:=PositionProperty(C!.properties,x->x[1]="characteristic");
C!.properties[pos][2]:=-1/2;
return C;

end);
###################################################
###################################################

