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
local n,k,S,R,b,i,x,dels,len,B,One, M, newbound, BoundaryC,dimension;;

for len in [0..1000] do   ##Sloppy!
if K!.dimension(1+len)=0 then break; fi;
od;

B:=List([0..len],i->[1..K!.dimension(i)]);

dels:=[];
for n in [0..len] do
for k in [1..K!.dimension(n)] do
S:=K!.stabilizer(n,k);
R:=ResolutionFiniteGroup(S,1);
for i in [1..R!.dimension(1)] do
b:=R!.boundary(1,i);
for x in b do
if not R!.elts[x[2]] in K!.elts then
Add(K!.elts,R!.elts[x[2]]);
fi;
if K!.action(n,k,Position(K!.elts,R!.elts[x[2]]))=-1 then
AddSet(dels,[n,k]);
fi;
od;
od;
od;
od;

for x in dels do
RemoveSet(B[x[1]+1],x[2]);
od;

#######################################################
dimension:=function(n)
if n>len then return 0; fi;
return Length(B[n+1]);
end;
#######################################################

#######################################################
newbound:=function(n,k)
local b, x, bnd, p;
b:=K!.boundary(n,B[n+1][k]);
bnd:=[];
for x in b do
p:=Position(B[n],AbsInt(x[1]));
if not p=fail then
Add(bnd,[SignInt(x[1])*p,x[2]]);
fi;
od;

return bnd;
end;
#######################################################

One:=1;
M:=[1..len];

#####################################################################
BoundaryC:=function(n,k)
local
        row, Mt, i, j, x, sum;

if n <0 then return false; fi;
if n=0 then return [0]; fi;
if dimension(n-1)=0 then return [0]; fi;

if M[n]=n then
Mt:=[];

for i in [1..dimension(n-1)] do
row:=[];
   for j in [1..dimension(n)] do
   sum:=0;
          for x in newbound(n,j) do
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
            dimension:=dimension,
            boundary:=BoundaryC,
            properties:=
            [["length",len],
             ["connected",true],
             ["type", "chainComplex"],
             ["characteristic", -1/2]]));


end);
###################################################
###################################################

