######################################################
######################################################
InstallMethod(DualComplex,
"Hom_{ZG}(-,ZG) for  group resolutions",
[IsHapResolution],
function(R)
return DualComplex_FreeResolution(R);
end);
######################################################
######################################################

######################################################
######################################################
InstallOtherMethod(DualComplex,
"Hom_{Z}(-,Z) for non-free group resolutions",
[IsHapNonFreeResolution],
function(R)
return DualComplex_NonFreeResolution(R);
end);
######################################################
######################################################

######################################################
######################################################
InstallGlobalFunction(DualComplex_FreeResolution,
function(K)
local Dimension, Boundary, Stabilizer, Action, dim, BoundaryRec, 
Boundaries, 
id, b, k, n, x, y, g, s, bnd, pos, Y,dm ,cells, B;

###########################
dim:=0;
while K!.dimension(dim)>0 and dim<Length(K) do
dim:=dim+1;
od;
dim:=dim-1;
###########################

###########################
Dimension:=function(n);
if n<0 or n>dim then return 0; fi;
return K!.dimension(dim-n);
end;
###########################

###########################
BoundaryRec:=List([0..dim],i->[]);
for n in [1..dim] do
for k in [1..Dimension(n)] do
BoundaryRec[n][k]:=[];
od;
od;

for n in [1..dim] do
for k in [1..K!.dimension(n)] do
for x in K!.boundary(n,k) do


g:=K!.elts[x[2]]^-1;

pos:=Position(K!.elts,g);
if pos = fail then Add(K!.elts,g); pos:=Length(K!.elts); fi;
y:=[SignInt(x[1])*k,pos];


if not y in BoundaryRec[dim-n+1][AbsInt(x[1])] and not [-y[1],y[2]] in BoundaryRec[dim-n+1][AbsInt(x[1])] then
AddSet(BoundaryRec[dim-n+1][AbsInt(x[1])],y);
fi;
od;
od;
od;
###########################

###########################
Boundary:=function(n,k);
return BoundaryRec[n][k];
end;
###########################


return Objectify(HapResolution,
            rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=fail,
            elts:=K!.elts,
            group:=K!.group,
            properties:=
            [["length",EvaluateProperty(K,"length")],
             ["characteristic",EvaluateProperty(K,"characteristic")],
             ["type","resolution"],
             ["resolutionOfModule",true],
             ["reduced",false]]  ));

end);
######################################################
######################################################

######################################################
######################################################
InstallGlobalFunction(DualComplex_NonFreeResolution,
function(K)
local Dimension, Boundary, Stabilizer, Action, dim, BoundaryRec,
StandardForm, Boundaries,
id, i, b, k, n, x, y, g, s, bnd, pos, poss;

###########################
StandardForm:=function(n,k,g)
local h,s,pos;
h:= Minimum(g*Elements(K!.stabilizer(n,k)));
s:=h^-1*g;
pos:=Position(K!.elts,s);
if pos=fail then Add(K!.elts,s); pos:=Length(K!.elts); fi;
s:=K!.action(n,k,pos);
return [h,s];
end;
###########################

###########################
dim:=0;
while K!.dimension(dim)>0 do
   dim:=dim+1;
od;
dim:=dim-1;
###########################

###########################
Dimension:=function(n);
if n<0 or n>dim then return 0; fi;
return K!.dimension(dim-n);
end;
###########################

###########################
BoundaryRec:=List([0..dim],i->[]);
for n in [1..dim] do
for k in [1..Dimension(n)] do
BoundaryRec[n][k]:=[];
od;
od;

for n in [1..dim] do
for k in [1..K!.dimension(n)] do
for x in K!.boundary(n,k) do
for s in K!.stabilizer(n-1,AbsInt(x[1])) do

   g:=s*K!.elts[x[2]]^-1;
   g:=StandardForm(n,k,g);
poss:=Position(K!.elts,s);
if poss=fail then Add(K!.elts,s); poss:=Length(K!.elts); fi;
   pos:=Position(K!.elts,g[1]);
   if pos = fail then Add(K!.elts,g[1]); pos:=Length(K!.elts); fi;
   y:=[K!.action(n-1,AbsInt(x[1]),poss)*g[2]*SignInt(x[1])*k,pos];

   if not y in BoundaryRec[dim-n+1][AbsInt(x[1])] and 
      not [-y[1],y[2]] in BoundaryRec[dim-n+1][AbsInt(x[1])] then
   AddSet(BoundaryRec[dim-n+1][AbsInt(x[1])],y);
   fi;

od;
od;
od;
od;
###########################

###########################
Boundary:=function(n,k);
if k>0 then
return BoundaryRec[n][k];
else 
return BoundaryRec[n][-k];
fi;
end;
###########################

###########################
Stabilizer:=function(n,k);
return K!.stabilizer(dim-n,k);
end;
###########################

###########################
Action:=function(n,k,g);
return K!.action(dim-n,k,g);
end;
###########################

return Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=fail,
            elts:=K!.elts,
            group:=K!.group,
            stabilizer:=Stabilizer,
            action:=Action,
            properties:=
            [["length",EvaluateProperty(K,"length")],
             ["characteristic",EvaluateProperty(K,"characteristic")],
             ["type","resolution"],
             ["resolutionOfModule",true],
             ["reduced",false]]  ));

end);
######################################################
######################################################

