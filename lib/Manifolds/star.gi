######################################
######################################
InstallGlobalFunction(RemoveStar,
function(K,V)
local maxK, maxSt, max, St, v, L;
# Inputs a simplicial complex K and set V of vertices.
# Returns K - (the union of the stars of the vertices in V).
maxSt:=[];
for v in V do
St:=VertexStar(K,v);
L:=MaximalSimplicesOfSimplicialComplex(St);
L:=List(L,s->SortedList(s));
if Length(Intersection(Flat(L),Flat(maxSt)))>0 then
Print("Oh oh\n"); fi;
Append(maxSt,L);
od;

maxK:=MaximalSimplicesOfSimplicialComplex(K);
maxK:=List(maxK,s->SortedList(s));
max:=Difference(maxK,maxSt);
return MaximalSimplicesToSimplicialComplex(max);

end);
######################################
######################################



######################################
######################################
InstallGlobalFunction(NonManifoldVertices,
function(YY)
local K, KK, Y, L, dim, mx, mn,v,issphr, sing, b;
# For a suitable regular CW-complex YY, or simplicial complex YY, of 
# dimension 3 this function returns the list of vertices whose
# links are not spheres. If K is not suitable then the function
# returns fail. "Suitable" means that YY - (union of stars of the vertices)
# is a manifold with boundary, where barycentric subdivision might
# need to be applied to YY before removal to ensure that the stars are 
# disjoint.

if not (IsHapRegularCWComplex(YY) or IsHapSimplicialComplex(YY)) then
Print("Argument must be a regular CW or a simlicial complex.\n");
return fail;
fi;

if IsHapRegularCWComplex(YY) then
if not IsPureRegularCWComplex(YY) then return fail; fi;
K:=BarycentricSubdivision(YY);
Y:=RegularCWComplex(K);
else
K:=YY;
Y:=RegularCWComplex(K);
if not IsPureRegularCWComplex(Y) then return fail; fi;
fi;

#So K is now a pure simplicial complex and Y is its RegularCWComplex
#representation.

dim:=Dimension(Y);
if dim<>3 then return fail; fi;
L:=List(Y!.coboundaries[dim],x->x[1]);
mx:=Maximum(L);
mn:=Minimum(L);
if not (mx=2 and mn=2) then return fail; fi;


#####################################
issphr:=function(W)
local L,F,H;
L:=CriticalCells(W);
if [0,2]=SortedList(List(L,x->x[1])) then return true; fi;
H:=List([0..3],i->Homology(W,i));
if not H=[[0],[],[],[0]] then return false; fi;
F:=FundamentalGroup(W);
if GeneratorsOfGroup(F)=[] then return true; fi;
if Order(F)=1 then return true; fi;
return false;
end;
#####################################

sing:=[];
for v in K!.vertices do
L:=VertexLink(K,v);
if not IsClosedManifold(L) then return fail; fi;
L:=RegularCWComplex(L);
#L:=CriticalCells(L);
#if not [0,2]=SortedList(List(L,x->x[1])) then return false; fi;
b:=issphr(L);
if not b then Add(sing,v); fi;
od;
return sing;
end);
######################################
######################################

######################################
######################################
InstallGlobalFunction(ThreeManifoldWithBoundary,
function(YY)
local Y,K,V,L,S,bool,v;

if not IsHapRegularCWComplex(YY) then
Print("Argument must be a regular CW-complex.\n");
return fail;
fi;

if IsHapRegularCWComplex(YY) then
if not IsPureRegularCWComplex(YY) then 
Print("Complex is not pure.\n");
return fail; fi;
K:=BarycentricSubdivision(YY);
Y:=RegularCWComplex(K);
else
K:=YY;
Y:=RegularCWComplex(K);
if not IsPureRegularCWComplex(Y) then
Print("Complex is not pure.\n");
return fail; fi;
fi;

V:=NonManifoldVertices(K);
if Length(V)=0 then return YY; fi;

L:=[];
bool:=false;
for v in V do
S:=VertexStar(K,v);
if Size(Intersection(L,S!.vertices))>0 then bool:=true;
break; fi;
Append(L,S!.vertices);
od;

if bool then
K:=BarycentricSubdivision(Y);
fi;

return BarycentricallySimplifiedComplex(RegularCWComplex(RemoveStar(K,V)));

end);
#######################################
######################################
