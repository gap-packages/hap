######################################
######################################
InstallGlobalFunction(IsClosedManifold,
function(YY)
local K, KK, Y, L, dim, mx, mn,v,issphr;
#For a regular CW-complex X, or simplicial complex X, of dimension <4
# this function returns true/false depending on whether X is a closed
# manifold. If dim(X)>3 then the function returns false or fail. 

if not (IsHapRegularCWComplex(YY) or IsHapSimplicialComplex(YY)) then
Print("Argument must be a regular CW or a simlicial complex.\n");
return fail;
fi;

if IsHapRegularCWComplex(YY) then
if not IsPureRegularCWComplex(YY) then return false; fi;
K:=BarycentricSubdivision(YY);
Y:=RegularCWComplex(K);
else
K:=YY;
Y:=RegularCWComplex(K);
if not IsPureRegularCWComplex(Y) then return false; fi;
fi;

#So K is now a pure simplicial complex and Y is its RegularCWComplex
#representation.

dim:=Dimension(Y);
L:=List(Y!.coboundaries[dim],x->x[1]);
mx:=Maximum(L);
mn:=Minimum(L);
if not (mx=2 and mn=2) then return false; fi;

if dim>3 then 
Print("This function is fully implemented only for spaces of dimension <4.\n");
return fail;
fi;


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

if dim=1 or dim=2 then return true; fi;
for v in K!.vertices do
L:=VertexLink(K,v);
if not IsClosedManifold(L) then return false; fi;
L:=RegularCWComplex(L);
#L:=CriticalCells(L);
#if not [0,2]=SortedList(List(L,x->x[1])) then return false; fi;
if not issphr(L) then return false; fi;
od;

return true;
end);
######################################
######################################

