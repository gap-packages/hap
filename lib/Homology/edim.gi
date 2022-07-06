#####################################################################
InstallGlobalFunction(HomologyPrimePart,
function(C,n,p)
local
         M2, Smith, Divs,  Dimension, Boundary,
                i;

if n <0 then return false; fi;

Dimension:=C!.dimension;
Boundary:=C!.boundary;

#######################
M2:=[];
for i in [1..Dimension(n+1)] do
M2[i]:=Boundary(n+1,i);
od;

ConvertToMatrixRep(M2);
Divs:=ElementaryDivisorsPPartRk(M2,p);
M2:=0;
Smith:=[];
for i in [1..Length(Divs)-1] do
Append(Smith,ListWithIdenticalEntries(Divs[i]-Divs[i+1] ,p^i));
od;

#######################

return Smith ;

end);
#####################################################################
#####################################################################


#####################################################################
InstallGlobalFunction(CohomologyPrimePart,
function(C,N,p)
local
 M2, Divs, Smith, Dimension, Boundary, i,row,n;

n:=N-1;

if N <0 then return [ ]; fi;

Dimension:=C!.dimension;
Boundary:=C!.boundary;

M2:=[];

if n=-1 then M2:=[List([1..Dimension(0)],i->0)];
else
for i in [1..Dimension(n)] do
M2[i]:=Boundary(n,i);
od;
fi;

Divs:=ElementaryDivisorsPPartRk(M2,p);
Smith:=[];
for i in [1..Length(Divs)-1] do
Append(Smith,ListWithIdenticalEntries(Divs[i]-Divs[i+1] ,p^i));
od;


return Smith;

end);
#####################################################################
#####################################################################

