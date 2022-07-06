#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(IntegralCohomologyOfCochainComplex,
function(C,N)
local  
	M1, M2, 
	DimKerd1,  
	 Smith, TorsionCoefficients,
	Dimension, Boundary,
	i,row,n;

n:=N-1;

if N <0 then return [ ]; fi;
if C!.dimension(N)=0 then return []; fi; #ADDED FEB 2014

Dimension:=C!.dimension;
Boundary:=C!.boundary;
M1:=[];
M2:=[];

if n=-1 then M2:=[List([1..Dimension(0)],i->0)];
else
for i in [1..Dimension(n)] do
M2[i]:=Boundary(n,i);
od;
fi;

for i in [1..Dimension(n+1)] do
M1[i]:=Boundary(n+1,i);
od;

if Length(M1)=0 then
DimKerd1:=0;
else
DimKerd1:=Length(M1)-Rank(M1);
fi;

Smith:= SmithNormalFormIntegerMat(M2);

TorsionCoefficients:=[];

for i in [1..DimKerd1] do
	if i<=Length(Smith) then
	TorsionCoefficients[i]:=Smith[i][i];
	else 
	TorsionCoefficients[i]:=0;
	fi;
od;

return Filtered(TorsionCoefficients, i-> not (i=1));

end);
#####################################################################
#####################################################################



