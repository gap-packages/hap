#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(IntegralHomologyOfChainComplex,
function(C,n)
local  
	M1, M2, DimKerd1,  Smith, TorsionCoefficients, Dimension, Boundary,
	i;

if n <0 then return false; fi;

Dimension:=C!.dimension;
Boundary:=C!.boundary;


########################
if n=0 then 
DimKerd1:=Dimension(n);

else
M1:=[];

for i in [1..Dimension(n)] do
M1[i]:=Boundary(n,i);
od;

if Dimension(n-1)=0 then M1:=List(M1,x->[0]); fi;

ConvertToMatrixRep(M1);
if Length(M1)=0 then 
DimKerd1:=0;
else
DimKerd1:=Length(M1)-Rank(M1);
fi;
M1:=0;

fi;
#######################

#######################
M2:=[];
for i in [1..Dimension(n+1)] do
M2[i]:=Boundary(n+1,i);
od;

if Dimension(n)=0 then M2:=List(M2,x->[0]); fi;

ConvertToMatrixRep(M2);
Smith:= SmithNormalFormIntegerMat(M2);
M2:=0;
#######################

TorsionCoefficients:=[];

for i in [1..DimKerd1] do
	if i<=Length(Smith) then
	TorsionCoefficients[i]:=Smith[i][i];
	else 
	TorsionCoefficients[i]:=0;
	fi;
od;

return Filtered(TorsionCoefficients, i-> not (i=1)) ;

end);
#####################################################################
#####################################################################



