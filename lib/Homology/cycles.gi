#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(CR_IntegralCycleToClass,
function(arg)
local  
	R, n, Dimension, Boundary, 
	M1, M2, row, 
	dim, BasisKerd1, BasisImaged2, Rels,
	Smith, SmithRecord, TorsionCoefficients,
	ColMat, InvColMat, 
	RemoveRowsMat, InsertRowsList, 
	CycleToClass,
	i, j, x, sum;

R:=arg[1];
n:=arg[2];
Dimension:=R!.dimension;
Boundary:=R!.boundary;

if n <0 then return false; fi;
if n=0 then return [0]; fi;

	################CONSTRUCT BOUNDARY MATRICES M1 AND M2########
M1:=[];
M2:=[];

for i in [1..Dimension(n-1)] do
row:=[];
        for j in [1..Dimension(n)] do
        sum:=0;
                for x in Boundary(n,j) do
                if AbsoluteValue(x[1])=i then
                sum := sum + SignInt(x[1]);
                fi;
                od;
        row[j]:=sum;
        od;
M1[i]:=row;
od;

for i in [1..Dimension(n)] do
row:=[];
        for j in [1..Dimension(n+1)] do
        sum:=0;
                for x in Boundary(n+1,j) do
                if AbsoluteValue(x[1])=i then
                sum := sum + SignInt(x[1]);
                fi;
                od;
        row[j]:=sum;
       	od;
M2[i]:=row;
od;

	################MATRICES M1 AND M2 CONSTRUCTED###############

BasisKerd1:=LLLReducedBasis(TransposedMat(M1),"linearcomb").relations;
BasisImaged2:=LLLReducedBasis(TransposedMat(M2)).basis;
dim:=Length(BasisImaged2);

Rels:=[];
for i in [1..dim] do
        Rels[i]:=SolutionMat(BasisKerd1,BasisImaged2[i]);
od;

	
SmithRecord:= SmithNormalFormIntegerMatTransforms(Rels);
Smith:=SmithRecord.normal;
ColMat:=TransposedMat(SmithRecord.coltrans);
InvColMat:=Inverse(ColMat);

TorsionCoefficients:=[];
for i in [1..Length(BasisKerd1)] do
if i<=Length(Smith) then
    TorsionCoefficients[i]:=Smith[i][i];
else
    TorsionCoefficients[i]:=0;
fi;
od;

InsertRowsList:=[];
RemoveRowsMat:=IdentityMat(Length(TorsionCoefficients));
for i in [1..Length(BasisKerd1)] do
	if TorsionCoefficients[i]=1 then 
	RemoveRowsMat[i]:=47;
	Append(InsertRowsList,[i]);
	fi;
od;
RemoveRowsMat:=Filtered(RemoveRowsMat,r->not (r=47));
if Length(RemoveRowsMat)=0 then
TorsionCoefficients:=[];
else
TorsionCoefficients:= RemoveRowsMat*TorsionCoefficients;
fi;

#####################################################################
CycleToClass:=function(v)
local u, i;

u:=SolutionMat(BasisKerd1,v);
u:=ColMat*u;
u:=RemoveRowsMat*u;

for i in [1..Length(u)] do
u[i]:=u[i]  mod TorsionCoefficients[i];
od;

return u;
end;
#####################################################################




return 	
	 	CycleToClass ;


end);
#####################################################################
#####################################################################


