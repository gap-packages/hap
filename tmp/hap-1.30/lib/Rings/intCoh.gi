#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(CR_CocyclesAndCoboundaries,
function(arg)
local  
	R, n, toggle, Dimension, Boundary, 
	M1, M2, row, 
	dim, BasisKerd1, BasisImaged2, Rels,
	Smith, SmithRecord, TorsionCoefficients,
	ColMat, InvColMat, 
	RemoveRowsMat, InsertRowsList, 
	CycleToClass, ClassToCycle,
	i, j, x, sum;

R:=arg[1];
n:=arg[2];
if Length(arg)>2 then toggle := arg[3]; else toggle := false; fi;
Dimension:=R!.dimension;
Boundary:=R!.boundary;

if n <0 then return false; fi;
if n=0 then return [0]; fi;

	################CONSTRUCT BOUNDARY MATRICES M1 AND M2########
M1:=[];
M2:=[];

for i in [1..Dimension(n)] do
row:=[];
        for j in [1..Dimension(n-1)] do
        sum:=0;
                for x in Boundary(n,i) do
                if AbsoluteValue(x[1])=j then
                sum := sum + SignInt(x[1]);
                fi;
                od;
        row[j]:=sum;
        od;
M1[i]:=row;
od;

if Dimension(n+1)>0 then
for i in [1..Dimension(n+1)] do
row:=[];
        for j in [1..Dimension(n)] do
        sum:=0;
                for x in Boundary(n+1,i) do
                if AbsoluteValue(x[1])=j then
                sum := sum + SignInt(x[1]);
                fi;
                od;
        row[j]:=sum;
       	od;
M2[i]:=row;
od;

else

row:=[];
for j in [1..Dimension(n)] do
row[j]:=0;
od;
M2[1]:=row;
fi;
	################MATRICES M1 AND M2 CONSTRUCTED###############

BasisKerd1:=LLLReducedBasis(TransposedMat(M2),"linearcomb").relations;
BasisImaged2:=LLLReducedBasis(TransposedMat(M1)).basis;
dim:=Length(BasisImaged2);

Rels:=[];
for i in [1..dim] do
        Rels[i]:=SolutionMat(BasisKerd1,BasisImaged2[i]);
od;

if Length(Rels)=0 and Length(BasisKerd1)>0 then
Append(Rels,[List([1..Length(BasisKerd1)],x->0)]);
fi;				#CHECK THE MATHS HERE!

if toggle=false then 
return rec(
		cocyclesBasis:=BasisKerd1, 
		boundariesCoefficients:=Rels,
		torsionCoefficients:=fail,
		cocycleToClass:=fail,
		classToCocycle:=fail );
fi;
	################STOP HERE IF TOGGLE=FALSE####################


SmithRecord:= SmithNormalFormIntegerMatTransforms(Rels);
Smith:=SmithRecord.normal;
ColMat:=TransposedMat(SmithRecord.coltrans);
InvColMat:=Inverse(ColMat);	#Only valid for finite groups

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
if TorsionCoefficients[i]>0 then
u[i]:=u[i]  mod TorsionCoefficients[i];
fi;
od;

return u;
end;
#####################################################################

#####################################################################
ClassToCycle:=function(u)
local v,w, i, temp;

for i in [1..Length(u)] do
if TorsionCoefficients[i]>0 then
u[i]:=u[i] mod TorsionCoefficients[i];     
fi;
od;

v :=[];
temp:=0;
for i in [1..Length(BasisKerd1)] do
if i in InsertRowsList then v[i]:=0; 
else 
temp:=temp+1;
v[i] := u[ temp ];
fi;
od;

v:=InvColMat*v;

w:=[];
for i in [1..Dimension(n)] do
w[i]:=0;
od;

for i in [1..Length(v)] do
w:=w + v[i]*BasisKerd1[i];
od;

return w;
end;
#####################################################################

return 	rec(
		cocyclesBasis:=BasisKerd1,
	 	boundariesCoefficients:=Rels,
	 	torsionCoefficients:=TorsionCoefficients,
	 	cocycleToClass:=CycleToClass,        
	 	classToCocycle:=ClassToCycle );

end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(CR_IntegralCohomology,
function(R,n)
local A, i, Smith, TorsionCoefficients;

A:=CR_CocyclesAndCoboundaries(R,n);
Smith:= SmithNormalFormIntegerMat(A.boundariesCoefficients);

TorsionCoefficients:=[];

for i in [1..Length(A.cocyclesBasis)] do
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


#####################################################################
#####################################################################
InstallGlobalFunction(CR_IntegralClassToCocycle,
function(arg)
local
	R, u, n, A,
	i;

R:=arg[1];
u:=arg[2];
n:=arg[3];
if Length(arg)>3 then 
A:=arg[4];
else
A:=CR_CocyclesAndCoboundaries(R,n,true);
fi;

return A.classToCocycle(u);
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(CR_IntegralCocycleToClass,
function(arg)
local
        R, v, n, A,
        i;

R:=arg[1];
v:=arg[2];
n:=arg[3];
if Length(arg)>3 then
A:=arg[4];
else
A:=CR_CocyclesAndCoboundaries(R,n,true);
fi;

return A.cocycleToClass(v);
end);
#####################################################################
#####################################################################

