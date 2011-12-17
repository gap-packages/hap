#######################################################
InstallGlobalFunction(IsLieAlgebraHomomorphism,
function(f)

local 
BasSour,BasRan,i,j,k,m,n,l,a,b,DimSour,DimRan,Vector1,Vector2,TabSour,TabRan,count,Vec;

if not IsMapping(f) then return false; fi;
if not IsLieAlgebra(Source(f)) then  return false; fi;
if not IsLieAlgebra(Range(f)) then  return false; fi;

BasSour:=Basis(Source(f));
BasRan:=Basis(Range(f));
DimSour:=Length(BasSour);
DimRan:=Length(Basis(Range(f)));
TabSour:=StructureConstantsTable(BasSour);
TabRan:=StructureConstantsTable(BasRan);


if not IsLeftModuleHomomorphism(f) then
	return false;
else
	count:=0;
	for i in [1..DimSour-1] do
	for j in [i+1..DimSour] do
		Vector1:=0*BasSour[1];
		Vector2:=List([1..DimRan], r->0);
		for k in [1..Length(TabSour[i][j][1])] do
			Vector1:=Vector1+BasSour[TabSour[i][j][1][k]]*TabSour[i][j][2][k];
		od;
		Vector1:=Image(f,Vector1);
		Vector1:=Coefficients(BasRan,Vector1);

		a:=Coefficients(BasRan,Image(f,BasSour[i]));
		b:=Coefficients(BasRan,Image(f,BasSour[j]));
		for m in [1..DimRan] do
		for n in [1..DimRan] do
			Vec:=List([1..DimRan], r->0);
			for l in [1..Length(TabRan[m][n][1])] do
				Vec[TabRan[m][n][1][l]]:=TabRan[m][n][2][l];
			od;
			Vector2:=Vector2+a[m]*b[n]*Vec;
		od;od;
		if not Vector1=Vector2 then
			count:=1;
		fi;
	od;od;
	if count=1 then
		return false;
	else
		return true;
	fi;
fi;

end);
#######################################################
