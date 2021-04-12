#(C) Pablo Fernandez Ascariz, 2005-2006

#####################################################################
#####################################################################
#####################################################################
InstallGlobalFunction(LowerCentralSeriesLieAlgebra,
function(X)

local LieAlgebr,LieMap;

#####################################################################
#####################################################################
LieAlgebr:=function(G)

local LCSeries,Lgth,i,AbInvariants,j,aux,QuotLCS,Dim,DimLieAlgebra,Generat,n,ITPair,Bracket,SCTabl,L,NatHom,Homs;

if not IsPcpGroup(G) then 
Print("The group must be a pcp group.");
return fail; fi;

########################################################################
NatHom:=function(i)

local NH;

if IsPolycyclicGroup(G) then
	NH:=NaturalHomomorphism(LCSeries[i],LCSeries[i+1]);
else
	NH:=NaturalHomomorphismByNormalSubgroupNC(LCSeries[i],LCSeries[i+1]);
fi;

return NH;

end;
########################################################################

LCSeries:=LowerCentralSeries(G);
Lgth:=Length(LCSeries);
Dim:=[];
Generat:=[];
Homs:=List([1..Lgth-1],i->NatHom(i));
QuotLCS:=List([1..Lgth-1],i->Range(Homs[i]));
aux:=AbelianInvariants(QuotLCS[1])[1];

for i in [1..Lgth-1] do
	AbInvariants:=AbelianInvariants(QuotLCS[i]);
	Dim[i]:=Length(AbInvariants);
	for j in [1..Length(AbInvariants)] do;
		if not AbInvariants[j]=aux then
			return "The abelian invariants are not all the same";
		fi;
	od;
	Generat[i]:=GeneratorsOfGroup(QuotLCS[i]);
	if not Length(Generat[i])=Dim[i] then
		return "Error finding generators of abelian group";
	fi;
od;

DimLieAlgebra:=Sum([1..Lgth-1],i->Dim[i]);

########################################################################
ITPair:=function(n)

local s,suma,Term,BasisEl;

if n>DimLieAlgebra then
	return [Lgth,1];
fi;

suma:=0;

for s in [1..Lgth-1] do
	suma:=suma+Dim[s];
	if n<=suma then
		Term:=s;
		break;
	fi;
od;

suma:=suma-Dim[Term];
for s in [1..Dim[Term]] do
	if n=suma+s then
		BasisEl:=s;
		break;
	fi;
od;

return [Term,BasisEl];

end;
#######################################################################

#######################################################################
SCTabl:=function(i,j)

local x,y,Bracket,Commut,k,aux,F,Homomor,a,b;

Bracket:=[];
for k in [1..2*DimLieAlgebra] do
	if IsOddInt(k)=true then
		Bracket[k]:=0;
	else
		Bracket[k]:=1;
	fi;
od;
if ITPair(i)[1]+ITPair(j)[1]>Lgth-1 then
	return Bracket;
fi;
if ITPair(i+j)[1]>Lgth-1 then
	return Bracket;
fi;
aux:=Sum([1..ITPair(i)[1]+ITPair(j)[1]-1],k->Dim[k]);
x:=ITPair(i);
y:=ITPair(j);
a:=PreImagesRepresentative(Homs[x[1]],Generat[x[1]][x[2]]);
b:=PreImagesRepresentative(Homs[y[1]],Generat[y[1]][y[2]]);

Commut:=a*b*Inverse(a)*Inverse(b);
Commut:=Image(Homs[ITPair(i)[1]+ITPair(j)[1]],Commut);
F:=FreeGroup(Length(Generat[ITPair(i)[1]+ITPair(j)[1]]));
Homomor:=GroupHomomorphismByImagesNC(F,QuotLCS[ITPair(i)[1]+ITPair(j)[1]],GeneratorsOfGroup(F),Generat[ITPair(i)[1]+ITPair(j)[1]]);
Commut:=PreImagesRepresentative(Homomor,Commut);
Commut:=LetterRepAssocWord(Commut);

for k in [1..Length(Commut)] do
	Bracket[2*(aux+AbsInt(Commut[k]))-1]:=Bracket[2*(aux+AbsInt(Commut[k]))-1]+SignInt(Commut[k]);
	Bracket[2*(aux+AbsInt(Commut[k]))]:=aux+AbsInt(Commut[k]);
od;

return Bracket;

end;
#####################################################################

L:=EmptySCTable(DimLieAlgebra,0,"antisymmetric");
for i in [1..DimLieAlgebra] do
for j in [i..DimLieAlgebra] do
	SetEntrySCTable(L,i,j,SCTabl(i,j));
od;od;

if aux=0 then
	L:=LieAlgebraByStructureConstants(Integers,L);
else
	L:=LieAlgebraByStructureConstants(GF(aux),L);
fi;

return L;

end;
#####################################################################
#####################################################################






####################################################################
####################################################################
LieMap:=function(f)

local Map,Sour,Ran,BasisSour,i,LCSeriesSour,LgthSour,HomsSour,QuotLCSSour,GeneratSour,DimSour,ITPairSour,Imag,
LCSeriesRan,LgthRan,HomsRan,QuotLCSRan,GeneratRan,DimRan,BasisRan;

Sour:=LowerCentralSeriesLieAlgebra(Source(f));
Ran:=LowerCentralSeriesLieAlgebra(Range(f));
BasisSour:=Basis(Sour);
BasisRan:=Basis(Ran);
LCSeriesSour:=LowerCentralSeries(Source(f));
LCSeriesRan:=LowerCentralSeries(Range(f));
LgthSour:=Length(LCSeriesSour);
LgthRan:=Length(LCSeriesRan);
HomsSour:=List([1..LgthSour-1],i->(NaturalHomomorphism(LCSeriesSour[i],LCSeriesSour[i+1])));
HomsRan:=List([1..LgthRan-1],i->(NaturalHomomorphism(LCSeriesRan[i],LCSeriesRan[i+1])));
QuotLCSSour:=List([1..LgthSour-1],i->Range(HomsSour[i]));
QuotLCSRan:=List([1..LgthRan-1],i->Range(HomsRan[i]));
GeneratSour:=[];
GeneratRan:=[];
DimSour:=[];
DimRan:=[];
Imag:=[];
for i in [1..LgthSour-1] do
	GeneratSour[i]:=GeneratorsOfGroup(QuotLCSSour[i]);
	DimSour[i]:=Length(GeneratSour[i]);
od;
for i in [1..LgthRan-1] do
	GeneratRan[i]:=GeneratorsOfGroup(QuotLCSRan[i]);
	DimRan[i]:=Length(GeneratRan[i]);
od;

########################################################################
ITPairSour:=function(n)

local s,suma,Term,BasisEl;

suma:=0;

for s in [1..LgthSour-1] do
	suma:=suma+DimSour[s];
	if n<=suma then
		Term:=s;
		break;
	fi;
od;

suma:=suma-DimSour[Term];
for s in [1..DimSour[Term]] do
	if n=suma+s then
		BasisEl:=s;
		break;
	fi;
od;

return [Term,BasisEl];

end;
#######################################################################

####################################################################
Map:=function(n)

local Preim,Img,aux,F,Homomor,k,coef,map,sum,j;

#if n>Length(Basis(Sour)) then
#	return fail;
#fi;
if ITPairSour(n)[1]>LgthRan-1 then
	return 0*BasisRan[1];
fi;

sum:=0;
aux:=ITPairSour(n);
map:=One(Ran!.LeftActingDomain);
coef:=List([1..DimRan[aux[1]]],i->0);
Preim:=PreImagesRepresentative(HomsSour[aux[1]],GeneratSour[aux[1]][aux[2]]);
Img:=Image(f,Preim);

Img:=Image(HomsRan[aux[1]],Img);
F:=FreeGroup(DimRan[aux[1]]);
Homomor:=GroupHomomorphismByImagesNC(F,QuotLCSRan[aux[1]],GeneratorsOfGroup(F),GeneratRan[aux[1]]);
Img:=PreImagesRepresentative(Homomor,Img);
Img:=LetterRepAssocWord(Img);

if Img=[] then
	return 0*BasisRan[1];
fi;

for k in [1..Length(Img)] do
	coef[AbsInt(Img[k])]:=coef[AbsInt(Img[k])]+SignInt(Img[k]);
od;

for k in [1..aux[1]-1] do
	sum:=sum+DimRan[k];
od;

for k in [1..DimRan[aux[1]]] do
if not coef[k]=0 then
for j in [1..coef[k]] do
	map:=map*BasisRan[sum+k];
od;fi;od;

return map;

end;
###################################################################

for i in [1..Length(Basis(Sour))] do
	Imag[i]:=Map(i);
od;
return LeftModuleHomomorphismByImagesNC(Sour,Ran,BasisSour,Imag);
end;
###################################################################
###################################################################

if IsGroup(X)=true then
	return LieAlgebr(X);
else if IsGroupHomomorphism(X)=true then
	return LieMap(X);
else
	return fail;
fi;fi;
end);
###################################################################
###################################################################
###################################################################


