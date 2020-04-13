#(C) Pablo Fernandez Ascariz, 2005-2006

#####################################################################
#####################################################################
#####################################################################
InstallGlobalFunction(ChevalleyEilenbergComplex,
function(L,x)

local ChevalleyComplex,ChevalleyMap;

# L: Lie algebra or Lie algebras morphism
# x: Length of the chain complex/map to be calculated


###################################################################
###################################################################
ChevalleyComplex:=function(A,s)

local Sctab,B,Dim,n,i,d,j,Boundary,bound,r,Comb,ITT,TTI,ONE,Charact;

B:=Basis(A);
Sctab:=StructureConstantsTable(B);
d:=Length(B);
Comb:=List([1..s],r->Combinations([1..d],r));
#################################################################

################################################################
Charact:=function(M)
if not IsFinite(M!.LeftActingDomain) then 
	if Name(M!.LeftActingDomain)=Name(Integers) then 
		return 0;
	fi;
	if  Name(M!.LeftActingDomain)=Name(Rationals) then 
		return -(1/2);
	fi;
else
 ONE:=One(M!.LeftActingDomain);
	return Characteristic(M!.LeftActingDomain); 
fi;
end;
################################################################

###############################################################
ITT:=function(j,n);
return StructuralCopy(Comb[n][j]);
end;
###############################################################

###############################################################
TTI:=function(n);
return Position(Comb[Length(n)],SSortedList(n));
end;
###############################################################

################################################################
Dim:=function(n);
if n>s then return fail; else
return Binomial( d, n ); fi;
end;
###############################################################

###############################################################
Boundary:=function(n,j)
#Boundary returns the image of the j-th generator of C(n)

local a,b,x,p,q,R,m,Q;

if n>s then
	return fail;
fi;
if j>Binomial(d,n) then
	return fail;
fi;

bound:=List([1..Binomial(d,n-1)], i->0);
Q:=ITT(j,n);

for a in [1..n-1] do
for b in [a+1..n] do
	p:=Length(Sctab[Q[a]][Q[b]][1]);
	for m in [1..p] do
		R:=StructuralCopy(Q);
		Add(R,Sctab[Q[a]][Q[b]][1][m],1);
		Remove(R,a+1);
		Remove(R,b);
		if IsSSortedList(SortedList(R))=true then
			q:=Position(SortedList(R),R[1]);
			bound[TTI(R)]:=bound[TTI(R)]+(-1)^(a+b+q-1)*Sctab[Q[a]][Q[b]][2][m];
		fi;		
	od;
od;
od;

if IsPrimeInt(Charact(A)) then 
return bound*ONE;
else
return bound; fi;

end;
#############################################################

return Objectify(HapChainComplex,rec(dimension:=Dim,
           boundary:=Boundary,
           properties:=
                [["length",s],
                 ["reduced",true],
                 ["type","chainComplex"],
                 ["characteristic",Charact(A)]]));
end;
#############################################################
#############################################################





#############################################################
#############################################################
ChevalleyMap:=function(A,s)

local Sour,Ran,DimSour,DimRan,Charact,Map,Comb,CombSour,ITT,TTI,j,n,r,BasisSour,a,BasisRan,Applicat,SourComp,RanComp,Chara,ONE,Arr;

Sour:=Source(A);
Ran:=Range(A);
BasisSour:=Basis(Sour);
BasisRan:=Basis(Ran);
DimSour:=Length(BasisSour);
DimRan:=Length(BasisRan);
CombSour:=List([1..s],r->Combinations([1..DimSour],r));
Comb:=List([1..s],r->Combinations([1..DimRan],r));
SourComp:=ChevalleyEilenbergComplex(Sour,s);
RanComp:=ChevalleyEilenbergComplex(Ran,s);
Arr:=List([1..s],r->Arrangements([1..DimRan],r));
############################################################

################################################################
Chara:=function(M)
if not IsFinite(M!.LeftActingDomain) then 
	if Name(M!.LeftActingDomain)=Name(Integers) then 
		return 0;
	fi;
	if  Name(M!.LeftActingDomain)=Name(Rationals) then 
		return -(1/2);
	fi;
else
 	ONE:=One(M!.LeftActingDomain);
	return Characteristic(M!.LeftActingDomain); 
fi;
end;
################################################################

###############################################################
ITT:=function(j,n);
return StructuralCopy(CombSour[n][j]);
end;
###############################################################

###############################################################
TTI:=function(n);
return Position(Comb[Length(n)],SSortedList(n));
end;
###############################################################

#############################################################
Map:=function(n,i)

local Mapping,j,k,coef,count,perm;

if n>s then
	return fail;
fi;

Mapping:=List([1..Binomial(DimRan,n)], i->0);
for j in Arr[n] do
	coef:=1;
	count:=0;
	for k in j do
		count:=count+1;
		a:=BasisSour[ITT(i,n)[count]];
		coef:=Coefficients(BasisRan,Image(A,a))[k]*coef;
	od;
	perm:=SortingPerm(j);
	Mapping[TTI(SortedList(j))]:=Mapping[TTI(SortedList(j))]+SignPerm(perm)*coef;
od;
return Mapping;
end;
############################################################

############################################################
Applicat:=function(v,n)
local w;

if not Length(v)=SourComp!.dimension(n) then
	return fail;
fi;
w:=Sum([1..SourComp!.dimension(n)], i->v[i]*Map(n,i));
return w;

end;
############################################################

return Objectify(HapChainMap,rec(source:=SourComp,
		target:=RanComp, 
		mapping:=Applicat,
		properties:=
	     	     [["type","chainMap"],
	     		["characteristic",Chara(Sour)]]));
end;
############################################################
###########################################################

#if IsLieAlgebra(L)=true then
if IsLeftModule(L)=true then
	return ChevalleyComplex(L,x);
else if IsLeftModuleHomomorphism(L) then
	return ChevalleyMap(L,x);
else
return fail;
fi;fi;
end);

#####################################################################
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(LieAlgebraHomology,
function(A,n);

if not IsLeftModule(A) then return fail; fi;
return
Homology(ChevalleyEilenbergComplex(A,n+1),n);

end);
#####################################################################
#####################################################################
