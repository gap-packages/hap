#(C) Pablo Fernandez Ascariz, 2005-2006

#####################################################################
#####################################################################
#####################################################################
InstallGlobalFunction(LeibnizComplex,
function(arg)

local Leibnizcomplex,LeibnizMap, L, x, sparse;

# L: Leibniz algebra or Leibniz algebra morphism
# x: Length of the chain complex/map to be calculated
L:=arg[1];
x:=arg[2];
if Length(arg)=3 then sparse:=arg[3]; else sparse:=false; fi;

###################################################################
###################################################################
Leibnizcomplex:=function(A,s)

local Sctab,B,Dim,n,i,d,j,Boundary,bound,r,Tup,ITT,TTI,ONE,Charact;

B:=Basis(A);
Sctab:=StructureConstantsTable(B);
d:=Length(B);
Tup:=List([1..s],r->Tuples([1..d],r));
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
return StructuralCopy(Tup[n][j]);
end;
###############################################################

###############################################################
TTI:=function(n);
return Position(Tup[Length(n)],n);
end;
###############################################################

################################################################
Dim:=function(n);
if n>s then return fail; else
return d^n; fi;
end;
###############################################################

###############################################################
Boundary:=function(n,j)
#Boundary returns the image of the j-th generator of C(n)

local a,b,x,p,q,R,m,Q;

if n>s then
	return fail;
fi;
if j>d^n then
	return fail;
fi;

bound:=List([1..d^(n-1)], i->0);
Q:=ITT(j,n);

for a in [1..n-1] do
for b in [a+1..n] do
	p:=Length(Sctab[Q[a]][Q[b]][1]);
	for m in [1..p] do
		R:=StructuralCopy(Q);
		Add(R,Sctab[Q[a]][Q[b]][1][m],a);
		Remove(R,a+1);
		Remove(R,b);
		bound[TTI(R)]:=bound[TTI(R)]+((-1)^b)*Sctab[Q[a]][Q[b]][2][m];		
	od;
od;
od;

if IsPrimeInt(Charact(A)) then 
return bound*ONE;
else
return bound; fi;

end;
#############################################################

if not sparse then
return Objectify(HapChainComplex,rec(dimension:=Dim,
           boundary:=Boundary,
           properties:=
                [["length",s],
                 ["reduced",true],
                 ["type","chainComplex"],
                 ["characteristic",Charact(A)]]));
fi;


###############################################################
Boundary:=function(n,j)
#Boundary returns the image of the j-th generator of C(n)

local a,b,x,p,q,R,m,Q,pos;

if n>s then
        return fail;
fi;
if j>d^n then
        return fail;
fi;

bound:=[];
Q:=ITT(j,n);

for a in [1..n-1] do
for b in [a+1..n] do
        p:=Length(Sctab[Q[a]][Q[b]][1]);
        for m in [1..p] do
                R:=StructuralCopy(Q);
                Add(R,Sctab[Q[a]][Q[b]][1][m],a);
                Remove(R,a+1);
                Remove(R,b);
                #bound[TTI(R)]:=bound[TTI(R)]+((-1)^b)*Sctab[Q[a]][Q[b]][2][m];
                pos:=Position(bound,y->y[1]=TTI(R));
                if IsInt(pos) then 
                bound[pos][2]:=bound[pos][2] + ((-1)^b)*Sctab[Q[a]][Q[b]][2][m];
                else
		Add(bound,[TTI(R),((-1)^b)*Sctab[Q[a]][Q[b]][2][m] ]  );
		fi;
        od;
od;
od;

if IsPrimeInt(Charact(A)) then
Apply(bound, y->[y[1],y[2]*ONE]);
return bound;
else
return bound; fi;

end;
#############################################################

return Objectify(HapSparseChainComplex,rec(dimension:=Dim,
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
LeibnizMap:=function(A,s)

local Sour,Ran,DimSour,DimRan,Charact,Map,Tup,TupSour,ITT,TTI,j,n,r,BasisSour,a,BasisRan,Applicat,SourComp,RanComp,Chara,ONE,Arr;

Sour:=Source(A);
Ran:=Range(A);
BasisSour:=Basis(Sour);
BasisRan:=Basis(Ran);
DimSour:=Length(BasisSour);
DimRan:=Length(BasisRan);
TupSour:=List([1..s],r->Tuples([1..DimSour],r));
Tup:=List([1..s],r->Tuples([1..DimRan],r));
SourComp:=Leibnizcomplex(Sour,s);
RanComp:=Leibnizcomplex(Ran,s);
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
return StructuralCopy(TupSour[n][j]);
end;
###############################################################

###############################################################
TTI:=function(n);
return Position(Tup[Length(n)],n);
end;
###############################################################

#############################################################
Map:=function(n,i)

local Mapping,j,k,coef,count,perm;

if n>s then
	return fail;
fi;

Mapping:=List([1..DimRan^n], i->0);
for j in Tup[n] do
	coef:=1;
	count:=0;
	for k in j do
		count:=count+1;
		a:=BasisSour[ITT(i,n)[count]];
		coef:=Coefficients(BasisRan,Image(A,a))[k]*coef;
	od;
	Mapping[TTI(j)]:=Mapping[TTI(j)]+coef;
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

if IsLeftModule(L)=true then
	return Leibnizcomplex(L,x);
else if IsLeftModuleHomomorphism(L) then
	return LeibnizMap(L,x);
else
return fail;
fi;fi;
end);

#####################################################################
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(LeibnizAlgebraHomology,
function(A,n);

if not IsLeftModule(A)=true then
if not IsLeftModuleHomomorphism(A)=true then
	return fail;
fi;fi;
return Homology(ContractedComplex(LeibnizComplex(A,n+1)),n);

end);
#####################################################################
#####################################################################
