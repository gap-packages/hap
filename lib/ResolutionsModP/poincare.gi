#(C) Graham Ellis, 2005-2006

HapConstantPolRing:=PolynomialRing(Rationals,1);
#Sloppy!!! Should tidy this up.

#####################################################################
#####################################################################
InstallGlobalFunction(PoincareSeries,
function(arg)
local
	G,dim,
	R,Rdim,
	L, 
	TrialPoly,
	Dcoeffs,
	Ncoeffs,
	PolRing, x,
	k,bool,A,AA,i;

PolRing:=HapConstantPolRing;
x:=GeneratorsOfAlgebra(PolRing)[2];

if Length(arg)=1 then
if IsGroup(arg[1]) then
if IsPrimePowerInt(Order(arg[1])) then
return PoincareSeriesApproximation(arg[1],x);
else
Print("The group is not a p-group.\n");
return fail;
fi;
fi;
if IsAlgebra(arg[1]) then A:=arg[1]; fi;
if not IsBound(A!.degree) then Print("The algebra does not seem to be graded.\n"); return fail; fi;
A:=List(Basis(A),a->A!.degree(a));
AA:=[];
for i in [0..Maximum(A)] do
AA[i+1]:=Size(Filtered(A,a->a=i));
od;
return PoincareSeries(AA,Length(AA));
fi;

bool:=true;
dim:=arg[2];

if IsGroup(arg[1]) then
if IsPrimePowerInt(Order(arg[1])) then
G:=arg[1];
#R:=ResolutionPrimePowerGroup(G,dim);
Rdim:=RankPrimeHomology(G,dim);
bool:=false;
fi;
fi;

if IsHapResolution(arg[1]) then
R:=arg[1];
Rdim:=R!.dimension;
G:=R!.group;
if IsPrimeInt(EvaluateProperty(R,"characteristic")) then
bool:=false;fi;
fi;

L:=[];
if IsList(arg[1]) then
L:=arg[1];
bool:=false;
fi;

if bool then return fail; fi;

if Length(L)=0 then
L:=List([0..dim],i->Rdim(i));
fi;

####################################################################
TrialPoly:=function(L,k)
local M,i;

M:=[];
for i in [1..dim-k] do
Add(M,L{[i..i+k]});
od;
M:=TransposedMat(M);

return NullspaceMat(M);
end;
####################################################################

Dcoeffs:=[];
for k in [1..dim] do
Dcoeffs:=TrialPoly(L,dim-k);
if Length(Dcoeffs)=1 then break; fi;
od;

if not Length(Dcoeffs)=1 then 
#Print("Poincare series was not computed. \n");
return fail; fi;

Dcoeffs:=Reversed(Dcoeffs[1]);
Ncoeffs:=[];
for k in [1..Length(Dcoeffs)-1] do
Ncoeffs[k]:=Reversed(Dcoeffs{[1..k]})*L{[1..k]};
od;

return ValuePol(Ncoeffs,x)/ValuePol(Dcoeffs,x);
end);
#####################################################################
#####################################################################



#####################################################################
#####################################################################
InstallGlobalFunction(PoincareSeriesApproximation,
function(arg)
local
	G,kk,H,Hilbert,factos,x;

G:=arg[1];x:=arg[2];
factos:=DirectFactorsOfGroup(G);
if Length(factos)>1 then ; 
factos:=List(factos,g-> PoincareSeries(g));
return Product(factos);
else

return RankPrimeHomology(G,-1);

###########################################################
#Hilbert:=function(G,n,k)
#local P,R,i,B,PolRing,H;

#R:=ResolutionPrimePowerGroup(G,n);
#P:=[];
#for i in [1..k] do
#P[i]:=PoincareSeries(R,n+1-i);
#od;

#B:=List([1..k-1],i-> P[i] = P[i+1]);

#if not false in B then return P[1]; 
#else return fail; fi;
#end;
##########################################################

	
#kk:=10;
#H:=Hilbert(G,kk,2);
#while H=fail do
#kk:=kk+3;H:=Hilbert(G,kk,3);
#od;

#return H;
fi;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(ExpansionOfRationalFunction,
function(p,deg)
local
 	k, 
	Ncoeffs, 
	Dcoeffs, 
	Expansion;

if not (IsRationalFunction(p) and IsInt(deg)) then
Print("The input must be a rational function follwed by a positive integer. \n");
return fail;
fi;

Ncoeffs:=NumeratorOfRationalFunction(p);
Ncoeffs:=CoefficientsOfUnivariatePolynomial(Ncoeffs);
Ncoeffs:=ShallowCopy(Ncoeffs);

Dcoeffs:=DenominatorOfRationalFunction(p);
Dcoeffs:=CoefficientsOfUnivariatePolynomial(Dcoeffs);
Dcoeffs:=ShallowCopy(Dcoeffs);

#for k in [1..deg-Length(Ncoeffs)+1] do
for k in [1..deg+2] do
Append(Ncoeffs,[0]);
Append(Dcoeffs,[0]);
od;

Expansion:=[];
Expansion[1]:=Ncoeffs[1]/Dcoeffs[1];

for k in [2..deg+1] do
Expansion[k]:= 
(Ncoeffs[k] - Dcoeffs{[2..k]}*Reversed(Expansion{[1..k-1]}))/Dcoeffs[1];
od;
 
return Expansion;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
if true then
InstallGlobalFunction(PoincareSeriesPrimePart,
function(G,p,n)
local
	P,R,L,F;

P:=SylowSubgroup(G,p);
#R:=ResolutionFiniteGroup(P,n+1,false,p);
R:=ResolutionPrimePowerGroup(P,n+1);

####################################################################
F:=function(R);
return TensorWithIntegersModP(R,p);
end;
####################################################################

L:=List([1..n],i->Length(PrimePartDerivedFunctor(G,R,F,i)));
L:=Reversed(L);
Add(L,1);
L:=Reversed(L);

return(PoincareSeries(L,n));

end);
fi;
#####################################################################
#####################################################################
