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
	L, Ex, 
	TrialPoly,
	Dcoeffs,
	Ncoeffs,
	PolRing, krull,
        D,x,f,k,bool,A,AA,i;
####################
#Inputs one of
#      : group G of prime-power order
#      : graded algebra A
#      : minimal resolution R of prime-power group and integer n>0
#      : group G of prime-power order and integer n>0
#      : list L of integers equal to the homology ranks in low degrees and an integer n>0 
#returns a rational function f(x)=N(x)/D(x) whose series is guaranteed to 
#agree with the cohomology ring Poincares series in degrees < dim. 
####################

PolRing:=HapConstantPolRing;
x:=GeneratorsOfAlgebra(PolRing)[2];

#####One Input Variable###########################################
if Length(arg)=1 then
   if IsGroup(arg[1]) then
      if IsPrimePowerInt(Order(arg[1])) then
         return PoincareSeriesApproximation(arg[1],2);
      else
         Print("The group is not a p-group.\n");
         return fail;
      fi;
   fi;
   if IsAlgebra(arg[1]) then A:=arg[1]; 
   fi;
   if not IsBound(A!.degree) then 
      Print("The algebra does not seem to be graded.\n"); return fail; 
   fi;
   A:=List(Basis(A),a->A!.degree(a));
   AA:=[];
   for i in [0..Maximum(A)] do
      AA[i+1]:=Size(Filtered(A,a->a=i));
   od;
   return PoincareSeries(AA,Length(AA)-1);
fi;
##################################################################

#####Two Input Variables###########################################
bool:=true;
dim:=arg[2];

if IsGroup(arg[1]) then
   if IsPrimePowerInt(Order(arg[1])) then
         return PoincareSeriesApproximation(arg[1],arg[2]);
   else
         Print("The group is not a p-group.\n");
         return fail;
   fi;

#      G:=arg[1];
#      Rdim:=RankPrimeHomology(G,dim);
#      krull:=Prank(G);
#      bool:=false;
fi;

if IsHapResolution(arg[1]) then
   R:=arg[1];
   Rdim:=R!.dimension;
   G:=R!.group;
   krull:=Prank(G);
   if IsPrimeInt(EvaluateProperty(R,"characteristic")) then
      bool:=false;
   fi;
fi;

L:=[];
if IsList(arg[1]) then
   L:=arg[1];
   bool:=false;
fi;

if bool then return fail; 
fi;

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
   if Length(Dcoeffs)=1 then break; 
   fi;
od;

if not Length(Dcoeffs)=1 then 
#   Print("Try inputting more degrees of a resolution/algebra/sequence. \n");
   return fail; 
fi;

Dcoeffs:=Reversed(Dcoeffs[1]);
Ncoeffs:=[];
for k in [1..Length(Dcoeffs)-1] do
   Ncoeffs[k]:=Reversed(Dcoeffs{[1..k]})*L{[1..k]};
od;

f:=ValuePol(Ncoeffs,x)/ValuePol(Dcoeffs,x);

Ex:=ExpansionOfRationalFunction(f,1000); #Check at least that the first 500 
                                         #terms are integers and the initial 
                                         #terms agree with the input data and
                                         #krull is the number of poles (x-1) 
                                         #in the denominator of f (in cases 
                                         #where krull is available.
if (not fail=PositionProperty(Ex,x->not IsInt(x)) ) 
   or
   ( not L{[1..dim]}=Ex{[1..dim]}) 
   then 
#      Print("Try inputting more degrees of a resolution/algebra/sequence. \n");
      return fail; 
fi;  
if IsBound(krull) then
   D:=Factors(DenominatorOfRationalFunction(f));
   D:=Filtered(D,x->x^2 = (IndeterminateOfUnivariateRationalFunction(f)-1)^2 );
   if not krull = Length(D) then
#      Print("Try inputting more degrees of a resolution/algebra/sequence. \n");
      return fail;
   fi;
fi;
#    Print("The series is guaranteed correct for group cohomology in degrees < ",dim,"\n");
return f;
end);
#####################################################################
#####################################################################



#####################################################################
#####################################################################
InstallGlobalFunction(PoincareSeriesApproximation,
function(arg)
local
	G,kk,H,Hilbert,factos,guarantee;

G:=arg[1];
factos:=DirectFactorsOfGroup(G);
if Length(factos)>1 then  
   if Length(arg)=1 then
      factos:=List(factos,g-> RankPrimeHomology(g,-1,5));
      guarantee:=Minimum(List(factos,x->x[2]));
      factos:=List(factos,x->x[1]);
   else
      factos:=List(factos,g-> RankPrimeHomology(g,-1,arg[2]));
      guarantee:=Minimum(List(factos,x->x[2]));
      factos:=List(factos,x->x[1]);
   fi;
Print("The series is guaranteed correct for group cohomology in degrees < ",guarantee,"\n");
return Product(factos);
else

if Length(arg)=1 then
factos:= RankPrimeHomology(G,-1);
else
factos:= RankPrimeHomology(G,-1,arg[2]);
fi;

Print("The series is guaranteed correct for group cohomology in degrees < ",factos[2],"\n");
return factos[1];

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

#L:=List([1..n],i->Length(PrimePartDerivedFunctor(G,R,F,i)));
L:=List([1..n],i->Length(PrimePartDerivedFunctorViaSubgroupChain(G,R,F,i)));
L:=Reversed(L);
Add(L,1);
L:=Reversed(L);

Print("The series is guaranteed correct for group cohomology in degrees < ",n+1,"\n");
return(PoincareSeries(L,n));

end);
fi;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(ModPCohomologyPresentationBounds,
function(G)
local D,irr,S,SS,i,g;

if IsCyclic(G) then 
if Order(G)=2 then return rec( generators_degree_bound:=1,
      relators_degree_bound:=0      );; fi;
return rec( generators_degree_bound:=1,
      relators_degree_bound:=2      );;;
fi;

D:=DirectFactorsOfGroup(G);

if Length(D)=1 then
################################################
irr:=IrreducibleRepresentations(G);

for i in [1..Length(irr)] do
   for S in Combinations(irr,i) do
   if Order(Intersection(List(S,r->Kernel(r)))) = 1 then
      SS:=List(S,r->( Length( One( Range(r) ) ) )^2  );
      return rec( generators_degree_bound:=Sum(SS),
      relators_degree_bound:=2*Sum(SS)      );
   fi;
   od;
od;
########################################
fi;

SS:=[];
for g in D do
Add(SS,(ModPCohomologyPresentationBounds(g)).generators_degree_bound);
od;

return rec( generators_degree_bound:=Maximum(SS),
      relators_degree_bound:=2*Maximum(SS)      );;
end);
#####################################################################
#####################################################################
