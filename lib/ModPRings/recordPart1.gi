#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(ModPCohomologyRing_part_1,
function(arg)
local
	G,R,lngth,A,arg3,Degree,
	GeneratorLIFTS,
	Lift,
	prime,
	GactVec,
	MT,
	eltsG,
	pp,
	one,zero,
	IntPairList,
	IntToPair, IntToPairModified,
	PairToInt, PairToIntModified,
	ComposeGens,
	ComposeGensAugmented,
	SCT,
	dim,
	SolutionMatBoundaryMatrices,
	RngGens,
	NiceBas,
	bastmp,
	gensSupport,
	LstWIE,LstWIEWW,
	EltMats,
	M,n,i,j,k,g,mm;

if IsHapResolution(arg[1]) then R:=arg[1]; 
G:=R!.group; 
else
G:=arg[1];
R:=ResolutionPrimePowerGroup(G,arg[2]); fi;

pp:=Order(G);

if not 
	(
	IsPrimePowerInt(pp)
	and
	EvaluateProperty(R,"isMinimal")=true 
	and
	EvaluateProperty(R,"characteristic")>0
	)
	then 
Print("The resolution must be minimal, for a prime-power group, and over a finite field. \n");
return fail;
fi;

if (IsHapResolution(arg[1]) and Length(arg)=1)
or (IsGroup(arg[1]) and Length(arg)=2) then
arg3:="medium";
else
arg3:=arg[Length(arg)];
fi;

prime:=Factors(pp)[1];
one:=Identity(GF(prime));
zero:=0*one;


eltsG:=Elements(G);
MT:=MultiplicationTable(eltsG);

IntPairList:=[];
for n in [1..Length(R)] do
for i in [1..R!.dimension(n)] do
Append(IntPairList,[[n,i]]);
od;
od;
dim:=Length(IntPairList);


#####################################################################
#GactVec:=function(g,v) 
#local k,q,h,C;

#C:=ListWithIdenticalEntries(Length(v),one);

#k:=0;
#for q in [0..(-1+Length(v)/pp)] do

#for h in [1..pp] do
#C[k+MT[g][h]]:=v[k+h]; 
#od;
#k:=k+pp;
#od;

#ConvertToVectorRepNC(C,prime);

#return C;

#end;
#####################################################################

####################################################################
EltMats:=[];
for g in [1..pp] do
M:=R!.actMat(g,IdentityMat(Size(G))*one);
M:=TransposedMat(M);
ConvertToMatrixRepNC(M);
Add(EltMats,M);
od;
####################################################################

#####################################################################
GactVec:=function(g,v)
local q,C;
C:=[];
 
return
Concatenation(List(
[0..-1+Length(v)/pp],q->v{[q*pp+1..(q+1)*pp]}*EltMats[g]));

end;
####################################################################

#####################################################################
IntToPair:=function(k);
return IntPairList[k]; 
end;
#####################################################################

#####################################################################
IntToPairModified:=function(k);
if k = 1 then return [0,1]; 
else           
return IntPairList[k-1]; 
fi;
end;
#####################################################################

#####################################################################
PairToInt:=function(x);
return Position(IntPairList,x);
end;
#####################################################################

#####################################################################
PairToIntModified:=function(x);
if x = [0,1] then  return 1;       
else                      
return Position(IntPairList,x)+1;  
fi;
end;
#####################################################################



#####################################################################
SolutionMatBoundaryMatrices:=R!.solutionMatBoundaryMatrices;
#####################################################################

GeneratorLIFTS:=[1..Length(R)];
GeneratorLIFTS:=List(GeneratorLIFTS,n->[1..Dimension(R)(n)]);

	#The cohomology generators are represented by pairs [n,i] where n
	#lies between 1 and Length(R), and i lies between 1 and 
	#Dimension(R)(n).
	
	#GeneratorLifts[n,i] is a list [F[0],F[1],...,F[Length(R)-n]
	#where F[m] represents a ZG-homomorphism R_{n+m}-->R_{m}. Here F[m] is
	#actually a list [w_1,...,w_d] of words in R_m, with w_j the
	#image of the j-th generator of R_{n+m}. The ZG-homomorphism is 
	#induced by the standard cocycle representing the generator [n,i].

LstWIE:=ListWithIdenticalEntries(pp,zero);
LstWIEWW:=[];
for i in [1..Length(R)] do
LstWIEWW[i]:=ListWithIdenticalEntries(pp*Dimension(R)(i-1),zero);
ConvertToVectorRep(LstWIEWW[i],prime);
od;

#####################################################################
Lift:=function(arg) 	#This function calculates ALL GeneratorLifts[n,i]
local			#OR just lifts up to deg=arg[3].
	n,i,
	FM,j,W,m,S,x,deg,t;

n:=arg[1];
i:=arg[2];
if Length(arg)>2 then deg:=arg[3]; 
else
deg:=Minimum(Length(R)-n,Int(Length(R)/2),n);
fi;

if  IsInt(GeneratorLIFTS[n][i]) then   #
GeneratorLIFTS[n][i]:=[];
fi;    					    #

if not IsBound(GeneratorLIFTS[n][i][1]) then
FM:=[]; #m=0
for j in [1..Dimension(R)(n)] do
if j=i then 
FM[j]:=StructuralCopy(LstWIE); FM[j][1]:=one;
else 
FM[j]:=StructuralCopy(LstWIE);
fi;
ConvertToVectorRep(FM[j],prime);
od;
GeneratorLIFTS[n][i][1]:=FM;
fi;

for m in [1..deg] do
if not IsBound(GeneratorLIFTS[n][i][m+1]) then

FM:=[];

for j in [1..Dimension(R)(n+m)] do
W:=StructuralCopy(LstWIEWW[m]);

for x in R!.boundary(n+m,j) do
if prime=2 then
S:=GeneratorLIFTS[n][i][m][AbsInt(x[1])];
else
t:=SignInt(x[1])*one;
S:=t*GeneratorLIFTS[n][i][m][AbsInt(x[1])];
fi;
W:=W+GactVec(x[2],S);

od;

FM[j]:=SolutionMatBoundaryMatrices(m,W);

od;
GeneratorLIFTS[n][i][m+1]:=FM;


fi;
od;
end;
#####################################################################

if arg3="medium" or arg3="low" then
#####################################################################
ComposeGens:=function(x,y)
local xx,yy,prdct,prod,j,z;


if y>x then prod:=ComposeGens(y,x);
if IsOddInt(IntToPair(x)[1]*IntToPair(y)[1]) then  #graded commutativity!!
for j in [1..Length(prod)/2] do
prod[2*j-1]:=-prod[2*j-1];
od;
fi;
return prod;
fi;

xx:=IntToPair(x);
yy:=IntToPair(y);

if xx[1]+yy[1]>Length(R) then return []; fi;

prdct:=[];
for j in [1..R!.dimension(xx[1]+yy[1])] do
prdct[j]:=Sum(GeneratorLIFTS[xx[1]][xx[2]][yy[1]+1][j]{[(yy[2]-1)*pp+1..yy[2]*pp]});

od;

prdct:= prdct*one ;
prod:=[];

j:=PairToInt([xx[1]+yy[1],1])-1;
for z in [1..Length(prdct)] do
if not prdct[z]=zero then
Append(prod,[prdct[z],j+z]);
fi;
od;

return prod; 
end;
#####################################################################
fi;
if arg3="high" then
#####################################################################
ComposeGens:=function(x,y)
local xx,yy,prdct,prod,j,z;


#if y>x then prod:=ComposeGens(y,x);
#if IsOddInt(IntToPair(x)[1]*IntToPair(y)[1]) then  #graded commutativity!!
#for j in [1..Length(prod)/2] do
#prod[2*j-1]:=-prod[2*j-1];
#od;
#fi;
#return prod;
#fi;

xx:=IntToPair(x);
yy:=IntToPair(y);

if xx[1]+yy[1]>Length(R) then return []; fi;

prdct:=[];
for j in [1..R!.dimension(xx[1]+yy[1])] do
prdct[j]:=Sum(GeneratorLIFTS[xx[1]][xx[2]][yy[1]+1][j]{[(yy[2]-1)*pp+1..yy[2]*pp]});

od;

prdct:= prdct*one ;
prod:=[];

j:=PairToInt([xx[1]+yy[1],1])-1;
for z in [1..Length(prdct)] do
if not prdct[z]=zero then
Append(prod,[prdct[z],j+z]);
fi;
od;

return prod;
end;
#####################################################################
fi;
#####################################################################
ComposeGensAugmented:=function(i,j) #ComposeGens and all other above
				    #functions do not incorporate 
local prod,k,yy;		    #H^0(G,Z_p). ComposeGensAugmented
				    #does so, and thus involves a 
				    #clumsy dimension shift.
if i>1 and j>1 then
prod:=ComposeGens(i-1,j-1);
for k in [1..Length(prod)/2] do
prod[2*k]:=prod[2*k]+1;
od;
return prod;
fi;

if i=1 then return [1,j]; fi;
if j=1 then return [1,i]; fi;
end;
#####################################################################

#####################################################################
Degree:=function(x)
local i, w, bas;
# returns the highest degree of a non-zero coefficient of x

if IsZero(x) then return 0; fi;
i:=Position(GeneratorsOfAlgebra(A),x);

if i=1 then return 0; fi;

if not i=fail then return IntToPair(i-1)[1]; fi;

bas:=Basis(A);
w:=Coefficients(bas,x);
w:=Filtered([2..Length(bas)],i->not IsZero(w[i]));
w:=List(w,i->IntPairList[i-1][1]);
return Maximum(w);
end;
#####################################################################


SCT:=EmptySCTable(dim+1,Zero(GF(prime)));

#################################################
if arg3="high" then


for i in [1..R!.dimension(1)] do
Lift(1,i,Length(R)-1);
od;

for i in [1..Dimension(R)(1)+1] do
for j in [1..dim+1] do
SetEntrySCTable(SCT,i,j,ComposeGensAugmented(i,j));
od;
od;

A:=AlgebraByStructureConstants(GF(prime),SCT);
A!.degree:=Degree;
A!.intToPair:=IntToPair;
A!.intToPairModified:=IntToPairModified;
A!.pairToInt:=PairToInt;
A!.pairToIntModified:=PairToIntModified;
RngGens:=ModPRingGenerators(A);
if Length(RngGens)=Dimension(R)(1)+1 then
return A; fi;

RngGens:=RngGens{[Dimension(R)(1)+2..Length(RngGens)]};
RngGens:=List(RngGens,x->Position(Basis(A),x));



for k in RngGens do
n:=IntToPair(k-1);
Lift(n[1],n[2],Length(R)-n[1]); ###
od;

for i in RngGens do
for j in [1..dim+1] do
SetEntrySCTable(SCT,i,j,ComposeGensAugmented(i,j));
od;
od;

A:=AlgebraByStructureConstants(GF(prime),SCT);
A!.degree:=Degree;
A!.intToPair:=IntToPair;
A!.intToPairModified:=IntToPairModified;
A!.pairToInt:=PairToInt;
A!.pairToIntModified:=PairToIntModified;

return A;
fi;
#################################################

#################################################
if arg3="medium" then

for n in [1..Length(R)-1] do
for i in [1..R!.dimension(n)] do
Lift(n,i);
od;
od;

for i in [1..dim+1] do
for j in [1..dim+1] do
SetEntrySCTable(SCT,i,j,ComposeGensAugmented(i,j));
od;
od;

A:=AlgebraByStructureConstants(GF(prime),SCT);

A!.degree:=Degree;
A!.intToPair:=IntToPair;
A!.intToPairModified:=IntToPairModified;
A!.pairToInt:=PairToInt;
A!.pairToIntModified:=PairToIntModified;

return A;

fi;
#################################################

#################################################
if arg3="low" then

for n in [1..Length(R)-1] do
for i in [1..R!.dimension(n)] do
Lift(n,i,1);
od;
od;

for i in [1..Dimension(R)(1)+1] do
for j in [1..dim+1] do
SetEntrySCTable(SCT,i,j,ComposeGensAugmented(i,j));
SetEntrySCTable(SCT,j,i,ComposeGensAugmented(j,i));
od;
od;

A:=AlgebraByStructureConstants(GF(prime),SCT);
A!.degree:=Degree;
A!.intToPair:=IntToPair;
A!.intToPairModified:=IntToPairModified;
A!.pairToInt:=PairToInt;
A!.pairToIntModified:=PairToIntModified;
RngGens:=ModPRingGenerators(A);
if Length(RngGens)=Dimension(R)(1)+1 then
return A; fi;

RngGens:=RngGens{[Dimension(R)(1)+2..Length(RngGens)]};
bastmp:=Basis(Subalgebra(A,RngGens));
#bastmp:=Basis(Ideal(A,RngGens)); #There is a mathematical problem here.
				  #The resulting algebra is NOT associative
				  #when "Subalgebra" is used. But using "Ideal"
				  #slows the algorithm (though yields the right
				  #answer!
gensSupport:=[];		
for i in bastmp do
Append(gensSupport,
Filtered([1..Dimension(A)],x->not IsZero(Coefficients(Basis(A),i)[x]))
);
od;
gensSupport:=SSortedList(gensSupport);

mm:=Maximum(List(RngGens,x->Degree(x)));###

RngGens:=List(RngGens,x->Position(Basis(A),x));

#for k in RngGens do
for k in gensSupport do ###
n:=IntToPair(k-1);
#Lift(n[1],n[2]);
Lift(n[1],n[2],Minimum(Length(R)-n[1],Int(Length(R)/2),n[1],mm)); ### 
od;

#for i in RngGens do
#for j in RngGens do
for i in gensSupport do ###
for j in gensSupport do  ###
SetEntrySCTable(SCT,i,j,ComposeGensAugmented(i,j));
SetEntrySCTable(SCT,j,i,ComposeGensAugmented(j,i));
od;
od;

A:=AlgebraByStructureConstants(GF(prime),SCT);
A!.degree:=Degree;
A!.intToPair:=IntToPair;
A!.intToPairModified:=IntToPairModified;
A!.pairToInt:=PairToInt;
A!.pairToIntModified:=PairToIntModified;
A!.gensSupport:=gensSupport;

return A;
fi;
#################################################


end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(ModPRingGeneratorsAlt,
function(A)
local S, gens, gensA, x, dim, dimgensA;
#I need to think about this. The following procedure is inefficient.


S:=GeneratorsOfAlgebra(A);
dim:=Dimension(A);
gens:=[S[1]];
gensA:=Subalgebra(A,gens,"basis");
dimgensA:= Dimension(gensA);

for x in S do

if not x in gensA then
Append(gens,[x]);
gensA:=SubalgebraNC(A,gens);
dimgensA:=Dimension(gensA);

fi;

if dimgensA=dim then 
return gens; fi;

od;


end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(ModPRingGenerators,
function(A)
local S, gradedgens, deg, mingens, i,j,n, vecs,V,x,y;

if "mingens" in NamesOfComponents(A) then return A!.mingens;fi;


S:=GeneratorsOfAlgebra(A);
S:=Filtered(S,a->not IsZero(a));  #Added October 2023
deg:=A!.degree(S[Length(S)]);
gradedgens:=List([1..1+deg],i->[]);

for i in [1..Length(S)] do
Add(gradedgens[1+A!.degree(S[i])],S[i]);
od;

#gradedgens[i+1] contains generators of degree i.

mingens:=Concatenation(gradedgens[1+0],gradedgens[1+1]);

for n in [2..deg] do
vecs:=[];

   for i in [1..Int(n/2)] do
 	for x in gradedgens[1+i] do
	for y in gradedgens[n+1-i] do
	Add(vecs,x*y);
	od;
	od;
   od;
	V:=SubspaceNC(A,vecs);
	for x in gradedgens[n+1] do
	if not x in V then Add(mingens,x); 
	Add(vecs,x);
	V:=SubspaceNC(A,vecs);
	fi;
	od;

od;

A!.mingens:=mingens;

return mingens;
end);
#####################################################################
#####################################################################

