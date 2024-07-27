##########################################################
##########################################################
InstallGlobalFunction(QuadraticNumberField,
function(d) 
local F;
F:=Field(Sqrt(d));
F!.bianchiInteger:=d;
Setter(IsQuadraticNumberField)(F,true);
SetName(F,Concatenation("Q(Sqrt(", String(d), "))"  ));
return F;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(RingOfIntegers,
"Ring of integers of quadratic number fields",
[IsNumberField],
function(F)
local d,D,R;
if not IsBound(F!.bianchiInteger) then TryNextMethod(); fi;
d:=F!.bianchiInteger;

if d mod 4 =1 then D:=(1+Sqrt(d))/2; fi;
if d mod 4 =2 or d mod 4 = 3 then D:=Sqrt(d); fi;
R:=Ring(D);
R!.bianchiInteger:=d;
Setter(IsRingOfQuadraticIntegers)(R,true);
Setter(AssociatedNumberField)(R,F);
SetName(R,Concatenation("O(",Name(F),")"));

return R;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Norm,
"Norm of an element in a ring of quadratic integers",
[IsRingOfQuadraticIntegers,IsCyclotomic],
function(R,x)
return Norm(AssociatedNumberField(R),x);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Trace,
"Trace of an element in a ring of quadratic integers",
[IsRingOfQuadraticIntegers,IsObject],
function(R,x)
return Trace(AssociatedNumberField(R),x);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(IN,
"for ring of quadratic integers",
[IsCyclotomic,IsRingOfQuadraticIntegers and IsRing],
function(x,R)
local p; 

p:=PartsOfQuadraticInteger(R,x);
if p=fail then return false; fi;

if IsInt(p[1]) and IsInt(p[2]) then return true;
else return false;
fi;

  
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(QuadraticIdeal,
"ideal in a ring of quadratic integers",
[IsRing and IsRingOfQuadraticIntegers, IsCyclotomic],
function(R,x);
return QuadraticIdeal(R,[x]);
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallOtherMethod(QuadraticIdeal,
"ideal in a ring of quadratic integers",
[IsRing and IsRingOfQuadraticIntegers, IsList],
function(R,X)
local I, D, N, g, p, q, G, A;

if IsCyclotomic(X) then
I:=TwoSidedIdeal(R,[X]);
else
I:=TwoSidedIdeal(R,X);
fi;
Setter(AssociatedRing)(I,R);
Setter(IsIdealOfQuadraticIntegers)(I,true);
Setter(IsRingOfQuadraticIntegers)(I,false);

D:=GeneratorsOfRing(AssociatedRing(I));
D:=D[1];

G:=GeneratorsOfTwoSidedIdeal(I);
A:=[];
for g in G do
p:=PartsOfQuadraticInteger(R,g);
q:=PartsOfQuadraticInteger(R,g*D);
Add(A,p); Add(A,q);
od;

N:=HermiteNormalFormIntegerMat(A); #rows are new basis of I
Setter(NormOfIdeal)(I,N[1][1]*N[2][2]);
#if R!.bianchiInteger mod 4 = 1 then D:=2*(D-(1/2)); fi;
I!.hermiteBasis:=[N,D];
SetName(I,Concatenation("ideal of norm ", String(Norm(I))," in ",Name(R)) );
return I;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(IN,
"for ideal in a ring of quadratic integers",
[IsCyclotomic,IsIdealOfQuadraticIntegers and IsRing ],
function(x,I)
local g;

if x mod I = 0 then return true; else return false; fi;

##################This code is never used#########
#g:=GeneratorsOfTwoSidedIdeal(I)[1]; 
#if not 
#IsInt(  Norm(AssociatedRing(I),x) / Norm(I) )
#then return false; fi;
#return x*g^-1 in AssociatedRing(I);
##################################################
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(MOD,
"for a cyclotomic  and an ideal in a ring of quadratic integers",
[IsCyclotomic,IsIdealOfQuadraticIntegers and IsRing],
function(x,I)   
local g, N, D, p, int;

############This code is never used################
#for g in RightTransversal(I) do
#if g-x in I then return g; fi;
#od;
###################################################

##############################
int:=function(y);
if y>0 or IsInt(y) then return Int(y); 
else return Int(y)-1; fi;
end;
##############################

N:=I!.hermiteBasis[1];
D:=I!.hermiteBasis[2];

p:=PartsOfQuadraticInteger(AssociatedRing(I),x);
p:=p-int(p[1]/N[1][1])*N[1];
p:=p-int(p[2]/N[2][2])*N[2];

return p[1]+p[2]*D;

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(InverseOp,
"Inverse of an element mod an ideal",
[IsIdealOfQuadraticIntegers,IsCyclotomic],  #NEED TO IMPROVE THIS!!!
function(I,x)
local T,y;

T:=RightTransversal(I);
for y in T do
if x*y-1 in I then return y; fi;
od;
return fail;
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(Norm,
"Norm of an  ideal in a ring of quadratic integers",
[IsIdealOfQuadraticIntegers],
function(I)
local m;
return I!.NormOfIdeal;
end);
##########################################################
##########################################################


DeclareOperation("RightTransversal_alt",[IsIdealOfQuadraticIntegers]);
##########################################################
##########################################################
InstallOtherMethod(RightTransversal,
"Coset representatives of an ideal in a ring of quadratic integers",
[IsIdealOfQuadraticIntegers],
function(I)
local cosetreps,R,d,D,N,i,j;

if IsBound(I!.RightTransversal) then return I!.RightTransversal; fi;
R:=AssociatedRing(I);
d:=R!.bianchiInteger;
D:=I!.hermiteBasis[2];;
N:=I!.hermiteBasis[1];

cosetreps:=[];
for i in [0..N[1][1]-1] do
for j in [0..N[2][2]-1] do
Add(cosetreps,i+j*D);
od;od;

I!.RightTransversal:=cosetreps;
return cosetreps;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(RightTransversal_alt,
"For test purposes only: Coset representatives of an ideal in a ring of quadratic integers",
[IsIdealOfQuadraticIntegers],
function(I)
local cosetreps,leaves,newleaves,isinsetModI,norm,R,A,H,P,g,a,b,x,y,i,j,D;

if IsBound(I!.RightTransversal) then return I!.RightTransversal; fi;

norm:=Norm(I);
if IsPrimeInt(norm) then
#   I!.RightTransversal:=[0..norm-1];
#   return I!.RightTransversal;
fi;

D:=GeneratorsOfRing(AssociatedRing(I));
D:=D[1];
g:=GeneratorsOfTwoSidedIdeal(I);
g:=g[1];
R:=AssociatedRing(I);

cosetreps:=[0];
leaves:=[0];

#################################
isinsetModI:=function(S,x)
local s;
for s in S do
if (s-x) in I then return true; fi;
od;
return false;

end;
#################################

while Size(leaves)>0 do
newleaves:=[];
for x in leaves do
y:=x+1; 
if not isinsetModI(cosetreps,y) then
Add(cosetreps,y);  
Add(newleaves,y);
fi;
y:=x+D;
if not isinsetModI(cosetreps,y) then
Add(cosetreps,y); 
Add(newleaves,y);
fi;
od;
leaves:=newleaves;
od;

I!.RightTransversal:=cosetreps;
return cosetreps;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Discriminant,
"Discriminant of a quadratic number field",
[IsNumberField],
function(Q);
if Q!.bianchiInteger mod 4 = 1 then return Q!.bianchiInteger;
else return 4*Q!.bianchiInteger; fi;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Discriminant,
"Discriminant of the quadratic number field associated to a ring of integers",
[IsRingOfQuadraticIntegers],
function(R);
return Discriminant(AssociatedNumberField(R));
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Units,
"Units of a ring of quadratic integers",
[IsRingOfQuadraticIntegers],
1000000,#Hmm!
function(R) local d;
d:=R!.bianchiInteger;
if d>=0 then 
  Print("Not implemented for d>0.\n");
  return fail; 
fi;

if d=-1 then return [ [1,1] , [Sqrt(-1),-Sqrt(-1)] ]; fi;
if d=-3 then return [ [1,1] , [(1+Sqrt(-3))/2,(1-Sqrt(-3))/2] ]; fi;
return [[1,1]];

end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallGlobalFunction(PartsOfQuadraticInteger,
function(R,x)
local q,a,b,d,r;
#This will aslo work for a non-integer.

d:=R!.bianchiInteger;

#a:=RealPart(x);
#b:=ImaginaryPart(x)/Sqrt(AbsInt(d));
#if d mod 4 = 1 then a:=a-b; b:=2*b;fi;
#return [a,b];

#################IGNORE EVERYTHING ABOVE 

if IsCycInt(x) then
q:=Quadratic(x);

a:=q.a/q.d;
b:=q.b/q.d;


if d mod 4 = 1 then a:=a-b; b:=2*b;fi;
return [a,b];
fi;

r:=Sqrt(d);
a:=Trace(R,x)/2;
b:=(x-a)/r;

if not d mod 4 = 1 then
  return [a,b];
else
  return [a-b,2*b];
fi;


end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(SL2QuadraticIntegers,
function(arg)
local X,R,S,d,gens,G;

X:=arg[1];

if IsInt(X) then
   R:=ResolutionSL2QuadraticIntegers(X,2);
   gens:=Generators(R,true);
   G:=Group(gens);
   SetName(G,Concatenation("SL(2,O(Q(",String(X),"))"));
   return G;
fi;
if IsRing(X) then
  if not IsBound(X!.associatedIdeal) then return fail; 
  else S:=X!.associatedIdeal;  S:=S!.AssociatedRing;
  fi;
  if not IsBound(S!.bianchiInteger) then return fail; 
  else d:=S!.bianchiInteger; fi;

   R:=ResolutionSL2QuadraticIntegers(d,2);
   gens:=Generators(R,true);
   G:=Group(gens*One(S));
   SetName(G,Concatenation("SL(2,", Name(X),")"));
   return G;
  
fi;
end);
##########################################################
##########################################################




