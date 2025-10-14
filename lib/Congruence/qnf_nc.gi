DeclareFilter("QuadraticNF_NC");
DeclareProperty("IsQuadraticNumberField_NC",QuadraticNF_NC);
DeclareGlobalFunction("QuadraticNumberField_NC");
DeclareFilter("RingOfQuadraticIntegers_NC");
DeclareProperty("IsRingOfQuadraticIntegers_NC",RingOfQuadraticIntegers_NC);
DeclareAttribute("AssociatedQuadraticNumberField",IsQuadraticNumberField_NC);
DeclareFilter("IdealOfQuadraticIntegers_NC");
DeclareProperty("IsIdealOfQuadraticIntegers_NC",IdealOfQuadraticIntegers_NC);


##########################################################
##########################################################
InstallGlobalFunction(QuadraticNumberField_NC,
function(d)
local F;
F:=Field(HAPSqrt(d));
F!.Chatacteristic:=0;
F!.bianchiInteger:=d;
Setter(IsQuadraticNumberField_NC)(F,true);
SetName(F,Concatenation("Q(Sqrt(", String(d), "))"  ));
return F;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(RingOfIntegers,
"Ring of integers of quadratic number fields (non-cyclotomic)",
[IsQuadraticNumberField_NC],
function(F)
local d,D,R;
d:=F!.bianchiInteger;

if d mod 4 =1 then D:=(1+HAPSqrt(d))/2; fi;
if d mod 4 =2 or d mod 4 = 3 then D:=HAPSqrt(d); fi;
R:=Ring(D);
R!.bianchiInteger:=d;
Setter(IsRingOfQuadraticIntegers_NC)(R,true);
Setter(AssociatedQuadraticNumberField)(R,F);
SetName(R,Concatenation("O(",Name(F),")"));

return R;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Norm,
"Norm of an element in a ring of quadratic integers",
[IsQuadraticNumberField_NC,IsHapQuadraticNumber],
1000000,
function(R,x)
if not R!.bianchiInteger=x!.bianchiInteger then return fail; fi;
return x!.rational^2 -x!.bianchiInteger*x!.irrational^2;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Norm,
"Norm of an element in a ring of quadratic integers",
[IsRingOfQuadraticIntegers_NC,IsHapQuadraticNumber],
function(R,x)
return Norm(AssociatedQuadraticNumberField(R),x);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Trace,
"Trace of an element in a field of quadratic integers",
[IsQuadraticNumberField_NC,IsHapQuadraticNumber],
1000000,
function(R,x)
if not R!.bianchiInteger=x!.bianchiInteger then return fail; fi;
return 2*x!.rational;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Trace,
"Trace of an element in a ring of quadratic integers",
[IsRingOfQuadraticIntegers_NC,IsHapQuadraticNumber],
function(R,x)
if not R!.bianchiInteger=x!.bianchiInteger then return fail; fi;
return Trace(AssociatedQuadraticNumberField(R),x);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(IN,
"for ring of quadratic integers",
[IsHapQuadraticNumber,IsRingOfQuadraticIntegers_NC and IsRing],
function(x,R)
local p;

if not R!.bianchiInteger=x!.bianchiInteger then return false; fi;

if R!.bianchiInteger mod 4 = 1 then
if IsInt(2*x!.rational) and IsInt(x!.rational-x!.irrational) then return true; 
fi;
fi;

if IsInt(x!.rational) and IsInt(x!.irrational) then return true; fi;

return false;

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(QuadraticIdeal,
"ideal in a ring of quadratic integers",
[IsRing and IsRingOfQuadraticIntegers_NC, IsHapQuadraticNumber],
function(R,x);
return QuadraticIdeal(R,[x]);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(QuadraticIdeal,
"ideal in a ring of quadratic integers",
[IsRing and IsRingOfQuadraticIntegers_NC, IsList],
function(R,X)
local I, D, N, g, p, q, G, A, gD;

I:=TwoSidedIdeal(R,X);

Setter(AssociatedRing)(I,R);
Setter(IsIdealOfQuadraticIntegers_NC)(I,true);
Setter(IsRingOfQuadraticIntegers_NC)(I,false);

D:=GeneratorsOfRing(AssociatedRing(I));
D:=D[1];

G:=GeneratorsOfTwoSidedIdeal(I);
A:=[];
if R!.bianchiInteger mod 4 <> 1 then
for g in G do
p:=[g!.rational,g!.irrational];
gD:=g*D;
q:=[gD!.rational,gD!.irrational];;
Add(A,p); Add(A,q);
od;

else

for g in G do
p:=[g!.rational-g!.irrational,2*g!.irrational];
gD:=g*D;
q:=[gD!.rational-gD!.irrational,2*gD!.irrational];;
Add(A,p); Add(A,q);
od;


fi;

N:=HermiteNormalFormIntegerMat(A); #rows are new basis of I
Setter(NormOfIdeal)(I,N[1][1]*N[2][2]);

I!.hermiteBasis:=[N,D];
SetName(I,Concatenation("ideal of norm ", String(Norm(I))," in ",Name(R)) );
return I;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Norm,
"Norm of an  ideal in a ring of quadratic integers",
[IsIdealOfQuadraticIntegers_NC],
function(I)
local m;
return I!.NormOfIdeal;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(MOD,
"for a quadratic number  and an ideal in a ring of quadratic integers",
[IsHapQuadraticNumber,IsIdealOfQuadraticIntegers_NC and IsRing],
function(x,I)
local g, N, D, p;


N:=I!.hermiteBasis[1];
D:=I!.hermiteBasis[2];

#p:=PartsOfQuadraticInteger(AssociatedRing(I),x);
if not x!.bianchiInteger=AssociatedRing(I)!.bianchiInteger then return fail; fi;
if x!.bianchiInteger mod 4=1 then
p:=[x!.rational-x!.irrational,2*x!.irrational];
else
p:=[x!.rational,x!.irrational];
fi;

p:=p-Hap_int(p[1]/N[1][1])*N[1];
p:=p-Hap_int(p[2]/N[2][2])*N[2];

return p[1]+p[2]*D;

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(IN,
"for ideal in a ring of quadratic integers",
[IsHapQuadraticNumber,IsIdealOfQuadraticIntegers_NC and IsRing ],
function(x,I)
local g;

if x mod I = 0 then return true; else return false; fi;

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(RightTransversal,
"Coset representatives of an ideal in a ring of quadratic integers",
[IsIdealOfQuadraticIntegers_NC],
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
InstallOtherMethod(Discriminant,
"Discriminant of a quadratic number field",
[IsQuadraticNumberField_NC],
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
[IsRingOfQuadraticIntegers_NC],
function(R);
return Discriminant(AssociatedQuadraticNumberField(R));
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Units,
"Units of a ring of quadratic integers",
[IsRingOfQuadraticIntegers_NC],
1000000,#Hmm!
function(R) local d;
d:=R!.bianchiInteger;
if d>=0 then
  Print("Not implemented for d>0.\n");
  return fail;
fi;


if d=-1 then return [ [1,1]*One(HAPSqrt(d)) , [HAPSqrt(-1),-HAPSqrt(-1)] ]; fi;
if d=-3 then return [ [1,1]*One(HAPSqrt(d)) , [(1+HAPSqrt(-3))/2,(1-HAPSqrt(-3))/2] ]; fi;
return [[1,1]*One(HAPSqrt(d))];

end);
##########################################################
##########################################################

#################################
InstallOtherMethod(IsPrime,
"for deals in a ring of quadratic integers",
[IsIdealOfQuadraticIntegers_NC],
1000000,#Again, there must be a better way
function(I)
local N, p;
N:=Norm(I);
if IsPrimeInt(N) then return true;fi;
if N=1 then return false; fi;
N:=Factors(N);
if Length(N)>2 then return false; fi;
if Length(N)=1 then return true; fi;
p:=N[1]; #So Norm(I)=p^2.
if QuadraticCharacter(AssociatedRing(I),p)=-1 then return true;
else return false; fi;
end);
##################################

##########################################################
##########################################################
InstallOtherMethod(InverseOp,
"Inverse of an element mod an ideal",
[IsIdealOfQuadraticIntegers_NC,IsHapQuadraticNumber],  #NEED TO IMPROVE THIS!!!
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

