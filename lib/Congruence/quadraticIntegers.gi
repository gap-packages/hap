##########################################################
##########################################################
InstallGlobalFunction(QuadraticNumberField,
function(d) 
local F;
F:=Field(Sqrt(d));
F!.bianchiInteger:=d;
Setter(IsQuadraticNumberField)(F,true);
SetName(F,Concatenation("Rationals( Sqrt(", String(d), ") )"  ));
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
Setter(AssociatedField)(R,F);
SetName(R,Concatenation("ring of integers of ",Name(F)));

return R;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Norm,
"Norm of an element in a ring of quadratic integers",
[IsRingOfQuadraticIntegers,IsObject],
function(R,x)
return Norm(AssociatedField(R),x);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Trace,
"Trace of an element in a ring of quadratic integers",
[IsRingOfQuadraticIntegers,IsObject],
function(R,x)
return Trace(AssociatedField(R),x);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(IN,
"for ring of quadratic integers",
[IsObject,IsRingOfQuadraticIntegers and IsRing],
function(x,R)
local a,b;
a:=Norm(R,x);
b:=Trace(R,x);
return IsInt(a) and IsInt(b);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(PrincipalIdeal,
"principal ideal in a ring of quadratic integers",
[IsRing and IsRingOfQuadraticIntegers,IsObject],
function(R,x)
local I;
I:=TwoSidedIdeal(R,[x]);
Setter(AssociatedRing)(I,R);
Setter(IsPrincipalIdeal)(I,true);
SetName(I,Concatenation("principal ideal in ",Name(R)) );
return I;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(IN,
"for principal ideal in a ring of quadratic integers",
[IsObject,IsPrincipalIdeal and IsRing],
function(x,I)
local g;
g:=GeneratorsOfTwoSidedIdeal(I)[1];
return x*g^-1 in AssociatedRing(I);
end);
##########################################################
##########################################################



