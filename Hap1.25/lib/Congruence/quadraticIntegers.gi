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
local a,b; 
b:=Trace(R,x);
if not IsInt(b) then return false; fi;
a:=Norm(R,x);
if not IsInt(a) then return false; 
else return true; fi;;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(PrincipalIdeal,
"principal ideal in a ring of quadratic integers",
[IsRing and IsRingOfQuadraticIntegers,IsCyclotomic],
function(R,x)
local I;

I:=TwoSidedIdeal(R,[x]);
Setter(AssociatedRing)(I,R);
Setter(IsPrincipalIdeal)(I,true);
Setter(IsRingOfQuadraticIntegers)(I,false);
Setter(NormOfPrincipalIdealGenerator)(I,Norm(R,x));
SetName(I,Concatenation("principal ideal of norm ", String(Norm(I))," in ",Name(R)) );
return I;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(IN,
"for principal ideal in a ring of quadratic integers",
[IsCyclotomic,IsPrincipalIdeal and IsRing ],
function(x,I)
local g;

g:=GeneratorsOfTwoSidedIdeal(I)[1]; 
return x*g^-1 in AssociatedRing(I);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(MOD,
"for a cyclotomic  and a principal ideal in a ring of quadratic integers",
[IsCyclotomic,IsPrincipalIdeal and IsRing],
function(x,I)   #THERE IS A BETTER METHOD!!!
local n,i,g,norm ;


if x in I then return 0; fi;

if IsBound(I!.RightTransversal) then
for g in I!.RightTransversal do
if g-x in I then return g; fi;
od;
fi;


norm:=NormOfPrincipalIdealGenerator(I);
g:=GeneratorsOfRing(AssociatedRing(I));
g:=g[1];
n:=0;
while true do
n:=n+1;
for i in [0..Maximum(n,norm)] do
if (x - i - (n-i)*g) in I then 
return i + (n-i)*g;
fi;
od;
od;

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(InverseOp,
"Inverse of an element mod an ideal",
[IsPrincipalIdeal,IsCyclotomic],
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
"Norm of a principal ideal in a ring of quadratic integers",
[IsPrincipalIdeal],
function(I)
local m;
m:=GeneratorsOfTwoSidedIdeal(I);
m:=m[1];
return Norm(AssociatedRing(I),m);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(RightTransversal,
"Coset representatives of an ideal in a ring of quadratic integers",
[IsPrincipalIdeal],
function(I)
local cosetreps,leaves,newleaves,isinset,norm,x,y,i,j,D;

if IsBound(I!.RightTransversal) then return I!.RightTransversal; fi;

D:=GeneratorsOfRing(AssociatedRing(I));
D:=D[1];
norm:=NormOfPrincipalIdealGenerator(I);
cosetreps:=[0];
leaves:=[0];

#################################
isinset:=function(S,x)
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
if not isinset(cosetreps,y) then
Add(cosetreps,y);  
Add(newleaves,y);
fi;
y:=x+D; 
if not isinset(cosetreps,y) then
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



