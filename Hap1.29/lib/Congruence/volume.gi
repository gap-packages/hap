##################################
InstallGlobalFunction(QuadraticCharacter,
function(Q,m);
return Jacobi(Discriminant(Q),m);
end);
##################################

##################################
InstallGlobalFunction(Lfunction,
function(arg)
local Q,s,t,L,n;
Q:=arg[1];
s:=arg[2];
if Length(arg)=3 then t:=arg[3]; else t:=1000; fi;
L:=0;
for n in [1..t] do
L:=L+ QuadraticCharacter(Q,n)*n^(-s);
od;
return Int(10^8*L)/10^8;
end);
##################################

##################################
InstallOtherMethod(IsomorphismFpGroup,
"for HAP SL2Z subgroups",
[IsHapSL2ZSubgroup],
1000000,#Again, there must be a better way
function(H)
local QH, epi, ElementToWord, iso;
QH:=HAP_RightTransversalSL2ZSubgroups(H,false);
iso:=IsomorphismFpGroup(QH);
epi:=QH!.epimorphism;
ElementToWord:=QH!.ElementToWord;
return GroupHomomorphismByFunction(H, Image(iso), x->
                Image(iso,Image(epi,ElementToWord(x)) ) );


end);
##################################

#################################
InstallOtherMethod(AbelianInvariants,
"for HAP SL2Z subgroups",
[IsHapSL2ZSubgroup],
1000000,#Again, there must be a better way
function(H);
return AbelianInvariants(HAP_RightTransversalSL2ZSubgroups(H,false));
end);
##################################

#################################
InstallOtherMethod(IsPrime,
"for HAP SL2Z subgroups",
[IsIdealOfQuadraticIntegers],
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

