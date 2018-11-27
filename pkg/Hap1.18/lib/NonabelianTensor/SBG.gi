#(C) Graham Ellis 2005-06

#####################################################################
InstallGlobalFunction(ThirdHomotopyGroupOfSuspensionB_alt,
function(G)
local H2,Primes, p, S,M,N,ans,T,K,KS,L,x,y;

H2:=GroupHomology(G,2);
ans:=[];
Primes:=SSortedList(Factors(Product(AbelianInvariants(G))));
Primes:=Filtered(Primes, i->not i=1);

for p in Primes do
S:=SylowSubgroup(G,p);
T:=NonabelianTensorSquare(S);
K:=Intersection(DerivedSubgroup(G),S);
K:=Elements(RightTransversal(K,DerivedSubgroup(S)));
KS:=Elements(RightTransversal(S,DerivedSubgroup(S)));
N:=[];
for x in K do
for y in KS do
Append(N,[T.pairing(x,y),T.pairing(y,x)]);
od;
od;
N:=Group(N);


L:=[];
for x in KS do
Append(L,[T.pairing(x,x)]);
od;
L:=Group(L);
M:=AbelianInvariants(L/Intersection(L,N));    
Append(ans,M);
od;

return  [H2,ans];
end);
#####################################################################

#####################################################################
InstallGlobalFunction(NonabelianSymmetricKernel_alt,
function(G)
local H2,Primes, p, S,M,N,ans,T,K,KS,L,x,y;

H2:=GroupHomology(G,2);
ans:=[];
Primes:=SSortedList(Factors(Product(AbelianInvariants(G))));
Primes:=Filtered(Primes, i->not i=1);

for p in Primes do
S:=SylowSubgroup(G,p);
T:=NonabelianSymmetricSquare(S);
K:=Intersection(DerivedSubgroup(G),S);
K:=Elements(RightTransversal(K,DerivedSubgroup(S)));
KS:=Elements(RightTransversal(S,DerivedSubgroup(S)));
N:=[];
for x in K do
for y in KS do
Append(N,[T.pairing(x,y),T.pairing(y,x)]);
od;
od;
N:=Group(N);


L:=[];
for x in KS do
Append(L,[T.pairing(x,x)]);
od;
L:=Group(L);
M:=AbelianInvariants(L/Intersection(L,N));
Append(ans,M);
od;

return  [H2,ans];
end);
#####################################################################

