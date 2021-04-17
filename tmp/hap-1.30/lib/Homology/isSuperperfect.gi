#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(IsSuperperfect,
function(G)
local
	C,P,primes,R,prm,F;

F:= function(X);
return TensorWithIntegersModP(X,prm);
end;


if Length(AbelianInvariants(G))>0 then return false; fi;
if not IsFinite(G) then X:=GroupHomology(G,2);
	if Length(X)>0 then return false; else return true; fi;
fi;
primes:=(Factors(Order(G)/Product(SSortedList(Factors(Order(G))))));
primes:=SSortedList(Reversed(primes));

for  prm in primes do
P:=SylowSubgroup(G,prm);
if Order(P)<=3200 then
R:=ResolutionFiniteGroup(P,3);
else
R:=ResolutionNormalSeries(BigStepLCS(P,9),3);
fi;
C:=F(R);
#if Homology(C,2)>0 then
if Length(PrimePartDerivedFunctor(G,R,F,2))>0 then return 
false; fi;
#fi;
od;


return true;

end);
#####################################################################

