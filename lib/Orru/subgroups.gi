###################################################################
###################################################################
InstallGlobalFunction(HAP_GenericCongruenceSubgroup,
function(family,integer,ring,ideal)
local R, n, one, mat, type, G, I, fam;

# This function returns a template for a congruence subgroup G in some 
# classical group such as GL(n,R), SL(n,R), Sp_n(R) over a ring of integers R.
# The subgroup G is defined using some condition involving an ideal I<=R. 
# The input "family" is a string such as "SL", "GL", Sp". In this code we
# denote any such family by "SL".
# 
# Various components such as "membership, "tree", "cosetPos" will ineed to be 
# implemented using a specific secondary function that depends on the specific 
# family and ring.
#
# The *fast* method will use a formula to first implement "cosetPos" and then
# compute the coset "tree", The *slow* method will use a general brute force 
# algorithm to first comput the "tree" and then implement "cosetPos".

R:=ring;
n:=integer;
if IsInt(ideal) then I:=IdealByGenerators(R,[ideal]); 
else
    I:=ideal; 
fi;
if not IsBound(I) then return fail; 
fi;

one:=One(R);
mat:=one*IdentityMat(n);
fam:=Concatenation(family,"(",String(n),", ",Name(R)," )");

    type := NewType( FamilyObj([ mat ]),
                     IsGroup and
                     IsAttributeStoringRep and
                     IsFinitelyGeneratedGroup and
                     IsMatrixGroup and
                     IsHapCongruenceSubgroup);

G:=rec(
    ringOfIntegers:=R,		#Matrices in the subgroup G lie in the ring R.
    level:=I,                	#This ideal I<=R is used to define G. 
    fam:=fam,			#The string "SL(n,R)".
    dimension:=n,
    membership:= fail,          #true if a matrix g lies in G.
    membershipLight:=fail,      #true if a matrix g in SL(n,R) lies in G.
    gens:=fail,			#"Nice" generating set for SL(n,R).
    tree:= fail,   		#Coset tree of G with respect to gens.
    generators:= fail,		#Generating set for G. 
    index:=fail,   		#Index of G in SL(n,R).
    cosetRep:= fail,            #CosetRep(g) represents g*G for g in SL(n,R).
    cosetPos:= fail,	       	#cosetPos(g) is the position of the coset g*G. 
    ugrp:= Group(mat),		#The trivial (or vertex stabilizer) group.
    name:="Congruence subgroup");

ObjectifyWithAttributes(	G, type,
				DimensionOfMatrixGroup, n
							);

return G;
end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallMethod( HAPCongruenceSubgroupGamma0,
"for integer n and integer m",
[IsInt, IsInt],
function(n,m)
local G,sl,membership,membershipLight,CosetRep, CosetPos;
if not (n=2 and m>0) then TryNextMethod(); fi;

#The following implements G=Gamm0(m) < SL(2,Z)

sl:=SL(2,Integers);
G:=HAP_GenericCongruenceSubgroup("SL",2,Integers,m);

###################################################
membership:=function(g)
if not g in sl then return false; fi;
if not g[2][1] mod m = 0  then return false;
else return true; fi;
end;
###################################################
###################################################
membershipLight:=function(g)
if not g[2][1] mod m = 0  then return false;
else return true; fi;
end;
###################################################

G!.membership:=membership;
G!.membershipLight:=membershipLight;
G!.level:=m;

G!.ugrp:=Group([[1,0],[0,1]]);
G!.name:="CongruenceSubgroupGamma0";
if m=1 then
G!.index:=m;
else
G!.index:=m*Product(List(SSortedList(Factors(m)), p->1+1/p));
fi;

if IsPrimeInt(m) then   #NEEDS TO BE EXTENDED TO NONE PRIMES. THAT IS, NEED 
			#TO IMPLEMENT THE PROJECTIVE LINE P^1(Z/mZ)
###########################################
CosetPos:=function(g)
if g[1][1] mod m =0 then return m+1; fi;
return 1 +((g[2][1]*g[1][1]^-1) mod m);
end;
###########################################

###########################################
CosetRep:=function(g)
if g[1][1] mod m=0 then return [[0,-1],[1,0]]; fi;
return [[1,0],[(g[2][1]*g[1][1]^-1) mod m,1]];
end;
###########################################
G!.cosetRep:=CosetRep;
G!.cosetPos:=CosetPos;
fi;

G:=ObjectifyWithAttributes(G, TypeObj(G),
IsIntegerMatrixGroup, true,
IsFinite, false);
return G;

end);
###################################################################
###################################################################

###################################################################
###################################################################
#THERE SHOULD BE JUST ONE IMPLEMENTATION FOR ALL CASES. 
#NEED TO ADJUST THIS.
InstallMethod(HAPCongruenceSubgroupTree,
"Coset tree for congruence subgroup",
[IsHapCongruenceSubgroup],
function(G)
if not (G!.dimension=2 and G!.ringOfIntegers=Integers) then TryNextMethod(); fi;
HAP_SL2ZSubgroupTree_fast(G);
end);
###################################################################
###################################################################

###################################################################
###################################################################

InstallOtherMethod( RightTransversal,
"Right transversal of congruence subgroup G in SL(2,Z)",
[IsMatrixGroup, IsHapCongruenceSubgroup],
100000,
function(G,H)
if not (H!.dimension=2 and Name(G)="SL(2,Integers)") then TryNextMethod(); fi;
HAPCongruenceSubgroupTree(H);
return HAP_TransversalCongruenceSubgroups(G,H);
end);
###################################################################
###################################################################


