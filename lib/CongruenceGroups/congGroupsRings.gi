#############################################################################
##
## This file contains constructors of congruence subgroups for Bianchi groups 
## and the like,  and methods for general congruence subgroups. Faster methods 
## for specific groups of interest are installed elsewhere.

#############################################################################
##
## CongruenceSubgroupGamma0( < GG >, < n >)
InstallOtherMethod( CongruenceSubgroupGamma0, "for integer matrix group and positive integer",
    [ IsMatrixGroup, IsIdealOfQuadraticIntegers ],
    function(GG,I)
    local type, G, membership, ambientMembership,
    CosetPos,CosetRep,CosetReps,R,one,zero;
    type := NewType( FamilyObj([One(GG)]),
                     IsGroup and
                     IsAttributeStoringRep and
                     IsFinitelyGeneratedGroup and
                     IsMatrixGroup and
                     IsHAPCongruenceSubgroup);
    G := rec();
    ObjectifyWithAttributes( G, type,
      DimensionOfMatrixGroup, DimensionOfMatrixGroup(GG),
      AmbientGroupOfCongruenceSubgroup, GG,
      OneImmutable, One(GG),
      IsIntegerMatrixGroup, IsIntegerMatrixGroup(GG),
      IsFinite, IsFinite(GG),
      StabilizerSubgroup, Group(One(GG)),
      LevelOfCongruenceSubgroup, I,  
      IsHAPPrincipalCongruenceSubgroup, false,
      IsHAPCongruenceSubgroupGamma0, true,
      IsHAPCongruenceSubgroupGamma1, false,
      IsHAPConjugatedCongruenceSubgroup, false,
      Characteristic,0
      );
    G!.bianchiInteger:=AssociatedRing(I)!.bianchiInteger;     

    membership:=function(g)
        local i;

        if not g in GG then   #WARNING: This will not work. Must fix it.
            return false;
        fi;

        ## Gamma0: first column is (*,0,...,0) mod n.
        for i in [2..DimensionOfMatrixGroup(GG)] do
            if not g[i][1] mod I = 0 then
                return false;
            fi;
        od;

        return true;
    end;
    ###################################################
    ###################################################
    ambientMembership:=function(g)
        local i;

        for i in [2..DimensionOfMatrixGroup(GG)] do
            if not g[i][1] mod I = 0 then
                return false;
            fi;
        od;

        return true;
    end;

    G!.membership:=membership;
    G!.ambientMembership:=ambientMembership;

if HasQuadraticArithmetic(I) then
   one:=One(HAPSqrt(AssociatedRing(I)!.bianchiInteger));
   zero:=Zero(one);
   SetQuadraticArithmetic(G,true);
else
   one:=One(I!.GeneratorsOfTwoSidedIdeal[1]);
   zero:=Zero(one);
fi;

    ##################################################
    ##
    ## Remaining components to be computed in other functions.
    #ProjectiveSpace(G);
    if DimensionOfMatrixGroup(GG)<>2 or not IsPrime(I) then
    AmbientTransversal(G);
    AmbientPosition(G); #A generic method will be used to construct the
    AmbientRepresentation(G); #these two functions except for cases with a
                         #bespoke implementation.
                         #cosetPos and poscan are used interchangably.
    return G;
    fi;

    ##################################################
    ##
    ## We'll at least get the 2x2 case with I prime working well.
    ## But this section needs to be rewritten for the nxn case 
    ## and maybe also the non-prime ideal case.
    if DimensionOfMatrixGroup(GG)=2 and IsPrime(I) then

R:=RightTransversal(I);
###########################################
CosetPos:=function(g)
if g[1][1] mod I =zero then return 1+Norm(I); fi;
return Position(R, ((g[2][1]*InverseOp(I,g[1][1])) mod I));
end;
###########################################

###########################################
CosetRep:=function(g)
if g[1][1] mod I=zero then return one*[[0,-1],[1,0]]; fi;
return one*[[1,0],[(g[2][1]*InverseOp(g[1][1])) mod I,one]];
end;
###########################################

CosetReps:=List([1..Norm(I)],i->[one*[1,0],one*[R[i],1]]);
Add(CosetReps,one*[[0,-1],[1,0]]);
G!.AmbientRepresentation:=CosetRep;
#G!.cosetReps:=CosetReps;
G!.AmbientPosition:=CosetPos;
G!.AmbientTransversal:=List(CosetReps,x->x^-1);
return G;

    fi;
        end);


