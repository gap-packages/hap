##########################################################################
##
## Methods for congruence subgroups Gamma0 of SL(m,Z) for arbitrary m.
## They use the generic HAP_FiniteProjectiveSpaceIntegers(dim, n), whose
## RepOf component is the canonical-representative *function*.
##
## Installed with rank -1 so that the bespoke methods in sl3methods.gi and
## sl4methods.gi are still preferred for m = 3, 4.

##########################################################################
##
## ProjectiveSpace( <G> )
##
InstallMethod(ProjectiveSpace,
     "Projective space P^(m-1)(Z/nZ) for Gamma0 in SL(m,Z)",
     [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ], -1,
     function(G)
     local n, m;
        n := LevelOfCongruenceSubgroup(G);
        m := DimensionOfMatrixGroup(G);
        return FiniteProjectiveSpace(m-1, n);
     end);

##########################################################################
##
## AmbientPosition( <G> )
##
## Returns a function cosetPos(g) giving the position of the coset gG in
## the ambient group.
InstallMethod(AmbientPosition,
"Returns cosetPos(g) function for Gamma0 in SL(m,Z)",
[ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ], -1,
    function(G)
        local cosetPos, canonicalRep, n, m, ProjSpace;

        n := LevelOfCongruenceSubgroup(G);
        m := DimensionOfMatrixGroup(G);
        ProjSpace := ProjectiveSpace(G);

        canonicalRep := function(g)
            return ProjSpace.RepOf(List([1..m], i -> g[i][1] mod n));
        end;

        cosetPos := function(g)
            return Position(ProjSpace.Reps, canonicalRep(g));
        end;

        return cosetPos;
    end);

##########################################################################
##
## HAP_MatrixInSLn_Hermite( <v> )
##
## For a unimodular integer vector v (mod n), returns a matrix in SL(m,Z)
## whose first column is v. Shared by AmbientRepresentation and
## AmbientTransversal below.
InstallGlobalFunction(HAP_MatrixInSLn_Hermite,
function(v)
    local Q;
    Q := HermiteNormalFormIntegerMatTransform(List(v, x -> [x]));
    return Inverse(Q!.rowtrans);
end);

##########################################################################
##
## AmbientRepresentation( <G> )
##
## Returns a function, cosetRep(g), giving a canonical representative of the
## coset gG in the ambient group.
InstallMethod(AmbientRepresentation,
"Returns cosetRep(g) function for Gamma0 in SL(m,Z)",
[ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ], -1,
     function(G)
        local cosetRep, ProjSpace, cosetPos;

        ProjSpace := ProjectiveSpace(G);
        cosetPos := AmbientPosition(G);

        cosetRep := function(g)
            return HAP_MatrixInSLn_Hermite(ProjSpace.Reps[cosetPos(g)]);
        end;

        return cosetRep;
     end);

##########################################################################
##
## AmbientTransversal( <G> )
##
## Right transversal for a congruence subgroup G in its ambient group GG
InstallMethod(AmbientTransversal,
"Right transversal for Gamma0 in SL(m,Z)",
[ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ], -1,
function(G)
    local GG, poscan, cosetPos, transversal, ProjSpace;

    ProjSpace := ProjectiveSpace(G);
    GG := AmbientGroupOfCongruenceSubgroup(G);
    cosetPos := AmbientPosition(G);

    poscan := function(g)
        return cosetPos(g^-1);
    end;

    transversal := List(ProjSpace.Reps, v -> HAP_MatrixInSLn_Hermite(v)^-1);

    return Objectify( NewType( FamilyObj( GG ),
                IsHapRightTransversalSLnZSubgroup and IsList and
                IsDuplicateFreeList and IsAttributeStoringRep ),
                rec( group := GG,
                     subgroup := G,
                     cosets:=transversal,
                     poscan:=poscan 
                ));
end);

