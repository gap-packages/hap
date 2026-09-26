##########################################################################
##
## Methods for 3x3 congruence subgroups of SL3

##########################################################################
##
## ProjectiveSpace( <G> )
##

InstallMethod(ProjectiveSpace,
     "Projective space",
     [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
     function(G)
     local n;
        if DimensionOfMatrixGroup(G)<>3 then TryNextMethod(); fi;

        n := LevelOfCongruenceSubgroup(G);

        return FiniteProjectivePlane(n);
     end);
##########################################################################
##
## AmbientPosition( <G> )
##
## Returns a function cosetPos(g) giving the position of the coset gG in 
## the ambient group. 
InstallMethod(AmbientPosition,
"Returns cosetPos(g) function for the congruence subgroup G",
[ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
    function(G)
        local cosetPos, canonicalRep, n, ProjPlane;
        if DimensionOfMatrixGroup(G) <> 3 then
            TryNextMethod();
        fi;

        n := LevelOfCongruenceSubgroup(G);

        canonicalRep := function(g)
            local x, y, z, d_x, q_x, a, d_y, q_y, y_0, b, d, q, z_0;

            x := g[1][1] mod n;
            y := g[2][1] mod n;
            z := g[3][1] mod n;

            d_x := Gcd(x,n);
            q_x := n/d_x;
            a := Gcdex(x/d_x, q_x).coeff1 mod q_x;

            d_y := Gcd(y,n);
            q_y := n/d_y;
            y_0 := a*(y/d_y) mod Gcd(q_y,q_x);

            while not Gcd(y_0, q_y) = 1 do
                y_0 := y_0 + Gcd(q_y,q_x);
            od;

            b := ChineseRem([q_x,q_y],[a, Gcdex(y/d_y, q_y).coeff1*y_0 mod q_y]);

            d := Gcd(d_x,d_y);
            q := n/d;

            z_0 := b*z mod q;
            while not Gcd(z_0, d) = 1 do
                z_0 := z_0 + q;
            od;

            return [d_x mod n, d_y*y_0 mod n, z_0];
        end;
        
        ProjPlane := ProjectiveSpace(G);

        cosetPos := function(g)
            return Position(ProjPlane.Reps, canonicalRep(g));
        end;

        return cosetPos;
    end);
##########################################################################
##
## AmbientRepresentation( <G> )
##
## Returns a function, cosetRep(g), giving a canonical representative of the 
## coset gG in the ambient group. 
InstallMethod(AmbientRepresentation,
"Returns cosetRep(g) function for the congruence subgroup G",
[ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
     function(G)
        local MatrixInSL3_Hermite, cosetOfInt, cosetRep, n, ProjPlane, cosetPos;
        if DimensionOfMatrixGroup(G) <> 3 then
            TryNextMethod();
        fi;

        n := LevelOfCongruenceSubgroup(G);

        MatrixInSL3_Hermite := function(v)
            local Herm;
            Herm := HermiteNormalFormIntegerMatTransform([[v[1]],[v[2]],[v[3]]]);
            return Inverse(Herm!.rowtrans);
        end;

        ProjPlane := ProjectiveSpace(G);

        cosetOfInt:=function(i)
            local x,y,z;
            x := ProjPlane.Reps[i][1];
            y := ProjPlane.Reps[i][2];
            z := ProjPlane.Reps[i][3];

            return MatrixInSL3_Hermite([x,y,z]);
        end;

        cosetPos := AmbientPosition(G);

        cosetRep:=function(g);
            return cosetOfInt(cosetPos(g));
        end;

        return cosetRep;
     end);

     ##
## AmbientTransversal( <G> )
##
## Right transversal for a congruence subgroup G in its ambient group GG
InstallMethod(AmbientTransversal,
"Right transversal for a congruence subgroup G in its ambient group",
[ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
function(G)
    local n, GG, poscan, cosetPos, transversal, ProjPlane, cosetOfInt, MatrixInSL3_Hermite;

    if DimensionOfMatrixGroup(G) <> 3 then
        TryNextMethod();
    fi;

    n := LevelOfCongruenceSubgroup(G);
    ProjPlane := ProjectiveSpace(G);
    
    GG:=AmbientGroupOfCongruenceSubgroup(G);
    cosetPos:=AmbientPosition(G);

    MatrixInSL3_Hermite := function(v)
        local Herm;
        Herm := HermiteNormalFormIntegerMatTransform([[v[1]],[v[2]],[v[3]]]);
        return Inverse(Herm!.rowtrans);
    end;

    cosetOfInt:=function(i)
        local x,y,z;
        x := ProjPlane.Reps[i][1];
        y := ProjPlane.Reps[i][2];
        z := ProjPlane.Reps[i][3];
        return MatrixInSL3_Hermite([x,y,z]);
        end;

    poscan := function(g)
        return cosetPos(g^-1);   
    end;

    transversal := List([1..Length(ProjPlane.Reps)],i->cosetOfInt(i)^-1);

    return Objectify( NewType( FamilyObj( GG ),
                IsHapRightTransversalSLnZSubgroup and IsList and  #SL2???
                IsDuplicateFreeList and IsAttributeStoringRep ),
                rec( group := GG,
                     subgroup := G,
                     cosets:=transversal,
                     poscan:=poscan 
                ));
end);
