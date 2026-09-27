##########################################################################
##
## Methods for 4x4 congruence subgroups of SL4

##########################################################################
##
## ProjectiveSpace( <G> )
##

InstallMethod(ProjectiveSpace,
     "Projective space",
     [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
     function(G)
     local n;
        if DimensionOfMatrixGroup(G)<>4 then TryNextMethod(); fi;

        n := LevelOfCongruenceSubgroup(G);

        return FiniteProjectiveSpace3(n);
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
        local cosetPos, canonicalRep, n, ProjSpace;
        if DimensionOfMatrixGroup(G) <> 4 then
            TryNextMethod();
        fi;

        n := LevelOfCongruenceSubgroup(G);

        canonicalRep := function(g)
             local x, y, z, t, d_x, q_x, a, d_y, q_y, y_0, b, d_z, q_z, z_0, c, t_0, d, q;

            x := g[1][1] mod n;
            y := g[2][1] mod n;
            z := g[3][1] mod n;
            t := g[4][1] mod n;

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

            d_z := Gcd(z,n);
            q_z := n/d_z;
            z_0 := b*(z/d_z) mod Gcd(LcmInt(q_x,q_y),q_z);

            while not Gcd(z_0, q_z) = 1 do
                z_0 := z_0 + Gcd(LcmInt(q_x,q_y),q_z);
            od;

            c := ChineseRem([LcmInt(q_x,q_y),q_z],[b, Gcdex(z/d_z, q_z).coeff1*z_0 mod q_z]);

            d := Gcd(d_x,d_y,d_z);
            q := n/d;

            t_0 := c*t mod q;
            while not Gcd(t_0, d) = 1 do
                t_0 := t_0 + q;
            od;

            return [d_x mod n, d_y*y_0 mod n, d_z*z_0 mod n, t_0];
        end;
        
        ProjSpace := ProjectiveSpace(G);

        cosetPos := function(g)
            return Position(ProjSpace.Reps, canonicalRep(g));
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
        local MatrixInSL4_Hermite, cosetOfInt, cosetRep, n, ProjSpace, cosetPos;
        if DimensionOfMatrixGroup(G) <> 4 then
            TryNextMethod();
        fi;

        n := LevelOfCongruenceSubgroup(G);

        MatrixInSL4_Hermite := function(v)
            local Herm;
            Herm := HermiteNormalFormIntegerMatTransform([[v[1]],[v[2]],[v[3]],[v[4]]]);
            return Inverse(Herm!.rowtrans);
        end;

        ProjSpace := ProjectiveSpace(G);

        cosetOfInt:=function(i)
            local x,y,z,t;
            x := ProjSpace.Reps[i][1];
            y := ProjSpace.Reps[i][2];
            z := ProjSpace.Reps[i][3];
            t := ProjSpace.Reps[i][4];

            return MatrixInSL4_Hermite([x,y,z,t]);
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
    local n, GG, poscan, cosetPos, transversal, ProjSpace, cosetOfInt, MatrixInSL4_Hermite;

    if DimensionOfMatrixGroup(G) <> 4 then
        TryNextMethod();
    fi;

    n := LevelOfCongruenceSubgroup(G);
    ProjSpace := ProjectiveSpace(G);
    
    GG:=AmbientGroupOfCongruenceSubgroup(G);
    cosetPos:=AmbientPosition(G);

    MatrixInSL4_Hermite := function(v)
        local Herm;
        Herm := HermiteNormalFormIntegerMatTransform([[v[1]],[v[2]],[v[3]],[v[4]]]);
        return Inverse(Herm!.rowtrans);
    end;

    cosetOfInt:=function(i)
        local x,y,z,t;
        x := ProjSpace.Reps[i][1];
        y := ProjSpace.Reps[i][2];
        z := ProjSpace.Reps[i][3];
        t := ProjSpace.Reps[i][4];
        return MatrixInSL4_Hermite([x,y,z,t]);
        end;

    poscan := function(g)
        return cosetPos(g^-1);   
    end;

    transversal := List([1..Length(ProjSpace.Reps)],i->cosetOfInt(i)^-1);

    return Objectify( NewType( FamilyObj( GG ),
                IsHapRightTransversalSLnZSubgroup and IsList and  #SL2???
                IsDuplicateFreeList and IsAttributeStoringRep ),
                rec( group := GG,
                     subgroup := G,
                     cosets:=transversal,
                     poscan:=poscan 
                ));
end);
