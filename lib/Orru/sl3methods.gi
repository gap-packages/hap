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
        local cosetPos, n, ProjPlane, U, UU;
        if DimensionOfMatrixGroup(G) <> 3 then
            TryNextMethod();
        fi;

        n := LevelOfCongruenceSubgroup(G);
    
        ProjPlane := ProjectiveSpace(G);

        if not IsBound(G!.Units) then
	UU := Units(Integers mod n); #Graham
        G!.Units:=List(U,Int);  #Graham
        fi;
        
	cosetPos := function(g)
            local v, vv, u, w, p;
            v := [g[1][1], g[2][1], g[3][1]];
            vv := List(v, x -> x mod n);
         
            for u in G!.Units do  #Graham
                w := List(vv, x -> (u*x) mod n); #Graham
                p:=Position(ProjPlane.Reps,w);   #Graham
                if not p=fail then return p; fi; #Graham
            od;
        end;

        return cosetPos;
    end);
##########################################################################
##
## AmbientRepresentation( <G> )
##
## Returns a function cosetPos(g) giving a canonical rpresentative of the 
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
                IsHapRightTransversalSL2ZSubgroup and IsList and  #SL2???
                IsDuplicateFreeList and IsAttributeStoringRep ),
                rec( group := GG,
                     subgroup := G,
                     cosets:=transversal,
                     poscan:=poscan 
                ));
end);
