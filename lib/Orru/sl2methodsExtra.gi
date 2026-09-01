# methods for power of primes

InstallMethod(AmbientPosition,
    "Returns cosetPos(g) function for the congruence subgroup G",
    [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
    function(G)
        local cosetPos, canonicalRep, n, countRep;

        if DimensionOfMatrixGroup(G) <> 2 then
            TryNextMethod();
        fi;

        n := LevelOfCongruenceSubgroup(G);

        if not IsPrimePowerInt(n) then
            TryNextMethod();
        fi;

        if IsPrime(n) then
            TryNextMethod();
        fi;

        canonicalRep := function(g)
            local v, vv, U, d, dd, x, y;
            v := [g[1][1], g[2][1]];
            vv := List(v, x -> x mod n);
            U := Filtered([0..n],i->Gcd(i,n)=1);
            if vv[1] mod n = 0 then
                return [0,1];
            elif vv[1] mod n in U then
                return [1,(Inverse(vv[1]) mod n)*vv[2] mod n];
            else
                d := Gcd(vv[1],n);
                dd := n/d;
                x := vv[1]/d;
                y := vv[2]/x mod dd;
                while not Gcd(d,y) = 1 do
                    y := y + dd;
                od;
                return [d, y];
            fi;
        end;

        countRep := function(m)
            local p, e, i, countp;

            p := Set(Factors(m))[1];
            e := Length(Factors(m));

            countp := [1,m];

            for i in [2..e] do
                Add(countp, p^(e-i)*(p-1));
            od;

            return countp;
        end;

        cosetPos := function(g)
            local w, count, e, U;
            w := canonicalRep(g);

            count := countRep(n);

            if w[1] = 0 then
                return 1;
            elif w[1] = 1 then
                U := [0..n-1];
                return 1 + Position(U,w[2]);
            else
                e := Length(Factors(w[1]));
                U := Filtered([0..n/w[1]], i-> Gcd(i,n/w[1]) = 1);
                return Sum(count{[1..1 + e]}) + Position(U, w[2]);
            fi;
        end;

        return cosetPos;
    end);

InstallMethod(AmbientPosition,
    "Returns cosetPos(g) function for the congruence subgroup G",
    [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
    function(G)
        local cosetPos, canonicalRep, n, ProjLine, U;

        if DimensionOfMatrixGroup(G) <> 2 then
            TryNextMethod();
        fi;

        n := LevelOfCongruenceSubgroup(G);

        if IsPrime(n) then
            TryNextMethod();
        fi;

        if IsPrimePowerInt(n) then
            TryNextMethod();
        fi;

        ProjLine := ProjectiveSpace(G);

        U := Filtered([0..n],i -> Gcd(i,n) = 1);

        canonicalRep := function(g)
            local v, vv, d, dd, x, y;
            v := [g[1][1], g[2][1]];
            vv := List(v, x -> x mod n);
            if vv[1] mod n = 0 then
                return [0,1];
            elif (vv[1] mod n) in U then
                return [1,(Inverse(vv[1]) mod n)*vv[2] mod n];
            else
                d := Gcd(vv[1],n);
                dd := n/d;
                x := vv[1]/d;
                y := vv[2]/x mod dd;
                while not Gcd(d,y) = 1 do
                    y := y + dd;
                od;
                return [d, y];
            fi;
        end;

        cosetPos := function(g)
            local w;
            w := canonicalRep(g);
            return Position(ProjLine,w);
        end;

        return cosetPos;
    end);

InstallMethod(AmbientRepresentation,
    "Returns cosetPos(g) function for the congruence subgroup G",
    [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
    function(G)
        local cosetOfInt, cosetRep, n, ProjLine, cosetPos;

        if DimensionOfMatrixGroup(G) <> 2 then
            TryNextMethod();
        fi;

        n := LevelOfCongruenceSubgroup(G);

        if IsPrime(n) then
            TryNextMethod();
        fi;

        ProjLine := ProjectiveSpace(G);
        
        cosetOfInt := function(i)
            local a, c, b, d, gg;
            a := ProjLine[i][1];
            c := ProjLine[i][2];
            if a = 0 then
                return [[0,-1],[1,0]];
            fi;
            gg := Gcdex(a,c);
            b := -gg.coeff2;
            d :=  gg.coeff1;
            return [[a,b],[c,d]];
        end;

        cosetPos := AmbientPosition(G);

        cosetRep:=function(g);
            return cosetOfInt(cosetPos(g));
        end;

        return cosetRep;
    end);

    ##########################################################################
##
## AmbientTransversal( <G> )
##
## Right transversal for a congruence subgroup G in its ambient group GG
     InstallMethod(AmbientTransversal,
     "Right transversal for a congruence subgroup G in its ambient group",
     [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
     function(G)
        local n, GG, poscan, cosetPos, transversal, ProjLine, cosetOfInt;
        if DimensionOfMatrixGroup(G) <> 2 then
            TryNextMethod();
        fi;

        n:=LevelOfCongruenceSubgroup(G);

        if IsPrime(n) then
            TryNextMethod();
        fi;

        ProjLine := ProjectiveSpace(G);
        
        GG:=AmbientGroupOfCongruenceSubgroup(G);

        cosetPos:=AmbientPosition(G);

        cosetOfInt := function(i)
            local a, c, b, d, gg;
            a := ProjLine[i][1];
            c := ProjLine[i][2];
            if a = 0 then
                return [[0,-1],[1,0]];
            fi;
            gg := Gcdex(a,c);
            b := -gg.coeff2;
            d :=  gg.coeff1;
            return [[a,b],[c,d]];
        end;

        poscan := function(g)
            return cosetPos(g^-1);
        end;

        transversal := List([1..Length(ProjLine)],i->cosetOfInt(i)^-1);

        return Objectify( NewType( FamilyObj( GG ),
                    IsHapRightTransversalSLnZSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
                    rec( group := GG,
                         subgroup := G,
                         cosets:=transversal,
                         poscan:=poscan 
                    ));
     end);
