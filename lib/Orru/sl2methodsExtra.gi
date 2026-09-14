InstallMethod(AmbientPosition,
    "Returns cosetPos(g) function for the congruence subgroup G",
    [ IsIntegerMatrixGroup and IsHAPCongruenceSubgroupGamma0 ],
    function(G)
        local cosetPos, canonicalRep, n, countRep, count, offset, rank, divs, U, i, q, d, j, e;

        if DimensionOfMatrixGroup(G) <> 2 then
            TryNextMethod();
        fi;

        n := LevelOfCongruenceSubgroup(G);

        if IsPrime(n) then
            TryNextMethod();
        fi;

        countRep := function(m)
            local d, q, count;

            count := [1];

            for d in DivisorsInt(m) do
                q := m/d;
                Add(count, q*Phi(Gcd(d,q))/Gcd(d,q));
            od;

            return count;
        end;

        canonicalRep := function(g)
            local v, vv, d, dd, x, y;
            v := [g[1][1], g[2][1]];
            vv := List(v, x -> x mod n);
            if vv[1] mod n = 0 then
                return [0,1];
            elif Gcd(vv[1] mod n, n) = 1 then
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

        count := countRep(n);
        divs := DivisorsInt(n);

        offset := [];
        rank := [];

        for e in [2..Length(divs)-1] do
            d := divs[e];
            q := n / d;

            U := Unique(Filtered([1..n], i -> Gcd(i,d) = 1) mod q);

            offset[e] := Sum(count{[1..e]});
            rank[e] := [];

            for j in [1..Length(U)] do
                rank[e][U[j] + 1] := j;   # residue 0 goes in GAP position 1
            od;
        od;

        cosetPos := function(g)
            local w, e, U, q;

            w := canonicalRep(g);
            
            if w[1] = 0 then
                return 1;
            elif w[1] = 1 then
                U := [0..n-1];
                return 1 + Position(U,w[2]);
            else
                e := Position(divs, w[1]);
                q := n / w[1];
                return offset[e] + rank[e][(w[2] mod q) + 1];
            fi;
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
