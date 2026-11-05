# B_n(p) for SL(n,Z[1/p]), using the conventions of buildin.g.
# After LoadPackage("hap"); B := BruhatTitsBuildingSL(6,2);

# Exact membership in a cell stabilizer; no enumeration of infinite groups.
InstallMethod(\in, [IsMatrix, IsGroup and HasBTBuildingLatticeBases],
    function(g,h)
        local bases, n;
        bases := BTBuildingLatticeBases(h); n := Length(bases[1]);
        return Length(g)=n and ForAll(g,r -> Length(r)=n)
            and ForAll(Concatenation(g),IsRat) and DeterminantMat(g)=1
            and ForAll(bases,a -> ForAll(Concatenation(a*g*a^-1),IsInt));
    end);

InstallGlobalFunction(BruhatTitsBuildingSL, function(n,p)
    local one, elementary, pairs, gens, a, cells, rotation, powers,
          kernel, torus, flagCache, stabilizers, principalGenerators,
          sl2EdgeGenerators, flagGenerators, dimension, checkCell,
          boundary, stabilizer, u, v, i, g, diagonal;
    if not IsInt(n) or n < 2 then Error("n must be an integer >= 2"); fi;
    if not IsInt(p) or not IsPrimeInt(p) then Error("p must be prime"); fi;
    one := IdentityMat(n);
    elementary := function(i,j,x)
        local e;
        e := StructuralCopy(one); e[i][j] := x; return e;
    end;
    pairs := Filtered(Tuples([1..n],2),ij -> ij[1]<>ij[2]);
    gens := List(pairs,ij -> elementary(ij[1],ij[2],1));

    # Row bases a_i=diag(1,...,1,p,...,p), with i trailing p's.
    # Degree d representatives are the (d+1)-subsets, oriented increasingly.
    a := List([0..n-1],i -> DiagonalMat(
        Concatenation(List([1..n-i],j -> 1),List([1..i],j -> p))));
    cells := List([1..n],k -> Combinations([1..n],k));
    rotation := NullMat(n,n); rotation[1][n] := p;
    for i in [2..n] do rotation[i][i-1] := 1; od;
    powers := List([0..n-1],i -> rotation^i);
    kernel := fail; flagCache := rec();
    stabilizers := List([1..n],i -> []);

    # Generators of Gamma_n(p), n>=3. Sury--Venkataramana's formula,
    # as stated in Detinko--Flannery--Hulpke, Proposition 1.10:
    # https://arxiv.org/pdf/1906.10423 (see BUILDING_GENERAL.md).
    principalGenerators := function()
        local sigma, upper, ij, c, k;
        sigma := [one];
        for ij in Combinations([1..n],2) do
            c := StructuralCopy(one);
            c[ij[1]] := one[ij[2]]; c[ij[2]] := one[ij[1]];
            Add(sigma,c);                    # coordinate transposition
        od;
        for k in [1..n-1] do
            c := StructuralCopy(one);
            c[k][k] := -1; c[k+1][k+1] := -1; c[k+1][k] := 1;
            Add(sigma,c);
        od;
        upper := List(Combinations([1..n],2),ij -> elementary(ij[1],ij[2],p));
        return Set(Concatenation(List(sigma,c -> List(upper,t -> c^-1*t*c))));
    end;

    # Lifts of generators of the determinant-one diagonal torus mod p.
    # The 2x2 block has determinant 1 over Z, not merely modulo p.
    torus := [];
    if p > 2 and n > 2 then
        u := PrimitiveRootMod(p); v := PowerModInt(u,-1,p^2);
        for i in [1..n-1] do
            g := StructuralCopy(one);
            g[i][i] := u; g[i][n] := p;
            g[n][i] := (u*v-1)/p; g[n][n] := v;
            Add(torus,g);
        od;
    fi;

    # Rank 2: Schreier generators for the row-line stabilizer in SL(2,Z).
    # Transversal rows represent [1,t] (0<=t<p) and [0,1].
    sl2EdgeGenerators := function()
        local reps, result, r, s, x, pos;
        reps := List([0..p-1],t -> [[1,t],[0,1]]);
        Add(reps,[[0,1],[-1,0]]); result := [];
        for r in reps do
            for s in gens do
                x := r*s;
                if x[1][1] mod p=0 then pos := p+1;
                else pos := (x[1][2]*PowerModInt(x[1][1],-1,p)) mod p+1; fi;
                x := x*reps[pos]^-1;
                if x<>one then AddSet(result,x); fi;
            od;
        od;
        return result;
    end;

    # Full preimage of a standard row-flag stabilizer in SL(n,F_p).
    # cuts lists the dimensions of the prefix subspaces in that flag.
    flagGenerators := function(cuts)
        local key, allowed;
        if IsEmpty(cuts) then return gens; fi;
        key := String(cuts);
        if not IsBound(flagCache.(key)) then
            if n=2 then flagCache.(key) := sl2EdgeGenerators();
            else
                if kernel=fail then kernel := principalGenerators(); fi;
                allowed := Filtered(pairs,ij ->
                    ForAll(cuts,d -> not (ij[1]<=d and d<ij[2])));
                flagCache.(key) := Union(kernel,torus,
                    List(allowed,ij -> elementary(ij[1],ij[2],1)));
            fi;
        fi;
        return flagCache.(key);
    end;

    dimension := function(d)
        if IsInt(d) and d in [0..n-1] then return Length(cells[d+1]); fi;
        return 0;
    end;
    checkCell := function(d,k)
        if not IsInt(k) or k=0 or AbsInt(k)>dimension(d) then
            Error("Cell index out of range");
        fi;
    end;
    boundary := function(d,k)
        local c;
        checkCell(d,k);
        if d=0 then return []; fi;
        c := cells[d+1][AbsInt(k)];
        return Reversed(List([1..Length(c)],j ->
            [SignInt(k)*(-1)^(j-1)*Position(cells[d],Difference(c,[c[j]])),1]));
    end;
    stabilizer := function(d,k)
        local c, shift, cuts, r, h;
        checkCell(d,k); k := AbsInt(k);
        if not IsBound(stabilizers[d+1][k]) then
            c := cells[d+1][k]; shift := c[1]-1;
            # Rotate the first vertex to v0, compute its mod-p flag,
            # then conjugate the integral stabilizer back to the given cell.
            cuts := Set(List(c{[2..Length(c)]},j -> n-(j-c[1])));
            r := powers[shift+1];
            h := Group(List(flagGenerators(cuts),g -> r^-1*g*r));
            SetBTBuildingLatticeBases(h,a{c}); SetIsFinite(h,false);
            stabilizers[d+1][k] := h;
        fi;
        return stabilizers[d+1][k];
    end;
    diagonal := List([1..n],i -> 1); diagonal[1] := p; diagonal[2] := 1/p;
    return Objectify(HapNonFreeResolution,rec(
        group := Group(Concatenation(gens,[DiagonalMat(diagonal)])),
        elts := [one], dimension := dimension, boundary := boundary,
        stabilizer := stabilizer,
        action := function(d,k,g) return 1; end,
        homotopy := fail,
        cellRepresentatives := List(cells,level -> List(level,c -> a{c})),
        properties := [["length",n-1],["characteristic",0],["type","resolution"]]
    ));
end);
