if LoadPackage("hap") <> true then Error("HAP is required"); fi;

ReadPackage("hap", "tst/Wunden/building_sl3_reference.g");;                 # Regression comparison; file unchanged.

TestBuildingGeneral := function()
    local check, parameters, n, p, b, old, d, k, c, h, g, t, f, sums,
          integral, image, cuts, blocks, ends, expected, glOrder,
          start, elapsed, a, kernelRows, kernelMatrices, j, x;
    check := function(ok,message)
        if not ok then Error(message); fi;
    end;
    glOrder := function(m,q)
        return Product([0..m-1],i -> q^m-q^i);
    end;
    old := BruhatTitsBuildingSL3ZHalf();
    for parameters in [[2,2],[2,3],[2,5],[2,7],[3,2],[3,3],[3,5],
                       [4,2],[4,3],[5,2],[5,3],[6,2],[6,3],[6,5]] do
        n := parameters[1]; p := parameters[2]; start := Runtime();
        b := BruhatTitsBuildingSL(n,p);
        check(IsHapNonFreeResolution(b) and Length(b)=n-1,"HAP type/length");
        check(List([0..n-1],b!.dimension)=List([1..n],k -> Binomial(n,k)),
            "Cell orbit counts");
        check(b!.dimension(-1)=0 and b!.dimension(n)=0,"Dimension outside range");
        for d in [0..n-1] do
            for k in [1..b!.dimension(d)] do
                c := b!.cellRepresentatives[d+1][k];
                check(b!.boundary(d,-k)=List(b!.boundary(d,k),t -> [-t[1],t[2]]),
                    "Signed boundary");
                if d>=2 then
                    sums := List([1..b!.dimension(d-2)],i -> 0);
                    for t in b!.boundary(d,k) do
                        for f in b!.boundary(d-1,t[1]) do
                            sums[AbsInt(f[1])] := sums[AbsInt(f[1])]+SignInt(f[1]);
                        od;
                    od;
                    check(ForAll(sums,x -> x=0),"Boundary squared");
                fi;
                if d=1 then
                    check(Sum(b!.boundary(d,k),t -> SignInt(t[1]))=0,"Augmentation");
                fi;
                h := b!.stabilizer(d,k);
                check(IsIdenticalObj(h,b!.stabilizer(d,-k)),"Stabilizer cache");
                check(not IsFinite(h) and IdentityMat(n) in h,"Infinite stabilizer");
                for g in GeneratorsOfGroup(h) do
                    check(DeterminantMat(g)=1 and
                        ForAll(c,a -> ForAll(Concatenation(a^-1*g*a),IsInt)),
                        "Generators preserve every lattice");
                    for t in b!.boundary(d,k) do
                        check(g in b!.stabilizer(d-1,t[1]),"Face stabilizer inclusion");
                    od;
                od;
                if n=3 and p=2 then
                    check(b!.boundary(d,k)=old!.boundary(d,k),"Original boundaries");
                    check(c=old!.cellRepresentatives[d+1][k],"Original representatives");
                    check(ForAll(GeneratorsOfGroup(h),g -> TransposedMat(g) in old!.stabilizer(d,k))
                        and ForAll(GeneratorsOfGroup(old!.stabilizer(d,k)),g -> TransposedMat(g) in h),
                        "Original stabilizers, both inclusions");
                fi;
                # Independently check orders of the finite parabolic images.
                # Normalize at the first vertex, then use the block-size formula.
                if n<=4 or d=n-1 then
                    a := c[1];
                    integral := List(GeneratorsOfGroup(h),g -> a^-1*g*a);
                    image := Group(List(integral,g -> g*One(GF(p))));
                    cuts := [];
                    for j in [2..Length(c)] do
                        x := c[j]*a^-1;
                        # Non-unit entries p mark the complement of a subspace.
                        Add(cuts,Number([1..n],i -> x[i][i]=1));
                    od;
                    ends := Concatenation([0],Set(cuts),[n]);
                    blocks := List([2..Length(ends)],i -> ends[i]-ends[i-1]);
                    expected := p^((n^2-Sum(blocks,i -> i^2))/2)
                        *Product(blocks,i -> glOrder(i,p))/(p-1);
                    check(Size(image)=expected,"Finite parabolic image order");
                fi;
            od;
        od;
        # Check that the congruence generators visible in the chamber have
        # the full trace-zero image Gamma(p)/Gamma(p^2), independently of
        # the installed infinite-group membership method.
        if n>=3 then
            kernelMatrices := Filtered(GeneratorsOfGroup(b!.stabilizer(n-1,1)),
                g -> ForAll(Concatenation(g-IdentityMat(n)),x -> x mod p=0));
            kernelRows := List(kernelMatrices,
                g -> Concatenation((g-IdentityMat(n))/p)*One(GF(p)));
            check(RankMat(kernelRows)=n^2-1,"Full first congruence layer");
        fi;
        # Negative membership: a rational determinant-one matrix with a
        # forbidden denominator, and one that fails a chamber congruence.
        g := IdentityMat(n); g[1][n] := 1/(p+1);
        check(not g in b!.stabilizer(0,1),"Reject other denominators");
        g := IdentityMat(n); g[n][1] := 1;
        check(g in b!.stabilizer(0,1) and not g in b!.stabilizer(n-1,1),
            "Reject chamber non-member");
        elapsed := Runtime()-start;
        Print("n=",n,", p=",p,": all ",2^n-1," stabilizers checked (",elapsed," ms).\n");
    od;
    Print("General building tests passed.\n");
end;;
TestBuildingGeneral();
