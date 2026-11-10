InstallGlobalFunction(HAP_Gamma0_SL3ZHeckeCosetReps,
    function(n,p)
    local Matrices, a, b, c, M;

    Matrices := [];

    for a in [1..p] do
        for b in [1..p] do
            M := [[1,0,0],[0,1,0],[(a-1)*n, b-1, p]];
            Matrices := Concatenation(Matrices, [M]);
        od;
    od;

    for c in [1..p] do
        M := [[1,0,0],[(c-1)*n,p,0],[0,0,1]];
        Matrices := Concatenation(Matrices, [M]);
    od;
    
    M := [[p,0,0],[0,1,0],[0,0,1]];
    Matrices := Concatenation(Matrices, [M]);

    return Matrices;
end);

InstallGlobalFunction(HAP_UnimodularDecomposition,
function(M)
local SumOfThree, CheckDecomposition, Dec, Check, NewDec;

    SumOfThree := function(Q)
        local ReducedQ, m, R, ReducedQmodm, K, x, v, Q_1, Q_2, Q_3, Dec, x_2;

        if DeterminantIntMat(Q) = 1 then
            return [Q];
        fi;    

        ReducedQ := TransposedMat(HermiteNormalFormIntegerMat(TransposedMat(Q)));

        if AbsInt(ReducedQ[1][1]) <> 1 then
            m := ReducedQ[1][1];
        elif AbsInt(ReducedQ[2][2]) <> 1 then
            m := ReducedQ[2][2];
        else
            m := ReducedQ[3][3];
        fi;

        if IsPrime(m) or m = 1 then
            R := Integers mod m;

            ReducedQmodm := List(ReducedQ, row -> List(row, x -> x mod m)) * One(R);

            K := NullspaceMat(ReducedQmodm);

            x := List(K[1], k -> Int(k));

            v := (x[1]/m)*Q[1]+(x[2]/m)*Q[2]+(x[3]/m)*Q[3];

            Q_1 := ShallowCopy(Q);
            Q_1[1] := v;

            Q_2 := ShallowCopy(Q);
            Q_2[2] := v;

            Q_3 := ShallowCopy(Q);
            Q_3[3] := v;

            Dec := [Q_1, Q_2, Q_3];
            Dec := Filtered(Dec, D -> DeterminantIntMat(D) <> 0);

            return Dec;
        else
            #Print("m is not prime\n");
            if AbsInt(ReducedQ[1][1]) <> 1 then
                x := [1,0,0];
            elif AbsInt(ReducedQ[2][2]) <> 1 then
                x := [(m-ReducedQ[2,1]) mod m,1,0];
            else
                x_2 := (m-ReducedQ[3][2]) mod m;
                x := [(m-ReducedQ[3][1]-ReducedQ[2][1]*x_2) mod m,x_2,1];
            fi;
            v := (x[1]/m)*Q[1]+(x[2]/m)*Q[2]+(x[3]/m)*Q[3];

            Q_1 := ShallowCopy(Q);
            Q_1[1] := v;

            Q_2 := ShallowCopy(Q);
            Q_2[2] := v;

            Q_3 := ShallowCopy(Q);
            Q_3[3] := v;

            Dec := [Q_1, Q_2, Q_3];
            Dec := Filtered(Dec, D -> DeterminantIntMat(D) <> 0);

            return Dec;
        fi;
    end;

    CheckDecomposition := function(Dec)
        local i;

        for i in [1..Length(Dec)] do
            if AbsInt(DeterminantIntMat(Dec[i])) <> 1 then
            return rec(
                check := false,
                pos := i
            );
            fi;
        od;
        return rec(
            check := true,
            pos := 0
        );
    end;

    Dec := SumOfThree(M);

    Check := CheckDecomposition(Dec);

    while not Check.check do
        NewDec := SumOfThree(Dec[Check.pos]);
        if Check.pos = 1 then
            Dec := Concatenation(NewDec, Dec{[2..Length(Dec)]});
        else
            Dec := Concatenation(Dec{[1..Check.pos - 1]}, NewDec, Dec{[Check.pos+1.. Length(Dec)]});
        fi;
        Check := CheckDecomposition(Dec);
    od;

    return Dec;
end);

InstallGlobalFunction(HAP_Gamma0_SL3ZHeckeLift,
function(ProjPlane,n,p,S)

    local npts, Betas, A, i, Q_x, beta, U, y;

    npts := Length(ProjPlane.Reps);
    Betas := HAP_Gamma0_SL3ZHeckeCosetReps(n, p);
    A := [];

    for i in S.survivors do
        A[i] := ListWithIdenticalEntries(npts, 0);
        Q_x := HAP_MatrixInSLn_Hermite(ProjPlane.Reps[i]);
        for beta in Betas do
            for U in HAP_UnimodularDecomposition(Q_x*beta) do
                y := Position(ProjPlane.Reps,ProjPlane.RepOf(TransposedMat(U)[1]));
                A[i][y] := A[i][y] + 1;
            od;
        od;
    od;

    return A;
end);

InstallGlobalFunction(Gamma0_SL3ZHeckeOperator,
function(n,p)

    local ProjPlane, A, M, S, T, i, g, gA, free, torsion, Tfree, Ttorsion, torsionDivisors;

    ProjPlane := FiniteProjectiveSpace(2,n);
    M := Gamma0_SL3ZModularSymbolRelations(ProjPlane,n);
    S := SparseSmithNormalForm(M,true);
    A := HAP_Gamma0_SL3ZHeckeLift(ProjPlane,n,p,S);

    T := [];
    
    for g in S.generators do
        gA := Sum(S.survivors, x -> g[x] * A[x]);
        Add(T, S.Coordinates(gA));
    od; 

    torsionDivisors := Filtered(S.divisors, d -> d <> 0);
    free := Filtered([1..Length(S.divisors)], i -> S.divisors[i] = 0);
    torsion := Filtered([1..Length(torsionDivisors)], i -> torsionDivisors[i] <> 0);

    Tfree := T{free}{free};
    Ttorsion := T{torsion}{torsion};

    return rec(Tfree := Tfree, Ttorsion := Ttorsion);
end);