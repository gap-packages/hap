# Smith normal form of a sparse integer matrix.
#
# Rows are relations, columns are generators (as in ThirdHomologySparseMatrix),
# so coker = Z^cols / (row span) and the SNF describes that quotient.
#
# Method:
#   1. Sparse elimination of unit pivots. If a relation r has coefficient
#      e = +-1 on column b, then [b] = -e*(rest of r) in the cokernel: we record
#      this substitution, add multiples of r to the other relations to clear
#      column b, and delete r and b. Each step contributes one 1 to the SNF.
#   2. Dense SmithNormalFormIntegerMatTransforms on what is left (the
#      "survivor" columns), which is where any torsion lives.
#
# The full transforms P, Q are dense (rows x rows, cols x cols), so they are
# not returned. With transforms = true you get instead what is needed to work
# in the cokernel:
#   generators  : the generators g_j of the nontrivial summands Z/d_j
#                 (d_j <> 1), as dense vectors in Z^cols (rows of Q^-1);
#   divisors    : the d_j (0 = free summand);
#   Coordinates : v in Z^cols  ->  coordinates of the class of v in the g_j,
#                 reduced mod d_j (this is v*Q restricted to the g_j).

# SparseSmithNormalForm(S)  or  SparseSmithNormalForm(S, true)
# S is a HAP sparse matrix (S!.mat, S!.rows, S!.cols).
InstallGlobalFunction(SparseSmithNormalForm,
function(arg)
    local S, transforms, ncols, rels, cobounds, subs, alive, changed, order,
          k, row, cands, i, b, e, bnd, j, pos, c, x, survivors, posOf,
          smallM, v, nsurv, SNF, Qs, Qsinv, smallDiag, smallRank, nelim,
          idx, divisors, generators, Coordinates, result, r, SNF_NormalizeRow, SNF_AddRows;

          # Sort a sparse row [[col,val],...], add up repeated columns, drop zeros.
    SNF_NormalizeRow := function(row)
        local out, x;

        out := [];
        for x in SortedList(row) do
            if Length(out) > 0 and out[Length(out)][1] = x[1] then
                out[Length(out)][2] := out[Length(out)][2] + x[2];
            else
                Add(out, [x[1], x[2]]);
            fi;
        od;
        return Filtered(out, x -> x[2] <> 0);
    end;

    # Sum of two normalized sparse rows.
    SNF_AddRows := function(B, Y)
        local U, pB, pY, lB, lY;

        lB := Length(B); lY := Length(Y);
        U := []; pB := 1; pY := 1;
        while pB <= lB and pY <= lY do
            if B[pB][1] < Y[pY][1] then
                Add(U, B[pB]); pB := pB + 1;
            elif B[pB][1] > Y[pY][1] then
                Add(U, Y[pY]); pY := pY + 1;
            else
                if B[pB][2] + Y[pY][2] <> 0 then
                    Add(U, [B[pB][1], B[pB][2] + Y[pY][2]]);
                fi;
                pB := pB + 1; pY := pY + 1;
            fi;
        od;
        Append(U, B{[pB..lB]});
        Append(U, Y{[pY..lY]});
        return U;
    end;

    S := arg[1];
    transforms := Length(arg) > 1 and arg[2] = true;
    ncols := S!.cols;

    rels := List(Filtered(S!.mat, r -> Length(r) > 0), SNF_NormalizeRow);
    rels := Filtered(rels, r -> Length(r) > 0);

    # cobounds[x] = relations in which column x (may) occur.
    cobounds := List([1..ncols], x -> []);
    for k in [1..Length(rels)] do
        for x in rels[k] do
            AddSet(cobounds[x[1]], k);
        od;
    od;

    subs := [];
    alive := ListWithIdenticalEntries(ncols, true);

    ########## 1. Sparse elimination of unit pivots ##########
    changed := true;
    while changed do
        changed := false;
        # Short relations first: less fill-in.
        order := Filtered([1..Length(rels)], k -> rels[k] <> 0);
        SortBy(order, k -> Length(rels[k]));
        for k in order do
            row := rels[k];
            if row = 0 then continue; fi;
            cands := Filtered([1..Length(row)], i -> AbsInt(row[i][2]) = 1);
            if cands = [] then continue; fi;
            # Among the unit entries, eliminate the column occurring least often.
            i := cands[PositionMinimum(List(cands, i -> Length(cobounds[row[i][1]])))];

            # row = e*[b] + rest = 0 with e = +-1, so [b] = -e*rest.
            b := row[i][1];
            e := row[i][2];
            bnd := List(Concatenation(row{[1..i-1]}, row{[i+1..Length(row)]}),
                        x -> [x[1], -e*x[2]]);
            Add(subs, [b, bnd]);
            alive[b] := false;
            rels[k] := 0;
            changed := true;

            # Substitute [b] = bnd in every other relation containing b.
            for j in cobounds[b] do
                if rels[j] = 0 then continue; fi;
                pos := PositionProperty(rels[j], x -> x[1] = b);
                if pos = fail then continue; fi;
                c := rels[j][pos][2];
                rels[j] := SNF_AddRows(
                    Concatenation(rels[j]{[1..pos-1]}, rels[j]{[pos+1..Length(rels[j])]}),
                    List(bnd, x -> [x[1], c*x[2]]));
                for x in bnd do
                    AddSet(cobounds[x[1]], j);
                od;
                if rels[j] = [] then rels[j] := 0; fi;
            od;
            cobounds[b] := [];
        od;
    od;
    nelim := Length(subs);

    ########## 2. Dense SNF on the survivor columns ##########
    survivors := Filtered([1..ncols], x -> alive[x]);
    nsurv := Length(survivors);
    posOf := [];
    for i in [1..nsurv] do
        posOf[survivors[i]] := i;
    od;

    smallM := [];
    for row in Filtered(rels, r -> r <> 0) do
        v := ListWithIdenticalEntries(nsurv, 0);
        for x in row do
            v[posOf[x[1]]] := x[2];
        od;
        Add(smallM, v);
    od;

    for i in [1..Length(smallM)] do
        r := smallM[i];
        r := SignInt(First(r, x -> x <> 0))*r;
        smallM[i] := r;
    od;

    smallM := Set(smallM);

    if nsurv = 0 or smallM = [] then
        Qs := IdentityMat(nsurv);
        smallDiag := [];
        smallRank := 0;
    else
        SNF := SmithNormalFormIntegerMatTransforms(smallM);
        Qs := SNF.coltrans;
        smallRank := SNF.rank;
        smallDiag := List([1..smallRank], j -> AbsInt(SNF.normal[j][j]));
    fi;

    result := rec(
        rows := S!.rows,
        cols := ncols,
        rank := nelim + smallRank,
        # Nonzero diagonal entries of the SNF, in divisibility order.
        elementaryDivisors := Concatenation(ListWithIdenticalEntries(nelim, 1), smallDiag),
        torsion := Filtered(smallDiag, d -> d <> 1),
        freeRank := ncols - nelim - smallRank,
        survivors := survivors,
        smallMatrix := smallM
    );

    if not transforms then
        return result;
    fi;

    ########## Generators and coordinates of the cokernel ##########
    # divisors[j] for the survivor-level Smith columns j = 1..nsurv.
    divisors := List([1..nsurv], j -> 0);
    for j in [1..smallRank] do
        divisors[j] := smallDiag[j];
    od;
    idx := Filtered([1..nsurv], j -> divisors[j] <> 1);

    if nsurv = 0 then Qsinv := []; else Qsinv := Qs^-1; fi;
    generators := [];
    for j in idx do
        v := ListWithIdenticalEntries(ncols, 0);
        for i in [1..nsurv] do
            v[survivors[i]] := Qsinv[j][i];
        od;
        Add(generators, v);
    od;

    Coordinates := function(w)
        local s, c, x, u, i;

        w := ShallowCopy(w);
        # Substitutions in elimination order: bnd of the j-th elimination
        # only involves survivors and columns eliminated later.
        for s in subs do
            c := w[s[1]];
            if c <> 0 then
                w[s[1]] := 0;
                for x in s[2] do
                    w[x[1]] := w[x[1]] + c*x[2];
                od;
            fi;
        od;
        if nsurv = 0 then return []; fi;
        u := (w{survivors}*Qs){idx};
        for i in [1..Length(idx)] do
            if divisors[idx[i]] <> 0 then
                u[i] := u[i] mod divisors[idx[i]];
            fi;
        od;
        return u;
    end;

    result.divisors := divisors{idx};
    result.generators := generators;
    result.Coordinates := Coordinates;
    return result;
end);
