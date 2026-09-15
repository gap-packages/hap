#############################################################################
##
##  ResolutionFiniteGroup.g
##
##  Experimental performance-oriented rewrite of HAP's
##  ResolutionFiniteGroup.
##
##  Main changes:
##    * O(1) completion test via RemainingZeros
##    * dictionary lookup for group element -> index
##    * reduced allocation in Differential / Contraction
##    * single-pass FirstZero
##    * event-driven consequence propagation using residual counts
##
##  Usage:
##
##      LoadPackage("hap");
##      Read("ResolutionFiniteGroup.g");
##
##      G := SmallGroup(64,134);
##      R := ResolutionFiniteGroup(G,5);
##
##  The function is deliberately named ResolutionFiniteGroupFast so that
##  it can coexist with HAP's ResolutionFiniteGroup for benchmarking and
##  regression testing.
##
##  IMPORTANT:
##  The event-driven collapse order can differ from HAP's original scan
##  order.  The resulting resolution should therefore be validated on the
##  groups/degrees relevant to your application, especially with
##  tietze = true.
##
#############################################################################

InstallGlobalFunction(ResolutionFiniteGroup,
function(arg)

local
    R, Gens, K, tietze, G, AllElts, Elts, ExtendedElts, N, MT,
    ElementIndex, ActionIndex, InverseIndex,
    ChangeSign, MaxComplex, Dimension, Boundary, PseudoBoundary,
    ContractionMatrix, Contraction, Homotopy,
    ComputedContractions, Differential,
    FirstZero, NextResTerm,
    Spheres, DiffLengths, InitComputedContractions,
    RemainingZeros,
    saveSpace, Charact, AlgebraicRed,
    ExtendRes, Extendible,
    IncidenceByFace, Residual,
    ConsequenceQueue, ConsequenceQueueHead, ConsequenceQueued,
    InitConsequenceEngine, RegisterBoundaryRow,
    EnqueueConsequence, MarkMCEntry, ProcessConsequences,
    UniqueZeroFace, AddSphere,
    i, iso,
    AbsInt, SignInt;

AbsInt := AbsInt_HAP;
SignInt := SignInt_HAP;

#############################################################################
##
##  Arguments
##
#############################################################################

if IsGroup(arg[1]) then
    Gens := GeneratorsOfGroup(arg[1]);
    if Length(Gens) = 0 then
        Gens := [Identity(arg[1])];
    fi;
else
    Gens := arg[1];
fi;

Gens := SSortedList(StructuralCopy(Gens));
K := StructuralCopy(arg[2]);

if Length(arg) > 2 then
    tietze := arg[3];
else
    tietze := false;
fi;

Charact := 0;
if Length(arg) > 3 and IsInt(arg[4]) then
    Charact := arg[4];
fi;

if Length(arg) > 4 then
    if arg[5] = "extendible" then
        Extendible := true;
        saveSpace := false;
    else
        saveSpace := arg[5];
        Extendible := false;
    fi;
else
    saveSpace := false;
    Extendible := false;
fi;

G := GroupWithGenerators(Gens);
N := Order(G);

#############################################################################
##
##  Matrix groups: compute via a permutation representation.
##
#############################################################################

if IsMatrixGroup(G) then
    iso := IsomorphismPermGroup(G);

    if Length(arg) = 2 then
        R := ResolutionFiniteGroup(Image(iso, G), K);
    elif Length(arg) = 3 then
        R := ResolutionFiniteGroup(Image(iso, G), K, tietze);
    elif Length(arg) = 4 then
        R := ResolutionFiniteGroup(
            Image(iso, G), K, tietze, Charact
        );
    else
        R := ResolutionFiniteGroup(
            Image(iso, G), K, tietze, Charact, arg[5]
        );
    fi;

    R!.elts := List(R!.elts, x -> PreImageElm(iso, x));
    R!.group := G;

    return R;
fi;

#############################################################################
##
##  Elements and fast element -> index lookup.
##
##  Keep AllElts in its native enumeration order for use as the dictionary
##  domain, while Elts is reordered so that the identity is at position 1,
##  matching HAP's convention.
##
#############################################################################

AllElts := Elements(G);
Elts := ShallowCopy(AllElts);

if Elts[1] <> Identity(G) then
    i := Position(Elts, Identity(G));
    Elts[i] := Elts[1];
    Elts[1] := Identity(G);
fi;

RemoveSet(Gens, Identity(G));

ElementIndex := NewDictionary(false, true, AllElts);

for i in [1..N] do
    AddDictionary(ElementIndex, Elts[i], i);
od;

ExtendedElts :=
    List(Gens, g -> LookupDictionary(ElementIndex, g));

Append(ExtendedElts, [1..N]);

InverseIndex :=
    List(
        [1..N],
        g -> LookupDictionary(
            ElementIndex,
            Inverse(Elts[g])
        )
    );

#############################################################################
##
##  Algebraic reduction.
##
#############################################################################

if Charact = 0 then
    AlgebraicRed := AlgebraicReduction;
else
    AlgebraicRed := function(w)
        return AlgebraicReduction(w, Charact);
    end;
fi;

#############################################################################
##
##  Group multiplication lookup.
##
#############################################################################

if N <= 5060 then
    MT := MultiplicationTable(Elts);

    ActionIndex := function(g, h)
        return MT[g][h];
    end;
else
    MT := fail;

    ActionIndex := function(g, h)
        return LookupDictionary(
            ElementIndex,
            Elts[g] * Elts[h]
        );
    end;
fi;

#############################################################################
##
##  Sign helper.
##
#############################################################################

ChangeSign := function(j, b)
local r, x;

    if j > 0 then
        return b;
    fi;

    r := [];

    for x in b do
        Add(r, [-x[1], x[2]]);
    od;

    return r;
end;

#############################################################################
##
##  Resolution state.
##
#############################################################################

MaxComplex := [
    [ListWithIdenticalEntries(N, 0)]
];

MaxComplex[1][1][1] := 1;

ContractionMatrix := [];
ComputedContractions := [];

RemainingZeros := 0;

#############################################################################
##
##  Dimension and boundary.
##
#############################################################################

Dimension := function(i)

    if i < 0 then
        return 0;
    elif i = 0 then
        return 1;
    fi;

    return Length(PseudoBoundary[i]);
end;

PseudoBoundary := [];

Boundary := function(i, j)

    if i <= 0 then
        return [];
    fi;

    return ChangeSign(
        j,
        PseudoBoundary[i][AbsInt(j)]
    );
end;

#############################################################################
##
##  Differential-length based choice of the next unmatched cell.
##
##  This avoids building temporary lists of every zero position and every
##  corresponding differential length.
##
#############################################################################

FirstZero := function(MC, degree)
local
    j, g, len, best, bestlen;

    best := fail;
    bestlen := fail;

    for j in [1..Length(MC)] do
        for g in [1..N] do

            if MC[j][g] = 0 then

                len := DiffLengths[j][g];

                if len = 0 then
                    len := Length(
                        Differential(degree, [j, g])
                    );

                    DiffLengths[j][g] := len;
                fi;

                if len = 1 then
                    return [j, g];
                fi;

                if best = fail or len < bestlen then
                    best := [j, g];
                    bestlen := len;
                fi;

            fi;

        od;
    od;

    return best;
end;

#############################################################################
##
##  Contraction.
##
#############################################################################

Contraction := function(i, x)
local
    m, c, e, y, factor, ax1, z, cached;

    if i < 1 then
        return [[-1, 1]];
    fi;

    ax1 := AbsInt(x[1]);

    cached := ComputedContractions[i][ax1][x[2]];

    if cached <> 0 then
        return ChangeSign(x[1], cached);
    fi;

    z := [ax1, x[2]];

    if ContractionMatrix[i][z[1]][z[2]] = 1 then
        return [];
    fi;

    m := ContractionMatrix[i][z[1]][z[2]];

    c := [[
        -SignInt(x[1]) * m[1],
        m[2]
    ]];

    factor :=
        -SignInt(x[1]) * SignInt(m[1]);

    for e in PseudoBoundary[i][AbsInt(m[1])] do

        y := [
            factor * e[1],
            ActionIndex(m[2], e[2])
        ];

        if AbsInt(y[1]) <> ax1 or y[2] <> x[2] then
            Append(c, Contraction(i, y));
        fi;

    od;

    if i < K or Extendible then
        ComputedContractions[i][z[1]][z[2]] :=
            ChangeSign(x[1], c);
    fi;

    return c;
end;

Homotopy := function(i, p)

    if i < 0 then
        return fail;
    fi;

    return ChangeSign(-1, Contraction(i + 1, p));
end;

#############################################################################
##
##  Differential.
##
##  Process translated boundary cells directly rather than constructing a
##  complete translated copy first.
##
#############################################################################

Differential := function(i, p)
local
    j, k, Diff, e, x, s;

    j := p[1];
    k := p[2];

    Diff := [p];

    if i = 1 then
        Add(Diff, [-1, 1]);
        return Diff;
    fi;

    s := SignInt(j);

    for e in PseudoBoundary[i - 1][AbsInt(j)] do

        x := [
            s * e[1],
            ActionIndex(k, e[2])
        ];

        Append(
            Diff,
            Contraction(i - 1, x)
        );

    od;

    return Diff;
end;

#############################################################################
##
##  Contraction cache.
##
#############################################################################

InitComputedContractions := function(i)
local d;

    d := Dimension(i - 1);

    ComputedContractions[i] :=
        List(
            [1..d],
            x -> ListWithIdenticalEntries(N, 0)
        );
end;

#############################################################################
##
##  Incremental consequence engine.
##
##  For every translated i-cell [j,g], Residual[j][g] is the number of
##  boundary terms which still correspond to zero entries of MC.
##
##  IncidenceByFace[a] contains [j,t] whenever the boundary of the base
##  i-cell j contains a term [+-a,t].
##
##  When MC[a][h] changes 0 -> 1, [j,g] is affected exactly when
##
##      g * t = h,
##
##  hence
##
##      g = h * t^-1.
##
#############################################################################

InitConsequenceEngine := function(MC)
local a;

    IncidenceByFace := [];

    for a in [1..Length(MC)] do
        IncidenceByFace[a] := [];
    od;

    Residual := [];
    ConsequenceQueued := [];

    ConsequenceQueue := [];
    ConsequenceQueueHead := 1;
end;

EnqueueConsequence := function(j, g)

    if not ConsequenceQueued[j][g] then
        Add(ConsequenceQueue, [j, g]);
        ConsequenceQueued[j][g] := true;
    fi;

end;

#############################################################################
##
##  Return the unique currently-unmatched boundary face of [j,g].
##
#############################################################################

UniqueZeroFace := function(i, MC, j, g)
local e, h;

    for e in PseudoBoundary[i][j] do

        h := ActionIndex(g, e[2]);

        if MC[AbsInt(e[1])][h] = 0 then
            return [e[1], h];
        fi;

    od;

    return fail;
end;

#############################################################################
##
##  Record a translated sphere when Tietze reduction is enabled.
##
#############################################################################

AddSphere := function(i, j, g)

    if tietze then

        Add(
            Spheres,
            List(
                PseudoBoundary[i][j],
                e -> [
                    e[1],
                    ActionIndex(g, e[2])
                ]
            )
        );

    fi;

end;

#############################################################################
##
##  Mark one MC entry and update only affected translated cells.
##
#############################################################################

MarkMCEntry := function(i, MC, p)
local
    a, h,
    incidence,
    j, t, g,
    r;

    a := AbsInt(p[1]);
    h := p[2];

    if MC[a][h] = 1 then
        return;
    fi;

    MC[a][h] := 1;
    RemainingZeros := RemainingZeros - 1;

    for incidence in IncidenceByFace[a] do

        j := incidence[1];
        t := incidence[2];

        g := ActionIndex(
            h,
            InverseIndex[t]
        );

        Residual[j][g] := Residual[j][g] - 1;
        r := Residual[j][g];

        if r = 1 then
            EnqueueConsequence(j, g);
        elif r = 0 then
            AddSphere(i, j, g);
        fi;

    od;

end;

#############################################################################
##
##  Register one newly-created boundary row.
##
##  PseudoBoundary[i] grows during NextResTerm, so incidences are registered
##  incrementally.  The new row's translated residual counts are initialized
##  against the CURRENT MC state.
##
#############################################################################

RegisterBoundaryRow := function(i, MC, j)
local
    e,
    a, t,
    g, h,
    r;

    for e in PseudoBoundary[i][j] do

        a := AbsInt(e[1]);
        t := e[2];

        Add(
            IncidenceByFace[a],
            [j, t]
        );

    od;

    Residual[j] :=
        ListWithIdenticalEntries(N, 0);

    ConsequenceQueued[j] :=
        ListWithIdenticalEntries(N, false);

    for g in [1..N] do

        r := 0;

        for e in PseudoBoundary[i][j] do

            h := ActionIndex(g, e[2]);

            if MC[AbsInt(e[1])][h] = 0 then
                r := r + 1;
            fi;

        od;

        Residual[j][g] := r;

        if r = 1 then
            EnqueueConsequence(j, g);
        elif r = 0 then
            AddSphere(i, j, g);
        fi;

    od;

end;

#############################################################################
##
##  Drain the consequence queue.
##
##  Queue entries may become stale while waiting, so Residual is checked
##  immediately before processing.
##
#############################################################################

ProcessConsequences := function(i, MC)
local
    x,
    j, g,
    p,
    a, h;

    while ConsequenceQueueHead <= Length(ConsequenceQueue) do

        x := ConsequenceQueue[ConsequenceQueueHead];
        ConsequenceQueueHead :=
            ConsequenceQueueHead + 1;

        j := x[1];
        g := x[2];

        ConsequenceQueued[j][g] := false;

        if Residual[j][g] = 1 then

            p := UniqueZeroFace(i, MC, j, g);

            if p <> fail then

                a := AbsInt(p[1]);
                h := p[2];

                if i < K or Extendible then

                    ContractionMatrix[i][a][h] :=
                        [
                            SignInt(p[1]) * j,
                            g
                        ];

                fi;

                MaxComplex[i + 1][j][g] := 1;

                MarkMCEntry(i, MC, p);

            fi;

        fi;

    od;

    ConsequenceQueue := [];
    ConsequenceQueueHead := 1;

end;

#############################################################################
##
##  Compute the next resolution term.
##
#############################################################################

NextResTerm := function(i)
local
    ii,
    p,
    MC,
    l,
    Diff,
    j;

    PseudoBoundary[i] := [];
    MaxComplex[i + 1] := [];

    ContractionMatrix[i] :=
        ShallowCopy(MaxComplex[i]);

    MC :=
        StructuralCopy(MaxComplex[i]);

    Spheres := [];

    DiffLengths :=
        List(
            [1..Length(MC)],
            x -> ListWithIdenticalEntries(N, 0)
        );

    InitComputedContractions(i);

    l := ListWithIdenticalEntries(N, 0);
    l[1] := 1;

    RemainingZeros := 0;

    for ii in [1..Length(MC)] do
        RemainingZeros :=
            RemainingZeros +
            Number(
                MC[ii],
                x -> x = 0
            );
    od;

    InitConsequenceEngine(MC);

    while RemainingZeros > 0 do

        p := FirstZero(MC, i);

        if tietze then

            Diff :=
                TietzeReduction(
                    Spheres,
                    AlgebraicRed(
                        Differential(i, p)
                    )
                );

        else

            Diff :=
                AlgebraicRed(
                    Differential(i, p)
                );

        fi;

        Add(
            PseudoBoundary[i],
            Diff
        );

        Add(
            MaxComplex[i + 1],
            ShallowCopy(l)
        );

        j := Length(PseudoBoundary[i]);

        if i < K or Extendible then
            ContractionMatrix[i][p[1]][p[2]] :=
                [j, 1];
        fi;

        # Mark the explicitly selected cell.  This updates all previously
        # registered boundary rows.
        MarkMCEntry(
            i,
            MC,
            p
        );

        # Register the new boundary row against the current MC state.
        RegisterBoundaryRow(
            i,
            MC,
            j
        );

        # Propagate all newly available elementary collapses.
        ProcessConsequences(
            i,
            MC
        );

    od;

    DiffLengths := 0;
    Residual := 0;
    IncidenceByFace := 0;
    ConsequenceQueued := 0;
    ConsequenceQueue := 0;
    MC := 0;

end;

#############################################################################
##
##  Resolution object.
##
#############################################################################

R := Objectify(
    HapResolution,
    rec(
        dimension := Dimension,
        boundary := Boundary,
        homotopy := Homotopy,
        elts := Elts,
        group := G,
        vectorField := ContractionMatrix,

        properties := [
            ["length", 0],
            ["reduced", true],
            ["type", "resolution"],
            ["characteristic", Charact]
        ]
    )
);

#############################################################################
##
##  Extension.
##
#############################################################################

ExtendRes := function()
local i;

    i := R!.properties[1][2] + 1;

    NextResTerm(i);

    R!.properties[1][2] := i;

    MaxComplex[i] := [];

    if i > 1 and saveSpace then
        InitComputedContractions(i - 1);
    fi;
end;

if Extendible then
    R!.extend := ExtendRes;
fi;

for i in [1..K] do
    ExtendRes();
od;

return R;

end);
