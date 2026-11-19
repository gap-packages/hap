# Pull membership back to the original group, which may have a specialized
# membership method for an infinite matrix group.
InstallMethod(\in,
    [IsMultiplicativeElementWithInverse,
     IsGroup and HasNonFreeConjugationData],
    function(k, K)
        local data;
        data := NonFreeConjugationData(K);
        if not IsIdenticalObj(FamilyObj(k), FamilyObj(One(K))) then
            return false;
        fi;
        return k^data.conjugator in data.original;
    end);

InstallGlobalFunction(ConjugateNonFreeResolution, function(R, g)
    local inverse, conjugateGroup, elts, S;
    if not IsHapNonFreeResolution(R) then
        Error("Expected a HapNonFreeResolution");
    fi;
    inverse := g^-1;
    conjugateGroup := function(H)
        local K;
        K := Group(List(GeneratorsOfGroup(H), h -> h^inverse), One(H));
        SetNonFreeConjugationData(K, rec(original := H, conjugator := g));
        if HasIsFinite(H) then SetIsFinite(K, IsFinite(H)); fi;
        return K;
    end;

    # The source can append elements as it computes boundaries or homotopies.
    # Keep the indices shared and conjugate elements only when accessed.
    elts := ListToPseudoList([]);
    SetIsPseudoListWithFunction(elts, true);
    elts!.eltfun := i -> R!.elts[i]^inverse;
    elts!.lnthfun := function() return Length(R!.elts); end;
    elts!.posfun := k -> Position(R!.elts, k^g);
    elts!.addfun := function(k) Add(R!.elts, k^g); end;

    S := rec(
        group := conjugateGroup(R!.group),
        elts := elts,
        dimension := R!.dimension,
        boundary := R!.boundary,
        stabilizer := function(n, i)
            return conjugateGroup(R!.stabilizer(n, i));
        end,
        action := R!.action,
        properties := ShallowCopy(R!.properties)
    );
    if IsBound(R!.homotopy) then S.homotopy := R!.homotopy; fi;
    return Objectify(HapNonFreeResolution, S);
end);
