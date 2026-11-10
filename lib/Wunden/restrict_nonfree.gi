# Restrict a HAP non-free resolution to a finite-index subgroup.
# After LoadPackage("hap"); S := RestrictNonFreeResolution(R,K);

# Exact membership in K intersect t*H*t^-1, even when H is infinite.
InstallMethod(\in,
    [IsMultiplicativeElementWithInverse,
     IsGroup and HasRestrictionCellStabilizerData],
    function(g,h)
        local data;
        data := RestrictionCellStabilizerData(h);
        return g in data.subgroup and
            data.translate^-1*g*data.translate in data.original;
    end);

InstallGlobalFunction(RestrictNonFreeResolution, function(R,K)
    local G, one, length, top, T, cosetPosition, q, g, x,
          data, records, elts, originalIndex, subgroupIndex, prepare,
          dimension, checkCell, fromOriginal, toOriginal, boundary,
          stabilizer, action, homotopy, boundaryCache, RT;
    if not IsHapNonFreeResolution(R) or not IsGroup(K) then
        Error("Expected a HapNonFreeResolution and a subgroup");
    fi;
    G := R!.group; one := One(G); length := Length(R);
    if one <> One(K) then
        Error("R!.group and K must use the same element representation; ",
            "their identity elements differ. Transport K or R first.");
    fi;
    if not IsInt(length) or length < 0 then
        Error("The input must have a finite nonnegative recorded length");
    fi;
    # Enumerate LEFT cosets x*K using membership in K. Finite index and
    # terminating membership tests in K are preconditions.
    if HasIsFinite(G) and IsFinite(G) then
        if not IsSubgroup(G,K) then Error("K must be a subgroup of R!.group"); fi;
    fi;
#    T := [one];
#    ## cosetPosition := g -> PositionProperty(T,t -> t^-1*g in K);
#if IsBound(K!.ambientMembership) then                         # Originally 
#  cosetPosition:=function(g)                                  # changed to
#  local g1;                                                   # this faster
#  g1:=g^-1;                                                   # version before
#  return PositionProperty(T, t->K!.ambientMembership(g1*t));  # realizing that
#  end;                                                        # this is just a
#else                                                          # left 
#  cosetPosition:=function(g)                                  # transversal
#  local g1;                                                   #
#  g1:=g^-1;                                                   #
#  return PositionProperty(T, t->g1*t in K);                   #
#  end;                                                        #
#fi;                                                           #
#    q := 1;
#    while q <= Length(T) do
#        for g in GeneratorsOfGroup(G) do
#            x := g*T[q];
#            if cosetPosition(x)=fail then Add(T,x); fi;
#        od;
#        q := q+1;
#    od;

#The next three lines replace (and speed up) the above commented lines.
RT:=RightTransversal(G,K);
T:=List(RT,x->x^-1);
cosetPosition:=function(g); return PositionCanonical(RT,g^-1); end;

    top := length;
    while top > 0 and R!.dimension(top)=0 do top := top-1; od;
    data := List([0..top],n -> []);
    records := List([0..top],n -> []);
    boundaryCache := List([0..top],n -> []);
    elts := [one];

    originalIndex := function(g)
        local pos;
        #pos := Position(R!.elts,g); #Changed this tp focus on speed
        #if pos=fail then            #rather than on memory
        Add(R!.elts,g); 
        #pos := Length(R!.elts); fi; #
        return Length(R!.elts);
        #return pos;                 #
    end;
    subgroupIndex := function(g)
        local pos;
        #pos := Position(elts,g);    #Changed this to focus on speed
        #if pos=fail then            #rather than on memory
        Add(elts,g); 
        #pos := Length(elts); fi;    #
        return Length(elts);
        #return pos;                 #
    end;

    # Orbits of H on the LEFT cosets G/K. For a base coset x*K,
    # witness[q]*x*K = T[q]*K; the new cell is x^-1*e_i.
    prepare := function(n)
        local i, H, generators, owner, witness, base, orbit, cursor, q,
              next, s, t, z, h, loops, number, entry;
        if n<0 or n>top or not IsEmpty(data[n+1]) then return; fi;
        for i in [1..R!.dimension(n)] do
            H := R!.stabilizer(n,i); generators := GeneratorsOfGroup(H);
            owner := []; witness := [];
            for base in [1..Length(T)] do
                if not IsBound(owner[base]) then
                    number := Length(records[n+1])+1; t := T[base]^-1;
                    owner[base] := number; witness[base] := one;
                    orbit := [base]; cursor := 1; loops := [];
                    while cursor <= Length(orbit) do
                        q := orbit[cursor];
                        for s in generators do
                            next := cosetPosition(s*T[q]);
                            if not IsBound(owner[next]) then
                                owner[next] := number;
                                witness[next] := s*witness[q];
                                Add(orbit,next);
                            fi;
                            z := witness[next]^-1*s*witness[q];
                            if z<>one then AddSet(loops,t*z*t^-1); fi;
                        od;
                        cursor := cursor+1;
                    od;
                    # Schreier's lemma: these generate K intersect t*H*t^-1.
                    h := SubgroupNC(K,loops);
                    SetRestrictionCellStabilizerData(h,
                        rec(subgroup := K, original := H, translate := t));
                    if HasIsFinite(H) and not IsFinite(H) then
                        SetIsFinite(h,false);
                    fi;
                    entry := rec(originalCell := i, translate := t, stabilizer := h);
                    Add(records[n+1],entry);
                fi;
            od;
            data[n+1][i] := rec(owner := owner, witness := witness);
        od;
    end;
    dimension := function(n)
        if not IsInt(n) or n<0 or n>top then return 0; fi;
        prepare(n); return Length(records[n+1]);
    end;
    checkCell := function(n,k)
        if not IsInt(k) or k=0 or AbsInt(k)>dimension(n) then
            Error("Cell index out of range");
        fi;
    end;

    # Rewrite y*e_i as sign*h*(t*e_i), where y=h*t*s and s fixes e_i.
    # Only the stabilizer's orientation character is needed here.
    fromOriginal := function(n,word)
        local result, term, i, y, q, info, k, entry, s, h, sign;
        if IsEmpty(word) then return []; fi;
        if n<0 or n>top then Error("Degree out of range"); fi;
        prepare(n); result := [];
        for term in word do
            i := AbsInt(term[1]); y := R!.elts[term[2]];
            info := data[n+1][i]; q := cosetPosition(y^-1);
            k := info.owner[q]; entry := records[n+1][k];
            s := info.witness[q]^-1; h := y*s^-1*entry.translate^-1;
            sign := R!.action(n,i,originalIndex(s));
            Add(result,[SignInt(term[1])*sign*k,subgroupIndex(h)]);
        od;
        return result;
    end;
    toOriginal := function(n,word)
        local result, term, entry;
        result := [];
        for term in word do
            checkCell(n,term[1]); entry := records[n+1][AbsInt(term[1])];
            Add(result,[SignInt(term[1])*entry.originalCell,
                originalIndex(elts[term[2]]*entry.translate)]);
        od;
        return result;
    end;
    boundary := function(n,k)
        local entry, word;
        checkCell(n,k);
        if n=0 then return []; fi;
        if not IsBound(boundaryCache[n+1][AbsInt(k)]) then
            entry := records[n+1][AbsInt(k)];
            word := List(R!.boundary(n,entry.originalCell),term ->
                [term[1],originalIndex(entry.translate*R!.elts[term[2]])]);
            boundaryCache[n+1][AbsInt(k)] := fromOriginal(n-1,word);
        fi;
        return List(boundaryCache[n+1][AbsInt(k)],
            term -> [SignInt(k)*term[1],term[2]]);
    end;
    stabilizer := function(n,k)
        checkCell(n,k); return records[n+1][AbsInt(k)].stabilizer;
    end;
    action := function(n,k,j)
        local entry;
        checkCell(n,k); entry := records[n+1][AbsInt(k)];
        return R!.action(n,entry.originalCell,
            originalIndex(entry.translate^-1*elts[j]*entry.translate));
    end;
    homotopy := fail;
    if IsFunction(R!.homotopy) then
        homotopy := function(n,e)
            if n=-1 then return fromOriginal(0,R!.homotopy(-1,e)); fi;
            return fromOriginal(n+1,R!.homotopy(n,toOriginal(n,[e])[1]));
        end;
    fi;
    return Objectify(HapNonFreeResolution,rec(
        group := K, elts := elts, dimension := dimension, boundary := boundary,
        stabilizer := stabilizer, action := action, homotopy := homotopy,
        originalResolution := R, cosetRepresentatives := T,
        restrictionCells := records, toOriginal := toOriginal,
        fromOriginal := fromOriginal,
        properties := [["length",length],
            ["characteristic",EvaluateProperty(R,"characteristic")],
            ["reduced",false],["type","resolution"]]
    ));
end);
