#Replace the k-star of the e-th n-dimensional cell in L by a single cell
InstallGlobalFunction(Replace2, function(L,e,n,k)
local perm, i, j, d, l, topdim, matched, bndnk1, cnt, orbits, bool, StabU, sign, Elts, g, g1, lg, poslg, c, cell, pos,
newbnddup, newbnd, bnd, bndcell, h, instar, candidates, cellstar, incell, Boundary, 
Stab, Dim, action, oldpos, oldindex, RemapBoundary, oldi, signs, found,
kstar, ReduceTermOld, NormalizeBoundaryOld, matches,
coefficients, original, expanded, coefficient, rep, twist;

Elts := ShallowCopy(L!.elts);
kstar:=kStar(L,e,n,k);

pos := function(g)
local p;
    p := Position(Elts, g);
    if p = fail then
        Add(Elts, g);
        return Length(Elts);
    fi;
    return p;
end;

oldpos := function(g)
local p;
    p := Position(L!.elts, g);
    if p = fail then
        Add(L!.elts, g);
        return Length(L!.elts);
    fi;
    return p;
end;

#ki 09/07
## ReduceTermOld(degree,j,h): j is an OLD (pre-perm) cell index,
## h is an index into Elts. Returns [rep, eps] where rep is the smallest
## index such that Elts[h] and Elts[rep] lie in the same left coset of
## Stab(degree,j).
ReduceTermOld := function(degree, j, h)
local stab, p, t;
    stab := L!.stabilizer(degree, j);
    for p in [1..Length(Elts)] do
        t := Inverse(Elts[h]) * Elts[p];
        if t in stab then
            return [ p, L!.action(degree, j, oldpos(t)) ];
        fi;
    od;
    Error("ReduceTermOld: no representative found (should be unreachable)");
end;

## NormalizeBoundaryOld(bnd,degree): takes a raw list of [±j,h] terms (old
## indices), groups them by (j, canonical left-coset rep), sums signed
## multiplicities (folding in the eps twist), drops classes that cancel
## to zero, and expands surviving classes back into repeated [±j,rep]
## entries. This is what correctly cancels a coset that shows up twice
## with opposite signs -- which the old same-signed-value dedup could
## not do.
NormalizeBoundaryOld := function(bnd, degree)
local out, term, j, sgn, h, red, rep, eps, entry, matched, result, cnt;
    out := [];
    for term in bnd do
        j   := AbsoluteValue(term[1]);
        sgn := SignInt(term[1]);
        h   := term[2];
        red := ReduceTermOld(degree, j, h);
        rep := red[1];
        eps := red[2];
        matched := false;
        for entry in out do
            if entry[1] = j and entry[2] = rep then
                entry[3] := entry[3] + sgn*eps;
                matched := true;
                break;
            fi;
        od;
        if not matched then
            Add(out, [ j, rep, sgn*eps ]);
        fi;
    od;
    out := Filtered(out, x -> x[3] <> 0);
    result := [];
    for entry in out do
        for cnt in [1..AbsoluteValue(entry[3])] do
            Add(result, [ SignInt(entry[3]) * entry[1], entry[2] ]);
        od;
    od;
    return result;
end;
#09/07 ging bis hier

perm:=[];  ##perm[i][j] will be the new position of
           ##the old j-th cell of dimension n+i-1;
for i in [1..k+1] do
    if kstar[i]=[] then break;
    fi;
    perm[i]:=[];
    cnt:=0;
    orbits:=List(kstar[i],x->x[1]);
    for j in [1..L!.dimension(n+i-1)] do
        if j in orbits then cnt:=cnt+1;
        else perm[i][j]:=j-cnt;
        fi;
    od;
od;
#nrremovedtopcells=cnt
topdim:=L!.dimension(n+k)-cnt+1;

#ki #???Fehler: d hat keine Funktion
oldindex := function(p, i, d)
local q;
    q := Position(p, i);
    if q = fail then
        Error("No old index corresponding to new index ",
              i, " in dimension ", d);
    fi;
    return q;
end;

signs:=SimplifyLattice(L,kstar,n)[2];

newbnddup:=[];
#signs:=SimplifyLattice(L,kstar,n)[2];
for c in [1..Length(kstar[k+1])] do
    cell:=kstar[k+1][c];
    bndcell:=L!.boundary(n+k,cell[1]);
    for bnd in bndcell do
        h:=cell[2]*Elts[bnd[2]];
        instar:=false;
        candidates:=Filtered(kstar[k],y->y[1]=AbsoluteValue(bnd[1]));
        for cellstar in candidates do
            if Inverse(h)*cellstar[2] in cellstar[3] then #changed from h*Inverse(cellstar[2]) to the way it is now
                instar:=true;
                break;
            fi;
        od;
        if not instar then
            #sign:=L!.action(n+k,cell[1],pos(cell[2]));
            #sign:=SignInt(SimplifyLattice(L,kstar,n)[2][c]);
            Add(newbnddup,[signs[c]*bnd[1],pos(h)]);
        fi;
    od;
od;
#Print(newbnddup, "\n");

# newbnd:=[];
# for cell in newbnddup do
#     bool:=true;
#     for incell in Filtered(newbnd,y->y[1]=cell[1]) do
#         if Inverse(Elts[incell[2]])*Elts[cell[2]] in L!.stabilizer(n+k-1,AbsoluteValue(cell[1])) then
#         # changed from if Elts[incell[2]]*Inverse(Elts[cell[2]]) in L!.stabilizer(n+k-1,AbsoluteValue(cell[1])) then
#             bool:=false;
#             break;
#         fi;
#     od;
#     if bool then
#         Add(newbnd,cell);
#     fi;
# od;
newbnd := NormalizeBoundaryOld(newbnddup,n+k-1); #ki


newbnd:=List(newbnd, 
    cell -> [SignInt(cell[1]) * perm[k][AbsoluteValue(cell[1])],cell[2]]);

#ki
RemapBoundary := function(b, p, d, i)
local term, out, newidx;
    out := [];
    for term in b do
        newidx := p[AbsoluteValue(term[1])];
        if newidx = fail then
            Error("Boundary remap failed at d=", d, " i=", i,
                  " old face=", AbsoluteValue(term[1]));
        fi;
        Add(out, [ SignInt(term[1]) * newidx, term[2] ]);
    od;
    return out;
end;


Boundary := function(d,i)
local oldi, b;

    if d<n or d>n+k+1 then
        return L!.boundary(d,i);

    elif d=n then
        # boundaries in degree n go to degree n-1, which was not changed
        oldi := oldindex(perm[d-n+1], i, d);
        return L!.boundary(d, oldi);

    elif d>n and d<n+k then
        oldi := oldindex(perm[d-n+1], i, d);
        b := L!.boundary(d, oldi);
        return RemapBoundary(b, perm[d-n], d, i);

    elif d=n+k then
        if AbsoluteValue(i)<topdim then
            oldi := oldindex(perm[d-n+1], AbsoluteValue(i), d);
            b := L!.boundary(d, oldi);
            b := RemapBoundary(b, perm[d-n], d, AbsoluteValue(i));
            if i<0 then
                b := List(b, t -> [-t[1],t[2]]);
            fi;
            return b;
        elif AbsoluteValue(i)=topdim then
            if i<0 then
                return List(newbnd, t -> [-t[1],t[2]]);
            fi;
            return newbnd;
        fi;

    elif d=n+k+1 then
        if i<0 then
            return List(bndnk1[-i], t -> [-t[1],t[2]]);
        fi;
        return bndnk1[i];
    fi;
end;

StabU:=Stabstar(L,e,n,k);
Stab:=function(d,i)
    if d<n or d>n+k+1 then
        return L!.stabilizer(d,i);
    elif d>=n and d<n+k then
        return L!.stabilizer(d,Position(perm[d-n+1],AbsoluteValue(i)));
    elif d=n+k then
        if AbsoluteValue(i)<topdim then
            return L!.stabilizer(d,Position(perm[d-n+1],AbsoluteValue(i)));
        elif AbsoluteValue(i)=topdim then
            return StabU;
        fi;
    elif d=n+k+1 then
        return L!.stabilizer(d,i);
    fi;
end; 



Dim:=function(d)
local p;
    if d<n or d>=n+k+1 then
        return L!.dimension(d);
    elif d>=n and d<n+k then
        p := Compacted(perm[d-n+1]);
        if Length(p)=0 then
            return 0;
        fi;
        return Maximum(p);
    elif d=n+k then
        p := Compacted(perm[d-n+1]);
        if Length(p)=0 then
            return 1;
        fi;
        return Maximum(p)+1;
    fi;
end;

#changes start here



action := function(d,i,j)
local g, h, representative, cell, cell2;
    g := Elts[j];

    if d = n+k and AbsoluteValue(i) = Dim(d) then
        if not (g in StabU) then
            # FreeGResolution also asks for signs outside the stabilizer.
            # Use the same coset decomposition as HAP's ExtendScalars:
            # g = representative * h, and evaluate the character on h.
            representative:=CanonicalRightCountableCosetElement(StabU,Inverse(g))^-1;
            g:=Inverse(representative)*g;
        fi;
        # Read the orientation from a constituent top cell. This also works
        # when the new cell's cellular boundary is zero.
        cell := kstar[k+1][1];
        for cell2 in kstar[k+1] do
            if cell2[1] = cell[1] then
                h := Inverse(cell2[2])*g*cell[2];
                if h in cell[3] then
                    return signs[1]*signs[Position(kstar[k+1],cell2)]*
                           L!.action(n+k,cell[1],oldpos(h));
                fi;
            fi;
        od;
        Error("Replace2: star stabilizer does not preserve its top cells");

    elif d >= n and d <= n+k then
        return L!.action(d, oldindex(perm[d-n+1], AbsoluteValue(i), d), oldpos(g));

    else
        return L!.action(d, AbsoluteValue(i), oldpos(g));
    fi;
end;


# An attaching boundary contains the entire oriented chain of each star
# translate, not one copy of the new cell per constituent face. Collect
# coefficients modulo the star stabilizer (including its orientation),
# then divide by the number of constituent top cells. Preserve genuine
# multiplicities, and verify by expanding back into the old cell module.
# This must follow the definition of action, used for orientation twists.
orbits:=List(kstar[k+1],x->x[1]);
bndnk1:=[];
for j in [1..L!.dimension(n+k+1)] do
    bndnk1[j]:=[];
    coefficients:=[];
    original:=[];
    for cell in L!.boundary(n+k+1,j) do
        if AbsoluteValue(cell[1]) in orbits then
            Add(original,cell);
            c:=Position(orbits,AbsoluteValue(cell[1]));
            lg:=L!.elts[cell[2]]*Inverse(kstar[k+1][c][2]);
            poslg:=pos(lg);
            rep:=First([1..poslg],p->Inverse(lg)*Elts[p] in StabU);
            twist:=action(n+k,topdim,pos(Inverse(lg)*Elts[rep]));
            if not IsBound(coefficients[rep]) then coefficients[rep]:=0; fi;
            coefficients[rep]:=coefficients[rep]+SignInt(cell[1])*signs[c]*twist;
        else
            Add(bndnk1[j],[SignInt(cell[1])*perm[k+1][AbsoluteValue(cell[1])],cell[2]]);
        fi;
    od;
    expanded:=[];
    for rep in [1..Length(coefficients)] do
        if IsBound(coefficients[rep]) and coefficients[rep]<>0 then
            coefficient:=coefficients[rep]/Length(kstar[k+1]);
            if not IsInt(coefficient) then
                Error("Replace2: attaching boundary is not a whole star chain in dimension ",
                      n+k+1," cell ",j);
            fi;
            for cnt in [1..AbsoluteValue(coefficient)] do
                Add(bndnk1[j],[SignInt(coefficient)*topdim,rep]);
                for c in [1..Length(kstar[k+1])] do
                    Add(expanded,[SignInt(coefficient)*signs[c]*kstar[k+1][c][1],
                                  pos(Elts[rep]*kstar[k+1][c][2])]);
                od;
            od;
        fi;
    od;
    if SortedList(NormalizeBoundaryOld(original,n+k)) <>
       SortedList(NormalizeBoundaryOld(expanded,n+k)) then
        Error("Replace2: attaching boundary does not expand to the original star chain in dimension ",
              n+k+1," cell ",j);
    fi;
od;


return Objectify(HapNonFreeResolution,
        rec(
            dimension := Dim,
            boundary := Boundary,
            stabilizer := Stab,
            action := action,
            elts := Elts,
            group := L!.group,
            homotopy := fail,
            properties := [
                ["characteristic",0],
                ["type","resolution"],
                ["length",1000]
            ]
        )
    );
end);
###############################

