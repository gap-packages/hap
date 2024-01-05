
# RegularCWMapToCWSubcomplex
################################################################################
############ Input: a strictly cellular map f: Y -> X of regular ###############
################### CW-complexes ###############################################
################################################################################
########### Output: the CW-subcomplex [X,S] corresponding to the ###############
################### image of f #################################################
################################################################################
InstallGlobalFunction(
    RegularCWMapToCWSubcomplex,
    function(f)
    local
        dim, src_bnd, S, i, j;

    dim:=EvaluateProperty(f!.source,"dimension")*1;
    src_bnd:=f!.source!.boundaries*1;
    S:=List([0..dim],x->[]);

    for i in [1..dim+1] do
        for j in [1..Length(src_bnd[i])] do
            Add(S[i],f!.mapping(i-1,j));
        od;
    od;
    
    return [ShallowCopy(f!.target),S];
end);

# CWSubcomplexToRegularCWMap
################################################################################
############ Input: a CW-subcomplex [X,S] ######################################
################################################################################
########### Output: the corresponding inclusion map f: Y -> X ##################
################################################################################
InstallGlobalFunction(
    CWSubcomplexToRegularCWMap,
    function(YS)
    local
        map, src, i, j, trg_cell;

    map:={i,j}->YS[2][i+1][j];

    src:=List([1..Length(Filtered(YS[2],y->y<>[]))+1],x->[]);
    src[1]:=List(YS[2][1],x->[1,0]);

    for i in [2..Length(src)-1] do
        for j in [1..Length(YS[2][i])] do
            trg_cell:=YS[1]!.boundaries[i][YS[2][i][j]]*1;
            trg_cell:=trg_cell{[2..trg_cell[1]+1]};
            trg_cell:=List(trg_cell,x->Position(YS[2][i-1],x));
            Add(trg_cell,Length(trg_cell),1);
            Add(src[i],trg_cell*1);
        od;
    od;

    return Objectify(
        HapRegularCWMap,
        rec(
            source:=RegularCWComplex(src),
            target:=ShallowCopy(YS[1]),
            mapping:=map
        )
    );
end);

# IntersectionCWSubcomplex
################################################################################
############ Input: two CW-subcomplexes [X,S], [X,S'] ##########################
################################################################################
########### Output: the CW-subcomplex [X,S''] corresponding to #################
################### their intersection #########################################
################################################################################
InstallGlobalFunction(
    IntersectionCWSubcomplex,
    function(XS_1,XS_2)
    local
        max, xs_1, xs_2, i, j;
    
    max:=Maximum(Length(XS_1[2]),Length(XS_2[2]));
    xs_1:=ShallowCopy(XS_1[2]);
    xs_2:=ShallowCopy(XS_2[2]);

    for i in [xs_1,xs_2] do
        if Length(i)<max then
            for j in [Length(i)+1..max] do
                Add(i,[]);
            od;
        fi;
    od;

    return[
        ShallowCopy(XS_1[1]),
        List(
            [1..max],
            x->Intersection(xs_1[x],xs_2[x])
        )
    ];
end);

# PathComponentsCWSubcomplex
################################################################################
############ Input: a CW-subcomplex [X,S] ######################################
################################################################################
########### Output: a list of CW-subcomplexes [[X,S1],...,[X,Sn]] ##############
################### arising as the path components of [X,S] ####################
################################################################################
InstallGlobalFunction(
    PathComponentsCWSubcomplex,
    function(XS)
    local   
        ccs, i, j, k, cell, int, l;

    ccs:=List(
        XS[2][1]*1,
        x->Concatenation(
            [[x]],
            List([2..Length(XS[2])],y->[])
        )
    );
    for i in [1..Length(ccs)] do
        for j in [2..Length(ccs[i])] do
            for k in [1..Length(XS[2][j])] do
                cell:=XS[1]!.boundaries[j][XS[2][j][k]]*1;
                cell:=cell{[2..cell[1]+1]};
                int:=Intersection(cell,ccs[i][j-1]);
                if int<>[] then
                    for l in [1..Length(cell)] do
                        Add(ccs[i][j-1],cell[l]);
                    od;
                    Add(ccs[i][j],XS[2][j][k]);
                fi;
            od;
        od;
    od;
    for i in [1..Length(ccs)] do
        for j in Filtered([1..Length(ccs)],x->x<>i) do
            if ccs[i]<>[] and ccs[j]<>[] then
                for k in [1..Length(ccs[i])] do
                    if ccs[j]<>[] then
                        if Intersection(ccs[i][k],ccs[j][k])<>[] then
                            for l in [1..Length(ccs[i])] do
                                ccs[i][l]:=Set(
                                    Concatenation(
                                        ccs[i][l],
                                        ccs[j][l]
                                    )
                                );
                            od;
                            ccs[j]:=[];
                        fi;
                    fi;
                od;
            fi;
        od;
    od;
    ccs:=Filtered(List(ccs,x->List(x,y->Set(y))),z->z<>[]);

    return List(ccs,x->[ShallowCopy(XS[1]),x]);
end);

# ClosureCWCell
################################################################################
############ Input: a CW-complex Y and two integers k>=0 & i>1 #################
################################################################################
########### Output: a CW-subcomplex [Y,S] corresponding to the #################
################### topological closure of the ith k-cell of Y #################
################################################################################
InstallGlobalFunction(
    ClosureCWCell,
    function(Y,k,i)
    local
        complex, clsr, l, m, bnd, n;

    complex:=Y!.boundaries*1;

    clsr:=List([1..k+1],x->[]);
    Add(clsr[k+1],i);
    
    for l in Reversed([2..k+1]) do
        for m in [1..Length(clsr[l])] do
            bnd:=[];
            for n in [2..Length(complex[l][clsr[l][m]])] do
                Add(bnd,complex[l][clsr[l][m]][n]);
            od;
            clsr[l-1]:=Union(clsr[l-1],bnd);
        od;
    od;

    return [RegularCWComplex(complex),clsr];
end);

# HAP_KK_AddCell
################################################################################
############ Input: a cell complex B, an integer k>=0 and lists of #############
################### positive integers b and c ##################################
################################################################################
########### Output: the cell complex B with an added k-cell whose ##############
################### boundary is specified by b and whose coboundary ############
################### is specified by c ##########################################
################################################################################
InstallGlobalFunction(
    HAP_KK_AddCell,
    function(B,k,b,c)
    local
        i;

    Add(b,Length(b),1);
    Add(B[k+1],b);
    
    for i in [1..Length(c)] do
        Add(B[k+2][c[i]],Length(B[k+1]));
        B[k+2][c[i]][1]:=B[k+2][c[i]][1]+1;
    od;
end);

# BarycentricallySubdivideCell
################################################################################
############ Input: an inclusion of cell complexes f:Z->Y and two integers #####
################### k>0 & n>=0 #################################################
################################################################################
########### Output: the regular CW-map f':Z'->Y' corresponding to the same #####
################### spaces but where the kth n-cell of Y has been ##############
################### barycentrically subdivided #################################
################################################################################
InstallGlobalFunction(
    BarycentricallySubdivideCell,
    function(f,n,k)

    local
        inc, Y, clsr, complex,
        Subdivide, i, j;

    inc:=RegularCWMapToCWSubcomplex(ShallowCopy(f));
    Y:=inc[1];

    clsr:=ClosureCWCell(Y,n,k);
    complex:=Y!.boundaries*1;
    clsr:=clsr[2]*1;

    Subdivide:=function(a,b)
    # subdivides the bth a-cell into as many a-cells
    # as there are in its barycentric subdivision
        local
            sub_clsr, i, j, l, ints, bnd_ints,
            bnd, pre_len, m, pre_len_bnd;

        sub_clsr:=ClosureCWCell(RegularCWComplex(complex),a,b)[2];

        Add(complex[1],[1,0]); # the barycentre
        if Length(inc[2])>=a+1 then
            if b in inc[2][a+1] then
                Add(inc[2][1],Length(complex[1]));
            fi;
        fi;

        pre_len_bnd:=List([1..a],x->Length(complex[x]));

        for i in [1..a] do
            pre_len:=Length(complex[i+1])*1;
            for j in [1..Length(sub_clsr[i])] do
                if i=1 then # edge case is dealt with separately
                    if a=1 and j=1 then # we overwrite the bth a-cell with
                    # the first a-cell of its barycentric subdivision
                        complex[2][b]:=[
                            2,
                            sub_clsr[1][1],
                            Length(complex[1])
                        ];
                    else
                        Add(
                            complex[2],
                            [
                                2,
                                sub_clsr[1][j],
                                Length(complex[1])
                            ]
                        );
                        if Length(inc[2])>=a+1 then
                            if b in inc[2][a+1] then
                                Add(inc[2][2],Length(complex[2]));
                            fi;
                        fi;
                    fi;
                else
                    bnd:=[];
                    ints:=List([pre_len_bnd[i]+1..Length(complex[i])]);
                    bnd_ints:=List(
                        ints,
                        x->complex[i][x]{[2..complex[i][x][1]+1]}
                    );
                    bnd_ints:=List(
                        bnd_ints,
                        x->Intersection(
                            x,
                            complex[i][sub_clsr[i][j]]{
                                [2..complex[i][sub_clsr[i][j]][1]+1]
                                }
                        )<>[]
                    );
                    bnd:=Concatenation(
                        [sub_clsr[i][j]],
                        Filtered(ints,x->bnd_ints[Position(ints,x)]=true)
                    );
                    Add(bnd,Length(bnd),1);
                    if i=a and j=1 then # overwrite as before
                        complex[a+1][b]:=bnd;
                    else
                        Add(complex[i+1],bnd);

                        if Length(inc[2])>=a+1 then
                            if b in inc[2][a+1] then
                                Add(inc[2][i+1],Length(complex[i+1]));
                            fi;
                        fi;
                    fi;
                fi;

                if i=a and a<Dimension(Y) then # rewrite the coboundary of
                # the replaced a-cell to contain all the a-cells in its
                # barycentric subdivision
                    for l in [2..Length(Y!.coboundaries[a+1][b])] do
                        for m in [pre_len+1..Length(complex[a+1])] do
                            Add(
                                complex[a+2][Y!.coboundaries[a+1][b][l]],
                                m
                            );
                        od;
                        complex[a+2][Y!.coboundaries[a+1][b][l]][1]:=Length(
                            complex[a+2][Y!.coboundaries[a+1][b][l]]
                        )-1;
                    od;
                fi;
            od;
        od;  
    end;

    for i in [2..Length(clsr)] do # inductively apply Subdivide to obtain the
        for j in [1..Length(clsr[i])] do # barycentric subdivision of the
            Subdivide(i-1,clsr[i][j]); # kth n-cell of Y as required
        od;
    od;

    return CWSubcomplexToRegularCWMap(
        [
            RegularCWComplex(complex),
            inc[2]
        ]
    );
end);

# SubdivideCell
################################################################################
############ Input: an inclusion of cell complexes f:Z->Y and two integers #####
################### k>0 & n>=0 #################################################
################################################################################
########### Output: the regular CW-map f':Z'->Y' corresponding to the same #####
################### spaces but where the kth n-cell of Y has been ##############
################### subdivided into as many n-cells as there are (n-1)-cells ###
################### in the boundary of that cell ###############################
################################################################################
InstallGlobalFunction(
    SubdivideCell,
    function(f,n,k)
    local
        sub, Y, closure, plus1,
        i, j, bnd, Last;
Last:=function(L); return L[Length(L)];end;

    sub:=RegularCWMapToCWSubcomplex(ShallowCopy(f));
    Y:=sub[1]; # the actual CW-complex
    sub:=sub[2]*1; # the indexing of the subcomplex
    closure:=ClosureCWCell(Y,n,k)[2];
    Y:=Y!.boundaries*1;

    Add(Y[1],[1,0]); # the barycentre of the kth n-cell
    plus1:=List([1..Length(closure)],x->[]);
    # this will associate an x-cell in the closure to
    # the resulting (x+1)-cell in the subdivision

    for i in [1..Length(closure)-1] do
        for j in [1..Length(closure[i])] do
            if i=1 then
                Add(Y[2],[2,closure[i][j],Length(Y[1])]);
                Add(plus1[1],Length(Y[2]));
            else
                bnd:=Y[i][closure[i][j]];
                bnd:=bnd{[2..bnd[1]+1]};
                bnd:=List(bnd,x->plus1[i-1][Position(closure[i-1],x)]);
                Add(bnd,closure[i][j],1);
                Add(bnd,Length(bnd),1);
                if i=Length(closure)-1 and j=1 then
                    Y[i+1][closure[i+1][1]]:=bnd;
                    Add(plus1[i],closure[i+1][1]);
                else
                    Add(Y[i+1],bnd);
                    Add(plus1[i],Length(Y[i+1]));
                fi;
            fi;
        od;
    od;
    for i in [1..Length(Y[Length(closure)+1])] do
        if Last(closure)[1] in
            Y[Length(closure)+1][i]{[2..Y[Length(closure)+1][i][1]+1]} then
                Append(Y[Length(closure)+1][i],plus1[Length(closure)]);
                Unbind(Y[Length(closure)+1][i][1]);
                Y[Length(closure)+1][i]:=Set(Y[Length(closure)+1][i]);
                Add(Y[Length(closure)+1][i],Length(Y[Length(closure)+1][i]),1);
        fi;
    od;
    return CWSubcomplexToRegularCWMap(
        [
            RegularCWComplex(Y),
            sub
        ]
    );
end);

# RegularCWComplexComplement
################################################################################
############ Input: an inclusion f: Z -> Y of regular CW-complexes #############
################################################################################
########### Output: let N(Z) denote an open tubular neighbourhood of ###########
################### Z in Y. this will output an inclusion B -> C where C #######
################### is homeomorphic to Y \ N(Z) and B is homeomorphic to the ###
################### boundary of Y \ N(Z) #######################################
################################################################################
InstallGlobalFunction(
    RegularCWComplexComplement,
    function(arg...)
    local
        f, check, subdiv, details, Y, B, IsInternal, count, total,
        path_comp, cobound_subcomplex, cbnd, i, j, clsr, int, crit,
        bary, IsSubpathComponent, ext_cell_2_f_notation,
        f_notation_2_ext_cell, ext_cells, k, ext_cell_bnd, e_n_bar,
        ext_cell_cbnd;

    if Length(arg)=1 then
    	arg:=[arg[1],"some","basic",false];
    fi;
    f:=ShallowCopy(arg[1]);
    check:=arg[2]; # which cells we check for the contractible closures:
    # "some" or "all"
    subdiv:=arg[3]; # the subdivision we use:
    # "basic", "barycentric" or "none"
    details:=arg[4]; # do/don't suppress progress bars and other output
    
    Y:=RegularCWMapToCWSubcomplex(f);
    for i in [1..Dimension(Y[1])-Length(Y[2])+2] do
        Add(Y[2],[]);
    od;
    B:=List([1..Dimension(Y[1])+1],x->[]);

    IsInternal:=function(n,k) # is the kth n-cell of Y in Y\Z?
        if n+1>Length(Y[2]) then
            return true;
        elif k in Y[2][n+1] then
            return false;
        fi;
        return true;
    end;

    if details then
        Print("Testing contractibility...\n");
    fi;
    bary:=[];
    path_comp:=List([1..Length(Y[1]!.boundaries)-1],x->[]);

    if check="some" then
    # we check only those cells that are `close' to the subcomplex
    # i.e. those cells lying within the coboundaries of the coboundaries
    # of (...) the cells of the subcomplex
        cobound_subcomplex:=List(
            [1..Length(Y[2])],
            x->List(
                [1..Length(Y[2][x])],
                y->Y[1]!.coboundaries[x][Y[2][x][y]]{
                    [2..Y[1]!.coboundaries[x][Y[2][x][y]][1]+1]
                }
            )
        );
        cobound_subcomplex:=List(cobound_subcomplex,Concatenation);
        cobound_subcomplex:=List(cobound_subcomplex,Set);
        cobound_subcomplex:=List(
            [1..Length(cobound_subcomplex)],
            x->Filtered(cobound_subcomplex[x],y->not y in Y[2][x+1])
        );
        Add(cobound_subcomplex,[],1); # keep indexing w/out 0-cells
        for i in [1..Length(cobound_subcomplex)-1] do
            for j in [1..Length(cobound_subcomplex[i])] do
                cbnd:=Y[1]!.coboundaries[i][cobound_subcomplex[i][j]];
                for k in [2..cbnd[1]+1] do
                    Add(cobound_subcomplex[i+1],cbnd[k]);
                od;
            od;
        od;
        cobound_subcomplex:=List(cobound_subcomplex,Set);
        for i in [1..Length(cobound_subcomplex)-1] do   
            cobound_subcomplex[i]:=Filtered(
                cobound_subcomplex[i],
                x->not x in Y[2][i+1]
            );
        od;
    fi;

    count:=0;
    total:=Sum(List(Y[1]!.boundaries{[1..Length(Y[1]!.boundaries)]},Length));

    # list of all of the path components of the intersection between
    # the closure of each internal cell and the subcomplex Z < Y
    # note: entries may be an empty list if they correspond to an
    #       empty intersection or they may not be assigned a
    #       value at all if the associated cell lies in Z
    for i in [1..Length(Y[1]!.boundaries)] do
        for j in [1..Length(Y[1]!.boundaries[i])] do
            count:=count+1;
            if IsInternal(i-1,j) then
                Add(B[i],Y[1]!.boundaries[i][j]);
                # output will contain all internal cells
                if check="some" then
                    if j in cobound_subcomplex[i] then
                        clsr:=ClosureCWCell(Y[1],i-1,j);
                        int:=IntersectionCWSubcomplex(clsr,Y);
                        path_comp[i][j]:=PathComponentsCWSubcomplex(int);
                    else
                        path_comp[i][j]:=[];
                    fi;
                else
                    clsr:=ClosureCWCell(Y[1],i-1,j);
                    int:=IntersectionCWSubcomplex(clsr,Y);
                    path_comp[i][j]:=PathComponentsCWSubcomplex(int);
                fi; 
                # CONTRACTIBILITY TEST
                # must be a non-empty subcomplex of dimension > 0
                if details then
                    Print(count," out of ",total," cells tested.","\r");
                fi;
                if subdiv<>"none" then
                    # test the contractibility of each path component
                    # if we find anything non-contractible,
                    # subdivide the problematic cells and restart
                    for k in [1..Length(path_comp[i][j])] do
                        if path_comp[i][j][k][2]<>List(path_comp[i][j][k][2],x->[])
                        and path_comp[i][j][k][2][2]<>[] then
                            crit:=CWSubcomplexToRegularCWMap(path_comp[i][j][k]);
                            crit:=Source(crit);
                            crit:=CriticalCellsOfRegularCWComplex(crit);
                            if not Length(crit)=1 then
                                Add(bary,[i-1,j]);
                            fi;
                        fi;
                    od;
                fi;
            else
                Add(B[i],"*"); # temporary entry to keep correct indexing
            fi;
        od;
    od;
    bary:=Set(bary);

    if bary=[] then
        if details then
            Print("\nThe input is compatible with this algorithm.\n");
        fi;
    else
        if details then
            Print("\nSubdividing ",Length(bary)," cell(s):\n");
        fi;
        for i in [1..Length(bary)] do
            if subdiv="basic" then
                f:=SubdivideCell(
                    f,
                    bary[i][1],
                    bary[i][2]
                );
            elif subdiv="barycentric" then
                f:=BarycentricallySubdivideCell(
                    f,
                    bary[i][1],
                    bary[i][2]
                );
            fi;
            if details then
                Print(Int(100*i/Length(bary)),"\% complete. \r");
            fi;
        od;
        Print("\n");
        return RegularCWComplexComplement(f,check,subdiv,details); # really bad
    fi;

    for i in [1..Length(B)] do
        for j in [1..Length(B[i])] do
            if i>1 and B[i][j]<>"*" then
                B[i][j]:=Filtered(
                        B[i][j]{[2..B[i][j][1]+1]},
                        x->"*"<>B[i-1][x]
                    );
                Add(B[i][j],Length(B[i][j]),1);
            fi;
        od;
    od;
    # at this point, B corresponds to the cell complex Y\Z

    IsSubpathComponent:=function(super,sub)
        local
            i;

        for i in [1..Length(sub[2])-1] do
            if not IsSubset(super[2][i],sub[2][i]) then
                return false;
            fi;
        od;

        return true;
    end;

    ext_cell_2_f_notation:=NewDictionary([],true);
    f_notation_2_ext_cell:=NewDictionary([],true);
    # takes a pair [n,k] corresponding to the kth n-cell (external) of B and
    # returns its associated "f notation" i.e. a triple [n+1,k',A] corresponding
    # to the k'th (n+1)-cell (internal) whose closure intersected with Z in the
    # path component A

    ext_cells:=List([1..Length(B)],x->[]);

    for i in [2..Length(path_comp)] do
        for j in [1..Length(path_comp[i])] do
            if IsBound(path_comp[i][j]) then
                if path_comp[i][j]<>[] then
                    # we add as many external (i-1)-cells to B as
                    # there are path components in path_comp[i][j]
                    for k in [1..Length(path_comp[i][j])] do
                        if i=2 then
                            ext_cell_bnd:=[0];
                        else
                            e_n_bar:=ClosureCWCell(Y[1],i-1,j)[2];
                            ext_cell_bnd:=List(
                                ext_cells[i-2],
                                x->LookupDictionary(
                                    ext_cell_2_f_notation,
                                    [i-3,x]
                                )
                            );
                            ext_cell_bnd:=Filtered(
                                ext_cell_bnd,
                                x->
                                    x[2] in e_n_bar[i-1]
                                    and
                                    IsSubpathComponent(path_comp[i][j][k],x[3])
                            );
                            ext_cell_bnd:=List(
                                ext_cell_bnd,
                                x->LookupDictionary(
                                    f_notation_2_ext_cell,
                                    x
                                )[2]
                            );
                        fi;
                        ext_cell_cbnd:=[j];

                        HAP_KK_AddCell(
                            B,
                            i-2,
                            ext_cell_bnd,
                            ext_cell_cbnd
                        );
                        Add(ext_cells[i-1],Length(B[i-1]));
                        AddDictionary(
                            ext_cell_2_f_notation,
                            [i-2,Length(B[i-1])],
                            [i-1,j,path_comp[i][j][k]]
                        );
                        AddDictionary(
                            f_notation_2_ext_cell,
                            [i-1,j,path_comp[i][j][k]],
                            [i-2,Length(B[i-1])]
                        );
                    od;
                fi;
            fi;
        od;
    od;

    # a final reindexing of B and removal of "*"
    # entries that once corresponded to cells of Z
    for i in [2..Length(B)] do
        for j in [1..Length(B[i])] do
            for k in [2..Length(B[i][j])] do
                B[i][j][k]:=
                B[i][j][k]-
                Length(
                    Filtered(
                        B[i-1]{[1..B[i][j][k]]},
                        x->x="*"
                    )
                );
            od;
        od;
    od;
    B:=List(B,x->Filtered(x,y->y<>"*"));
    Add(B,[]);

    return RegularCWComplex(B);
end);

# SequentialRegularCWComplexComplement
InstallGlobalFunction(
    SequentialRegularCWComplexComplement,
    function(arg...)
    local
        subdiv, method, details, map, sub, clsr, seq, tub, i;

    if Length(arg)=1 then
        arg:=[arg[1],"some","basic",false];
    fi;
    method:=arg[2];
    subdiv:=arg[3];
    details:=arg[4];

    map:=RegularCWMapToCWSubcomplex(arg[1]);
    sub:=SortedList(map[2][Length(map[2])]);
    sub:=List(sub,x->x-Position(sub,x)+1);
    clsr:=ClosureCWCell(map[1],2,sub[1])[2];;
    seq:=CWSubcomplexToRegularCWMap([map[1],clsr]);;
    tub:=RegularCWComplexComplement(seq,method,subdiv,details);
    for i in [2..Length(sub)] do
        clsr:=ClosureCWCell(tub,2,sub[i])[2];;
        seq:=CWSubcomplexToRegularCWMap([tub,clsr]);;
        tub:=RegularCWComplexComplement(seq,method,subdiv,details);
    od;
    
    return tub;
end);

# LiftColouredSurface
################################################################################
############ Input: an inclusion i:Y->X of regular CW-complexes ################
################### with component object i!.colour(n,k) returning an integer ##
################### in the range [a,b] associated to the kth n-cell of Y #######
################################################################################
########### Output: an inclusion of regular CW-complexes i~:Y~->Xx[a,b] ########
################### where Y~ is the lifted subcomplex of Xx[a,b] as ############
################### specified by i!.colour #####################################
################################################################################
InstallGlobalFunction(
    LiftColouredSurface,
    function(f)
    local
        src, trg, lens, colours,
        cbnd4_1cells, cbnd4_1cells_bnd,
        i, j, clr, closure, k,
        min, max, l_colours,
        cell, col1, col2, int,
        bnd, bbnd, l, cobnd4,
        l_src, prods;

    trg:=RegularCWMapToCWSubcomplex(f);
    src:=trg[2]*1;
    trg:=trg[1];
    lens:=List([1..Length(trg!.boundaries)-1],x->Length(trg!.boundaries[x]));
    colours:=List([1..Length(src)],x->List([1..trg!.nrCells(x-1)],y->[]));
    cbnd4_1cells:=[]; # list of all 1-cells of src whose coboundary consists of 4 2-cells
    cbnd4_1cells_bnd:=[]; # the boundary of the above 1-cells
    for i in [1..Source(f)!.nrCells(1)] do
        if Source(f)!.coboundaries[2][i][1]=4 then
            Add(cbnd4_1cells,src[2][i]);
            for j in [2,3] do
                Add(
                    cbnd4_1cells_bnd,
                    src[1][Source(f)!.boundaries[2][i][j]]
                );
            od;
        fi;
    od;
    cbnd4_1cells_bnd:=Set(cbnd4_1cells_bnd);
    # we want to compute the inclusion src* -> trg x I where
    # I is the interval [min,max] and where min and max are the
    # smallest/largest integers that f!.colour assigns to 2-cells
    #       we begin by associating to each cell e^n of src a list
    # of integers e.g. [-1,2] indicating that e^n x {-1} and
    # e^n x {2} will be in src x I
    for i in [1..Length(src[Length(src)])] do
        clr:=f!.colour(Length(src)-1,src[Length(src)][i]);
        closure:=ClosureCWCell(trg,Length(src)-1,src[Length(src)][i])[2];
        for j in [1..Length(closure)] do
            for k in [1..Length(closure[j])] do
                colours[j][closure[j][k]]:=Set(
                    Concatenation(
                        colours[j][closure[j][k]],
                        clr
                    )
                );
            od;
        od;
    od;
    min:=Minimum(List(Filtered(colours[Length(colours)],x->x<>[]),Minimum));
    # ^smallest & \/largest 'colours'
    max:=Maximum(List(Filtered(colours[Length(colours)],x->x<>[]),Maximum));
    # make a copy of trg for each colour in [min-1,min+1]
    trg:=trg!.boundaries*1; Add(trg,[]);
    l_colours:=List(
        [1..Length(trg)],
        x->List([1..Length(trg[x])],
            y->[[x-1,y],min-1]
        )
    );
    for i in [min..max+1] do
        for j in [1..lens[1]] do
            Add(trg[1],[1,0]);
            Add(l_colours[1],[[0,j],i]);
        od;
    od;
    for i in [1..Length([min..max+1])] do
        for j in [1..Length(lens)-1] do
            for k in [1..lens[j+1]] do
                cell:=trg[j+1][k]*1;
                cell:=Concatenation(
                    [cell[1]],
                    cell{[2..Length(cell)]}+lens[j]*i
                );
                Add(trg[j+1],cell);
                Add(l_colours[j+1],[[j,k],[min..max+1][i]]);
            od;
        od;
    od;
    # `join' each Target(f) x {i-1} to Target(f) x {i}
    for i in [1..Length(lens)] do
        for j in [1..lens[i]] do
            cell:=[i-1,j];
            for k in [1..Length([min-1..max])] do
                col1:=[min-1..max][k];
                col2:=[min-1..max+1][k+1];
                int:=[col1,col2];
                bnd:=[]; # boundary of cell x int
                Add(
                    bnd,
                    Position(l_colours[i],[cell,col1])
                );
                Add(
                    bnd,
                    Position(l_colours[i],[cell,col2])
                );
                if i>1 then
                    bbnd:=trg[i][j]*1;
                    bbnd:=bbnd{[2..bbnd[1]+1]};
                    bbnd:=List(
                        bbnd,
                        x->[i-2,x]
                    );
                    for l in [1..Length(bbnd)] do
                        Add(
                            bnd,
                            Position(l_colours[i],[bbnd[l],int])
                        );
                    od;
                fi;
                bnd:=Set(bnd); Add(bnd,Length(bnd),1);
                Add(trg[i+1],bnd);
                Add(l_colours[i+1],[cell,int]);
            od;
        od;
    od;
    # identify the correct subcomplex of X x [min-1,max+1]
    cobnd4:=[]; # there are some 1-cells which shouldn't be lifted
    l_src:=List(src,x->[]);
    for i in [1..Length(colours)] do
        for j in [1..Length(colours[i])] do
            if colours[i][j]<>[] then
                prods:=[];
                if Length(colours[i][j])=1 then
                    Add(prods,colours[i][j][1]);
                else
                    int:=[Minimum(colours[i][j]),Maximum(colours[i][j])];
                    for k in [2..Length(int)] do
                        Add(prods,int[k-1]);
                        if not (i=2 and j in cbnd4_1cells) and
                            not (i=1 and j in cbnd4_1cells_bnd) then
                            Add(prods,[int[k-1],int[k]]);
                        fi;
                        Add(prods,int[k]);
                    od;
                    prods:=Set(prods);
                fi;
                for k in [1..Length(prods)] do
                    if IsInt(prods[k]) then
                        Add(
                            l_src[i],
                            Position(
                                l_colours[i],
                                [[i-1,j],prods[k]]
                            )
                        );
                    else
                        Add(
                            l_src[i+1],
                            Position(
                                l_colours[i+1],
                                [[i-1,j],prods[k]]
                            )
                        );
                    fi;
                od;
            fi;
        od;
    od;
    return CWSubcomplexToRegularCWMap(
        [
            RegularCWComplex(trg),
            l_src
        ]
    );
end);

# ViewArc2Presentation
################################################################################
############ Input: three lists a, b, c. a corresponds to an arc ###############
################### presentation. b's entries are either 0, 1 or -1 and they ###
################### determine whether a given crossing in the arc presentation #
################### is an intersection, an overcrossing or an undercrossing. c #
################### is a list whose entries are 1, 2, 3 or 4. they #############
################### correspond to an intersection going from blue to blue, #####
################### blue to red, red to blue or red to red (see manual) ########
################################################################################
########### Output: a png depicting the associated coloured arc diagram ########
################################################################################
InstallGlobalFunction(
    ViewArc2Presentation,
    function(l)
    local
        arc, cross, cols, AppendTo, PrintTo, file, filetxt, filepng,
        bin_arr, coord, res, colours, i, j, x, k, y, z, clr, tmpdir;

    arc:=List(l[1],x->[SignInt(x[1])*x[1],SignInt(x[2])*x[2]]);
    cross:=l[2]*1;
    cols:=l[3]*1;
    AppendTo:=HAP_AppendTo;
    PrintTo:=HAP_PrintTo;
    #file:="/tmp/HAPtmpImage";
    #filetxt:="/tmp/HAPtmpImage.txt";
tmpdir := DirectoryTemporary();;
file:=Filename( tmpdir , "HAPtmpImage" );
filetxt:=Filename( tmpdir , "HAPtmpImage.txt" );
filepng:=Filename( tmpdir , "HAPtmpImage.png" );

    # create a binary array from the arc presentation
    bin_arr:=Sum(PureCubicalKnot(arc)!.binaryArray);
    bin_arr:=TransposedMat(bin_arr); # graham has this backwards!
    bin_arr:=List(
        bin_arr{[2..Length(bin_arr)-3]},
        x->x{[2..Length(bin_arr[1])-3]}
    );
    coord:=Concatenation(
        List(
            [1..Length(bin_arr)],
            x->List(
                [1..Length(bin_arr)],
                y->[x,y]
            )
        )
    );
    coord:=Filtered(coord,x->bin_arr[x[1]][x[2]]=2);

    res:=String(Minimum(50*Length(bin_arr),950));

    PrintTo(
        filetxt,
        "# ImageMagick pixel enumeration: ",
        Length(bin_arr)+2,
        ",",
        Length(bin_arr[1])+2,
        ",255,RGB\n"
    );

    colours:=NewDictionary([0,""],true);
    colours!.entries:=[
        [0,"(255,255,255)"], # white
        [1,"(131,205,131)"], # green
        [2,"(131,205,131)"], # green
        [3,"(131,205,131)"], # green
        [4,"(205,131,131)"], # red
        [5,"(131,131,205)"], # blue
    ];
    
    for i in [1..Length(bin_arr)] do
        for j in [1..Length(bin_arr)] do
            if bin_arr[i][j]=3 then
            # remove vertical bars that have -entry in l[1]
                x:=-1*(Int(j/3)+1);
                if x in Concatenation(l[1]) then
                    for k in [i+1..Length(bin_arr)] do
                        if bin_arr[k][j]=1 then
                            bin_arr[k][j]:=0;
                        fi;
                    od;
                fi;
            elif bin_arr[i][j]=2 then
            # check crossing type and adjust surrounding pixels
                x:=Position(coord,[i,j]);
                y:=cross[x];
                if y=-1 then
                    bin_arr[i][j-1]:=0;
                    bin_arr[i][j+1]:=0;
                elif y=1 then
                    bin_arr[i-1][j]:=0;
                    bin_arr[i+1][j]:=0;
                elif y=0 then
                    z:=cols[Position(Positions(cross,0),x)];
                    if z=1 then
                        bin_arr[i][j-1]:=5;
                        bin_arr[i][j+1]:=5;
                    elif z=2 then
                        bin_arr[i][j-1]:=5;
                        bin_arr[i][j+1]:=4;
                    elif z=3 then
                        bin_arr[i][j-1]:=4;
                        bin_arr[i][j+1]:=5;
                    elif z=4 then
                        bin_arr[i][j-1]:=4;
                        bin_arr[i][j+1]:=4;
                    fi;
                fi;
            fi;
        od;
    od;
    # image gets cropped if i dont do this
    for i in [1..2] do
        Add(bin_arr,bin_arr[1]*0);
    od;
    for j in [1..Length(bin_arr)] do
        bin_arr[j]:=Concatenation(bin_arr[j],[0,0]);
    od;
    #return bin_arr;
    # swap each entry of the binary area with an rgb value
    for i in [1..Length(bin_arr)] do
        for j in [1..Length(bin_arr[i])] do
            clr:=LookupDictionary(colours,bin_arr[i][j]);
            AppendTo(filetxt,j,",",i,": ",clr,"\n");
        od;
    od;

    # convert the binary array to an upscaled 500x500 png
    Exec(
        Concatenation("convert ",filetxt," ","-scale ",res,"x",res," ",file,".png")
    );
    # delete the old txt file
    Exec(
        Concatenation("rm ",filetxt)
    );
    # display the image
    #Exec(
    #    Concatenation(DISPLAY_PATH," ","/tmp/HAPtmpImage.png")
    #);
    Exec(
        Concatenation(DISPLAY_PATH," ",filepng)
    );
    # delete the image on window close
    #Exec(
    #    Concatenation("rm  ","/tmp/HAPtmpImage.png")
    #);
    Exec(
        Concatenation("rm  ",filepng)
    );
end);

# Tube
################################################################################
############ Input: an arc 2-presentation ######################################
################################################################################
########### Output: a regular CW-complex corresponding to the complement #######
################### of the knotted surface associated to that arc ##############
################### 2-presentation via the Tube map ############################
################################################################################
InstallGlobalFunction(
    Tube,
    function(arc)
    local
        ribbon;

    arc:=KinkArc2Presentation(arc);

    ribbon:=ArcDiagramToTubularSurface(arc);
    ribbon:=LiftColouredSurface(ribbon);
    ribbon:=SequentialRegularCWComplexComplement(ribbon,"some","none",false);

    return ribbon;
end);

# KinkArc2Presentation
################################################################################
############ Input: an arc 2-presentation ######################################
################################################################################
########### Output: an arc 2-presentation with no two crossings occuring #######
################### on the same row/column #####################################
################################################################################
InstallGlobalFunction(
    KinkArc2Presentation,
    function(arc)
    local
        mat, i, x, j, Twist, z, coord, keep, arc_;

    mat:=List([1..Length(arc[1])],x->[1..Length(arc[1])]*0);
    for i in [1..Length(arc[1])] do
        mat[Length(mat)-i+1][AbsoluteValue(arc[1][i][1])]:=3;
        mat[Length(mat)-i+1][AbsoluteValue(arc[1][i][2])]:=3;
    od;
    for x in [1,2] do
        if x=2 then
            mat:=MutableTransposedMat(mat);
        fi;
        for i in [1..Length(mat)] do
            for j in [2..Length(mat)-1] do
                if 3 in mat[i]{[1..j-1]} and
                3 in mat[i]{[j+1..Length(mat)]} then
                    mat[i][j]:=mat[i][j]+1;
                fi;
            od;
        od;
    od;
    mat:=MutableTransposedMat(mat);
    for i in Set(Concatenation(arc[1])) do
        if i<0 then
            for j in [2..Length(mat)-1] do
                if mat[j][-i] in [1,2] then
                    mat[j][-i]:=0;
                fi;
            od;
        fi;
    od;
    # mat models the arc 2-presentation, even should it correspond to
    # (a) knotted sphere(s)

    Twist:=function(i,j)
        local
            x, k, l, left, right, up, down;

        # first create a new row and column in mat 
        for x in [1..Length(mat)] do
            Add(mat[x],4,j);
        od;
        Add(
            mat,
            Concatenation(
                List([1..j-1],y->4),
                [3],
                mat[i]{[j+1..Length(mat)+1]}
            ),
            i+1
        );
        mat[i]:=Concatenation(
            mat[i]{[1..j-1]},
            [3],
            List([j+1..Length(mat)],y->4)
        );

        # now determine the values of the 2n+1 new entries
        # where n is the length of mat before Twist
        for k in [1..Length(mat)] do
            for l in [1..Length(mat[k])] do
                if mat[k][l]=4 then
                    if k in [1,Length(mat)] then
                    # the entry is either in the top/bottom row
                    # so there is no need to check above/below it
                        left:=mat[k]{[1..l-1]};
                        right:=mat[k]{[l+1..Length(mat)]};
                        if 3 in left and 3 in right then
                            if 0 in mat[k]{[l-1,l+1]} then
                            # this is to check if this row corresponds
                            # to an omitted row in the arc 2-pres.
                                mat[k][l]:=0;
                            else
                                mat[k][l]:=1;
                            fi;
                        else
                            mat[k][l]:=0;
                        fi;
                    elif l in [1,Length(mat)] then
                    # the entry is in the first/last column and
                    # we don't need to check to the left/right of it
                        up:=mat{[1..k-1]}[l];
                        down:=mat{[k+1..Length(mat)]}[l];
                        if 3 in up and 3 in down then
                            if 0 in mat{[k,k+1]}[l] then
                                mat[k][l]:=0;
                            else
                                mat[k][l]:=1;
                            fi;
                        else
                            mat[k][l]:=0;
                        fi;
                    else
                    # we are required to check above/below/left/right
                    # of the entry to determine its value
                        left:=mat[k]{[1..l-1]};
                        right:=mat[k]{[l+1..Length(mat)]};
                        up:=mat{[1..k-1]}[l];
                        down:=mat{[k+1..Length(mat)]}[l];
                        if 3 in left and 3 in right then
                            if 0 in mat[k]{[l-1,l+1]} then
                                mat[k][l]:=0;
                            else
                                mat[k][l]:=1;
                            fi;
                        elif 3 in up and 3 in down then
                            if 0 in mat{[k,k+1]}[l] then
                                mat[k][l]:=0;
                            else
                                mat[k][l]:=1;
                            fi;
                        else
                            mat[k][l]:=0;
                        fi;
                    fi;
                fi;
            od;
        od;

        return mat;
    end;

    for z in [1,2] do
        if z=2 then
            mat:=MutableTransposedMat(mat);
        fi;
        for i in [1..Sum(List(mat,y->Maximum(0,Length(Positions(y,2))-1)))] do
        # number of times we need to perform Twist on the rows of mat
            coord:=Filtered(
                Concatenation(
                    List(
                        [1..Length(mat)],
                        x->List(
                            [1..Length(mat)],
                            y->[x,y,mat[x][y]]
                        )
                    )
                ),
                x->x[3]=2
            );
            keep:=[];
            for j in [1..Length(coord)] do
                if coord[j][1] in List(coord,x->x[1]){[1..j-1]} then
                    Add(keep,j);
                fi;
            od;
            coord:=First(coord{keep});
            # the coordinate of the first instance of a crossing that
            # lies in a row where multiple crossings occur
            mat:=Twist(coord[1],coord[2]);
        od;
    od;
    mat:=MutableTransposedMat(mat);

    # before output we need to see if there are any columns that were
    # omitted. if so we need to find their new positions and reflect
    # that with a negative entry in arc[1]
    

    # now to recover the new arc 2-presentation from mat
    arc_:=ShallowCopy(arc);
    arc_[1]:=[];
    for i in Reversed([1..Length(mat)]) do
        Add(
            arc_[1],
            Filtered([1..Length(mat)],y->mat[i][y]=3)
        );
    od;

    return arc_;
end);

# NumberOfCrossingsInArc2Presentation
################################################################################
############ Input: an arc 2-presentation ######################################
################################################################################
########### Output: the number of crossings in the arc 2-presentation ##########
################################################################################
InstallGlobalFunction(
    NumberOfCrossingsInArc2Presentation,
    function(arc)
    local
        mat, i, x, j;

    mat:=List([1..Length(arc[1])],x->[1..Length(arc[1])]*0);
    for i in [1..Length(arc[1])] do
        mat[Length(mat)-i+1][AbsoluteValue(arc[1][i][1])]:=3;
        mat[Length(mat)-i+1][AbsoluteValue(arc[1][i][2])]:=3;
    od;
    for x in [1,2] do
        if x=2 then
            mat:=MutableTransposedMat(mat);
        fi;
        for i in [1..Length(mat)] do
            for j in [2..Length(mat)-1] do
                if 3 in mat[i]{[1..j-1]} and
                3 in mat[i]{[j+1..Length(mat)]} then
                    mat[i][j]:=mat[i][j]+1;
                fi;
            od;
        od;
    od;

    return Sum(List(mat,y->Length(Filtered(y,z->z=2))));
end);

# RandomArc2Presentation
################################################################################
############ Input: nothing ####################################################
################################################################################
########### Output: a random arc 2-presentation arising as the knot sum ########
################### of one or more prime knots on <12 crossings ################
################################################################################
InstallGlobalFunction(
    RandomArc2Presentation,
    function()
    local
        len, rand, prime, cross, final_arc, i, arc;
        
    len:=1;
    rand:=Random([0,1]);

    while rand=1 do
        len:=len+1;
        rand:=Random([0,1]);
    od;

    prime:=[0,0,1,1,2,3,7,21,49,165,552];

    final_arc:=Random([3..11]);
    final_arc:=[final_arc,Random([1..prime[final_arc]])];
    final_arc:=ArcPresentation(
        PureCubicalKnot(final_arc[1],final_arc[2])
    );
    
    for i in [2..len] do
        arc:=Random([3..11]);
        arc:=[arc,Random([1..prime[arc]])];
        final_arc:=ArcPresentation(
            KnotSum(
                PureCubicalKnot(final_arc),
                PureCubicalKnot(arc[1],arc[2])
            )
        );
    od;

    cross:=List(
        [1..NumberOfCrossingsInArc2Presentation([final_arc,[],[]])],
        x->Random([-1,0,1])
    );

    final_arc:=[
        final_arc,
        cross,
        List([1..Length(Filtered(cross,x->x=0))],y->Random([1..4]))    
    ];

    return final_arc;
end);
