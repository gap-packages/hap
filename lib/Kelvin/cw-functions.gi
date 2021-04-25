
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

# AddCell
################################################################################
############ Input: a cell complex B, an integer k>=0 and lists of #############
################### positive integers b and c ##################################
################################################################################
########### Output: the cell complex B with an added k-cell whose ##############
################### boundary is specified by b and whose coboundary ############
################### is specified by c ##########################################
################################################################################
InstallGlobalFunction(
    AddCell,
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
        i, j, bnd;

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
    for i in [1..Length(Y[Length(closure)])] do
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
        f, compatible, Y, B, IsInternal, count, total, path_comp,
        i, j, clsr, int, crit, bary, IsSubpathComponent,
        ext_cell_2_f_notation, f_notation_2_ext_cell,
        ext_cells, k, ext_cell_bnd, e_n_bar,
        ext_cell_cbnd;

    if Length(arg)=1 then
        f:=ShallowCopy(arg[1]);
        compatible:=false;
    elif Length(arg)=2 then
        f:=ShallowCopy(arg[1]);
        compatible:=arg[2];
    fi;

    Y:=RegularCWMapToCWSubcomplex(f);
    B:=List([1..Dimension(Y[1])+1],x->[]);
    
    # sub:=List([1..Dimension(Y[1])+1],x->[]);
    # keeps track of the indexing of those cells
    # which comprise the boundary of Y \ N(Z)

    IsInternal:=function(n,k) # is the kth n-cell of Y in Y\Z?
        if n+1>Length(Y[2]) then
            return true;
        elif k in Y[2][n+1] then
            return false;
        fi;
        return true;
    end;

    if not compatible then
        Print("Testing contractibility...\n");
    fi;
    count:=0;
    total:=Sum(List(Y[1]!.boundaries,Length));
    bary:=[];

    path_comp:=List([1..Length(Y[1]!.boundaries)-1],x->[]);
    # list of all of the path components of the intersection between
    # the closure of each internal cell and the subcomplex Z < Y
    # note: entries may be an empty list if they correspond to an
    #       empty intersection or they may not be assigned a
    #       value at all if the associated cell lies in Z
    for i in [1..Length(Y[1]!.boundaries)-1] do
        for j in [1..Length(Y[1]!.boundaries[i])] do
            count:=count+1;
            if IsInternal(i-1,j) then
                Add(B[i],Y[1]!.boundaries[i][j]);
                # output will contain all internal cells
                clsr:=ClosureCWCell(Y[1],i-1,j);
                int:=IntersectionCWSubcomplex(clsr,Y);
                
                path_comp[i][j]:=PathComponentsCWSubcomplex(int);

                # CONTRACIBILITY TEST
                # must be a non-empty subcomplex of dimension > 0
                if not compatible then
                    Print(count," out of ",total," cells tested.","\r");

                    # test the contractibility of each path component!
                    # if you find anything non contractible,
                    # subdivide the problem cells and restart
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
    
    if not compatible then
        if bary=[] then
            Print("\nThe input is compatible with this algorithm.\n");
        else
            Print("\nBarycentrically subdividing ",Length(bary)," cell(s):\n");
            for i in [1..Length(bary)] do
                f:=BarycentricallySubdivideCell(
                    f,
                    bary[i][1],
                    bary[i][2]
                );
                Print(Int(100*i/Length(bary)),"\% complete. \r");
            od;
            Print("\n");
            return RegularCWComplexComplement(f,true); # bypass compatibility check
        fi;
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

                        AddCell(
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

    # now to find the boundary of Y \ N(Z)
    # and form the inclusion map
    #for i in [3..Length(B)-1] do
    #    for j in [1..Length(B[i])] do
    #        for k in [2..Length(B[i][j])] do
    #            Add(sub[i-1],B[i][j][k]);
    #        od;
    #    od;
    #od;
    #sub:=List(
    #    List(
    #        sub,
    #        x->Filtered(
    #            x,
    #            y->Length(Positions(x,y))=1
    #        )
    #    ),
    #    Set
    #);
    #for i in Reversed([2..Length(sub)]) do
    #    for j in [1..Length(sub[i])] do
    #        for k in B[i][j]{[2..B[i][j][1]+1]} do
    #            Add(sub[i-1],k);
    #        od;
    #    od;
    #od;
    #sub:=List(sub,Set);

    return RegularCWComplex(B);#BoundaryMap(RegularCWComplex(B));
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

# ViewColouredArcDiagram
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
    ViewColouredArcDiagram,
    function(arc,cross,levels)
    local
        AppendTo, PrintTo, file, filetxt, ArcPresentationToArray,
        bin_arr, i, res, crs, crs_coord, j, colours, n, bar1, loc_3,
        bar2, nsew, clr;

    AppendTo:=HAP_AppendTo;
    PrintTo:=HAP_PrintTo;
    file:="/tmp/HAPtmpImage";
    filetxt:="/tmp/HAPtmpImage.txt";

    ArcPresentationToArray:=function(arc)
        local
            arr, i, pair, j, k;

        arr:=List([0..5*(Length(arc)-1)],x->List([0..5*(Length(arc)-1)])*0);

        for i in [1..Length(arc)] do
            pair:=arc[i]*1;
            arr[Length(arr)-5*(i-1)][1+5*(pair[1]-1)]:=3;
            arr[Length(arr)-5*(i-1)][1+5*(pair[2]-1)]:=3;
        od;

        for i in [1..Length(arr)] do
            for j in [1..Length(arr[i])] do
                if arr[i][j]=3 then
                    k:=j*1;
                    while arr[i][k+1]<>3 do
                        k:=k+1;
                        arr[i][k]:=1;
                    od;
                    break;
                fi;
            od;
        od;
        for i in [1..Length(arr)-1] do
            for j in [1..Length(arr[i])] do
                if arr[i][j]=3 then
                    if 3 in List(arr{[i+1..Length(arr)]},x->x[j]) then
                        k:=i*1;
                        while arr[k+1][j]<>3 do
                            k:=k+1;
                            arr[k][j]:=arr[k][j]+1;
                        od;
                    fi;
                fi;
            od;
        od;

        return arr;
    end;

    bin_arr:=FrameArray(ArcPresentationToArray(arc))-3;
    bin_arr:=List(bin_arr,x->Concatenation(x,[-3,-3]));
    for i in [1,2] do
        Add(bin_arr,bin_arr[1]*0-3);
    od;

    # resolution of the output, chosen to display on (my) laptop nicely
    res:=String(Minimum(50*Length(bin_arr),950));

    crs:=List(bin_arr,x->Positions(x,-1));
    crs_coord:=[]; # the co-ordinates of each crossing
    for i in [1..Length(crs)] do
        if crs[i]<>[] then
            for j in [1..Length(crs[i])] do
                Add(crs_coord,[i,crs[i][j]]);
            od;
        fi;
    od;

    #if Length(crs_coord)<>Length(cross) then
    #    Error("number of specified crossings is incorrect");
    #elif levels<>[] then
    #    if Length(Positions(cross,0))<Length(levels) then
    #        Error("all ",Length(Positions(cross,0))," 4-d crossing(s) must be accounted for");
    #    fi;
    #fi;

    PrintTo(
        filetxt,
        "# ImageMagick pixel enumeration: ",
        Length(bin_arr[1]),
        ",",
        Length(bin_arr),
        ",255,RGB\n"
    );

    colours:=NewDictionary([0,""],true);
    colours!.entries:=[
        [-3,"(255,255,255)"], # white
        [-2,"(0,255,255)"], # cyan
        [-1,"(0,0,255)"], # blue
        [0,"(0,255,0)"], # green
        [1,"(255,0,0)"], # red
        [2,"(255,255,0)"], # yellow
        [3,"(255,0,255)"] # magenta
    ];

    for i in [1..Length(bin_arr)] do
        for j in [1..Length(bin_arr[i])] do
            if bin_arr[i][j] in [-1,0] then
                bin_arr[i][j]:=3; # temp. colour for the corners/crossings
            elif bin_arr[i][j]=-2 then
                bin_arr[i][j]:=0; # everything should be green 
            fi; # unless there's a 4-d crossing involved
        od;
    od;

    n:=0;
    for i in [1..Length(crs_coord)] do
        if cross[i]=0 then
            n:=n+1;
            bar1:=List([1..crs_coord[i][1]-1]);
            loc_3:=Positions(List(bar1,x->bin_arr[x][crs_coord[i][2]]),3);
            bar1:=bar1{
                [loc_3[Length(loc_3)]+1..Length(bar1)]
            };
            bar2:=List([crs_coord[i][1]+1..Length(bin_arr)]);
            loc_3:=Positions(List(bar2,x->bin_arr[x][crs_coord[i][2]]),3);
            bar2:=bar2{
                [1..loc_3[1]-1]
            };
            for j in bar1 do
                if levels[n] in [1,2] then
                    bin_arr[j][crs_coord[i][2]]:=-1;
                else
                    bin_arr[j][crs_coord[i][2]]:=1;
                fi;
            od;
            for j in bar2 do
                if levels[n] in [1,3] then
                    bin_arr[j][crs_coord[i][2]]:=-1;
                else
                    bin_arr[j][crs_coord[i][2]]:=1;
                fi;
            od;
        elif cross[i]=-1 then
            j:=crs_coord[i]*1;
            bin_arr[j[1]][j[2]-1]:=-3;
            bin_arr[j[1]][j[2]+1]:=-3;
        fi;
    od;
    for i in [1..Length(crs_coord)] do
        if cross[i]=1 then
            j:=crs_coord[i]*1;
            bin_arr[j[1]-1][j[2]]:=-3;
            bin_arr[j[1]+1][j[2]]:=-3;
        fi;
    od;

    for i in [1..Length(bin_arr)] do
        for j in [1..Length(bin_arr[i])] do
            if bin_arr[i][j]=3 then
                nsew:=[
                    bin_arr[i-1][j],
                    bin_arr[i+1][j],
                    bin_arr[i][j+1],
                    bin_arr[i][j-1]
                ];
                if -3 in nsew then
                    if Length(Positions([nsew[1],nsew[2]],-3))=1 then # corner
                        if -1 in nsew then
                            bin_arr[i][j]:=-2;
                        elif 1 in nsew then
                            bin_arr[i][j]:=2;
                        else
                            bin_arr[i][j]:=0;
                        fi;
                    else # int. point
                        if nsew[1]=-3 then
                            bin_arr[i][j]:=0;
                        else
                            if nsew[1]=nsew[2] then
                                bin_arr[i][j]:=nsew[1]*1;
                            elif Set([nsew[1],nsew[2]])=Set([1,0]) then
                                bin_arr[i][j]:=2;
                            elif Set([nsew[1],nsew[2]])=Set([-1,0]) then
                                bin_arr[i][j]:=-2;
                            fi;
                        fi;
                    fi;
                else
                    bin_arr[i][j]:=0;
                fi;
            fi;
        od;
    od;

    # colour entries according to levels
    for i in [1..Length(bin_arr)] do
        for j in [1..Length(bin_arr[i])] do
            clr:=LookupDictionary(colours,bin_arr[i][j]);
            AppendTo(filetxt,j,",",i,": ",clr,"\n");
        od;
    od;

    # convert the binary array to an upscaled png
    Exec(
        Concatenation("convert ",filetxt," ","-scale ",res,"x",res," ",file,".png")
    );
    # delete the old txt file
    Exec(
        Concatenation("rm ",filetxt)
    );
    # display the image
    Exec(
        Concatenation(DISPLAY_PATH," ","/tmp/HAPtmpImage.png")
    );
    # delete the image on window close
    Exec(
        Concatenation("rm  ","/tmp/HAPtmpImage.png")
    );
end);
