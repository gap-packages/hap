# ArcPresentationToKnottedOneComplex_Alt
################################################################################
############ Input: an arc presentation ########################################
################################################################################
########### Output: an inclusion map of regular CW-complexes from ##############
################### a 1-dimensional link to a 3-dimensional ball ###############
################################################################################
InstallGlobalFunction(
    ArcPresentationToKnottedOneComplex_Alt,
    function(arc)
# this method differs than the one below by outputting an embedding which is
# always compatible with the RegularCWComplexComplement function (i.e. there are
# no cells which intersect with the subcomplex in a non-contractible way)
    local
        gn, grid;

    gn:=Length(arc); # grid number

    grid:=List([1..2*gn]);
end);

# ArcPresentationToKnottedOneComplex
################################################################################
############ Input: an arc presentation ########################################
################################################################################
########### Output: an inclusion map of regular CW-complexes from ##############
################### a 1-dimensional link to a 3-dimensional ball ###############
################################################################################
InstallGlobalFunction(
	ArcPresentationToKnottedOneComplex,
	function(arc)
    local
        gn, grid, i, kbnd, bnd, map, imap,
        IsIntersection, ints, ext_ints, j, 0c,
        B2Decomposition, B3Decomposition, embed;

    gn:=Length(arc); # the grid number
    
    grid:=List([1..5*gn],x->[1..5*gn]*0); # form a (3*gn) x gn
    for i in [0..gn-1] do # matrix from the arc presentation
        grid[5*(gn-i)-2][5*arc[i+1][1]-2]:=1;
        grid[5*(gn-i)-2][5*arc[i+1][2]-2]:=1;
    od;
    grid:=FrameArray(grid);

    kbnd:=List([1..3],x->[]); # boundary list of the knot
    bnd:=List([1..5],x->[]); # boundary list of the complement of the knot
    map:=List([1..2],x->[]); # inclusion map from kbnd to bnd
    imap:=List([1..4],x->[]); # inverse image of the above inclusion map

    IsIntersection:=function(i,j) # finds where crossings occur in the grid
        if grid[i][j]=0 then
            if 1 in grid[i]{[1..j]} then
                if 1 in grid[i]{[j..5*gn+2]} then
                    if 1 in List([1..i],x->grid[x][j]) then
                        if 1 in List([i..5*gn+2],x->grid[x][j]) then
                            return true;
                        fi;
                    fi;
                fi;
            fi;
        fi;

        return false;
    end;

    ints:=[]; # records the coordinates of each crossing
    ext_ints:=[]; # records those 0-cells which arise in intersection points
    for i in [1..5*gn+2] do
        for j in [1..5*gn+2] do
            if IsIntersection(i,j) then
                Add(ext_ints,[i-1,j]); Add(ext_ints,[i+1,j]);
                Add(ints,[i,j]);
                grid[i-1][j]:='*';
                grid[i][j]:='*';
                grid[i+1][j]:='*';
            fi;
        od;
    od;

    0c:=4; # label the entries of grid so that it models the 0-skeleton
    for i in [1..5*gn+2] do
        for j in [1..5*gn+2] do
            if grid[i][j]<>0 then
                grid[i][j]:=0c;
                0c:=0c+1;
            fi;
        od;
    od;
    0c:=0c+2;

    ints:=List(ints,x->grid[x[1],x[2]]);

    B2Decomposition:=function()
# takes what we have so far and uses it to form a regular CW-decomposition of
# the 2-ball with the appropriate inclusion map
        local i, j, hslice, vslice, 2SkeletonOfDisk, DuplicateDisk;

        kbnd[1]:=List([1..0c-6],x->[1,0]);
        bnd[1]:=List([1..0c],x->[1,0]);
        map[1]:=[1..0c-6]+3;
        imap[1]:=Concatenation([0,0,0],map[1]-3,[0,0,0]);

        Add(bnd[2],[2,1,2]); Add(bnd[2],[2,1,3]); Add(bnd[2],[2,2,0c-2]);
        Add(bnd[2],[2,3,0c-1]); Add(bnd[2],[2,0c-2,0c]);
        Add(bnd[2],[2,0c-1,0c]);
        Add(bnd[2],[2,3,4]); Add(bnd[2],[2,0c-3,0c-2]);
        # add some 1-cells to just the complement to stay regular
        # these act as a frame to the knot

        for i in [1..8] do
            Add(imap[2],0);
        od;

        for i in [1..5*gn] do # add the horizontal arcs of the knot first
            hslice:=[];
            for j in [1..5*gn] do
                if grid[i][j]<>0 and not [i,j] in ext_ints then
                    Add(hslice,grid[i][j]);
                fi;
            od;
            for j in [1..Length(hslice)-1] do
                Add(kbnd[2],[2,hslice[j]-3,hslice[j+1]-3]);
                Add(bnd[2],[2,hslice[j],hslice[j+1]]);
                Add(map[2],Length(bnd[2]));
                Add(imap[2],Length(kbnd[2]));
            od;
        od;
        for j in [1..5*gn] do # now add the vertical arcs
            vslice:=[];
            for i in [1..5*gn] do
                if grid[i][j]<>0 then
                    Add(vslice,grid[i][j]);
                fi;
            od;
            for i in [1..Length(vslice)-1] do
                Add(bnd[2],[2,vslice[i],vslice[i+1]]);
                if not(vslice[i] in ints or vslice[i+1] in ints) then
                    Add(kbnd[2],[2,vslice[i]-3,vslice[i+1]-3]); # introduce a 
                    Add(map[2],Length(bnd[2])); # gap in the knot to be 'lifted' 
                    Add(imap[2],Length(kbnd[2])); # and filled in later
                else
                    Add(imap[2],0);
                fi;
            od;
        od;

        2SkeletonOfDisk:=function(bnd)
            local
                ori, i, j, cell, top, rgt, btm, lft,
                Clockwise, path, FaceTrace;

            grid[1][1]:=1; grid[1][5*gn+2]:=2; grid[4][1]:=3;
            grid[5*gn-1][5*gn+2]:=0c-2; grid[5*gn+2][1]:=0c-1;
            grid[5*gn+2][5*gn+2]:=0c;

            ori:=List([1..0c],x->[1..4]*0);
            # each 0-cell will have the 0-cells N/E/S/W of it
            # listed in that order 

            for i in [1..5*gn+2] do
                for j in [1..5*gn+2] do
                    cell:=grid[i][j];
                    if cell<>0 then
                        top:=List([1..i-1],x->grid[x][j]);
                        top:=Filtered(top,x->x<>0);
                        if top<>[] then
                            ori[cell][1]:=Position
                            (
                                bnd[2],
                                Concatenation([2],Set([cell,top[Length(top)]]))
                            );
                        fi;
                        rgt:=grid[i]{[j+1..5*gn+2]};
                        rgt:=Filtered(rgt,x->x<>0);
                        if rgt<>[] then
                            ori[cell][2]:=Position
                            (
                                bnd[2],
                                Concatenation([2],Set([cell,rgt[1]]))
                            );
                        fi;
                        btm:=List([i+1..5*gn+2],x->grid[x][j]);
                        btm:=Filtered(btm,x->x<>0);
                        if btm<>[] then
                            ori[cell][3]:=Position(
                                bnd[2],
                                Concatenation([2],Set([cell,btm[1]]))
                            );
                        fi;
                        lft:=grid[i]{[1..j-1]};
                        lft:=Filtered(lft,x->x<>0);
                        if lft<>[] then
                            ori[cell][4]:=Position(
                                bnd[2],
                                Concatenation([2],Set([cell,lft[Length(lft)]]))
                            );
                        fi;
                        if [i,j] in ext_ints then # 0-cells in ext_ints
                            ori[cell][2]:=0; # never have edges from the 
                            ori[cell][4]:=0; # left or right of them
                        fi;
                    fi;
                od;
            od;

            
################################################################################
# repurposed from KnotComplement and KnotComplementWithBoundary
            FaceTrace:=function(path)
                local
                    unselected, sourceORtarget, x, ClockwiseTurn,
                    2cell, sORt, dir, e1, e0, i, bool;

                unselected:=Concatenation
                (
                    [1..Length(bnd[2])],
                    [7..Length(bnd[2])] # the first 6 1-cells occur just once
                );

                ClockwiseTurn:=function(p,e)
                    local f;
                    
                    f:=(Position(p,e) mod 4)+1;
                    while p[f]=0 do
                        f:=(f mod 4)+1;
                    od;
                    
                    return p[f];
                end;

                bool:=false;
                sourceORtarget:=List([1..Length(bnd[2])],y->[3,2]);
                x:=1;
                while unselected<>[] do
                    while (not x in unselected) and (not e1 in unselected) do
                        x:=x+1;
                    od;
                    2cell:=[x];
                    sORt:=sourceORtarget[x][Length(sourceORtarget[x])];
                    Unbind(sourceORtarget[x][Length(sourceORtarget[x])]);

                    dir:=path[bnd[2][x][sORt]];
                    e0:=bnd[2][x][sORt];
                    e1:=ClockwiseTurn(dir,x);
                    while e1<>x do
                        Add(2cell,e1);
                        e0:=Filtered(bnd[2][e1]{[2,3]},y->y<>e0)[1];
                        dir:=path[e0];
                        e1:=ClockwiseTurn(dir,e1);
                    od;
                    Add(2cell,Length(2cell),1);
                    if (not Set(2cell) in List(bnd[3],x->Set(x))) then
                        for i in Filtered(
                            2cell{[2..Length(2cell)]},
                            y->y in unselected
                        ) do
                            Unbind(unselected[Position(unselected,i)]);
                        od;
                        if not bool then # to save some checks
                            if Set(2cell)=[1,2,3,4,5,6] then
                                bool:=true;
                            else
                                Add(bnd[3],2cell);
                                Add(imap[3],0);
                            fi;
                        else
                            Add(bnd[3],2cell);
                            Add(imap[3],0);
                        fi;
                    fi;
                od;
################################################################################
                bnd[3]:=List(
                    bnd[3],
                    x->Concatenation(
                        [x[1]],
                        Set(x{[2..Length(x)]})
                    )
                ); # order the 2-cells nicely
            end;
            FaceTrace(ori);
        end;
        2SkeletonOfDisk(bnd);
    end;
    B2Decomposition();

    B3Decomposition:=function()
        local 
            k0, b0, k1, b1, b2, DuplicateDisk, b22, CrossI, 3c;

        k0:=Length(kbnd[1]); b0:=Length(bnd[1]);
        k1:=Length(kbnd[2]); b1:=Length(bnd[2]);
        b2:=Length(bnd[3]);

        DuplicateDisk:=function() # make a duplicate of everything
            local i, n, mult;

            for i in [1..Length(ext_ints)/2] do
                Add(kbnd[1],[1,0]);
                Add(kbnd[1],[1,0]);
                Add(kbnd[1],[1,0]);
            od;

            bnd[1]:=Concatenation(bnd[1],bnd[1]);

            imap[1]:=Concatenation(imap[1],imap[1]+Length(kbnd[1])-k0);
            for i in [b0+1..2*b0] do
                if imap[1][i]=Length(kbnd[1])-k0 then
                    imap[1][i]:=0;
                fi;
            od;

            n:=k0+1;
            mult:=1;
            for i in [1..b1] do
                Add(bnd[2],bnd[2][i]+[0,b0,b0]);
                if i>8 and not (bnd[2][i]-[0,3,3] in kbnd[2]) then
                    Add(
                        kbnd[2],
                        Concatenation
                        (
                            [2],
                            [n,n+1]
                        )
                    );
                    n:=n+1+Int(mult/2);

                    Add(map[1],bnd[2][i][2]+b0);
                    if Int(mult/2)=1 then 
                        mult:=1;
                        Add(map[1],bnd[2][i][3]+b0);
                    else
                        mult:=2;
                    fi;

                    Add(map[2],Length(bnd[2]));
                    Add(imap[2],Length(map[2]));
                else
                    Add(imap[2],0);
                fi;
            od;

            for i in [1..b2] do
                Add(
                    bnd[3],
                    Concatenation(
                        [bnd[3][i][1]],
                        bnd[3][i]{[2..Length(bnd[3][i])]}+b1
                    )
                );
                Add(imap[3],0);
            od;

            return Length(bnd[3]);
        end;
        b22:=DuplicateDisk();

        CrossI:=function() # connect the two disks 
            local l, i, j, n, bool, 3cell;
            # each n-cell gives rise to an (n+1)-cell

            l:=[];
            for i in [1..5*gn+2] do
                for j in [1..5*gn+2] do
                    if grid[i][j]<>0 then
                        Add(l,grid[i][j]);
                    fi;
                od;
            od;

            n:=k0+1;
            bool:=false;
            for i in l do
                Add(bnd[2],[2,i,i+b0]);
                if i in List(ext_ints,x->grid[x[1]][x[2]]) then
                    Add(kbnd[2],[2,i-3,n]);
                    if not bool then
                        n:=n+2; bool:=true;
                    else
                        n:=n+1; bool:=false;
                    fi;
                    Add(map[2],Length(bnd[2]));
                    Add(imap[2],Length(map[2]));
                else
                    Add(imap[2],0);
                fi;
            od;

            for i in [1..b1] do
                Add(
                    bnd[3],
                    [
                        4,
                        i,
                        i+b1,
                        Position(
                            bnd[2],
                            [
                                2,
                                bnd[2][i][2],
                                bnd[2][i+b1][2]
                            ]
                        ),
                        Position(
                            bnd[2],
                            [
                                2,
                                bnd[2][i][3],
                                bnd[2][i+b1][3]
                            ]
                        )
                    ]
                );
                Add(imap[4],[i,Length(bnd[3])]); # not exactly the inverse
            od; # image any more, but this list is used directly below

            for i in [1..b2] do
                3cell:=Concatenation(
                    [
                        i,
                        i+b2
                    ],
                    List(
                        bnd[3][i]{[2..Length(bnd[3][i])]},
                        x->imap[4]
                            [
                                Position(
                                    List(imap[4],y->y[1]),
                                    x
                                )
                            ][2]
                    )
                );
                Add(3cell,Length(3cell),1);
                Add(bnd[4],3cell);
            od;
        end;
        CrossI();

        # add a cap to B3 /// not sure if necessary
        Add(bnd[3],[6,1,2,3,4,5,6]);
        Add(bnd[3],[6,b1+1,b1+2,b1+3,b1+4,b1+5,b1+6]);
        3c:=Concatenation([Length(bnd[3])-1],[1..b2]);
        Add(3c,Length(3c),1);
        Add(bnd[4],3c);
        3c:=Concatenation([Length(bnd[3])],[b2+1..b22]);
        Add(3c,Length(3c),1);
        Add(bnd[4],3c);

    end;
    B3Decomposition();

    embed:={n,k}->map[n+1][k];

    #return [map,grid,kbnd,bnd];
    return Objectify(
        HapRegularCWMap,
        rec(
            source:=RegularCWComplex(kbnd),
            target:=RegularCWComplex(bnd),
            mapping:=embed,
			grid:=grid
        )
    );
end);

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

# RegularCWComplexComplement
################################################################################
############ Input: an inclusion f: Z -> Y of regular CW-complexes #############
################################################################################
########### Output: let N(Z) denote an open tubular neighbourhood of ###########
################### Z in Y. this will output an inclusion B -> C where C #######
################### is homeomorphic to Y \ N(Z) and B is homeomorphic to the ###
################### boundary of N(Z) ###########################################
################################################################################
InstallGlobalFunction(
    RegularCWComplexComplement,
    function(arg...)
    local
        f, compatible, Y, B, IsInternal, count, total, path_comp,
        i, j, clsr, int, hom_lst, bary, IsSubpathComponent,
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
                    # if you find anything non contractible, barycentrically
                    # subdivide the problem cells and restart with optional
                    # 'true' cheat code 8-)
                        for k in [1..Length(path_comp[i][j])] do
                            if path_comp[i][j][k][2]<>List(path_comp[i][j][k][2],x->[]) and path_comp[i][j][k][2][2]<>[] then
                                hom_lst:=List(
                                    [1..Length(path_comp[i][j][k][2])-1],
                                    x->Homology(
                                        Source(
                                            CWSubcomplexToRegularCWMap(
                                                path_comp[i][j][k]
                                            )
                                        ),
                                        x
                                    )
                                );
                                if hom_lst<>List(hom_lst,x->[]) then
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
    
    for i in [1..Length(bary)] do
        if bary[i][2] in Concatenation(List(Filtered(bary,x->x[1]=bary[i][1]+1),x->Y[1]!.boundaries[x[1]+1][x[2]])) then
            Unbind(bary[i]);
        fi;
    od;
    bary:=Set(bary);

    if not compatible then
        if bary=[] then
            Print("\nThe input is compatible with this algorithm.\n");
        else
            Print("\n",Length(bary)," cell(s) will be barycentrically subdivided.\n");
            for i in [1..Length(bary)] do
                f:=BarycentricallySubdivideCell(
                    f,
                    bary[i][1],
                    bary[i][2]
                );
            od;
            return RegularCWComplexComplement(f,true);
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

    return RegularCWComplex(B);
end);
