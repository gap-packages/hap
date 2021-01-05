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

    return [map,grid,kbnd,bnd];
    return Objectify(
        HapRegularCWMap,
        rec(
            source:=RegularCWComplex(kbnd),
            target:=RegularCWComplex(bnd),
            mapping:=embed
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
end);

# Closure
################################################################################
############ Input: a CW-complex X and two integers k>=0 & i>1 #################
################################################################################
########### Output: a CW-subcomplex [X,S] corresponding to the #################
################### topological closure of the ith k-cell of X #################
################################################################################
InstallGlobalFunction(
    Closure,
    function(X,k,i)
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
end);

# ComplementTubularNeighbourhood
################################################################################
############ Input: an inclusion f: Y -> X of regular CW-complexes #############
################################################################################
########### Output: let N(Y) denote an open tubular neighbourhood of ###########
################### Y in X. this will output an inclusion B -> C where C #######
################### is homeomorphic to X \ N(Y) and B is homeomorphic to the ###
################### boundary of N(Y) ###########################################
################################################################################
InstallGlobalFunction(
    ComplementTubularNeighbourhood,
    function(f)
end);
