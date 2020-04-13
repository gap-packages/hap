################################################################################
################### Regular CW-decomposition of knot complements ###############
################################################################################
############ Input: an arc presentation of some knot. ##########################
################################################################################
########### Output: a regular CW-complex representing the complement of ########
################### the input knot/link. #######################################
################################################################################
InstallGlobalFunction(
    KnotComplement,
    function(arg...)
    local
        rand, D, len, signless, PuncturedDisk,
        P, grid, loop_correction, PuncturedTube, i;

    if Length(arg)>1
        then
        rand:=true;
        Print("Random 2-cell selection is enabled.\n");
    else
        rand:=false;
    fi;

    D:=arg[1];
    len:=Length(D);
    signless:=List(D,x->[AbsInt(x[1]),AbsInt(x[2])]);

    PuncturedDisk:=function(D)
        local
            grid, i, IsIntersection,
            CornerConfiguration, bound,
            bigGrid, GridFill, j, tick,
            hslice, vslice, k, 0c, Orient,
            path, FaceTrace, traced_bound, cgrid;

        grid:=List([1..len],x->List([1..len],y->0));
        for i in [1..len]
            do
            grid[len-i+1][D[i][1]]:=1;
            grid[len-i+1][D[i][2]]:=1;
        od;

        IsIntersection:=function(i,j)

            if grid[i][j]=0
                then
                if 1 in grid[i]{[1..j]}
                    then
                    if 1 in grid[i]{[j..len]}
                        then
                        if 1 in List([1..i],x->grid[x][j])
                            then
                            if 1 in List([i..len],x->grid[x][j])
                                then
                                return true;
                            fi;
                        fi;
                    fi;
                fi;
            fi;

        return false;
        end;

        CornerConfiguration:=function(i,j);

            if grid[i][j]=1
                then
                if Size(Positions(grid[i]{[j..len]},1))=2
                    then
                    if Size(Positions(List([i..len],x->grid[x][j]),1))=2
                        then # Corner type 1, i.e : __
                        return 1; #                |
                    elif Size(Positions(List([1..i],x->grid[x][j]),1))=2
                        then # Corner type 3, i.e :
                        return 3; #                |__
                    fi;
                elif Size(Positions(grid[i]{[1..j]},1))=2
                    then
                    if Size(Positions(List([i..len],x->grid[x][j]),1))=2
                        then # Corner type 2, i.e : __
                        return 2; #                   |
                    elif Size(Positions(List([1..i],x->grid[x][j]),1))=2
                        then # Corner type 4, i.e :
                        return 4; #                 __|
                    fi;
                fi;
            fi;

            return 0;
        end;

        bound:=[[],[],[],[],[]];
        bigGrid:=List([1..2*len],x->List([1..2*len],y->0));

        GridFill:=function(c,i,j);
# places an * at each point where a 0-cell is to be added to bigGrid
            if c=1 or c=4
                then
                bigGrid[(2*i)-1][(2*j)-1]:='*';
                bigGrid[2*i][2*j]:='*';
            elif c=2 or c=3
                then
                bigGrid[(2*i)-1][2*j]:='*';
                bigGrid[2*i][(2*j)-1]:='*';
            fi;
        end;

        for i in [1..len]
            do # loop through bigGrid and add temporary *s
            for j in [1..len]
                do
                if IsIntersection(i,j)
                    then # four 0-cells at an intersection
                    bigGrid[(2*i)-1][(2*j)-1]:='*';
                    bigGrid[(2*i)-1][2*j]:='*';
                    bigGrid[2*i][(2*j)-1]:='*';
                    bigGrid[2*i][2*j]:='*';
                elif grid[i][j]=1
                    then # two 0-cells at the endpoints of each horizontal bar
                    GridFill(CornerConfiguration(i,j),i,j);
                fi;
            od;
        od;

        tick:=2;
        for i in [1..2*len]
            do # number the 0-cells row-by-row
            for j in [1..2*len]
                do
                if bigGrid[i][j]='*'
                    then
                    bigGrid[i][j]:=tick;
                    tick:=tick+1;
                fi;
            od;
        od;

        for i in [1..2*len]
            do # connect all 0-cells that lie in the same
            hslice:=[]; # horizontal/vertical 'slice'
            vslice:=[];
            for j in [1..2*len]
                do
                if not bigGrid[i][j]=0
                    then
                    Add(hslice,bigGrid[i][j]);
                fi;
                if not bigGrid[j][i]=0
                    then
                    Add(vslice,bigGrid[j][i]);
                fi;
            od;
            for k in [1..Length(hslice)]
                do
                if Length(hslice)>k
                    then
                    Add(bound[2],[2,hslice[k],hslice[k+1]]);
                fi;
            od;
            for k in [1..Length(vslice)]
                do
                if Length(vslice)>k
                    then
                    Add(bound[2],[2,vslice[k],vslice[k+1]]);
                fi;
            od;
        od;

        for i in [1..len]
            do # add the looping 1-cells to the 1-skeleton
            for j in [1..len]
                do
                if CornerConfiguration(i,j) in [1,4]
                    then
                    Add(bound[2],[
                        2,
                        bigGrid[(2*i)-1][(2*j)-1],
                        bigGrid[2*i][2*j]
                        ]
                    );
                    Add(bound[2],[
                        2,
                        bigGrid[(2*i)-1][(2*j)-1],
                        bigGrid[2*i][2*j]
                        ]
                    );
                elif CornerConfiguration(i,j) in [2,3]
                    then
                    Add(bound[2],[
                        2,
                        bigGrid[(2*i)-1][2*j],
                        bigGrid[2*i][(2*j)-1]
                        ]
                    );
                    Add(bound[2],[
                        2,
                        bigGrid[(2*i)-1][2*j],
                        bigGrid[2*i][(2*j)-1]
                        ]
                    );
                fi;
            od;
        od;

        0c:=Maximum(List(bigGrid,x->Maximum(x)));
        for i in [1..0c+1]
            do
            Add(bound[1],[1,0]);
        od;

        Add(bound[2],[2,1,2]); # connect the central component to the
        Add(bound[2],[2,Length(bound[1])-1,Length(bound[1])]); # circumference
        Add(bound[2],[2,1,Length(bound[1])]); # of the disk
        Add(bound[2],[2,1,Length(bound[1])]);

        bigGrid:=FrameArray(bigGrid);
        bigGrid[1][2]:=1; # Adds the first and last 0-cells to bigGrid
        bigGrid[Length(bigGrid)][Length(bigGrid[1])-1]:=0c+1;

        Orient:=function(bound)
# traces the 1-skeleton in a clockwise walk to yield the 2-cells
            local 
                unchosen, neighbours, i, j,
                Clockwise;

            unchosen:=List(ShallowCopy(bound[2]),x->[x[2],x[3]]);
            neighbours:=List(ShallowCopy(bound[1]),x->[]);

            for i in [1..Length(bound[1])]
                do
                for j in [1..Length(unchosen)]
                    do
                    if i in unchosen[j]
                        then
                        Add(neighbours[i],j);
                    fi;
                od;
            od;

            Clockwise:=function(neighbours)
                local # orders clockwise the neighbours of each 0-cell
                    oriented, first0, last0,
                    i, j, x, k, l, posi, posx;

                oriented:=List(neighbours,x->List([1..12],y->"pass"));
                first0:=SortedList(neighbours[1]);
                last0:=SortedList(neighbours[Length(neighbours)]);

                oriented[1][7]:=first0[1];
                oriented[1][6]:=first0[3];
                oriented[1][8]:=first0[2];
            # these two orderings are always fixed;
            # they correspond to the circumferential edges
                oriented[Length(oriented)][1]:=last0[1]; 
                oriented[Length(oriented)][2]:=last0[3];
                oriented[Length(oriented)][12]:=last0[2];

                for i in [2..Length(neighbours)-1]
                    do # excludes the 1st and last 0-cells
                    for j in [1..Length(neighbours[i])]
                        do # x is a neighbouring 0-cell to i
                        x:=bound[2][neighbours[i][j]];
                        x:=Filtered(x{[2,3]},y->y<>i)[1];
                        for k in [1..Length(bigGrid)]
                            do
                            for l in [1..Length(bigGrid[1])]
                                do
                                if i=bigGrid[k][l]
                                    then
                                    posi:=[k,l];
                                fi;
                                if x=bigGrid[k][l]
                                    then
                                    posx:=[k,l];
                                fi;
                            od;
                        od;
                        # below are the checks for orientation,
                        # there are 12 in total (two for each diagonal):
                        # _\\|//_
                        #  //|\\
                        # ! ugly code warning !
                        if posi[1]>posx[1]
                            then
                            if posi[2]=posx[2]
                                then
                                oriented[i][1]:=neighbours[i][j];
                            elif posi[2]<posx[2]
                                then
                                if oriented[i][2]="pass"
                                    then # *always assigns the upper loop first*
                                    oriented[i][2]:=neighbours[i][j];
                                else
                                    oriented[i][3]:=neighbours[i][j];
                                fi;
                            elif posi[2]>posx[2]
                                then
                                if oriented[i][12]="pass"
                                    then
                                    oriented[i][12]:=neighbours[i][j];
                                else
                                    oriented[i][11]:=neighbours[i][j];
                                fi;
                            fi;
                        elif posi[1]=posx[1]
                            then
                            if posi[2]<posx[2]
                                then
                                oriented[i][4]:=neighbours[i][j];
                            elif posi[2]>posx[2]
                                then
                                oriented[i][10]:=neighbours[i][j];
                            fi;
                        elif posi[1]<posx[1]
                            then
                            if posi[2]=posx[2]
                                then
                                oriented[i][7]:=neighbours[i][j];
                            elif posi[2]<posx[2]
                                then
                                if oriented[i][5]="pass"
                                    then
                                    oriented[i][5]:=neighbours[i][j];
                                else
                                    oriented[i][6]:=neighbours[i][j];
                                fi;
                            elif posi[2]>posx[2]
                                then
                                if oriented[i][9]="pass"
                                    then
                                    oriented[i][9]:=neighbours[i][j];
                                else
                                    oriented[i][8]:=neighbours[i][j];
                                fi;
                            fi;
                        fi;
                    od;
                od;
                
                return oriented;
            end;

            return Clockwise(neighbours);
        end;

        path:=Orient(bound);
        # this is an ordered list of the neighbours of each 1-cell

        FaceTrace:=function(path)
            local
                unselectedEdges, sourceORtarget, faceloops,
                x, ClockwiseTurn, IsLoop, loop_correction, edge,
                2nd_loop, 2cell, sORt, ori, e1, e0, i;

            unselectedEdges:=List([1..Length(bound[2])-2]);
            unselectedEdges:=Concatenation(unselectedEdges,unselectedEdges);
            Add(unselectedEdges,Length(bound[2])-1);
            Add(unselectedEdges,Length(bound[2]));
# list of two of each edge except for the circumferential edges

            ClockwiseTurn:=function(p,e)
# inputs the orientation list of a node and the number of an edge in that list,
# outputs the next edge after a clockwise turn
                local
                    f;
                
                f:=(Position(p,e) mod 12)+1;
                while p[f]="pass"
                    do
                    f:=(f mod 12)+1;
                od;
                
                return p[f];
            end;

            ############ ADDED 15/10/19 ############
            IsLoop:=function(n)

                if Length(Positions(bound[2],bound[2][n]))=2 then
                    return true;
                else
                    return false;
                fi;

            end;

            loop_correction:=List(bound[2],x->0);
            ########################################

            sourceORtarget:=List([1..Length(bound[2])],y->[3,2]);
            x:=1;
            while unselectedEdges<>[]
                do # main loop, locates all 2-cells
                if rand
                    then
                    x:=Random([1..Length(bound[2])]); # select a random edge
                fi;
                while (not x in unselectedEdges) and (not e1 in unselectedEdges)
                    do # reselect edge if it already has two
                    if rand # 2-cells in its coboundary
                        then
                        x:=Random([1..Length(bound[2])]);
                    else
                        x:=x+1;
                    fi;
                od;

                2cell:=[x]; # the 2-cell begins with just x in its boundary
                if rand
                    then
                    sORt:=Random([2,3]);
                else
                    sORt:=sourceORtarget[x][Length(sourceORtarget[x])];
                    Unbind(sourceORtarget[x][Length(sourceORtarget[x])]);
                fi;
                ori:=path[bound[2][x][sORt]]; # the orientation of x's target
                e0:=bound[2][x][sORt];
                e1:=ClockwiseTurn(ori,x); # next edge to travel along
                while e1<>x
                    do
                    Add(2cell,e1);
                    e0:=Filtered(bound[2][e1]{[2,3]},y->y<>e0)[1]; # e1's target
                    ori:=path[e0];
                    e1:=ClockwiseTurn(ori,e1);
                od;
                Add(2cell,Length(2cell),1);
                if (not Set(2cell) in List(bound[3],x->Set(x)))
                    then
                    for i in Filtered(2cell{[2..Length(2cell)]},
                    y->y in unselectedEdges)
                        do
                        Unbind(unselectedEdges[Position(unselectedEdges,i)]);
                    od;
                    Add(bound[3],2cell);
                fi;

                ############ ADDED 15/10/19 ############
                # orders any loops that are present in the 2cell by the order
                # in which they were selected (doesn't include redundant 2cells
                # which are filtered out after the main while loop)
                if 2cell[1]<>2 then
                    faceloops:=Filtered(2cell{[2..Length(2cell)]},IsLoop);
                    if faceloops<>[] then
                        for edge in faceloops do
                            if loop_correction[edge]=0 then
                                loop_correction[edge]:=1;
                                2nd_loop:=Filtered(
                                    Positions(
                                        bound[2],
                                        bound[2][edge]
                                    ),
                                    y->y<>edge
                                )[1];
                                loop_correction[2nd_loop]:=2;
                            fi;
                        od;
                    fi;
                fi;
                ########################################
                        
            od;

            bound[3]:=Filtered(bound[3],y->y[1]<>2);
            return [bound,loop_correction];
        end;

        cgrid:=grid*0; # this is needed at the very end when
        for i in [1..Length(grid)] # patching the tubes together
            do
            for j in [1..Length(grid)]
                do
                cgrid[i][j]:=CornerConfiguration(i,j);
            od;
        od;

        traced_bound:=FaceTrace(path);

        return [traced_bound[1],cgrid,traced_bound[2]];
    end;

    P:=PuncturedDisk(D);
    grid:=P[2];
    loop_correction:=P[3];
    P:=P[1];

    PuncturedTube:=function(bound)
        local
            l0, l1, l2, DuplicateDisk,
            JoinDisks, Patch, prepatch,
            postpatch, Cap;

        l0:=Length(bound[1]);
        l1:=Length(bound[2]);
        l2:=Length(bound[3]);

        DuplicateDisk:=function(bound)
            local # creates a disjoint copy of the punctured
                i, edges2, faces2; # disk and concatenates the two

            for i in [1..l0]
                do
                Add(bound[1],[1,0]);
            od;

            edges2:=List(ShallowCopy(bound[2]),x->[2,x[2]+l0,x[3]+l0]);
            bound[2]:=Concatenation(bound[2],edges2);

            faces2:=List(ShallowCopy(bound[3]),
            x->Concatenation([x[1]],x{[2..Length(x)]}+l1));
            bound[3]:=Concatenation(bound[3],faces2);

            return bound;
        end;

        bound:=DuplicateDisk(bound);
        loop_correction:=Concatenation(loop_correction,loop_correction);

        JoinDisks:=function(bound)
# patch together the two disks via 1-cells, 2-cells & 3-cells
# mathematically speaking, form the space P x I where I is the unit interval
            local
                i, x, y, 3cell;

            for i in [1..l0]
                do # connect the 2 disks by 1-cells
                Add(bound[2],[2,i,l0+i]);
            od;

            for i in [1..l1]
                do # for each base 1-cell, form a 2-cell
                x:=bound[2][i][2];
                y:=bound[2][i][3];
                Add(bound[3],[4,i,l1+i,(l1*2)+x,(l1*2)+y]);
            od;

            for i in [1..l2]
                do # form a 3-cell from each 2-cell in the base disk
                x:=List(bound[3][i]{[2..Length(bound[3][i])]},y->y+(2*l2));
                3cell:=Concatenation([i,l2+i],x);
                Add(3cell,Length(3cell),1);
                Add(bound[4],3cell);
            od;

            return bound;
        end;

        bound:=JoinDisks(bound);

        prepatch:=Length(bound[3]);
        postpatch:=0;

        Patch:=function(bound)
            local # close the tubes to complete the construction
                loops, horizontals, verticals, i,
                lst, htube, h1, h2, vtube, x,
                cycle, loopless;

            loops:=Filtered(
                [1..l1-4],
                x->Length(Positions(bound[2],bound[2][x]))>1 and
                bound[2][x][2]<>1
            );

            horizontals:=Filtered(
                [1..l1-4],
                x->bound[2][x][2]=bound[2][x][3]-1
            );

            verticals:=Filtered(
                [1..l1-4],
                x->not (x in loops or x in horizontals)
            );
            verticals:=verticals+l1;

            for i in [1..Length(loops)/4]
                do
                lst:=[0,0];
                if 1 in grid[i]
                    then # check for corner configuration, important in deciding
                    lst[1]:=2; # which loop to use in the 2-cell (top/bottom)
                fi;
                if 2 in grid[i]
                    then
                    lst[2]:=4;
                fi;
                if 3 in grid[i]
                    then
                    lst[1]:=1;
                fi;
                if 4 in grid[i]
                    then
                    lst[2]:=3;
                fi;

                htube:=loops{lst+4*(i-1)};
                h1:=Filtered(
                    horizontals,
                    x->bound[2][x][2] in
                    [bound[2][htube[1]][2]..bound[2][htube[2]][2]]
                );
                h2:=Filtered(
                    horizontals,
                    x->bound[2][x][2] in
                    [bound[2][htube[1]][3]..bound[2][htube[2]][3]]
                );
                htube:=Concatenation(htube,h1,h2);
                Add(htube,Length(htube),1);
                Add(bound[3],htube);
            od;

            postpatch:=Length(bound[3]);

            loops:=loops+l1;

            vtube:=[];
            Add(vtube,verticals[1]);
            x:=bound[2][verticals[1]][3];
            cycle:=0;
            loopless:=[];
            for i in [2..Length(verticals)]
                do
                if bound[2][verticals[i]][2]=x
                    then
                    Add(vtube,verticals[i]);
                    x:=bound[2][verticals[i]][3];
                else
                    cycle:=cycle+1;
                    if cycle=2
                        then
                        cycle:=0;
                        Add(loopless,vtube);
                        vtube:=[];
                    fi;
                    x:=bound[2][verticals[i]][3];
                    Add(vtube,verticals[i]);
                fi;
            od;
            Add(loopless,vtube);

            for i in loopless
                do
                Add(i,Filtered(
                    loops,
                    y->bound[2][i[1]][2] in bound[2][y]{[2,3]})[1]
                );
                Add(i,Filtered(
                    loops,
                    y->bound[2][i[Length(i)-1]][3] in bound[2][y]{[2,3]})[2]
                );
                Add(i,Length(i),1);
                Add(bound[3],i);
            od;

            return bound;
        end;

        bound:=Patch(bound);

        Cap:=function(bound)
            local
                bottom, btm, top, tp,
                i, x, j, k;

            Add(bound[3],[2,l1-1,l1]); # the upper and lower
            Add(bound[3],[2,(2*l1)-1,2*l1]); # domes

            bottom:=[1..l2];
            btm:=[];
            for i in bound[3]{[prepatch+1..postpatch]}
                do
                x:=(Length(i)-3)/2;
                for j in [4..3+x]
                    do
                    for k in bottom
                        do
                        if 
                        i[j] in bound[3][k]{[2..Length(bound[3][k])]} and
                        i[j+x] in bound[3][k]{[2..Length(bound[3][k])]}
                            then
                            Add(btm,k);
                        fi;
                    od;
                od;
            od;

            bottom:=Difference(bottom,btm);
            bottom:=Concatenation(
                bottom, # all base 2-cells without the overlap
                [prepatch+1..postpatch], # the patch 2-cells enclosing the tubes
                [Length(bound[3])-1] # the dome
            );
            Add(bottom,Length(bottom),1);

            top:=[l2+1..2*l2];
            tp:=[];
            for i in bound[3]{[postpatch+1..Length(bound[3])-2]}
                do
                x:=(Length(i)-3)/2;
                for j in [2..1+x]
                    do
                    for k in top
                        do
                        if
                        i[j] in bound[3][k]{[2..Length(bound[3][k])]} and
                        i[j+x] in bound[3][k]{[2..Length(bound[3][k])]}
                            then
                            Add(tp,k);
                        fi;
                    od;
                od;
            od;

            top:=Difference(top,tp);
            top:=Concatenation(
                top,
                [postpatch+1..Length(bound[3])-2],
                [Length(bound[3])]
            );
            Add(top,Length(top),1);

            Add(bound[4],bottom);
            Add(bound[4],top);

            return bound;
        end;

        return Cap(bound);
    end;

    P:=PuncturedTube(P);
    P:=RegularCWComplex(P);

    for i in [1..Length(P!.boundaries[2])-Length(loop_correction)] do
        Add(loop_correction,0);
    od;

    P!.loopCorrection:=loop_correction;

    return P;
end);

################################################################################
################### An inclusion map for regular CW-structures on knots ########
################################################################################
############ Input: an arc presentation of some knot. ##########################
################################################################################
########### Output: an inclusion map of regular CW-complexes: from the #########
################### boundary of some knot to its complement. ###################
################################################################################
InstallGlobalFunction(
    KnotComplementWithBoundary,
    function(arc)
    local
        comp, RegularCWKnot, knot, hcorrection,
        threshold, inclusion, iota, inv_mapping;

    comp:=KnotComplement(arc);

    RegularCWKnot:=function(arc)
        local
            D, len, signless, HollowTubes, max,
            bigGrid, correction, threshold, hcorrection, TubeJoiner;

        D:=arc;
        len:=Length(D);
        signless:=List(D,x->[AbsInt(x[1]),AbsInt(x[2])]);

        HollowTubes:=function(D)
            local
                grid, i, IsIntersection,
                CornerConfiguration, bound,
                bigGrid, GridFill, j, tick, correction, hcorrection,
                hslice1, hslice2, l1, l2, l3,
                max, vslice1, vslice2, threshold;

            grid:=List([1..len],x->List([1..len],y->0));
            for i in [1..len]
                do
                grid[len-i+1][D[i][1]]:=1;
                grid[len-i+1][D[i][2]]:=1;
            od;

            IsIntersection:=function(i,j)

                if grid[i][j]=0
                    then
                    if 1 in grid[i]{[1..j]}
                        then
                        if 1 in grid[i]{[j..len]}
                            then
                            if 1 in List([1..i],x->grid[x][j])
                                then
                                if 1 in List([i..len],x->grid[x][j])
                                    then
                                    return true;
                                fi;
                            fi;
                        fi;
                    fi;
                fi;

            return false;
            end;

            CornerConfiguration:=function(i,j);

                if grid[i][j]=1
                    then
                    if Size(Positions(grid[i]{[j..len]},1))=2
                        then
                        if Size(Positions(List([i..len],x->grid[x][j]),1))=2
                            then # Corner type 1, i.e : __
                            return 1; #                |
                        elif Size(Positions(List([1..i],x->grid[x][j]),1))=2
                            then # Corner type 3, i.e :
                            return 3; #                |__
                        fi;
                    elif Size(Positions(grid[i]{[1..j]},1))=2
                        then
                        if Size(Positions(List([i..len],x->grid[x][j]),1))=2
                            then # Corner type 2, i.e : __
                            return 2; #                   |
                        elif Size(Positions(List([1..i],x->grid[x][j]),1))=2
                            then # Corner type 4, i.e :
                            return 4; #                 __|
                        fi;
                    fi;
                fi;

                return 0;
            end;

            bound:=[[],[],[],[],[]];
            bigGrid:=List([1..2*len],x->List([1..2*len],y->0));

            GridFill:=function(c,i,j);
                if c=1 or c=4
                    then
                    bigGrid[(2*i)-1][(2*j)-1]:='*';
                    bigGrid[2*i][2*j]:='*';
                elif c=2 or c=3
                    then
                    bigGrid[(2*i)-1][2*j]:='*';
                    bigGrid[2*i][(2*j)-1]:='*';
                fi;
            end;

            for i in [1..len]
                do
                for j in [1..len]
                    do
                    if IsIntersection(i,j)
                        then
                        bigGrid[(2*i)-1][(2*j)-1]:='*';
                        bigGrid[(2*i)-1][2*j]:='*';
                        bigGrid[2*i][(2*j)-1]:='*';
                        bigGrid[2*i][2*j]:='*';
                    elif grid[i][j]=1
                        then
                        GridFill(CornerConfiguration(i,j),i,j);
                    fi;
                od;
            od;

            tick:=1;
            for i in [1..2*len]
                do
                for j in [1..2*len]
                    do
                    if bigGrid[i][j]='*'
                        then
                        bigGrid[i][j]:=tick;
                        tick:=tick+1;
                    fi;
                od;
            od;

    # UPDATE: needed to account for configuration of corners at the end-step
    # (i.e. when matching the loops of one layer to the other).
    # There are sometimes disparities in the ordering on 1-cells from
    # left-to-right vs. when ordering from top-to-bottom. Not realising this was
    # causing 111 of the pre-stored knots in HAP to yield incorrect
    # CW-decomposition.
            correction:=[];
            hcorrection:=[];
            for i in [1..len] do
                for j in [1..len] do
                    if CornerConfiguration(i,j)<>0 then
                        if CornerConfiguration(i,j) in [1,4] then
                            Add(correction,1);
                            Add(correction,-1);
                            Add(hcorrection,2);
                            Add(hcorrection,1);
                        else
                            Add(correction,0);
                            Add(correction,0);
                            Add(hcorrection,1);
                            Add(hcorrection,2);
                        fi;
                    fi;
                od;
            od;

            ### add the 0, 1 & 2-cells ###
            ########## to bound ##########
            for i in [1..2*Maximum(bigGrid[Length(bigGrid)])] do
                Add(bound[1],[1,0]);
            od;

            for i in [1..len] do # add the 'horizontal' 2-cells
                hslice1:=Filtered(bigGrid[2*i-1],x->x<>0);
                hslice2:=Filtered(bigGrid[2*i],x->x<>0);
                l2:=Length(bound[2]);
                for j in [1..Length(hslice1)-1] do
                    Add(
                        bound[2],
                        [2,hslice1[j],hslice2[j]]
                    );
                    if j=1 then
                        Add(
                            bound[2],
                            [2,hslice1[j],hslice2[j]]
                        );
                    fi;
                    if j<>1 then
                        l1:=Length(bound[2]);
                        Add(
                            bound[3],
                            [4,l1-3,l1-2,l1-1,l1]
                        );
                    fi;
                    Add(
                        bound[2],
                        Concatenation([2],hslice1{[j,j+1]})
                    );
                    Add(
                        bound[2],
                        Concatenation([2],hslice2{[j,j+1]})
                    );
                    if j=Length(hslice1)-1 then
                        Add(
                            bound[2],
                            [2,hslice1[j+1],hslice2[j+1]]
                        );
                        l1:=Length(bound[2]);
                        Add(
                            bound[2],
                            [2,hslice1[j+1],hslice2[j+1]]
                        );
                        Add(
                            bound[3],
                            [4,l1-3,l1-2,l1-1,l1]
                        );
                    fi;
                od;
                l3:=Concatenation(
                    [l2+1],
                    Filtered(
                        [l2+3..Length(bound[2])-2],
                        x->AbsInt(bound[2][x][2]-bound[2][x][3])=1
                    ),
                    [Length(bound[2])]
                );
                Add(l3,Length(l3),1);
                Add(bound[3],l3);
            od;

            max:=Maximum(bigGrid[Length(bigGrid)]);
            for i in [1..2*len] do
                for j in [1..2*len] do
                    if bigGrid[i][j]<>0 then
                        bigGrid[i][j]:=bigGrid[i][j]+max;
                    fi;
                od;
            od;

            for i in [1..len] do # add the 'vertical' 2-cells
                vslice1:=Filtered(List([1..2*len],x->bigGrid[x][2*i-1]),x->x<>0);
                vslice2:=Filtered(List([1..2*len],x->bigGrid[x][2*i]),x->x<>0);
                l2:=Length(bound[2]);
                for j in [1..Length(vslice1)-1] do
                    Add(
                        bound[2],
                        Concatenation(
                            [2],
                            Set(
                                [vslice1[j],
                                vslice2[j]]
                            )
                        )
                    );
                    if j=1 then
                        Add(
                            bound[2],
                            Concatenation(
                                [2],
                                Set(
                                    [vslice1[j],
                                    vslice2[j]]
                                )
                            )
                        );
                    fi;
                    if j<>1 then
                        l1:=Length(bound[2]);
                        Add(
                            bound[3],
                            [4,l1-3,l1-2,l1-1,l1]
                        );
                    fi;
                    Add(
                        bound[2],
                        Concatenation([2],Set(vslice1{[j,j+1]}))
                    );
                    Add(
                        bound[2],
                        Concatenation([2],Set(vslice2{[j,j+1]}))
                    );
                    if j=Length(vslice1)-1 then
                        Add(
                            bound[2],
                            Concatenation(
                                [2],
                                Set(
                                    [vslice1[j+1],
                                    vslice2[j+1]]
                                )
                            )
                        );
                        l1:=Length(bound[2]);
                        Add(
                            bound[2],
                            Concatenation(
                                [2],
                                Set(
                                    [vslice1[j+1],
                                    vslice2[j+1]]
                                )
                            )
                        );
                        Add(
                            bound[3],
                            [4,l1-3,l1-2,l1-1,l1]
                        );
                    fi;
                od;
                l3:=Concatenation(
                    [l2+1],
                    Filtered(
                        [l2+3..Length(bound[2])-2],
                        x-> not
                        (bound[2][x][2] in vslice1 and bound[2][x][3] in vslice2)
                        and
                        not
                        (bound[2][x][2] in vslice2 and bound[2][x][3] in vslice1)
                    ),
                    [Length(bound[2])]
                );
                Add(l3,Length(l3),1);
                Add(bound[3],l3);
            od;

            threshold:=Length(bound[2])/2;
            ##############################

            return [max,bound,bigGrid,correction,threshold,hcorrection];
        end;

        D:=HollowTubes(D);
        max:=D[1];
        bigGrid:=D[3];
        correction:=D[4];
        threshold:=D[5];
        hcorrection:=D[6];
        D:=D[2];

        TubeJoiner:=function(D)
            local
                loops, size, hloops, vloops,
                i, l;

            loops:=Filtered([1..Length(D[2])],x->Length(Positions(D[2],D[2][x]))=2);
            size:=Length(loops)/2;
            hloops:=loops{[1..size]};
            vloops:=List([1..size],x->Position(D[2],D[2][hloops[x]]+[0,max,max]));
            for i in [1..size] do
                if i mod 2 = 0 then
                    vloops[i]:=vloops[i]+1;
                fi;
            od;
            vloops:=vloops+correction;

            for i in [1..size/2] do
                Add(D[2],[2,D[2][hloops[2*i]][2],D[2][vloops[2*i]][2]]);
                Add(D[2],[2,D[2][hloops[2*i]][3],D[2][vloops[2*i]][3]]);
                l:=Length(D[2]);

                Add(D[3],[4,hloops[2*i-1],vloops[2*i-1],l-1,l]);
                Add(D[3],[4,hloops[2*i],vloops[2*i],l-1,l]);
            od;

            return D;
        end;

        D:=TubeJoiner(D);
        D:=RegularCWComplex(D);
        D!.grid:=bigGrid;
        D!.arcPresentation:=arc;

        return [D,hcorrection,threshold];
    end;

    knot:=RegularCWKnot(arc);
    hcorrection:=knot[2];
    threshold:=knot[3];
    knot:=knot[1];

    inclusion:=function(bound)
        local
            bound1, bound2, len,
            1c1, 1c2, 2c2, HorizontalIndex, inc;

        bound1:=bound[1];
        bound2:=bound[2];

        len:=Length(bound1[1])/2;

        1c1:=bound1[2]*1;
        1c1:=List(
            1c1,
            x->Concatenation(
                [2],
                [x[2]+1+2*Int((x[2]-1)/len),
                x[3]+1+2*Int((x[3]-1)/len)]
            )
        );

        1c2:=bound2[2]*1;
        1c2:=List(
            1c2,
            x->Concatenation(
                [2],
                Set(
                    [x[2],x[3]]
                )
            )
        );

        2c2:=List(bound2[3],x->Concatenation([x[1]],Set(x{[2..Length(x)]})));

        HorizontalIndex:=function(n)

            return hcorrection[
                Position(
                    Filtered(
                        [1..threshold],
                        x->Length(
                            Positions(
                                bound1[2],
                                bound1[2][x]
                            )
                        )=2
                    ),n
                )
            ];
        end;

        inc:=function(n,k)
            local
                ind, 2cell;

            if n=0 then
                return k+1+2*Int((k-1)/len);
            elif n=1 then
                ind:=1;
                if k in [1..threshold] then
                    if k>1 and 1c1[k-1]=1c1[k] then
                        ind:=HorizontalIndex(k);
                    elif 1c1[k]=1c1[k+1] then
                        ind:=HorizontalIndex(k);
                    fi;
                elif k>threshold then
                    if 1c1[k-1]=1c1[k] then
                        ind:=2;
                    elif k<Length(1c1) and 1c1[k]=1c1[k+1] then
                        ind:=1;
                    fi;
                fi;
                return Positions(1c2,1c1[k])[ind];
            elif n=2 then
                2cell:=List(bound1[3][k]{[2..Length(bound1[3][k])]},x->inc(1,x));
                2cell:=Concatenation([Length(2cell)],Set(2cell));
                return Position(2c2,2cell);
            else
                return fail;
            fi;
        end;

        return inc;
    end;

    iota:=inclusion([knot!.boundaries,comp!.boundaries]);

    inv_mapping:=[[],[],[]];
    inv_mapping[1]:=List([1..knot!.nrCells(0)],x->iota(0,x));
    inv_mapping[2]:=List([1..knot!.nrCells(1)],x->iota(1,x));
    inv_mapping[3]:=List([1..knot!.nrCells(2)],x->iota(2,x));

    return Objectify(
        HapRegularCWMap,
        rec(
            source:=knot,
            target:=comp,
            mapping:=iota,
            properties:=[
                ["image",inv_mapping]
            ]
        )
    );
end);
################################################################################