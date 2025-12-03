
# KnotComplement
################################################################################
############ Input: an arc presentation ########################################
################################################################################
########### Output: a regular CW-complex representing the complement of ########
################### the input link #############################################
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

# KnotComplementWithBoundary
################################################################################
############ Input: an arc presentation of some link K #########################
################################################################################
########### Output: an inclusion of CW-complexes f: b(K) -> B^3 \ K ############
################### where B^3 is homeomorphic to the 3-ball and b(K) ###########
################### denotes the boundary of an open tubular neighbourhood ######
################### of K in B^3 ################################################
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

# ArcDiagramToTubularSurface
################################################################################
############ Input: a list [arc,crs,col] where arc is an arc presentation, #####
################### crs is a list whose entries are -1, 0 or 1 and whose #######
################### length is the number of crossings in arc and col is ########
################### a list whose entries are 1, 2, 3 or 4 corresponding to #####
################### colours as detailed below ##################################
################################################################################
########### Output: an inclusion of regular CW-complexes from a self- ##########
################### intersecting tubular surface into the 3-ball ###############
################################################################################
InstallGlobalFunction(
    ArcDiagramToTubularSurface,
    function(arc)
    local
        prs, crs, clr, grd, i, IsIntersection,
        CornerConfiguration, GRD, crossings, j,
        k, nr0cells, bnd, sub, hbars, hslice, cell,
        int, max, vbars, vslice, loops1, loops2,
        unchosen, neighbours, Clockwise,
        path, unselectedEdges, sourceORtarget,
        faceloops, x, ClockwiseTurn, 2cell, sORt,
        ori, e1, e0, loops, present_loops, vertices,
        check, l0, l1, l2, l1_, l2_, IsEdgeInDuplicate,
        copy1, hbars2, vbars2, copy2, 3cell, colour, lcap, reg,
        closure, ucap, l1__, l2__, path_comp, pipes, HorizontalOrVertical,
        l, AboveBelow0Cell, IntersectingCylinders, pos, colour_,Last;
Last:=function(L); return L[Length(L)]; end;
    if IsList(arc[1][1]) then
        prs:=arc[1]*1;
        crs:=arc[2]*1;
#                -1 |     +1 |      0 |
# crossing types: --|-- or ----- or --+--
#                   |        |        |
        if Length(arc)=3 then
            clr:=arc[3]*1;
# colours: 1; bgb, 2; bgr, 3; rgb, 4; rgr
        fi;
    else
        prs:=arc*1;
    fi;

# (i) the 0-skeleton of the disk
####################################################################################
    grd:=List([1..Length(prs)],x->[1..Length(prs)]*0);
    for i in [0..Length(prs)-1] do
        grd[Length(prs)-i][prs[i+1][1]]:=1;
        grd[Length(prs)-i][prs[i+1][2]]:=1;
    od;

    IsIntersection:=function(i,j)
        if grd[i][j]=0 and
            1 in grd[i]{[1..j]} and
                1 in grd[i]{[j..Length(prs)]} and
                    1 in List([1..i],x->grd[x][j]) and
                        1 in List([i..Length(prs)],x->grd[x][j]) then
            return true;
        fi;
        return false;
    end;
    CornerConfiguration:=function(i,j);
        if grd[i][j]=1 then
            if Size(Positions(grd[i]{[j..Length(prs)]},1))=2 then
                if Size(Positions(List([i..Length(prs)],x->grd[x][j]),1))=2 then
                        # Corner type 1, i.e :  __
                    return 1; #                |
                elif Size(Positions(List([1..i],x->grd[x][j]),1))=2 then
                        # Corner type 3, i.e :
                    return 3; #                |__
                fi;
            elif Size(Positions(grd[i]{[1..j]},1))=2 then
                if Size(Positions(List([i..Length(prs)],x->grd[x][j]),1))=2 then
                        # Corner type 2, i.e :  __
                    return 2; #                   |
                elif Size(Positions(List([1..i],x->grd[x][j]),1))=2 then
                        # Corner type 4, i.e :
                    return 4; #                 __|
                fi;
            fi;
        fi;
        return 0;
    end;
    
    GRD:=List([1..4*Length(prs)],x->[1..4*Length(prs)]*0);
    crossings:=[];
# quadruple the size of grd to allow for the 0-skeleton
# to be displayed nicely without overlap in an array
    for i in [1..Length(prs)] do
        for j in [1..Length(prs)] do
            if CornerConfiguration(i,j) in [1,4] then
                GRD[4*i-3][4*j-3]:=1;
                GRD[4*i][4*j]:=1;
            elif CornerConfiguration(i,j) in [2,3] then
                GRD[4*i-3][4*j]:=1;
                GRD[4*i][4*j-3]:=1;
            elif IsIntersection(i,j) then
                for k in [0,3] do
                    GRD[4*i-3][4*j-3+k]:=1;
                    GRD[4*i][4*j-3+k]:=1;
                    Add(crossings,[4*i-3,4*j-3+k]);
                    Add(crossings,[4*i,4*j-3+k]);
                od;
            fi;
        od;
    od;
# label the 0-cells row by row
    nr0cells:=2;
    for i in [1..4*Length(prs)] do
        for j in [1..4*Length(prs)] do
            if GRD[i][j]=1 then
                GRD[i][j]:=nr0cells;
                nr0cells:=nr0cells+1;
            fi;
        od;
    od;
    crossings:=List(crossings,x->GRD[x[1]][x[2]]);
    crossings:=List([1..Length(crossings)/4],x->List([1..4]+4*x-4,y->crossings[y]));
    GRD:=FrameArray(GRD);
    GRD[1][1]:=1;
    GRD[4*Length(prs)+2][4*Length(prs)+2]:=nr0cells;

    bnd:=List([1..5],x->[]); # eventual boundary list of the 3-ball containing
    sub:=List([[],[],[]]); # the boundary of a knotted surface as a subcomplex
    bnd[1]:=List([1..nr0cells],x->[1,0]);
    sub[1]:=[2..nr0cells-1];
    if IsBound(crs) then
        for i in [1..Length(crs)] do
            if crs[i]=0 then
                sub[1]:=Difference(sub[1],crossings[i]);
            fi;
        od;
    fi;
####################################################################################

# (ii) the 1-skeleton of the disk
####################################################################################
    # add the horizontal 1-cells
    hbars:=[];
    for i in [2..Length(GRD)-1] do
        hslice:=Filtered(GRD[i],x->x<>0);
        if hslice<>hslice*0 then
            Add(hbars,hslice);
        fi;
        for j in [1..Length(hslice)-1] do
            cell:=[2,hslice[j],hslice[j+1]];
            Add(bnd[2],cell);
            int:=List(crossings,x->Length(Intersection(cell{[2,3]},x)));
            max:=PositionMaximum(int);
            if IsBound(crs) then
                if int[max]=1 and
                    crs[max]=-1 and
                        not Length(bnd[2]) in sub[2] then
                            Add(sub[2],Length(bnd[2]));
                elif int[max]=2 and
                    crs[max]<>0 and
                        not Length(bnd[2]) in sub[2] then
                            Add(sub[2],Length(bnd[2])); 
                fi;
            else
                if max<>fail then
                    if int[max]=2 then #  horizontal 1-cells default to the top
                        Add(sub[2],Length(bnd[2])); # so they're not all included here
                    fi;
                fi;
            fi;
        od;
    od;
    hbars:=List([1..Length(hbars)/2],x->Concatenation(hbars[2*x-1],hbars[2*x]));

    # add the vertical 1-cells
    vbars:=[];
    for i in TransposedMat(GRD){[2..Length(GRD)-1]} do
        vslice:=Filtered(i,x->x<>0);
        if vslice<>vslice*0 then
            Add(vbars,vslice);
        fi;
        for j in [1..Length(vslice)-1] do
            cell:=[2,vslice[j],vslice[j+1]];
            Add(bnd[2],cell);
            int:=List(crossings,x->Length(Intersection(cell{[2,3]},x)));
            max:=PositionMaximum(int);
            if IsBound(crs) then
                if int[max]=1 then
                    if crs[max]=1 and
                        not Length(bnd[2]) in sub[2] then
                            Add(sub[2],Length(bnd[2]));
                    fi;
                elif int[max]=2 then
                    if crs[max]<>0 and
                        not Length(bnd[2]) in sub[2] then
                            Add(sub[2],Length(bnd[2]));
                    fi;
                else
                    Add(sub[2],Length(bnd[2]));
                fi;
            else
                Add(sub[2],Length(bnd[2])); # vertical 1-cells default to the bottom
            fi;
        od;
    od;
    vbars:=List([1..Length(vbars)/2],x->Concatenation(vbars[2*x-1],vbars[2*x]));

    # add the loops
    for i in [2,6..Length(GRD)-4] do
        loops1:=Filtered(GRD[i],x->x<>0);
        loops2:=Filtered(GRD[i+3],x->x<>0);
        for j in [1,2] do
            Add(bnd[2],[2,loops1[1],loops2[1]]);
            Add(sub[2],Length(bnd[2])); # loops always in subcomplex... for now
            Add(bnd[2],[2,loops1[Length(loops1)],loops2[Length(loops2)]]);
            Add(sub[2],Length(bnd[2]));
        od;
    od;

    # add the remaining four 1-cells to keep things regular
    Add(bnd[2],[2,1,2]); Add(bnd[2],[2,nr0cells-1,nr0cells]);
    Add(bnd[2],[2,1,nr0cells]); Add(bnd[2],[2,1,nr0cells]);
####################################################################################

# (iii) the 2-skeleton of the disk
####################################################################################
    unchosen:=List(bnd[2],x->x{[2,3]});
    neighbours:=List(bnd[1],x->[]);

    for i in [1..Length(bnd[1])] do
        for j in [1..Length(unchosen)] do
            if i in unchosen[j] then
                Add(neighbours[i],j);
            fi;
        od;
    od;

    Clockwise:=function(neighbours)
        local
            oriented, first0, last0,
            i, j, x, k, l, posi, posx;

        oriented:=List(neighbours,x->List([1..12],y->"pass"));
        first0:=SortedList(neighbours[1]);
        last0:=SortedList(neighbours[Length(neighbours)]);

        oriented[1][7]:=first0[1];
        oriented[1][6]:=first0[3];
        oriented[1][8]:=first0[2];
        oriented[Length(oriented)][1]:=last0[1]; 
        oriented[Length(oriented)][2]:=last0[3];
        oriented[Length(oriented)][12]:=last0[2];

        for i in [2..Length(neighbours)-1] do
            for j in [1..Length(neighbours[i])] do
                x:=bnd[2][neighbours[i][j]];
                x:=Filtered(x{[2,3]},y->y<>i)[1];
                for k in [1..Length(GRD)] do
                    for l in [1..Length(GRD[1])] do
                        if i=GRD[k][l] then
                            posi:=[k,l];
                        fi;
                        if x=GRD[k][l]then
                            posx:=[k,l];
                        fi;
                    od;
                od;
                # _\\|//_
                #  //|\\
                if posi[1]>posx[1] then
                    if posi[2]=posx[2] then
                        oriented[i][1]:=neighbours[i][j];
                    elif posi[2]<posx[2] then
                        if oriented[i][2]="pass" then
                            oriented[i][2]:=neighbours[i][j];
                        else
                            oriented[i][3]:=neighbours[i][j];
                        fi;
                    elif posi[2]>posx[2] then
                        if oriented[i][12]="pass" then
                            oriented[i][12]:=neighbours[i][j];
                        else
                            oriented[i][11]:=neighbours[i][j];
                        fi;
                    fi;
                elif posi[1]=posx[1] then
                    if posi[2]<posx[2] then
                        oriented[i][4]:=neighbours[i][j];
                    elif posi[2]>posx[2] then
                        oriented[i][10]:=neighbours[i][j];
                    fi;
                elif posi[1]<posx[1] then
                    if posi[2]=posx[2] then
                        oriented[i][7]:=neighbours[i][j];
                    elif posi[2]<posx[2] then
                        if oriented[i][5]="pass" then
                            oriented[i][5]:=neighbours[i][j];
                        else
                            oriented[i][6]:=neighbours[i][j];
                        fi;
                    elif posi[2]>posx[2] then
                        if oriented[i][9]="pass" then
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

    path:=Clockwise(neighbours);

    unselectedEdges:=List([1..Length(bnd[2])-2]);
    unselectedEdges:=Concatenation(unselectedEdges,unselectedEdges); # use Append() instead
    Add(unselectedEdges,Length(bnd[2])-1);
    Add(unselectedEdges,Length(bnd[2]));

    ClockwiseTurn:=function(p,e)
        local f;
        f:=(Position(p,e) mod 12)+1;
        while p[f]="pass" do
            f:=(f mod 12)+1;
        od;
        return p[f];
    end;

    sourceORtarget:=List([1..Length(bnd[2])],y->[3,2]);
    x:=1;
    while unselectedEdges<>[] do
        while (not x in unselectedEdges) and (not e1 in unselectedEdges) do
            x:=x+1;
        od;
        2cell:=[x];
        sORt:=sourceORtarget[x][Length(sourceORtarget[x])];
        Unbind(sourceORtarget[x][Length(sourceORtarget[x])]);

        ori:=path[bnd[2][x][sORt]];
        e0:=bnd[2][x][sORt];
        e1:=ClockwiseTurn(ori,x);
        while e1<>x do
            Add(2cell,e1);
            e0:=Filtered(bnd[2][e1]{[2,3]},y->y<>e0)[1];
            ori:=path[e0];
            e1:=ClockwiseTurn(ori,e1);
        od;
        Add(2cell,Length(2cell),1);
        if not Set(2cell) in List(bnd[3],Set) then
            for i in Filtered(2cell{[2..Length(2cell)]},y->y in unselectedEdges) do
                Unbind(unselectedEdges[Position(unselectedEdges,i)]);
            od;
            Add(bnd[3],2cell);
            if Length(Intersection(2cell{[2..2cell[1]+1]},sub[2]))=2cell[1] and
                2cell[1]<>2 then
                    Add(sub[3],Length(bnd[3]));
            fi;
        fi;
    od;
    bnd[3]:=List(bnd[3],x->Concatenation([x[1]],Set(x{[2..Length(x)]})));
####################################################################################

# (iv) the direct product of the disk with [0,1]
####################################################################################
# make a duplicate of the disk
    l0:=Length(bnd[1]);
    l1:=Length(bnd[2]); l1_:=Length(sub[2]);
    l2:=Length(bnd[3]); l2_:=Length(sub[3]);

    IsEdgeInDuplicate:=function(k)
        local x;

        if Length(Positions(bnd[2],bnd[2][k]))=2 and
            bnd[2][k]<>[2,l0+1,2*l0] then
                return true;
        elif not Position(bnd[2],bnd[2][k]-[0,l0,l0]) in sub[2] and
            not bnd[2][k]-[0,l0,l0] in [[2,1,2],[2,1,l0],[2,l0-1,l0]] then
                return true;
        else
            x:=List(
                List(crossings,y->y+l0),
                z->Length(Intersection(z,bnd[2][k]{[2,3]}))=2
            );
            if true in x then
                if IsBound(crs) then
                    if crs[Position(x,true)] in [1,-1] then
                        return true;
                    fi;   
                else
                    return true;
                fi;
            fi;                 
        fi;
        return false;
    end;

    loops:=Filtered(bnd[2]*1,x->Length(Positions(bnd[2],x))=2);
    loops:=Set(Concatenation(List(loops,x->x{[2,3]})));
    bnd[1]:=Concatenation(bnd[1],bnd[1]);
    sub[1]:=Concatenation(sub[1],[l0+2..2*l0-1]); 
# sub contains all duplicate 0-cells except for those in the frame... for now

    copy1:=List(bnd[2],x->x+[0,l0,l0]);
    bnd[2]:=Concatenation(bnd[2],copy1);
    for i in [l1+1..2*l1] do
        if IsEdgeInDuplicate(i) then
            Add(sub[2],i);
        fi;
    od;

    hbars2:=List(hbars,x->x+l0);
    vbars2:=List(vbars,x->x+l0);
    copy2:=List(bnd[3],x->Concatenation([x[1]],List(x{[2..Length(x)]},y->y+l1)));
    bnd[3]:=Concatenation(bnd[3],copy2);
    for i in [l2..2*l2] do
        cell:=bnd[3][i]*1;
        x:=[]; # all 0-cells in a given 2-cell
        for j in [2..Length(cell)] do
            for k in [2,3] do
                Add(x,bnd[2][cell[j]][k]);
            od;
        od;
        x:=Set(x);
        if Length(Intersection(cell{[2..cell[1]+1]},sub[2]))=cell[1] and
            (
                true in List(hbars2,y->IsSubset(y,x)) or
                true in List(vbars2,y->IsSubset(y,x))
            ) and
                cell[1]<>2 then
                    Add(sub[3],i);
        fi;
    od;
    l1__:=Length(sub[2]); l2__:=Length(sub[3]);

# 13-02-21 remove any path components of sub that consist of two 0-cells and two 1-cells only
    path_comp:=PathComponentsCWSubcomplex([RegularCWComplex(bnd),sub]);
    for i in [1..Length(path_comp)] do
        if Length(path_comp[i][2][1])=2 and
            Length(path_comp[i][2][2])=2 then
                for j in [1,2] do
                    sub[j]:=Difference(sub[j],path_comp[i][2][j]);
                od;
        fi;
    od;

# join the original disk to the copy
# each n-cell of the disk yields an (n+1)-cell which connects it
# to its duplicate
    for i in [1..l0] do
        Add(bnd[2],[2,i,i+l0]);
        if i in loops and
            i in sub[1] and
                i+l0 in sub[1] and
                    not i in [1,l0] then
                        Add(sub[2],Length(bnd[2]));
        fi;
    od;

    for i in [1..l1] do
        Add(
            bnd[3],
            [
                4,
                i,
                Position(bnd[2],[2,bnd[2][i][2],bnd[2][i+l1][2]]),
                Position(bnd[2],[2,bnd[2][i][3],bnd[2][i+l1][3]]),
                i+l1
            ]
        );
        if i in sub[2] and
            Length(Positions(bnd[2],bnd[2][i]))=2 and
                i+l1 in sub[2] then
                Add(sub[3],Length(bnd[3]));
        fi;       
    od;

    for i in [1..l2] do
        3cell:=[];
        # 13-02-21 plug the corner holes
            if bnd[3][i][1]=2 and
                bnd[3][i][2] in sub[2] and
                    not bnd[3][i][2]+l1 in sub[2] then
                        Add(sub[3],i);
            fi;
            if bnd[3][i][1]=2 and
                not bnd[3][i][2] in sub[2] and
                    bnd[3][i][2]+l1 in sub[2] then
                        Add(sub[3],i+l2);
            fi;
        for j in bnd[3][i]{[2..Length(bnd[3][i])]} do
            Add(
                3cell,
                Position(
                    bnd[3],
                    [
                        4,
                        j,
                        Position(bnd[2],[2,bnd[2][j][2],bnd[2][j+l1][2]]),
                        Position(bnd[2],[2,bnd[2][j][3],bnd[2][j+l1][3]]),
                        j+l1
                    ]
                )
            );
        od;
        Add(3cell,i);
        Add(3cell,i+l2);
        Add(bnd[4],Concatenation([bnd[3][i][1]+2],3cell));
    od;

    for i in [3,4] do
        bnd[i]:=List(bnd[i],x->Concatenation([x[1]],Set(x{[2..Length(x)]})));
    od;
####################################################################################

# (v) join the open ended 'pipes' at either end of B^2xI according to crs
####################################################################################
    colour:=List([1..4],x->[]);

    if not IsBound(crs) then
        crs:=List([1..Length(crossings)],x->1);
    fi;

# (a) start with the lower caps, they're more straight forward #####################
    lcap:=[];
    reg:=RegularCWComplex(bnd);
    for i in [1..l2] do
        if bnd[3][i][1]<>2 and not i in sub[3] then
            Add(lcap,i);
        elif bnd[3][i][1]=2 and not i in sub[3] then
            closure:=ClosureCWCell(reg,2,i);
            if not IsSubset(sub[2],closure[2][2]) then
                Add(lcap,i);
            fi;
        fi;
    od;

    pipes:=[]; # the 0,1 and 2-cells of each 'pipe' in the lower dome which will join
    # together to form the (intersecting) tubular surface
    for i in [1..l2] do
        if i in sub[3] then
            Add(
                pipes,
                [
                    [i],
                    bnd[3][i]{[2..bnd[3][i][1]+1]},
                    Set(
                        Concatenation(
                            List(
                                bnd[3][i]{[2..bnd[3][i][1]+1]},
                                x->bnd[2][x]{[2,3]}
                            )
                        )
                    )
                ]
            );
        fi;
    od;
    pipes:=Filtered(pipes,x->Length(x[3])>2);
    for i in [1..Length(pipes)] do
        for j in [i+1..Length(pipes)] do
            if Intersection(
                Intersection(
                    pipes[i][3],pipes[j][3]
                ),
                Concatenation(crossings)
            )<>[] then
                Append(pipes[i][1],pipes[j][1]);
                Append(pipes[i][2],pipes[j][2]);
                Append(pipes[i][3],pipes[j][3]);
                pipes[j][3]:=[];
            fi;
        od;
    od;
    Apply(pipes,x->[Set(x[1]),Set(x[2]), Set(x[3])]);
    pipes:=Filtered(pipes,x->x[3]<>[]);
    for i in [1..Length(pipes)] do # either end of each pipe has two 2-cells which
    # are currently absent from pipe[i][1], this step will add them
        for j in Filtered(bnd[3],x->x[1]=2) do
            int:=Intersection(pipes[i][2],j{[2,3]});
            if int<>[] then
                Add(pipes[i][1],Position(bnd[3],j));
            fi;
        od;
    od;
    for i in [1..Length(pipes)] do # if two pipes' 0-skeletons intersect, add a
        for j in [i+1..Length(pipes)] do # new 1-cell whose boundary is their intersection
            int:=Intersection(pipes[i][3],pipes[j][3]);
            if int<>[] then
                Add(int,2,1);
                Add(bnd[2],int);
                Add(sub[2],Length(bnd[2]));
            fi;
        od;
    od;
    for i in [1..Length(pipes)] do # replace the 1-cells that occur twice with
        for j in [1..Length(pipes[i][2])] do # either the new cell we just added or the
        # other occurance of that cell which is not in pipes[i][2]
            pos:=Positions(
                bnd[2],
                [
                    2,
                    bnd[2][pipes[i][2][j]][2],
                    bnd[2][pipes[i][2][j]][3]
                ]
            );
            if Length(pos)>=2 then
                cell:=Last(Filtered(pos,x->x<>pipes[i][2][j]));
                pipes[i][2][j]:=cell;
            fi;
        od;
    od;
    HorizontalOrVertical:=function(l)
        if l in hbars then
            return "horizontal";
        fi;
        return "vertical";
    end;
    for i in [1..Length(pipes)] do # if a pipe contains 1-cells from crossings then
        # filter out the horizontal/vertical 1-cells
        # should said pipe be vertical/horizontal itself
        if Intersection(pipes[i][3],Concatenation(crossings))<>[] then
            ori:=HorizontalOrVertical(pipes[i][3]);
            l:=Filtered(crossings,x->Intersection(x,pipes[i][3])<>[])[1];
            if ori="horizontal" then
                pipes[i][2]:=Difference(
                    pipes[i][2],
                    [
                        Position(bnd[2],[2,l[1],l[2]]),
                        Position(bnd[2],[2,l[3],l[4]])
                    ]
                );
            else
                pipes[i][2]:=Difference(
                    pipes[i][2],
                    [
                        Position(bnd[2],[2,l[1],l[3]]),
                        Position(bnd[2],[2,l[2],l[4]])
                    ]
                );
            fi;
        fi;
    od;
    for i in [1..Length(pipes)] do # add new 2-cells to bnd[3] whose boundary is
    # the pipes[i][2] 1-cells
        cell:=Set(pipes[i][2]);
        Add(cell,Length(cell),1);
        Add(bnd[3],cell);
        Add(sub[3],Length(bnd[3]));
        Add(lcap,Length(bnd[3]));
        Add(pipes[i][1],Length(bnd[3]));
    od;
    for i in [1..Length(pipes)] do # find where pipes intersect in their 2-skeleta
    # use this to form the 3-cells comprising the interiors of the pipes
        for j in [i+1..Length(pipes)] do
            if Intersection(pipes[i][1],pipes[j][1])<>[] then
                Append(pipes[i][1],pipes[j][1]);
                pipes[j][1]:=[];
            fi;
        od;
    od;
    for i in [1..Length(pipes)] do
        if pipes[i][1]<>[] then
            cell:=Set(pipes[i][1]);
            Add(cell,Length(cell),1);
            Add(bnd[4],cell);
        fi;
    od;

# (b) now for the upper caps, 0 in crs leads to a very elaborate CW-structure #######
# 13-02-21 this function was added to help alter the intersection structure below
# after discovering it was incorrect
    AboveBelow0Cell:=function(n,str)
        local n_, i, j, l;

        if n>l0 then
            n_:=n-l0;
        else
            n_:=n*1;
        fi;

        i:=Position(List(GRD,x->n_ in x),true);
        j:=Position(GRD[i],n_);
        
        if str="above" then
            l:=Reversed([1..i-1]);
        else
            l:=[i+1..Length(GRD)];
        fi;

        for k in l do
            if GRD[k][j]<>0 then
                if n>l0 then
                    return GRD[k][j]+l0;
                else
                    return GRD[k][j];
                fi;
            fi;
        od;
    end;
    
    ucap:=[];
    reg:=RegularCWComplex(bnd);
    for i in [l2+1..2*l2] do
        if bnd[3][i][1]<>2 and not i in sub[3] then
            Add(ucap,i);
        elif bnd[3][i][1]=2 and not i in sub[3] then
            closure:=ClosureCWCell(reg,2,i);
            if not IsSubset(sub[2],closure[2][2]) then
                Add(ucap,i);
            fi;
        fi;
    od;

    pipes:=[]; # the 0,1 and 2-cells of each 'pipe' in the upper dome which will join
    # together to form the (intersecting) tubular surface
    for i in [l2+1..2*l2] do
        if i in sub[3] then
            Add(
                pipes,
                [
                    [i],
                    bnd[3][i]{[2..bnd[3][i][1]+1]},
                    Set(
                        Concatenation(
                            List(
                                bnd[3][i]{[2..bnd[3][i][1]+1]},
                                x->bnd[2][x]{[2,3]}
                            )
                        )
                    )
                ]
            );
        fi;
    od;
    pipes:=Filtered(pipes,x->Length(x[3])>2);
    
    IntersectingCylinders:=function(a,b,c,d) ######################################
        local
            n, i, m, j, l, a_abv,
            b_abv, c_blw, d_blw,
            del_2_cell_1, del_2_cell_2,
            del_2_cell_1_i, del_2_cell_2_i,
            dif, f_clr, f_clr1, f_clr2;
# attaches to a 0 crossing some additional regular CW-structure
# to allow for a self-intersection to occur
        n:=1*Length(bnd[1])+1;
        for i in [1..8] do # 0-skeleton of intersection
            Add(bnd[1],[1,0]);
            Add(sub[1],Length(bnd[1]));
        od;
        # 1-skeleton of intersection
        m:=1*Length(bnd[2])+1;
        Add(bnd[2],[2,a,n]); Add(sub[2],Length(bnd[2])); # m
        Add(bnd[2],[2,b,n+1]); Add(sub[2],Length(bnd[2])); # m+1
        Add(bnd[2],[2,c,n+6]); Add(sub[2],Length(bnd[2])); # m+2
        Add(bnd[2],[2,d,n+7]); Add(sub[2],Length(bnd[2])); # m+3
        for i in [0..3] do
            for j in [1,2] do
                Add(bnd[2],[2,n+2*i,n+1+2*i]); # m+4, m+5, m+6, m+7, m+10, m+11, m+14, m+15
                # top first, then bottom (refer to drawing)
                Add(sub[2],Length(bnd[2]));
            od;
            if i>0 then
                Add(bnd[2],[2,n+2*i-2,n+2*i]); Add(sub[2],Length(bnd[2])); # m+8, m+12, m+16
                Add(bnd[2],[2,n+2*i-1,n+2*i+1]); Add(sub[2],Length(bnd[2])); # m+9, m+13, m+17
            fi;
        od;

        # 13-02-21 need to remove cells from sub to accommodate the
        # new self-intersection structure
        a_abv:=AboveBelow0Cell(a,"above");
        b_abv:=AboveBelow0Cell(b,"above");
        c_blw:=AboveBelow0Cell(c,"below");
        d_blw:=AboveBelow0Cell(d,"below");
        del_2_cell_1:=[
            Position(bnd[2],[2,a_abv,a]),
            Position(bnd[2],[2,b_abv,b]),
            Position(bnd[2],[2,a,b])
        ];
        del_2_cell_1_i:=Position(
            List(
                bnd[3],
                x->Length(
                    Intersection(
                        x{[2..x[1]+1]},
                        del_2_cell_1
                    )
                )=3
            ),
            true
        );
        del_2_cell_2:=[
            Position(bnd[2],[2,c,c_blw]),
            Position(bnd[2],[2,d,d_blw]),
            Position(bnd[2],[2,c,d])
        ];
        del_2_cell_2_i:=Position(
            List(
                bnd[3],
                x->Length(
                    Intersection(
                        x{[2..x[1]+1]},
                        del_2_cell_2
                    )
                )=3
            ),
            true
        );
        sub[3]:=Difference(sub[3],[del_2_cell_1_i,del_2_cell_2_i]);
        sub[2]:=Difference(
            sub[2],
            [
                Position(bnd[2],[2,a_abv,a]),
                Position(bnd[2],[2,b_abv,b]),
                Position(bnd[2],[2,c,c_blw]),
                Position(bnd[2],[2,d,d_blw])
            ]
        );

        # 14-02-21 this new structure has 4 additional 1-cells
        Add(bnd[2],[2,a_abv,n]); Add(sub[2],Length(bnd[2])); # m+18
        Add(bnd[2],[2,b_abv,n+1]); Add(sub[2],Length(bnd[2])); # m+19
        Add(bnd[2],[2,n+6,c_blw]); Add(sub[2],Length(bnd[2])); # m+20
        Add(bnd[2],[2,n+7,d_blw]); Add(sub[2],Length(bnd[2])); # m+21

        # 2-skeleton of intersection
        l:=1*Length(bnd[3])+1;
        Add( # l
            bnd[3],
            [
                4,
                Position(bnd[2],[2,a,b]),
                m,
                m+1,
                m+5
            ]
        );
        Add(sub[3],Length(bnd[3]));
        Add( # l+1
            bnd[3],
            [
                4,
                Position(bnd[2],[2,c,d]),
                m+2,
                m+3,
                m+15
            ]
        );
        Add(sub[3],Length(bnd[3]));
        # these 2-cells are those which should be coloured #########################
        if IsBound(clr) then                                                      ##
            pos:=Position(List(crossings,x->Set(x)+l0),Set([a,b,c,d]));           ##
            pos:=pos-Length(Filtered(crs{[1..pos-1]},x->x<>0));                   ##
        fi;                                                                       ##
        Add(bnd[3],[4,m+4,m+6,m+8,m+9]); # l+2                                    ##
        Add(sub[3],Length(bnd[3]));                                               ##              
        Add(bnd[3],[4,m+5,m+7,m+8,m+9]); # l+3                                    ##
        Add(sub[3],Length(bnd[3]));                                               ## 
        for i in [0,1] do                                                         ##
            Add(bnd[3],[4,m+6+4*i,m+10+4*i,m+12+4*i,m+13+4*i]); # l+4, l+6        ##  
            Add(sub[3],Length(bnd[3]));                                           ##
            Add(bnd[3],[4,m+7+4*i,m+11+4*i,m+12+4*i,m+13+4*i]); # l+5, l+7        ##
            Add(sub[3],Length(bnd[3]));                                           ##
        od;                                                                       ##
        if IsBound(clr) then                                                      ##
            if clr[pos]=1 then                                                    ##
                colour[3][Length(bnd[3])-5]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])-4]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])-1]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])]:=[-1]; # blue                           ##
            elif clr[pos]=2 then                                                  ##
                colour[3][Length(bnd[3])-5]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])-4]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])-1]:=[1]; # red                           ##
                colour[3][Length(bnd[3])]:=[1]; # red                             ##
            elif clr[pos]=3 then                                                  ##
                colour[3][Length(bnd[3])-5]:=[1]; # red                           ##
                colour[3][Length(bnd[3])-4]:=[1]; # red                           ##
                colour[3][Length(bnd[3])-1]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])]:=[-1]; # blue                           ##
            else                                                                  ##
                colour[3][Length(bnd[3])-5]:=[1]; # red                           ##
                colour[3][Length(bnd[3])-4]:=[1]; # red                           ##
                colour[3][Length(bnd[3])-1]:=[1]; # red                           ##
                colour[3][Length(bnd[3])]:=[1]; # red                             ##
            fi;                                                                   ##
        fi;                                                                       ##
        # 14-02-21 there are four more 2-cells                                    ##
        dif:=Difference(bnd[3][del_2_cell_1_i]{[2..5]},del_2_cell_1)[1];          ##
        Add( # l+8 bottom                                                         ##
            bnd[3],                                                               ##
            [                                                                     ##
                4,                                                                ##
                dif,                                                              ##
                m+5,                                                              ##
                m+18,                                                             ##
                m+19                                                              ##
            ]                                                                     ##
        );                                                                        ##
        Add(sub[3],Length(bnd[3]));                                               ##
        dif:=Difference(bnd[3][del_2_cell_2_i]{[2..5]},del_2_cell_2)[1];          ##
        Add( # l+9 bottom                                                         ##
            bnd[3],                                                               ##
            [                                                                     ##
                4,                                                                ##
                dif,                                                              ##
                m+15,                                                             ##
                m+20,                                                             ##
                m+21                                                              ##
            ]                                                                     ##
        );                                                                        ##
        Add(sub[3],Length(bnd[3]));                                               ##
        if IsBound(clr) then                                                      ##
            if clr[pos]=1 then                                                    ##
                colour[3][Length(bnd[3])-1]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])]:=[-1]; # blue                           ##
            elif clr[pos]=2 then                                                  ##
                colour[3][Length(bnd[3])-1]:=[-1]; # blue                         ##
                colour[3][Length(bnd[3])]:=[1]; # red                             ##
            elif clr[pos]=3 then                                                  ##
                colour[3][Length(bnd[3])-1]:=[1]; # red                           ##
                colour[3][Length(bnd[3])]:=[-1]; # blue                           ##
            else                                                                  ##
                colour[3][Length(bnd[3])-1]:=[1]; # red                           ##
                colour[3][Length(bnd[3])]:=[1]; # red                             ##
            fi;                                                                   ##
        fi;                                                                       ##
        ############################################################################
        Add(bnd[3],[2,m+4,m+5]); # l+10
        Add(sub[3],Length(bnd[3]));
        Add(bnd[3],[2,m+14,m+15]); # l+11
        Add(sub[3],Length(bnd[3]));
        # from this point onwards, cells added are only present in bnd, not sub
        Add( # l+12
            bnd[3],
            [
                6,
                Position(bnd[2],[2,a,c]),
                m,
                m+2,
                m+8,
                m+12,
                m+16
            ]
        );
        Add( # l+13
            bnd[3],
            [
                6,
                Position(bnd[2],[2,b,d]),
                m+1,
                m+3,
                m+9,
                m+13,
                m+17
            ]
        );
        # 14-02-21 another few 2-cells to bnd only
        Add( # l+14
            bnd[3],
            [
                3,
                Position(bnd[2],[2,a_abv,a]),
                m,
                m+18
            ]
        );
        Add( # l+15
            bnd[3],
            [
                3,
                Position(bnd[2],[2,b_abv,b]),
                m+1,
                m+19
            ]
        );
        Add( # l+16
            bnd[3],
            [
                3,
                Position(bnd[2],[2,c,c_blw]),
                m+2,
                m+20
            ]
        );
        Add( # l+17
            bnd[3],
            [
                3,
                Position(bnd[2],[2,d,d_blw]),
                m+3,
                m+21
            ]
        );
        # 19/04/21 found some 2-cells missing from ucap causing the chain complex
        # boundary matrices to detect that this wasn't a regular cw-complex
        for i in [0..3] do Add(ucap,l+14+i); od; 
        # 3-skeleton of intersection
        Add(
            bnd[4],
            [
                8,
                l+2,
                l+3,
                l+4,
                l+5,
                l+6,
                l+7,
                l+10,
                l+11
            ]
        );
        Add(
            bnd[4],
            [
                8,
                Position(
                    List(bnd[3],Set),
                    Set(
                        [
                            4,
                            Position(bnd[2],[2,a,b]),
                            Position(bnd[2],[2,c,d]),
                            Position(bnd[2],[2,a,c]),
                            Position(bnd[2],[2,b,d])
                        ]
                    )
                ),
                l,
                l+1,
                l+3,
                l+5,
                l+7,
                l+12,
                l+13
            ]
        );
        # 14-02-21 need two more 3-cells to finish off the structure
        Add(
            bnd[4],
            [
                5,
                del_2_cell_1_i,
                l,
                l+8,
                l+14,
                l+15
            ]
        );
        Add(
            bnd[4],
            [
                5,
                del_2_cell_2_i,
                l+1,
                l+9,
                l+16,
                l+17
            ]
        );

        if IsBound(clr) then
            f_clr:=clr[pos]*1;
            if f_clr=1 then
                f_clr1:=-1;
                f_clr2:=-1;
            elif f_clr=2 then
                f_clr1:=-1;
                f_clr2:=1;
            elif f_clr=3 then
                f_clr1:=1;
                f_clr2:=-1;
            else
                f_clr1:=1;
                f_clr2:=1;
            fi;
        else
            f_clr1:=0;
            f_clr2:=0;
        fi;

        for i in [1..Length(pipes)] do
            if IsBound(pipes[i][1][1]) then
                if pipes[i][1][1]=del_2_cell_1_i then
                    pipes[i]:=[
                        [
                            Position(
                                bnd[3],
                                Concatenation(
                                    [2],
                                    Set(
                                        Positions(
                                            bnd[2],
                                            Concatenation(
                                                [2],
                                                Set(
                                                    [a_abv,b_abv]
                                                )
                                            )
                                        )
                                    )
                                )
                            ),
                            l+8,
                            l+10
                        ],
                        [
                            Filtered(
                                pipes[i][2],
                                x->Length(
                                    Positions(
                                        bnd[2],
                                        bnd[2][x]
                                    )
                                )=2
                            )[1],
                            m+4,
                            m+18,
                            m+19
                        ],
                        [
                            a_abv,
                            b_abv,
                            n,
                            n+1
                        ],
                        '*', # just to mark this pipe as the vertical part of an intersection
                        f_clr1
                    ];
                elif pipes[i][1][1]=del_2_cell_2_i then
                    pipes[i]:=[
                        [
                            Position(
                                bnd[3],
                                Concatenation(
                                    [2],
                                    Set(
                                        Positions(
                                            bnd[2],
                                            Concatenation(
                                                [2],
                                                Set(
                                                    [c_blw,d_blw]
                                                )
                                            )
                                        )
                                    )
                                )
                            ),
                            l+9,
                            l+11
                        ],
                        [
                            Filtered(
                                pipes[i][2],
                                x->Length(
                                    Positions(
                                        bnd[2],
                                        bnd[2][x]
                                    )
                                )=2
                            )[1],
                            m+14,
                            m+20,
                            m+21
                        ],
                        [
                            c_blw,
                            d_blw,
                            n+6,
                            n+7
                        ],
                        '*',
                        f_clr2
                    ];
                fi;
            fi;
        od;
        for i in [1..Length(pipes)] do
            if not IsBound(pipes[i][4]) then
                int:=Set(
                    Intersection(
                        pipes[i][3],
                        [a,b,c,d]
                    )
                );
                if Length(int)=4 then
                    if crs[Position(List(crossings+l0,Set),int)]=0 then
                        pipes[i]:=[
                            [
                                l+2,
                                l+4,
                                l+6,
                                l+12,
                                l+13
                            ],
                            [
                                m,
                                m+1,
                                m+2,
                                m+3,
                                m+4,
                                m+14
                            ],
                            [
                                a,
                                b,
                                c,
                                d,
                                n,
                                n+1,
                                n+6,
                                n+7
                            ],
                        ];
                    fi;
                fi;
            fi;
        od;
    end; ##########################################################################

    for i in [1..Length(crs)] do
        if crs[i]=0 then
            IntersectingCylinders(
                crossings[i][1]+l0,
                crossings[i][3]+l0,
                crossings[i][2]+l0,
                crossings[i][4]+l0
            );
        fi;
    od;
    for i in [1..Length(pipes)] do # ignoring all self-intersections, join all
    # pipes which intersect at a common crossing point
        for j in [i+1..Length(pipes)] do
            if not IsBound(pipes[i][4]) and not IsBound(pipes[j][4]) then
                int:=Intersection(pipes[i][3],pipes[j][3]);
                if Length(int)>=2 and Intersection(int,Concatenation(crossings)+l0)=int then
                    Append(pipes[i][1],pipes[j][1]);
                    pipes[i][1]:=Set(pipes[i][1]); pipes[j][1]:=[];
                    Append(pipes[i][2],pipes[j][2]);
                    pipes[i][2]:=Set(pipes[i][2]); pipes[j][2]:=[];
                    Append(pipes[i][3],pipes[j][3]);
                    pipes[i][3]:=Set(pipes[i][3]); pipes[j][3]:=[];
                fi;
            fi;
        od;
    od;
    pipes:=Filtered(pipes,x->x<>[[],[],[]]);
    for i in [1..Length(pipes)] do # add the degree two 2-cells to the ends of each pipe
    # except for the vertical pipes at each intersection (they already have them)
        if not IsBound(pipes[i][4]) then
            for j in Filtered(bnd[3]{[l2+1..2*l2]},x->x[1]=2) do
                int:=Intersection(pipes[i][2],j{[2,3]});
                if int<>[] then
                    Add(pipes[i][1],Position(bnd[3],j));
                fi;
            od;
        fi;
    od;
    for i in [1..Length(pipes)] do # if two pipes' 0-skeletons intersect (at somewhere other
    # than a crossing), add a new 1-cell whose boundary is their intersection
        for j in [i+1..Length(pipes)] do
            int:=Intersection(pipes[i][3],pipes[j][3]);
            if int<>[] then
                if Intersection(int,Concatenation(crossings)+l0)=[] and
                    int<[2*l0,2*l0] then
                        Add(int,2,1);
                        Add(bnd[2],int);
                        Add(sub[2],Length(bnd[2]));
                fi;
            fi;
        od;
    od;
    for i in [1..Length(pipes)] do # replace the 1-cells (added pre IntersectingCylinders)
    # that occur twice with either the new cell we just added or the other occurance of
        for j in [1..Length(pipes[i][2])] do # that cell which is not in pipes[i][2]
            pos:=Positions(
                bnd[2],
                [
                    2,
                    bnd[2][pipes[i][2][j]][2],
                    bnd[2][pipes[i][2][j]][3]
                ]
            );
            if Length(pos)>=2 and
                [
                    bnd[2][pipes[i][2][j]][2],
                    bnd[2][pipes[i][2][j]][3]
                ]<[2*l0,2*l0] then
                    cell:=Last(Filtered(pos,x->x<>pipes[i][2][j]));
                    pipes[i][2][j]:=cell;
            fi;
        od;
    od;
    HorizontalOrVertical:=function(l)
        if l in hbars+l0 then
            return "horizontal";
        fi;
        return "vertical";
    end;
    for i in [1..Length(pipes)] do # if a pipe contains 1-cells from crossings then
        # filter out the horizontal/vertical 1-cells
        # should said pipe be vertical/horizontal itself
        # however if a pipe contains a self intersection remove both the horizontal
        # AND vertical 1-cells
        int:=Intersection(pipes[i][3],Concatenation(crossings)+l0);
        if Length(int)=4 then
            pos:=Position(List(crossings,Set)+l0,Set(int));
            if crs[pos]=0 then
                pipes[i][2]:=Difference(
                    pipes[i][2],
                    [
                        Position(bnd[2],[2,crossings[pos][1]+l0,crossings[pos][2]+l0]),
                        Position(bnd[2],[2,crossings[pos][3]+l0,crossings[pos][4]+l0]),
                        Position(bnd[2],[2,crossings[pos][1]+l0,crossings[pos][3]+l0]),
                        Position(bnd[2],[2,crossings[pos][2]+l0,crossings[pos][4]+l0])
                    ]
                );
            else
                ori:=HorizontalOrVertical(pipes[i][3]);
                if ori="horizontal" then
                    pipes[i][2]:=Difference(
                    pipes[i][2],
                    [
                        Position(bnd[2],[2,crossings[pos][1]+l0,crossings[pos][2]+l0]),
                        Position(bnd[2],[2,crossings[pos][3]+l0,crossings[pos][4]+l0]),
                    ]
                );
                else
                    pipes[i][2]:=Difference(
                    pipes[i][2],
                    [
                        Position(bnd[2],[2,crossings[pos][1]+l0,crossings[pos][3]+l0]),
                        Position(bnd[2],[2,crossings[pos][2]+l0,crossings[pos][4]+l0])
                    ]
                );
                fi;
            fi;
        fi;
    od;
    for i in [1..Length(pipes)] do # add new 2-cells to bnd[3] whose boundary is
    # the pipes[i][2] 1-cells
        cell:=Set(pipes[i][2]);
        Add(cell,Length(cell),1);
        Add(bnd[3],cell);
        Add(sub[3],Length(bnd[3]));
        Add(ucap,Length(bnd[3]));
        Add(pipes[i][1],Length(bnd[3]));
        if IsBound(pipes[i][5]) then
            if IsBound(clr) then
                colour[3][Length(bnd[3])]:=[pipes[i][5]];
            fi;
        fi;
    od;
    for i in [1..Length(pipes)] do # find where pipes intersect in their 2-skeleta
    # use this to form the 3-cells comprising the interiors of the pipes
        for j in [i+1..Length(pipes)] do
            if Intersection(pipes[i][1],pipes[j][1])<>[] then
                Append(pipes[i][1],pipes[j][1]);
                pipes[j][1]:=[];
            fi;
        od;
    od;
    for i in [1..Length(pipes)] do
        if pipes[i][1]<>[] then
            cell:=Set(pipes[i][1]);
            Add(cell,Length(cell),1);
            Add(bnd[4],cell);
        fi;
    od;
###################################################################################

# add a cap to both ends of D x [0,1] 
####################################################################################
    Add(
        bnd[3],
        [
            2,
            Positions(bnd[2],[2,1,l0])[1],
            Positions(bnd[2],[2,1,l0])[2]
        ]
    );
    Add(lcap,Length(bnd[3]));
    lcap:=Set(lcap);
    Add(lcap,Length(lcap),1);
    Add(bnd[4],lcap);
    Add(
        bnd[3],
        [
            2,
            Positions(bnd[2],[2,l0+1,2*l0])[1],
            Positions(bnd[2],[2,l0+1,2*l0])[2]
        ]
    );
    Add(ucap,Length(bnd[3]));
    ucap:=Set(ucap);
    Add(ucap,Length(ucap),1);
    Add(bnd[4],ucap);
####################################################################################

# add colour
####################################################################################
    for i in [1..Length(bnd[3])] do
        if not IsBound(colour[3][i]) then
            colour[3][i]:=[0];
        fi;
    od;
    colour_:={n,k}->colour[n+1][k];
####################################################################################
    x:=CWSubcomplexToRegularCWMap([RegularCWComplex(bnd),sub]);
    x!.colour:=colour_;
    return x;
end);

