InstallGlobalFunction(SimplifyLattice, function(L,kstar,n) #n=dimension of the center of the kstar
local Changed, b, d, j, m, h, g, g1, g2, dim, face, cell, signs, faceorbit, Stab, topcells, involvedcell1, involvedcell2,
ind1, ind2, sign1, sign2, incidence1, incidence2, mergecell1, mergecell2, FLgroups, FLcodim1, permFL, FL, cob;

FL:=ForgetG(L,kstar,n);
permFL:=StructuralCopy(FL);
FL:=StructuralCopy(FL);
FLcodim1:=[1..Length(FL[Length(FL)-2])];
FLgroups:=List([1..Length(FL[Length(FL)-1])], i -> [i]);
dim:=Length(FL)-1;

repeat
    Changed:=false;
    #topcells:=
    for d in [1 .. Length(FL)-2] do
        m:=1;
        # JoinCells can remove an entire dimension during this pass.
        while IsBound(FL[d]) and m <= Length(FL[d]) do
            if FL[d][m][1] > 0 and ShouldRemove(FL,d,m) then
                if d=Length(FL)-2 then #for top dimensional cells, we have to keep track of signs
                    mergecell1:=Cobound(FL)[d][m][2];
                    mergecell2:=Cobound(FL)[d][m][3]; #these are the two cells which will be merged
                    #adjust signs
                    b:=FLcodim1[m]; #old index of the cell that will be deleted
                    ind1 := fail;
                    involvedcell1 := fail;
                    ind2 := fail;
                    involvedcell2 := fail; 
                    for j in FLgroups[mergecell1] do
                        if b in permFL[dim][AbsoluteValue(j)]{[2..Length(permFL[dim][AbsoluteValue(j)])]} then
                            #involvedcell1:=permFL[dim][j];
                            ind1:=AbsoluteValue(j);
                            sign1:=SignInt(j);
                            break;
                        fi;
                    od;
                    for j in FLgroups[mergecell2] do
                        if b in permFL[dim][AbsoluteValue(j)]{[2..Length(permFL[dim][AbsoluteValue(j)])]} then
                            #involvedcell2:=permFL[dim][j];
                            ind2:=AbsoluteValue(j);
                            sign2:=SignInt(j);
                            break;
                        fi;
                    od;
                    #Print("Merging the cells ", FLgroups[mergecell1], " and ", FLgroups[mergecell2]," along the face ",b,"\n", " the components involved are ",ind1," and ",ind2,"\n");
                    if ind1 = fail then
                        Error("Could not find face ", b, " in group ", FLgroups[mergecell1]);
                    fi;
                    face:=kstar[dim-1][b];
                    faceorbit:=face[1];
                    g:=face[2];
                    involvedcell1:=kstar[dim][ind1][1];
                    g1:=kstar[dim][ind1][2];
                    involvedcell2:=kstar[dim][ind2][1];
                    g2:=kstar[dim][ind2][2];
                    for cell in Filtered(L!.boundary(n+dim-1,involvedcell1),p->AbsoluteValue(p[1])=faceorbit) do
                        h:=L!.elts[cell[2]];
                        Stab:=kstar[dim-1][b][3];
                        if Inverse(h)*Inverse(g1)*g in Stab then
                            incidence1:=SignInt(cell[1]);
                            break;
                        fi;
                    od;
                    for cell in Filtered(L!.boundary(n+dim-1,involvedcell2),p->AbsoluteValue(p[1])=faceorbit) do
                        h:=L!.elts[cell[2]];
                        Stab:=kstar[dim-1][b][3];
                        if Inverse(h)*Inverse(g2)*g in Stab then
                            incidence2:=SignInt(cell[1]);
                            break;
                        fi;
                    od;
                    if sign1*incidence1*sign2*incidence2 = 1 then
                        FLgroups[mergecell2]:=List(FLgroups[mergecell2], y -> -y);
                    fi;
                    Append(FLgroups[mergecell1],FLgroups[mergecell2]);
                    Remove(FLgroups,mergecell2);
                    Remove(FLcodim1,m);
                fi;

                # NEW: if this JoinCells call will delete a cell from the
                # codim-1 level (FL[Length(FL)-2]) via its d2-side removal,
                # keep FLcodim1 in sync.
                if d = Length(FL)-3 then
                    cob := Cobound(FL)[d][m];
                    FL:=JoinCells(FL,d,m);
                    Remove(FLcodim1, cob[3]);
                else
                    FL:=JoinCells(FL,d,m);
                fi;
                
                # FL:=JoinCells(FL,d,m);
                Changed:=true;               
                m:=1;
            else
                m:=m+1;
            fi;

        od;
    od;
until not Changed;

signs := List([1..Length(FLgroups[1])], x -> 0);
for j in FLgroups[1] do
    signs[AbsoluteValue(j)] := SignInt(j);
od;
return [FL,signs];
end);

InstallGlobalFunction(HomeomorphicToDisk, function(L,e,n,k)
local kstar, B, simplified;

B:=false;
kstar:=kStar(L,e,n,k);
# A zero cellular boundary can leave no incident cells at this level.
if kstar[k+1]=[] then return false; fi;
simplified:=SimplifyLattice(L,kstar,n)[1];
if Length(simplified)=2 and Length(simplified[1])=1 then
    B:=true;
fi;
return B;
end);
############################################

