InstallGlobalFunction(EntireBoundary, function(bnd,d,n)
local V, U, D, tmp, v , cells;

if d=0 then return []; fi;

D:=d+1;
V:=StructuralCopy(bnd[D][n]);
V:=V{[2..Length(V)]};
V:=SSortedList(V);
cells:=List(V,i->[D-2,i]);
D:=D-1;

while D>1 do
tmp:=[];
for v in V do
U:=StructuralCopy(bnd[D][v]);
U:=U{[2..Length(U)]};
Append(tmp,U);
Append(cells,List(U,i->[D-2,i]));   #Changed D-1 to D-2  (May 2022)
od;
tmp:=SSortedList(tmp);
cells:=SSortedList(cells);
V:=tmp;
D:=D-1;
od;

return cells;
end);

InstallGlobalFunction(Cobound, function(bnd)
local dim, cobnd, b, k, i, j, n;

dim:=PositionProperty(bnd,IsEmpty)-2;
cobnd:=[];; #Cobnd[n+1] contains the info on n-cells.
for n in [0..dim] do
  #k:=2*(n+1)+1;#k:=1+2^(n+1);
  cobnd[n+1]:=List(bnd[n+1],i->[0]);
  for j in [1..Length(bnd[n+2])] do
    b:=bnd[n+2][j];
    k:=Length(b);
    for i in b{[2..k]} do
      cobnd[n+1][i][1]:=cobnd[n+1][i][1]+1;
      Add(cobnd[n+1][i],j);
    od;
  od;
od;
cobnd[dim+1]:=List(bnd[dim+1],a->[0]);
return cobnd;
end);

InstallGlobalFunction(ShouldRemove, function(FL,d1,n)
local d2, cob, V1, V2, V3, bnd, cobnd, Y;

bnd := FL;

#Y := RegularCWComplex(FL);
cobnd := Cobound(bnd);

if n > Length(bnd[d1]) then
    return false;
fi;

if not (cobnd[d1][n][1] = 2 and bnd[d1][n][1] > 0) then
    return false;
fi;

d2 := d1+1;
cob := 1*cobnd[d1][n];

if Length(cob) < 3 then
    return false;
fi;

if d2 > Length(cobnd) then
    return false;
fi;

if cob[2] > Length(cobnd[d2]) or cob[3] > Length(cobnd[d2]) then
    return false;
fi;

if not SortedList(cobnd[d2][cob[2]]) = SortedList(cobnd[d2][cob[3]]) then
    return false;
fi;


# Print("d1=", d1, " n=", n, "\n");
# Print("d2=", d2, "\n");
# Print("cob=", cob, "\n");
# Print("Length(bnd)=", Length(bnd), "\n");
# Print("Length(bnd[d1])=", Length(bnd[d1]), "\n");
# Print("Length(bnd[d2])=", Length(bnd[d2]), "\n");

V1:=EntireBoundary(bnd,d1-1,n);
V2:=EntireBoundary(bnd,d1,cobnd[d1][n][2]);
V3:=EntireBoundary(bnd,d1,cobnd[d1][n][3]);

if not 1+Size(V1) = Size(Intersection(V2,V3)) then
    return false;
fi;

return true;
end);


InstallGlobalFunction(JoinCells, function(FL,d1,n)
local bnd, cobnd, cob, d, d2, d3, m, t, poss, s, pos, u, perm, cnt, F, x, i;

d2:=d1+1;
d3:=d1+2;
bnd:=FL;
cobnd:=Cobound(bnd);
#cobnd:=Y!.coboundaries;
cob:=1*cobnd[d1][n];
##
##REDUCE COBOUNDARIES OF BOUNDARIES OF nTH CELL
if d1>1 then
d:=d1-1;;
for m in bnd[d1][n]{[2..Length(bnd[d1][n])]} do
    t:=1*cobnd[d][m]{[2..Length(cobnd[d][m])]};;
    poss:=Position(t,n);;
    Remove(t,poss);;
    cobnd[d][m]:=Concatenation([Length(t)],t);;
od;;
fi;;
#Print("cbnd ",cobnd,"\n");
##
##COBOUNDARIES OF BOUNDARIES REDUCED

##
##REMOVE nTH CELL, ITS COBOUNDARY, AND ADJUST ITS PRESENCE IN BOUNDARIES
##OF ITS COBOUNDARIES
bnd[d1][n]:=[0];;
cobnd[d1][n]:=[0];;
#Print("cobnd",cobnd,"\n");
#Print("bnd",bnd,"\n");

s:=bnd[d2][cob[2]];;
s:=s{[2..Length(s)]};;
pos:=Position(s,n);;
Remove(s,pos);;


t:=bnd[d2][cob[3]];;
t:=t{[2..Length(t)]};;
pos:=Position(t,n);;
Remove(t,pos);;


u:=SortedList(Concatenation(s,t));;

if IsEmpty(u) then
    bnd[d2][cob[2]]:=[1,0];;
else
    bnd[d2][cob[2]]:=Concatenation([Length(u)],u);;
fi;;

##
##nTH CELL AND ITS PRESENCE REMOVED


##
##FOR SECOND CELL OF DIMENSION n+1  REDUCE THE BOUNDARIES OF ITS 
##COBOUNDARIES

for m in cobnd[d2][cob[3]]{[2..Length(cobnd[d2][cob[3]])]} do
t:=bnd[d3][m]{[2..Length(bnd[d3][m])]};;

pos:=Position(t,cob[3]);;
Remove(t,pos);;
bnd[d3][m]:=Concatenation([Length(t)],t);;
od;;


##
##SECOND CELL  BOUNDARIES OF ITS COBOUNDARIES REDUCED

##
##REMOVE PRESENCE OF SECOND CELL IN COBOUNDARIES OF ITS BOUNDARIES
for m in bnd[d2][cob[3]]{[2..Length(bnd[d2][cob[3]])]} do
if cobnd[d1][m][1]>0 then

t:=cobnd[d1][m]{[2..Length(cobnd[d1][m])]};;

pos:=Position(t,cob[3]);;
t[pos]:=cob[2];;
t:=SSortedList(t);;
cobnd[d1][m]:=Concatenation([Length(t)],t);;

fi;;
od;;


bnd[d2][cob[3]]:=[0];;

cobnd[d2][cob[3]]:=[0];;

perm:=[];  ##perm[d][i] will be the new position of
           ##the old i-th cell of dimension d;

for d in [1..Length(bnd)] do
perm[d]:=[];
cnt:=0;
for n in [1..Length(bnd[d])] do
if bnd[d][n][1]=0 then cnt:=cnt+1;
else
perm[d][n]:=n-cnt;
fi;
od;
od;



for d in [1..Length(bnd)] do
F:=1*Filtered([1..Length(bnd[d])],i->not bnd[d][i][1]=0);
bnd[d]:=1*((bnd[d]){F});

if d>1 then
d1:=d-1;
for x in bnd[d] do
    if x[2]>0 then
        for i in [2..Length(x)] do
            x[i]:=1*perm[d1][x[i]];
        od;
    fi;
od;
fi;
od;

#if all lower-dimensional cells have been deleted, move all the remaining cells downwards
while Length(bnd) > 0 and IsEmpty(bnd[1]) do
    Remove(bnd, 1);
od;

return bnd;
end);
#################################

