

########################################################
########################################################
InstallGlobalFunction(RegularCWComplexWithRemovedCell,
function(Y,n,k)
local B, A;

if not n=Dimension(Y) then 
Print("Only top dimensional cells can be removed.\n");
return fail;
fi;

B:=1*Y!.boundaries;
A:=Filtered( [1..Length(B[n+1])],i->not i=k);
B[n+1]:=B[n+1]{A};
return RegularCWComplex(B);
end);
########################################################
########################################################

########################################################
########################################################
InstallGlobalFunction(RegularCWComplex_AttachCellDestructive,
function(arg)
local Y, n, b, B, CB, e, NrCells, F, k, L1, L2, s,
d, d1, d2, bnd, i, j, bool, bb, t, L, S, T, orien;
##
##WARNING: This function unbinds/assigns fail to all components except 
##         Y!.nrCells, Y!.boundaries, Y!.coboundaries, Y!.orientation.
##         It inputs a CW-complex Y, a dimension n, a list b of 
##         boundary n-1-cells of anew n-cell.
##         It returns the number e of the new n-cell. 

Y:=arg[1];
n:=arg[2];
if Length(arg)>2 then b:=1*arg[3]; fi;

######KEEP BASIC DATA of Y AND DELETE THE REST
######

###########################
NrCells:=function(n);
return Length(Y!.boundaries[n+1]);
end;
###########################

Y!.nrCells:=NrCells;
Y!.vectorField:=fail;
Y!.inverseVectorField:=fail;
Y!.criticalCells:=fail;
L1:=[ "coboundaries", "properties", "vectorField", "boundaries", 
  "nrCells", "orientation", "inverseVectorField", "criticalCells" ];
L2:=NamesOfComponents(Y);
for k in L2 do
if not k in L1 then s:=EvalString(Concatenation("Y!.",k)); 
Unbind(s); fi;
od;
if n>Dimension(Y) then
for k in Y!.properties do
if k[1]="dimension" then k[2]:=n; break; fi;
od;
Add(Y!.boundaries,[]);
Add(Y!.coboundaries,[]);
Add(Y!.orientation,[]);
fi;

########
########NON-BASIC DATA HAS NOW BEEN REMOVED

B:=Y!.boundaries;;
CB:=Y!.coboundaries;;
e:=Length(B[n+1])+1;

if n=0 then
Add(B[1],[1,0]);
else
Add(B[n+1], Concatenation( [Length(b)] , 1*SortedList(b)));
for k in b do
CB[n][k][1]:=CB[n][k][1]+1; Add(CB[n][k],e);
od;
fi;
CB[n+1][e]:=[0];

if n=0 then Y!.orientation[1][e]:=[1]; fi;
if n=1 then Y!.orientation[2][e]:=[1, -1]; fi;
if n>1 then Y!.orientation[n+1][e]:=0*[1..Length(Y!.boundaries[n+1][e])-1];

#####################
#####################MODIFIED FROM ORIENTREGULARCWCOMPLEX
d:=n+1;
bnd:=Y!.boundaries;
orien:=Y!.orientation;
  d1:=d-1;
  d2:=d-2;
  #for i in [1..Length(bnd[d])] do
  i:=e;
    b:=bnd[d][i]{[2..Length(bnd[d][i])]};
    bb:=[];
    for j in [1..Length(b)] do
      Add(bb, bnd[d1][b[j]]{[2..Length(bnd[d1][b[j]])]}   );
    od;
    orien[d][i][1]:=1;
    S:=[1..Length(b)];
    T:=[1..Length(b)];
    for s in [2..Length(b)] do
    Unbind(S[s]);
    od;
    Unbind(T[1]);

    while 0 in orien[d][i] do
    ###############################
    bool:=false;
    for s in S do
    for t in T do
      L:=Intersection(bb[s], bb[t]);
      if Length(L)>0 then
        S[t]:=t;
        Unbind(T[t]);
        bool:=true;
        if orien[d][i][s]*orien[d1][b[s]][Position(bb[s],L[1])]=
           orien[d1][b[t]][Position(bb[t],L[1])]
           then orien[d][i][t]:=-1;
           else
           orien[d][i][t]:=1;
        fi;
        break;
      fi;
    od;
    od;
    ###############################
    od;
  #od;


#####################END OF MODIFIED CODE
#####################

fi;
return e;
end);
########################################################
########################################################

########################################################
########################################################
RegularCWSphere:=function(n)
local S,k;

S:=[];
S[1]:=[ [1,0],[1,0] ];
for k in [1..n] do
Add(S, [ [2,1,2],[2,1,2] ]);
od;
Add(S,[]);
return RegularCWComplex(S);

end;
########################################################
########################################################

########################################################
########################################################
RegularCWClosedSurface:=function(n)
local S;

S:=ClosedSurface(n);
S:=RegularCWComplex(S);
S:=SimplifiedComplex(S);  #SHOULD RE-CODE THIS TO AVOID UNNECESSARY COMPUTATIONS
return S;

end;
########################################################
########################################################

########################################################
########################################################
RegularCWDiscreteSpace:=function(n)
local S, bnd;

bnd:=[[],[],[]];
bnd[1]:=List([1..n],i->[1,0]);
S:=RegularCWComplex(bnd);
return S;

end;
########################################################
########################################################

########################################################
########################################################
SphericalKnotComplementWithBoundary:=function(AP)
local K,f;

K:=SphericalKnotComplement(AP);;
f:=BoundaryPairOfPureRegularCWComplex(K);;

return f;
end;
########################################################
########################################################
RegularCWComplex_WedgeSum:=function(X,Y,uu,vv)
local  bndX, bndY, bnd, perm0,  perm, u, v, x, y, k;

u:=uu; v:=vv;
bndX:=1*X!.boundaries;
bndY:=1*Y!.boundaries;
for k in [1..Dimension(X)-Dimension(Y)] do
Add(bndY,[]);
od;
for k in [1..Dimension(Y)-Dimension(X)] do
Add(bndX,[]);
od;

######################
perm0:=function(i);
if i=v then return u; fi;
if i<v then return i+X!.nrCells(0); fi;
return i+X!.nrCells(0)-1;
end;
######################

######################
perm:=function(k,i);
return i+X!.nrCells(k);
end;
######################

bnd:=[];
bnd[1]:=List([1..X!.nrCells(0)+Y!.nrCells(0)-1],i->[1,0]);
bnd[2]:=bndX[2];
for x in bndY[2] do
y:=Concatenation( [x[1]], List(x{[2..x[1]+1]},i->perm0(i)));
Add(bnd[2],1*y);
od;
for k in [2..Maximum(Dimension(X),Dimension(Y))] do
bnd[k+1]:=bndX[k+1];
   for x in bndY[k+1] do
   y:=Concatenation( [x[1]], List(x{[2..x[1]+1]},i->perm(k-1,i)));
   Add(bnd[k+1],y);
od;

od;

Add(bnd,[]);
return RegularCWComplex(bnd);
end;
########################################################
########################################################

########################################################
########################################################
RegularCWComplex_DisjointUnion:=function(X,Y)
local  bndX, bndY, bnd, perm, u, v, x, y, k;

bndX:=1*X!.boundaries;
bndY:=1*Y!.boundaries;
for k in [1..Dimension(X)-Dimension(Y)] do
Add(bndY,[]);
od;
for k in [1..Dimension(Y)-Dimension(X)] do
Add(bndX,[]);
od;

######################
perm:=function(k,i);
return i+X!.nrCells(k);
end;
######################

bnd:=[];
bnd[1]:=List([1..X!.nrCells(0)+Y!.nrCells(0)],i->[1,0]);
#bnd[2]:=bndX[2];
#for x in bndY[2] do
#y:=Concatenation( [x[1]], List(x{[2..x[1]+1]},i->perm0(i)));
#Add(bnd[2],1*y);
#od;
for k in [1..Maximum(Dimension(X),Dimension(Y))] do
bnd[k+1]:=bndX[k+1];
   for x in bndY[k+1] do
   y:=Concatenation( [x[1]], List(x{[2..x[1]+1]},i->perm(k-1,i)));
   Add(bnd[k+1],y);
od;

od;

Add(bnd,[]);
return RegularCWComplex(bnd);
end;
########################################################
########################################################


########################################################
########################################################
RegularCWComplexWithAttachedRelatorCells :=function(arg)
local YY,G,m,P, gpath, split,  loops, vpairs, bool, i, j, x, vertices,
      e0,e1,v,u, L, pairs1, pairs2,nodes,extra1cells,pos,
      wedge, T, gg,  g,Y;

YY:=arg[1];
G:=arg[2];

Y:=RegularCWComplex(1*YY!.boundaries);
loops:=G!.loops;

for m in [3..Length(arg)] do

gg:=arg[m];
g:=ExtRepOfObj(gg);

gpath:=[];;
## gpath is the sequence of signed edges representing the group element g.

for i in [1..Length(g)/2] do
if g[2*i]>0 then
   for j in [1..g[2*i]] do Append(gpath,loops[g[2*i-1]]); od;
else
   for j in [1..-g[2*i]] do Append(gpath,-1*Reversed(loops[g[2*i-1]])); od;
fi;
od;

## The following removes any [ ...,k,-k, ...] occurences in gpath as
## these are not needed and their removal yields a minor efficiency gain.

bool:=true;
while bool do
bool:=false;
for i in [1..Length(gpath)-1] do
if  gpath[i]+gpath[i+1]=0 then gpath[i]:=0; gpath[i+1]:=0; bool:=true; fi;
od;
gpath:=Filtered(gpath,a-> not a=0);
od;

## vpairs is the list of boundary vertex pairs of the edges in gpath

vpairs:=[];
for x in gpath do
Add(vpairs,Y!.boundaries[2][AbsInt(x)]{[2,3]});
od;

## Next we'll order each pair in vpairs to achieve the form
## [...,[u,v],[v,w],[w,x],...]

if Length(gpath)>1 then
if vpairs[1][1] in vpairs[2] then vpairs[1]:=Reversed(vpairs[1]); fi;
for i in [2..Length(vpairs)] do
if not vpairs[i][1] = vpairs[i-1][2] then vpairs[i]:=Reversed(vpairs[i]); fi;
od;
fi;

## pairs1[i] contains the source vertex of edge i in gpath
## pairs2[i] contains the target vertex of edge i in gpath
pairs1:=List(vpairs,x->x[1]);
pairs2:=List(vpairs,x->x[2]);

## split is a decomposition of gpath into a sequence of simply connected paths.

split:=[];
vertices:=[pairs1[1],pairs2[1]]; #This list is used to make sure a
                                 #vertex is not hit twice.
P:=[gpath[1]];
for i in [2..Length(gpath)] do
if pairs2[i] in vertices then
    Add(split,List(P,AbsInt));
    P:=[gpath[i]];
    vertices:=[pairs1[i],pairs2[i]];
else
    Add(P,gpath[i]);
    Add(vertices,pairs1[i]);
    Add(vertices,pairs2[i]);
fi;
od;
Add(split,List(P,AbsInt));


if Length(split)=1 then
RegularCWComplex_AttachCellDestructive(Y,2,split[1]);
return Y;
fi;

## Attach one 0-cell
e0:=RegularCWComplex_AttachCellDestructive(Y,0);

## Attach one 1-cell for each initial vertex in the paths in split
L:=[];
pos:=1;
for i in [1..Length(split)] do
Add(L,pos); pos:=pos+Length(split[i]);
od;
nodes:=List(L,i->pairs1[i]);
extra1cells:=List(nodes,n->RegularCWComplex_AttachCellDestructive(Y,1,SortedList([e0,n])));

## Attach one "wedge" 2-cell for each term in the list split
wedge:=[];
for i in [1..Length(split)] do
    P:=1*split[i];
    u:=extra1cells[i];
    if i<Length(split) then
        v:=extra1cells[i+1];
    else
        v:=extra1cells[1];
    fi;
    Add(P,u);
    Add(P,v);
wedge[i]:=RegularCWComplex_AttachCellDestructive(Y,2,SortedList(P));
od;

od;
Y:=SimplifiedComplex(Y);
return Y;
end;
########################################################
########################################################

#########################################
#########################################
InstallMethod(Suspension,
"Suspension of regular CW complex",
[IsHapRegularCWComplex],
function(Y)
local B, n;;

B:=[];
B[1]:= [[1,0],[1,0]];
B[2]:=List([1..Y!.nrCells(0)], i->[2,1,2]);

for n in [1..Dimension(Y)] do
B[n+2]:=1*Y!.boundaries[n+1];
od;
Add(B,[]);

return RegularCWComplex(B);
end);
############################################
############################################

#########################################
#########################################
InstallOtherMethod(Suspension,
"n-fold suspension of regular CW complex",
[IsHapRegularCWComplex,IsInt],
function(Y,n) local S, i;

if n=0 then return Y; fi;
S:=Y;
for i in [1..n] do
S:=Suspension(S);
od;
return S;
end);
##########################################
##########################################

##########################################
#########################################
InstallOtherMethod(Suspension,
"n-fold suspension of regular CW complex",
[IsHapSimplicialComplex,IsInt],
function(Y,n) local S, i;

if n=0 then return Y; fi;
S:=Y;
for i in [1..n] do
S:=Suspension(S);
od;
return S;
end);
##########################################
##########################################


##########################################
#########################################
InstallOtherMethod(Suspension,
"Suspension of simplicial complex",
[IsHapSimplicialComplex],
function(K)
local M, S, top, bot, m;;

M:=MaximalSimplicesOfSimplicialComplex(IntegerSimplicialComplex(K));
bot:=Maximum(Flat(M))+1;
top:=bot+1;
S:=[];

for m in M do
Add(S,m);
Add(S,Concatenation(m,[bot]));
Add(S,Concatenation(m,[top]));
od;
return SimplicialComplex(S);
end);
############################################
############################################


############################################
############################################
InstallGlobalFunction(Suspension_alt,
function(Y)
local SY, B, top, bot, dim, dims,dimss, k,n,x,bnd;;

dim:=Dimension(Y);
dims:=List([0..dim],i->Y!.nrCells(i));
B:=1*Y!.boundaries;
Add(B,[]);
Add(B[1],[1,0]);
Add(B[1],[1,0]);
top:=Length(B[1]);
bot:=top-1;

for k in [1..dims[1]] do
bnd:=[2,k,top];
Add(B[2],bnd);
od;
for k in [1..dims[1]] do
bnd:=[2,k,bot];
Add(B[2],bnd);
od;

for n in [1..dim] do
for k in [1..dims[n+1]] do
x:=1*Y!.boundaries[n+1][k];
bnd:=[x[1]+1,k];
x:=x{[2..Length(x)]}+dims[n+1];
Append(bnd,SortedList(x));
Add(B[n+2],bnd);
od;
od;

dimss:=List([1..dim+1],i->Length(B[i]));
dimss[2]:=dimss[2]-dims[1];

for n in [1..dim] do
for k in [1..dims[n+1]] do
x:=1*Y!.boundaries[n+1][k];
bnd:=[x[1]+1,k];
x:=x{[2..Length(x)]}+dimss[n+1];
Append(bnd,SortedList(x));
Add(B[n+2],bnd);
od;
od;

SY:=RegularCWComplex(B);
return SY;
end);
#########################################
#########################################

