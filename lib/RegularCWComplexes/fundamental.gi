#(C) Graham Ellis

##########################################################
##########################################################
InstallGlobalFunction(FundamentalGroupOfRegularCWComplex,
function(arg)
local P,Y,base,e,bool, b, vertices,edges,F, G, FhomG, r,x,w, gens, rels, 
      cells, 0cells,1cells, 2cells, 2boundaries, deform, EdgeToWord,
      EdgeToLoop, VertexToPath, loops, BOOL, homotopyOrientation,R,S,SS,i,ii,jj,      bnd,sn, y, j,
      indx, A, B, fn,LL;

Y:=arg[1];
base:=1;
BOOL:=true;

if Length(arg)>1 then 
if IsInt(arg[2]) then base:=arg[2]; fi;
if IsString(arg[2]) then BOOL:=false;  fi;
fi;

if Length(arg)>2 then 
if IsInt(arg[3]) then base:=arg[3];  fi;
if IsString(arg[3]) then BOOL:=false;  fi;
fi;



if Dimension(Y)<4 then    
     cells:=CriticalCellsOfRegularCWComplex(Y);
#    cells:=CocriticalCellsOfRegularCWComplex(Y,4);
else
   if IsBound(Y!.allcocriticalcells) then 
      cells:=CocriticalCellsOfRegularCWComplex(Y,Y!.allcocriticalcells);
   else
      cells:=CocriticalCellsOfRegularCWComplex(Y,3);
   fi;
fi;
Y!.criticalCells:=cells;

###########################################################
P:=SSortedList(List(Y!.orientation[2],Sum)); # changed 1 to 2 (Nov 2019)
if P=[0] then Y!.homotopyOrientation:=Y!.orientation{[1,2,3]};
else 
P:=TruncatedRegularCWComplex(Y,2);;
P!.orientation:=fail;
OrientRegularCWComplex(P);
Y!.homotopyOrientation:=P!.orientation;
fi;
Unbind(P);
###########################################################
0cells:=Filtered(cells,x->x[1]=0);
Apply(0cells,x->x[2]);
1cells:=Filtered(cells,x->x[1]=1);
Apply(1cells,x->x[2]);
2cells:=Filtered(cells,x->x[1]=2);
Apply(2cells,x->x[2]);
2boundaries:=1*List(2cells,x->[Y!.boundaries[3][x],Y!.homotopyOrientation[3][x]]);
Apply(2boundaries,x->[x[1]{[2..Length(x[1])]},x[2]]);


##NEED TO ORDER EACH BOUNDARY  #CHANGED 27/11/2018
for i in [1..Length(2boundaries)] do
R:=1*2boundaries[i];
S:=[1]; #SS:=[1];
indx:=[2..Length(R[1])];
for ii in [1..Length(R[1])-1] do
A:=Y!.boundaries[2][R[1][S[ii]]];
A:=A{[2,3]};
for jj in indx do
B:=Y!.boundaries[2][R[1][jj]];
B:=B{[2,3]};
if Length(Intersection(A,B))>0 then
Add(S,jj);
RemoveSet(indx,jj);
break;
fi;
od;
od;

2boundaries[i]:=[2boundaries[i][1]{S},2boundaries[i][2]{S}];
od;
##BOUNDARIES ORDERED

for j in [1..Length(2boundaries)] do
sn:=[];
bnd:=2boundaries[j][1];
###### This is a really thoughtless way to get the signs right!!!
###### And it is also wasteful of time.
A:=[Y!.boundaries[2][bnd[1]]{[2,3]}];
for i in [2..Length(bnd)]  do
x:=Y!.boundaries[2][bnd[i]]{[2,3]}; y:=A[Length(A)];;
if x[1] in y then sn[i]:=Y!.homotopyOrientation[2][bnd[i]][1]; Add(A,x{[1,2]}); 
else sn[i]:=-Y!.homotopyOrientation[2][bnd[i]][1]; Add(A,x{[2,1]}); fi;
od;
if A[1][1] in A[2] then A[1]:=A[1]{[2,1]}; sn[1]:=-Y!.homotopyOrientation[2][bnd[1]][1];
else sn[1]:=Y!.homotopyOrientation[2][bnd[1]][1]; fi;
2boundaries[j][2]:=sn;
######
######
########################################
od;




Apply(2boundaries,x->List([1..Length(x[1])],i->x[1][i]*x[2][i]));



deform:=ChainComplex(Y)!.homotopicalDeform;

Apply(2boundaries,x->Flat(List(x,a->deform(1,a))));

vertices:=[deform(0,base)];
edges:=[];
###################################
###################################
if not Length(0cells)=1 then 

bool:=true;
while bool do
bool:=false;

for e in 1cells do
b:=Y!.boundaries[2][e];
b:=b{[2,3]};
Apply(b,x->deform(0,x));

if b[1] in vertices and not b[2] in vertices
then Add(edges,e); Add(vertices,b[2]); bool:=true;
fi;
if b[2] in vertices and not b[1] in vertices
then Add(edges,e); Add(vertices,b[1]); bool:=true;
fi;
od;

od;

1cells:=Difference(1cells,edges);

1cells:=Filtered(1cells,e->deform(0,Y!.boundaries[2][e][2]) in vertices);
2cells:=Filtered(2cells,e->deform(1,Y!.boundaries[3][e][2]) in 1cells);
fi;
###################################
###################################


F:=FreeGroup(Length(1cells));
gens:=GeneratorsOfGroup(F);
if Length(gens)=0 then return F; fi;
rels:=[];
for r in 2boundaries do

w:=Identity(F);
for x in r do
if (not AbsInt(x) in edges) and deform(0,Y!.boundaries[2][AbsInt(x)][2]) in vertices then
w:=w*gens[Position(1cells,AbsInt(x))]^(SignInt(x));
fi;
od;

Add(rels,w);
od;


if BOOL then 
P:=PresentationFpGroup(F/rels);
TzOptions(P).printLevel:=0;
SimplifyPresentation(P);; 
G:=FpGroupPresentation(P);

else

G:=F/rels;
FhomG:=GroupHomomorphismByImagesNC(F,G,GeneratorsOfGroup(F),GeneratorsOfGroup(G));
fi;

##############################################
EdgeToWord:=function(e)
local r, x, w;

#r:=Flat(deform(1,e));   ###changed
r:=deform(1,e);


w:=Identity(F);
for x in r do
if (not AbsInt(x) in edges) and deform(0,Y!.boundaries[2][AbsInt(x)][2]) in vertices then
w:=w*gens[Position(1cells,AbsInt(x))]^(SignInt(x)) ; 
fi;
od;

return Image(FhomG,w);

end;
##############################################

if not BOOL then
G!.edgeToWord:=EdgeToWord;
fi;

loops:=StructuralCopy(1cells);

########################
VertexToPath:=function(v)
local path, e, pos;


path:=[];

while true do
if [v] in vertices then return path; 
else
e:=Y!.inverseVectorField[1][v];
w:=Y!.boundaries[2][e];
w:=w{[2,3]};
pos:=Position(w,v);
if pos=2 then v:=w[1]; Add(path,-e); else v:=w[2]; Add(path,e); fi;
fi;
od;

end;
########################

########################
EdgeToLoop:=function(e)
local loop, b;

b:=Y!.boundaries[2][e];
loop:=-Reversed(VertexToPath(b[2]));
Add(loop,e);
Append(loop,VertexToPath(b[3]));
return loop;
end;
########################

if Length(arg)>2 then
Apply(loops,EdgeToLoop);
G!.loops:=loops;
fi;

G!.deform:=deform;

return G;

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(FundamentalGroup,
"for regular CW-complexes",
[IsHapRegularCWComplex],
function(Y)
local F;
F:= FundamentalGroupOfRegularCWComplex(Y);
return F;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(FundamentalGroupWithPathReps,
"for regular CW-complexes",
[IsHapRegularCWComplex],
function(Y)
local F;
F:= FundamentalGroupOfRegularCWComplex(Y,1,"nosimplify");
return F;
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallMethod(FundamentalGroup,
"for regular CW-complex",
[IsHapRegularCWComplex,IsInt],
function(Y,n)
local bool,F;
F:= FundamentalGroupOfRegularCWComplex(Y,n);
return F;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(FundamentalGroup,
"for simplicial complexes",
[IsHapSimplicialComplex],
function(K)
local Y,c;
if Dimension(K)=2 then
return FundamentalGroupSimplicialTwoComplex(K);
fi;
Y:=SimplicialComplexToRegularCWComplex(K,3);;
c:=CocriticalCellsOfRegularCWComplex(Y,3);
return FundamentalGroup(Y);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(FundamentalGroup,
"for  pure cubical complexes",
[IsHapPureCubicalComplex],
function(M)
local Y,c;
Y:=CubicalComplexToRegularCWComplex(M,3);;
if Dimension(Y)<4 then 
c:=CriticalCellsOfRegularCWComplex(Y);
else
c:=CocriticalCellsOfRegularCWComplex(Y,3);
fi;
return FundamentalGroup(Y);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(FundamentalGroup,
"for  pure permutahedral complexes",
[IsHapPurePermutahedralComplex],
function(M)
local Y,c;
Y:=RegularCWComplex(M);;
if Dimension(Y)<4 then
c:=CriticalCellsOfRegularCWComplex(Y);
else
c:=CocriticalCellsOfRegularCWComplex(Y,3);
fi;
return FundamentalGroup(Y);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(FundamentalGroup,
"for  pure Regular CW-Maps",
[IsHapRegularCWMap],
function(map);
return FundamentalGroupOfRegularCWMap(map);
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallOtherMethod(FundamentalGroup,
"for  pure Regular CW-Maps with specified base-point",
[IsHapRegularCWMap,IsInt],
function(map,base);
return FundamentalGroupOfRegularCWMap(map,base);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(FundamentalGroup,
"for cubical complexes",
[IsHapCubicalComplex],
function(M)
local Y,c;
Y:=CubicalComplexToRegularCWComplex(M,3);;
if Dimension(Y)<4 then
c:=CriticalCellsOfRegularCWComplex(Y);
else
c:=CocriticalCellsOfRegularCWComplex(Y,3);
fi;
return FundamentalGroup(Y);
end);
##########################################################
##########################################################



#################################################
#################################################
InstallGlobalFunction(BoundaryPairOfPureRegularCWComplex,
function(Y)
local B, map, perm,invperm, x, pm, cnt;

B:=BoundaryOfPureRegularCWComplex(Y);
perm:=B!.perm;
invperm:=List([1..Length(perm)],i->[]);
for x in [1..Length(perm)] do
pm:=perm[x];
cnt:=0;
while cnt<Length(pm) do
cnt:=cnt+1;
if IsBound(pm[cnt]) then invperm[x][pm[cnt]]:=cnt; fi;
od;
od;

#########################
map:=function(n,i);
return 1*invperm[n+1][i];
end;
#########################

return Objectify(HapRegularCWMap,
       rec(
           source:=B,
           target:=Y,
           mapping:=map));
end);
#################################################
#################################################

#################################################
#################################################
InstallOtherMethod(Source,
"Source of a RegularCWMap",
[IsHapRegularCWMap],
function(map)
return map!.source;
end);
#################################################
#################################################

#################################################
#################################################
InstallOtherMethod(Target,
"Target of a RegularCWMap",
[IsHapRegularCWMap],
function(map)
return map!.target;
end);
#################################################
#################################################


#################################################
#################################################
InstallGlobalFunction(FundamentalGroupOfRegularCWMap,
function(arg)
local map, pntS, pntT,GS, GT, S, T, mapfn, loops,gensS, x, w;

map:=arg[1];
S:=Source(map);
T:=Target(map);
mapfn:=map!.mapping;

if Length(arg)>1 then pntS:=arg[2]; else pntS:=1; fi;
pntT:=mapfn(0,pntS);

GS:=FundamentalGroupOfRegularCWComplex(S,pntS,"nosimplify");
GT:=FundamentalGroupOfRegularCWComplex(T,pntT,"nosimplify");

gensS:=GeneratorsOfGroup(GS);

if Length(gensS)=0 then return
GroupHomomorphismByImagesNC(Group(Identity(GT)),GT,[Identity(GT)],[Identity(GT)]); fi;

loops:=[];
for x in GS!.loops do
w:= List(x,i->SignInt(i)*mapfn(1,AbsInt(i))) ;

#Apply(w,i->GT!.edgeToWord(AbsInt(i))^SignInt(i));
#change:
Apply(w,GT!.edgeToWord);

Add(loops, Product(w));  
od;

return GroupHomomorphismByImagesNC(GS,GT,gensS,loops);;
end);
#################################################
#################################################


############################################
InstallGlobalFunction(CriticalBoundaryCells,
function(Y,n,k)
local C, Deform, B, BB, j;

if n<=0 then return []; fi;

C:=ChainComplexOfRegularCWComplexWithVectorField(Y);
Deform:=C!.deform;

B:=Y!.boundaries[n+1,k];
B:=B{[2..B[1]+1]};
BB:=[];

for j in [1..Length(B)] do
Append(BB,1*Deform(n-1,Y!.orientation[n+1][k][j]*B[j]));
od;


return BB;

end);
############################################

