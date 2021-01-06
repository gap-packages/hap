
####################################################
InstallGlobalFunction(LiftedRegularCWMap,
function(f,p)
local Y, B, W, WB, WBcells, WBcellsinv, Bcells, WBbnd,  
WBorient, cnt, ff, n, i;

#       ff
#    WB ---> W                 ff is the lift of f.
#            |
#            | p
#       f    v
#    B ----> Y
Y:=Target(f);
B:=Source(f);
W:=Source(p);
if not Y=Target(p) then return fail; fi;

Bcells:=[];
for n in [1..Length(B!.boundaries)] do
Bcells[n]:=SSortedList(List([1..B!.nrCells(n-1)], i->f!.mapping(n-1,i)));
od;

WBcells:=[];
WBcellsinv:=[];
for n in [1..Length(W!.boundaries)] do
WBcells[n]:=[];
WBcellsinv[n]:=[];
cnt:=0;
for i in [1..Length(W!.boundaries[n])] do
if p!.mapping(n-1,i) in Bcells[n] then Add(WBcells[n], i); 
cnt:=cnt+1;WBcellsinv[n][i]:=cnt; fi;
od;
od;

##########################
ff:=function(n,i);
return WBcells[n+1][i];
end;
##########################

WBbnd:=[];
WBbnd[1]:=List(WBcells[1],i->[1,0]);
for n in [2..Length(WBcells)] do
WBbnd[n]:=List(WBcells[n],i->1*W!.boundaries[n][i]);
WBbnd[n]:=List(WBbnd[n], x-> Concatenation([x[1]],
List(x{[2..1+x[1]]},i->WBcellsinv[n-1][i])
));
od;

WB:=RegularCWComplex(WBbnd);
OrientRegularCWComplex(W);
WBorient:=[];
for n in [1..1+Dimension(WB)] do
WBorient[n]:=List([1..WB!.nrCells(n-1)],i->W!.orientation[n][ff(n-1,i)]);
od;
WB!.orientation:=WBorient;

return Objectify(HapRegularCWMap,
       rec(
           source:=WB,
           target:=W,
           mapping:=ff));

end);
####################################################


###############################################
InstallGlobalFunction(FirstHomologyCoveringCokernels,
function(f,n)
local Y, U, G, L, M, p, x, y;

Y:=Target(f);
U:=UniversalCover(Y);
G:=U!.group;
L:=LowIndexSubgroupsFpGroup(G,n);;
L:=Filtered(L,H->Index(G,H)=n);;
L:=List(L, H->EquivariantCWComplexToRegularCWMap(U,H));
L:=List(L,p->LiftedRegularCWMap(f,p));

M:=[];
for p in L do
CocriticalCellsOfRegularCWComplex(Source(p),Dimension(Source(p)));
for y in CocriticalCellsOfRegularCWComplex(Source(p),0) do
Add(M,[p,y]);
od;od;


M:=List(M,x->FundamentalGroup(x[1],x[2][2]));
M:=List(M,h-> AbelianInvariants(Target(h)/GeneratorsOfGroup(Image(h))) );


return SortedList(M);
end);
###############################################

