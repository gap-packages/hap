#M:=PureCubicalComplex([[1,1,1],[1,0,1],[1,1,1]]);; T:=DirectProductOfPureCubicalComplexes(M,M);; T:=DirectProductOfPureCubicalComplexes(T,M);; Y:=CubicalComplexToRegularCWSpace(T);; F:=FundamentalGroup(Y);; RelatorsOfFpGroup(F);


##########################################################
##########################################################
InstallGlobalFunction(FundamentalGroupOfRegularCWSpace,
function(arg)
local P,Y,base,e,bool, b, vertices,edges,F, r,x,w, gens, rels, 
      cells, 0cells,1cells, 2cells, 2boundaries, deform;

Y:=arg[1];
if Length(arg)>1 then base:=arg[2]; else base:=1; fi;

cells:=CriticalCellsOfRegularCWSpace(Y);
0cells:=Filtered(cells,x->x[1]=0);
Apply(0cells,x->x[2]);
1cells:=Filtered(cells,x->x[1]=1);
Apply(1cells,x->x[2]);
2cells:=Filtered(cells,x->x[1]=2);
Apply(2cells,x->x[2]);
2boundaries:=List(2cells,x->[Y!.boundaries[3][x],Y!.orientation[3][x]]);
Apply(2boundaries,x->[x[1]{[2..Length(x[1])]},x[2]]);
Apply(2boundaries,x->List([1..Length(x[1])],i->x[1][i]*x[2][i]));
deform:=ChainComplex(Y)!.deform;
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
fi;
###################################
###################################

F:=FreeGroup(Length(1cells));
gens:=GeneratorsOfGroup(F);
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

P:=PresentationFpGroup(F/rels);
SimplifyPresentation(P);;
return FpGroupPresentation(P);

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(FundamentalGroup,
"for regular CW-spaces",
[IsHapRegularCWSpace],
function(Y)
return FundamentalGroupOfRegularCWSpace(Y);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(FundamentalGroup,
"for regular CW-spaces",
[IsHapRegularCWSpace,IsInt],
function(Y,n)
return FundamentalGroupOfRegularCWSpace(Y,n);
end);
##########################################################
##########################################################


