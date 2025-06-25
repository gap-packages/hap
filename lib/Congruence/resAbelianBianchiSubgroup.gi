
#######################################################
#######################################################
InstallGlobalFunction(ResolutionAbelianBianchiSubgroup,
function(G,n)
local R, A, A1, A2, D, p1, p2, f1, f2, f, gensA1, gensA2, gens, L, LL;

R:=ResolutionAbelianGroup([2,0,0],n);
A:=R!.group;
D:=A!.DirectProductInfo;
p1:=D!.projections[1];
p2:=D!.projections[2];
A1:=D!.groups[1];
A2:=D!.groups[2];
gensA1:=GeneratorsOfGroup(A1);
gensA2:=GeneratorsOfGroup(A2);
gens:=GeneratorsOfGroup(G);
f1:=GroupHomomorphismByImagesNC(A1,G,gensA1,gens{[3]});
f2:=GroupHomomorphismByImagesNC(A2,G,gensA2,gens{[1,2]});

f:=GroupHomomorphismByFunction(A,G,x-> Image(f1,Image(p1,x))*Image(f2,Image(p2,x)));

L:=R!.elts;
LL:=List([1..1000],i->Image(f,L[i])); #SLOPPY! This should be a lazy list
R!.elts:=LL;
R!.group:=G;

return R;
end);
#######################################################
#######################################################

