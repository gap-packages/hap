
#######################################################
#######################################################
InstallGlobalFunction(WeakCommutativityGroup,
function(GG)
local G, epi, F, gensF, relsG, FF, n, FhomFF1, FhomFF2, gensFF, relsX, 
       g, x, i, j, y, iso, Chi;


iso:=IsomorphismFpGroup(GG);
G:=Image(iso);
F:=FreeGroupOfFpGroup(G);
gensF:=GeneratorsOfGroup(F);
epi:=GroupHomomorphismByImagesNC(F,G,gensF, GeneratorsOfGroup(G));
relsG:=RelatorsOfFpGroup(G);

n:=Length(gensF);
FF:=FreeGroup(2*n);
gensFF:=GeneratorsOfGroup(FF);

FhomFF1:=GroupHomomorphismByImages(F,FF,gensF,gensFF{[1..n]});
FhomFF2:=GroupHomomorphismByImages(F,FF,gensF,gensFF{[n+1..2*n]});

relsX:=List(relsG,x->Image(FhomFF1,x));
relsX:=Concatenation(relsX,List(relsG,x->Image(FhomFF2,x)) );

for g in GG do
x:=ImagesRepresentative(iso,g);
x:=PreImagesRepresentative(epi,x);
Add(relsX, Comm(Image(FhomFF1,x), Image(FhomFF2,x) ) );
od;


Chi:= FF/relsX;
Chi:=SimplifiedFpGroup(Chi);
#if IsPGroup(GG) then Chi:=Image(NqEpimorphismNilpotentQuotient(Chi)); fi;
if IsNilpotentGroup(GG) then Chi:=Image(NqEpimorphismNilpotentQuotient(Chi)); fi;

return Chi;
end);
#######################################################
#######################################################


#######################################################
#######################################################
InstallGlobalFunction(SymmetricCommutativityGroup,
function(GG)
local G, epi, F, gensF, relsG, FF, n, FhomFF1, FhomFF2, gensFF, relsX,
       g, h, x, i, j, y, iso, Chi;


iso:=IsomorphismFpGroup(Group(MinimalGeneratingSet(GG)));
G:=Image(iso);
F:=FreeGroupOfFpGroup(G);
gensF:=GeneratorsOfGroup(F);
epi:=GroupHomomorphismByImagesNC(F,G,gensF, GeneratorsOfGroup(G));
relsG:=RelatorsOfFpGroup(G);

n:=Length(gensF);
FF:=FreeGroup(2*n);
gensFF:=GeneratorsOfGroup(FF);

FhomFF1:=GroupHomomorphismByImages(F,FF,gensF,gensFF{[1..n]});
FhomFF2:=GroupHomomorphismByImages(F,FF,gensF,gensFF{[n+1..2*n]});

relsX:=List(relsG,x->Image(FhomFF1,x));
relsX:=Concatenation(relsX,List(relsG,x->Image(FhomFF2,x)) );

for g in GG do
for h in GG do
x:=ImagesRepresentative(iso,g);
x:=PreImagesRepresentative(epi,x);
y:=ImagesRepresentative(iso,h);
y:=PreImagesRepresentative(epi,y);
Add(relsX, Comm(Image(FhomFF1,x), Image(FhomFF2,y) ) * Comm(Image(FhomFF1,y), Image(FhomFF2,x) ) );
od;od;


Chi:= FF/relsX;
#Chi:=SimplifiedFpGroup(Chi);

return Chi;
end);
#######################################################
#######################################################


