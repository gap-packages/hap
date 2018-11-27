###############################################
###############################################
InstallGlobalFunction(HomToGModule_hom,
function(F,A)
local R,C,grps,maps, properties, HomF, lnth, n, B,S, gensS, 
      pid, imS, hom, fn, act, gouts;

R:=Source(F);
C:=HomToGModule(R,A);
lnth:=EvaluateProperty(R,"length");

################################
grps:=[];gouts:=[];
for n in [0..lnth-1] do
B:=C!.boundary(n);
if n=0 then S:=Source(B); 
gouts[1]:=S;
grps[1]:=[S!.ActedGroup,S!.OuterAction]; fi;
S:=Target(B); 
gouts[n+2]:=S;
grps[n+2]:=[S!.ActedGroup,S!.OuterAction];
od;
################################

pid:=Position(R!.elts,One(R!.group));
################################
fn:=function(S,act,x,n)   ##This is a sloppy piece of code!
local w,y,i,fwrd, u;
w:=One(S);
for i in [1..R!.dimension(n)] do
fwrd:=F!.mapping([[i,pid]],n);
for y in fwrd do
u:=act(R!.elts[y[2]],x);
u:=Image(Projection(S,AbsInt(y[1])),u^SignInt(y[1]));
u:=Image(Embedding(S,i),u);
w:=u*w;
od;
od;
return w;
end;
################################

################################
maps:=[];
for n in [1..Length(grps)-1] do
S:=grps[n][1]; act:=grps[n][2];
gensS:=GeneratorsOfGroup(S);
imS:=List(gensS,x->fn(S,act,x,n-1));
hom:=GroupHomomorphismByImages(S,S,gensS,imS);
hom:=GOuterGroupHomomorphism(gouts[n],gouts[n],hom);
maps[n]:=hom;
od;
################################

properties:=[["length",lnth],["type","GcomplexMap"]];

HomF:=rec(source:=C,target:=C,map:=maps,properties:=properties);

return Objectify(HapGComplexMap,HomF);
end);
###############################################
###############################################

###############################################
###############################################
InstallOtherMethod(HomToGModule ,
        "for GOuterGroup endomorphisms",
        [IsHapEquivariantChainMap, IsGOuterGroup],
    function(F,A)
    return HomToGModule_hom(F,A);
end);
###############################################
###############################################

###############################################
###############################################
InstallOtherMethod(Source ,
        "for GComplex map",
        [IsHapGComplexMap],
    function(F)
    return F!.source;
end);
###############################################
###############################################

###############################################
###############################################
InstallOtherMethod(Target ,
        "for GComplex map",
        [IsHapGComplexMap],
    function(F)
    return F!.target;
end);
###############################################
###############################################



###############################################
###############################################
InstallGlobalFunction(CohomologyModule_Gmap,
function(F,n)
local C,D,HC,HD,Hom,hom,NatC, NatD, natc, natd, ind, c, d, imc, gensc,x;

C:=Source(F);
D:=Target(F);
HC:=CohomologyModule(C,n);
HD:=CohomologyModule(D,n);
Hom:=F!.map[n+1];
hom:=Hom!.Mapping;
NatC:=HC!.nat;
NatD:=HD!.nat;
natc:=NatC!.Mapping;
natd:=NatD!.Mapping;

c:=Target(natc);
d:=Target(natd);
gensc:=GeneratorsOfGroup(c);
imc:=[];

for x in gensc do
Add(imc,Image(natd,Image(hom,PreImagesRepresentative(natc,x)))  );
od;

ind:=GroupHomomorphismByImages(c,d,gensc,imc);

return GOuterGroupHomomorphism(HC,HD,ind);
end);
###############################################
###############################################

###############################################
###############################################
InstallOtherMethod(CohomologyModule ,
        "for GComplex maps",
        [IsHapGComplexMap, IsInt],
    function(F,n)
    return CohomologyModule_Gmap(F,n);
end);
###############################################
###############################################


###############################################
###############################################
InstallGlobalFunction(CommonEndomorphisms,
function(G, L)
local aut, L1, G1, GisoG1, G1isoG, i;

#Inputs a group G and a list of endomorphisms a_i:Hi-->Hi
#where  H1, H2, ... are all isomorphic to G.
#Returns the corresponding list of endomorphisms a_i:G-->G

L1:=[];
for i in [1..Length(L)] do
G1:=Source(L[i]);
GisoG1:=IsomorphismGroups(G,G1);
G1isoG:=InverseGeneralMapping(GisoG1);
Add(L1,GisoG1*L[i]*G1isoG);
od;

return L1;

end);
###############################################
###############################################

###############################################
###############################################
InstallGlobalFunction(CohomologyModule_AsAutModule,
function(C,aut,n)
local R,A, gens, gens1, H, CH, AUT, hom; 

if not IsBound(C!.resolution) then
Print("No resolution is associated to the cocomplex.\n");
return fail;
fi;

CH:=CohomologyModule(C,n);
H:=CH!.ActedGroup;;

R:=C!.resolution;
A:=C!.module;
gens1:=SmallGeneratingSet(aut);;
gens:=List(gens1,h-> EquivariantChainMap(R,R,h));;
gens:=List(gens,F->HomToGModule(F,A));;
gens:=List(gens,a->CohomologyModule(a,n));;
gens:=List(gens,x->x!.Mapping);;
gens:=CommonEndomorphisms(H,gens);;
AUT:=Group(gens);;

hom:=GroupHomomorphismByImages(aut,AUT,gens1, List(gens,x->x^-1));

CH!.ActingGroup:=aut;

#############################
CH!.OuterAction:=function(g,a);
return Image(Image(hom,g),a);
end;
#############################

return CH; 
end);
###############################################
###############################################

###############################################
###############################################
InstallOtherMethod(CohomologyModule ,
        "for GComplex maps",
        [IsHapGCocomplex, IsGroup, IsInt],
    function(C,aut,n)
    return CohomologyModule_AsAutModule(C, aut,n);
end);
###############################################
###############################################


