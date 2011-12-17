#(C) Graham Ellis, October 2005

#####################################################################
InstallGlobalFunction(NonabelianTensorSquare,
function(AG)
local
	gensAG, NiceGensAG, 
	G, gensG, relsG, 
	BG, GhomBG, BG1homF, BG2homF,
	F, relsT, gensF, gensF1, gensF2,
	AF, FhomAF,
	AGhomG, G1homF, G2homF, AG1homF, AG2homF,
	SF, gensSF, gensSFG, FhomSF, AFhomSF, AG1homSF, AG2homSF, SFhomAG,
	TensorSquare, delta,
	Trans,
	CrossedPairing,
	i,v,w,x,y,z;

# AG and SF are groups whose elements are essentially enumerated. AG is 
# isomorphic to G and to BG. SF is equal to F/relsT and AF. Two isomorphic 
# copies of AG lie inside SF, and the homomorphisms AG1homSF, AG2homSF 
# identify the two copies. delta is the commutator map from TensorSquare to AG.
# The homomorphisms GhomBG, AGhomG, FhomSF, FhomAF, AFhomSF are all 
# isomorphisms. The relationship between the groups is summarized in the 
# following diagrams:   AG->G->BG->F->AF->SF and SF->AG.

gensAG:=GeneratorsOfGroup(AG);
AGhomG:=IsomorphismFpGroupByGenerators(AG,gensAG);	
G:=Image(AGhomG);
gensG:=FreeGeneratorsOfFpGroup(G);
relsG:=RelatorsOfFpGroup(G);
BG:=Group(gensG);
GhomBG:=GroupHomomorphismByImagesNC(G,BG, GeneratorsOfGroup(G),gensG);
			#I hope GhomBG really is the identity map!

F:=FreeGroup(2*Length(gensG));
gensF:=GeneratorsOfGroup(F); gensF1:=[]; gensF2:=[];
for i in [1..Length(gensG)] do
Append(gensF1,[gensF[i]]);
Append(gensF2,[gensF[Length(gensG)+i]]);
od;

BG1homF:=GroupHomomorphismByImagesNC(BG,F,gensG,gensF1);
G1homF:=GroupHomomorphismByFunction(G,F,x->Image(BG1homF,Image(GhomBG,x)));
BG2homF:=GroupHomomorphismByImagesNC(BG,F,gensG,gensF2);
G2homF:=GroupHomomorphismByFunction(G,F,x->Image(BG2homF,Image(GhomBG,x)));
AG1homF:=GroupHomomorphismByFunction(AG,F,g->Image(G1homF,Image(AGhomG,g)));
AG2homF:=GroupHomomorphismByFunction(AG,F,g->Image(G2homF,Image(AGhomG,g)));

NiceGensAG:=List(UpperCentralSeries(AG),x->GeneratorsOfGroup(x));
NiceGensAG[1]:=[Identity(AG)];
NiceGensAG:=Flat(NiceGensAG);
Trans:=RightTransversal(AG,Group(NiceGensAG));
Append(NiceGensAG,Elements(Trans));

relsT:=[];
for x in relsG do
Append(relsT,[Image(BG1homF,x), Image(BG2homF,x)]);
od;

for z in GeneratorsOfGroup(AG) do
for x in NiceGensAG do
for y in NiceGensAG do
v:=Comm(Image(AG1homF,x),Image(AG2homF,y))^Image(AG1homF,z) ;
w:=Comm(Image(AG2homF,y^z),Image(AG1homF,x^z) );
Append(relsT,[v*w]);
v:=Comm(Image(AG1homF,x),Image(AG2homF,y))^Image(AG2homF,z);
Append(relsT,[v*w]);
od;
od;
od;

#####################################################################IF
if not IsNilpotent(AG) then

AF:=F/relsT;
FhomAF:=
GroupHomomorphismByImages(F,AF,GeneratorsOfGroup(F),GeneratorsOfGroup(AF));

AFhomSF:=IsomorphismSimplifiedFpGroup(AF);
SF:=Image(AFhomSF);
FhomSF:=
GroupHomomorphismByFunction(F,SF,x->Image(AFhomSF,Image(FhomAF,x)) );


else

AF:=F/relsT;
FhomAF:=
GroupHomomorphismByImages(F,AF,GeneratorsOfGroup(F),GeneratorsOfGroup(AF));

AFhomSF:=EpimorphismNilpotentQuotient(AF);
SF:=Image(AFhomSF);
FhomSF:=
GroupHomomorphismByFunction(F,SF,x->Image(AFhomSF,Image(FhomAF,x)) );

fi;
#####################################################################FI

AG1homSF:=GroupHomomorphismByFunction(AG,SF,x->Image(FhomSF,Image(AG1homF,x)));
AG2homSF:=GroupHomomorphismByFunction(AG,SF,x->Image(FhomSF,Image(AG2homF,x)));

TensorSquare:=Intersection(
NormalClosure(SF,Group(List(GeneratorsOfGroup(AG),x->Image(AG1homSF,x)))),
NormalClosure(SF,Group(List(GeneratorsOfGroup(AG),x->Image(AG2homSF,x))))
);


gensSF:=List(gensF,x->Image(FhomSF,x));
gensSFG:=[];
for i in [1..Length(gensAG)] do
Append(gensSFG,[gensAG[i]]);
od;
for i in [1..Length(gensAG)] do
Append(gensSFG,[gensAG[i]]);
od;

SFhomAG:=GroupHomomorphismByImages(SF,AG,gensSF,gensSFG);

delta:=GroupHomomorphismByFunction(TensorSquare,AG,x->Image(SFhomAG,x));

#####################################################################
CrossedPairing:=function(x,y)

return Comm(Image(AG1homSF,x), Image(AG2homSF,y));

end;
#####################################################################

return rec(homomorphism:=delta, pairing:=CrossedPairing);
end);
#####################################################################

#####################################################################
InstallGlobalFunction(ThirdHomotopyGroupOfSuspensionB,
function(G);

return AbelianInvariants(Kernel(
			NonabelianTensorSquare(G).homomorphism));

end);
#####################################################################

