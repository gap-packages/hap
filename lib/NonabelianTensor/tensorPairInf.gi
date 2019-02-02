#(C) Graham Ellis, October 2005

#####################################################################
InstallGlobalFunction(NonabelianTensorProduct_Inf,
function(arg)
local
	AG, AH, SizeOrList,
	gensAG, NiceGensAG,  
        gensAH, NiceGensAH,
	G,G1, gensG, relsG, 
        H, H1, gensH, relsH,
	BG, GhomBG, BG1homF, BG2homF,
        BH, HhomBH, 
	F, relsT, gensF, gensF1, gensF2,
	AF, FhomAF,
	AGhomG, AHhomH, G1homF, G2homF, AG1homF, AG2homF,
	SF, gensSF, gensSFG, FhomSF, AFhomSF, AG1homSF, AG2homSF, SFhomAG,
	AFhomSSF,SSF,gensSF2,SSFhomSF,
	TensorSquare, delta,
	Trans,
	CrossedPairing, 
	AGhomBG, 
	BGhomAG, BHhomAH,
	i,v,w,x,y,z;

AG:=arg[1];
AH:=arg[2];

# AG and SF are groups whose elements are essentially enumerated. AG is 
# isomorphic to G and to BG. SF is equal to F/relsT and AF. Two isomorphic 
# copies of AG lie inside SF, and the homomorphisms AG1homSF, AG2homSF 
# identify the two copies. delta is the commutator map from TensorSquare to AG.
# The homomorphisms GhomBG, AGhomG, FhomSF, FhomAF, AFhomSF are all 
# isomorphisms. The relationship between the groups is summarized in the 
# following diagrams:   AG->G->BG->F->AF->SF and SF->AG.

gensAG:=GeneratorsOfGroup(AG);
NiceGensAG:=gensAG;
AGhomG:=IsomorphismFpGroup(AG);
G:=Range(AGhomG);

gensAH:=GeneratorsOfGroup(AH);
NiceGensAH:=gensAH;
AHhomH:=IsomorphismFpGroup(AH);
H:=Range(AHhomH);

gensG:=FreeGeneratorsOfFpGroup(G);
relsG:=RelatorsOfFpGroup(G);
BG:=FreeGroupOfFpGroup(G);
BGhomAG:=GroupHomomorphismByImagesNC(BG,AG, GeneratorsOfGroup(BG),gensAG);

gensH:=FreeGeneratorsOfFpGroup(H);
relsH:=RelatorsOfFpGroup(H);
BH:=FreeGroupOfFpGroup(H);
BHhomAH:=GroupHomomorphismByImagesNC(BH,AH, GeneratorsOfGroup(BH),gensAH);

F:=FreeGroup(Length(gensG)+Length(gensH));
gensF:=GeneratorsOfGroup(F); gensF1:=[]; gensF2:=[];
for i in [1..Length(gensG)] do
Append(gensF1,[gensF[i]]);
od;
for i in [1..Length(gensH)] do
Append(gensF2,[gensF[Length(gensG)+i]]);
od;

BG1homF:=GroupHomomorphismByImagesNC(BG,F,gensG,gensF1);
BG2homF:=GroupHomomorphismByImagesNC(BH,F,gensH,gensF2);
AG1homF:=GroupHomomorphismByFunction(AG,F,g->Image(BG1homF,PreImagesRepresentative(BGhomAG,g)));
AG2homF:=GroupHomomorphismByFunction(AH,F,g->Image(BG2homF,PreImagesRepresentative(BHhomAH,g)));


relsT:=[];
for x in relsG do
Append(relsT,[Image(BG1homF,x)]);
od;
for x in relsH do
Append(relsT,[Image(BG2homF,x)]);
od;


for z in NiceGensAG do
for x in gensAG do
for y in gensAH do
v:=Comm(Image(AG1homF,x),Image(AG2homF,y))^Image(AG1homF,z) ;
w:=Comm(Image(AG2homF,y^z),Image(AG1homF,x^z) );
Append(relsT,[v*w]);
od;
od;
od;

for z in NiceGensAH do
for x in gensAG do
for y in gensAH do
w:=Comm(Image(AG2homF,y^z),Image(AG1homF,x^z) );
v:=Comm(Image(AG1homF,x),Image(AG2homF,y))^Image(AG2homF,z);
Append(relsT,[v*w]);
od;
od;
od;

relsT:=SSortedList(relsT);

#####################################################################IF
AF:=F/relsT;
FhomAF:=
GroupHomomorphismByImagesNC(F,AF,GeneratorsOfGroup(F),GeneratorsOfGroup(AF));

AFhomSSF:=IsomorphismSimplifiedFpGroup(AF);
SSF:=Image(AFhomSSF);

	    SSFhomSF:=HAP_NqEpimorphismNilpotentQuotient(SSF);

SF:=Range(SSFhomSF);

gensSF2:=List(GeneratorsOfGroup(AF),x->Image(SSFhomSF,Image(AFhomSSF,x)));

AFhomSF:=GroupHomomorphismByImagesNC(AF,SF,GeneratorsOfGroup(AF),gensSF2);

FhomSF:=
GroupHomomorphismByFunction(F,SF,x->Image(AFhomSF,Image(FhomAF,x)) );

#####################################################################FI

AG1homSF:=GroupHomomorphismByFunction(AG,SF,x->Image(FhomSF,Image(AG1homF,x)));
AG2homSF:=GroupHomomorphismByFunction(AH,SF,x->Image(FhomSF,Image(AG2homF,x)));


#TensorSquare:=Intersection(
#NormalClosure(SF,Group(List(GeneratorsOfGroup(AG),x->Image(AG1homSF,x)))),
#NormalClosure(SF,Group(List(GeneratorsOfGroup(AG),x->Image(AG2homSF,x))))
#);

TensorSquare:=CommutatorSubgroup(
Group(List(GeneratorsOfGroup(AG),x->Image(AG1homSF,x))),
Group(List(GeneratorsOfGroup(AH),x->Image(AG2homSF,x))))
;



gensSF:=List(gensF,x->Image(FhomSF,x));
gensSFG:=[];
for i in [1..Length(gensAG)] do
Append(gensSFG,[gensAG[i]]);
od;
for i in [1..Length(gensAH)] do
Append(gensSFG,[gensAH[i]]);
od;

SFhomAG:=GroupHomomorphismByImagesNC(SF,AG,gensSF,gensSFG);

delta:=GroupHomomorphismByImagesNC(TensorSquare,AG,
GeneratorsOfGroup(TensorSquare),
List(GeneratorsOfGroup(TensorSquare),x->Image(SFhomAG,x)));

#####################################################################
CrossedPairing:=function(x,y)

return Comm(Image(AG1homSF,x), Image(AG2homSF,y));

end;
#####################################################################

return rec(homomorphism:=delta, pairing:=CrossedPairing);
end);
#####################################################################


