#(C) Graham Ellis, October 2005

#####################################################################
InstallGlobalFunction(NonabelianSymmetricSquare_inf,
function(arg)
local
	AG, SizeOrList,
	gensAG, NiceGensAG,  
	G,G1, gensG, relsG, 
	BG, GhomBG, BG1homF, BG2homF,
	F, relsT, gensF, gensF1, gensF2,
	AF, FhomAF,
	AGhomG, G1homF, G2homF, AG1homF, AG2homF,
	SF, gensSF, gensSFG, FhomSF, AFhomSF, AG1homSF, AG2homSF, SFhomAG,
	AFhomSSF,SSF,gensSF2,SSFhomSF,
	SymmetricSquare, delta,
	Trans,
	CrossedPairing, 
	AGhomBG,
	BGhomAG,
	i,v,w,x,y,z;






AG:=arg[1];

# AG and SF are groups whose elements are essentially enumerated. AG is 
# isomorphic to G and to BG. SF is equal to F/relsT and AF. Two isomorphic 
# copies of AG lie inside SF, and the homomorphisms AG1homSF, AG2homSF 
# identify the two copies. delta is the commutator map from SymmetricSquare to AG.
# The homomorphisms GhomBG, AGhomG, FhomSF, FhomAF, AFhomSF are all 
# isomorphisms. The relationship between the groups is summarized in the 
# following diagrams:   AG->G->BG->F->AF->SF and SF->AG.

gensAG:=GeneratorsOfGroup(AG);
NiceGensAG:=gensAG;
AGhomG:=IsomorphismFpGroup(AG);
G:=Range(AGhomG);


gensG:=FreeGeneratorsOfFpGroup(G);
relsG:=RelatorsOfFpGroup(G);
BG:=FreeGroupOfFpGroup(G);
			#I hope GhomBG really is the identity map!
BGhomAG:=GroupHomomorphismByImagesNC(BG,AG, GeneratorsOfGroup(BG),gensAG);
F:=FreeGroup(2*Length(gensG));
gensF:=GeneratorsOfGroup(F); gensF1:=[]; gensF2:=[];
for i in [1..Length(gensG)] do
Append(gensF1,[gensF[i]]);
Append(gensF2,[gensF[Length(gensG)+i]]);
od;

BG1homF:=GroupHomomorphismByImagesNC(BG,F,gensG,gensF1);
BG2homF:=GroupHomomorphismByImagesNC(BG,F,gensG,gensF2);
AG1homF:=GroupHomomorphismByFunction(AG,F,g->Image(BG1homF,PreImagesRepresentative(BGhomAG,g)));
AG2homF:=GroupHomomorphismByFunction(AG,F,g->Image(BG2homF,PreImagesRepresentative(BGhomAG,g)));


relsT:=[];
for x in relsG do
Append(relsT,[Image(BG1homF,x), Image(BG2homF,x)]);
od;

#for z in GeneratorsOfGroup(AG) do
#for x in NiceGensAG do
#for y in NiceGensAG do
for z in NiceGensAG do
for x in gensAG do
for y in gensAG do
v:=Comm(Image(AG1homF,x),Image(AG2homF,y))^Image(AG1homF,z) ;
w:=Comm(Image(AG2homF,y^z),Image(AG1homF,x^z) );
Append(relsT,[v*w]);
v:=Comm(Image(AG1homF,x),Image(AG2homF,y))^Image(AG2homF,z);
Append(relsT,[v*w]);
od;
od;
od;

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
AG2homSF:=GroupHomomorphismByFunction(AG,SF,x->Image(FhomSF,Image(AG2homF,x)));



SymmetricSquare:=CommutatorSubgroup(
Group(List(GeneratorsOfGroup(AG),x->Image(AG1homSF,x))),
Group(List(GeneratorsOfGroup(AG),x->Image(AG2homSF,x))))
;



gensSF:=List(gensF,x->Image(FhomSF,x));
gensSFG:=[];
for i in [1..Length(gensAG)] do
Append(gensSFG,[gensAG[i]]);
od;
for i in [1..Length(gensAG)] do
Append(gensSFG,[gensAG[i]]);
od;

SFhomAG:=GroupHomomorphismByImagesNC(SF,AG,gensSF,gensSFG);

delta:=GroupHomomorphismByImagesNC(SymmetricSquare,AG,
GeneratorsOfGroup(SymmetricSquare),
List(GeneratorsOfGroup(SymmetricSquare),x->Image(SFhomAG,x)));

#####################################################################
CrossedPairing:=function(x,y)

return Comm(Image(AG1homSF,x), Image(AG2homSF,y));

end;
#####################################################################

return rec(homomorphism:=delta, pairing:=CrossedPairing);
end);
#####################################################################


