#(C) Graham Ellis, October 2005

#####################################################################
InstallMethod(NonabelianExteriorProduct,
"Exterior produce of a group AG and normal subgroup AH",
[IsGroup,IsGroup],
function(AG,AH)
local
	gensAG, NiceGensAG, 
	G, gensG, relsG, 
	BG, GhomBG, BG1homF, 
	gensAH, NiceGensAH,
        H, gensH, relsH,
        BH, HhomBH, AHhomH,
	BG2homF,
	F, relsT, gensF, gensF1, gensF2,
	AF, FhomAF,
	AGhomG, G1homF, G2homF, AG1homF, AG2homF,
	SF, gensSF, gensSFG, FhomSF, AFhomSF, AG1homSF, AG2homSF, SFhomAG,
	ExteriorProduct, delta,
	Trans,
	CrossedPairing,
	i,v,w,x,y,z;

# This function is an adaption of the NonabelianTensorSquare() function and
# so names of variables are not always the most sensible. For instance, G1
# stands for G, and G2 stands for H.

# The group H is a normal subgroup of G.

# AG and SF are groups whose elements are essentially enumerated. AG is 
# isomorphic to G and to BG. SF is equal to F/relsT and AF. Two isomorphic 
# copies of AG lie inside SF, and the homomorphisms AG1homSF, AG2homSF 
# identify the two copies. delta is the commutator map from TensorSquare to AG.
# The homomorphisms GhomBG, AGhomG, FhomSF, FhomAF, AFhomSF are all 
# isomorphisms. The relationship between the groups is summarized in the 
# following diagrams:   AG->G->BG->F->AF->SF and SF->AG.

if not IsSubgroup(AG,AH) then
Print("The second argument should be a normal subgroup of the first.\n");
return fail;
fi;

if not IsNormal(AG,AH) then
Print("The second argument should be a normal subgroup of the first.\n");
return fail;
fi;


gensAG:=GeneratorsOfGroup(AG);
AGhomG:=IsomorphismFpGroupByGenerators(AG,gensAG);	
G:=Image(AGhomG);
gensG:=FreeGeneratorsOfFpGroup(G);
relsG:=RelatorsOfFpGroup(G);
BG:=Group(gensG);
GhomBG:=GroupHomomorphismByImagesNC(G,BG, GeneratorsOfGroup(G),gensG);
			#I hope GhomBG really is the identity map!


gensAH:=GeneratorsOfGroup(AH);
AHhomH:=IsomorphismFpGroupByGenerators(AH,gensAH);
H:=Image(AHhomH);
gensH:=FreeGeneratorsOfFpGroup(H);
relsH:=RelatorsOfFpGroup(H);
BH:=Group(gensH);
HhomBH:=GroupHomomorphismByImagesNC(H,BH, GeneratorsOfGroup(H),gensH);
                        #I hope HhomBH really is the identity map!

F:=FreeGroup(Length(gensG)+Length(gensH));
gensF:=GeneratorsOfGroup(F); gensF1:=[]; gensF2:=[];
for i in [1..Length(gensG)] do
Append(gensF1,[gensF[i]]);
od;
for i in [1..Length(gensH)] do
Append(gensF2,[gensF[Length(gensG)+i]]);
od;

BG1homF:=GroupHomomorphismByImagesNC(BG,F,gensG,gensF1);
G1homF:=GroupHomomorphismByFunction(G,F,x->Image(BG1homF,Image(GhomBG,x)));
BG2homF:=GroupHomomorphismByImagesNC(BH,F,gensH,gensF2);
G2homF:=GroupHomomorphismByFunction(H,F,x->Image(BG2homF,Image(HhomBH,x)));
AG1homF:=GroupHomomorphismByFunction(AG,F,g->Image(G1homF,Image(AGhomG,g)));
AG2homF:=GroupHomomorphismByFunction(AH,F,g->Image(G2homF,Image(AHhomH,g)));

if IsSolvable(AG) then #
NiceGensAG:=Pcgs(AG); #Need to check the maths here!
         else          #

NiceGensAG:=List(UpperCentralSeries(AG),x->GeneratorsOfGroup(x));
NiceGensAG[1]:=[Identity(AG)];
NiceGensAG:=Flat(NiceGensAG);
Trans:=RightTransversal(AG,Group(NiceGensAG));
Append(NiceGensAG,Elements(Trans));
fi;

if IsSolvable(AH) then #
NiceGensAH:=Pcgs(AH); #Need to check the maths here!
         else          #

NiceGensAH:=
List(UpperCentralSeries(AG),x->GeneratorsOfGroup(x));
NiceGensAH[1]:=[Identity(AH)];
NiceGensAH:=Flat(NiceGensAH);
NiceGensAH:=Filtered(NiceGensAH,x-> (x in AH));
Trans:=RightTransversal(AH,Group(NiceGensAH));
Append(NiceGensAH,Elements(Trans));
fi;


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


for y in AH do
v:=Comm(Image(AG1homF,y),Image(AG2homF,y));
Append(relsT,[v]);
od;

#####################################################################IF
if not IsSolvable(AG) then

AF:=F/relsT;
FhomAF:=
GroupHomomorphismByImagesNC(F,AF,GeneratorsOfGroup(F),GeneratorsOfGroup(AF));

AFhomSF:=IsomorphismSimplifiedFpGroup(AF);
SF:=Image(AFhomSF);
FhomSF:=
GroupHomomorphismByFunction(F,SF,x->Image(AFhomSF,Image(FhomAF,x)) );


else

AF:=F/relsT;
FhomAF:=
GroupHomomorphismByImagesNC(F,AF,GeneratorsOfGroup(F),GeneratorsOfGroup(AF));

if IsNilpotent(AG) then
AFhomSF:=EpimorphismNilpotentQuotient(AF);
else
AFhomSF:=EpimorphismSolvableQuotient(AF, SSortedList(Factors(Order(G))));
fi;

SF:=Image(AFhomSF);
FhomSF:=
GroupHomomorphismByFunction(F,SF,x->Image(AFhomSF,Image(FhomAF,x)) );

fi;
#####################################################################FI


AG1homSF:=GroupHomomorphismByFunction(AG,SF,x->Image(FhomSF,Image(AG1homF,x)));
AG2homSF:=GroupHomomorphismByFunction(AH,SF,x->Image(FhomSF,Image(AG2homF,x)));

ExteriorProduct:=Intersection(
NormalClosure(SF,Group(List(GeneratorsOfGroup(AG),x->Image(AG1homSF,x)))),
NormalClosure(SF,Group(List(GeneratorsOfGroup(AH),x->Image(AG2homSF,x))))
);


gensSF:=List(gensF,x->Image(FhomSF,x));
gensSFG:=[];
for i in [1..Length(gensAG)] do
Append(gensSFG,[gensAG[i]]);
od;
for i in [1..Length(gensAH)] do
Append(gensSFG,[gensAH[i]]);
od;

SFhomAG:=GroupHomomorphismByImagesNC(SF,AG,gensSF,gensSFG);

delta:=GroupHomomorphismByFunction(ExteriorProduct,AG,x->Image(SFhomAG,x));

#####################################################################
CrossedPairing:=function(x,y)

return Comm(Image(AG1homSF,x), Image(AG2homSF,y));

end;
#####################################################################

return rec(homomorphism:=delta, pairing:=CrossedPairing);
end);
#####################################################################


#####################################################################
InstallGlobalFunction(RelativeSchurMultiplier,
function(G,N);

if G=N  then return GroupHomology(G,2); fi;

return AbelianInvariants(Kernel
          (NonabelianExteriorProduct(G,N).homomorphism));

end);
#####################################################################

#####################################################################
InstallGlobalFunction(EpiCentre,
function(arg)
local 		G,N,gensG, Epi, Pairing,toggle,x,z;

G:=arg[1];
if Length(arg)>1 then N:=arg[2]; else 

if Order(Centre(G))=1 then return Centre(G); fi;
if IsNilpotent(G) and LoadPackage("nq")=true then
return UpperEpicentralSeries(G,1); fi;

N:=G; fi;

if IsTrivial(Centre(N)) then return Centre(N); fi;

gensG:=GeneratorsOfGroup(G);

Pairing:=NonabelianExteriorProduct(G,N).pairing;

Epi:=[];

for z in Center(N) do
toggle:=true;
for x in gensG do
if Order(Pairing(x,z))>1 then toggle:=false; break; fi;
od;
if toggle then Append(Epi,[z]); fi;
od;

return Group(Epi);
end);
#####################################################################
