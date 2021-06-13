#(C) Graham Ellis, October 2005


#####################################################################
InstallGlobalFunction(Nil3TensorSquare,
function(arg)
local
	AG, SizeOrList,
	gensAG, NiceGensAG,  
	G, gensG, relsG, 
	BG, GhomBG, BG1homF, BG2homF,
	F, relsT, gensF, gensF1, gensF2,
	AF, FhomAF,
	AGhomG, G1homF, G2homF, AG1homF, AG2homF,
	SF, gensSF, gensSFG, FhomSF, AFhomSF, AG1homSF, AG2homSF, SFhomAG,
	AFhomSSF,SSF,gensSF2,SSFhomSF,
	TensorSquare, delta,
	Trans,
	CrossedPairing, Action,
	UpperBound,
	Todd,i,v,w,x,y,z,t;

if not IsFinite(arg[1]) then return NonabelianTensorSquare_inf(arg[1]); fi;

Todd:=16;	#Use Todd-Coxeter if Order(G)<Todd and G is not nilpotent.
#####################################################################
UpperBound:=function(AG)
local Facts, p,P,hom,bnd;

Facts:=SSortedList(Factors(Order(AG)));
bnd:=1;

for p in Facts do
P:=SylowSubgroup(AG,p);
hom:=NonabelianTensorSquare(P).homomorphism;
bnd:=bnd*Order(Source(hom))/Order(DerivedSubgroup(P));
od;

return bnd*Order(DerivedSubgroup(AG))*Order(AG)^2;
end;
#####################################################################





AG:=arg[1];
if Length(arg)>1 then SizeOrList:=arg[2]*Order(AG)^2; 
else 
	if not IsSolvable(AG) then SizeOrList:=0;
	else
    		if not IsNilpotent(AG) and Size(AG)>Todd 
		          then SizeOrList:=UpperBound(AG); fi;
	   	if not IsNilpotent(AG) and Size(AG)<=Todd then SizeOrList:=0;fi;
	    	if IsNilpotent(AG) then SizeOrList:=-1; fi;
	fi;
fi;

# AG and SF are groups whose elements are essentially enumerated. AG is 
# isomorphic to G and to BG. SF is equal to F/relsT and AF. Two isomorphic 
# copies of AG lie inside SF, and the homomorphisms AG1homSF, AG2homSF 
# identify the two copies. delta is the commutator map from TensorSquare to AG.
# The homomorphisms GhomBG, AGhomG, FhomSF, FhomAF, AFhomSF are all 
# isomorphisms. The relationship between the groups is summarized in the 
# following diagrams:   AG->G->BG->F->AF->SF and SF->AG.

gensAG:=ReduceGenerators(GeneratorsOfGroup(AG),AG);
AGhomG:=IsomorphismFpGroupByGenerators(AG,gensAG);
G:=Range(AGhomG);

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

	if IsSolvable(AG) then 
	    NiceGensAG:=Pcgs(AG);
	else
	NiceGensAG:=List(UpperCentralSeries(AG),x->GeneratorsOfGroup(x));
	NiceGensAG[1]:=[Identity(AG)];
	NiceGensAG:=Flat(NiceGensAG);
	Trans:=RightTransversal(AG,Group(NiceGensAG));
	Append(NiceGensAG,Elements(Trans));
	fi;

relsT:=[];
for x in relsG do
Append(relsT,[Image(BG1homF,x), Image(BG2homF,x)]);
od;

for z in NiceGensAG do
for x in gensAG do
for y in gensAG do
for t in NiceGensAG do
v:=Comm(Comm(Image(AG1homF,x),Image(AG2homF,y)),Image(AG1homF,t))^Image(AG1homF,z) ;
w:=Comm(Image(AG1homF,t^z),Comm(Image(AG2homF,x^z),Image(AG1homF,y^z)));
Append(relsT,[v*w]);
v:=Comm(Comm(Image(AG1homF,x),Image(AG2homF,y)),Image(AG2homF,t))^Image(AG1homF,z) ;
w:=Comm(Image(AG2homF,t^z),Comm(Image(AG2homF,x^z),Image(AG1homF,y^z)));
Append(relsT,[v*w]);
v:=Comm(Comm(Image(AG1homF,x),Image(AG2homF,y)),Image(AG1homF,t))^Image(AG2homF,z) ;
w:=Comm(Image(AG1homF,t^z),Comm(Image(AG2homF,x^z),Image(AG1homF,y^z)));
Append(relsT,[v*w]);
v:=Comm(Comm(Image(AG1homF,x),Image(AG2homF,y)),Image(AG2homF,t))^Image(AG2homF,z) ;
w:=Comm(Image(AG2homF,t^z),Comm(Image(AG2homF,x^z),Image(AG1homF,y^z)));
Append(relsT,[v*w]);
od;
od;
od;
od;

#####################################################################IF
if SizeOrList=0 then

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

AFhomSSF:=IsomorphismSimplifiedFpGroup(AF);
SSF:=Image(AFhomSSF);

	if SizeOrList=-1 then 		#if nilpotent
	    SSFhomSF:=EpimorphismNilpotentQuotient(SSF);
	#SSFhomSF:=IsomorphismPcGroup(SSF);
	else				#if solvable and big
	SSFhomSF:=EpimorphismSolvableQuotient(SSF,SizeOrList); 
	fi;

SF:=Range(SSFhomSF);

gensSF2:=List(GeneratorsOfGroup(AF),x->Image(SSFhomSF,Image(AFhomSSF,x)));

AFhomSF:=GroupHomomorphismByImagesNC(AF,SF,GeneratorsOfGroup(AF),gensSF2);

FhomSF:=
GroupHomomorphismByFunction(F,SF,x->Image(AFhomSF,Image(FhomAF,x)) );

fi;
#####################################################################FI

AG1homSF:=GroupHomomorphismByFunction(AG,SF,x->Image(FhomSF,Image(AG1homF,x)));
AG2homSF:=GroupHomomorphismByFunction(AG,SF,x->Image(FhomSF,Image(AG2homF,x)));

#TensorSquare:=NormalIntersection(
#NormalClosure(SF,Group(List(GeneratorsOfGroup(AG),x->Image(AG1homSF,x)))),
#NormalClosure(SF,Group(List(GeneratorsOfGroup(AG),x->Image(AG2homSF,x))))
#);

TensorSquare:=CommutatorSubgroup( Group(List(GeneratorsOfGroup(AG),x->Image(AG1homSF,x))),   Group(List(GeneratorsOfGroup(AG),x->Image(AG2homSF,x))) );

#Print(StructureDescription(SF/TensorSquare),"\n");
gensSF:=List(gensF,x->Image(FhomSF,x));
gensSFG:=[];
for i in [1..Length(gensAG)] do
Append(gensSFG,[gensAG[i]]);
od;
for i in [1..Length(gensAG)] do
Append(gensSFG,[gensAG[i]]);
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

#####################################################################
Action:=function(x,t)

return Image(AG1homSF,x) * t * Image(AG1homSF,x)^-1;

end;
#####################################################################


return rec(homomorphism:=delta, pairing:=CrossedPairing, action:=Action,nu:=SF);
end);
#####################################################################


