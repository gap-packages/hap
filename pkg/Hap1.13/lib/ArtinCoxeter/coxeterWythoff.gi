#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(CoxeterWythoffComplex,
function(D,B,K)
local
	Dimension,
	Boundary,       
	Contraction,	#not yet used
	EltsG, EltsG1,	
	Vertices,
	G, gensG, Glist, G1,	#G is the Artin group of D. We treat G as a
 	GhomG1,	gensG1,		#free group. However, we output a copy G1 of
				#G with relators. 
	W, Wgens, GhomW,	#W is the Coxeter group of D
	ResGens,
	BoundaryCoeff,
	PseudoBoundary,
	BoundaryRecord,
	CoxeterDiagramInterval,
	CoxeterDiagramIsBlocking,
	CoxeterDiagramIsEssential,
	NonEssentialSets,
	EssentialSet,
	Vcomplement,
	Triv,
	AhomW,WhomWP,WP,WPev,AhomWP,EltsWP,A,StabilizerSubgroup,Action,
	EvenStabGroup,
	m, n,k,i,c,x,U, S, SD;

###########################
if not CoxeterDiagramIsSpherical(D) then
Print("This function is only implemented for finite Coxeter groups.\n");
return fail;
fi;
###########################

Vertices:=CoxeterDiagramVertices(D);

#######################
if not IsSubset(Vertices,B) then
Print("The specified vertices do not all lie in the vertex set of the Coxeter diagram. \n"); return fail; fi;
#######################



#######CREATE FUNCTIONS FOR FINDING ESSENTIAL SETS############
#######                                           ############
#######(The efficiency of these could easily be improved)#####

##############################################################
CoxeterDiagramInterval:=function(D,u,v)
#Inputs a diagram D and two vertices u,v.
#Outputs the list of vertices in some path from u to v.
local Paths, NewPaths, P, M, Vertices,i,j;

if u=v then return [u]; fi;

Vertices:=CoxeterDiagramVertices(D);
M:=CoxeterDiagramMatrix(D);

Paths:=[[u]];
NewPaths:=[];

while true do
for P in Paths do
for i in Vertices do
if M[P[Length(P)]][i]>2 and not i in P then
if i=v then return Concatenation(P,[i]);fi;
Add(NewPaths,Concatenation(P,[i])); fi;
od;
od;
Paths:=NewPaths; NewPaths:=[];
od;

end;
#############################################################

#############################################################
CoxeterDiagramIsBlocking:=function(D,V,U1,U)
local u,v;

for u in U do
for v in V do
if Length(Intersection(
CoxeterDiagramInterval(D,u,v),
U1))=0
then return false; fi;
od;
od;

return true;
end;
#############################################################

NonEssentialSets:=[];
for i in [1..Length(Vertices)-1] do
NonEssentialSets[i]:=[];
for c in Combinations(Vertices,i+1) do
for x in c do
U:=Filtered(c,j->not j=x);
if CoxeterDiagramIsBlocking(D,B,U,c) then
Add(NonEssentialSets[i],U); fi;
od;
od;
od;

#Print(NonEssentialSets,"\n");

##############################################################
CoxeterDiagramIsEssential:=function(T);
if Length(T)=Length(Vertices) then return true; fi;
if Length(T)=0 then return true; fi;
if T in NonEssentialSets[Length(T)] then return false; fi;
return true;
end;
############################################################## 

##############################################################
EssentialSet:=function(U)
local u,U1;;

U1:=SSortedList(U);

for u in U do
RemoveSet(U1,u);
if not CoxeterDiagramIsBlocking(D,B,U1,U) then
AddSet(U1,u); fi;
od;

return  U1;
end;
##############################################################

##############################################################
Vcomplement:=function(S);
return Filtered(Vertices,x->not x in S);
end;
##############################################################

######ESSENTIAL SET FUNCTIONS NOW CREATED#####################
######                                   #####################
##############################################################


Glist:=CoxeterDiagramFpArtinGroup(D);
G1:=Glist[1]/Glist[2];        	#Take care for this not to cause Knuth-Bendix 
gensG1:=GeneratorsOfGroup(G1);	#to start up later on!
G:=Glist[1];
gensG:=GeneratorsOfGroup(G);
GhomG1:=GroupHomomorphismByImagesNC(G,G1,gensG,gensG1);
EltsG:=[];
EltsG1:=[];


ResGens:=[];
ResGens[1]:=[[]];
for n in [1..K] do
ResGens[n+1]:=[];
for S in Combinations(Vertices,n) do
if  CoxeterDiagramIsEssential(Vcomplement(S)) then AddSet(ResGens[n+1],S); fi;
od;
od;

#####################################################################
Dimension:=function(n);

if n=0 then return 1;
else return Length(ResGens[n+1]); fi;

end;
#####################################################################

BoundaryRecord:=[];
for n in [1..K] do
BoundaryRecord[n]:=[];
for m in [1..Dimension(n)] do
BoundaryRecord[n][m]:=true;
od;
od;

#####################################################################
BoundaryCoeff:=function(S,T)	#S is a set of vertices generating a
				#finite Coxeter group WS. T is a
				#subset of S, and WT is the corresponding
				#subgroup of WS.
local 	SD, WS, gensWS, 
	WT, gensWT,
	Trans, 
	WShomG, Ggens,
	x,y;

SD:=CoxeterSubDiagram(D,S);
WS:=CoxeterDiagramFpCoxeterGroup(SD);
WS:=WS[1]/WS[2];
Ggens:=List(S,x->gensG[Position(Vertices,x)]);
gensWS:=GeneratorsOfGroup(WS);
WShomG:=GroupHomomorphismByImagesNC(WS,G,gensWS,Ggens);
gensWT:=List(T,x->gensWS[Position(S,x)]);
if Length(T)>0 then WT:=Group(gensWT);
else WT:=Group(Identity(WS)); fi;

Trans:=List(Elements(RightTransversal(WS,WT)),x->x^-1);

for x in Trans do
y:=Image(WShomG,x);
if not y in EltsG then Append(EltsG,[y]); 
y:=Image(GhomG1,y); 
Append(EltsG1,[y]);
fi;
od;

return List(Trans,x->Image(WShomG,x));
end;
#####################################################################

#####################################################################
PseudoBoundary:=function(S)	#S is a subset of vertices with finite
				#Coxeter group WS.
local T, bndry, a;

bndry:=[];
for T in Combinations(S,Length(S)-1) do
if  CoxeterDiagramIsEssential(Vcomplement(T)) then
a:=Difference(S,T)[1];
Append(bndry,[  [T,BoundaryCoeff(S,T),Position(S,a)]  ]);
fi;
od;

return bndry;
end;
#####################################################################

#####################################################################
Boundary:=function(n,kk)
local B, B1, FreeGWord, x, y, k;

#n:=AbsoluteValue(m);
if n<1 then return 0; fi;

k:=AbsoluteValue(kk);

if not BoundaryRecord[n][k]=true then 
if kk>0 then return BoundaryRecord[n][k]; 
else return NegateWord(BoundaryRecord[n][k]);fi;
fi;

B:=PseudoBoundary(ResGens[n+1][k]);
B1:=List(B,x->[Position(ResGens[n],x[1]),
	List(x[2],y->(-1)^(Length(y)+x[3])*Position(EltsG,y))  ]);
FreeGWord:=[];
for x in B1 do
for y in x[2] do
Append(FreeGWord,[ [SignInt(y)*x[1],AbsoluteValue(y)] ]);
od;
od;

BoundaryRecord[n][k]:=FreeGWord; 
if kk>0 then return FreeGWord;
else return NegateWord(FreeGWord); fi; 
end;
#####################################################################

####From here on we cut and past from CoxeterComplex and so need
#a change of notation!!
A:=G1;


for n in [1..K] do
for k in [1..Dimension(n)] do
i:=Boundary(n,k);
od; od;

W:=CoxeterDiagramFpCoxeterGroup(D);
W:=W[1]/W[2];
AhomW:=GroupHomomorphismByImagesNC(A,W,GeneratorsOfGroup(A),GeneratorsOfGroup(W));
WhomWP:=IsomorphismPermGroup(W);
WP:=Image(WhomWP);
WPev:=EvenSubgroup(WP);

AhomWP:=GroupHomomorphismByFunction(A,WP,x->Image(WhomWP,Image(AhomW,x)));

EltsWP:=List(EltsG1,x->Image(AhomWP,x));
EltsWP:=Concatenation(
EltsWP,
Filtered(Elements(WP),x->not x in EltsWP));

###############################################################
StabilizerSubgroup:=function(n,k)
local G,U;

U:=Filtered(Vertices,a->not a in ResGens[n+1][k]);
G:=List( Vcomplement(EssentialSet(U)), x->GeneratorsOfGroup(WP)[x]);
if Length(G)>0 then return Group(G);fi;
return Group(());
end;
###############################################################

###############################################################
EvenStabGroup:=function(n,k)
local G,U,V,x,y;

U:=Filtered(Vertices,a->not a in ResGens[n+1][k]);
V:=Vcomplement(EssentialSet(U));
U:=Difference(V,ResGens[n+1][k]); 
if Length(ResGens[n+1][k])>0 then 
V:=[];
for x in ResGens[n+1][k] do
for y in ResGens[n+1][k] do
Add(V,GeneratorsOfGroup(WP)[x]*GeneratorsOfGroup(WP)[y]);
od;
od;
V:=Concatenation(V,List(U,x->GeneratorsOfGroup(WP)[x]));
return Group(V);fi;

return Group(());;

end;
###############################################################

###############################################################
# This describes how the group WP acts on the orientation.
Action:=function(n,k,g);
if n=0 then return 1; fi;
if
EltsWP[g] in WPev  then return 1;
#EltsWP[g] in EvenStabGroup(n,AbsInt(k))  then return 1;
else return -1; fi;  #THIS IS WRONG - VERY WRONG! FIX IT ASAP!
end;
###############################################################


return Objectify(HapNonFreeResolution,
	    rec(
	    dimension:=Dimension,
	    boundary:=Boundary,
	    homotopy:=fail,
	    elts:=EltsWP,
	    group:=WP,
	    stabilizer:=StabilizerSubgroup,
	    action:=Action,
	    properties:=
	    [["length",n],
	     ["characteristic",0],
	     ["type","resolution"],
	     ["reduced",true]]  ));
end);
#####################################################################

