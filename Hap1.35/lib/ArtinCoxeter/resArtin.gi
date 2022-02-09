#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionArtinGroup,
function(D,K)
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
	m, n, i, S, SD,R;

Vertices:=CoxeterDiagramVertices(D);
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
SD:=CoxeterSubDiagram(D,S);
if CoxeterDiagramIsSpherical(SD) then AddSet(ResGens[n+1],S); fi;
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
if not y in EltsG then Add(EltsG,y); 
y:=Image(GhomG1,y); 
Add(EltsG1,y);
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
a:=Difference(S,T)[1];
Append(bndry,[  [T,BoundaryCoeff(S,T),Position(S,a)]  ]);
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


R:=         Objectify(HapResolution,
	    rec(
	    dimension:=Dimension,
	    boundary:=Boundary,
	    homotopy:=fail,
	    elts:=EltsG1,
	    group:=G1,
	    resGens:=ResGens,
	    properties:=
	    [["length",n],
	     ["characteristic",0],
	     ["type","resolution"],
	     ["reduced",true]]  ));
for n in [1..Length(R)] do
for i in [1..R!.dimension(n)] do
R!.boundary(n,i);
od;od;

return R;

end);
#####################################################################

