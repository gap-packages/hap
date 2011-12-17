#(C) Graham Ellis 2005-2006

#####################################################################
InstallGlobalFunction(PolytopalComplex,
function(G,StartVector)
local
	PG,
	GG,
	Action,
	VertexToVector,
	FaceToVertices,
	Hasse,
	p,x,
	CreatPoints,
	Points,
	Dimension,
	Boundary,
	lngth,
	StabilizerSubgroup,
	VectorToGroupElt;

PG:=PolytopalGenerators(G,StartVector);
lngth:=Length(PG.hasseDiagram);
Points:=[];
GG:=Filtered(Elements(G),x->not x=Identity(G));

#####################################################################
Dimension:=function(k);
if k<0 then return 0; fi;
if k=0 then return 1; fi;
return Length(PG.hasseDiagram[k]);
end;
#####################################################################

#########################CREATE POINTS###############################
CreatPoints:=function()
local x,w,D,i,v;

D:=Length(StartVector);

if IsPermGroup(G) then

for x in G do
w:=[];
	for i in [1..D] do
        Append(w,[StartVector[i^x]]);
        od;
Append(Points, [w]);
od;

else

for x in G do
w:=x*StartVector;
Append(Points, [w]);
od;
fi;

return Points;
end;
#####################################################################

CreatPoints();
										
#####################################################################
Action:=function(g,V)
local i,gV;
gV:=[];

for i in [1..Length(V)] do
gV[i]:=V[i^g];
od;

return gV;
end;
#####################################################################

#####################################################################
VertexToVector:=function(v);
return Action(PG.generators[v+1],StartVector) - StartVector;
end;
#####################################################################

#####################################################################
VectorToGroupElt:=function(v)  #This is the clumsiest possible implementation!
local g;
for g in G do
if Action(g,StartVector)=v then return g; fi;
od;
end;
#####################################################################

#####################################################################
FaceToVertices:=function(F)
local W,v,w,bool,V;
V:=[];
W:=BaseOrthogonalSpaceMat(List(F,x->VertexToVector(x)));

for p in Points do
bool:=true;
for w in W do
if not (p - StartVector)*w=0 then bool :=false; break; fi;
od;
if bool then Append(V,[p]); fi;
od;

return V;
end;
#####################################################################

Hasse:=[];
for x in PG.hasseDiagram do
Append(Hasse,[List(x,y->FaceToVertices(y))     ]);
od;

#####################################################################
StabilizerSubgroup:=function(k,n)
local S,T,verts,StabGroup,x;
if k=0 then return VectorStabilizer(G,StartVector); fi;
if k=lngth then return G; fi;

StabGroup:=[];
S:=Hasse[k][n];
T:=List(S, i->VectorToGroupElt(i));
for x in T do
if Length(Intersection(x*T,T))=Length(T) then Append(StabGroup,[x]); fi;
od;

StabGroup:=Concatenation(StabGroup,
	GeneratorsOfGroup(VectorStabilizer(G,StartVector)));
StabGroup:=ReduceGenerators(StabGroup,Group(StabGroup));
return Group(StabGroup);
end;
#####################################################################

return rec(
            dimension:=Dimension,
#            boundary:=Boundary,
            homotopy:=fail,
            elts:=Elements(G),
            group:=G,
	    stabilizer:=StabilizerSubgroup,
            properties:=
             [["type","nonFreeResolution"],
              ["length",lngth],
              ["characteristic", 0] ]);

end);
