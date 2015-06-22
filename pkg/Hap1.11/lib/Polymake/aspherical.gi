#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(IsAspherical,
function(F,Frels)
local
	Vertices,
	CollEdges,
	Edges,
	RelatorToVertices,
	RelatorToEdges,
	BoundaryMat,
	BoundaryKer,
	Vector,
	EqualitiesMat,
	EqualitiesVec,
	InEqualities,
	InEqSizes,
	Polymake,
        tmpdir, tmpin, tmpodir, tmpout,
	R,Redges,V, i,x,row,n,m,M;

tmpdir := DirectoryTemporary();;
tmpin:=Filename( tmpdir , "tmpIn.log" );
tmpodir := DirectoryTemporary();;
tmpout:=Filename( tmpdir , "tmpOut.log" );


Vertices:=[];
CollEdges:=[];
Edges:=[];

#####################################################################
RelatorToVertices:=function(R);
return SSortedList(LetterRepAssocWord(R));
end;
#####################################################################

for R in Frels do
UniteSet(Vertices, RelatorToVertices(R));
UniteSet(Vertices, -1*RelatorToVertices(R));
od;

#####################################################################
RelatorToEdges:=function(R)
local Rvertices, Redges, u,v,i;

Rvertices:=LetterRepAssocWord(R);
Redges:=[];

for i in [1..(Length(Rvertices)-1)] do
u:=Rvertices[i];
v:=-Rvertices[i+1];
Append(Redges,[  SortedList([u,v])  ]);
od;

u:=Rvertices[Length(Rvertices)];
v:=-Rvertices[1];
Append(Redges,[  SortedList([u,v])  ]);

return Redges;
end;
####################################################################

for R in Frels do
Append(CollEdges, RelatorToEdges(R));
od;

CollEdges:=Collected(CollEdges);
Edges:=List(CollEdges,x->x[1]);

BoundaryMat:=[];
for x in Edges do
row:=List([1..Length(Vertices)],i->0);
row[Position(Vertices,x[1])]:=1;
row[Position(Vertices,x[2])]:=1;
Append(BoundaryMat,[row]);
od;

#EQUATIONS
EqualitiesMat:=[];
EqualitiesVec:=[];

for R in Frels do
M:=[];
	Redges:=RelatorToEdges(R);
	n:=Length(LetterRepAssocWord(R));
	row:=List([1..Length(Edges)+1],i->0);

	for x in Redges do
	m:=Position(Edges,x);
	row[m]:=row[m]+1;
	od;
	row[Length(Edges)+1]:=n-2;

for x in [1..Length(Edges)] do
Append(M,[row[x]]);
od;
Append(EqualitiesMat,[M]);
Append(EqualitiesVec,[n-2]);
od;

#Vector:=SolutionMat(TransposedMat(EqualitiesMat),EqualitiesVec);

#INEQUALITIES
InEqualities:=NullspaceModQ(BoundaryMat,2);
V:=List([1..Length(InEqualities[1])],i->0);
RemoveSet(InEqualities,V);


InEqSizes:=List(InEqualities,x->Sum(x));
for i in [1..2*LogInt(Length(InEqualities)+1,2)] do
m:=Position(InEqSizes,Minimum(InEqSizes));
InEqSizes[m]:=Maximum(InEqSizes)+1;
V:=InEqualities[m];
if not V=0 then
for x in [1..Length(InEqualities)] do
	if V*InEqualities[x]=Sum(V) and (not x=m) then
	InEqualities[x]:=0; fi;
od;
fi;
od;
InEqualities:=Filtered(InEqualities,x->not x=0);

for x in CollEdges do
if x[2]>1 then
V:=List([1..Length(Edges)],i->0);
V[Position(Edges,x[1])]:=2;

Append(InEqualities,[V]);
fi;
od;

#####################################################################
Polymake:=function()
local V, i, x, input;

AppendTo(tmpin,"EQUATIONS","\n");

for i in [1..Length(EqualitiesMat)] do
AppendTo(tmpin,-EqualitiesVec[i]," ");
for x in [1..Length(EqualitiesMat[1])] do
AppendTo(tmpin,EqualitiesMat[i][x]," ");
od;
AppendTo(tmpin,"\n");
od;

AppendTo(tmpin,"\n","INEQUALITIES","\n");

for i in [1..Length(InEqualities)] do
AppendTo(tmpin,-2," ");
for x in [1..Length(InEqualities[1])] do
AppendTo(tmpin,InEqualities[i][x]," ");
od;
AppendTo(tmpin,"\n");
od;

for i in [1..Length(InEqualities[1])] do
AppendTo(tmpin,0," ");
for x in [1..Length(InEqualities[1])] do
if i=x then AppendTo(tmpin,1," ");
else  AppendTo(tmpin,0," "); fi;
od;
AppendTo(tmpin,"\n");
od;

Exec(Concatenation(POLYMAKE_PATH, tmpin, " FEASIBLE > ",tmpout));
Exec(Concatenation("rm ",tmpin));
input := InputTextFile(tmpout);
x:=ReadLine(input);
x:=Int(Chomp(ReadLine(input)));
Exec(Concatenation("rm ",tmpout));

return x;
end;
#####################################################################

x:= Polymake();
if x=1 then
Print("Presentation is aspherical.\n\n"); return true;
else
Print("Test inconclusive.\n\n");
return fail; fi;

end);
#####################################################################
#####################################################################

