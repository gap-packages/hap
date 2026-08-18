#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(IsAspherical,
function(arg)
local
        F,Frels,G,
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
	TestFeasible,
        L,r,
	R,Redges,V,ii,bol, sm, i,x,row,n,m,M;

if Length(arg)=2 then F:=arg[1]; Frels:=arg[2]; fi;
if Length(arg)=1 then
G:=arg[1];
F:=FreeGroupOfFpGroup(G); Frels:=RelatorsOfFpGroup(G); fi;
Frels:=List(Frels,r->CyclicallyReducedWord(r));
Vertices:=[];
CollEdges:=[];
Edges:=[];

#####################################################################
RelatorToVertices:=function(R)
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
InEqualities:=NullspaceModQ(BoundaryMat,2);  #THIS IS CLUMSY!!

InEqualities:=Filtered(InEqualities,r->not IsZero(r));  


        ###MODIFIED AUG 2016
#SortBy(InEqualities,Sum);
#L:=InEqualities;
#for ii in [1..Minimum(1000,Length(InEqualities))] do
##ANF WHY 1000??
#sm:=Sum(InEqualities[ii]);
#L:=Filtered(L,r->r*InEqualities[ii]<sm or r=InEqualities[ii]);
#od;
#InEqualities:=L;
        ###MODIFICATION STOPS HERE


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
TestFeasible:=function()
local poly, n, equations, inequalities;

n:=Length(InEqualities[1]);

equations:=List([1..Length(EqualitiesMat)],
                i->Concatenation([-EqualitiesVec[i]],EqualitiesMat[i]));

# every generator is asked for at least 2, and for at least 0 individually
inequalities:=Concatenation(
        List(InEqualities,row->Concatenation([-2],row)),
        List(IdentityMat(n),row->Concatenation([0],row)));

poly:=CreatePolymakeObject();
if Length(equations)>0 then
    # polymaking 0.8 has no method taking a matrix for a section it does not
    # know by name, so EQUATIONS has to go through a string; the rest of this
    # uses the higher level methods, which both versions provide
    AppendToPolymakeObject(poly,
            ConvertMatrixToPolymakeString("EQUATIONS",equations));
fi;
AppendInequalitiesToPolymakeObject(poly,inequalities);

return Polymake(poly,"FEASIBLE");
end;
#####################################################################

x:= TestFeasible();

if x=true then
Print("Presentation is aspherical.\n\n"); return true;
else
Print("Test inconclusive.\n\n");
return fail; fi;

end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallMethod(StarGraph, "For FPGroups", [IsFpGroup], StarGraphAttr);
InstallMethod(StarGraphAttr,
"For FpGroups",
[IsFpGroup],
function(G)
local
        F,Frels,Vertices,
        CollEdges,
        Edges,
        RelatorToVertices,
        RelatorToEdges,
        A,R,i,j,x;


F:=FreeGroupOfFpGroup(G);;
Frels:=RelatorsOfFpGroup(G);;
Frels:=List(Frels,r->CyclicallyReducedWord(r));
Vertices:=[];
CollEdges:=[];
Edges:=[];

#####################################################################
RelatorToVertices:=function(R)
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

A:=NullMat(Length(Vertices),Length(Vertices));
for x in CollEdges do
i:=Position(Vertices,x[1][1]);
j:=Position(Vertices,x[1][2]);
A[i][j]:=Minimum(2,A[i][j]+x[2]);A[j][i]:=Minimum(2,A[j][i]+x[2]);
od;
return IncidenceMatrixToGraph(A);

end);
#####################################################################
#####################################################################
