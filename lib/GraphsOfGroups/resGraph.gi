#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionGraphOfGroups,
function(arg)
local

	D,N,L,
        PositionName,
	Vertices, Edges, EdgeGroups,
	ColouredVertices, ComplementaryEdges,
	resVertices, resEdges, resEdgeGroups,
	FpHoms, F,Frels,PhomF, FundGroup,
	inc,

	Dimension,
	Elts,
	PseudoBoundary,
	FillPseudoBoundary,
	Boundary,
	hom, quohom,
	pos,e,f,g,x,P,R,S,T;
	
####################################READ INPUT DATA##################
D:=arg[1];
N:=arg[2];
if Length(arg)>2 then L:=arg[3];
else L:=[];
fi;
########################################DATA READ####################

if not GraphOfGroupsTest(D) then 
Print("The list D does not represent a Graph of Groups \n");
return fail;
fi;

#####################################################################
PositionName:=function(L,x);
return PositionProperty(L,n->Name(n)=Name(x));
end;
#####################################################################

###########EXTRACT VERTICES ABD EDGES OF GRAPH OF GROUPS ############
Vertices:=[];
Edges:=[];

for x in D do
if IsGroup(x) then Append(Vertices,[x]); fi;
if IsList(x) then Append(Edges,[x]); fi;
od;

EdgeGroups:=List(Edges,e->Source(e[1]));
##########EXTRACTED##################################################

####################FIND MAXIMAL TREE IN GRAPH OF GROUPS#############
ColouredVertices:=[];
ComplementaryEdges:=[];

for e in Edges do
if Name(Range(e[1])) in ColouredVertices
and Name(Range(e[2])) in ColouredVertices
or Name(Range(e[1]))=Name(Range(e[2])) 
then
Append(ComplementaryEdges,[[Name(Range(e[1])),Name(Range(e[2]))]]);
fi;
Append(ColouredVertices,[Name(Range(e[1])),Name(Range(e[2]))]);
od;



#########################TREE FOUND##################################

##########CREATE GRAPH OF RESOLUTIONS################################
resVertices:=[];
for x in Vertices do
pos:=Position(List(L,r->r!.group),x);
if pos=fail then
if Order(x)<121 then           #This suggests I should improve
R:=ResolutionFiniteGroup(x,N); #ResolutionGenericGroup() !!!
else
R:=ResolutionGenericGroup(x,N);
fi;
else
R:=L[pos];
fi;
SetName(R!.group,Name(x));
Append(resVertices,[R]);
od;

resEdgeGroups:=[];
for x in EdgeGroups do
pos:=Position(List(L,r->r!.group),x);
if pos=fail then
#R:= ResolutionFiniteGroup(x,N-1);
R:= ResolutionGenericGroup(x,N-1);
else
R:=L[pos];
fi;
SetName(R!.group,Name(x));
Append(resEdgeGroups,[R]);
od;


resEdges:=[];
for x in Edges do
R:=resVertices[PositionName(Vertices,Range(x[1]))];
S:=resEdgeGroups[PositionName(EdgeGroups,Source(x[1]))];
T:=resVertices[PositionName(Vertices,Range(x[2]))];
e:=[EquivariantChainMap(S,R,x[1]),EquivariantChainMap(S,T,x[2])];
Append(resEdges,[e]);
od;
##########CREATED####################################################

##########CREATE HOMOMORPHISM TO FP GROUP############################
F:=FreeGroup(Length(ComplementaryEdges)+
		Sum(List(resVertices,R->R!.dimension(1))));


FpHoms:=[];
Frels:=[];
x:=0;

for R in resVertices do
P:=PresentationOfResolution(R);

PhomF:=GroupHomomorphismByImagesNC(P.freeGroup,F,
GeneratorsOfGroup(P.freeGroup),
GeneratorsOfGroup(F){[x+1..x+R!.dimension(1)]});

Append(Frels, List(P.relators,r->Image(PhomF,r)));

Append(FpHoms,
[GroupHomomorphismByImagesNC(R!.group,F,
List([1..R!.dimension(1)],i->R!.elts[R!.boundary(1,i)[2][2]]), #CHANGED March 2021
GeneratorsOfGroup(F){[x+1..x+R!.dimension(1)]}
)]);
x:=x+R!.dimension(1);
od;


pos:=Sum(List(resVertices,R->R!.dimension(1)));

for e in Edges do
if not [Name(Range(e[1])),Name(Range(e[2]))] in ComplementaryEdges 
and 
not [Name(Range(e[2])),Name(Range(e[1]))] in ComplementaryEdges
then


for x in GeneratorsOfGroup(Source(e[1])) do
Append(Frels,[
Image(FpHoms[PositionName(Vertices,Range(e[1]))],
Image(e[1],x))*
Image(FpHoms[PositionName(Vertices,Range(e[2]))],
Image(e[2],x))^-1
]);
od;

else

pos:=pos+1;

for x in GeneratorsOfGroup(Source(e[1])) do
Append(Frels,[
Image(FpHoms[PositionName(Vertices,Range(e[1]))],
Image(e[1],x))*GeneratorsOfGroup(F)[pos]*
Image(FpHoms[PositionName(Vertices,Range(e[2]))],
Image(e[2],x))^-1*GeneratorsOfGroup(F)[pos]^-1
]);
od;


fi;
od;

FundGroup:=F/Frels;  					     ###Added Aug 2013
quohom:= 
GroupHomomorphismByImages(F,FundGroup,                       ###
                            GeneratorsOfGroup(F),            ###
                            GeneratorsOfGroup(FundGroup));   ###
FundGroup:=F;                                                ###

x:=0;
FpHoms:=[];
Elts:=[];
for R in resVertices do
P:=PresentationOfResolution(R);
Append(FpHoms,
[GroupHomomorphismByImagesNC(R!.group,FundGroup,
List([1..R!.dimension(1)],i->R!.elts[R!.boundary(1,i)[2][2]]), ###CHANGED March 2021
GeneratorsOfGroup(FundGroup){[x+1..x+R!.dimension(1)]}
)]);
Append(Elts,Image(FpHoms[Length(FpHoms)],R!.elts));

x:=x+R!.dimension(1);
od;
##########DREATED####################################################

#####################################################################
Dimension:=function(k)
local dim,R;

if k<0 then return 0; fi;

dim:=0;
for R in resVertices do
dim:=dim+R!.dimension(k);
od;
for R in resEdgeGroups do
dim:=dim+R!.dimension(k-1);
od;

return dim;
end;
#####################################################################

#####################################################################
inc:=function(x,n)	#x is the vertex or edge group/resolution and
local pos,i,increment;	#n is the dimension in the TOTAL resolution

increment:=0;
pos:=PositionName(Vertices,x);
if not pos=fail then

for i in [1..pos-1] do
increment:=increment + resVertices[i]!.dimension(n);
od;
return increment;

else

for i in [1..Length(resVertices)] do
increment:=increment + resVertices[i]!.dimension(n);
od;

pos:=PositionName(EdgeGroups,x);

for i in [1..pos-1] do
increment:=increment + resEdgeGroups[i]!.dimension(n-1);
od;
return increment;

fi;

end;
#####################################################################

PseudoBoundary:=[];

#####################################################################
FillPseudoBoundary:=function()
local b,f,g,i,k,x,hom,hom1,hom2,R,dim,dim1,dim2,bnd;

for k in [1..N] do
PseudoBoundary[k]:=[];


	for R in resVertices do
	x:=PositionName(List(FpHoms,f->Source(f)),R!.group);
	hom:=FpHoms[x]; 
	dim:=inc(R!.group,k-1);

	for i in [1..R!.dimension(k)] do
	Append(PseudoBoundary[k],[
	List(R!.boundary(k,i),a->
	[
	SignInt(a[1])*(AbsoluteValue(a[1])+dim),
	Position(Elts,Image(hom,R!.elts[a[2]]))
	])]);
	od;
	od;
	
	#######################
	
	for e in Edges do
	hom:=e[1];
	x:=PositionName(EdgeGroups,Source(e[1]));
	R:=Source(resEdges[x][1]);
	dim:=inc(R!.group,k-1);
	S:=Target(resEdges[x][1]);
	T:=Target(resEdges[x][2]);
	f:=resEdges[x][1];
	g:=resEdges[x][2];

	dim1:=inc(S!.group,k-1);
	x:=PositionName(List(FpHoms,f->Source(f)),S!.group);
	hom1:=FpHoms[x];

	dim2:=inc(T!.group,k-1);
	x:=PositionName(List(FpHoms,f->Source(f)),T!.group);
        hom2:=FpHoms[x];


	for i in [1..R!.dimension(k-1)] do
	bnd:=[];
	
	b:=List( f!.mapping([[i,1]],k-1),a->
	[
	SignInt(a[1])*(AbsoluteValue(a[1])+dim1),
	Position(Elts,Image(hom1,S!.elts[a[2]]))
	]);
if IsOddInt(k) then
Append(bnd,b);
else
Append(bnd,NegateWord(b));
fi;

	 b:=List( g!.mapping([[i,1]],k-1),a->
	 [
         SignInt(a[1])*(AbsoluteValue(a[1])+dim2),
#        Position(Elts,Image(hom1,T!.elts[a[2]]))  #Typo found Jan 2012
         Position(Elts,Image(hom2,T!.elts[a[2]]))
         ]);
if IsOddInt(k) then
Append(bnd,NegateWord(b));
else
Append(bnd,b);
fi;
	b:=List(R!.boundary(k-1,i),a->
	[
	SignInt(a[1])*(AbsoluteValue(a[1])+dim),
        Position(Elts,Image(hom1,Image(hom,R!.elts[a[2]])))
       ]);
	
	Append(bnd,b);

	Append(PseudoBoundary[k],[bnd]);
	od;
	od;

od;
end;
#####################################################################

FillPseudoBoundary();

#####################################################################
Boundary:=function(k,i);
if k<1 then return 0; fi;
if i>0 then return PseudoBoundary[k][i];fi;
return NegateWord(PseudoBoundary[k][-i]);
end;
####################################################################

Apply(Elts,x->Image(quohom,x));

return Objectify(HapResolution,
                rec(
                dimension:=Dimension,
                boundary:=Boundary,
                homotopy:=fail,
                elts:=Elts,
                group:=Range(quohom),
                properties:=
                [["length",N],
                 ["reduced",true],
                 ["type","resolution"],
                 ["characteristic",0]  ]));


end);
#####################################################################
