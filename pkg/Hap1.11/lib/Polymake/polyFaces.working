#(C) Graham Ellis 2005-2006
#RT:=0;
#####################################################################
InstallGlobalFunction(PolytopalComplex,
function(arg)
local
	G,StartVector, Gev,
	PG,
	Action,
	VertexToVector, VVRecord,
	FaceToVertices,
	Hasse,
	p,x,n,i,
	Points,
	Dimension,
	Boundary,
	lngth,
	StabilizerSubgroup,
	StabilizerRecord,
        StabilizerBasisRecord,
        StabilizerBasis,
	BoundaryComponent,
	EltsG,
	PseudoBoundary,
	OrbitReps,
	StabSum,
	StabAction,
	VSGS;

G:=arg[1];
if IsPermGroup(G) then
G:=Image(PermToMatrixGroup(G));
fi;

StartVector:=arg[2];
PG:=PolytopalGenerators(G,StartVector);

if Length(arg)>2 then lngth:=arg[3]; else lngth:=Length(PG.hasseDiagram); fi;
Points:=[];
EltsG:=Elements(G);
VSGS:=VectorStabilizer(G,StartVector);

#####################################################################
Dimension:=function(k);
if k<0 then return 0; fi;
if k=0 then return 1; fi;
return Length(Hasse[k]);
end;
#####################################################################

#####################################################################
Action:=function(g,V) ;
return g*V;    
end;
#####################################################################


#########################CREATE POINTS###############################
for x in G do
Add(Points, Action(x,StartVector));
od;
Points:=SSortedList(Points);
#####################################################################


#####################################################################
VertexToVector:=function(v);
return Action(PG.generators[v+1],StartVector) - StartVector;
end;
#####################################################################


#####################################################################
FaceToVertices:=function(F)
local W,v,w,V;
#RT:=RT-Runtime();
V:=[];
W:=BaseOrthogonalSpaceMat(List(F,x->VertexToVector(x)));

if W=[] then
W:=[Points[1]*0];
fi;

for p in Points do

if IsZero((p - StartVector)*TransposedMat(W)) then Add(V,p); fi;
od;
#RT:=RT+Runtime();
return V;
end;
#####################################################################


Hasse:=[];
for x in [1..lngth] do
Append(Hasse,[List(PG.hasseDiagram[x],y->FaceToVertices(y))     ]);
od;


#####################################################################
OrbitReps:=function(L)  #L=Hasse[i]
local g,R,S, T,Reps,bool,count;

Reps:=[];
for S in L do
bool:=true;
count:=0;
for g in G do
count:=count+1;
T:=List(S,x->Action(g,x));
for R in Reps do
if Length(T)=Length(Intersection(T,R)) then
bool:=false; break; fi;
od;
if bool =false then break;fi;
if count=Order(G) then Add(Reps,SSortedList(S)); fi;
od;
od;

return Reps;
end;
#####################################################################


Hasse:=List(Hasse,x->OrbitReps(x));




StabilizerRecord:=List([1..lngth],i->[1..Dimension(i)]);
StabilizerBasisRecord:=List([1..lngth],i->[1..Dimension(i)]);


#####################################################################
StabilizerSubgroup:=function(kk,nn)
local S,T,xT,verts,StabGroup,x,k,n;

k:=AbsInt(kk);
n:=AbsInt(nn);

if k=0 then return VSGS; fi;

if not IsInt(StabilizerRecord[k][n]) then
return StabilizerRecord[k][n]; fi;

if k=Length(PG.hasseDiagram) then return G; fi;

StabGroup:=Group(One(G));
S:=Hasse[k][n];

T:=StructuralCopy(S);
for x in G do
if not x in StabGroup then
xT:=List(T,a->Action(x,a));
if Length(Intersection(xT,T))=Length(T) then 
StabGroup:=  GeneratorsOfGroup(StabGroup); ;
StabGroup:=Concatenation(StabGroup,[x]);
StabGroup:=Group(StabGroup);
fi;
fi;
od;


StabGroup:=ReduceGenerators(GeneratorsOfGroup(StabGroup),StabGroup);
if Length(StabGroup)=0 then StabGroup:=[Identity(G)]; fi;

StabilizerRecord[k][n]:=Group(StabGroup);
return StabilizerRecord[k][n];
end;
#####################################################################

#####################################################################
StabilizerBasis:=function(kk,nn)
local S,T,verts,bas,CG,x,k,n;

k:=AbsInt(kk);
n:=AbsInt(nn);

if k=0 then return []; fi;

if not IsInt(StabilizerBasisRecord[k][n]) then
return StabilizerBasisRecord[k][n]; fi;

if k=Length(PG.hasseDiagram) then return IdentityMat(Length(StartVector)); fi;

S:=Hasse[k][n];

CG:=Sum(S)/Length(S);
bas:=List(S,x->x-CG);
bas:=SemiEchelonMat(bas).vectors;

StabilizerBasisRecord[k][n]:=bas;
return StabilizerBasisRecord[k][n];
end;
#####################################################################



StabSum:=List([1..lngth],k->
Sum(List([1..Dimension(k-1)],j->Order(StabilizerSubgroup(k-1,j)))-1));



#####################################################################
BoundaryComponent:=function(k,m,n)  	#Let Fm be the m-th face in 
					#dimension k, and Fn the n-th
					#face in dimension k-1. Return 
					#a list [g1,...,gd] of the elements
					#gi in G such that gi.Fn lies in the
					#boundary of Fm. The list is maximal
					#with respect to the property that
					#gi*gj^-1 is not in the stabilizer
					#of Fn.
local 	Fm,Fn, Stab, Component, CompCpy,test,
	g, gFn;

Fm:=Hasse[k][m];
if k>1 then Fn:=Hasse[k-1][n];
else Fn:=[StartVector]; fi;
Stab:=StabilizerSubgroup(k-1,n);

Component:=[];
CompCpy:=[];

for g in G do
gFn:=SSortedList(List(Fn,x->Action(g,x)));
if Size(gFn) = Size(Intersection(gFn,Fm)) then
if not gFn in CompCpy then
Add(Component,g); 
Add(CompCpy,gFn);
fi;
fi;
od;

return SSortedList(Component);
end;
#####################################################################

PseudoBoundary:=List([1..lngth],i->[1..Dimension(i)]);

#####################################################################
Boundary:=function(k,mm)
local b,bb,x,n, bnd,signedbnd,bndbnd,tmp,m;

m:=AbsoluteValue(mm);

if not IsInt(PseudoBoundary[k][m]) then 
if mm>0 then return PseudoBoundary[k][m]; 
else return NegateWord(PseudoBoundary[k][m]);fi;
fi;

bnd:=[];
for n in [1..Dimension(k-1)] do
tmp:=BoundaryComponent(k,m,n);
tmp:=List(tmp, x->Position(EltsG,x));    ##########
tmp:=List(tmp, x->[n,x]);
Append(bnd,tmp);
od;


#PseudoBoundary[k][m]:=bnd;  ##WARNING: REMOVE THIS
#return Boundary(k,mm);


	######Inserting the signs#########
if k=1 then
bnd[1][1]:=-bnd[1][1];
fi;

if k>1  then
#################################
bndbnd:=[];
bnd:=SortedList(bnd);
signedbnd:=[bnd[1]]; 
x:=bnd[1]; 
b:=Boundary(k-1,x[1]);
b:=List(b,y->[y[1], Position(EltsG,EltsG[x[2]]*EltsG[y[2]])]);

b:=List(b,y->[StabAction(k-2,y[1],y[2])*y[1], CanonicalRightCosetElement(StabilizerSubgroup(k-2,AbsInt(y[1])), EltsG[y[2]]^-1)^-1]);

if StabAction(k-1,x[1],x[2])<0 then
b:=NegateWord(b);
fi;

Append(bndbnd,b);
RemoveSet(bnd,bnd[1]);

while Length(bnd)>0 do
x:=Random(bnd);
b:=Boundary(k-1,x[1]);
b:=List(b,y->[y[1], Position(EltsG,EltsG[x[2]]*EltsG[y[2]])]);

b:=List(b,y->[StabAction(k-2,y[1],y[2])*y[1], CanonicalRightCosetElement(StabilizerSubgroup(k-2,AbsInt(y[1])), EltsG[y[2]]^-1)^-1]);

if StabAction(k-1,x[1],x[2])<0 then
b:=NegateWord(b);
fi;


if Length(Intersection(b,bndbnd))>0 then
Append(signedbnd, [[-x[1],x[2]]]);
Append(bndbnd,NegateWord(b));
RemoveSet(bnd,x);
else

if Length(Intersection(NegateWord(b),bndbnd))>0 then
Append(signedbnd, [[x[1],x[2]]]);
Append(bndbnd,b);
RemoveSet(bnd,x);
fi;

fi;
od;

bnd:=signedbnd;
########################################
fi;
	######Signs inserted##############


PseudoBoundary[k][m]:=bnd;

if mm>0 then return PseudoBoundary[k][m];
else return NegateWord(PseudoBoundary[k][m]);fi;


end;

###############################################################
# This describes how the group G acts on the orientation.
StabAction:=function(n,k,h)
local bas, Gbas, mat,id,r,u,H; 

#n:=AbsInt(nn);

if n=0 then return 1; fi;

H:=StabilizerSubgroup(n,k);

id:=CanonicalRightCosetElement(H,Identity(H));
r:=CanonicalRightCosetElement(H,EltsG[h]^-1);
r:=id^-1*r;
u:=r*EltsG[h];



bas:=StabilizerBasis(n,k);
Gbas:=List(bas,V->Action(u,V));
mat:=List(Gbas, b->SolutionMat(bas,b));

return SignInt(Determinant(mat));
end;
###############################################################

for n in [1..lngth] do
for i in [1..Dimension(n)] do
Boundary(n,i);
od;od;




#EltsG:=List(EltsG,x->x^-1);

#####################################################################
return Objectify(HapNonFreeResolution,
	   rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=fail,
            elts:=EltsG,
            group:=G,
	    stabilizer:=StabilizerSubgroup,
            basis:=StabilizerBasis,
	    action:=StabAction,
	    hasse:=Hasse,
            properties:=
             [["type","resolution"],
              ["length",lngth],
              ["characteristic", 0] ]));

end);
#####################################################################

