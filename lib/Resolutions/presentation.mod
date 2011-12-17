#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(PresentationOfResolution,
function(R)
local  
		Dimension,
		Boundary,
		Elts, Mult, Inv,
		F, Frels, Fgens,
		gens, rels,
		FirstBoundaryHomomorphism,
		Boundary2Relator, Relator2Word,
		start,
		b, r, x;

if not (IsHapResolution(R) or IsHapNonFreeResolution(R)) then
Print("This function must be applied to a resolution. \n");
return fail;
fi;

if not EvaluateProperty(R,"reduced")=true then
if R!.dimension(0)>1 then
Print("This function must be applied to a REDUCED resolution. \n");
return fail; fi;
fi;

if not EvaluateProperty(R,"characteristic")=0 then
Print("This function only works in characteristic 0. \n");
return fail;
fi;


Dimension:=R!.dimension;
Elts:=R!.elts;

#####################################################################
Mult:=function(g,h);
return Position(Elts,Elts[g]*Elts[h]);
end;
#####################################################################

#####################################################################
Inv:=function(g);
return Position(Elts,Elts[g]^-1);
end;
#####################################################################

#####################################################################
Boundary:=function(n,k)
local b;
b:= R!.boundary(n,k);
return List(b,x->[x[1],Mult(Inv(b[1][2]),x[2])]);

end;
#####################################################################

F:=FreeGroup(Dimension(1));
Fgens:=GeneratorsOfGroup(F);
Frels:=[];
gens:=[];

start:=List([1..Dimension(2)],x->List(Boundary(2,x),y->y[2]));
start:=SortedList(Intersection(start))[1];

#####################################################################
FirstBoundaryHomomorphism:=function(x)
local r;
r:=Boundary(1,x[1]);
r:=List(r,y->Mult(x[2],y[2]));
if x[1]>0 then return r;
else return Reversed(r); fi;
end;
#####################################################################

#####################################################################
Boundary2Relator:=function(b)
local c, rel, w;

b:=SortedList(AlgebraicReduction(b));
rel:=[start];

while Length(b)>0 do
	for x in b do
	w:=FirstBoundaryHomomorphism(x);
	if w[1]= rel[Length(rel)] then
	Append(rel, [w[2]]); RemoveSet(b,x); break; 
	else
	   if w[2]= rel[Length(rel)] then
	   Append(rel, [w[1]]); RemoveSet(b,x); break;
	   fi;
	fi;
	od;
od;

return rel;
end;
#####################################################################

for r in [1..Dimension(2)] do
Append(Frels,[Boundary2Relator(Boundary(2,r))]);
od;



for r in Frels do
if (not Inv(r[2]) in gens) then AddSet(gens,r[2]);fi;
if (not Inv(r[Length(r)-1]) in gens) then AddSet(gens,r[Length(r)-1]);fi;
od;

#####################################################################
Relator2Word:=function(r)
local w,v,g,h;
w:=Identity(F);

for v in [2..Length(r)] do
	for g in gens do
	if Mult(r[v-1],g)=r[v] then h:=Position(gens,g); break; fi;
	if Mult(r[v-1],Inv(g))=r[v] then h:=-Position(gens,g); break; fi;
	od;
w:=w*Fgens[AbsoluteValue(h)]^SignInt(h);
od;

return w;
end;
#####################################################################

rels:= List(Frels,g->Relator2Word(g));

return rec( 
		freeGroup:=F,
		relators:=rels );
end);
#####################################################################

