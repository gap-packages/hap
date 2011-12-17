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
		b, r, x;

Dimension:=R.dimension;
Boundary:=R.boundary;
Elts:=R.elts;
F:=FreeGroup(Dimension(1));
Fgens:=GeneratorsOfGroup(F);
Frels:=[];
gens:=[];

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

b:=SortedList(b);
rel:=[1];  

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

