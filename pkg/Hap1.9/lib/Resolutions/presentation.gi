#(C) Graham Ellis, 2005-2006

#####################################################################
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
Boundary:=R!.boundary;
Elts:=R!.elts;
F:=FreeGroup(Dimension(1));
Fgens:=GeneratorsOfGroup(F);
Frels:=[];
gens:=[];

#####################################################################
Mult:=function(g,h)
local pos;
pos:=Position(Elts,Elts[g]*Elts[h]);
if pos=fail then Add(Elts,Elts[g]*Elts[h]); 
pos:=Length(Elts);
fi;
return pos;
end;
#####################################################################

#####################################################################
Inv:=function(g)
local pos;
pos:= Position(Elts,Elts[g]^-1);
if pos=fail then Add(Elts,Elts[g]^-1);
pos:=Length(Elts);
fi;
return pos;

end;
#####################################################################


start:=List([1..Dimension(2)],x->List(Boundary(2,x),y->y[2]));

if Length(Intersection(start))=0 then
###############################################
Boundary:=function(n,k)
local w;
w:=Inv(R!.boundary(n,k)[1][2]) ;
return
List(R!.boundary(n,k),x->[x[1],Mult(w,x[2])]);
end;
##############################################

start:=List([1..Dimension(2)],x->List(Boundary(2,x),y->y[2]));
fi;

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

b:=SortedList(AlgebraicReduction(b));  #I'm assuming a regular CW-space.
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
		relators:=rels,
                gens:=gens );
end);
#####################################################################
#####################################################################



#####################################################################
#####################################################################
InstallGlobalFunction(ResolutionToResolutionOfFpGroup,
function(R)
local
         P,FF,epi,F,G,FFhomF,x,Iso,OldElts,gensFF,tmp;

OldElts:=StructuralCopy(R!.elts);
P:=PresentationOfResolution(R);
F:=P.freeGroup/P.relators;
G:=Group(R!.elts{P.gens});
epi:=EpimorphismFromFreeGroup(G);;
FF:=Source(epi);;
gensFF:=GeneratorsOfGroup(FF);
tmp:=List(gensFF,x->x^-1);
gensFF:=Concatenation(gensFF,tmp);

FFhomF:=GroupHomomorphismByImagesNC(FF,F,
                   GeneratorsOfGroup(FF), GeneratorsOfGroup(F));



####################################
Iso:=function(L)
local t,x,pos,cnt,w,nodup,newelts,Tips,NewTips,elts;

cnt:=0;
nodup:=[];
newelts:=[];
Tips:=[Identity(FF)];
NewTips:=[];
elts:=[];

while Length(Tips)>0 do
for t in Tips do 
for x in gensFF do
w:=Image(epi,t*x);
if not w in elts then Add(elts,w); Add(NewTips,t*x);
######
pos:=Position(L,Image(epi,t*x));
if IsInt(pos) then 
  if not IsBound(nodup[pos]) then 
  newelts[pos]:=Image(FFhomF,t*x); cnt:=cnt+1; 
  nodup[pos]:=pos;
  fi;
fi;
if cnt=Length(L) then break; fi;

######
fi;
od;
if cnt=Length(L) then break; fi;
od;
if cnt=Length(L) then break; fi;
Tips:=NewTips;
NewTips:=[];
od;

return newelts;
end;
####################################

R!.elts:=Iso(OldElts);
R!.group:=F;
R!.isomorphism:=Iso;

end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(PresentationOfResolution_alt,
function(R)
local
                Dimension,
                Boundary,
                Elts, Mult, Inv,
                F, Frels, Fgens,
                gens, rels,
                GENS,
                FirstBoundaryHomomorphism,
                Boundary2Relator, Relator2Word,
                i,b, r, x;

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
Boundary:=R!.boundary;
Elts:=R!.elts;
F:=FreeGroup(Dimension(1));
Fgens:=GeneratorsOfGroup(F);
Frels:=[];
gens:=[];
GENS:=[];

#####################################################################
Mult:=function(g,h)
local pos;
pos:=Position(Elts,Elts[g]*Elts[h]);
if pos=fail then Add(Elts,Elts[g]*Elts[h]);
pos:=Length(Elts);
fi;
return pos;
end;
#####################################################################

#####################################################################
Inv:=function(g)
local pos;
pos:= Position(Elts,Elts[g]^-1);
if pos=fail then Add(Elts,Elts[g]^-1);
pos:=Length(Elts);
fi;
return pos;

end;
#####################################################################

for i in [1..Dimension(1)] do
b:=Boundary(1,i);
x:=R!.elts[b[1][2]] * R!.elts[b[2][2]]^-1;
if b[1][1]<0 then x:=x^-1; fi;
Add(GENS,x);
GENS:=Reversed(GENS);
od;

gens:=List(GENS,x->Position(R!.elts,x));
return gens;

#####################################################################
FirstBoundaryHomomorphism:=function(i)
local r ;
r:=Boundary(1,i);
r:=List(r,w->[w[1],w[2],Mult(x[2],w[2])]);
if x[1]>0 then return r;
else return Reversed(r); fi;
end;
#####################################################################

#####################################################################
Boundary2Relator:=function(b)
local c, rel, w;

end;
#####################################################################

for r in [1..Dimension(2)] do
Append(Frels,[Boundary2Relator(Boundary(2,r))]);
od;


for r in Frels do
if (not Inv(r[2]) in gens) then AddSet(gens,r[2]);fi;
if (not Inv(r[Length(r)-1]) in gens) then AddSet(gens,r[Length(r)-1]);fi;
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
                relators:=rels,
                gens:=gens );
end);
#####################################################################
#####################################################################









