#(C) Graham Ellis, 2005-2006

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
		FirstBoundaryHomomorphism,
		Boundary2Relator, Relator2Word,
		start,
		b, r, x;

if not 
(IsHapResolution(R) or IsHapNonFreeResolution(R)
or IsHapEquivariantCWComplex(R)) then
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
InstallGlobalFunction(PresentationOfResolution,
function(R)
local
                Dimension,
                Boundary,
                Elts, Mult, Inv,
                HGens,HRels,HRels1,
                Tree,Verts, NTree, NNTree, VertElts,
                index, cnt, gens, Gens, F,
                src,trg,
                alpha, len,
                B,i,a, b,bb,x,y,g,lst,bool;

if not (IsHapResolution(R) or IsHapEquivariantCWComplex(R)) then
Print("This function must be applied to a resolution. \n");
return fail;
fi;

if IsHapResolution(R) and not EvaluateProperty(R,"characteristic")=0 then
Print("This function only works in characteristic 0. \n");
return fail;
fi;

Dimension:=R!.dimension;
Boundary:=R!.boundary;
Elts:=R!.elts;

######Let's make sure the first two boundaries of R are computed.
i:=List([1..Dimension(1)],i->Boundary(1,i));
i:=List([1..Dimension(2)],i->Boundary(2,i));
i:=0;

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

###############################
HGens:=List([1..Dimension(1)],i->Boundary(1,i));
###############################

###############################
Tree:=[];
Verts:=[1];
while Length(Verts)<Dimension(0) do
for x in HGens do
if AbsInt(x[1][1]) in Verts and not AbsInt(x[2][1]) in Verts then
Add(Verts, AbsInt(x[2][1])); 
Add(Tree,[x[1],x[2],Position(HGens,x)]);
fi;
if AbsInt(x[2][1]) in Verts and not AbsInt(x[1][1]) in Verts then
Add(Verts, AbsInt(x[1][1])); 
Add(Tree,[x[1],x[2],Position(HGens,x)]); 
fi;
od;
od;
NTree:=1*Tree;
Tree:=List(Tree,x->x[3]);
Verts:=SSortedList(Verts);
###############################



NTree:=SSortedList(NTree);
NNTree:=[];

if Length(NTree)>0 then
a:=NTree[1];
RemoveSet(NTree,a);
AddSet(NNTree,a);

while Length(NTree)>0 do
a:=1*Random(NTree);

for x in NNTree do

if AbsInt(x[1][1])=AbsInt(a[1][1]) then
RemoveSet(NTree,a);
g:=a[1][2]*1;
a[1][2]:=x[1][2]*1; a[2][2]:= Mult( x[1][2], Mult(Inv(g),a[2][2]))*1;
AddSet(NNTree,a);

break; fi;

if AbsInt(x[1][1])=AbsInt(a[2][1]) then
RemoveSet(NTree,a);
g:=a[2][2]*1;
a[2][2]:=x[1][2]*1; a[1][2]:= Mult( x[1][2], Mult(Inv(g),a[1][2]))*1;
AddSet(NNTree,a);

break; fi;

if AbsInt(x[2][1])=AbsInt(a[1][1]) then
RemoveSet(NTree,a);
g:=a[1][2]*1;
a[1][2]:=x[2][2]*1; a[2][2]:= Mult( x[2][2], Mult(Inv(g),a[2][2]))*1;
AddSet(NNTree,a);

break; fi;

if AbsInt(x[2][1])=AbsInt(a[2][1]) then
RemoveSet(NTree,a);
g:=a[2][2]*1;
a[2][2]:=x[2][2]*1; a[1][2]:= Mult( x[2][2], Mult(Inv(g),a[1][2]))*1;
AddSet(NNTree,a);

break; fi;

od;

od;
fi;

NTree:=NNTree;




VertElts:=[];
if Length(NTree)=0 then VertElts[1]:=1;fi;
for x in NTree do
VertElts[AbsInt(x[1][1])]:=x[1][2];
VertElts[AbsInt(x[2][1])]:=x[2][2];
od;



##########################
#HRels will be a list of relators given as ordered edges of a polygon.
HRels:=[];
for i in [1..Dimension(2)] do
b:=[];

for x in StructuralCopy(Boundary(2,i)) do
y:=StructuralCopy(HGens[AbsInt(x[1])]);
if SignInt(x[1])<0 then y[1][1]:=-y[1][1]; y[2][1]:=-y[2][1]; fi;
y[1][2]:=Mult(x[2],y[1][2]);
y[2][2]:=Mult(x[2],y[2][2]);
Add(b,[y[1],y[2],x[1]]);
od;

B:=1*b;
bb:=[b[1]]; b:=Difference(b,bb);
src:=bb[1][1];
trg:=bb[1][2];
while Length(b)>0 do
   lst:=StructuralCopy(trg);
   lst:=StructuralCopy([-lst[1],lst[2]]);
      for x in b do
      if x[1]=lst then src:=x[1]; trg:=x[2];  break; fi; 
      if x[2]=lst then src:=x[2]; trg:=x[1];  break; fi;
      od;
   Add(bb,x); b:=Difference(b,bb);
od;
Add(HRels,bb);
od;
###################################


cnt:=0;
index:=[];
Gens:=[];
for i in [1..Length(HGens)] do
if not i in Tree then cnt:=cnt+1; index[i]:=cnt; Add(Gens,HGens[i]);fi;
od;

#Gens contains the information needed to return gens
Gens:=List(Gens,x->
[   [x[1][1],VertElts[AbsInt(x[1][1])]]  ,  
    [x[2][1],Mult( Mult(VertElts[AbsInt(x[1][1])],Inv(x[1][2])) ,x[2][2])]  ]   );
Gens:=List(Gens,x-> Mult(x[2][2], Inv(VertElts[AbsInt(x[2][1])])));



HRels1:=[];
for x in HRels do
a:=Filtered(x,a->not AbsInt(a[3]) in Tree);
a:=List(a,i->SignInt(i[3])*index[AbsInt(i[3])]);
Add(HRels1,a);
od;

alpha:=["x","y","z","w","v","u","t","s","r","q","p"];
len:=Length(HGens)-Length(Tree);
if len<12 then
F:=FreeGroup(alpha{[1..len]});
else
F:=FreeGroup(Length(HGens)-Length(Tree));
fi;
gens:=GeneratorsOfGroup(F);
HRels:=[];

for x in HRels1 do
a:=Identity(F);
for i in x do
a:=a*gens[AbsInt(i)]^SignInt(i);
od;
Add(HRels,a);
od;

return rec(
                freeGroup:=F,
                relators:=HRels,
                gens:=Gens
          );
end);
#####################################################################
#####################################################################









