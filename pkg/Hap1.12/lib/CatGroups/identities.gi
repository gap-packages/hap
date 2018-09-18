#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(IdentityAmongRelatorsDisplay,
function(arg)
local  
		R,idnum,
		Dimension,
		Boundary,
		Elts, Mult, Inv,
		Frels, rels,
		Fgens,gens,
		FirstBoundaryHomomorphism,
		Boundary2Relator, 
		start,
		ActWord,
		idnt,
		Edges,
		Color,
		COLOURS,
  		tmpDir,tmpInlog,tmpIn2log,tmpIngif,
                AppendTo, PrintTo,
		b, r, x,i,X;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

tmpDir:=DirectoryTemporary();
tmpInlog:=Filename(tmpDir,"tmpInlog");
tmpIn2log:=Filename(tmpDir,"tmpIn2log");
tmpIngif:=Filename(tmpDir,"tmpIngif");

R:=arg[1];
idnum:=arg[2];

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
Frels:=[];
start:=List([1..Dimension(2)],x->List(Boundary(2,x),y->y[2]));
start:=SortedList(Intersection(start))[1];
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

COLOURS:=["blue","red","green","yellow","brown","black"];

#####################################################################
Color:=function(r)
local g,h;

        for g in gens do
        if Mult(r[1],g)=r[2] then h:=Position(gens,g); break; fi;
        if Mult(r[1],Inv(g))=r[2] then h:=Position(gens,g); break; fi;
        od;

h:=(h-1) mod 6;
h:=h+1;
return COLOURS[h];
end;
#####################################################################


#####################################################################

idnt:=[];

for b in Boundary(3,idnum) do
x:=Boundary(2,AbsInt(b[1]));
x:=Boundary2Relator(x);
if SignInt(b[1])=-1 then x:=Reversed(x); fi;
x:=List(x,t->Mult(b[2],t));
Add(idnt,x);
od;

Edges:=[];
for b in idnt do
for x in [1..Length(b)-1] do
if b[x]<b[x+1] then r:=[b[x],b[x+1]];
else r:=[b[x+1],b[x]]; fi;
#AddSet(Edges,r);
Add(Edges,r);
od;
od;

Edges:=SortedList(Edges);
Edges:=List([1..Length(Edges)/2],i->Edges[2*i]);


################ WRITE TO TMPIN.LOG #################################

AppendTo(tmpInlog," graph G { \n size=\"4,4\" \n subgraph cluster0 {\n node [shape=ellipse, width=.2,height=.2,fixedsize=true,style=filled, color=gray35,label=\"\"] \n edge [style=\"setlinewidth(2)\"] \n");

for x in Edges do
AppendTo(tmpInlog,x[1]," -- ", x[2], "[color=",Color(x), "] \n");
od;

####
X:=FreeGroup(Length(gens));
X:=GeneratorsOfGroup(X);
AppendTo(tmpInlog," }\n subgraph cluster1 {\n  node [shape=box, width=2,height=1,fixedsize=true,style=filled, color=white,fillcolor=white] \n ");

if Maximum(List(X,x->Length(String(x))))<20 then
for i in [1..Length(X)] do
AppendTo(tmpInlog,-i,"  [fontcolor= ",COLOURS[i],",label=\"", X[i],"\" ] \n");
od;
fi;

####

AppendTo(tmpInlog,"}\n }\n");

############### WRITTEN ############################################
Exec(Concatenation(NEATO_PATH,"-Tgif ",tmpInlog," > ",tmpIngif));

if Length(arg)=2 then
Exec(Concatenation(DISPLAY_PATH,tmpIngif));
Exec(Concatenation("rm ",tmpInlog,"; rm ",tmpIngif));

else

AppendTo(tmpIn2log, "Browser=",arg[3],"\n");
AppendTo(tmpIn2log,"$Browser ",tmpIngif);
Exec(Concatenation("chmod a+x ",tmpIn2log,"; ",tmpIn2log));
Exec(Concatenation("rm ",tmpInlog,"; rm ",tmpIngif,"; rm ",tmpIn2log,";"));
fi;


end);
#####################################################################


