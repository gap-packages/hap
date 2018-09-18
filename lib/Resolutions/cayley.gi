#(C) Graham Ellis 2005-2006

#################################################
InstallGlobalFunction(CayleyGraphOfGroup,
function(G,A)
local M, g, h, Elts;

M:=IdentityMat(Order(G))*0;
Elts:=Elements(G);

for g in [1..Order(G)] do
for h in [1..Order(G)] do
if Elts[g]*Elts[h]^-1 in A then M[g][h]:=1; fi;
if Elts[h]*Elts[g]^-1 in A then M[g][h]:=1; fi;

od;
od;

return IncidenceMatrixToGraph(M);

end);
#################################################


#####################################################################
InstallGlobalFunction(CayleyGraphOfGroupDisplay,
function(arg)
local G,X,Elts,M,i,j,COLOURS,tmpDir,tmpInlog,tmpIngif,tmpIn2log,
AppendTo, PrintTo;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

tmpDir:=DirectoryTemporary();

tmpInlog:=Filename(tmpDir,"tmpIn.log");
tmpIngif:=Filename(tmpDir,"tmpIn.gif");
tmpIn2log:=Filename(tmpDir,"tmpIn2.log");


COLOURS:=["blue","red","green","yellow","brown","black"];

G:=arg[1];
X:=arg[2];
if IsGroup(G) then
	Elts:=Elements(G);
	else
	Elts:=G;
fi;

if Length(X)>6 then
Print("There are too many generators.\n");
return fail;
fi;



M:=List([1..Length(Elts)],i->[]);

for i in [1..Length(Elts)] do
for j in [1..Length(Elts)] do
M[i][j]:=Position(X,Elts[i]^-1*Elts[j]);
od;
od;

################ WRITE TO TMPIN.LOG #################################

AppendTo(tmpInlog," graph G { \n size=\"4,4\" \n subgraph cluster0 {\n node [shape=ellipse, width=.2,height=.2,fixedsize=true,style=filled, color=gray35,label=\"\"] \n edge [style=\"setlinewidth(2)\"] \n");

for i in [1..Length(Elts)] do
for j in [1..Length(Elts)] do

if  M[i][j]=fail then
AppendTo(tmpInlog,i, " \n");
else
AppendTo(tmpInlog,i," -- ", j, "[color=",COLOURS[M[i][j]],"] \n");
fi;

od;od;

AppendTo(tmpInlog," }\n subgraph cluster1 {\n  node [shape=box, width=2,height=1,fixedsize=true,style=filled, color=white,fillcolor=white] \n ");

if Maximum(List(X,x->Length(String(x))))<20 then
for i in [1..Length(X)] do
AppendTo(tmpInlog,-i,"  [fontcolor= ",COLOURS[i],",label=\"", X[i],"\" ] \n");
od;
fi;


AppendTo(tmpInlog,"}\n }\n");
############### WRITTEN ############################################
Exec(Concatenation(NEATO_PATH,"-Tgif ", tmpInlog," > ", tmpIngif));

if Length(arg)=2 then
Exec(Concatenation(DISPLAY_PATH, tmpIngif));
Exec(Concatenation("rm ",tmpInlog,"; rm ",tmpIngif));

else

AppendTo(tmpIn2log, "Browser=",arg[3],"\n");
AppendTo(tmpIn2log,"$Browser ", tmpIngif);
Exec(Concatenation("chmod a+x ",tmpIn2log," ; ", tmpIn2log));
Exec(Concatenation("rm ",tmpInlog," ; rm ",tmpIngif,"; rm ",tmpIn2log,";"));
fi;


end);
#####################################################################


