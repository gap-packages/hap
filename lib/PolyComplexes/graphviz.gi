
#####################################################################
InstallGlobalFunction(SymmetricMatDisplay,
function(arg)
local AppendTo, PrintTo, M,t,X,V,i,j,COLOURS,tmpDir,tmpInlog,tmpIngif,tmpIn2log, mx;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

M:=StructuralCopy(arg[1]);
t:=arg[2];
COLOURS:=["blue","red","green","yellow","brown","black"];
V:=List([1..Length(M)],i->Length(COLOURS));


if Length(arg)>2 then
X:=arg[3];
if Length(X)>6 then
Print("There are too many colours required.\n");
return fail;
fi;
for i in [1..Length(X)] do
for j in X[i] do
V[j]:=i;
od;
od;
fi;

tmpDir:=DirectoryTemporary();

tmpInlog:=Filename(tmpDir,"tmpIn.log");
tmpIngif:=Filename(tmpDir,"tmpIn.gif");
tmpIn2log:=Filename(tmpDir,"tmpIn2.log");

#mx:=Maximum(Maximum(M));;
#for i in [1..Length(M)] do
#for j in [1..Length(M)] do
#if not i=j then
#M[i][j]:=1+Int(M[i][j]*100/mx);
#M[i][j]:=1+Int(M[i][j]);
#fi;
#od;od;


################ WRITE TO TMPIN.LOG #################################

#AppendTo(tmpInlog," graph G { \n size=\"4,4\" \n subgraph cluster0 {\n node [shape=ellipse, width=.2,height=.2,fixedsize=true,style=filled, color=gray35,label=\"\"] \n edge [style=\"setlinewidth(2)\"] \n");

AppendTo(tmpInlog," graph G { \n size=\"4,4\" \n subgraph cluster0 {\n node [shape=point, style=filled, label=\"\"] \n ");

for i in [1..Length(M)] do
AppendTo(tmpInlog,i, " [color=",COLOURS[V[i]],"]  \n");
od;

for i in [1..Length(M)] do
for j in [i+1..Length(M)] do

if M[i][j]<=t then
AppendTo(tmpInlog,i," -- ", j, "[len=",M[i][j]," color= ", COLOURS[Minimum(V[i],V[j])], "] \n");
else
AppendTo(tmpInlog,i," -- ", j, "[len=",M[i][j]," color=transparent] \n");
fi;

od;od;



AppendTo(tmpInlog," }\n subgraph cluster1 {\n  node [shape=box, width=2,height=1,fixedsize=true,style=filled, color=white,fillcolor=white] \n ");



AppendTo(tmpInlog,"}\n }\n");
############### WRITTEN ############################################
Exec(Concatenation(NEATO_PATH,"-Tgif ", tmpInlog," > ", tmpIngif));

Exec(Concatenation(DISPLAY_PATH, tmpIngif));
Exec(Concatenation("rm ",tmpInlog,"; rm ",tmpIngif,"; rmdir ",Filename(tmpDir,"")));


end);
#####################################################################

