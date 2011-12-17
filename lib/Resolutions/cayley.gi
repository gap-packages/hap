#(C) Graham Ellis 2005-2006

#####################################################################
InstallGlobalFunction(CayleyGraphDisplay,
function(arg)
local G,X,Elts,M,i,j,COLOURS;

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

AppendTo("/tmp/tmpIn.log"," graph G { \n size=\"4,4\" \n subgraph cluster0 {\n node [shape=ellipse, width=.2,height=.2,fixedsize=true,style=filled, color=gray35,label=\"\"] \n edge [style=\"setlinewidth(2)\"] \n");

for i in [1..Length(Elts)] do
for j in [1..Length(Elts)] do

if  M[i][j]=fail then
AppendTo("/tmp/tmpIn.log",i, " \n");
else
AppendTo("/tmp/tmpIn.log",i," -- ", j, "[color=",COLOURS[M[i][j]],"] \n");
fi;

od;od;

AppendTo("/tmp/tmpIn.log"," }\n subgraph cluster1 {\n  node [shape=box, width=2,height=1,fixedsize=true,style=filled, color=white,fillcolor=white] \n ");

if Maximum(List(X,x->Length(String(x))))<20 then
for i in [1..Length(X)] do
AppendTo("/tmp/tmpIn.log",-i,"  [fontcolor= ",COLOURS[i],",label=\"", X[i],"\" ] \n");
od;
fi;


AppendTo("/tmp/tmpIn.log","}\n }\n");
############### WRITTEN ############################################
Exec("neato -Tgif /tmp/tmpIn.log > /tmp/tmpIn.gif");

if Length(arg)=2 then
Exec("mozilla /tmp/tmpIn.gif");
Exec("rm /tmp/tmpIn.log; rm /tmp/tmpIn.gif");

else

AppendTo("/tmp/tmpIn2.log", "Browser=",arg[3],"\n");
AppendTo("/tmp/tmpIn2.log","$Browser /tmp/tmpIn.gif");
Exec("chmod a+x /tmp/tmpIn2.log; /tmp/tmpIn2.log");
Exec("rm /tmp/tmpIn.log; rm /tmp/tmpIn.gif; rm /tmp/tmpIn2.log;");
fi;


end);
#####################################################################


