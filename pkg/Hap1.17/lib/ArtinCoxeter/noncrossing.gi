#Graham Ellis

###############################################################
###############################################################
InstallGlobalFunction(NoncrossingPartitionsLatticeDisplay,
function(D)
local W,S,T,WLevel,Level,PreviousLevels,
      AppendTo, PrintTo,
      CoxeterElement,n,EltsW,NCP,tmpDir,Loggeddot,Loggedpng,x,w,y,wx,k;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

tmpDir:=DirectoryTemporary();
Loggeddot:=Filename(tmpDir,"Logged.dot");
Loggedpng:=Filename(tmpDir,"Logged.png");


W:=CoxeterDiagramFpCoxeterGroup(D);
W:=W[1]/W[2];
W:=Image(IsomorphismPermGroup(W));
S:=GeneratorsOfGroup(W);

T:=[];
for x in W do
for y in S do
AddSet(T,x*y*x^-1);
od;
od;

CoxeterElement:=Product(S);
S:=T;


##############################################
WLevel:=[];
WLevel[1]:=[Identity(W)];

n:=1;
while Length(WLevel[n])>0 do
WLevel[n+1]:=[];
PreviousLevels:=SSortedList(Flat(WLevel));  #UGHHHH!

for x in S do
for w in WLevel[n] do
wx:=w*x;
if not wx in PreviousLevels and not wx in WLevel[n+1] then 
Add( WLevel[n+1],wx);fi;
od;
od;
WLevel[n+1]:=SSortedList(WLevel[n+1]);
n:=n+1;
od;
##############################################

##############################################
Level:=function(x)
local k;
for k in [1..Length(WLevel)] do
if x in WLevel[k] then return k; fi;
od;
end;
##############################################

EltsW:=Elements(W);

NCP:=[CoxeterElement];

while not Identity(W) in NCP do
for w in NCP do
for x in S do
if Level(w*x) < Level(w) then AddSet(NCP,w*x);fi;
od;
od;
od;


PrintTo(Loggeddot,"digraph { \n node [shape=circle, style=filled, color=blue] \n "); 


for w in NCP do
for x in S do
k:=w*x;
if Level(w)<Level(k) and k in NCP then
AppendTo(Loggeddot, Position(EltsW,k), "->", Position(EltsW,w),";\n");
fi;
od;od;
AppendTo(Loggeddot,"}\n");

Exec(Concatenation("dot -Tpng ",Loggeddot," -o ",Loggedpng));
Exec(Concatenation("rm ",Loggeddot));
Exec(Concatenation(DISPLAY_PATH,Loggedpng));
Exec(Concatenation("rm ",Loggedpng));

end);
#####################################################################
#####################################################################
