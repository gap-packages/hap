#(C) Graham Ellis 2009

#########################################################
#########################################################
InstallGlobalFunction(SymmetricMatrixToIncidenceMatrix,
function(arg)
local
   	S,t,M,SortedM,len,i,j,T; 

S:=arg[1];
t:=arg[2];

M:=StructuralCopy(S);

if Length(arg)=2 then
len:=Length(M);
for i in [1..len] do
for j in [1..len] do
if M[i][j]>t then M[i][j]:=0;
else M[i][j]:=1;
fi;
od;
od;
fi;

if Length(arg)>2 then
T:=arg[3];
SortedM:=List(M,row->SSortedList(SortedList(row){[1..T]}));

len:=Length(M);
for i in [1..len] do
for j in [1..len] do
if S[i][j]>t  then M[i][j]:=0; 
else
if S[i][j] in SortedM[i] then M[i][j]:=1; M[j][i]:=1;  
else M[i][j]:=0; M[j][i]:=0; fi;
fi;
od;
od;
fi;


for i in [1..len] do
M[i][i]:=0;
od;

return M;

end);
#########################################################
#########################################################

#########################################################
#########################################################
InstallGlobalFunction(IncidenceMatrixToGraph,
function(M)
	local vertices;

vertices:=Length(M);

return 
 Objectify(HapGraph,
                rec(
                incidenceMatrix:=StructuralCopy(M),
                properties:=
                [
                ["numberofvertices",vertices]
                ]));


end);
#########################################################
#########################################################

#########################################################
#########################################################
InstallGlobalFunction(SymmetricMatrixToGraph,
function(arg)
local S,t,T;

S:=arg[1];
t:=arg[2];

if Length(arg)=2 then
return 
IncidenceMatrixToGraph(SymmetricMatrixToIncidenceMatrix(S,t));
fi;

T:=arg[3];
return
IncidenceMatrixToGraph(SymmetricMatrixToIncidenceMatrix(S,t,T));
end);
#########################################################
#########################################################

#########################################################
#########################################################
InstallGlobalFunction(PathComponentsOfGraph,
function(G,n)
local  MaximalConnectedGraph, BettiZero, COLOUR,singletons, FindOne, bool, c, i, j, v, M, A;

M:=StructuralCopy(G!.incidenceMatrix);

if not "pathComponents" in NamesOfComponents(G) then
################
################
################
singletons:=[];
for i in [1..Length(M)] do
if Sum(M[i])=0 then Add(singletons,i); fi;
od;


##################
FindOne:=function(M)
local i,j;

for i in [1..Length(M)] do
for j in [1..Length(M)] do
if M[i][j]=1 then return [i,j];fi;
od;od;
return false;

end;
##################

################################################
MaximalConnectedGraph:=function(M,v)
local i,j,colour,len,leaves,newleaves,vertices;
#Find a maximal tree containing the free edge v.
#Colour the edges of the tree.

len:=Length(M);
colour:=Maximum(List(M,x->Maximum(x)))+1;
M[v[1]][v[2]]:=colour;
M[v[2]][v[1]]:=colour;
leaves:=v;
vertices:=SSortedList(v);

while Size(leaves)>0 do
newleaves:=[];
for i in leaves do
for j in [1..len] do
if M[i][j]=1 then 
M[i][j]:=colour; M[j][i]:=colour; 
   if not j in vertices then
       AddSet(vertices,j);
       Add(newleaves,j);
   fi;
fi;
od;
od;
leaves:=newleaves;
od;

return colour;
end;
################################################

COLOUR:=1;
v:=FindOne(M);
while not v=false do
COLOUR:=MaximalConnectedGraph(M,v);
v:=FindOne(M);
od;

BettiZero:=Length(singletons)+COLOUR-1;
G!.bettiZero:=BettiZero;
G!.pathComponents:=M;
G!.singletons:=singletons;
####################
####################
####################
fi;

if n=0 then
return G!.bettiZero;
fi;
if n>G!.bettiZero then
Print("The are not so many path components!\n");
return fail;
fi;
if n<=Length(G!.singletons) then return
IncidenceMatrixToGraph([[0]]);
fi;

M:=G!.pathComponents;
c:=n-Length(G!.singletons)+1;
A:=Filtered(M,row->c in row);
A:=MutableCopyMat(TransposedMat(A));
A:=Filtered(A,row->c in row);
for i in [1..Length(A)] do
for j in [1..Length(A)] do
if not A[i][j]=0 then A[i][j]:=1; fi;
od;od;

return IncidenceMatrixToGraph(A);

end);
#########################################################
#########################################################





###########################################
###########################################
ReadBioData:=function(file)
local L,M,fun,T,t,it;

it:=InputTextFile(file);
Sleep(1);
L:=ReadLine(it);;
M:=[];

#################
fun:=function(x);
if x='\t' then return ',';
else return x; fi;
end;
#################

while true do
L:=ReadLine(it);;
if L=fail then break; fi;
t:=Position(L,'\t');
L:=L{[t..Length(L)-1]};
T:=List(L,fun);
T[1]:='[';
Add(T,']');
T:=Filtered(T,a->not a='.');
T:=EvalString(T);   
Add(M,T);
od;


return M;

end;
#######################################
#######################################


#######################################
#######################################
InstallGlobalFunction(VectorsToSymmetricMatrix,
function(arg)
local M, Distance,S,i,j;

#############################
M:=arg[1];
if Length(arg)=1 then
	Distance:=function(u,v);
	return Sum(List(u-v,x->AbsInt(x)));
	end;
else Distance:=arg[2];
fi;
#############################

S:=[];
for i in [1..Length(M)] do
S[i]:=[];
for j in [1..Length(M)] do
S[i][j]:=Distance(M[i],M[j]);
od;
od;

return S;
end);
#######################################
#######################################






#####################################################################
InstallGlobalFunction(GraphDisplay,
function(arg)
local G,X,Elts,M,i,j,COLOURS,tmpDir,tmpInlog,tmpIngif,tmpIn2log;

tmpDir:=DirectoryTemporary();

tmpInlog:=Filename(tmpDir,"tmpIn.log");
tmpIngif:=Filename(tmpDir,"tmpIn.gif");
tmpIn2log:=Filename(tmpDir,"tmpIn2.log");

G:=arg[1];

M:=G!.incidenceMatrix;


################ WRITE TO TMPIN.LOG #################################

AppendTo(tmpInlog," graph G { \n size=\"4,4\" \n subgraph cluster0 {\n node [shape=ellipse, width=.2,height=.2,fixedsize=true,style=filled, color=gray35,label=\"\"] \n edge [style=\"setlinewidth(2)\"] \n");

for i in [1..Length(M)] do
for j in [i+1..Length(M)] do

if  M[i][j]=0 then
AppendTo(tmpInlog,i, " \n");
else
AppendTo(tmpInlog,i," -- ", j, " \n");
fi;

od;od;

AppendTo(tmpInlog," }\n subgraph cluster1 {\n  node [shape=box, width=2,height=1,fixedsize=true,style=filled, color=white,fillcolor=white] \n ");



AppendTo(tmpInlog,"}\n }\n");
############### WRITTEN ############################################
Exec(Concatenation(NEATO_PATH,"-Tgif ", tmpInlog," > ", tmpIngif));

if Length(arg)=1 then
Exec(Concatenation(DISPLAY_PATH, tmpIngif));
Exec(Concatenation("rm ",tmpInlog,"; rm ",tmpIngif));

else

AppendTo(tmpIn2log, "Browser=",arg[2],"\n");
AppendTo(tmpIn2log,"$Browser ", tmpIngif);
Exec(Concatenation("chmod a+x ",tmpIn2log," ; ", tmpIn2log));
Exec(Concatenation("rm ",tmpInlog," ; rm ",tmpIngif,"; rm ",tmpIn2log,";"));
fi;


end);
#####################################################################

#######################################
#######################################
InstallGlobalFunction(SimplicialNerveOfGraph,
function(G,dim)
local  A, Vertices, NrSimplices, Simplices, SimplicesLst, EnumeratedSimplex,
       bool, s, VL,x, y, d, i,j,k,l,m,n,RedundantFaces;

A:=G!.incidenceMatrix;

Vertices:=[1..EvaluateProperty(G,"numberofvertices")];
VL:=Length(Vertices);

SimplicesLst:=[];

##########################
if dim>=0 then
SimplicesLst[1]:=List(Vertices,x->[x]); 
fi;

if dim>=1 then
SimplicesLst[2]:=[];
for i in Vertices do
for j in [i+1..VL] do
if A[i][j]>0 then Add(SimplicesLst[2],[i,j]); fi;
od;od;
SimplicesLst[2]:=SSortedList(SimplicesLst[2]);
fi;

if dim>=2 then
SimplicesLst[3]:=[];
#for i in Vertices do
#for j in [i+1..VL] do
for s in SimplicesLst[2] do
i:=s[1];j:=s[2];
for k in [j+1..VL] do
if A[i][j]>0 and A[i][k]>0 and A[j][k]>0
then Add(SimplicesLst[3],[i,j,k]); fi;
od;od;#od;
SimplicesLst[3]:=SSortedList(SimplicesLst[3]);
fi;

if dim>=3 then
SimplicesLst[4]:=[];
#for i in Vertices do
#for j in [i+1..VL] do
#for k in [j+1..VL] do
for s in SimplicesLst[3] do
i:=s[1];j:=s[2];k:=s[3];
for l in [k+1..VL] do
if A[i][j]>0 and A[i][k]>0 and A[i][l]>0 and
A[j][k]>0 and A[j][l]>0 and A[k][l]>0
then Add(SimplicesLst[4],[i,j,k,l]); fi;
od;od;#od;od;
SimplicesLst[4]:=SSortedList(SimplicesLst[4]);
fi;

if dim>=4 then
SimplicesLst[5]:=[];
#for i in Vertices do
#for j in [i+1..VL] do
#for k in [j+1..VL] do
#for l in [k+1..VL] do
for s in SimplicesLst[4] do
i:=s[1];j:=s[2];k:=s[3];l:=s[4];
for m in [l+1..VL] do
if A[i][j]>0 and A[i][k]>0 and A[i][l]>0 and A[i][m]>0
and A[j][k]>0 and A[j][l]>0 and A[j][m]>0 
and A[k][l]>0 and A[k][m]>0
and A[l][m]>0
then Add(SimplicesLst[5],[i,j,k,l,m]); fi;
od;od;#od;od;od;
SimplicesLst[5]:=SSortedList(SimplicesLst[5]);
fi;

if dim>=5 then
SimplicesLst[6]:=[];
#for i in Vertices do
#for j in [i+1..VL] do
#for k in [j+1..VL] do
#for l in [k+1..VL] do
#for m in [l+1..VL] do
for s in SimplicesLst[5] do
i:=s[1];j:=s[2];k:=s[3];l:=s[4];m:=s[5];
for n in [m+1..VL] do
if A[i][j]>0 and A[i][k]>0 and A[i][l]>0 and A[i][m]>0 and A[i][n]>0
and A[j][k]>0 and A[j][l]>0 and A[j][m]>0  and A[j][n]>0
and A[k][l]>0 and A[k][m]>0 and   A[k][n]>0
and A[m][m]>0 and   A[l][n]>0
and   A[m][n]>0
then Add(SimplicesLst[6],[i,j,k,l,m,n]); fi;
od;od;#od;od;od;od;
SimplicesLst[6]:=SSortedList(SimplicesLst[6]);
fi;



for d in [6..dim] do
SimplicesLst[d+1]:=[];
for y in Combinations(Vertices,d+1) do
bool:=true;
for x in Combinations(y,2) do
if A[x[1]][x[2]]=0 then bool:=false; break; fi;
od;
if bool then Add(SimplicesLst[d],y); fi;
od;
SimplicesLst[d+1]:=SSortedList(SimplicesLst[d+1]);
od;
##########################

##########################
##########################
if dim=2 then 
RedundantFaces:=[];

for i in [1..Length(A)] do
for j in [i+1..Length(A)] do
if not A[i][j]=0 then
for k in [j+1..Length(A)] do
if not A[i][j]*A[j][k]*A[i][k]=0 then
for l in [k+1..Length(A)] do
if not A[i][j]*A[i][k]*A[i][l]*A[j][k]*A[j][l]*A[k][l] = 0 then
Add(RedundantFaces,[i,j,k]);
break;
fi;
od;fi;od;fi;od;od;

SimplicesLst[dim]:=Difference(SimplicesLst[dim],RedundantFaces);
SimplicesLst[dim]:=SSortedList(SimplicesLst[dim]);
fi;
###########################
###########################

##########################
##########################
if dim=3 then
RedundantFaces:=[];

for i in [1..Length(A)] do
for j in [i+1..Length(A)] do
if not A[i][j]=0 then
for k in [j+1..Length(A)] do
if not A[i][j]*A[j][k]*A[i][k]=0 then
for l in [k+1..Length(A)] do
if not A[i][j]*A[i][k]*A[i][l]*A[j][k]*A[j][l]*A[k][l] = 0 then
for m in [l+1..Length(A)] do
if not A[i][j]*A[i][k]*A[i][l]*A[i][m]*A[j][k]*A[j][l]*A[j][m]*A[k][l]*A[k][m]*A[l][m] = 0 then
Add(RedundantFaces,[i,j,k,l]);
break;
fi;
od;fi;od;fi;od;fi;od;od;

SimplicesLst[dim]:=Difference(SimplicesLst[dim],RedundantFaces);
SimplicesLst[dim]:=SSortedList(SimplicesLst[dim]);
fi;
###########################
###########################


##########################
Simplices:=function(d,k);
return SimplicesLst[d+1][k];
end;
#########################


##########################
NrSimplices:=function(d);
return Length(SimplicesLst[d+1]);
end;
#########################

#########################
EnumeratedSimplex:=function(v);
return Position(SimplicesLst[Length(v)],v);
end;
#########################

Add(SimplicesLst,[]);

return
Objectify(HapSimplicialComplex,
           rec(
           vertices:=Vertices,
           nrSimplices:=NrSimplices,
           simplices:=Simplices,
           simplicesLst:=SimplicesLst,
           enumeratedSimplex:=EnumeratedSimplex,
           properties:=[
           ["dimension",dim]]
           ));



end);
#####################################################################

####################################################
InstallGlobalFunction(DensityMat,
function(S,k,T)
local A,s,p,F,B;

A:=[];

for s in S do
Add(A,SortedList(s)[k]);
od;

p:=SortedList(A)[T];

F:=Filtered([1..Length(S)],i->A[i]<=p);

B:=StructuralCopy(S);
B:=B{F};
B:=TransposedMat(B);
B:=MutableCopyMat(B);
B:=B{F};

return B;
end);
####################################################

