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
InstallGlobalFunction(SymmetricMatrixToFilteredGraph,
function(arg)
local S,T,M,A,  i,j;

S:=arg[1];
T:=arg[2];
if Length(arg)=3 then M:=arg[3]; else
M:=Maximum(Flat(S));
fi;

A:=StructuralCopy(S);;
#M:=Maximum(Maximum(A));

for i in [1..Length(A)] do
for j in [1..Length(A)] do
A[i][j]:=Int(T*A[i][j]/M);
if A[i][j]>T then A[i][j]:=0; fi;
od;
A[i][i]:=1;
od;

return Objectify(HapFilteredGraph,
                rec(
		incidenceMatrix:=A,
		filtrationLength:=T,
                properties:=
                [
                ["numberofvertices",Length(A)]
                ]));
 
end);
#########################################################
#########################################################

#########################################################
#########################################################
InstallGlobalFunction(PermGroupToFilteredGraph,
function(H,Metric)
local N,S,F,i,j, GroupToSymmetricMat;

if not IsPermGroup(H) then
Print("This function must be applied to a permutation group.\n");
return fail;
fi;

N:=Maximum(MovedPoints(H));

##########################################
GroupToSymmetricMat:=function(G)
local S,x,y, elts,F;

G:=G^Random(G);
elts:=Elements(G);;
S:=[];

for x in [1..Length(elts)] do
S[x]:=[];
for y in [1..Length(elts)] do
S[x][y]:=Metric(elts[x],elts[y],N);
od;od;

return S;
end;

##########################################



S:=GroupToSymmetricMat(H);
F:=SSortedList(Flat(S));;

for i in [1..Length(S)] do
for j in [1..Length(S)] do
S[i][j]:=PositionSorted(F,S[i][j]);
od;
od;

return SymmetricMatrixToFilteredGraph(S,Length(F),Length(F));


end);
#########################################################
#########################################################



#########################################################
#########################################################
InstallGlobalFunction(FiltrationTermOfGraph,
function(G,t)
local
	A,i,j;

A:=G!.incidenceMatrix*1;;

for i in [1..Length(A)] do
for j in [1..Length(A[1])] do
if A[i][j]<=t then A[i][j]:=1; else A[i][j]:=0; fi;
od;
A[i][i]:=0;
od;

return IncidenceMatrixToGraph(A);
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
local AppendTo,
G,X,Elts,M,i,j,mm,COLOURS,tmpDir,tmpInlog,tmpIngif,tmpIn2log, cs,  s;

AppendTo:=HAP_AppendTo;

tmpDir:=DirectoryTemporary();

tmpInlog:=Filename(tmpDir,"tmpIn.log");
tmpIngif:=Filename(tmpDir,"tmpIn.gif");
tmpIn2log:=Filename(tmpDir,"tmpIn2.log");

G:=arg[1];

M:=G!.incidenceMatrix;

if IsBound(G!.clustersizes) then
################ WRITE TO TMPIN.LOG #################################
cs:=G!.clustersizes;
s:=Sum(cs);
s:=[1..s]*(4/s);
s:=(1/5)+s;
s:=1.0*s;

AppendTo(tmpInlog," graph G { \n size=\"4,4\" \n subgraph cluster0 {\n node [shape=ellipse, width=.2,height=.2,fixedsize=true,style=filled, color=gray35,label=\"\"] \n edge [style=\"setlinewidth(2)\"] \n");

for i in [1..Length(M)] do
AppendTo(tmpInlog,i, "[width=",s[cs[i]],",height=",s[cs[i]],"] \n");
od;

for i in [1..Length(M)] do
for j in [i+1..Length(M)] do

if not M[i][j]=0 then
AppendTo(tmpInlog,i," -- ", j, " \n");
fi;

od;od;

AppendTo(tmpInlog," }\n subgraph cluster1 {\n  node [shape=box, width=2,height=1,fixedsize=true,style=filled, color=white,fillcolor=white] \n ");



AppendTo(tmpInlog,"}\n }\n");
############### WRITTEN ############################################
else
################ WRITE TO TMPIN.LOG #################################

AppendTo(tmpInlog," graph G { \n size=\"4,4\" \n subgraph cluster0 {\n node [shape=ellipse, width=.2,height=.2,fixedsize=true,style=filled, color=gray35,label=\"\"] \n edge [style=\"setlinewidth(2)\"] \n");

for i in [1..Length(M)] do
AppendTo(tmpInlog,i, " \n");
od;

for i in [1..Length(M)] do
for j in [i+1..Length(M)] do

#if not M[i][j]=0 then
#AppendTo(tmpInlog,i," -- ", j, " \n");
#fi;

for mm in [1..M[i][j]] do
AppendTo(tmpInlog,i," -- ", j, " \n");
od;

od;od;

AppendTo(tmpInlog," }\n subgraph cluster1 {\n  node [shape=box, width=2,height=1,fixedsize=true,style=filled, color=white,fillcolor=white] \n ");



AppendTo(tmpInlog,"}\n }\n");
############### WRITTEN ############################################
fi;

if IsBound(G!.dot) then
Exec(Concatenation(DOT_PATH," -Tgif -Gstart=rand -Gepsilon=.000001 ", tmpInlog," > ", tmpIngif));
else
Exec(Concatenation(NEATO_PATH," -Tgif -Gstart=rand ", tmpInlog," > ", tmpIngif));
fi;

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

#####################################################################
InstallValue(ViewGraph,GraphDisplay);
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
if Length(SimplicesLst[2])=0 then dim:=0; fi;
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
if Length(SimplicesLst[3])=0 then dim:=1; fi;
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
if Length(SimplicesLst[4])=0 then dim:=2; fi;
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
if Length(SimplicesLst[5])=0 then dim:=3; fi;
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
and A[l][m]>0 and   A[l][n]>0
and   A[m][n]>0
then Add(SimplicesLst[6],[i,j,k,l,m,n]); fi;
od;od;#od;od;od;od;
SimplicesLst[6]:=SSortedList(SimplicesLst[6]);
if Length(SimplicesLst[6])=0 then dim:=4; fi;
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
if Length(SimplicesLst[d+1])=0 then dim:=d-1; fi;
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
return PositionSet(SimplicesLst[Length(v)],v);
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


##########################################################
##########################################################
InstallGlobalFunction(BarCodeOfSymmetricMatrix,
#function(M, N, m)
function(arg)
local M, N, m, F, G, D, t, step, B, bettis, i, j;

M:=arg[1];
F:=Maximum(Flat(M));
if Length(arg)>1 then N:=arg[2];
else N:=F; fi;
if Length(arg)>2 then m:=arg[3];
else m:=0; fi;

step:= Int(F/N);

t:=0;
bettis:=[];

while t<F do
G:=SymmetricMatrixToGraph(M,t);
Add(bettis,RipsHomology(G,0,2));
t:=t+step;
od;

B:=NullMat(Length(bettis),Length(bettis));
for j in [1..Length(bettis)] do
for i in [1..j] do
B[i][j]:=bettis[j];
od;
od;

return B;

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(BarCodeOfFilteredPureCubicalComplex,
function(arg)
local M, T, F, m, t, bettis, B, i, j;

M:=arg[1];
if IsBound(M!.filtrationLength) then
F:=M!.filtrationLength;
else
F:=Maximum(Flat(M!.filtration));
fi;
if Length(arg)>2 then m:=arg[3];
else m:=0; fi;

t:=1;
bettis:=[];

while t<F do
T:=FiltrationTerm(M,t);;
Add(bettis,PathComponentOfPureCubicalComplex(T,0));
t:=t+1;
od;

B:=NullMat(Length(bettis),Length(bettis));
for j in [1..Length(bettis)] do
for i in [1..j] do
B[i][j]:=bettis[j];
od;
od;

return B;


end);
##########################################################
##########################################################


#######################################
#######################################
InstallGlobalFunction(SimplicialNerveOfFilteredGraph,
function(G,dim)
local  A, T, Vertices, NrSimplices, Simplices, SimplicesLst, EnumeratedSimplex,
       bool, s, VL,x, y, d, i,j,k,l,m,n,FilteredDims,FilteredDimension, t, mx,mn;

A:=G!.incidenceMatrix;
T:=G!.filtrationLength;

Vertices:=[1..EvaluateProperty(G,"numberofvertices")];
VL:=Length(Vertices);

SimplicesLst:=List([1..dim+1],i->List([1..T],j->[]));
FilteredDims:=[];

##########################
if dim>=0 then
for i in Vertices do
Add(SimplicesLst[1][A[i][i]],[i]);
od;
for t in [1..T] do
SimplicesLst[1][t]:=SSortedList(SimplicesLst[1][t]);
od;
FilteredDims[1]:=List(SimplicesLst[1],x->Length(x));
FilteredDims[1]:=List([1..T],t->Sum(FilteredDims[1]{[1..t]}));
SimplicesLst[1]:=Concatenation(SimplicesLst[1]);
fi;

if dim>=1 then
for i in Vertices do
for j in [i+1..VL] do
if A[i][j]>0 then Add(SimplicesLst[2][A[i][j]],[i,j]); fi;
od;od;
for t in [1..T] do
SimplicesLst[2][t]:=SSortedList(SimplicesLst[2][t]);
od;
FilteredDims[2]:=List(SimplicesLst[2],x->Length(x));
FilteredDims[2]:=List([1..T],t->Sum(FilteredDims[2]{[1..t]}));
SimplicesLst[2]:=Concatenation(SimplicesLst[2]);
fi;

if dim>=2 then
for s in SimplicesLst[2] do
i:=s[1];j:=s[2];
for k in [j+1..VL] do
mn:=Minimum([ A[i][j],A[i][k],A[j][k] ]);
if mn>0 then 
mx:=Maximum([ A[i][j],A[i][k],A[j][k] ]);
Add(SimplicesLst[3][mx],[i,j,k]); fi;
od;od;
for t in [1..T] do
SimplicesLst[3][t]:=SSortedList(SimplicesLst[3][t]);
od;
FilteredDims[3]:=List(SimplicesLst[3],x->Length(x));
FilteredDims[3]:=List([1..T],t->Sum(FilteredDims[3]{[1..t]}));
SimplicesLst[3]:=Concatenation(SimplicesLst[3]);

fi;

if dim>=3 then
for s in SimplicesLst[3] do
i:=s[1];j:=s[2];k:=s[3];
for l in [k+1..VL] do
mn:=Minimum([ A[i][j], A[i][k], A[i][l], A[j][k], A[j][l], A[k][l]]);
if mn>0 then 
mx:=Maximum([ A[i][j], A[i][k], A[i][l], A[j][k], A[j][l], A[k][l]]);
Add(SimplicesLst[4][mx],[i,j,k,l]); fi;
od;od;
for t in [1..T] do
SimplicesLst[4][t]:=SSortedList(SimplicesLst[4][t]);
od;
FilteredDims[4]:=List(SimplicesLst[4],x->Length(x));
FilteredDims[4]:=List([1..T],t->Sum(FilteredDims[4]{[1..t]}));
SimplicesLst[4]:=Concatenation(SimplicesLst[4]);

fi;

if dim>=4 then
for s in SimplicesLst[4] do
i:=s[1];j:=s[2];k:=s[3];l:=s[4];
for m in [l+1..VL] do
mn:=Minimum([ A[i][j], A[i][k], A[i][l], A[i][m], A[j][k], A[j][l], A[j][m], A[k][l], A[k][m],  A[l][m]]);
if mn>0 then 
mx:=Maximum([ A[i][j], A[i][k], A[i][l], A[i][m], A[j][k], A[j][l], A[j][m], A[k][l], A[k][m],  A[l][m]]);
Add(SimplicesLst[5][mx],[i,j,k,l,m]); fi;
od;od;
for t in [1..T] do
SimplicesLst[5][t]:=SSortedList(SimplicesLst[5][t]);
od;
FilteredDims[5]:=List(SimplicesLst[5],x->Length(x));
FilteredDims[5]:=List([1..T],t->Sum(FilteredDims[5]{[1..t]}));
SimplicesLst[5]:=Concatenation(SimplicesLst[5]);

fi;

if dim>=5 then
for s in SimplicesLst[5] do
i:=s[1];j:=s[2];k:=s[3];l:=s[4];m:=s[5];
for n in [m+1..VL] do
mn:=Minimum([ A[i][j], A[i][k], A[i][l], A[i][m], A[i][n], A[j][k], A[j][l], A[j][m], A[j][n], A[k][l], A[k][m],   A[k][n], A[m][m], A[l][n], A[m][n]]);
if mn>0 then 
mx:=Maximum([ A[i][j], A[i][k], A[i][l], A[i][m], A[i][n], A[j][k], A[j][l], A[j][m], A[j][n], A[k][l], A[k][m],   A[k][n], A[m][m], A[l][n], A[m][n]]);
Add(SimplicesLst[6][mx],[i,j,k,l,m,n]); fi;
od;od;
for t in [1..T] do
SimplicesLst[6][t]:=SSortedList(SimplicesLst[6][t]);
od;
FilteredDims[6]:=List(SimplicesLst[6],x->Length(x));
FilteredDims[6]:=List([1..T],t->Sum(FilteredDims[6]{[1..t]}));
SimplicesLst[6]:=Concatenation(SimplicesLst[6]);

fi;



for d in [6..dim] do
SimplicesLst[d+1]:=[];
for y in Combinations(Vertices,d+1) do
mn:=0;
mx:=0;
for x in Combinations(y,2) do
if A[x[1]][x[2]]=0 then mn:=0; break;
else
mx:=Maximum(bool,A[x[1]][x[2]]); fi;
od;
if mx>0 then Add(SimplicesLst[d][mx],y); fi;
od;
for t in [1..T] do
SimplicesLst[d+1][t]:=SSortedList(SimplicesLst[d+1][t]);
od;
FilteredDims[d+1]:=List(SimplicesLst[d+1],x->Length(x));
FilteredDims[d+1]:=List([1..T],t->Sum(FilteredDims[d+1]{[1..t]}));
SimplicesLst[d+1]:=Concatenation(SimplicesLst[d+1]);

od;
##########################
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
EnumeratedSimplex:=function(v)
local s,t,n,i,pos;
n:=Length(v);
s:=1;
for i in [1..T] do
t:=FilteredDimension(i,n-1);
pos:=
PositionSet(SimplicesLst[n]{[s..t]},v);
if IsInt(pos) then return pos;
else s:=t+1; fi;
od;
end;
#########################

#########################
EnumeratedSimplex:=function(v);
return Position(SimplicesLst[Length(v)],v);
end;
#########################


##########################
FilteredDimension:=function(t,d);
return FilteredDims[d+1][t];
end;
#########################


Add(SimplicesLst,[]);

return
Objectify(HapFilteredSimplicialComplex,
           rec(
           vertices:=Vertices,
           nrSimplices:=NrSimplices,
           simplices:=Simplices,
           filteredDimension:=FilteredDimension,
           simplicesLst:=SimplicesLst,
           enumeratedSimplex:=EnumeratedSimplex,
           filtrationLength:=T,
           properties:=[
           ["dimension",dim]]
           ));



end);
#####################################################################


#####################################################################
#####################################################################
InstallMethod(CliqueComplex,
"Clique complex of simplicial graph",
[IsHapSimplicialComplex, IsInt],
function(K,n)
local G;

if Dimension(K)<=1 then
G:=GraphOfSimplicialComplex(K);
return SimplicialNerveOfGraph(G,n);
fi;

if Dimension(K)=2 then
return SimplicialNerveOfTwoComplex(K,n);
fi;

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(CliqueComplex,
"Clique complex of  graph",
[IsHapGraph, IsInt],
function(K,n);

return SimplicialNerveOfGraph(K,n);

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(CliqueComplex,
"Clique complex of filtered graph",
[IsHapFilteredGraph, IsInt],
function(K,n);

return SimplicialNerveOfFilteredGraph(K,n);

end);
#####################################################################
#####################################################################


########################################################
########################################################
NerveOfCover:=function(C,N)
local S, V, cnt, i, x, tog,en;

V:=[1..Length(C)];

S:=List(V,i->[i]);

cnt:=1; tog:=true;
while cnt<=N and tog=true do
cnt:=cnt+1;tog:=false;
#for x in Combinations(V,cnt) do #CHANGED 06/02/2019
en:=IteratorOfCombinations(V,cnt);
for x in en do
if Length(Intersection(List(x,xi->C[xi])))>0 then Add(S,x); tog:=true; fi;
od;
od;

S:=List(S,x->SSortedList(x));
S:=SSortedList(S);
return SimplicialComplex(S);

end;
########################################################
########################################################


########################################################
########################################################
InstallGlobalFunction(Mapper_alt,
function(S,f,N,t,tol)
local M,m,s,I,x,y,i,j,c,G,A,Clusters,Mapper,vert,K,W;
##IntervalPreCover  function of N, f, S ,tol
#N number of divisions of image
# f  function
# S  distance matrix
#tol clustering parameter
# toverlap parameter
M:=Maximum(Flat(S));;
m:=Minimum(Filtered(Flat(S),x->not x=0));;
s:=(M-m)/(N);
I:=[];
for i in [1..N] do
I[i]:=[];
x:=(i-1)*s-t*s/100;
y:=i*s+t*s/100;
for j in [1..Length(S)] do
if f(j)<=y and f(j)>=x then Add(I[i],j); fi;
od;
od;
##PreCover Done


##Mapper Complex should output the nerve of a collection of sets
G:=[]; #tol:=5;
for i in [1..N] do
A:=List(I[i],j->S[j]*1);
A:=List(A,a->a{I[i]}*1);
G[i]:=SymmetricMatrixToGraph(A,Maximum(Flat(A))/tol);
od;

Clusters:=[];
for i in [1..N] do
Clusters[i]:=[];
vert:=EvaluateProperty(G[i],"numberofvertices");
for c in 1+[1..PathComponentsOfGraph(G[i],0)] do
K:=Filtered([1..vert],k->c in G[i]!.pathComponents[k]);
K:=List(K,a->I[i][a]);
if Length(K)>0 then Add(Clusters[i],K); fi;
od;
for c in G[i]!.singletons do
Add(Clusters[i],[c]);
od;
od;


W:=Concatenation(Clusters);
Mapper:=NerveOfCover(W,20);
Mapper!.clustersizes:=List(W,Size);
return Mapper;
end);
########################################################
########################################################

###########################################
###########################################
#Mapper:=function(S,dx,f,dz,P,r,cluster)
#Mapper_alt:=function(S,f,N,t,tol)
InstallGlobalFunction(Mapper,
function(arg)
local  L, i, s , B, N, t, tol, f, S, dx, dz, cluster, r, P;

###########################
if Length(arg)=5 then
S:=arg[1]; f:=arg[2]; N:=arg[3]; t:=arg[4]; tol:=arg[5];
return Mapper_alt(S,f,N,t,tol);
fi;
###########################

S:=arg[1]; dx:=arg[2]; f:=arg[3]; dz:=arg[4]; 
P:=arg[5]; r:=arg[6]; cluster:=arg[7];
if Length(arg)=8 then B:=arg[8]; else B:=5; fi; #CHANGED 06/02/2019
L:=List(P,x->[]);
for s in S do
for i in [1..Length(P)] do
if dz(f(s),P[i]) <=r then Add(L[i],s); fi;
od;
od;

L:=List(L,x->cluster(x));
L:=Concatenation(L);
L:=Filtered(L,x->Length(x)>0);

N:=NerveOfCover(L,B);
N!.clustersizes:=List(L,Size);
N!.clusters:=L;
return N;
end);
###########################################
###########################################

#############################################
#############################################
InstallGlobalFunction(VectorsToOneSkeleton,
function(S,epsilon,dx)
local Boundaries, i, j;

Boundaries:=[];;

Boundaries[1]:=List(S,x->[1,0]);

Boundaries[2]:=[];
for i in [1..Length(S)] do
for j in [i+1..Length(S)] do
if dx(S[i],S[j])<=epsilon then Add(Boundaries[2],[2,i,j]); fi;
od;
od;

if Length(Boundaries[2])>0 then
Boundaries[3]:=[];
fi;
return RegularCWComplex(Boundaries);

end);
#############################################
#############################################

