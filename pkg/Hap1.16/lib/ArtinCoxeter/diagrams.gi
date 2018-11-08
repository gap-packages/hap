#Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(CoxeterDiagramDisplay,
function(arg)
local AppendTo, PrintTo, D,V,M,i,j,tmpDir,tmpInlog,tmpIn2log,basicgif;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

tmpDir:=DirectoryTemporary();
tmpInlog:=Filename(tmpDir,"tmpIn.log");
tmpIn2log:=Filename(tmpDir,"tmpIn2.log");
basicgif:=Filename(tmpDir,"basic.gif");

D:=arg[1];
V:=CoxeterDiagramVertices(D);
M:=function(i,j);
return CoxeterDiagramMatrix(D,i,j);
end;

################ WRITE TO TMPIN.LOG #################################

AppendTo(tmpInlog," graph G { \n size=\"4,4\" \n node [shape=circle, style=filled, color=blue,label=\" \"] \n edge [style=\"setlinewidth(5)\"] \n");

for i in V do

AppendTo(tmpInlog,i, " \n");

for j in V do

if j>=i and M(i,j)=3 then 
AppendTo(tmpInlog,i," -- ", j, "\n");
fi;

if j>=i and (M(i,j)>3 or M(i,j)<2) then
AppendTo(tmpInlog,i," -- ", j, "[label=\" ",M(i,j), " \",fontsize=20 ] \n");
fi;
od;od;

if CoxeterDiagramIsSpherical(D) then

AppendTo(tmpInlog,"10000 [label=\" Spherical Coxeter\\n Diagram \", color=white, fontsize=20,fontcolor=red,width=1.5 ]  \n");

else
AppendTo(tmpInlog,"10000 [label=\" Non-spherical Coxeter\\n Diagram \", color=white, fontsize=20, fontcolor=red,width=1.5 ]  \n");

 fi;
AppendTo(tmpInlog,"} \n");
################ WRITTEN ############################################

Exec(Concatenation(NEATO_PATH,"-Tgif ",tmpInlog," > ",basicgif));

if Length(arg)=1 then
Exec(Concatenation(DISPLAY_PATH,basicgif));
Exec(Concatenation("rm ",tmpInlog,"; rm ",basicgif));

else
AppendTo(tmpIn2log, "Browser=",arg[2],"\n");
AppendTo(tmpIn2log,"$Browser ",basicgif);
Exec(Concatenation("chmod a+x ",tmpIn2log,"; ",tmpIn2log));
Exec(Concatenation("rm ",tmpInlog,"; rm ",basicgif,"; rm ",tmpIn2log,";"));
fi;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(CoxeterDiagramMatrix,
function(arg)              #Assume i<j
local i,j,r,s,L, D, i1,j1;

D:=arg[1];
############################################
if Length(arg)=1 then     #This is unnecessarily clumsy!
L:=[];
for i in CoxeterDiagramVertices(D) do
L[i]:=[];
for j in CoxeterDiagramVertices(D) do
L[i][j]:=CoxeterDiagramMatrix(D,i,j);
od;
od;

return L;
fi;
############################################

D:=arg[1];i1:=arg[2];j1:=arg[3];
if i1<j1 then i:=i1; j:=j1; else i:=j1;j:=i1;fi;

r:=PositionProperty(D,x->(x[1]=i));
if not r=fail then L:=ShallowCopy(D[r]); L[1]:=[L[1],2];
else return 2; fi;

s:=PositionProperty(L,x->(x[1]=j));
if s=fail then return 2;
else return L[s][2]; fi;

end);
####################################################################

#####################################################################
InstallGlobalFunction(CoxeterDiagramVertices,
function(D)
local Vertices, d, i;

Vertices:=[];
for d in D do
AddSet(Vertices,d[1]);
        for i in [2..Length(d)] do
        AddSet(Vertices,d[i][1]);
        od;
od;

return Vertices;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(CoxeterDiagramFpArtinGroup,
function(D)
local 	
	F, gensF, relsF,
	Relator,
	Vertices,
	d,i,j,n;

Vertices:=CoxeterDiagramVertices(D);

F:=FreeGroup(Length(Vertices));
gensF:=GeneratorsOfGroup(F);
relsF:=[];

	#############################################################
	Relator:=function(x,y,n)
	local l,r,N;

	l:=Identity(F);
	r:=Identity(F);
	N:=ShallowCopy(n);

	while N>0 do
	if IsOddInt(N) then 
	l:=l*y; r:=r*x;
	else
	l:=l*x; r:=r*y;
	fi;
	N:=N-1;
	od;

	if IsOddInt(n) then return r*l^-1; else return l*r^-1; fi;
	end;
	#############################################################

for i in [1..Length(Vertices)] do
for j in [i+1..Length(Vertices)] do
Append(relsF,
[Relator(gensF[i],gensF[j],CoxeterDiagramMatrix(D,Vertices[i],Vertices[j]))]);
od;
od;

return [F,relsF];
end);
#####################################################################

#####################################################################
InstallGlobalFunction(CoxeterDiagramFpCoxeterGroup,
function(D)
local L,x,F,relsF;

L:=CoxeterDiagramFpArtinGroup(D);
F:=L[1];
relsF:=L[2];

for x in GeneratorsOfGroup(F) do
Append(relsF,[x^2]);
od;

return [F,relsF];
end);
#####################################################################

#####################################################################
InstallGlobalFunction(CoxeterSubDiagram,
function(D,V)
local SD, n, u, v, singles, temp;

SD:=[];

for v in V do
SD[v]:=[v];
for u in V do
if u>v then
n:=CoxeterDiagramMatrix(D,u,v);
	if not n=2 then
	Append(SD[v],[[u,n]]);
	fi;
fi;
od;
od;

SD:=SortedList(SD);
singles:=Filtered(SD,x->Length(x)=1);

for u in singles do
temp:=ShallowCopy(SD);
RemoveSet(temp,u);
if Length(CoxeterDiagramVertices(temp))=Length(CoxeterDiagramVertices(SD)) then 
SD:=temp; fi;
od;

return SD;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(CoxeterDiagramComponents,
function(D)
local Components, RemainingVertices, C, Count, Vert, n, u, v, w;

Components:=[];
RemainingVertices:=CoxeterDiagramVertices(D);

while not Length(RemainingVertices)=0 do

Count:=0;
C:=CoxeterSubDiagram(D,[RemainingVertices[1]]);

	while Length(CoxeterDiagramVertices(C))>Count do
	Vert:=CoxeterDiagramVertices(C);
	Count:=Length(Vert);
	
	for v in CoxeterDiagramVertices(C) do
	for u in RemainingVertices do
	n:=CoxeterDiagramMatrix(D,v,u);
	if not n=2 then 
	AddSet(Vert,u);
	fi;
	od;
	od;

	C:=CoxeterSubDiagram(D,Vert);
	od;

C:=CoxeterSubDiagram(D,Vert);
Append(Components,[C]);
for v in CoxeterDiagramVertices(C) do
RemoveSet(RemainingVertices,v);
od;

od;

return Components;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(CoxeterDiagramDegree,
function(D,v)
local u, degree;

degree:=0;
for u in CoxeterDiagramVertices(D) do
if not CoxeterDiagramMatrix(D,u,v)=2 then degree:=degree+1; fi;
od;

return degree;
end);
#####################################################################



#####################################################################
InstallGlobalFunction(CoxeterDiagramIsSpherical,
function(D)
local C, Components, ComponentIsSpherical;

Components:=CoxeterDiagramComponents(D);

#####################################################################
ComponentIsSpherical:=function(C)
local
	Vertices,
	Verts,
	Degrees,
	Labels,
	LargeLabels,
	LargeDegrees,
	Leaves,
	e,n,u,v,w,N;

Vertices:=CoxeterDiagramVertices(C);
Degrees:=List(Vertices,x->CoxeterDiagramDegree(D,x));
LargeDegrees:=Filtered(Degrees,x->x>2);
Leaves:=Filtered(Degrees,x->x=1);

Labels:=[];
for u in Vertices do
for v in Vertices do
n:=CoxeterDiagramMatrix(D,u,v);
if (not n=2) and (u>v) then Append(Labels,[n]); fi;
od;
od;

LargeLabels:=Filtered(Labels,x->x>3);

if 0 in Labels then return false; fi;
if Length(Vertices)<3 then return true; fi;
if Length(LargeDegrees)>1 then return false; fi;
if Length(Leaves)<2 or Length(Leaves)>3 then return false; fi;
if Length(LargeLabels)>1 then return false; fi;
if Length(LargeDegrees)>0 and Length(LargeLabels)>0 then return false; fi;

  
if (not Length(LargeLabels)=0) then		 
	if Maximum(LargeLabels)>5 then return false; fi;
fi;
if (not Length(Degrees)=0) then
        if (Maximum(Degrees)>3) then return false; fi;
	fi;

# From this point on we can assume that:
# the Coxeter diagram is connected;
# there are more than two vertices;
# either two or three vertices have degree 1;
# no vertex has degree > 3;
# at most one vertex has degree 3 and, in this case, all labels are <4;
# no edge label is > 5 (and so no label is infinity=0);
# at most one edge label is >3.
# So now its just a case by case slog! 

# Case A0 (Forgot this in initia versions of HAP)
if Length(LargeDegrees)>0 then

Verts:=Filtered(Vertices,v->not v=
Vertices[Position(Degrees,LargeDegrees[1])]);
Verts:=List(
CoxeterDiagramComponents(CoxeterSubDiagram(D,Verts)),
x->Length(CoxeterDiagramVertices(x)));
Verts:=Filtered(Verts,x->x>1);
if Length(Verts)>2 then return false; fi;
fi;


# Case A
if Maximum(Labels)=3 and Maximum(Degrees)=2 then return true; fi;

# Case B
if Maximum(Labels)=4 and Maximum(Degrees)=2 then
	for u in Vertices do
	for v in Vertices do
	if CoxeterDiagramMatrix(D,u,v)=4 then e:=[u,v];fi;
	od;
	od;
	if (CoxeterDiagramDegree(D,e[1])=1) or (CoxeterDiagramDegree(D,e[2])=1) then
	return true;
	fi;
fi;

# Case D
if Maximum(Labels)=3 and Maximum(Degrees)=3 then
	for v in Vertices do
	if CoxeterDiagramDegree(D,v)=3 then w:=v; break; fi;
	od;
	e:=0;
	for v in Vertices do
	if CoxeterDiagramMatrix(D,w,v)=3 then e:=e+CoxeterDiagramDegree(D,v);fi;
	od;
	if e<5 then return true; fi;
fi;

# Case H3
if Length(Vertices)=3 and Maximum(Labels)=5 then return true; fi;

# Case H4
if Length(Vertices)=4 and Maximum(Labels)=5 then
	for u in Vertices do
        for v in Vertices do
	        if CoxeterDiagramMatrix(D,u,v)=5 then e:=[u,v];fi;
        od;
        od;
        if (CoxeterDiagramDegree(D,e[1])=1) or (CoxeterDiagramDegree(D,e[2])=1) then
        return true;
        fi;
fi;

# Case F4
if Length(Vertices)=4 and Maximum(Labels)=4 then
        for u in Vertices do
        for v in Vertices do
                if CoxeterDiagramMatrix(D,u,v)=4 then e:=[u,v];fi;
        od;
        od;
        if (CoxeterDiagramDegree(D,e[1])=2) and (CoxeterDiagramDegree(D,e[2])=2) then
        return true;
        fi;
fi;


# Cases E6, E7 and E8
if (Length(Vertices)=6 or Length(Vertices)=7 or Length(Vertices)=8)
and Maximum(Labels)=3 and Maximum(Degrees)=3 then

	for v in Vertices do
        if CoxeterDiagramDegree(D,v)=3 then w:=v; break; fi;
   	od;
	
        e:=[];
        
	for v in Vertices do
        if CoxeterDiagramMatrix(D,w,v)=3 then 
		Append(e,[v]);
	fi;
        od;
	
	e:=Filtered(e,x->CoxeterDiagramDegree(D,x)=2);
	if Length(e)<2 then return false; fi;

	for v in e do
	for u in Vertices do
	if CoxeterDiagramMatrix(D,u,v)=3 and CoxeterDiagramDegree(D,u)=1 then 
		return true;
	fi;
	od;
	od;
fi;

return false;
end;
#####################################################################

for C in Components do
if not ComponentIsSpherical(C) then return false; fi;
od;

return true;
end);
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(CoxeterDiagramMatCoxeterGroup,
function(D)
local M, cosine, G, i,j;

M:=CoxeterDiagramMatrix(D);
M:=M-IdentityMat(Length(M));

#################################
cosine := function(n)
if n>0 then
return (1/2)*( E(2*n) + E(2*n)^(2*n-1) );
elif n = 0 then  return 1; 
fi;
end;
#################################

G := [];
for i in [1..Length(M)] do
G[i] := IdentityMat(Length(M));
for j in [1..Length(M)] do
G[i][i][j] := G[i][i][j] + 2*cosine(M[i][j]);
od;
od;

return Group(G);
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallMethod(CoxeterMatrix,
"Coxeter matrix of Coxeter diagram",
[IsList],
function(L) local A;
A:= CoxeterDiagramMatrix(L);
A:=A-IdentityMat(Length(A));
return A;
end);
#####################################################################
#####################################################################

