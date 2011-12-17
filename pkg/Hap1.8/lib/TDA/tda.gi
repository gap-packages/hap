
#################################################################
InstallGlobalFunction(MatrixToChainComplex,
function(A)
local  row,col,i,j,x,Vertices, Edges, Faces,char,
       Boundary, BoundaryFace, BoundaryEdge, Dimension;

char:=0;
Vertices:=[];
Edges:=[];
Faces:=[];

row:=Length(A);
col:=Length(A[1]);

for i in [1..row] do
for j in [1..col] do
if A[i][j]=1 then
x:=[[2*(i-1),2*(j-1)],[2*(i-1),2*j],[2*i,2*(j-1)],[2*i,2*j]];
Append(Vertices,x);
x:=[[2*i-1,2*(j-1)],[2*i-1,2*j],[2*(i-1),2*j-1],[2*i,2*j-1]];
Append(Edges,x);
Add(Faces,[2*i-1,2*j-1]);
fi;
od;
od;

Vertices:=SSortedList(Vertices);
Edges:=SSortedList(Edges);
Faces:=SSortedList(Faces);;

######################
Dimension:=function(n);
if n = 0 then return Length(Vertices); fi;
if n = 1 then return Length(Edges); fi;
if n = 2 then return Length(Faces); fi;
return 0;
end;
######################

###
BoundaryFace:=function(e)
local Vec;
Vec:=List([1..Length(Edges)],i->0);

Vec[Position(Edges,[e[1],e[2]-1])]:=-1;  
Vec[Position(Edges,[e[1],e[2]+1])]:=1;
Vec[Position(Edges,[e[1]-1,e[2]])]:=1;
Vec[Position(Edges,[e[1]+1,e[2]])]:=-1;
return Vec;
end;
####

###
BoundaryEdge:=function(e)
local Vec;
Vec:=List([1..Length(Vertices)],i->0);

if IsEvenInt(e[1]) then
Vec[Position(Vertices,[e[1],e[2]-1])]:=-1;
Vec[Position(Vertices,[e[1],e[2]+1])]:= 1;
return Vec;
fi;

Vec[Position(Vertices,[e[1]-1,e[2]])]:=-1;
Vec[Position(Vertices,[e[1]+1,e[2]])]:= 1;
return Vec;

end;
###

####################
Boundary:=function(n,f);
if n=0 then return [0];fi;
if n=1 then return BoundaryEdge(Edges[f]); fi;
if n=2 then return BoundaryFace(Faces[f]); fi;
if n=3 then return ListWithIdenticalEntries([1..Dimension(2)],0); fi;
end;
###################


return
Objectify(HapChainComplex,
          rec(
                dimension:=Dimension,
                boundary:=Boundary,
                properties:=
                [["length",2],
                ["connected",false],
                ["type", "chainComplex"],
                ["characteristic",char]
                 ]));

end);
#################################################################
#################################################################

HAPAAA:=0;
#################################################################
#################################################################
InstallGlobalFunction(ReadImageFile,
function(file,threshold)
local i,j,prog,A;

prog:=Concatenation(GAP_ROOT_PATHS[1],"pkg/Hap1.8/lib/TDA/prog");

i:=Concatenation("convert ",file," /tmp/im.txt");
Exec(i);
i:=Concatenation("perl ",prog," /tmp/im.txt >/tmp/im.g");
Exec(i);

Read("/tmp/im.g");

for i in [1..Length(HAPAAA)] do
for j in [1..Length(HAPAAA[1])] do
if HAPAAA[i][j]<threshold then  HAPAAA[i][j]:=1; else  HAPAAA[i][j]:=0; fi;
od;od;

return HAPAAA;
end);
#################################################################
#################################################################


