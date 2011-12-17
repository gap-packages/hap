#####################################################################
#####################################################################
InstallGlobalFunction(ContractMatrix,
function(arg)
local
	ChooseRandomBox,
	Neighbours,
	MOVES,
	Contract,
	ContractImageMatrix,
	TOG1,
	TOG2;
A:=arg[1];
#####################################################################
ChooseRandomBox:=function(A)
local i,j, col, row;
			#Returns a random (=first!) black box
row:=Length(A);
col:=Length(A[1]);

for i in [1..row] do
for j in [1..col] do
if A[i][j]=1 then return [i,j]; fi;
od;
od;

return fail;

end;
#####################################################################

#####################################################################
Neighbours:=function(A,box)
local nbrs,left, top,right,bottom,cols,rows,topleft,topright,bottomleft,bottomright;

cols:=Length(A[1]);
rows:=Length(A);

nbrs:=[];

if box[2]=1 then left:=0; else left:=A[box[1]][box[2]-1]; fi;
if box[2]=cols then right:=0; else right:=A[box[1]][box[2]+1]; fi;
if box[1]=1 then top:=0; else top:=A[box[1]-1][box[2]]; fi;
if box[1]=rows then bottom:=0; else bottom:=A[box[1]+1][box[2]]; fi;

if box[2]=1 or box[1]=1 
then topleft:=0; else topleft:=A[box[1]-1][box[2]-1]; fi;

if box[2]=cols or box[1]=1 then topright:=0; else topright:=A[box[1]-1][box[2]+1]; fi;

if box[2]=1 or box[1]=rows then bottomleft:=0; else bottomleft:=A[box[1]+1][box[2]-1]; fi;

if box[2]=cols or box[1]=rows then bottomright:=0; else bottomright:=A[box[1]+1][box[2]+1]; fi;


return [left,topleft,top,topright,right,bottomright,bottom,bottomleft];

end;
######################################################################

MOVES:=[

[ 0, 0, 0, 0, 0, 1, 0, 0 ], 
[ 0, 0, 0, 0, 0, 0, 1, 0 ],
[ 0, 0, 0, 0, 0, 1, 1, 0 ], 
[ 0, 0, 0, 0, 0, 0, 0, 1 ],
[ 0, 0, 0, 0, 0, 0, 1, 1 ], 
[ 0, 0, 0, 0, 0, 1, 1, 1 ],
[ 0, 0, 0, 0, 1, 0, 0, 0 ], 
[ 0, 0, 0, 0, 1, 1, 0, 0 ],
[ 0, 0, 0, 0, 1, 0, 1, 0 ], 
[ 0, 0, 0, 0, 1, 1, 1, 0 ],
[ 0, 0, 0, 0, 1, 0, 1, 1 ], 
[ 0, 0, 0, 0, 1, 1, 1, 1 ],
[ 1, 0, 0, 0, 0, 0, 0, 0 ], 
[ 1, 0, 0, 0, 0, 0, 1, 0 ],
[ 1, 0, 0, 0, 0, 1, 1, 0 ], 
[ 1, 0, 0, 0, 0, 0, 0, 1 ],
[ 1, 0, 0, 0, 0, 0, 1, 1 ], 
[ 1, 0, 0, 0, 0, 1, 1, 1 ],
[ 1, 0, 0, 0, 1, 0, 1, 0 ], 
[ 1, 0, 0, 0, 1, 1, 1, 0 ],
[ 1, 0, 0, 0, 1, 0, 1, 1 ], 
[ 1, 0, 0, 0, 1, 1, 1, 1 ],
[ 0, 0, 0, 1, 0, 0, 0, 0 ], 
[ 0, 0, 0, 1, 1, 0, 0, 0 ],
[ 0, 0, 0, 1, 1, 1, 0, 0 ], 
[ 0, 0, 0, 1, 1, 0, 1, 0 ],
[ 0, 0, 0, 1, 1, 1, 1, 0 ], 
[ 0, 0, 0, 1, 1, 0, 1, 1 ],
[ 0, 0, 0, 1, 1, 1, 1, 1 ], 
[ 1, 0, 0, 1, 1, 0, 1, 0 ],
[ 1, 0, 0, 1, 1, 1, 1, 0 ], 
[ 1, 0, 0, 1, 1, 0, 1, 1 ],
[ 1, 0, 0, 1, 1, 1, 1, 1 ], 
[ 0, 0, 1, 0, 0, 0, 0, 0 ],
[ 0, 0, 1, 0, 1, 0, 0, 0 ], 
[ 0, 0, 1, 0, 1, 1, 0, 0 ],
[ 0, 0, 1, 0, 1, 0, 1, 0 ], 
[ 0, 0, 1, 0, 1, 1, 1, 0 ],
[ 0, 0, 1, 0, 1, 0, 1, 1 ], 
[ 0, 0, 1, 0, 1, 1, 1, 1 ],
[ 1, 0, 1, 0, 0, 0, 0, 0 ], 
[ 1, 0, 1, 0, 0, 0, 1, 0 ],
[ 1, 0, 1, 0, 0, 1, 1, 0 ], 
[ 1, 0, 1, 0, 0, 0, 0, 1 ],
[ 1, 0, 1, 0, 0, 0, 1, 1 ], 
[ 1, 0, 1, 0, 0, 1, 1, 1 ],
[ 1, 0, 1, 0, 1, 0, 0, 0 ], 
[ 1, 0, 1, 0, 1, 1, 0, 0 ],
[ 1, 0, 1, 0, 1, 0, 0, 1 ], 
[ 1, 0, 1, 0, 1, 1, 0, 1 ],
[ 0, 0, 1, 1, 0, 0, 0, 0 ], 
[ 0, 0, 1, 1, 1, 0, 0, 0 ],
[ 0, 0, 1, 1, 1, 1, 0, 0 ], 
[ 0, 0, 1, 1, 1, 0, 1, 0 ],
[ 0, 0, 1, 1, 1, 1, 1, 0 ], [ 0, 0, 1, 1, 1, 0, 1, 1 ],
[ 0, 0, 1, 1, 1, 1, 1, 1 ], [ 1, 0, 1, 1, 0, 0, 0, 0 ],
[ 1, 0, 1, 1, 0, 0, 1, 0 ], [ 1, 0, 1, 1, 0, 1, 1, 0 ],
[ 1, 0, 1, 1, 0, 0, 0, 1 ], [ 1, 0, 1, 1, 0, 0, 1, 1 ],
[ 1, 0, 1, 1, 0, 1, 1, 1 ], [ 1, 0, 1, 1, 1, 0, 0, 0 ],
[ 1, 0, 1, 1, 1, 1, 0, 0 ], [ 1, 0, 1, 1, 1, 0, 0, 1 ],
[ 1, 0, 1, 1, 1, 1, 0, 1 ], [ 0, 1, 0, 0, 0, 0, 0, 0 ],
[ 1, 1, 0, 0, 0, 0, 0, 0 ], [ 1, 1, 0, 0, 0, 0, 1, 0 ],
[ 1, 1, 0, 0, 0, 1, 1, 0 ], [ 1, 1, 0, 0, 0, 0, 0, 1 ],
[ 1, 1, 0, 0, 0, 0, 1, 1 ], [ 1, 1, 0, 0, 0, 1, 1, 1 ],
[ 1, 1, 0, 0, 1, 0, 1, 0 ], [ 1, 1, 0, 0, 1, 1, 1, 0 ],
[ 1, 1, 0, 0, 1, 0, 1, 1 ], [ 1, 1, 0, 0, 1, 1, 1, 1 ],
[ 1, 1, 0, 1, 1, 0, 1, 0 ], [ 1, 1, 0, 1, 1, 1, 1, 0 ],
[ 1, 1, 0, 1, 1, 0, 1, 1 ], [ 1, 1, 0, 1, 1, 1, 1, 1 ],
[ 0, 1, 1, 0, 0, 0, 0, 0 ], [ 0, 1, 1, 0, 1, 0, 0, 0 ],
[ 0, 1, 1, 0, 1, 1, 0, 0 ], [ 0, 1, 1, 0, 1, 0, 1, 0 ],
[ 0, 1, 1, 0, 1, 1, 1, 0 ], [ 0, 1, 1, 0, 1, 0, 1, 1 ],
[ 0, 1, 1, 0, 1, 1, 1, 1 ], [ 1, 1, 1, 0, 0, 0, 0, 0 ],
[ 1, 1, 1, 0, 0, 0, 1, 0 ], [ 1, 1, 1, 0, 0, 1, 1, 0 ],
[ 1, 1, 1, 0, 0, 0, 0, 1 ], [ 1, 1, 1, 0, 0, 0, 1, 1 ],
[ 1, 1, 1, 0, 0, 1, 1, 1 ], [ 1, 1, 1, 0, 1, 0, 0, 0 ],
[ 1, 1, 1, 0, 1, 1, 0, 0 ], [ 1, 1, 1, 0, 1, 0, 0, 1 ],
[ 1, 1, 1, 0, 1, 1, 0, 1 ], [ 0, 1, 1, 1, 0, 0, 0, 0 ],
[ 0, 1, 1, 1, 1, 0, 0, 0 ], [ 0, 1, 1, 1, 1, 1, 0, 0 ],
[ 0, 1, 1, 1, 1, 0, 1, 0 ], [ 0, 1, 1, 1, 1, 1, 1, 0 ],
[ 0, 1, 1, 1, 1, 0, 1, 1 ], [ 0, 1, 1, 1, 1, 1, 1, 1 ],
[ 1, 1, 1, 1, 0, 0, 0, 0 ], [ 1, 1, 1, 1, 0, 0, 1, 0 ],
[ 1, 1, 1, 1, 0, 1, 1, 0 ], [ 1, 1, 1, 1, 0, 0, 0, 1 ],
[ 1, 1, 1, 1, 0, 0, 1, 1 ], [ 1, 1, 1, 1, 0, 1, 1, 1 ],
[ 1, 1, 1, 1, 1, 0, 0, 0 ], [ 1, 1, 1, 1, 1, 1, 0, 0 ],
[ 1, 1, 1, 1, 1, 0, 0, 1 ], [ 1, 1, 1, 1, 1, 1, 0, 1 ]

];;
MOVES:=SSortedList(MOVES);

######################################################################
Contract:=function(A,box,Moves);

if Neighbours(A,box) in Moves
then A[box[1]][box[2]]:=0; return true;else 
return false; fi;
end;
######################################################################

######################################################################
ContractImageMatrix:=function(A,Moves)
local i,j,row,col,tog,ii,jj;
row:=Length(A);
col:=Length(A[1]);
tog:=false;

for i in [1..row] do
for j in [1..col] do
if A[i][j]=1 then
tog:=Contract(A,[i,j],Moves) or tog ;
fi;
od;
od;

for i in [1..row] do
for j in [1..col] do
ii:=row+1-i;jj:=col+1-j;
if A[ii][jj]=1 then 
tog:=Contract(A,[ii,jj],Moves) or tog;
fi;
od;
od;

return tog;
end;
######################################################################

TOG2:=false;
TOG1:=true;
while TOG1 do
TOG1:= ContractImageMatrix(A,MOVES);
if TOG1 then TOG2:=true;fi;
od;

return TOG2;

end);
######################################################################
######################################################################




######################################################################
######################################################################
InstallGlobalFunction(BettiNumbersOfMatrix,
function(AA)
local 
	A,B, i,j,
	Neighbours,
	rows,cols,firstBetti,
	secondBetti,box,
	colour,
	ChooseBox,
	neibs, newneibs,
	EdgeCount,
	VertexCount,
	Edges,
	Vertices,
	Faces,
	NeighbourValues,
	v;

A:=StructuralCopy(AA);
ContractMatrix(A);

rows:=Length(A);
cols:=Length(A[1]);
B:=StructuralCopy(A);

Edges:=0;
Vertices:=0;
Faces:=0;

#####################################################################
ChooseBox:=function(A)
local i,j, col, row;
#Returns a random (=first!) black box
row:=Length(A);
col:=Length(A[1]);
for i in [1..row] do
for j in [1..col] do
if A[i][j]=1 then return [i,j]; fi;
od;
od;
return fail;
end;
#####################################################################

#####################################################################
Neighbours:=function(B,box)
local nbs,i,j;
nbs:=[
[box[1],box[2]]-[1,1],
[box[1],box[2]]-[1,0],
[box[1],box[2]]-[1,-1],
[box[1],box[2]]-[0,1],
[box[1],box[2]]-[0,-1],
[box[1],box[2]]-[-1,1],
[box[1],box[2]]-[-1,0],
[box[1],box[2]]-[-1,-1]
];

if box[1]=1 or box[1]=rows or box[2]=1 or box[2]=cols then
nbs:=Filtered(nbs,x->x[1]>0 and x[1]<=rows and x[2]>0 and x[2]<=cols);
fi;
nbs:= Filtered(nbs,x->B[x[1]][x[2]]=1);
return nbs;
end;
####################################################################

#####################################################################
NeighbourValues:=function(A,box)
local nbrs,left, top,right,bottom,cols,rows,topleft,topright,bottomleft,bottomright;

cols:=Length(A[1]);
rows:=Length(A);

nbrs:=[];

if box[2]=1 then left:=0; else left:=A[box[1]][box[2]-1]; fi;
if box[2]=cols then right:=0; else right:=A[box[1]][box[2]+1]; fi;
if box[1]=1 then top:=0; else top:=A[box[1]-1][box[2]]; fi;
if box[1]=rows then bottom:=0; else bottom:=A[box[1]+1][box[2]]; fi;

if box[2]=1 or box[1]=1
then topleft:=0; else topleft:=A[box[1]-1][box[2]-1]; fi;

if box[2]=cols or box[1]=1 then topright:=0; else topright:=A[box[1]-1][box[2]+1]; fi;

if box[2]=1 or box[1]=rows then bottomleft:=0; else bottomleft:=A[box[1]+1][box[2]-1]; fi;

if box[2]=cols or box[1]=rows then bottomright:=0; else bottomright:=A[box[1]+1][box[2]+1]; fi;


return [left,topleft,top,topright,right,bottomright,bottom,bottomleft];

end;
######################################################################

####################################################################
EdgeCount:=function(v);
return 4-Sum(v{[1,3,5,7]})/2;
end;
####################################################################

####################################################################
VertexCount:=function(v);
return 
1/(1+v[1]+v[2]+v[3])
+1/(1+v[3]+v[4]+v[5])
+1/(1+v[5]+v[6]+v[7])
+1/(1+v[7]+v[8]+v[1]);
end;
####################################################################
box:=ChooseBox(B);
colour:=1;

while not box=fail do
   colour:=colour+1;
   B[box[1]][box[2]]:=colour;
    v:=NeighbourValues(A,box);
    Vertices:=Vertices+VertexCount(v);
    Edges:=Edges+EdgeCount(v);
    Faces:=Faces+1;

   neibs:=Neighbours(B,box);

	while Length(neibs)>0 do
	for box in neibs do
	B[box[1]][box[2]]:=colour;
	v:=NeighbourValues(A,box);
	Vertices:=Vertices+VertexCount(v);
	Edges:=Edges+EdgeCount(v);
	Faces:=Faces+1;
	od;
	newneibs:=Concatenation(List(neibs,x->Neighbours(B,x)));
	neibs:=SSortedList(newneibs);
	od;

   box:=ChooseBox(B);
od;

return  [colour-1, colour-Faces+Edges-Vertices-1];
end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(ThickenedMatrix,
function(arg)
local A,N,row,col,B,Neighbours,i,j,box;

A:=arg[1];
if Length(arg)>1 then N:=arg[2];
if N>1 then return ThickenedMatrix(ThickenedMatrix(A,N-1));fi;
fi;

B:=StructuralCopy(A);
row:=Length(B);
col:=Length(B[1]);

#####################################################################
Neighbours:=function(B,box)
local nbs,i,j;
nbs:=[
[box[1],box[2]]-[1,1],
[box[1],box[2]]-[1,0],
[box[1],box[2]]-[1,-1],
[box[1],box[2]]-[0,1],
[box[1],box[2]]-[0,-1],
[box[1],box[2]]-[-1,1],
[box[1],box[2]]-[-1,0],
[box[1],box[2]]-[-1,-1]
];

if box[1]=1 or box[1]=row or box[2]=1 or box[2]=col then
nbs:=Filtered(nbs,x->x[1]>0 and x[1]<=row and x[2]>0 and x[2]<=col);
fi;
return nbs;
end;
####################################################################

for i in [1..row] do
for j in [1..col] do
if A[i][j]=1 then
for box in Neighbours(A,[i,j]) do
B[box[1]][box[2]]:=1;
od;
fi;
od;od;

return B;
end);
######################################################################
######################################################################
