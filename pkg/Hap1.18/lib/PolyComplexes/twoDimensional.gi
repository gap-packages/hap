######################################################################
######################################################################
InstallGlobalFunction(BettinumbersOfPureCubicalComplex_dim_2,
function(AA)
local 
	A,B,i,j,
	Neighbours,
	rows,cols,firstBetti,
	box,
	colour,
	ChooseBox,
	neibs, newneibs,
	EdgeCount,
	VertexCount,
	Edges,
	Vertices,
	Faces,
	NeighbourValues,
	v,vv;

A:=FrameArray(AA!.binaryArray);
#if EvaluateProperty(AA,"contracted")=fail then
#ContractMatrix(A);
#fi;
B:=StructuralCopy(A);

rows:=Length(A);
cols:=Length(A[1]);

Edges:=0;
Vertices:=0;
Faces:=0;

#####################################################################
ChooseBox:=function(A)
local i,j, col, row;
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
local nbs;
nbs:=[
[box[1]-1,box[2]-1],
[box[1]-1,box[2]],
[box[1]-1,box[2]+1],
[box[1],box[2]-1],
[box[1],box[2]+1],
[box[1]+1,box[2]-1],
[box[1]+1,box[2]],
[box[1]+1,box[2]+1]
];

return Filtered(nbs,x->B[x[1]][x[2]]=1);
end;
####################################################################

#####################################################################
NeighbourValues:=function(A,box)
local left,top,right,bottom,topleft,topright,bottomleft,bottomright;

left:=A[box[1]][box[2]-1]; 
right:=A[box[1]][box[2]+1]; 
top:=A[box[1]-1][box[2]]; 
bottom:=A[box[1]+1][box[2]]; 
topleft:=A[box[1]-1][box[2]-1]; 
topright:=A[box[1]-1][box[2]+1]; 
bottomleft:=A[box[1]+1][box[2]-1]; 
bottomright:=A[box[1]+1][box[2]+1]; 

return [left,topleft,top,topright,right,bottomright,bottom,bottomleft];
end;
######################################################################

####################################################################
box:=ChooseBox(B);
colour:=1;

while not box=fail do
    colour:=colour+1;
    B[box[1]][box[2]]:=colour;
    vv:=NeighbourValues(A,box);
    Vertices:=Vertices+
    1/(1+vv[1]+vv[2]+vv[3]) +1/(1+vv[3]+vv[4]+vv[5])
    +1/(1+vv[5]+vv[6]+vv[7]) +1/(1+vv[7]+vv[8]+vv[1]);
    Edges:=Edges+4-Sum(vv{[1,3,5,7]})/2;
    Faces:=Faces+1;

    neibs:=Neighbours(B,box);

	while Length(neibs)>0 do
	for box in neibs do
	B[box[1]][box[2]]:=colour;
	v:=NeighbourValues(A,box);
	Vertices:=Vertices+
	1/(1+v[1]+v[2]+v[3]) +1/(1+v[3]+v[4]+v[5]) 
	+1/(1+v[5]+v[6]+v[7]) +1/(1+v[7]+v[8]+v[1]);
	Edges:=Edges+4-Sum(v{[1,3,5,7]})/2;
	Faces:=Faces+1;
	od;
	newneibs:=Concatenation(List(neibs,x->Neighbours(B,x)));
	neibs:=SSortedList(newneibs);
	od;

    box:=ChooseBox(B);
od;

return  [colour-1, colour-Faces+Edges-Vertices-1, 0];
end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(ContractMatrix,
function(A);
return HomotopyEquivalentSmallerSubMatrix(A,A*0);
end);
######################################################################
######################################################################

HAP_MOVES_DIM_2:=0;

#####################################################################
#####################################################################
InstallGlobalFunction(HomotopyEquivalentSmallerSubMatrix,
function(AA,SS)
local
        A,S,
	tog,
        col,
        row,
        i,j,ii,jj,
        left,top,right,bottom,topleft,topright,bottomleft,bottomright;

if HAP_MOVES_DIM_2=0 then 
ReadPackage("HAP","lib/PolyComplexes/hapMovesDim2B.txt"); 
fi;

A:=FrameArray(AA);
S:=FrameArray(SS);
row:=Length(A);
col:=Length(A[1]);

tog:=true;
######################################################
while tog do
tog:=false;

for i in [1..row] do
for j in [1..col] do
if A[i][j]=1 and S[i][j]=0 then
left:=A[i][j-1];
right:=A[i][j+1];
top:=A[i-1][j];
bottom:=A[i+1][j];
topleft:=A[i-1][j-1];
topright:=A[i-1][j+1];
bottomleft:=A[i+1][j-1];
bottomright:=A[i+1][j+1];
if 
HAP_MOVES_DIM_2[1+left+2*topleft+4*top+8*topright+16*right+32*bottomright+64*bottom+128*bottomleft]
then A[i][j]:=0; tog:=true;
fi;
fi;
od;
od;

if tog then
for ii in [1..row] do
for jj in [1..col] do
i:=row+1-ii;j:=col+1-jj;
if A[i][j]=1 and S[i][j]=0 then 
left:=A[i][j-1];
right:=A[i][j+1];
top:=A[i-1][j];
bottom:=A[i+1][j];
topleft:=A[i-1][j-1];
topright:=A[i-1][j+1];
bottomleft:=A[i+1][j-1];
bottomright:=A[i+1][j+1];
if 
HAP_MOVES_DIM_2[1+left+2*topleft+4*top+8*topright+16*right+32*bottomright+64*bottom+128*bottomleft]
then A[i][j]:=0; tog:=true;
fi;
fi;
od;
od;
fi;

od;
#############################################

return UnframeArray(A);
end);
######################################################################
######################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(ContractibleSubMatrix,
function(AA);
return HomotopyEquivalentLargerSubMatrix(AA,AA*0);
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(HomotopyEquivalentLargerSubMatrix,
function(AA,BB)
local
	A,B,
        tog,
        col,
        row,
	start,
        i,j,ii,jj,
        left,top,right,bottom,topleft,topright,bottomleft,bottomright;

if HAP_MOVES_DIM_2=0 then 
ReadPackage("HAP","lib/PolyComplexes/hapMovesDim2B.txt"); 
fi;

A:=FrameArray(AA);
B:=FrameArray(BB);
row:=Length(A);
col:=Length(A[1]);

##If B is empty then###############################
if ArraySum(B)=0 then

start:=0;

############################
for i in [1..row] do
for j in [1..col] do
if A[i][j]=1 then
start:=[i,j]; break; fi;
od;
if IsList(start) then break; fi;
od;
############################

#B:=A*0;
B[start[1]][start[2]]:=1;

fi;
##B is now probably non-empty######################

tog:=true;
while tog do
tog:=false;

for i in [1..row] do
for j in [1..col] do
if A[i][j]=1 and B[i][j]=0 then
left:=B[i][j-1];
right:=B[i][j+1];
top:=B[i-1][j];
bottom:=B[i+1][j];
topleft:=B[i-1][j-1];
topright:=B[i-1][j+1];
bottomleft:=B[i+1][j-1];
bottomright:=B[i+1][j+1];
if 
HAP_MOVES_DIM_2[1+left+2*topleft+4*top+8*topright+16*right+32*bottomright+64*bottom+128*bottomleft]
then B[i][j]:=1; tog:=true;
fi;
fi;
od;
od;

if tog then
for ii in [1..row] do
i:=row-ii+1;
for jj in [1..col] do
j:=col-jj+1;
if A[i][j]=1 and B[i][j]=0 then
left:=B[i][j-1];
right:=B[i][j+1];
top:=B[i-1][j];
bottom:=B[i+1][j];
topleft:=B[i-1][j-1];
topright:=B[i-1][j+1];
bottomleft:=B[i+1][j-1];
bottomright:=B[i+1][j+1];
if 
HAP_MOVES_DIM_2[1+left+2*topleft+4*top+8*topright+16*right+32*bottomright+64*bottom+128*bottomleft]
then B[i][j]:=1; tog:=true;
fi;
fi;
od;
od;
fi;

od;

return UnframeArray(B);
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(ThickenedPureCubicalComplex_dim2,
function(M)
local
	A,B,rows,cols,ball,x,i,j;

A:=FrameArray(M!.binaryArray);
B:=A*0;
rows:=Length(B);
cols:=Length(B[1]);
ball:=Cartesian([-1,0,1],[-1,0,1]);
#RemoveSet(ball,[0,0]);

for i in [2..rows-1] do
for j in [2..cols-1] do
if A[i][j]=1 then
  for x in ball do
  B[i+x[1]][j+x[2]]:=1;
  od;
fi;
od;
od;

B:=UnframeArray(B);
return PureCubicalComplex(B);

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(ThickenedHEPureCubicalComplex,
function(M)
local
        A,B,rows,cols,ball,x,i,j,left,right,top,bottom,topleft,topright,bottomleft, bottomright;

A:=FrameArray(M!.binaryArray);
A:=FrameArray(A);
B:=A*1;
rows:=Length(B);
cols:=Length(B[1]);
ball:=Cartesian([-1,0,1],[-1,0,1]);
#RemoveSet(ball,[0,0]);

if HAP_MOVES_DIM_2=0 then
ReadPackage("HAP","lib/PolyComplexes/hapMovesDim2B.txt");
fi;


for i in [2..rows-1] do
for j in [2..cols-1] do
if A[i][j]=1 then
  for x in ball do
  if B[i+x[1]][j+x[2]]=0 then
###
left:=B[i+x[1]][j-1+x[2]];
right:=B[i+x[1]][j+1+x[2]];
top:=B[i-1+x[1]][j+x[2]];
bottom:=B[i+1+x[1]][j+x[2]];
topleft:=B[i-1+x[1]][j-1+x[2]];
topright:=B[i-1+x[1]][j+1+x[2]];
bottomleft:=B[i+1+x[1]][j-1+x[2]];
bottomright:=B[i+1+x[1]][j+1+x[2]];
    if
    HAP_MOVES_DIM_2[1+left+2*topleft+4*top+8*topright+16*right+32*bottomright+64*bottom+128*bottomleft]
    then; B[i+x[1]][j+x[2]]:=1;
    fi;
###
  fi;
  od;
fi;
od;
od;

B:=UnframeArray(B);
B:=UnframeArray(B);
return PureCubicalComplex(B);

end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(ExcisedPureCubicalPair_dim_2,
function(M,S)
local
        A,B,C,D,rows,cols,ball,x,i,j,tog;

A:=FrameArray(M!.binaryArray);
B:=A*0;D:=A*0;
C:=FrameArray(S!.binaryArray);
rows:=Length(B);
cols:=Length(B[1]);
ball:=Cartesian([-1,0,1],[-1,0,1]);

for i in [2..rows-1] do
for j in [2..cols-1] do
if C[i][j]=1 then
tog:=false;
  for x in ball do
  if C[i+x[1]][j+x[2]]=0 then tog:=true; break;fi;
  od;
if tog then B[i][j]:=1; fi;
fi;
od;
od;

for i in [2..rows-1] do
for j in [2..cols-1] do
if A[i][j]=1 and ( C[i][j]=0 or B[i][j]=1) then D[i][j]:=1; fi;  
od;
od;

B:=UnframeArray(B);
D:=UnframeArray(D);
return [PureCubicalComplex(D), PureCubicalComplex(B)];

end);
#####################################################################
#####################################################################





