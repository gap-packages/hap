#C 2007 Graham Ellis

#####################################################################
#####################################################################
InstallGlobalFunction(DeleteAndContract,
function(A,box)
local
	Neighbours,
	NeighboursValues,
	MOVES,
	Contract,
	ContractImageMatrix,
	TOG1,
	TOG2,
	NeighboursOfBox,
	rows,cols;

A[box[1]][box[2]]:=0;
rows:=Length(A);
cols:=Length(A[1]);

#####################################################################
Neighbours:=function(B,box,rows,cols)
local nbs,i,j;
nbs:=[
[box[1],box[2]]-[1,0],
[box[1],box[2]]-[0,1],
[box[1],box[2]]-[0,-1],
[box[1],box[2]]-[-1,0],

#[box[1],box[2]]-[1,1],
#[box[1],box[2]]-[1,-1],
#[box[1],box[2]]-[-1,1],
#[box[1],box[2]]-[-1,-1]

];

if box[1]=1 or box[1]=rows or box[2]=1 or box[2]=cols then
nbs:=Filtered(nbs,x->x[1]>0 and x[1]<=rows and x[2]>0 and x[2]<=cols);
fi;
nbs:= Filtered(nbs,x->B[x[1]][x[2]]=1);
return nbs;
end;
####################################################################


#####################################################################
NeighboursValues:=function(A,box)
local nbrs,left, top,right,bottom,cols,rows,topleft,topright,bottomleft,bottomright;

cols:=Length(A[1]);
rows:=Length(A);

nbrs:=[];

if box[2]=1 then left:=0; else left:=StructuralCopy(A[box[1]][box[2]-1]); fi;
if box[2]=cols then right:=0; else right:=StructuralCopy(A[box[1]][box[2]+1]); fi;
if box[1]=1 then top:=0; else top:=StructuralCopy(A[box[1]-1][box[2]]); fi;
if box[1]=rows then bottom:=0; else bottom:=StructuralCopy(A[box[1]+1][box[2]]); fi;

if box[2]=1 or box[1]=1 
then topleft:=0; else topleft:=StructuralCopy(A[box[1]-1][box[2]-1]); fi;

if box[2]=cols or box[1]=1 then topright:=0; else topright:=
StructuralCopy(A[box[1]-1][box[2]+1]); fi;

if box[2]=1 or box[1]=rows then bottomleft:=0; else bottomleft:=
StructuralCopy(A[box[1]+1][box[2]-1]); fi;

if box[2]=cols or box[1]=rows then bottomright:=0; else bottomright:=
StructuralCopy(A[box[1]+1][box[2]+1]); fi;


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
Contract:=function(A,box,Moves,NeighboursOfBox)
local tmp;

if NeighboursValues(A,box) in Moves
then A[box[1]][box[2]]:=0; 
tmp:=Neighbours(A,box,rows,cols);
tmp:=Filtered(tmp,b->NeighboursValues(A,b) in Moves);
Append(NeighboursOfBox,tmp);
return true;else 
return false; fi;
end;
######################################################################

######################################################################
ContractImageMatrix:=function(A,Moves,NeighboursOfBox)
local i,j,row,col,tog,togII,ii,jj,box;
row:=Length(A);
col:=Length(A[1]);
tog:=false;

for box in NeighboursOfBox do

if A[box[1]][box[2]]=1 then
tog:=Contract(A,box,Moves,NeighboursOfBox) or tog; 
fi;
od;

return tog;
end;
######################################################################

NeighboursOfBox:=Neighbours(A,box,rows,cols);
TOG2:=false;
TOG1:=true;
while TOG1 do
TOG1:= ContractImageMatrix(A,MOVES,NeighboursOfBox);
if TOG1 then TOG2:=true;fi;
NeighboursOfBox:=Filtered(NeighboursOfBox,x->x=1);
od;

return TOG2;

end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(BoundaryOfMatrix,
function(AA)
local
	NeighboursValues,
	A,
	i,j,
	rows,cols,
	TOG;
A:=StructuralCopy(AA);
rows:=Length(A);
cols:=Length(A[1]);

#####################################################################
NeighboursValues:=function(A,box)
local nbrs,left, top,right,bottom,cols,rows,topleft,topright,bottomleft,bottomright;

cols:=Length(A[1]);
rows:=Length(A);

nbrs:=[];

if box[2]=1 then left:=0; else left:=StructuralCopy(A[box[1]][box[2]-1]); fi;
if box[2]=cols then right:=0; else right:=StructuralCopy(A[box[1]][box[2]+1]); fi;
if box[1]=1 then top:=0; else top:=StructuralCopy(A[box[1]-1][box[2]]); fi;
if box[1]=rows then bottom:=0; else bottom:=StructuralCopy(A[box[1]+1][box[2]]); fi;

if box[2]=1 or box[1]=1
then topleft:=0; else topleft:=StructuralCopy(A[box[1]-1][box[2]-1]); fi;

if box[2]=cols or box[1]=1 then topright:=0; else topright:=
StructuralCopy(A[box[1]-1][box[2]+1]); fi;

if box[2]=1 or box[1]=rows then bottomleft:=0; else bottomleft:=
StructuralCopy(A[box[1]+1][box[2]-1]); fi;

if box[2]=cols or box[1]=rows then bottomright:=0; else bottomright:=
StructuralCopy(A[box[1]+1][box[2]+1]); fi;


return [left,topleft,top,topright,right,bottomright,bottom,bottomleft];
end;
#####################################################################

TOG:=true;
while TOG do
TOG:=false;
	for i in [1..rows] do
	for j in [1..cols] do
		if A[i][j]=1 then
		if Sum(NeighboursValues(A,[i,j]))=8 then
		TOG:=true;
		DeleteAndContract(A,[i,j]);
		fi;fi;
	od;od;
od;

TOG:=true;
while TOG do
TOG:=false;
        for i in [1..rows] do
        for j in [1..cols] do
                if A[i][j]=1 then
                if Sum(NeighboursValues(AA,[i,j]))=8 then
                TOG:=true;
		DeleteAndContract(A,[i,j]);
		fi;fi;
        od;od;
od;

return A;

end);
######################################################################
######################################################################

######################################################################
######################################################################
InstallGlobalFunction(SingularityMatrix,
function(AA)		#More experimentation is needed to
			#find better values of r and s.
local
	SingularPoint,
	A,B,
	C,
	r,s,t,
	rows,
	cols,
	i,j,b,vv;

r:=5;;
s:=2;;
t:=20;;

###########################################################
SingularPoint:=function(A,box,r)
local i,x,y,v,vecs;
vecs:=[];

for x in [-r..r] do
if A[box[1]-r][box[2]+x]=1 then
for i in [-s..s] do
if	 A[box[1]+r+i][box[2]-x]=1 then return false; fi;
if	 A[box[1]+r][box[2]-x+i]=1 then return false; fi;
od;
fi;
od;

for y in [-r..r] do
if A[box[1]+y][box[2]+r]=1 then
for i in [-s..s] do
if	A[box[1]-y+i][box[2]-r]=1 then return false; fi;
if	A[box[1]-y][box[2]-r+i]=1 then return false; fi;
od;
fi;
od;

return true;
end;
###########################################################

A:=[]; #A will be be equal to AA with a border of zeros of width t.
vv:=ListWithIdenticalEntries(Length(AA[1])+2*t,0);
for i in [1..t] do
A[i]:=vv;
A[t+i+Length(AA)]:=vv;
od;
for i in [1..Length(AA)] do
A[t+i]:=Concatenation(List([1..t],a->0),AA[i],List([1..t],a->0));
od;

B:=BoundaryOfMatrix(A);

C:=0*B;

rows:=Length(B);
cols:=Length(B[1]);
for i in [1..rows] do
for j in [1..cols] do
if B[i][j]=1 then
if SingularPoint(B,[i,j],r)then C[i][j]:=1; fi;
fi;
od;od;

C:=ThickenedMatrix(C,4);;

return C;

end);
######################################################################
######################################################################
