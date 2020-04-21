
HAP_PERMMOVES_DIM_2:=0;

#####################################################################
#####################################################################
InstallGlobalFunction(HomotopyEquivalentSmallerSubPermMatrix,
function(AA,SS)
local
        A,S,
	tog,
        col,
        row,
        i,j,ii,jj,
        left,top,right,bottom,topright,bottomleft;

if HAP_PERMMOVES_DIM_2=0 then 
ReadPackage("HAP","lib/PolyComplexes/hapPermMovesDim2B.txt"); 
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
topright:=A[i-1][j+1];
bottomleft:=A[i+1][j-1];
if 
HAP_PERMMOVES_DIM_2[1+left+2*top+4*topright+8*right+16*bottom+32*bottomleft]
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
topright:=A[i-1][j+1];
bottomleft:=A[i+1][j-1];
if 
HAP_PERMMOVES_DIM_2[1+left+2*top+4*topright+8*right+16*bottom+32*bottomleft]
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



######################################################################
######################################################################
InstallGlobalFunction(ContractPermMatrix,
function(A);
return HomotopyEquivalentSmallerSubPermMatrix(A,A*0);
end);
######################################################################
######################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(HomotopyEquivalentLargerSubPermMatrix,
function(AA,BB)
local
        A,B,
        tog,
        col,
        row,
        start,
        i,j,ii,jj,
        left,top,right,bottom,topleft,topright,bottomleft,bottomright;

if HAP_PERMMOVES_DIM_2=0 then
ReadPackage("HAP","lib/PolyComplexes/hapPermMovesDim2B.txt");
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
topright:=B[i-1][j+1];
bottomleft:=B[i+1][j-1];
if
HAP_PERMMOVES_DIM_2[1+left+2*top+4*topright+8*right+16*bottom+32*bottomleft]
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
topright:=B[i-1][j+1];
bottomleft:=B[i+1][j-1];
if
HAP_PERMMOVES_DIM_2[1+left+2*top+4*topright+8*right+16*bottom+32*bottomleft]
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






