#########################################################
#########################################################
InstallGlobalFunction(HAP_Binlisttoint,
function(binlist)
local n,l;
n:=0;
for l in [1..26] do
	if binlist[l]=1 then
		n:=n+2^(l-1);
	fi;
od;
return n;
end);
#########################################################
#########################################################

HAP_MOVES_DIM_3:=0;

#########################################################
#########################################################
InstallGlobalFunction(HomotopyEquivalentSmallerSubArray3D,
function(AA,SS)
local
        A,S,
        tog,
        col,
        row,
        dep,
        i,j,ii,jj,k,kk,LTP,TP,RTP,LP,P,RP,LBP,BP,RBP,LT,T,
        RT,L,R,LB,B,RB,LTF,TF,RTF,LF,F,RF,LBF,BF,RBF, LST;

if HAP_MOVES_DIM_3=0 then
ReadPackage("HAP","lib/PolyComplexes/hapMovesDim3B.txt");
fi;

LST:=HAP_MOVES_DIM_3;
A:=FrameArray(AA);
S:=FrameArray(SS);
row:=Length(A);
col:=Length(A[1]);
dep:=Length(A[1][1]);
tog:=true;

######################################################
while tog do
tog:=false;

for i in [1..row] do
for j in [1..col] do
for k in [1..dep] do
if A[i][j][k]=1 and S[i][j][k]=0 then

LTP:=A[i-1][j-1][k-1];
TP:=A[i-1][j-1][k];
RTP:=A[i-1][j-1][k+1];
LP:=A[i-1][j][k-1];
P:=A[i-1][j][k];
RP:=A[i-1][j][k+1];
LBP:=A[i-1][j+1][k-1];
BP:=A[i-1][j+1][k];
RBP:=A[i-1][j+1][k+1];

LT:=A[i][j-1][k-1];
T:=A[i][j-1][k];
RT:=A[i][j-1][k+1];
L:=A[i][j][k-1];
R:=A[i][j][k+1];
LB:=A[i][j+1][k-1];
B:=A[i][j+1][k];
RB:=A[i][j+1][k+1];

LTF:=A[i+1][j-1][k-1];
TF:=A[i+1][j-1][k];
RTF:=A[i+1][j-1][k+1];
LF:=A[i+1][j][k-1];
F:=A[i+1][j][k];
RF:=A[i+1][j][k+1];
LBF:=A[i+1][j+1][k-1];
BF:=A[i+1][j+1][k];
RBF:=A[i+1][j+1][k+1];

if  LST[HAP_Binlisttoint([LTP,TP,RTP,LP,P,RP,LBP,BP,RBP,LT,T,
RT,L,R,LB,B,RB,LTF,TF,RTF,LF,F,RF,LBF,BF,RBF])+1]
then A[i][j][k]:=0; tog:=true;
fi;

fi;
od;
od;
od;

if tog then
for ii in [1..row] do
for jj in [1..col] do
for kk in [1..dep] do
i:=row+1-ii;j:=col+1-jj;k:=dep+1-kk;
if A[i][j][k]=1 and S[i][j][k]=0 then

LTP:=A[i-1][j-1][k-1];
TP:=A[i][j-1][k-1];
RTP:=A[i+1][j-1][k-1];
LP:=A[i-1][j][k-1];
P:=A[i][j][k-1];
RP:=A[i+1][j][k-1];
LBP:=A[i-1][j+1][k-1];
BP:=A[i][j+1][k-1];
RBP:=A[i+1][j+1][k-1];

LT:=A[i-1][j-1][k];
T:=A[i][j-1][k];
RT:=A[i+1][j-1][k];
L:=A[i-1][j][k];
R:=A[i+1][j][k];
LB:=A[i-1][j+1][k];
B:=A[i][j+1][k];
RB:=A[i+1][j+1][k];

LTF:=A[i-1][j-1][k+1];
TF:=A[i][j-1][k+1];
RTF:=A[i+1][j-1][k+1];
LF:=A[i-1][j][k+1];
F:=A[i][j][k+1];
RF:=A[i+1][j][k+1];
LBF:=A[i-1][j+1][k+1];
BF:=A[i][j+1][k+1];
RBF:=A[i+1][j+1][k+1];

if  LST[HAP_Binlisttoint([LTP,TP,RTP,LP,P,RP,LBP,BP,RBP,LT,T,RT,L,R,LB,B,
RB,LTF,TF,RTF,LF,F,RF,LBF,BF,RBF])+1]
then A[i][j][k]:=0; tog:=true;
fi;

fi;
od;
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
InstallGlobalFunction(HomotopyEquivalentLargerSubArray3D,
function(AA,BB)
local
        A,B,
        tog,
        col,
        row,
        dep,
        start,
        i,j,k,ii,jj,kk,
        LTP,TP,RTP,LP,RP,P,LBP,BP,RBP,T,RT,L,R,
        LB,RB,LTF,TF,RTF,LF,F,RF,LBF,BF,RBF,LT,Bk,
        LST ;

if HAP_MOVES_DIM_3=0 then
ReadPackage("HAP","lib/PolyComplexes/hapMovesDim3B.txt");
fi;
LST:=HAP_MOVES_DIM_3;

A:=FrameArray(AA);
B:=FrameArray(BB);
row:=Length(A);
col:=Length(A[1]);
dep:=Length(A[1][1]);

##If B is empty then###############################
if ArraySum(B)=0 then

start:=0;
############################
for k in [1..dep] do
for j in [1..col] do
for i in [1..row] do
if A[i][j][k]=1 then
start:=[i,j,k]; break; fi;
od;
if IsList(start) then break; fi;
od;
if IsList(start) then break; fi;
od;
############################

#B:=A*0;
B[start[1]][start[2]][start[3]]:=1;

fi;
##B is now probably non-empty######################

tog:=true;
while tog do
tog:=false;

for i in [1..row] do
for j in [1..col] do
for k in [1..dep] do
if A[i][j][k]=1 and B[i][j][k]=0 then

LTP:=B[i-1][j-1][k-1];
TP:=B[i-1][j-1][k];
RTP:=B[i-1][j-1][k+1];
LP:=B[i-1][j][k-1];
P:=B[i-1][j][k];
RP:=B[i-1][j][k+1];
LBP:=B[i-1][j+1][k-1];
BP:=B[i-1][j+1][k];
RBP:=B[i-1][j+1][k+1];

LT:=B[i][j-1][k-1];
T:=B[i][j-1][k];
RT:=B[i][j-1][k+1];
L:=B[i][j][k-1];
R:=B[i][j][k+1];
LB:=B[i][j+1][k-1];
Bk:=B[i][j+1][k];
RB:=B[i][j+1][k+1];

LTF:=B[i+1][j-1][k-1];
TF:=B[i+1][j-1][k];
RTF:=B[i+1][j-1][k+1];
LF:=B[i+1][j][k-1];
F:=B[i+1][j][k];
RF:=B[i+1][j][k+1];
LBF:=B[i+1][j+1][k-1];
BF:=B[i+1][j+1][k];
RBF:=B[i+1][j+1][k+1];

if  LST[HAP_Binlisttoint([LTP,TP,RTP,LP,P,RP,LBP,BP,RBP,LT,T,
RT,L,R,LB,Bk,RB,LTF,TF,RTF,LF,F,RF,LBF,BF,RBF])+1]
then B[i][j][k]:=1; tog:=true;
fi;

fi;
od;
od;
od;

if tog then
for ii in [1..row] do
i:=row-ii+1;

for jj in [1..col] do
j:=col-jj+1;

for kk in [1..dep] do
k:=dep-kk+1;



if A[i][j][k]=1 and B[i][j][k]=0 then

LTP:=B[i-1][j-1][k-1];
TP:=B[i-1][j-1][k];
RTP:=B[i-1][j-1][k+1];
LP:=B[i-1][j][k-1];
P:=B[i-1][j][k];
RP:=B[i-1][j][k+1];
LBP:=B[i-1][j+1][k-1];
BP:=B[i-1][j+1][k];
RBP:=B[i-1][j+1][k+1];

LT:=B[i][j-1][k-1];
T:=B[i][j-1][k];
RT:=B[i][j-1][k+1];
L:=B[i][j][k-1];
R:=B[i][j][k+1];
LB:=B[i][j+1][k-1];
Bk:=B[i][j+1][k];
RB:=B[i][j+1][k+1];

LTF:=B[i+1][j-1][k-1];
TF:=B[i+1][j-1][k];
RTF:=B[i+1][j-1][k+1];
LF:=B[i+1][j][k-1];
F:=B[i+1][j][k];
RF:=B[i+1][j][k+1];
LBF:=B[i+1][j+1][k-1];
BF:=B[i+1][j+1][k];
RBF:=B[i+1][j+1][k+1];

if  LST[HAP_Binlisttoint([LTP,TP,RTP,LP,P,RP,LBP,BP,RBP,LT,T,RT,L,R,LB,Bk,
RB,LTF,TF,RTF,LF,F,RF,LBF,BF,RBF])+1]
then B[i][j][k]:=1; tog:=true;
fi;

fi;
od;
od;
od;
fi;

od;

return UnframeArray(B);
end);
#####################################################################
#####################################################################

