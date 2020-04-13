#########################################################
#########################################################
InstallGlobalFunction(HAP_PermBinlisttoint,
function(binlist)
local n,l;
n:=1;
for l in [1..14] do
        if binlist[l]=1 then
                n:=n+2^(l-1);
        fi;
od;
return n;
end);
#########################################################
#########################################################


HAP_PERMMOVES_DIM_3:=0;

#########################################################
#########################################################
InstallGlobalFunction(HomotopyEquivalentSmallerSubPermArray3D,
function(AA,SS)
local
        A,S,
        tog,
        col,
        row,
        dep,
        i,j,ii,jj,k,kk,LTP,TP,RTP,LP,P,RP,LBP,BP,RBP,LT,T,
        RT,L,R,LB,B,RB,LTF,TF,RTF,LF,F,RF,LBF,BF,RBF, LST;

if HAP_PERMMOVES_DIM_3=0 then
ReadPackage("HAP","lib/PolyComplexes/hapPermMovesDim3B.txt");
fi;

LST:=HAP_PERMMOVES_DIM_3;
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

RTP:=A[i-1][j][k];
P:=A[i-1][j][k+1];
RP:=A[i-1][j+1][k];
BP:=A[i-1][j+1][k+1];

T:=A[i][j-1][k];
RT:=A[i][j-1][k+1];
L:=A[i][j][k-1];
R:=A[i][j][k+1];
LB:=A[i][j+1][k-1];
B:=A[i][j+1][k];

TF:=A[i+1][j-1][k-1];
LF:=A[i+1][j-1][k];
F:=A[i+1][j][k-1];
LBF:=A[i+1][j][k];

if  LST[HAP_PermBinlisttoint([RTP,P,RP,BP,T,RT,L,R,LB,B,TF,LF,F,LBF])]
then 
A[i][j][k]:=0; tog:=true;

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

RTP:=A[i-1][j][k];
P:=A[i-1][j][k+1];
RP:=A[i-1][j+1][k];
BP:=A[i-1][j+1][k+1];

T:=A[i][j-1][k];
RT:=A[i][j-1][k+1];
L:=A[i][j][k-1];
R:=A[i][j][k+1];
LB:=A[i][j+1][k-1];
B:=A[i][j+1][k];

TF:=A[i+1][j-1][k-1];
LF:=A[i+1][j-1][k];
F:=A[i+1][j][k-1];
LBF:=A[i+1][j][k];

if LST[HAP_PermBinlisttoint([RTP,P,RP,BP,T,RT,L,R,LB,B,TF,LF,F,LBF])] 
then 

A[i][j][k]:=0; tog:=true;

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
InstallGlobalFunction(HomotopyEquivalentLargerSubPermArray3D,
function(AA,BB)
local
        A,B,
        tog,
        col,
        row,
        dep,
        start,
        i,j,k,ii,jj,kk,
        RTP,PX,RP,BP,TX,RT,LX,RX,LB,BX,TF,LF,FX,LBF,
        Bk, LST ;

if HAP_PERMMOVES_DIM_3=0 then
ReadPackage("HAP","lib/PolyComplexes/hapPermMovesDim3B.txt");
fi;
LST:=HAP_PERMMOVES_DIM_3;

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

RTP:=B[i-1][j][k];
PX:=B[i-1][j][k+1];
RP:=B[i-1][j+1][k];
BP:=B[i-1][j+1][k+1];

TX:=B[i][j-1][k];
RT:=B[i][j-1][k+1];
LX:=B[i][j][k-1];
RX:=B[i][j][k+1];
LB:=B[i][j+1][k-1];
BX:=B[i][j+1][k];

TF:=B[i+1][j-1][k-1];
LF:=B[i+1][j-1][k];
FX:=B[i+1][j][k-1];
LBF:=B[i+1][j][k];

if LST[HAP_PermBinlisttoint([RTP,PX,RP,BP,TX,RT,LX,RX,LB,BX,TF,LF,FX,LBF])]
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

RTP:=B[i-1][j][k];
PX:=B[i-1][j][k+1];
RP:=B[i-1][j+1][k];
BP:=B[i-1][j+1][k+1];

TX:=B[i][j-1][k];
RT:=B[i][j-1][k+1];
LX:=B[i][j][k-1];
RX:=B[i][j][k+1];
LB:=B[i][j+1][k-1];
BX:=B[i][j+1][k];

TF:=B[i+1][j-1][k-1];
LF:=B[i+1][j-1][k];
FX:=B[i+1][j][k-1];
LBF:=B[i+1][j][k];

if LST[HAP_PermBinlisttoint([RTP,PX,RP,BP,TX,RT,LX,RX,LB,BX,TF,LF,FX,LBF])]
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

