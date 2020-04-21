
################################################################
################################################################
InstallGlobalFunction(HopfSatohSurface,
function()
local S,t,z,hbar,hbar2,vbar, a,b,c,d,A,B,C,D,circ,v,w,aa;


S:=[];;
t:=4;
z:=5;

#####
hbar:=function(V,L,t)
local  circ, i,c,x;
circ:=[[0,1,0],[0,-1,0],[0,0,1],[0,0,-1],[0,1,1],[0,-1,-1],[0,1,-1],[0,-1,1]];
for i in [0..L] do
for c in circ do
x:=V+[i,0,0]+c;
Add(S,[x[1],x[2],x[3],t]);
od;
od;
Add(S,[V[1],V[2],V[3],t]);
Add(S,[V[1]+L,V[2],V[3],t]);
S:=SSortedList(S);
end;
#####

#####
hbar2:=function(V,L,M,t,s)
local  circ, i,c,x;
circ:=[[0,1,0],[0,-1,0],[0,0,1],[0,0,-1],[0,1,1],[0,-1,-1],[0,1,-1],[0,-1,1]];
for i in [0..L] do
for c in circ do
x:=V+[i,0,0]+c;
Add(S,[x[1],x[2],x[3],t]);
od;
od;
for i in [L..M] do
for c in circ do
x:=V+[i,0,0]+c;
Add(S,[x[1],x[2],x[3],s]);
od;
od;
Add(S,[V[1],V[2],V[3],t]);
Add(S,[V[1]+M,V[2],V[3],s]);
S:=SSortedList(S);
end;
#####


#####
vbar:=function(V,L,t)
local  circ, j,c,i,x;
circ:=[
[-3,0,-3],[-2,0,-3],[-1,0,-3],[0,0,-3],[1,0,-3],[2,0,-3],[3,0,-3],
[-3,0,-2],[3,0,-2],
[-3,0,-1],[3,0,-1],
[-3,0,-0],[3,0,-0],
[-3,0,1],[3,0,1],
[-3,0,2],[3,0,2],
[-3,0,3],[-2,0,3],[-1,0,3],[0,0,3],[1,0,3],[2,0,3],[3,0,3]
];
for j in [0..L] do
for c in circ do
x:=V+[0,j,0]+c;
Add(S,[x[1],x[2],x[3],t]);
od;
od;
for i in [-3..3] do
for j in [-3..3] do
if not [i,j]=[0,0] then
x:=V+[i,0,j];
Add(S,[x[1],x[2],x[3],t]);
x:=V+[i,L,j];
Add(S,[x[1],x[2],x[3],t]);
fi;
od;od;
S:=SSortedList(S);
end;
#####

############# A
#           #
#     ##### # ##### B
#     #     #     #
##### # #####     # C
      #           #
      ############# D
#a     b     c     d    
a:=5;
b:=14;
c:=28+1;
d:=38;
A:=5;
B:=10;
C:=15;
D:=20;

hbar([a,A,z],c-a,t);
RemoveSet(S,[a+1,A+1,z,t]);
RemoveSet(S,[c-1,A+1,z,t]);
hbar2([b,B,z],c-b-1,d-b,t-2,t+2);
RemoveSet(S,[b+1,B+1,z,t-2]);
RemoveSet(S,[d-1,B+1,z,t+2]);
#####
circ:=[[1,0,0],[-1,0,0],[0,0,1],[0,0,-1],[1,0,1],[-1,0,-1],[1,0,-1],[-1,0,1]];
for v in circ do
w:=v+[b+1,B+1,z];
w:=[w[1],w[2],w[3],t-1];
AddSet(S,w);
w:=[w[1],w[2],w[3],t-2];
AddSet(S,w);
w:=[w[1],w[2],w[3],t];
AddSet(S,w);



w:=v+[d-1,B+1,z];
w:=[w[1],w[2],w[3],t+1];
AddSet(S,w);
w:=[w[1],w[2],w[3],t+2];
AddSet(S,w);
w:=[w[1],w[2],w[3],t];
AddSet(S,w);

od;

circ:=[[0,1,0],[0,-1,0],[0,0,1],[0,0,-1],[0,1,1],[0,-1,-1],[0,1,-1],[0,-1,1]];

for v in circ do


w:=v+[c-1,B,z];
w:=[w[1],w[2],w[3],t-1];
AddSet(S,w);
w:=[w[1],w[2],w[3],t];
AddSet(S,w);
w:=[w[1],w[2],w[3],t+1];
AddSet(S,w);
w:=[w[1],w[2],w[3],t+2];
AddSet(S,w);
w:=[w[1],w[2],w[3],t-2];
AddSet(S,w);







od;
#####

hbar([a,C,z],c-a,t-2);
RemoveSet(S,[a+1,C-1,z,t-2]);
RemoveSet(S,[c-1,C-1,z,t-2]);
#####
circ:=[[1,0,0],[-1,0,0],[0,0,1],[0,0,-1],[1,0,1],[-1,0,-1],[1,0,-1],[-1,0,1]];

for v in circ do
w:=v+[a+1,C-1,z];
w:=[w[1],w[2],w[3],t-1];
AddSet(S,w);
w:=[w[1],w[2],w[3],t-2];
AddSet(S,w);
w:=[w[1],w[2],w[3],t];
AddSet(S,w);


w:=v+[c-1,C-1,z];
w:=[w[1],w[2],w[3],t-1];
AddSet(S,w);
w:=[w[1],w[2],w[3],t-2];
AddSet(S,w);
w:=[w[1],w[2],w[3],t];
AddSet(S,w);


od;
#####
hbar([b,D,z],d-b,t);
RemoveSet(S,[b+1,D-1,z,t]);
RemoveSet(S,[d-1,D-1,z,t]);
vbar([a+1,A+2,z],C-A-4,t);
vbar([c-1,A+2,z],C-A-4,t);
vbar([b+1,B+2,z],D-B-4,t);
vbar([d-1,B+2,z],D-B-4,t);


A:=NullMat(69,45);;
A:=List([1..10],i->1*A);;
A:=List([1..8],i->1*A);;

for aa in S do
A[aa[4]][aa[3]][aa[1]][aa[2]]:=1;
od;

return PureCubicalComplex(A);
end);
################################################################
################################################################

