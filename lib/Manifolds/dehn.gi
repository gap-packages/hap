
#################################################################
#################################################################
InstallGlobalFunction(ThreeManifoldViaDehnSurgery,
function(K,pp,qq)
local start, edges, vertices, faces, mins, fn, TORUS, ffmod, fffmod,
      fmod, mymod, loop, B3, diff, BALL, MiddleEdge, LeftVert, RightVert, 
      bool, lloop, pos, G1, G, F1, F, loop2, loopvertices, Bvertices, LB,
      Lvertices,Tvertices, Rvertices, 
      cdba, dcab, bacd, abdc, CDBA, DCAB, BACD, ABDC,  
      DB, BD, CA, AC, db, bd, ac, dc, cd, ca, DC, CD, BA, AB, ab,
      ba, Vd, VD, VC, VB, Vc, Vb, Va, VA, T3, Pdcab, Pcdba, Pbacd,
      Pabdc, Pbd, Pca, Pac, Pcd, Pba, Pdb, Pdc, Pab, PVd, PVc, PVb, PVa,
      Vabdc, vabdc, sqboundary, BR1, BR0, BL1, BL0, TR1, TL0, TL1, TR0,
      xtrballs, xtrfaces, xtrloop, xtrFaces, xtrvertices, Faces, SGN,
      A,B,e,ee,f,i,k,L,LL,m,M,MB,p,P,q,Q,s,ss,S,T,v,W, x,Y,fg,epi,g,YY,ff,TT;

YY:=SphericalKnotComplement(K);;
ss:=100000;;
s:=Size(BoundaryOfPureRegularCWComplex(YY));
while s<ss do
   ss:=s;
   YY:=RegularCWComplex(BarycentricSubdivision(YY));
   YY:=SimplifiedComplex(YY);
   s:=Size(BoundaryOfPureRegularCWComplex(YY));
od;

if s>16 then Print("PROBLEM: Torus has ",s," cells\n"); fi;
ff:=BoundaryMap(YY);
TT:=Source(ff);
YY:=Target(ff);


B3:=List(TT!.boundaries[3],x->1*x{[2..5]});

abdc:=B3[1];   #I think this choice needs to be refined.

####################################################
#Let's figure out the "orientation"

fg:=FundamentalGroupWithPathReps(YY);
e:=fg!.edgeToWord(ff!.mapping(1,abdc[1]));
epi:=NqEpimorphismNilpotentQuotient(fg,1);
e:=Image(epi,e);
g:=Target(epi)/Group(e);
####################################################

if Order(g)=infinity then
p:=qq; q:=-pp;                   #WARNING:  This is fine for the trivial and the trefoil
else                             #          knots. Need to think more about arbitrary knots
p:=pp;q:=-(qq-p*(Order(g)+1));
fi;

SGN:=SignInt(p*q)<0;
p:=AbsInt(p); q:=AbsInt(q);
x:=Gcd(p,q);
p:=p/x; q:=q/x;

m:=p*q;
L:=[];
x:=[0,0];
start:=[0,0];
edges:=[];

while not x in L do
   Add(L,x);
   x:=[  (x[1]+p) mod m,
         (x[2]+q) mod m];
   if x[1]=0 and x[2]>0 then 
      Add(L,1*[m,x[2]]); 
      Add(edges,SortedList([start,1*[m,x[2]]]));
      start:=1*x;
   fi;
   if x[2]=0 and x[1]>0 then 
      Add(L,[x[1],m]); 
      Add(edges,SortedList([start,1*[x[1],m]]));
      start:=1*x;
      fi;
od;
      Add(L,[m,m]);
      Add(edges,SortedList([start,1*[m,m]]));


#GenLoop:=1*edges;  #This is the loop determining the torus homeomorphism
edges:=SSortedList(edges);

vertices:=[];;     #This will be the ordered list of vertices in the torus fd
faces:=[];
for x in edges do
   AddSet(vertices,x[1]);
   AddSet(vertices,x[2]);
od;
AddSet(vertices,[0,m]);
AddSet(vertices,[m,0]);
AddSet(vertices,[m,m]);


Bvertices:=Filtered(vertices,x->x[2]=0);
Tvertices:=Filtered(vertices,x->x[2]=m);
Lvertices:=Filtered(vertices,x->x[1]=0);
Rvertices:=Filtered(vertices,x->x[1]=m);

for x in Rvertices do
#AddSet(vertices, x+[p,0]);
AddSet(edges, SortedList(1*[x, x+[p,0]]));
AddSet(edges, SortedList(1*[x, x+[0,q]]));
AddSet(edges, SortedList(1*[x+[p,0], x+[p,q]]));
od;
for x in Tvertices do
#AddSet(vertices, 1*[x,x+[0,q]]);
AddSet(edges, SortedList(1*[x, x+[0,q]]));
AddSet(edges, SortedList(1*[x, x+[p,0]]));
AddSet(edges, SortedList(1*[x+[0,q],x+[p,q]]));
od;

for x in Lvertices do
AddSet(edges, SortedList(1*[x, x+[0,q]]));
od; 

for x in Bvertices do
AddSet(edges, SortedList(1*[x, x+[p,0]]));
od;


#AddSet(vertices,[m+p,m+q]);
#AddSet(edges,SortedList([[m,m] , [m+p,m]]));
AddSet(edges,SortedList([[m+p,m],[m+p,m+q]]));

vertices:=[];
for e in edges do
Append(vertices,1*e);
od;
vertices:=SSortedList(vertices);

Y:=RegularCWDiscreteSpace(Length(vertices));
for e in edges do
RegularCWComplex_AttachCellDestructive(Y,1,[Position(vertices,e[1]), Position(vertices,e[2])]);
od;

for i in [0..q-1] do
P:=[];
Add(P,Position(edges,[[i*p,m],[(i+1)*p,m]]));
Add(P,Position(edges,[[i*p,m+q],[(i+1)*p,m+q]]));
Add(P,Position(edges,[[i*p,m],[i*p,m+q]])  );
Add(P,Position(edges,[[(i+1)*p,m],[(i+1)*p,m+q]]));
Add(faces,P);
RegularCWComplex_AttachCellDestructive(Y,2,SortedList(P));
od;

for i in [0..p] do
P:=[];
Add(P,Position(edges,[[m,i*q],[m,(i+1)*q]]));
Add(P,Position(edges,[[m+p,i*q],[m+p,(i+1)*q]]));
Add(P,Position(edges,[[m,i*q],[m+p,i*q]])  );
Add(P,Position(edges,[[m,(i+1)*q],[m+p,(i+1)*q]]));
Add(faces,P);
RegularCWComplex_AttachCellDestructive(Y,2,SortedList(P));
od;


LL:=[];
for v in Lvertices do
L:=Filtered(edges, x->v in x);
for x in L do
if (not x[1][1]=x[2][1]) and (not x[1][2]=x[2][2]) then
Add(LL,x);
fi;
od;
od;

LB:=[];
for v in Bvertices do
L:=Filtered(edges, x->v in x);
for x in L do
if (not x[1][1]=x[2][1]) and (not x[1][2]=x[2][2]) then
Add(LB,x);
fi;
od;
od;

for i in [0..p-2] do
P:=[];
x:=[[0,i*q],[0,(i+1)*q]];
Add(P,x);
pos:=Position(Lvertices,x[1]); if pos=fail then Print("OOPS\,"); fi;
Add(P,LL[pos]);
pos:=Position(Lvertices,x[2]); if pos=fail then Print("OOPS\,"); fi;
Add(P,LL[pos]);
x:=SortedList([P[2][2], P[3][2]]);
Add(P,x);
P:=List(P,a->Position(edges,a));
Add(faces,P);
RegularCWComplex_AttachCellDestructive(Y,2,P);
od;
i:=p-1;
P:=[];
x:=[[0,i*q],[0,(i+1)*q]];
Add(P,x);
pos:=Position(Lvertices,x[1]); if pos=fail then Print("OOPS\,"); fi;
Add(P,LL[pos]);
x:=[[0,m],[p,m]];
Add(P,x);
P:=List(P,a->Position(edges,a));
Add(faces,P);
RegularCWComplex_AttachCellDestructive(Y,2,P);



for i in [0..q-2] do
P:=[];
x:=[[i*p,0],[(i+1)*p,0]];
Add(P,x);
pos:=Position(Bvertices,x[1]); if pos=fail then Print("OOPS\,"); fi;
Add(P,LB[pos]);
pos:=Position(Bvertices,x[2]); if pos=fail then Print("OOPS\,"); fi;
Add(P,LB[pos]);
x:=SortedList([P[2][2], P[3][2]]);
Add(P,x);
P:=List(P,a->Position(edges,a));
Add(faces,P);
RegularCWComplex_AttachCellDestructive(Y,2,P);
od;

P:=[];
x:=[[m-p,0],[m,0]];
Add(P,x);
pos:=Position(Bvertices,x[1]); if pos=fail then Print("OOPS\,"); fi;
Add(P,LB[pos]);
x:=[[m,0],[m,q]];
Add(P,x);
P:=List(P,a->Position(edges,a));
Add(faces,P);
RegularCWComplex_AttachCellDestructive(Y,2,P);

############################
vertices:=10*vertices;
edges:=10*edges;
edges:=List(edges,x->SortedList(x));
#faces stays uncganged
############################

loop:=[];
for x in edges do
if x[1][1]<>x[2][1] and x[1][2]<>x[2][2] then
Add(loop,x);
fi;
od;

loopvertices:=List(loop,x->1*x[2]);
for v in loopvertices do
if v[1]=10*m then Add(loop,[v,v+[10*p,0]]); fi;
if v[2]=10*m then Add(loop,[v,v+[0,10*q]]); fi;
od;
pos:=Position(loop,10*[[m,m],[m+p,m]]);
Remove(loop,pos);
Add(loop,[[m,0],[m+p,0]]*10);

loop2:=[];
if p<=q then                                 ##ADDED "=" May 2022
Add(loop2, [10*[m,m]-[0,q], 10*[m,m]+[p,0]]);
fi;
for x in loop do
if x[1][1]=0 and x[1][2]>0 and x[2][1]<10*m then
Add(loop2, [x[1]-[0,q], x[2]+[p,0]]); 
fi;
if x[1][1]=0 and x[1][2]>0 and x[2][1]=10*m then
Add(loop2, [x[1]-[0,q], x[2]+[0,-q]]);
if x[2][2]=10*m then Add(loop2, [x[2]+[0,-q], x[2]+[p,0]]); fi;###############
fi;

if x[1][2]=0 and x[1][1]<10*m and x[2][1]=10*m then
Add(loop2, [x[1]+[p,0], x[2]-[0,q]]); 
fi;

if x[1][2]=0 and x[1][1]<10*m and x[2][1]<10*m then
Add(loop2, [x[1]+[p,0], x[2]+[p,0]]);   #NEW
fi;

if x[1][2]=10*m and x[1][1]<10*m then
Add(loop2, [x[1]+[p,0], x[2]+[p,0]]); 
fi;
if x[1][2]=10*m and x[1][1]=10*m then
Add(loop2, [x[1]+[p,0], x[2]+[p,-q]]); 
fi;

if x[1][1]=10*m and x[1][2]>0  and x[1][2]<10*m then
Add(loop2, [x[1]+[0,-q], x[2]+[0,-q]]); 
fi;
od;

Add(loop2,[[m,m+q], [m+p,m+q]]*10);
Add(loop2,[10*[m,m+q]+[p,-q], 10*[m+p,m+q]-[0,q]]);
Add(loop2,[[0,m+q]*10-[0,q], [0,m+q]*10+[p,0]]);


for i in Concatenation([0..p-2],[p]) do
Add(loop2,[10*[0,i*q], 10*[0,(i+1)*q]-[0,q]]  );
od;
Add(loop2,[10*[0,(p-1)*q], 10*[0,p*q]]  );
Add(loop2,[10*[0,(p+1)*q]-[0,q], 10*[0,p*q+q]]  );

for i in [0..p-1] do
Add(loop2,[10*[m,i*q], 10*[m,(i+1)*q]-[0,q]]  );
od;

for i in Concatenation([0..p-2],[p]) do
Add(loop2,[10*[m+p,i*q], 10*[m+p,(i+1)*q]-[0,q]]  );
od;
Add(loop2,[10*[m+p,(p-1)*q], 10*[m+p,p*q]]  );
Add(loop2,[10*[m+p,(p+1)*q]-[0,q], 10*[m+p,p*q+q]]  );

for i in [0..q-1] do
Add(loop2,[10*[i*p,0]+[p,0], 10*[(i+1)*p,0]]  );
od;

for i in [0..q-1] do
Add(loop2,[10*[i*p,m+q]+[p,0], 10*[(i+1)*p,m+q]]  );
od;
Add(loop2,[10*[0,m+q],10*[0,m+q]+[p,0]]  );

for i in [1..q] do
Add(loop2,[10*[i*p,m]+[p,0], 10*[(i+1)*p,m]]  );
od;
Add(loop2,[10*[0,m], 10*[p,m]]  );

for i in [1..p-1] do
Add(loop2,[10*[0,i*q]+[0,-q], 10*[0,i*q]]  );
Add(loop2,[10*[m,i*q]+[0,-q], 10*[m,i*q]]  );
Add(loop2,[10*[m+p,i*q]+[0,-q], 10*[m+p,i*q]]  );
od;
i:=p; Add(loop2,[10*[m,i*q]+[0,-q], 10*[m,i*q]]  ); 


for i in [0..q-1] do
Add(loop2,[10*[i*p,0], 10*[i*p,0]+[p,0]]  );
od;
for i in [1..q] do
Add(loop2,[10*[i*p,m], 10*[i*p,m]+[p,0]]  );
od;
for i in [1..q-1] do
Add(loop2,[10*[i*p,m+q], 10*[i*p,m+q]+[p,0]]  );
od;




loop:=Concatenation(loop,loop2);
loop:=List(loop,x->SortedList(x));

###################################################################
###################################################################
##VISUAL CHECK THAT THE COMPLEX IS CORRECT
S:=[];

#for e in edges do
for e in loop do

if (not e[1][1]=e[2][1] ) and (not e[1][2]=e[2][2]) then
   x:=10*e[1];
   #while x[1]<= 10*m*10 and x[2]<=10*m*10 do
while x[1]<= 10*e[2][1] and x[2]<=10*e[2][2] do
      Add(S,x);
      x:=x+[p,q];
   od;
fi;
 
if e[1][1]=e[2][1] then 
for i in [e[1][2]..e[2][2]] do
Add(S,[10*e[2][1],10*i]);
od;
fi;

if e[1][2]=e[2][2] then
for i in [e[1][1]..e[2][1]] do
Add(S,[10*i,10*e[2][2]]);
od;
fi;
od;
#return S;
###################################################################
###################################################################

vertices:=[];
for e in loop do
Append(vertices,1*e);
od;
vertices:=SSortedList(vertices);

Y:=RegularCWDiscreteSpace(Length(vertices));
for e in loop do
RegularCWComplex_AttachCellDestructive(Y,1,[Position(vertices,e[1]), Position(vertices,e[2])]);
od;

TL0:=[];
TL1:=[];
TR0:=[];
TR1:=[];
BL0:=[];
BL1:=[];
BR0:=[];
BR1:=[];

P:=[];;
#Add(P,[10*[0,m+q]-[0,q],10*[0,m+q]]);
#Add(P,[10*[0,m+q],10*[0,m+q]+[p,0]]);
#Add(P,[10*[0,m+q]-[0,q],10*[0,m+q]+[p,0]]);
#Add(TL0,P);

P:=[];;
Add(P,[10*[0,m],10*[0,m+q]-[0,q]]);
Add(P,[10*[0,m+q]-[0,q], 10*[0,m+q]+[p,0]]);
Add(P,[10*[0,m+q]+[p,0], 10*[p,m+q]]);
Add(P,[10*[0,m], 10*[p,m]]);
Add(P,[10*[p,m], 10*[p,m+q]]);
Add(TL0,P);

for i in [1..q-1] do
P:=[];;
Add(P,[10*[i*p,m]+[p,0],10*[i*p,m+q]+[p,0]]);
Add(P,[10*[(i+1)*p,m],10*[(i+1)*p,m+q]]);
Add(P,[10*[i*p,m]+[p,0],10*[(i+1)*p,m]]);
Add(P,[10*[i*p,m+q]+[p,0],10*[(i+1)*p,m+q]]);
Add(TL0,P);
od;

for i in [1..q-1] do
P:=[];;
Add(P,[10*[i*p,m],10*[i*p,m+q]]);
Add(P,[10*[i*p,m]+[p,0],10*[i*p,m+q]+[p,0]]);
Add(P,[10*[i*p,m],10*[i*p,m]+[p,0]]);
Add(P,[10*[i*p,m+q],10*[i*p,m+q]+[p,0]]);
Add(TL1,P);
od;

P:=[];
Add(P,[10*[0,m+q]-[0,q],10*[0,m+q]]);
Add(P,[10*[0,m+q],10*[0,m+q]+[p,0]]);
Add(P,[10*[0,m+q]-[0,q],10*[0,m+q]+[p,0]]);
Add(TL1,P);



P:=[];;
Add(P,[10*[m,m]+[p,0],10*[m,m+q]+[p,-q]]);
Add(P,[10*[m,m]+[p,0],10*[m+p,m]]);
Add(P,[10*[m,m+q]+[p,-q],10*[m+p,m+q]+[0,-q]]);
Add(P,[10*[m+p,m],10*[m+p,m+q]+[0,-q]]);
Add(TR0,P);

P:=[];;
Add(P,[10*[m,m],10*[m,m+q]]);
Add(P,[10*[m,m],10*[m,m]+[p,0]]);
Add(P,[10*[m,m+q]+[p,-q],10*[m+p,m+q]+[0,-q]]);
Add(P,[10*[m,m+q],10*[m+p,m+q]]);
Add(P,[10*[m,m]+[p,0],10*[m,m+q]+[p,-q]]);
Add(P,[10*[m+p,m+q]+[0,-q],10*[m+p,m+q]]);
Add(TR1,P);

for i in [0..p-2] do
P:=[];
Add(P,[10*[m,i*q],10*[m+p,i*q]]);
Add(P,[10*[m,(i+1)*q]+[0,-q],10*[m+p,(i+1)*q]+[0,-q]]);
Add(P,[10*[m,i*q],10*[m,(i+1)*q]-[0,q]]);
Add(P,[10*[m+p,i*q],10*[m+p,(i+1)*q]-[0,q]]);
Add(BR0,P);
od;

P:=[];
i:=p-1;
Add(P,[10*[m,i*q],10*[m,(i+1)*q]-[0,q]]);
Add(P,[10*[m,(i+1)*q]-[0,q], 10*[m,(i+1)*q]+[p,0]]);
Add(P,[10*[m,(i+1)*q]+[p,0] , 10*[m+p,(i+1)*q]]);
Add(P,[10*[m,i*q],10*[m+p,i*q]]);
Add(P,[10*[m+p,i*q], 10*[m+p,(i+1)*q]]);
Add(BR0,P);

for i in [1..p-1] do
P:=[];
Add(P,[10*[m,i*q]-[0,q], 10*[m+p,i*q]-[0,q]]);
Add(P,[10*[m,i*q], 10*[m+p,i*q]]);
Add(P,[10*[m,i*q]-[0,q], 10*[m,i*q]]);
Add(P,[10*[m+p,i*q]-[0,q], 10*[m+p,i*q]]);
Add(BR1,P);
od;

P:=[];
Add(P,[10*[m,m], 10*[m,m]+[p,0]]);
Add(P,[10*[m,m]-[0,q], 10*[m,m]]);
Add(P,[10*[m,m]-[0,q], 10*[m,m]+[p,0]]);
Add(BR1,P);

F:=Filtered(loop,e->e[1][1]=0);
F:=Filtered(F,e->e[2][1]>0);  
F1:=List(F,e->e[1]);

for i in [0..p-2] do
P:=[];
A:=10*[0,i*q];
B:=10*[0,(i+1)*q]-[0,q];
Add(P, [A, B]);
pos:=Position(F1,A);
Add(P,F[pos]);
pos:=Position(F1,B);
Add(P,F[pos]);
Add(P,SortedList([P[2][2],P[3][2]]));
Add(BL0,P);
od;

P:=[];
A:=10*[0,m-q];
B:=10*[0,m];
pos:=Position(F1,A);
Add(P,F[pos]);
Add(P,[A,B]);
Add(P,[B,10*[p,m]]);
Add(BL0,P);

G:=Filtered(loop,e->e[1][2]=0);
G:=Filtered(G,e->e[2][2]>0);
G1:=List(G,e->e[1]);

for i in [0..q-2] do
P:=[];
A:=10*[i*p,0]+[p,0];
B:=10*[(i+1)*p,0];
Add(P,[A,B]);
pos:=Position(G1,A);
Add(P,G[pos]);
pos:=Position(G1,B);
Add(P,G[pos]);
Add(P,SortedList([P[2][2],P[3][2]]));
Add(BL0,P);
od;

P:=[];
A:=10*[m-p,0]+[p,0];
B:=10*[m,0];
Add(P,[A,B]);
pos:=Position(G1,A);
Add(P,G[pos]);
Add(P,[P[1][2],P[2][2]]);
Add(BL0,P);

for i in [1..p-1] do
P:=[];
A:=10*[0,i*q]-[0,q];
B:=10*[0,i*q];
Add(P,[A,B]);
pos:=Position(F1,A);
Add(P,F[pos]);
pos:=Position(F1,B);
Add(P,F[pos]);
Add(P, SortedList([P[2][2],P[3][2]]));
Add(BL1,P);
od;

for i in [0..q-1] do
P:=[];
A:=10*[i*p,0];
B:=10*[i*p,0]+[p,0];
Add(P,[A,B]);
pos:=Position(G1,A);
Add(P,G[pos]);
pos:=Position(G1,B);
Add(P,G[pos]);
Add(P, SortedList([P[2][2],P[3][2]]));
Add(BL1,P);
od;


faces:=Concatenation(TL0,TL1,TR0,TR1,BL0,BL1,BR0,BR1);

##########################
mymod:=function(x,m);
if x=infinity then return infinity; fi;
if x=-infinity then return -infinity; fi;
return x mod m;
end;
##########################

##########################
fmod:=function(x);
return [mymod( x[1] , 10*(m+p) ), mymod(x[2] , 10*(m+q) )];
end;
##########################
ffmod:=function(e);
return SortedList([ fmod(e[1]), fmod(e[2]) ]);
end;
##########################
fffmod:=function(P);
return SortedList(List(P,ffmod));
end;
##########################



vertices:=List(vertices,fmod);
vertices:=SSortedList(vertices);


lloop:=[];
for e in loop do
Add(lloop,ffmod(e));
od;
loop:=lloop;
 loop:=SSortedList(loop);


faces:=List(faces,fffmod);

Faces:=List(1*faces,P->List(P,e->Position(loop,1*e)));

Y:=RegularCWDiscreteSpace(Length(vertices));
for e in loop do
RegularCWComplex_AttachCellDestructive(Y,1,SortedList([Position(vertices,e[1]), Position(vertices,e[2])]));
od;

for P in Faces do
RegularCWComplex_AttachCellDestructive(Y,2,P);
od;


#Now add the thickened relator disk to Y=TORUS
LeftVert:=fmod([-infinity,-infinity]);
RightVert:=fmod([infinity,infinity]);
MiddleEdge:=ffmod([LeftVert,RightVert]);
xtrvertices:=[LeftVert,RightVert];
Append(vertices,xtrvertices);
xtrloop:=[MiddleEdge];
xtrfaces:=[];
xtrballs:=[];


for P in TL1 do
BALL:=[fffmod(P)];
if Size(P)=4 then
for e in P do
A:=e[1]; B:=e[2];
if A[2]<B[2] and (A[1] mod (10*p)=0) then
   Add(xtrloop,ffmod([A,LeftVert])); 
   Add(xtrloop,ffmod([B,LeftVert]));
   Q:=[];
   Add(Q, e);
   Add(Q, ffmod([A,LeftVert]));
   Add(Q, ffmod([B,LeftVert]));
   Add(xtrfaces,fffmod(Q));
   Add(BALL,fffmod(Q));
fi;
if A[2]<B[2] and not (A[1] mod (10*p)=0) then 
   Add(xtrloop,ffmod([A,RightVert]));  
   Add(xtrloop,ffmod([B,RightVert]));
   Q:=[];
   Add(Q, e);
   Add(Q, ffmod([A, RightVert]));
   Add(Q, ffmod([B,RightVert]));
   Add(xtrfaces,fffmod(Q));
   Add(BALL,fffmod(Q));
fi;
if A[2]=B[2]  then 
   ee:=SortedList(e);
   Q:=[e,MiddleEdge, [ee[1],LeftVert], [ee[2],RightVert]];
   Q:=fffmod(Q);
   Add(xtrfaces,Q);
   Add(BALL,fffmod(Q));
fi;
od;
fi;
if Size(P)=3 then
for e in P do
A:=e[1]; B:=e[2];
if A[1]=0 and B[1]=0 then 
   Add(xtrloop, ffmod([B,LeftVert]));   
   Add(xtrloop, ffmod([A,RightVert]));
   Q:=[e,[B,LeftVert],[A,RightVert],MiddleEdge];
   Q:=fffmod(Q);
   Add(xtrfaces,Q);
   Add(BALL,fffmod(Q));
fi;
if A[2]=10*(m+q) and B[2]=10*(m+q) then 
   Add(xtrloop, ffmod([A,LeftVert]));
   Add(xtrloop, ffmod([B,RightVert]));
   Q:=[e,[A,LeftVert],[B,RightVert],MiddleEdge];
   Q:=fffmod(Q);
   Add(xtrfaces,Q);
   Add(BALL,fffmod(Q));
fi;
if (not A[1]=B[1]) and (not A[2]=B[2])  then 
   Q:=[e,[A,RightVert],[B,RightVert]];
   Q:=fffmod(Q);
   Add(xtrfaces,Q);   ##########??
   Add(BALL,fffmod(Q));
fi;
od;
fi;
Add(xtrballs,BALL);
od;

for P in TR1 do
BALL:=[fffmod(P)];
for e in P do
A:=e[1]; B:=e[2];
bool:=true;
if (A[1]=10*m and B[1]=10*m) 
   or 
   (A[2]=10*(m+q) and B[2]=10*(m+q)) 
   then
   bool:=false;
   Add(xtrloop, ffmod([A,LeftVert]) );
   Add(xtrloop, ffmod([B,LeftVert]) );
   Q:=[e, ffmod([A,LeftVert]),  ffmod([B,LeftVert]) ];
   Q:=fffmod(Q);
   Add(xtrfaces,Q);
   Add(BALL,fffmod(Q));
fi;
if (A[1]=(10*m+p) and B[1]=(10*m+p) )
   or
   (A[2]=(10*m+9*q) and B[2]=(10*m+9*q) )
   then
   bool:=false;
   Add(xtrloop, ffmod([A,RightVert]) );
   Add(xtrloop, ffmod([B,RightVert]) );
   Q:=[e, ffmod([A,RightVert]),  ffmod([B,RightVert]) ];
   Q:=fffmod(Q);
   Add(xtrfaces,Q);
   Add(BALL,fffmod(Q));
fi;
if bool then
   if A[1]=(10*m)  then 
   Q:=[e, MiddleEdge, ffmod([A,LeftVert]),  ffmod([B,RightVert]) ];
   Q:=fffmod(Q);
   Add(xtrfaces,Q);   ###########??
   Add(BALL,fffmod(Q));

   else 

   Q:=[e, MiddleEdge, ffmod([B,LeftVert]),  ffmod([A,RightVert]) ];
   Q:=fffmod(Q);
   Add(xtrfaces,Q);  ################??
   Add(BALL,fffmod(Q));
   fi;
fi;
od;
Add(xtrballs,BALL);
od;


for P in BR1 do
BALL:=[fffmod(P)];
if Size(P)=4 then
for e in P do
A:=e[1]; B:=e[2];
if A[2]=B[2] and A[2] mod (10*q) =0 then
   Add(xtrloop,ffmod([A,LeftVert]));
   Add(xtrloop,ffmod([B,LeftVert]));
   Q:=[];
   Add(Q, e);
   Add(Q, ffmod([A, LeftVert]));
   Add(Q, ffmod([B,LeftVert]));
   Add(xtrfaces,fffmod(Q));
   Add(BALL,fffmod(Q));
fi;
if A[2]=B[2] and  not (A[2]) mod (10*q) =0 then
   Add(xtrloop,ffmod([A,RightVert]));
   Add(xtrloop,ffmod([B,RightVert]));
   Q:=[];
   Add(Q, e);
   Add(Q, ffmod([A, RightVert]));
   Add(Q, ffmod([B,RightVert]));
   Add(xtrfaces,fffmod(Q));
   Add(BALL,fffmod(Q));
fi;
if A[1]=B[1] then
   Q:=[e, MiddleEdge, [A,RightVert], [B,LeftVert]];
   Add(xtrfaces,fffmod(Q));
   Add(BALL,fffmod(Q));
fi;
od;
fi;

if Size(P)=3 then
for e in P do
A:=e[1]; B:=e[2];
if A[1]=10*m and B[1]=10*m then 
  Add(xtrloop, [A,RightVert]);
  Add(xtrloop, [B,LeftVert]); 
  Q:=[e,[A,RightVert], [B,LeftVert], MiddleEdge];
  Add(xtrfaces,fffmod(Q));
  Add(BALL,fffmod(Q));
fi;
if A[2]=10*m and B[2]=10*m then
  Add(xtrloop, [B,RightVert]);
  Add(xtrloop, [A,LeftVert]);
  Q:=[e,[B,RightVert], [A,LeftVert], MiddleEdge];
  Add(xtrfaces,fffmod(Q)); #######################??
  Add(BALL,fffmod(Q));
fi;
if (not A[1]=B[1]) and (not A[2]=B[2]) then
  Q:=[e,[A,RightVert], [B,RightVert]];
  Add(xtrfaces,fffmod(Q));  ######################??
  Add(BALL,fffmod(Q));
fi;
od;
fi;
Add(xtrballs,BALL);
od;

for P in BL1 do
BALL:=[fffmod(P)];
for e in P do
A:=e[1]; B:=e[2];
if (not A[1]=B[1]) and (not A[2]=B[2]) 
   and 
   (A[2] mod (10*q)=0) and (B[2] mod (10*q)=0) 
   and
   A[1] mod (10*p)=0 then     #NEW########################
   Add(xtrloop, ffmod([A,LeftVert]));
   Add(xtrloop, ffmod([B,LeftVert]));
   Q:=[e,[A,LeftVert], [B,LeftVert]];
   Add(xtrfaces,fffmod(Q)); 
   Add(BALL,fffmod(Q)); 
fi;
if (not A[1]=B[1]) and (not A[2]=B[2])
   and
   (A[2] mod (10*q)=0) and (B[2] mod (10*q)=0) 
   and
   A[1] mod (10*p)=p then     #NEW########################
   Add(xtrloop, ffmod([A,RightVert]));
   Add(xtrloop, ffmod([B,RightVert]));
   Q:=[e,[A,RightVert], [B,RightVert]];
   Add(xtrfaces,fffmod(Q));
   Add(BALL,fffmod(Q));
fi;

if (not A[1]=B[1]) and (not A[2]=B[2]) 
   and 
   (( not A[2] mod (10*q)=0) or  ( not B[2] mod (10*q)=0))
   then
   Add(xtrloop, ffmod([A,RightVert]));
   Add(xtrloop, ffmod([B,RightVert]));
   Q:=[e,[A,RightVert], [B,RightVert]];
   Add(xtrfaces,fffmod(Q)); ###############??
   Add(BALL,fffmod(Q)); 
fi;
if A[1]=B[1] and A[1]=0
   then
   Q:=[e,[A,RightVert], [B,LeftVert], MiddleEdge];
   Add(xtrfaces,fffmod(Q)); ###################??
   Add(BALL,fffmod(Q));
fi;
if A[1]=B[1] and A[1]=10*m
   then
   Q:=[e,[A,RightVert], [B,LeftVert], MiddleEdge];
   Add(xtrfaces,fffmod(Q)); ###########################??
   Add(BALL,fffmod(Q));
fi;
if A[2]=B[2] and A[2]=0
   then
   Q:=[e,[A,LeftVert], [B,RightVert], MiddleEdge];
   Add(xtrfaces,fffmod(Q)); ######################??
   Add(BALL,fffmod(Q));
fi;
if A[2]=B[2] and A[2]=10*m
   then
   Q:=[e,[A,LeftVert], [B,RightVert], MiddleEdge];
   Add(xtrfaces,fffmod(Q)); ##########################??
   Add(BALL,fffmod(Q));
fi;
od;
Add(xtrballs,BALL);
od;



vertices:=List(vertices,fmod);
xtrloop:=List(xtrloop,e->ffmod(e));
xtrloop:=SSortedList(xtrloop);
Append(loop,xtrloop);
xtrfaces:=List(xtrfaces,f->fffmod(f));
xtrfaces:=SSortedList(xtrfaces);
xtrFaces:=List(xtrfaces,P->List(P,e->Position(loop,e)));
Append(Faces,xtrFaces);
Append(faces,xtrfaces);
xtrballs:=List(xtrballs,b->SortedList(b));
xtrballs:=SSortedList(xtrballs);
xtrballs:=List(xtrballs,P->SortedList(List(P,e->Position(faces,e))));


for v in xtrvertices do
RegularCWComplex_AttachCellDestructive(Y,0);
od;


for e in xtrloop do
RegularCWComplex_AttachCellDestructive(Y,1,SortedList([Position(vertices,e[1]), Position(vertices,e[2])]));
od;


for P in xtrFaces do
RegularCWComplex_AttachCellDestructive(Y,2,P);
od;

for P in xtrballs do
RegularCWComplex_AttachCellDestructive(Y,3,P);
od;

TORUS:=Y;
#return TORUS;

##   A----------B-----------A
##   |          |           |
##   |          |           |
##   C----------D-----------C
##   |          |           |
##   |          |           |
##   A----------B-----------A

VA:=10*[0,m+q];;
VB:=10*[m,m+q];;
VC:=10*[0,m];;
VD:=10*[m,m];;
AB:=[];
for i in [0..q-1] do
e:= [10*[i*p,m+q], 10*[i*p,m+q]+[p,0]];
Add(AB,e);
f:=[10*[i*p,m+q]+[p,0], 10*[(i+1)*p,m+q]];
Add(AB,f);
od;
CD:=[ [10*[0,m],10*[p,m]]];
for i in [1..q-1] do
e:= [10*[i*p,m], 10*[i*p,m]+[p,0]];
Add(CD,e);
f:=[10*[i*p,m]+[p,0], 10*[(i+1)*p,m]];
Add(CD,f);
od;
BA:=[ [10*[m,m+q], 10*[m+p,m+q]] ];
DC:=[ [10*[m,m] , 10*[m,m]+[p,0]] , [10*[m,m]+[p,0], 10*[m+p,m]]  ];
CA:=[];
for i in [0..p-2] do
e:=[10*[0,i*q],10*[0,(i+1)*q]+[0,-q]];
Add(CA,e);
f:=[10*[0,(i+1)*q]+[0,-q],  10*[0,(i+1)*q]];
Add(CA,f);
od;
Add(CA,[10*[0,(p-1)*q],10*[0,m]]);

DB:=[];
for i in [0..p-1] do
e:=[10*[m,i*q],10*[m,(i+1)*q]+[0,-q]];
Add(DB,e);
f:=[10*[m,(i+1)*q]+[0,-q],  10*[m,(i+1)*q]];
Add(DB,f);
od;

AC:=[ [10*[0,m],10*[0,m+q]-[0,q]] , [ 10*[0,m+q]-[0,q], 10*[0,m+q]]     ];
BD:=[ [10*[m,m],10*[m,m+q]]  ];

ABDC:=Concatenation(TL0,TL1);
BACD:=Concatenation(TR0,TR1);
CDBA:=Concatenation(BL0,BL1);
DCAB:=Concatenation(BR0,BR1);

##
##The next peice of code is a continuation of the code at the beginning!  Sloppy!!
##
f:=ff;
T:=TT;
Y:=YY;

for P in B3{[2,3,4]} do
if Intersection(abdc,P)=[] then dcab:=P; break; fi;
od;

diff:=Difference(B3,[abdc,dcab]);
vabdc:=List(abdc,i->T!.boundaries[2][i]{[2,3]});
sqboundary:=1*vabdc;
Vabdc:=[vabdc[1]]; vabdc:=Difference(vabdc,Vabdc);
for P in vabdc do
if Size(Intersection(P,Vabdc[1]))>0 then Add(Vabdc,P); vabdc:=Difference(vabdc,Vabdc); break; fi;
od;
for P in vabdc do
if Size(Intersection(P,Vabdc[2]))>0 then Add(Vabdc,P); vabdc:=Difference(vabdc,Vabdc); break; fi;
od;
Add(Vabdc,vabdc[1]); 
for i in [2,3,4] do
if Vabdc[i-1][2]<>Vabdc[i][1] then Vabdc[i]:=Reversed(Vabdc[i]); fi;
od;
if Vabdc[1][2]<>Vabdc[2][1] then Vabdc[1]:=Reversed(Vabdc[1]); fi;

Va:=Vabdc[1][1];
Vb:=Vabdc[2][1];
Vd:=Vabdc[3][1];
Vc:=Vabdc[4][1];

ab:=abdc[Position(sqboundary,SortedList(Vabdc[1]))];
L:=Filtered([1..8],i->T!.boundaries[2][i]{[2,3]}=SortedList(Vabdc[1]));
ba:=Difference(L,[ab])[1];

bd:=abdc[Position(sqboundary,SortedList(Vabdc[2]))];
L:=Filtered([1..8],i->T!.boundaries[2][i]{[2,3]}=SortedList(Vabdc[2]));
db:=Difference(L,[bd])[1];

cd:=abdc[Position(sqboundary,SortedList(Vabdc[3]))];
L:=Filtered([1..8],i->T!.boundaries[2][i]{[2,3]}=SortedList(Vabdc[3]));
dc:=Difference(L,[cd])[1];

ac:=abdc[Position(sqboundary,SortedList(Vabdc[4]))];
L:=Filtered([1..8],i->T!.boundaries[2][i]{[2,3]}=SortedList(Vabdc[4]));
ca:=Difference(L,[ac])[1];

B3:=List(B3,x->SSortedList(x));
abdc:=Position(B3, SSortedList([ab,bd,cd,ac]));
bacd:=Position(B3, SSortedList([ba,ac,dc,bd]));
cdba:=Position(B3, SSortedList([cd,db,ab,ca]));
dcab:=Position(B3, SSortedList([dc,ca,ba,db]));

Va:=f!.mapping(0,Va);
Vb:=f!.mapping(0,Vb);
Vc:=f!.mapping(0,Vc);
Vd:=f!.mapping(0,Vd);
ab:=f!.mapping(1,ab);
ba:=f!.mapping(1,ba);
cd:=f!.mapping(1,cd);
dc:=f!.mapping(1,dc);
ac:=f!.mapping(1,ac);
ca:=f!.mapping(1,ca);
bd:=f!.mapping(1,bd);
db:=f!.mapping(1,db);
abdc:=f!.mapping(2,abdc);
bacd:=f!.mapping(2,bacd);
cdba:=f!.mapping(2,cdba);
dcab:=f!.mapping(2,dcab);

if SGN then
############################################
PVa:=Vb;
PVb:=Va;
PVc:=Vd;
PVd:=Vc;
Pab:=ab;
Pba:=ba;
Pcd:=cd;
Pdc:=dc;
Pac:=bd;
Pca:=db;
Pbd:=ac;
Pdb:=ca;
Pabdc:=abdc;
Pbacd:=bacd;
Pcdba:=cdba;
Pdcab:=dcab;

Va:=PVa;
Vb:=PVb;
Vc:=PVc;
Vd:=PVd;
ab:=Pab;
ba:=Pba;
cd:=Pcd;
dc:=Pdc;
ac:=Pac;
ca:=Pca;
bd:=Pbd;
db:=Pdb;
abdc:=Pabdc;
bacd:=Pbacd;
cdba:=Pcdba;
dcab:=Pdcab;
############################################
fi;




M:=RegularCWComplex_DisjointUnion(Y,TORUS);
MB:=1*M!.boundaries;



VA:=fmod(VA);
VB:=fmod(VB);
VC:=fmod(VC);
VD:=fmod(VD);
AB:=List(AB,ffmod);
AB:=List(AB,x->SortedList(x));
AB:=List(AB,e->Position(loop,e))+Y!.nrCells(1);
BA:=List(BA,ffmod);
BA:=List(BA,x->SortedList(x));
BA:=List(BA,e->Position(loop,e))+Y!.nrCells(1);
CD:=List(CD,ffmod);
CD:=List(CD,x->SortedList(x));
CD:=List(CD,e->Position(loop,e))+Y!.nrCells(1);
DC:=List(DC,ffmod);
DC:=List(DC,x->SortedList(x));
DC:=List(DC,e->Position(loop,e))+Y!.nrCells(1);
AC:=List(AC,ffmod);
AC:=List(AC,x->SortedList(x));
AC:=List(AC,e->Position(loop,e))+Y!.nrCells(1);
CA:=List(CA,ffmod);
CA:=List(CA,x->SortedList(x));
CA:=List(CA,e->Position(loop,e))+Y!.nrCells(1);
BD:=List(BD,ffmod);
BD:=List(BD,x->SortedList(x));
BD:=List(BD,e->Position(loop,e))+Y!.nrCells(1);
DB:=List(DB,ffmod);
DB:=List(DB,x->SortedList(x));
DB:=List(DB,e->Position(loop,e))+Y!.nrCells(1);
ABDC:=List(ABDC,fffmod);
BACD:=List(BACD,fffmod);
CDBA:=List(CDBA,fffmod);
DCAB:=List(DCAB,fffmod);

ABDC:=List(1*ABDC,P->SortedList(List(P,e->Position(loop,1*e))));
BACD:=List(1*BACD,P->SortedList(List(P,e->Position(loop,1*e))));
CDBA:=List(1*CDBA,P->SortedList(List(P,e->Position(loop,1*e))));
DCAB:=List(1*DCAB,P->SortedList(List(P,e->Position(loop,1*e))));
T3:=List(TORUS!.boundaries[3],x->x{[2..Length(x)]});;
ABDC:=List(ABDC,P->Position(T3,P));
BACD:=List(BACD,P->Position(T3,P));
CDBA:=List(CDBA,P->Position(T3,P));
DCAB:=List(DCAB,P->Position(T3,P));





#Reduce the number of 0-cells
MB[1]:=MB[1]{[1..Length(MB[1])-4]};
mins:=[Y!.nrCells(0)+Position(vertices,VA),
 Y!.nrCells(0)+Position(vertices,VB),
Y!.nrCells(0)+Position(vertices,VC),
Y!.nrCells(0)+Position(vertices,VD)];
mins:=SSortedList(mins);
#############################
fn:=function(k);
if k=Y!.nrCells(0)+Position(vertices,VA) then  return Va; fi;
if k=Y!.nrCells(0)+Position(vertices,VB) then  return Vb; fi;
if k=Y!.nrCells(0)+Position(vertices,VC) then  return Vc; fi;
if k=Y!.nrCells(0)+Position(vertices,VD) then  return Vd; fi;
if k<mins[1] then return k; fi;
if k<mins[2] then return k-1; fi;
if k<mins[3] then return k-2; fi;
if k<mins[4] then return k-3; fi;
return k-4;
end;
#############################

MB[2]:=List(MB[2],x->Concatenation( [[x[1]],List(x{[2..x[1]+1]},k->fn(k))] ));

#Reduce the number of 1-cells
mins:=SSortedList([ab, ba, cd, dc, ac, ca, bd, db]);
############################
fn:=function(k);
if  k=ab then return AB; fi; 
if  k=ba then return BA; fi; 
if  k=cd then return CD; fi; 
if  k=dc then return DC; fi; 
if  k=ac then return AC; fi; 
if  k=ca then return CA; fi; 
if  k=bd then return BD; fi; 
if  k=db then return DB; fi; 
return k;
end;
############################
MB[3]:=List(MB[3],x->Concatenation( [[x[1]],List(x{[2..x[1]+1]},k->fn(k))] ));
MB[3]:=List(MB[3],x->Flat(x));
MB[3]:=List(MB[3],x->Concatenation( [[Length(x)-1],x{[2..Length(x)]}] ));
############################
fn:=function(k);
if k<mins[1] then return k; fi;
if k<mins[2] then return k-1; fi;
if k<mins[3] then return k-2; fi;
if k<mins[4] then return k-3; fi;
if k<mins[5] then return k-4; fi;
if k<mins[6] then return k-5; fi;
if k<mins[7] then return k-6; fi;
if k<mins[8] then return k-7; fi;
return k-8;
end;
############################
MB[3]:=List(MB[3],x->Concatenation([x[1]],List(x{[2..Length(x)]},fn)));
L:=[1..Length(MB[2])];
L:=Difference(L,mins);
MB[2]:=MB[2]{L};

#Reduce the number of 2-cells;
mins:=SSortedList([abdc, bacd, cdba, dcab]);
############################
fn:=function(k);
#return k;
if k=abdc then return ABDC+Y!.nrCells(2); fi;
if k=bacd then return BACD+Y!.nrCells(2); fi;
if k=cdba then return CDBA+Y!.nrCells(2); fi;
if k=dcab then return DCAB+Y!.nrCells(2); fi;
return k;
end;
############################
MB[4]:=List(MB[4],x->Concatenation([x[1]],List(x{[2..Length(x)]},fn)));
MB[4]:=List(MB[4],x->Flat(x));
MB[4]:=List(MB[4],x->Concatenation([Length(x)-1],x{[2..Length(x)]}));

############################
fn:=function(k);
if k<mins[1] then return k; fi;
if k<mins[2] then return k-1; fi;
if k<mins[3] then return k-2; fi;
if k<mins[4] then return k-3; fi;
return k-4;
end;
############################
MB[4]:=List(MB[4],x->Concatenation([x[1]],List(x{[2..Length(x)]},fn)));
L:=[1..Length(MB[3])];
L:=Difference(L,mins);
MB[3]:=MB[3]{L};


W:=RegularCWComplex(MB);
P:=Filtered([1..W!.nrCells(2)],k->W!.coboundaries[3][k][1]=1);;
RegularCWComplex_AttachCellDestructive(W,3,P);
#return SimplifiedComplex(W);
return BarycentricallySimplifiedComplex(W);

end);
################################################################
################################################################


##########################################################
##########################################################
InstallGlobalFunction(HAP_AllHomomorphisms,
function(G,H)
local gensG, EltsH, S, L, x, h, s;
#G can be infinite but H must be finite (and very small).
#A VERY INEFFICIENT IMPLEMENTATION

gensG:=GeneratorsOfGroup(G);
S:=SymmetricGroup(Length(gensG));
EltsH:=Elements(H);
L:=[];

for x in Tuples(EltsH,Length(gensG)) do
h:=GroupHomomorphismByImages(G,H,gensG,x);
if not h=fail then Add(L,h);  fi;
od;
return L;

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(DijkgraafWittenInvariant,
function(W,G)
local R,C,F,homs, inv, h, A, B;

R:=ResolutionFiniteGroup(G,4);
C:=ChainComplexOfUniversalCover(W);
F:=C!.group;
homs:=HAP_AllHomomorphisms(F,G);
SetSize(F,infinity);

inv:=[];
for h in homs do
Add(inv,Homology(TensorWithIntegers(EquivariantChainMap(C,R,h)),3));
od;

A:=SortedList(List(inv,f->Image(f,GeneratorsOfGroup(Source(f))[1])));
B:=SortedList(List(inv,f->Image(f,GeneratorsOfGroup(Source(f))[1]^-1)));
return Minimum(A,B);
end);
##########################################################
##########################################################

