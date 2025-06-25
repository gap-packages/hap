

############################################
InstallGlobalFunction(HAP_BianchiRegularCWComplex,
function(OQ, UMP)
local Y,BND,D,U,US,GV,GE,GF, H, R, F, F2, K, Triples, Triples2, 
Pairs, Pairs2, Singletons, Points, BoundariesSingletons, 
BoundariesPairs, S, Heights, H1, H2, H3, TripleToPoint, OnLine,
fun,fun2, PairToLine,pos,J,i,p,q,xx,yy,d,x,y,bnd, ind;
#Inputs a ring OQ of imaginary quadratic integers and a list UMP of 
#unimodular pairs that "minimally covers" the fundamental region.  

R:=HAP_BianchiFundamentalRectangle(OQ);
xx:=R[2]-R[1];;
yy:=(R[4]-R[3])*HAPSqrt(OQ!.bianchiInteger);;

U:=[];
if OQ!.bianchiInteger mod 4 =1 then
D:=Cartesian(xx*[-1..1],(xx/2+yy)*[-1..1]);
else
D:=Cartesian(xx*[-1..1],yy*[-1..1]);
fi;
Apply(D,d->d[1]+d[2]);
##########
##########
for p in UMP do
for d in D do
   q:=p{[1..3]};
   q[2]:=q[2]+d*q[1];
   Append(q,UnimodularPairCoordinates(OQ,q));
   Add(U,UnimodularPairStandardForm(q));
od;
od;
##########
##########
U:=DuplicateFreeList(U);

K:=SimplicialComplexOfUnimodularPairs(OQ,U);
Points:=K!.Points;
S:=K!.simplicesLst[3];
Heights:=[];
for i in [1..Length(S)]  do
   if Points[i] = fail then Heights[i]:=-infinity; 
   else
      Heights[i]:= HAP_HeightOfPointOnSphere(OQ,Points[i],U[S[i][1]]);
   fi;
od;

########
fun:=function(v);
if v[1]<R[1] then return false; fi;
if v[1]>R[2] then return false; fi;
if v[2]!.irrational<R[3] then return false; fi;
if v[2]!.irrational>R[4] then return false; fi;
return true;
end;
########

########
fun2:=function(v)
local p ,J;
J:=Positions(Points,v);
J:=Maximum(Heights{J});
for p in U do
  if HAP_HeightOfPointOnSphere(OQ,v,p) > J then return false; fi;
od;
return v;
end;
########

##########################################
PairToLine:=function(y)
local L;
L:= UnimodularIntersectingLine(OQ,U[y[1]], U[y[2]]);
L[2]:=Minimum(L[2],-L[2]);    #a kind of normal form!
return L;
end;
##########################################

##########################################
OnLine:=function(L,v)
local w,u;
w:=[L[2][2],-L[2][1]];
u:=L[1]-v;
return u*w=0;
end;
##########################################

##########################################
TripleToPoint:=function(y);
return Points[Position(S,y)];
end;
##########################################

F2:=Filtered(Points,x->not x=fail);
F2:=List(F2,fun2);
F2:=Filtered(F2,x->not x=false);
F2:=DuplicateFreeList(F2);

F:=Filtered(F2,fun);
#F is a list of vertices (given as coordinates) in the fundamental recatangle.
#F2 involves a larger rectangle.

Triples:=List(F,f->Positions(Points,f));
Triples:=Concatenation(Triples);
Triples:=DuplicateFreeList(Triples);
Triples:=S{Triples};
Triples2:=List(F2,f->Positions(Points,f));
Triples2:=Concatenation(Triples2);
Triples2:=DuplicateFreeList(Triples2);
Triples2:=S{Triples2};

#Triples correspond (possibly in a many-to-one fashion) to vertices 
#in the fundamental rectangle. Triples2 involves as larger rectangle.

Singletons:=DuplicateFreeList(Flat(Triples)); 
#Singletons correspond bijectively to the faces intersecting the 
#fundamental rectangle. So U{Singletons} is the same as P!.unimodularPairs.


Pairs2:=List(Triples2,t->[t{[1,2]},t{[1,3]},t{[2,3]}]);
Pairs2:=Concatenation(Pairs2);
Pairs2:=List(Pairs2,p->SSortedList(p));
Pairs2:=Filtered(Pairs2,x->HAP_AreStrictlyIntersectingUnimodularPairs(OQ,U[x[1]],U[x[2]])  );
Pairs2:=DuplicateFreeList(Pairs2);
Pairs2:=Classify(Pairs2,PairToLine);
#An element of Pairs2 is an equivalence class of pairs.

BoundariesSingletons:=[];
for x in Singletons do
bnd:=[];
for y in Pairs2 do
if x in Flat(y) then
   Add(bnd,y[1]);
fi;
od;
Add(BoundariesSingletons,bnd);
od;

Pairs:=Concatenation(BoundariesSingletons);
Pairs:=DuplicateFreeList(Pairs);
Pairs:=Classify(Pairs,PairToLine);

BoundariesPairs:=[];
for x in Pairs do
bnd:=[];
for y in Triples2 do
if x[1][1] in y and x[1][2] in y then 
pos:=Position(S,y);
J:=Positions(Points,Points[pos]);
J:=Maximum(Heights{J});
if HAP_HeightOfPointOnSphere(OQ,Points[pos],U[x[1][1]]) = J then
   Add(bnd,Points[pos]); 
fi;
fi;
od;

Add(BoundariesPairs,Intersection(SortedList(bnd),F2));
od;

BoundariesPairs:=List(BoundariesPairs,x->DuplicateFreeList(x));
ind:=[1..Length(Pairs)];;
ind:=Filtered(ind,i->Length(BoundariesPairs[i])>1);
Pairs:=Pairs{ind};
BoundariesPairs:=BoundariesPairs{ind};


BND:=[];
GV:=Concatenation(BoundariesPairs);
GV:=DuplicateFreeList(GV);
BND[1]:=List([1..Length(GV)],i->[1,0]);
BND[2]:=[];

for x in BoundariesPairs do
bnd:=[Position(GV,x[1]),Position(GV,x[2])];
bnd:=SortedList(bnd);
bnd[3]:=bnd[2];
bnd[2]:=bnd[1];
bnd[1]:=2;
Add(BND[2],bnd);
od;


BND[3]:=[];
for x in BoundariesSingletons do
bnd:=List(x,i->PositionProperty(Pairs,x->i in x));
bnd:=Filtered(bnd,x->not fail = x);
bnd:=SortedList(bnd);
bnd:=Reversed(bnd);
Add(bnd,Length(bnd));
bnd:=Reversed(bnd);
Add(BND[3],bnd);
od;

BND[4]:=[];

Heights:=Heights{List(GV,p->Position(Points,p))};
#Heights:=List(GV,p->  Maximum(Heights{Positions(Points,p)}) );
for i in [1..Length(Heights)] do
Heights[i]:=(Heights[i])!.rational;
od;
#Points:=List([1..Length(GV)],i->[GV[i][1],GV[i][2],Sqrt(Heights[i])]);
Points:=List([1..Length(GV)],i->[GV[i][1],GV[i][2],Heights[i]]);
Y:=RegularCWComplex(BND);
Y!.points:=Points;
Y!.ring:=OQ;
Y!.sphereCentres:=List(U{Singletons},u->UnimodularPairCoordinates(OQ,u));

return Y;

end);
############################################


