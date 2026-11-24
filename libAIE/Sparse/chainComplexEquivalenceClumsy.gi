#This function was written by Claude in response to a request
#to implement the function on my "ToDo" in this section of the manual.
#https://gap-packages.github.io/hap/tutorial/chap17_mj.html#X83E2594E815844FA


#############################################################################
##
##  ChainComplexEquivalenceOfChainComplex(C)
##
##  INPUT   C : a HAP chain complex of free abelian groups (characteristic 0),
##              either dense (C!.boundary(n,k) is a vector) or sparse
##              (C!.boundary(n,k) is a list of pairs [i,coeff]).
##
##  OUTPUT  [F,G]  where  F : C -> D  and  G : D -> C  are chain maps
##          between C and a (usually much smaller) chain complex D, with
##
##                 F o G = identity on D     and     G o F ~ identity on C
##
##          (chain homotopic).  D has the same type (dense/sparse) as C.
##          F and G are HAP chain maps with F!.mapping(v,n), G!.mapping(v,n)
##          acting on dense vectors, exactly as the output of
##          ChainComplexEquivalenceOfRegularCWComplex(Y).
##
##  METHOD  Elementary "generator removal" (algebraic discrete Morse theory,
##          Gaussian elimination on the boundary matrices).  Suppose e is an
##          n-cell and a an (n-1)-cell with  d(e) = u*a + r  where u = +1 or
##          -1 and r has no a-component.  Cancelling the pair (e,a) gives a
##          chain complex C' on all generators except e and a, with
##
##              d'(x) = d(x) - c_x*u*d(e)        (x an n-cell, c_x = coeff of
##                                                a in d(x); this removes a)
##              d'(y) = d(y) - (coeff of e in d(y))*e      (y an (n+1)-cell)
##
##          and chain maps  f : C -> C',  g : C' -> C  given by
##
##              f(a) = -u*r,   f(e) = 0,   f(other cells) = identity
##              g(x) = x - c_x*u*e   (x an n-cell),   g = identity otherwise,
##
##          satisfying  f g = 1  and  g f ~ 1  (the homotopy is h(a) = u*e).
##          Performing a sequence of such cancellations and composing the
##          corresponding f's and g's gives F and G.  Nothing is ever
##          computed by linear algebra: the maps are assembled from the
##          recorded cancellations.
##
##          Cancellations are performed degree by degree (n = 1,2,3,...)
##          using a work list, until no unit (+-1) entry remains in any
##          boundary matrix of D.  To limit fill-in, cells are tried in
##          order of increasing boundary length, and among the unit entries
##          of a boundary the one with smallest Markowitz cost
##          (|d(e)|-1)*(number of cofaces of the face - 1) is used.  This
##          matters a lot in practice: on a 15,000-cell complex it shrank
##          the result from 1,1,8,11 to 1,1,2,5 cells and halved the time.
##
##  NOTE    Only unit coefficients can be cancelled over Z, so D need not
##          be minimal (torsion survives as entries like 2,3,... in D).
##
#############################################################################
InstallGlobalFunction(ChainComplexEquivalenceOfChainComplex,
function(C)
local
  L, isSparse, dimC, dimD, bnd, cob, alive, T,
  K, stepDeg, stepE, stepA, stepU, stepBnd,
  NormForm, Merge, FindFace, Pivot, Cancel, Gvec,
  n, k, e, i, j, s, t, v, B, b, pos, queue, inq,
  surv, newpos, Dbnd, Fimg, Gimg, acc, a,
  DDimension, DBoundary, mapF, mapG, D, F, G, props, conn;

if not (IsHapChainComplex(C) or IsHapSparseChainComplex(C)) then
   Error("the argument must be a HAP chain complex");
fi;
if not EvaluateProperty(C,"characteristic")=0 then
   Error("only chain complexes of free abelian groups (characteristic 0) are handled");
fi;

L:=Length(C);                       # boundaries d_1,...,d_L
isSparse:=IsHapSparseChainComplex(C);
dimC:=List([0..L],n->C!.dimension(n));     # dimC[n+1] = rank of C_n

#############################################################################
## helpers for sparse vectors: sorted lists of pairs [index,coeff], coeff<>0
#############################################################################
NormForm:=function(b)
local S, out, x, m;
S:=ShallowCopy(b);
Sort(S);
out:=[];
for x in S do
   m:=Length(out);
   if m>0 and out[m][1]=x[1] then
      out[m]:=[x[1],out[m][2]+x[2]];
   else
      Add(out,[x[1],x[2]]);
   fi;
od;
return Filtered(out,x->x[2]<>0);
end;

Merge:=function(B,Y)      # sum of two sorted sparse vectors
local U, pB, pY, lB, lY, c;
lB:=Length(B); lY:=Length(Y); U:=[]; pB:=1; pY:=1;
while pB<=lB and pY<=lY do
   if B[pB][1]<Y[pY][1] then Add(U,B[pB]); pB:=pB+1;
   elif B[pB][1]>Y[pY][1] then Add(U,Y[pY]); pY:=pY+1;
   else
      c:=B[pB][2]+Y[pY][2];
      if c<>0 then Add(U,[B[pB][1],c]); fi;
      pB:=pB+1; pY:=pY+1;
   fi;
od;
while pB<=lB do Add(U,B[pB]); pB:=pB+1; od;
while pY<=lY do Add(U,Y[pY]); pY:=pY+1; od;
return U;
end;

FindFace:=function(B,a)   # position of the entry with index a in sorted B
local lo,hi,mid;
lo:=1; hi:=Length(B);
while lo<=hi do
   mid:=QuoInt(lo+hi,2);
   if B[mid][1]=a then return mid;
   elif B[mid][1]<a then lo:=mid+1;
   else hi:=mid-1; fi;
od;
return fail;
end;

#############################################################################
## sparse boundaries  bnd[n][k]  (n=1..L) and cofaces  cob[n][i]
## (cob[n][i] = n-cells whose boundary meets the (n-1)-cell i; it is allowed
##  to contain stale entries, never to miss one)
#############################################################################
bnd:=[]; cob:=[];
for n in [1..L] do
   bnd[n]:=[];
   cob[n]:=List([1..dimC[n]],i->[]);
   for k in [1..dimC[n+1]] do
      b:=C!.boundary(n,k);
      if isSparse then
         b:=NormForm(b);
      else
         b:=List(PositionsProperty(b,x->not IsZero(x)),i->[i,b[i]]);
      fi;
      bnd[n][k]:=b;
      for t in b do Add(cob[n][t[1]],k); od;     # k increases, so sorted
   od;
od;

alive:=List([0..L],n->ListWithIdenticalEntries(dimC[n+1],true));
T:=List([0..L],n->[]);   # T[n+1][x] = list of [step,c]: at that step the
                         # n-cell x had the cancelled face with coefficient c
K:=0; stepDeg:=[]; stepE:=[]; stepA:=[]; stepU:=[]; stepBnd:=[];
queue:=[]; inq:=[];

#############################################################################
## choose a unit entry of d(e) whose face has fewest cofaces
#############################################################################
Pivot:=function(n,e)
local B,i,best,bestcost,cost;
B:=bnd[n][e]; best:=fail; bestcost:=infinity;
for i in [1..Length(B)] do
   if B[i][2]=1 or B[i][2]=-1 then
      cost:=(Length(B)-1)*(Length(cob[n][B[i][1]])-1);   # Markowitz cost
      if cost<bestcost then
         best:=i; bestcost:=cost;
         if cost=0 then break; fi;
      fi;
   fi;
od;
return best;
end;

#############################################################################
## cancel the pair (e,a), a = face number pos of the n-cell e
#############################################################################
Cancel:=function(n,e,pos)
local B,a,u,r,x,Bx,px,c,Y,y,By,py,t;
B:=bnd[n][e];
a:=B[pos][1]; u:=B[pos][2];
r:=Concatenation(B{[1..pos-1]},B{[pos+1..Length(B)]});      # d(e)=u*a+r
K:=K+1;
stepDeg[K]:=n; stepE[K]:=e; stepA[K]:=a; stepU[K]:=u;
stepBnd[K]:=List(r,t->[t[1],-u*t[2]]);                      # f(a) = -u*r

for x in cob[n][a] do
   if x<>e and alive[n+1][x] then
      Bx:=bnd[n][x];
      px:=FindFace(Bx,a);
      if px<>fail then
         c:=Bx[px][2];
         if not IsBound(T[n+1][x]) then T[n+1][x]:=[]; fi;
         Add(T[n+1][x],[K,c]);                             # g(x)=x-c*u*e
         Y:=List(r,t->[t[1],-c*u*t[2]]);                   # d'x = dx - c*u*de
         Bx:=Concatenation(Bx{[1..px-1]},Bx{[px+1..Length(Bx)]});
         bnd[n][x]:=Merge(Bx,Y);
         for t in Y do AddSet(cob[n][t[1]],x); od;
         if not inq[x] then Add(queue,x); inq[x]:=true; fi;
      fi;
   fi;
od;

if n<L then                      # remove e from the boundaries of (n+1)-cells
   for y in cob[n+1][e] do
      if alive[n+2][y] then
         By:=bnd[n+1][y];
         py:=FindFace(By,e);
         if py<>fail then
            bnd[n+1][y]:=Concatenation(By{[1..py-1]},By{[py+1..Length(By)]});
         fi;
      fi;
   od;
fi;

alive[n+1][e]:=false; alive[n][a]:=false;
bnd[n][e]:=[]; cob[n][a]:=[];
end;

#############################################################################
## the reduction, one degree at a time
#############################################################################
for n in [1..L] do
   queue:=[1..dimC[n+1]];
   SortBy(queue,k->Length(bnd[n][k]));      # cells with short boundaries first
   queue:=Reversed(queue);                  # (the queue is used as a stack)
   inq:=ListWithIdenticalEntries(dimC[n+1],true);
   while Length(queue)>0 do
      e:=Remove(queue); inq[e]:=false;
      if alive[n+1][e] then
         pos:=Pivot(n,e);
         if pos<>fail then Cancel(n,e,pos); fi;
      fi;
   od;
od;

#############################################################################
## the surviving generators and the boundary of D
#############################################################################
surv:=[]; newpos:=[]; dimD:=[];
for n in [0..L] do
   surv[n+1]:=Filtered([1..dimC[n+1]],k->alive[n+1][k]);
   newpos[n+1]:=ListWithIdenticalEntries(dimC[n+1],0);
   for i in [1..Length(surv[n+1])] do newpos[n+1][surv[n+1][i]]:=i; od;
   dimD[n+1]:=Length(surv[n+1]);
od;

Dbnd:=[];
for n in [1..L] do
   Dbnd[n]:=[];
   for j in [1..dimD[n+1]] do
      B:=bnd[n][surv[n+1][j]];
      b:=[];
      for t in B do
         if newpos[n][t[1]]=0 then
            Error("internal error: boundary meets a cancelled cell");
         fi;
         Add(b,[newpos[n][t[1]],t[2]]);
      od;
      Dbnd[n][j]:=b;
   od;
od;

DDimension:=function(n)
if n<0 or n>L then return 0; fi;
return dimD[n+1];
end;

if isSparse then
   DBoundary:=function(n,k)
   if n<1 or n>L then return []; fi;
   return Dbnd[n][k];
   end;
else
   DBoundary:=function(n,k)
   local w,t;
   if n<1 or n>L then return []; fi;
   w:=ListWithIdenticalEntries(dimD[n],0);
   for t in Dbnd[n][k] do w[t[1]]:=t[2]; od;
   return w;
   end;
fi;

props:=[["length",L],["type","chainComplex"],["characteristic",0]];
conn:=EvaluateProperty(C,"connected");
if conn<>fail then Add(props,["connected",conn]); fi;

if isSparse then
   D:=Objectify(HapSparseChainComplex,
        rec(dimension:=DDimension, boundary:=DBoundary, properties:=props));
else
   D:=Objectify(HapChainComplex,
        rec(dimension:=DDimension, boundary:=DBoundary, properties:=props));
fi;

#############################################################################
## F : C -> D.  Fimg[n+1][a] = image of a cancelled n-cell a (as a sparse
## vector in D_n); computed by running through the cancellations backwards.
## Surviving cells map to themselves, cells cancelled as "e" map to zero.
#############################################################################
Fimg:=List([0..L],n->[]);
for i in [K,K-1..1] do
   n:=stepDeg[i];
   acc:=[];
   for t in stepBnd[i] do
      j:=t[1];
      if alive[n][j] then
         Add(acc,[newpos[n][j],t[2]]);
      elif IsBound(Fimg[n][j]) then
         for s in Fimg[n][j] do Add(acc,[s[1],t[2]*s[2]]); od;
      fi;
   od;
   Fimg[n][stepA[i]]:=NormForm(acc);
od;

#############################################################################
## G : D -> C.  The image of a generator is computed by running through
## the cancellations of that degree in reverse order; coefficients on the
## cancelled "e" cells are pushed down through the recorded T-lists.
#############################################################################
Gvec:=function(n,s)
local Lam, pend, res, i, lam, x, t;
Lam:=[]; pend:=[]; res:=[[s,1]];
for t in T[n+1][s] do
   if not IsBound(Lam[t[1]]) then Lam[t[1]]:=0; AddSet(pend,t[1]); fi;
   Lam[t[1]]:=Lam[t[1]]-stepU[t[1]]*t[2];
od;
while Length(pend)>0 do
   i:=Remove(pend);                       # largest pending step
   lam:=Lam[i];
   if lam<>0 then
      x:=stepE[i];
      Add(res,[x,lam]);
      if IsBound(T[n+1][x]) then
         for t in T[n+1][x] do
            if not IsBound(Lam[t[1]]) then Lam[t[1]]:=0; AddSet(pend,t[1]); fi;
            Lam[t[1]]:=Lam[t[1]]-stepU[t[1]]*lam*t[2];
         od;
      fi;
   fi;
od;
return NormForm(res);
end;

Gimg:=List([0..L],n->[]);
for n in [0..L] do
   for j in [1..dimD[n+1]] do
      s:=surv[n+1][j];
      if not IsBound(T[n+1][s]) then T[n+1][s]:=[]; fi;
      Gimg[n+1][j]:=Gvec(n,s);
   od;
od;

mapF:=function(v,n)
local w,k,t,p;
if n<0 or n>L then return []; fi;
w:=ListWithIdenticalEntries(dimD[n+1],0);
for k in PositionsProperty(v,x->not IsZero(x)) do
   p:=newpos[n+1][k];
   if p>0 then
      w[p]:=w[p]+v[k];
   elif IsBound(Fimg[n+1][k]) then
      for t in Fimg[n+1][k] do w[t[1]]:=w[t[1]]+v[k]*t[2]; od;
   fi;
od;
return w;
end;

mapG:=function(v,n)
local w,j,t;
if n<0 or n>L then return []; fi;
w:=ListWithIdenticalEntries(dimC[n+1],0);
for j in PositionsProperty(v,x->not IsZero(x)) do
   for t in Gimg[n+1][j] do w[t[1]]:=w[t[1]]+v[j]*t[2]; od;
od;
return w;
end;

F:=Objectify(HapChainMap,
     rec(source:=C, target:=D, mapping:=mapF,
         properties:=[["characteristic",0],["type","chainMap"]]));
G:=Objectify(HapChainMap,
     rec(source:=D, target:=C, mapping:=mapG,
         properties:=[["characteristic",0],["type","chainMap"]]));

return [F,G];
end);
