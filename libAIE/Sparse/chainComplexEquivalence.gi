#############################################################################
##
## This function was produced by Claude AI in response to asking it to translate 
## the existing function ChainComplexEquivalentOfRegularCWComplex(Y) to a function
## ChainComplexEquivalenceOfChainComplex(C) that applies to an arbitrary chain 
## complex C of f.g. free abelian groups. The only difference is that coefficients
## in a boundary map can now be any integer rather than just -1, 0, 1. So when choosing 
## a pair (e^n,e^{n-1}) to eliminate we need to check additionally that e^{n-1} occurs
## with coefficient -1 or +1 and not -n or +n with n>1. Otherwise everything should be 
## the same.
##
##  ChainComplexEquivalenceOfChainComplex(C)
##
##  INPUT   C : a HAP chain complex of free abelian groups (characteristic 0),
##              dense (C!.boundary(n,k) is a vector) or sparse
##              (C!.boundary(n,k) is a list of pairs [i,coeff]).
##
##  OUTPUT  [F,G]  with F : C -> D and G : D -> C chain maps,
##                 F o G = identity on D  and  G o F ~ identity on C,
##          where D is the Morse complex of an acyclic matching on the basis
##          of C.  D has the same type (dense/sparse) as C.
##
##  METHOD  Discrete Morse theory, with the
##          matching found FIRST and the boundaries of the Morse complex
##          computed afterwards, lazily, by summing over gradient paths.
##          Nothing is ever rewritten.  This is the strategy of
##          ChainComplexEquivalenceOfRegularCWComplex(), generalised from
##          regular CW-complexes to arbitrary chain complexes.
##
##  1. THE MATCHING (greedy collapse, top dimension first).
##     The "support graph" of C has an edge e -> a whenever the basis
##     (n-1)-cell a occurs in the boundary of the basis n-cell e.  A cell is
##     *alive* until it is matched or declared critical.  An (n-1)-cell a
##     is a *free face* if exactly one alive n-cell e contains it.  If the
##     coefficient of a in d(e) is +1 or -1 we match (a,e) and kill both;
##     the other faces of e, and the faces of a, each lose one alive
##     coface, which may make them free.  When no free face with a unit
##     coefficient is left, the first alive cell of highest dimension is
##     declared critical (killed, not matched) and the process continues.
##     Only integer counters are kept: for each cell the NUMBER of alive
##     cofaces and the SUM of their indices (when the number is 1, the sum
##     is the index of the unique coface).
##
##     Acyclicity is automatic.  Number the pairs (a_1,e_1), (a_2,e_2), ...
##     in order of creation.  If a_j is a face of e_i with j<>i then j>i:
##     otherwise e_i, which is still alive when pair j is made, would be a
##     second coface of the free face a_j.  So along every gradient path the
##     pair numbers strictly increase and no cycle can occur.
##
##  2. THE MAPS (u_a = coefficient of a in d(e_a), +1 or -1).
##     F = phi : C -> D is given on a basis cell x of degree m by
##        phi(x) = x                                   x critical
##        phi(x) = 0                                   x the upper cell of a pair
##        phi(a) = - u_a * sum_{b<>a} [d(e_a):b] phi(b)   a the lower cell of a
##                                                     pair (a,e_a)
##     The boundary of D is  d_D(c) = phi(d c)  for critical c.
##     G : D -> C is  G(c) = c - sum_a u_a [d c:a] L(e_a)  (a over the faces
##     of c that are lower cells of pairs), where
##        L(e) = e - sum_{a' <> a_e} u_{a'} [d e:a'] L(e_{a'}).
##     G is evaluated by processing the upper cells reached in order of
##     their pair numbers, so every cell is visited once.
##     Everything is computed on demand and remembered (phi, d_D, G).
##
##  NOTES   D is NOT reduced further: it may still contain unit boundary
##          coefficients, and when a complex has no free faces (e.g. the
##          dense boundary matrices of a group resolution) almost nothing is
##          matched.  The method is meant for complexes that are close to
##          chain complexes of regular CW-complexes.
##
#############################################################################
InstallGlobalFunction(ChainComplexEquivalenceOfChainComplex,
function(C)
local
  L, isSparse, dimC, bnd, cnt, sumcof, status, stack, partner, unitc,
  partnerDown, ptime, K, ptr, crit, newpos, dimD, Phi, GImgRec, Wt, Stamp,
  stampval, DBndRec,
  NormForm, Free, RemoveCritical, EvalPhi, GImg, DBndSparse,
  n, m, k, b, t, c, found, props, conn, DDimension, DBoundary, mapF, mapG,
  D, F, G, w, s, j, qhead;

if not (IsHapChainComplex(C) or IsHapSparseChainComplex(C)) then
   Error("the argument must be a HAP chain complex");
fi;
if not EvaluateProperty(C,"characteristic")=0 then
   Error("only chain complexes of free abelian groups (characteristic 0) are handled");
fi;

L:=Length(C);
isSparse:=IsHapSparseChainComplex(C);
dimC:=List([0..L],n->C!.dimension(n));          # dimC[n+1] = rank of C_n

#############################################################################
## sorted sparse vectors [[index,coeff],...] without zero coefficients
#############################################################################
NormForm:=function(b)
local S, out, x, mm;
S:=ShallowCopy(b);
Sort(S);
out:=[];
for x in S do
   mm:=Length(out);
   if mm>0 and out[mm][1]=x[1] then
      out[mm]:=[x[1],out[mm][2]+x[2]];
   else
      Add(out,[x[1],x[2]]);
   fi;
od;
return Filtered(out,x->x[2]<>0);
end;

#############################################################################
## sparse boundaries bnd[n][k]; counters cnt[n][a], sumcof[n][a] over the
## n-cells containing the (n-1)-cell a
#############################################################################
bnd:=[]; cnt:=[]; sumcof:=[];
for n in [1..L] do
   bnd[n]:=[];
   cnt[n]:=ListWithIdenticalEntries(dimC[n],0);
   sumcof[n]:=ListWithIdenticalEntries(dimC[n],0);
   for k in [1..dimC[n+1]] do
      b:=C!.boundary(n,k);
      if isSparse then
         b:=NormForm(b);
      else
         b:=List(PositionsProperty(b,x->not IsZero(x)),i->[i,b[i]]);
      fi;
      bnd[n][k]:=b;
      for t in b do
         cnt[n][t[1]]:=cnt[n][t[1]]+1;
         sumcof[n][t[1]]:=sumcof[n][t[1]]+k;
      od;
   od;
od;

# status[n+1][x] for an n-cell x:  0 alive, 1 critical, 2 lower cell of a
# pair, 3 upper cell of a pair
status:=List([0..L],n->ListWithIdenticalEntries(dimC[n+1],0));
partner:=List([0..L],n->ListWithIdenticalEntries(dimC[n+1],0));      # lower -> upper
unitc:=List([0..L],n->ListWithIdenticalEntries(dimC[n+1],0));        # lower -> u
partnerDown:=List([0..L],n->ListWithIdenticalEntries(dimC[n+1],0));  # upper -> lower
ptime:=List([0..L],n->ListWithIdenticalEntries(dimC[n+1],0));        # upper -> pair number
K:=0;
stack:=[]; qhead:=ListWithIdenticalEntries(L,1);
for n in [1..L] do
   stack[n]:=Filtered([1..dimC[n]],a->cnt[n][a]=1);
od;

#############################################################################
## process the free faces of dimension n-1
#############################################################################
Free:=function(n)
local a,e,B,pos,c,t,b;
while qhead[n]<=Length(stack[n]) do
   a:=stack[n][qhead[n]]; qhead[n]:=qhead[n]+1;
   if status[n][a]=0 and cnt[n][a]=1 then
      e:=sumcof[n][a];
      B:=bnd[n][e];
      pos:=PositionSorted(B,[a]);
      c:=B[pos][2];
      if c=1 or c=-1 then
         K:=K+1;
         status[n][a]:=2; partner[n][a]:=e; unitc[n][a]:=c;
         status[n+1][e]:=3; partnerDown[n+1][e]:=a; ptime[n+1][e]:=K;
         for t in B do                              # faces of e lose e
            b:=t[1];
            cnt[n][b]:=cnt[n][b]-1; sumcof[n][b]:=sumcof[n][b]-e;
            if cnt[n][b]=1 and status[n][b]=0 then Add(stack[n],b); fi;
         od;
         if n>=2 then                               # faces of a lose a
            for t in bnd[n-1][a] do
               b:=t[1];
               cnt[n-1][b]:=cnt[n-1][b]-1; sumcof[n-1][b]:=sumcof[n-1][b]-a;
               if cnt[n-1][b]=1 and status[n-1][b]=0 then Add(stack[n-1],b); fi;
            od;
         fi;
      fi;
   fi;
od;
stack[n]:=[]; qhead[n]:=1;
end;

RemoveCritical:=function(n,c)       # declare the n-cell c critical
local t,b;
status[n+1][c]:=1;
if n>=1 then
   for t in bnd[n][c] do
      b:=t[1];
      cnt[n][b]:=cnt[n][b]-1; sumcof[n][b]:=sumcof[n][b]-c;
      if cnt[n][b]=1 and status[n][b]=0 then Add(stack[n],b); fi;
   od;
fi;
end;

#############################################################################
## the matching
#############################################################################
ptr:=ListWithIdenticalEntries(L+1,1);
for n in [L,L-1..1] do Free(n); od;
repeat
   found:=false;
   for n in [L,L-1..0] do
      while ptr[n+1]<=dimC[n+1] and status[n+1][ptr[n+1]]<>0 do
         ptr[n+1]:=ptr[n+1]+1;
      od;
      if ptr[n+1]<=dimC[n+1] then
         found:=true;
         RemoveCritical(n,ptr[n+1]);
         break;
      fi;
   od;
   if found then
      for m in [L,L-1..1] do Free(m); od;
   fi;
until not found;

crit:=[]; newpos:=[]; dimD:=[];
for n in [0..L] do
   crit[n+1]:=Filtered([1..dimC[n+1]],k->status[n+1][k]=1);
   newpos[n+1]:=ListWithIdenticalEntries(dimC[n+1],0);
   for j in [1..Length(crit[n+1])] do newpos[n+1][crit[n+1][j]]:=j; od;
   dimD[n+1]:=Length(crit[n+1]);
   if ForAny(status[n+1],x->x=0) then Error("internal error: undecided cell"); fi;
od;

#############################################################################
## phi : C -> D on a basis m-cell x (sparse vector in D_m), memoised.
## For a lower cell the recursion is run with an explicit stack.
#############################################################################
Phi:=List([0..L],n->[]);
EvalPhi:=function(m,x)
local st,y,e,B,t,b,ready,acc,u,coef,r;
if status[m+1][x]=1 then return [[newpos[m+1][x],1]]; fi;
if status[m+1][x]=3 then return []; fi;
if IsBound(Phi[m+1][x]) then return Phi[m+1][x]; fi;
st:=[x];
while Length(st)>0 do
   y:=st[Length(st)];
   if IsBound(Phi[m+1][y]) then
      Remove(st);
   else
      e:=partner[m+1][y];
      B:=bnd[m+1][e];
      ready:=true;
      for t in B do
         b:=t[1];
         if b<>y and status[m+1][b]=2 and not IsBound(Phi[m+1][b]) then
            Add(st,b); ready:=false;
         fi;
      od;
      if ready then
         u:=unitc[m+1][y];
         acc:=[];
         for t in B do
            b:=t[1];
            if b<>y then
               coef:=-u*t[2];
               if status[m+1][b]=1 then
                  Add(acc,[newpos[m+1][b],coef]);
               elif status[m+1][b]=2 then
                  for r in Phi[m+1][b] do Add(acc,[r[1],coef*r[2]]); od;
               fi;
            fi;
         od;
         Phi[m+1][y]:=NormForm(acc);
         Remove(st);
      fi;
   fi;
od;
return Phi[m+1][x];
end;

#############################################################################
## boundary of the Morse complex (sparse, memoised): d_D(c) = phi(d c)
#############################################################################
DBndRec:=List([0..L],n->[]);
DBndSparse:=function(n,j)
local c,acc,t,b,r;
if IsBound(DBndRec[n+1][j]) then return DBndRec[n+1][j]; fi;
c:=crit[n+1][j];
acc:=[];
for t in bnd[n][c] do
   b:=t[1];
   if status[n][b]=1 then
      Add(acc,[newpos[n][b],t[2]]);
   elif status[n][b]=2 then
      for r in EvalPhi(n-1,b) do Add(acc,[r[1],t[2]*r[2]]); od;
   fi;
od;
DBndRec[n+1][j]:=NormForm(acc);
return DBndRec[n+1][j];
end;

DDimension:=function(n)
if n<0 or n>L then return 0; fi;
return dimD[n+1];
end;

if isSparse then
   DBoundary:=function(n,k)
   if n<1 or n>L then return []; fi;
   return DBndSparse(n,k);
   end;
else
   DBoundary:=function(n,k)
   local w,t;
   if n<1 or n>L then return []; fi;
   w:=ListWithIdenticalEntries(dimD[n],0);
   for t in DBndSparse(n,k) do w[t[1]]:=t[2]; od;
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
## G : D -> C on a critical n-cell number j (sparse vector in C_n), memoised.
## The upper cells reachable from c along gradient paths are collected, then
## processed in order of their pair numbers, accumulating path weights.
#############################################################################
GImgRec:=List([0..L],n->[]);
Wt:=[]; Stamp:=[]; stampval:=0;
GImg:=function(n,j)
local c,cone,st,e,e2,ae,a,t,res,wt,coef;
if IsBound(GImgRec[n+1][j]) then return GImgRec[n+1][j]; fi;
c:=crit[n+1][j];
res:=[[c,1]];
if n>=1 then
   if not IsBound(Wt[n+1]) then
      Wt[n+1]:=ListWithIdenticalEntries(dimC[n+1],0);
      Stamp[n+1]:=ListWithIdenticalEntries(dimC[n+1],0);
   fi;
   stampval:=stampval+1;
   cone:=[]; st:=[];
   for t in bnd[n][c] do
      a:=t[1];
      if status[n][a]=2 then
         e:=partner[n][a];
         Wt[n+1][e]:=Wt[n+1][e]-unitc[n][a]*t[2];
         if Stamp[n+1][e]<>stampval then
            Stamp[n+1][e]:=stampval; Add(cone,e); Add(st,e);
         fi;
      fi;
   od;
   while Length(st)>0 do
      e:=Remove(st);
      ae:=partnerDown[n+1][e];
      for t in bnd[n][e] do
         a:=t[1];
         if a<>ae and status[n][a]=2 then
            e2:=partner[n][a];
            if Stamp[n+1][e2]<>stampval then
               Stamp[n+1][e2]:=stampval; Add(cone,e2); Add(st,e2);
            fi;
         fi;
      od;
   od;
   SortBy(cone,e->ptime[n+1][e]);
   for e in cone do
      wt:=Wt[n+1][e]; Wt[n+1][e]:=0;
      if wt<>0 then
         Add(res,[e,wt]);
         ae:=partnerDown[n+1][e];
         for t in bnd[n][e] do
            a:=t[1];
            if a<>ae and status[n][a]=2 then
               e2:=partner[n][a];
               Wt[n+1][e2]:=Wt[n+1][e2]-wt*unitc[n][a]*t[2];
            fi;
         od;
      fi;
   od;
fi;
GImgRec[n+1][j]:=res;
return res;
end;

mapF:=function(v,n)
local w,k,t,s;
if n<0 or n>L then return []; fi;
w:=ListWithIdenticalEntries(dimD[n+1],0);
for k in PositionsProperty(v,x->not IsZero(x)) do
   s:=status[n+1][k];
   if s=1 then
      w[newpos[n+1][k]]:=w[newpos[n+1][k]]+v[k];
   elif s=2 then
      for t in EvalPhi(n,k) do w[t[1]]:=w[t[1]]+v[k]*t[2]; od;
   fi;
od;
return w;
end;

mapG:=function(v,n)
local w,j,t;
if n<0 or n>L then return []; fi;
w:=ListWithIdenticalEntries(dimC[n+1],0);
for j in PositionsProperty(v,x->not IsZero(x)) do
   for t in GImg(n,j) do w[t[1]]:=w[t[1]]+v[j]*t[2]; od;
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
