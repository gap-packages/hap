##################################################
##################################################
InstallMethod(VertexStar,
"star of vertex in a simplicial complex",
[IsHapSimplicialComplex,IsInt],
function(K,v)
local L, S;
L:=Concatenation(K!.simplicesLst);
S:=Filtered(L, x->v in x);
return MaximalSimplicesToSimplicialComplex(S);
end);
##################################################
##################################################

##################################################
##################################################
InstallMethod(VertexLink,
"link of a vertex in a simplicial complex",
[IsHapSimplicialComplex,IsInt],
function(K,v)
local S,L;
S:=VertexStar(K,v);
L:=Concatenation(S!.simplicesLst);
L:=Filtered(L, x-> not v in x);
return MaximalSimplicesToSimplicialComplex(L);
end);
##################################################
##################################################


#####################################################################
#####################################################################
InstallOtherMethod(Dimension,
"Dimension of  simplicial complex",
[IsHapSimplicialComplex],
function(K) return EvaluateProperty(K,"dimension");
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Size,
"Size of  simplicial complex",
[IsHapSimplicialComplex],
function(K) local n,s;
s:=0;
for n in [0..Dimension(K)] do
s:=s+K!.nrSimplices(n);
od;
return s;
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallOtherMethod(EulerCharacteristic,
"Euler characteristic  of a simplicial complex",
[IsHapSimplicialComplex],
function(M);
return
Sum(List([0..Dimension(M)],i->((-1)^i)*M!.nrSimplices(i)));
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Homology,
"Integral homologies  of a simplicial complex",
[IsHapSimplicialComplex],
function(M);
return
List([0..Dimension(M)],n->Homology(M,n));
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Homology,
"Integral homologies  of a simplicial complex",
[IsHapSimplicialComplex,IsInt],
function(M,n) local P, PM, i, answer,C,A;
return Homology(RegularCWComplex(M),n);
#if Dimension(M)=2 and n=1 then
#return FirstHomologySimplicialTwoComplex(M);
#fi;
#P:=PathComponentsOfSimplicialComplex(M,0);
#if n=0 then return [1..P]*0;fi;
#answer:=[];
#for i in [1..P] do
#PM:=PathComponentsOfSimplicialComplex(M,i);
#PM:=SkeletonOfSimplicialComplex(PM,n+1);
#A:=ContractibleSubcomplexOfSimplicialComplex(PM);;
#C:=ChainComplexOfSimplicialPair(PM,A);
#answer:=Concatenation(answer,Homology(C,n));
#od;
#return answer;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Homology,
"Integral homologies  of a simplicial complex",
[IsHapSimplicialComplex,IsInt,IsInt],
function(M,n,p) local P, PM, i, answer,C,A;
#P:=PathComponentsOfGraph(GraphOfSimplicialComplex(M),0);
P:=PathComponentsOfSimplicialComplex(M,0);
if n=0 then return P;fi;
answer:=0;
for i in [1..P] do
PM:=PathComponentsOfSimplicialComplex(M,i);
PM:=SkeletonOfSimplicialComplex(PM,n+1);
A:=ContractibleSubcomplexOfSimplicialComplex(PM);;
C:=ChainComplexOfSimplicialPair(PM,A);
answer:=answer+Homology(TensorWithIntegersModP(C,p),n);
od;
return answer;
end);
#####################################################################
#####################################################################



#####################################################################
#####################################################################
InstallOtherMethod(Bettinumbers,
"Integral homologies  of a simplicial complex",
[IsHapSimplicialComplex],
function(M)
local H;
H:=Homology(M);
return List(H,x->Length(Filtered(x,a->a=0)));
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallOtherMethod(ChainComplex,
"Cellular chain complex of a simplicial complex",
[IsHapSimplicialComplex],
function(M)
return ChainComplexOfSimplicialComplex(M);;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(ChainComplexOfSimplicialComplex,
function(M)
local
        Dimsn,Boundary,nullvecs,i;

#############################
Dimsn:=function(n);
return M!.nrSimplices(n);
end;
#############################

nullvecs:=[];
for i in [0..Dimension(M)+1] do
nullvecs[i+1]:=List([1..M!.nrSimplices(i)],a->0);
od;

#############################
Boundary:=function(n,i)
local V,Vhat, j, bnd;
V:=M!.simplices(n,i);

#bnd:=List([1..M!.nrSimplices(n-1)],a->0);
bnd:=StructuralCopy(nullvecs[n]);

for j in [1..Length(V)] do
Vhat:=StructuralCopy(V);
RemoveSet(Vhat,V[j]);
if IsOddInt(j) then
bnd[M!.enumeratedSimplex(Vhat)]:=1;
else
bnd[M!.enumeratedSimplex(Vhat)]:=-1;
fi;
od;

return bnd;
end;
#############################


return
Objectify(HapChainComplex,
           rec(
           dimension:=Dimsn,
           boundary:=Boundary,
           properties:=[
           ["length",EvaluateProperty(M,"dimension")],
           ["type","chainComplex"],
           ["characteristic",0]]
           ));

end);
#################################################################
#################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(SparseChainComplexOfSimplicialComplex,
function(M)
local
        Dimsn,Boundary,i;

#############################
Dimsn:=function(n);
return M!.nrSimplices(n);
end;
#############################


#############################
Boundary:=function(n,i)
local V,Vhat, j, bnd;
V:=M!.simplices(n,i);

bnd:=[];

for j in [1..Length(V)] do
Vhat:=StructuralCopy(V);
RemoveSet(Vhat,V[j]);
if IsOddInt(j) then
Add(bnd,[M!.enumeratedSimplex(Vhat),1]);
else
Add(bnd,[M!.enumeratedSimplex(Vhat),-1]);
fi;
od;

return bnd;
end;
#############################


return
Objectify(HapSparseChainComplex,
           rec(
           dimension:=Dimsn,
           boundary:=Boundary,
           properties:=[
           ["length",EvaluateProperty(M,"dimension")],
           ["type","chainComplex"],
           ["characteristic",0]]
           ));

end);
#################################################################
#################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(SparseFilteredChainComplexOfFilteredSimplicialComplex,
function(M)
local
        Dimsn,Boundary,i;

#############################
Dimsn:=function(n);
return M!.nrSimplices(n);
end;
#############################


#############################
Boundary:=function(n,i)
local V,Vhat, j, bnd;
V:=M!.simplices(n,i);

bnd:=[];

for j in [1..Length(V)] do
Vhat:=StructuralCopy(V);
RemoveSet(Vhat,V[j]);
if IsOddInt(j) then
Add(bnd,[M!.enumeratedSimplex(Vhat),1]);
else
Add(bnd,[M!.enumeratedSimplex(Vhat),-1]);
fi;
od;

return bnd;
end;
#############################

return
Objectify(HapFilteredSparseChainComplex,
           rec(
           dimension:=Dimsn,
           boundary:=Boundary,
           filteredDimension:=M!.filteredDimension,
           properties:=[
           ["length",EvaluateProperty(M,"dimension")],
           ["type","chainComplex"],
           ["characteristic",0],
	   ["filtration_length",M!.filtrationLength]]
           ));

end);
#################################################################
#################################################################



#################################################################
#################################################################
InstallGlobalFunction(SimplicesToSimplicialComplex,
function(L)
local
        Vertices,
        NrSimplices,
        SimplicesLst,
        Simplices,
        EnumeratedSimplex,
        dim,
        n,v,x;

for n in [1..Length(L)] do
Apply(L[n],x->SSortedList(x));
od;
Apply(L,l->SSortedList(l));
dim:=PositionProperty(L,x->Size(x)=0);
if dim=fail then dim:=Length(L)-1; else dim:=dim-2; fi;
if Length(L)>0 then
Vertices:=Concatenation(L[1]);
else Vertices:=[];
fi;
Vertices:=SSortedList(Vertices);
SimplicesLst:=L;

#####################
NrSimplices:=function(n);
if n<0 or n>dim then return 0; fi;
return Length(SimplicesLst[n+1]);
end;
#####################

#####################
Simplices:=function(n,i);
return SimplicesLst[n+1][i];
end;
#####################


#####################
EnumeratedSimplex:=function(v);
return PositionSet(SimplicesLst[Length(v)],v);
end;
#####################


return
Objectify(HapSimplicialComplex,
           rec(
           vertices:=Vertices,
           nrSimplices:=NrSimplices,
           simplices:=Simplices,
           simplicesLst:=SimplicesLst,
           enumeratedSimplex:=EnumeratedSimplex,
           properties:=[
           ["dimension",dim]]
           ));

end);
#################################################################
#################################################################



#################################################################
#################################################################
InstallGlobalFunction(MaximalSimplicesToSimplicialComplex,
function(V)
local
	Vertices,
        NrSimplices,
	SimplicesLst,
	Simplices,
	EnumeratedSimplex,
	dim,
	n,v,x;

dim:=Maximum(List(V,v->Length(v)))-1;
Vertices:=Concatenation(V);
Vertices:=SSortedList(Vertices);
SimplicesLst:=List([1..dim+100],i->[]); #SLOPPY!!!
SimplicesLst[1]:=List(Vertices,x->[x]);

for v in V do
for n in [1..Length(v)] do
for x in Combinations(v,n) do
Add(SimplicesLst[n],x);
od;
od;
od;

for n in [1..Length(v)] do
SimplicesLst[n]:=SSortedList(SimplicesLst[n]);
od;

return SimplicesToSimplicialComplex(SimplicesLst);

end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(SkeletonOfSimplicialComplex,
function(S,k)
local
	T,p,prp;
T:=Objectify(HapSimplicialComplex, rec());
T!.vertices:=StructuralCopy(S!.vertices);
T!.nrSimplices:=StructuralCopy(S!.nrSimplices);
T!.simplices:=StructuralCopy(S!.simplices);
if IsBound(S!.simplicesLst) then
T!.simplicesLst:=StructuralCopy(S!.simplicesLst);
fi;
T!.enumeratedSimplex:=StructuralCopy(S!.enumeratedSimplex);
prp:=StructuralCopy(S!.properties);
p:=PositionProperty(prp,x->"dimension" in x);
prp[p][2]:=k;
T!.properties:=prp;

return T;
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(SimplicialMapNC,
function(K,L,F)
local
        M,SimpF;
#Assume K, L to be simplicial complexes and F to
#be a mapping of vertex sets. SimpF will be the
#mapping on simplexes.

###################################
SimpF:=function(s);
return List(s,F);
end;
###################################

M:=Objectify(HapSimplicialMap, rec());
M!.source:=K;
M!.target:=L;
M!.mapping:=SimpF;
M!.properties:=[];

return M;
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(SimplicialMap,
function(K,L,F)
local
        M,Lsimps,s,d,dim,l,Fs;

Lsimps:=[];
for d in [0..Dimension(K)] do
Lsimps[d+1]:=List([1..L!.nrSimplices(d)],i->SSortedList(L!.simplices(d,i)));
Lsimps[d+1]:=SSortedList(Lsimps[d+1]);
od;

for d in [0..Dimension(K)] do
for s in [1..K!.nrSimplices(d)] do
Fs:=List(K!.simplices(d,s),F);
Fs:=SSortedList(Fs);
if not Fs in Lsimps[Length(Fs)] then return fail; fi;
od;
od;

return SimplicialMapNC(K,L,F);
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(ChainMapOfSimplicialMap,
function(F)
	local K,L,Lsimps,zeros,CK,CL, map, d;

K:=Source(F);;
L:=Target(F);;
CK:=ChainComplex(K);
CL:=ChainComplex(L);


Lsimps:=[];
zeros:=[];
for d in [0..Dimension(L)] do
Lsimps[d+1]:=List([1..L!.nrSimplices(d)],i->SSortedList(L!.simplices(d,i)));
Lsimps[d+1]:=SSortedList(Lsimps[d+1]);
zeros[d+1]:=[1..L!.nrSimplices(d)]*0;
od;


###########################
map:=function(v,d)
local i,z,s;

z:=StructuralCopy(zeros[d+1]);

for i in [1..Length(v)] do
if not v[i]=0 then 
s:=SSortedList( F!.mapping( K!.simplices(d,i) ));
s:=Position(Lsimps[d+1],s);
z[s]:=z[s]+v[i];
fi;
od;

return z;
end;
###########################

return
Objectify(HapChainMap,
           rec(
           source:=CK,
           target:=CL,
           mapping:=map,
           properties:=[
           ["type","chainMap"],
           ["characteristic",0]]
           ));


end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(MaximalSimplicesOfSimplicialComplex,
function(K)
local
	MaxSimps, Simps, genSimps,dim,s,x,i;

dim:=EvaluateProperty(K,"dimension");
Simps:=1*K!.simplicesLst[dim+1];
MaxSimps:=Simps;

while dim>0 do
dim:=dim-1;
   genSimps:=[];
   for s in Simps do
   for x in s do
   Add(genSimps,SSortedList(Filtered(s,i->not i=x)));
   od;
   od;
genSimps:=SSortedList(genSimps);
Simps:=genSimps;
   for i in [1..K!.nrSimplices(dim)] do
   if not K!.simplices(dim,i) in genSimps then 
      Add(Simps,K!.simplices(dim,i));
      Add(MaxSimps,K!.simplices(dim,i)); 
   fi;
   od;
od;

return MaxSimps;
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(ContractSimplicialComplex,
function(K)
local M, i, j, k, mx, KK, MM, C, inv,A;

if IsList(K) then M:=1*K;
else M:=MaximalSimplicesOfSimplicialComplex(K);;
fi;

A:=NullMat(Length(M),Length(M));
for i in [1..Length(M)] do
   for j in [i+1..Length(M)] do
      A[i][j]:=Size(Intersection(M[i],M[j]));
      A[j][i]:=A[i][j];
   od;
od;
inv:=function(x); return x{[1..Length(x)-1]}; end;
MM:=[];
for i in [1..Length(M)] do
   mx:=0;
   for j in [1..Length(M)] do
      if not i=j then
         mx:= Maximum(A[i][j],mx);
      fi;
   od;
   mx:=Minimum(2+mx,Length(M[i]));
   C:=Combinations(M[i],mx);
   C:=Classify(C,inv);
   Apply(C,x->x[1]);
   Append(MM,C);
od;
MM:=SSortedList(MM);
KK:=MaximalSimplicesToSimplicialComplex(MM);
return KK;
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(ContractSimplicialComplex_alt,
function(K)
local RemoveApex,RemoveEdge,Vertices,G,A,LA,Faces,bool,dim,PC,iter,row,x,i,j,k;
#THIS FUNCTION COULD BE DELETED
G:=GraphOfSimplicialComplex(K);
Faces:=K!.simplicesLst;
dim:=PositionProperty(Faces,x->Size(x)=0);
if IsInt(dim) then dim:=dim-1; else dim:=Dimension(K); fi;
A:=G!.incidenceMatrix;
#PC:=PathComponentsOfGraph(G,0);
PC:=PathComponentsOfSimplicialComplex(K,0);
bool:=List(Faces,Length);;
Vertices:=StructuralCopy(K!.vertices);

################################################
##
RemoveApex:=function(i)
local x, y, cols, lc,D;

cols:=[];
for x in LA do
if not A[i][x]=0 then Add(cols,x); fi;
od;

if Length(cols)=0 or Length(cols)>dim then  return false; fi;

cols:=K!.vertices{cols};
AddSet(cols,K!.vertices[i]);
if not cols in Faces[Length(cols)] then return false; fi;

for x in LA do
A[i][x]:=0;
A[x][i]:=0;
od;

for x in [1..dim] do #####?
Faces[x]:=Filtered(Faces[x], a->not K!.vertices[i] in a);
od;

RemoveSet(Vertices,K!.vertices[i]);

return true;
end;
##
################################################

################################################
##
RemoveEdge:=function(i,j)
local cols, x, y,lc;


cols:=[];
for x in LA do
if (not A[i][x]=0) and (not A[j][x]=0) then Add(cols,x); fi;
od;

if Length(cols)=0  or Length(cols)>dim-2 then return false; fi;

cols:=K!.vertices{cols};
AddSet(cols,K!.vertices[i]);
AddSet(cols,K!.vertices[j]);
if not IsBound(Faces[Length(cols)]) then return false; fi;
if not cols in Faces[Length(cols)] then return false; fi;

A[i][j]:=0;
A[j][i]:=0;

for x in [2..dim+1] do #####?
Faces[x]:=Filtered(Faces[x], a->
not SSortedList([K!.vertices[i], K!.vertices[j]]) in a);
od;


return true;
end;
##
################################################


LA:=Filtered([1..Length(A)], i->Sum(A[i])>0);
iter:=Filtered(LA,i->Sum(A[i])<dim+1);

x:=true;
while x do
x:=false;
   for i in iter do
   if RemoveApex(i) then x:=true; fi;
   od;
od;

#K!.vertices:=Vertices;
#if bool=Sum(Sum(A)) then return false; fi;;;
#return true;

LA:=Filtered([1..Length(A)], i->Sum(A[i])>0);
   for i in [1..Length(LA)] do
   for j in [i+1..Length(LA)] do
   if RemoveEdge(LA[i],LA[j]) then x:=true;  fi;
   od;
   od;

if bool=List(Faces,Length) then return false; fi;;;
K!.vertices:=Vertices;
ContractSimplicialComplex(K);

return true;
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(ContractGraph,
function(G)
local RemoveApex,RemoveEdge,A,LA,bool,PC,row,x,i,j,k,RemoveNonCritical1And2Cells;

A:=G!.incidenceMatrix;
PC:=PathComponentsOfGraph(G,0);
Unbind(G!.pathComponents);
bool:=Sum(Sum(A));


###############################################
RemoveNonCritical1And2Cells:=function(A)
local bool, i,j,L;

bool:=true;
while bool do
bool:=false;
L:=Filtered([1..Length(A)],i->not Sum(A[i])=1);
if Length(L)<Length(A) then bool:=true; 
A:=A{L};
A:=List(A,row->row{L});
fi;

L:=[];
for i in [1..Length(A)-1] do
for j in [i+1..Length(A)] do
if A[i][j]=1 then
if A[i]*A[j]=1 then A[i][j]:=0; A[j][i]:=0; bool:=true;  fi;
fi;
od;
od;
od;

end;
#################################################


################################################
##
RemoveApex:=function(i)
local x, y, cols, lc,D;

cols:=[];
for x in LA do
if not A[i][x]=0 then Add(cols,x); fi;
od;

if Length(cols)=0 then return false; fi;

lc:=Length(cols);
for x in [1..lc] do
for y in [x+1..lc] do
if A[cols[x]][cols[y]]=0  then return false; fi;
od;
od;

for x in LA do
A[i][x]:=0;
A[x][i]:=0;
od;
return true;
end;
##
################################################

################################################
##
RemoveEdge:=function(i,j)
local cols, x, y,lc;

if A[i][j]=0 then return false; fi;

cols:=[];
for x in LA do
if (not A[i][x]=0) and (not A[j][x]=0) then Add(cols,x); fi;
od;

if Length(cols)=0  then return false; fi;

lc:=Length(cols);
for x in [1..lc] do
for y in [x+1..lc] do
if A[cols[x]][cols[y]]=0  then return false; fi;
od;
od;

A[i][j]:=0;
A[j][i]:=0;

return true;
end;
##
################################################ 

LA:=Filtered([1..Length(A)], i->Sum(A[i])>0);

RemoveNonCritical1And2Cells(A);

x:=true;
while x do
x:=false;
   for i in LA do
   if RemoveApex(i) then x:=true; fi;
   od;
od;

x:=true;
while x do
x:=false;
LA:=Filtered([1..Length(A)], i->Sum(A[i])>0);
   for i in [1..Length(LA)] do
   for j in [i+1..Length(LA)] do
   if RemoveEdge(LA[i],LA[j]) then x:=true;  fi;
   od;
   od;
od;

if bool=Sum(Sum(A)) then return false; fi;;;
ContractGraph(G);

#Now create the correct number of path components.
A:=Filtered(A,v->not IsZero(v));
A:=TransposedMatMutable(A);
A:=Filtered(A,v->not IsZero(v));
G!.incidenceMatrix:=A;
PC:=PC-PathComponentsOfGraph(G,0);

if Length(A)>0 then row:=0*A[1]; 
G!.singletons:=[];
for i in [1..PC] do
Add(A,0*row);
Add(G!.singletons,Length(A));
od;
A:=TransposedMatMutable(A);
for i in [1..PC] do
Add(A,0*A[1]);
od;
fi;

if  Length(A)=0 then row:=[0]; 
G!.singletons:=[];
for i in [1..PC] do
Add(A,0*row);
Add(G!.singletons,Length(A));
od;
A:=TransposedMatMutable(A);
for i in [1..PC-1] do #WHY PC-1???
Add(A,0*A[1]);
od;
fi;

G!.incidenceMatrix:=A;
Unbind(G!.pathComponents);
Unbind(G!.bettiZero);
k:=PositionProperty(G!.properties,x->"numberofvertices" in x);
G!.properties[k][2]:=Length(A);

return true;
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(RipsHomology,
function(arg)
local P,G,n,p,bool,PC,answer,i ;

G:=StructuralCopy(arg[1]);
n:=arg[2];
if Length(arg)>2 then p:=arg[3]; bool:=true; answer:=0;
 else bool:=false; answer:=[]; fi;

PC:=PathComponentsOfGraph(G,0);;

if n=0 and bool then return PC; fi;
if n=0 and not bool then return [1..PC]*0; fi;

for i in [1..PC] do
P:=PathComponentsOfGraph(G,i);
ContractGraph(P);

if bool then
  if n=1 then answer:=answer+Homology(TensorWithIntegersModP(
                    RipsChainComplex(P,n),p),n);
  else
     answer:=answer+Homology(SimplicialNerveOfGraph(P,n+1),n,p);
  fi;
else
  if n=1 then answer:=Concatenation(answer,Homology( RipsChainComplex(P,n),n));
  else 
  answer:=Concatenation(answer,Homology(SimplicialNerveOfGraph(P,n+1),n));
  fi;
fi;
od;

return answer;
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(ContractibleSubcomplexOfSimplicialComplex,
function(KK)
local L,K,Simplices, NrSimplices, EnumeratedSimplex,
      RemoveSimplexIfPoss, Link, RemovableSimplex, 
      Leaves,  G, HH,   
      Faces,dim,bool,d,k,s,t,tt,x;;

K:=PathComponentsOfSimplicialComplex(KK,1);
dim:=Dimension(K);
G:=GraphOfSimplicialComplex(K);
G:=G!.incidenceMatrix;

################
if IsBound(K!.simplicesLst) then
Faces:=StructuralCopy(K!.simplicesLst);
else
Faces:=[];
for d in [0..dim] do
Faces[d+1]:=[];
for k in [1..K!.nrSimplices(d)] do
Faces[d+1][k]:=K!.simplices(d,k);
od;
Faces[d+1]:=SSortedList(Faces[d+1]);
od;
fi;
################ 

#####################
RemoveSimplexIfPoss:=function(s)
local i,t,cnt,tt;
cnt:=0;
       for i in s do
       t:=StructuralCopy(s); RemoveSet(t,i);
       if not t in L[d-1] then cnt:=cnt+1; tt:=t; fi;
       if cnt>1 then break; fi;  
       od;

       if cnt=0 then Add(HH,s); else 
       if cnt=1 then Add(L[d],s); Add(HH,s); 
       Add(L[d-1],tt);
       return true;
       fi;
       fi;
return false;
end;
#####################


L:=[];
if Length(Faces[1])=0 then return SimplicesToSimplicialComplex(L); fi;

L[1]:=[Faces[1][1]];

for d in [2..dim+1] do
   L[d]:=[];
   HH:=[];
   bool:=true;
   while bool do
      bool:=false;
      Faces[d]:=Difference(Faces[d],HH);
      HH:=[];
      for s in Faces[d]  do 
         if RemoveSimplexIfPoss(s) then bool:=true;  fi;
      od;
   od;
od;

for d in [dim+2..dim+100] do ##SLOPPY!
L[d]:=[];
od;

return SimplicesToSimplicialComplex(L);

end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(GraphOfSimplicialComplex,
function(K)
local
	G,dim,V,d,e,x;


dim:=Length(K!.vertices);
V:=[];
for x in [1..dim] do
V[K!.vertices[x]]:=x;
od;
G:=NullMat(dim,dim);

for e in [1..K!.nrSimplices(1)] do
x:=K!.simplices(1,e);
G[V[x[1]]][V[x[2]]]:=1;
G[V[x[2]]][V[x[1]]]:=1;
od;

G:= IncidenceMatrixToGraph(G);
if IsBound(K!.clustersizes) then
G!.clustersizes:=K!.clustersizes;
fi;
return G;

end);
#################################################################
#################################################################

#############################################
#############################################
InstallGlobalFunction(PathComponentsOfSimplicialComplex,
function(KK,N)
local

      K,Colours, clr,leaves,i,v,e,newleaves,Edges, EDGES, L,nr;

K:=IntegerSimplicialComplex(KK);
L:=Maximum(Flat(K!.simplicesLst[1]));
Colours:=List([1..L],i->0);
EDGES:=List([1..L],i->[]);

for e in K!.simplicesLst[2] do
  Add(EDGES[e[1]],e[2]);
  Add(EDGES[e[2]],e[1]);
od;
clr:=0;

while Length(Filtered(Colours,c->c=0)) > L-KK!.nrSimplices(0) do
  clr:=clr+1;
  leaves:=[Position(Colours,0)];
  while Size(leaves)>0 do
    for v in leaves do
      Colours[v]:=clr;
    od;
    newleaves:=[];
    for v in leaves do
      Append(newleaves,
      Filtered(EDGES[v],w->Colours[w]=0));
    od;
    leaves:=Difference(Flat(newleaves), leaves);
  od;

od;

nr:=Length(SSortedList(Colours));

if N=0 then
return nr;
fi;

if N=1 and nr=1 then return IntegerSimplicialComplex(K); fi;


L:=StructuralCopy(K!.simplicesLst);
for i in [1..Length(L)] do
L[i]:=Filtered(L[i],x->not false in List(x,a->Colours[a]=N));
od;

return SimplicesToSimplicialComplex(L);;
end);
############################################
############################################

#################################################################
#################################################################
InstallGlobalFunction(PathComponentsOfSimplicialComplex_alt,
function(K,n)
local
        G,A,P,V,L,n1,i;

G:=GraphOfSimplicialComplex(K);
P:=PathComponentsOfGraph(G,0);
A:=G!.pathComponents;

if n=0 then return P; fi;

n1:=n+1;
L:=StructuralCopy(K!.simplicesLst);
V:=Filtered([1..Length(A)],x-> n1 in A[x]);
V:=StructuralCopy(K!.vertices{V});

for i in [1..Length(L)] do
L[i]:=Filtered(L[i],x->IsSubset(V,x));
od;
return SimplicesToSimplicialComplex(L);
end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(ChainComplexOfSimplicialPair,
function(K,A)
local
	Dimen,Boundary,Faces,nullvecs,dim,d,i;

dim:=Dimension(K);
Faces:=[];

for d in [1..dim+1] do
if A!.nrSimplices(d-1)>0 then
Faces[d]:=Difference(K!.simplicesLst[d],A!.simplicesLst[d]);
else
Faces[d]:=K!.simplicesLst[d];
fi;
od;

#############################
Dimen:=function(k);
if k<0 or k>dim then return 0; fi;
return Length(Faces[k+1]);
end;
#############################

nullvecs:=[];
for i in [0..dim] do
nullvecs[i+1]:=List([1..Maximum(Dimen(i),1)],a->0);
od;

#############################
Boundary:=function(n,i)
local V,Vhat, j, bnd,pos;
V:=Faces[n+1][i];

bnd:=StructuralCopy(nullvecs[n]);

for j in [1..Length(V)] do
Vhat:=StructuralCopy(V);
RemoveSet(Vhat,V[j]);
pos:=Position(Faces[n],Vhat);
if not pos=fail then

if IsOddInt(j) then
bnd[pos]:=1;
else
bnd[pos]:=-1;
fi;

fi;
od;

return bnd;
end;
#############################


return
Objectify(HapChainComplex,
           rec(
           dimension:=Dimen,
           boundary:=Boundary,
           properties:=[
           ["length",dim],
           ["type","chainComplex"],
           ["characteristic",0]]
           ));

end);
#################################################################
#################################################################

#############################################
#############################################
InstallGlobalFunction(IntegerSimplicialComplex,
function(K)
local

      Vertices, L, n;

if EvaluateProperty(K,"integer")=true then
  return K;
fi;

#if not false in List(K!.vertices,x->IsInt(x)) then
if K!.vertices=[1..Length(K!.vertices)] then     #Modified this May 2022
  Add(K!.properties,["integer",true]);
  return K;
fi;

Vertices:=K!.vertices;
L:=StructuralCopy(K!.simplicesLst);

for n in [1..Length(L)] do
L[n]:=List(L[n],x->List(x,v->Position(Vertices,v)));
L[n]:=SSortedList(L[n]);
od;

L:=SimplicesToSimplicialComplex(L);
Add(L!.properties,["integer",true]);
return L;

end);
############################################
############################################

#####################################################
#####################################################
InstallGlobalFunction(RandomSimplicialGraph,
function(N,p)
local num,den,P, K, i,j,k, Edges;

num:=NumeratorRat(p);
den:=DenominatorRat(p);

P:=[1..den];


Edges:=List([1..N],i->[i]);

for i in [1..N] do
for j in [i+1..N] do
if Random(P)<=num then
Add(Edges, [i,j]);
fi;
od;od;

K:=MaximalSimplicesToSimplicialComplex(Edges);

return K;
end);
#####################################################
#####################################################



#####################################################
#####################################################
InstallGlobalFunction(RandomSimplicialTwoComplex,
function(N,p)
local num,den,P, K, i,j,k, Faces;

num:=NumeratorRat(p);
den:=DenominatorRat(p);

P:=[1..den];


Faces:=[];

for i in [1..N] do
for j in [i+1..N] do
Add(Faces,[i,j]); #Here we actually add an edge and not a face!
od;od;


for i in [1..N] do
for j in [i+1..N] do
for k in [j+1..N] do
if Random(P)<=num then
Add(Faces, [i,j,k]);
fi;
od;od;od;

K:=MaximalSimplicesToSimplicialComplex(Faces);

return K;
end);
#####################################################
#####################################################


#####################################################
#####################################################
InstallGlobalFunction(FirstHomologySimplicialTwoComplex,
function(K)
local S, c, L, A, Y, x;

Y:=SimplicialComplexToRegularCWComplex(K);
c:=CriticalCellsOfRegularCWComplex(Y,1);
if Length(Filtered(c,x->x[1]=1))=0 then return []; fi;

Apply(K!.simplicesLst[3],x->SSortedList(x));
Apply(K!.simplicesLst[2],x->SSortedList(x));

S:=ContractibleSubcomplexOfSimplicialComplex(K);
L:=SSortedList(S!.simplicesLst[2]);
Apply(L,x->SSortedList(x));
A:=[];

for x in K!.simplicesLst[3] do
if x{[1,2]} in L and x{[1,3]} in L and x{[2,3]} in L then
Add(A,x);
fi;
od;

#Print("Length(A)= ",Length(A),"\n");
A:=Difference(A,S!.simplicesLst[3]);
#Print("Length(A)= ",Length(A),"\n");

L:=Difference(K!.simplicesLst[3],A);
L:=Concatenation(L,K!.simplicesLst[2],K!.simplicesLst[1]);
L:=MaximalSimplicesToSimplicialComplex(L);
Y:=SimplicialComplexToRegularCWComplex(L);

return Homology(Y,1);
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallGlobalFunction(FundamentalGroupSimplicialTwoComplex,
function(K)
local S, c, L, A, Y, x;

Y:=SimplicialComplexToRegularCWComplex(K);
c:=CriticalCellsOfRegularCWComplex(Y,1);
if Length(Filtered(c,x->x[1]=1))=0 then return
Image(IsomorphismFpGroup(Group(()))); fi;

Apply(K!.simplicesLst[3],x->SSortedList(x));
Apply(K!.simplicesLst[2],x->SSortedList(x));

S:=ContractibleSubcomplexOfSimplicialComplex(K);
L:=SSortedList(S!.simplicesLst[2]);
Apply(L,x->SSortedList(x));
A:=[];

for x in K!.simplicesLst[3] do
if x{[1,2]} in L and x{[1,3]} in L and x{[2,3]} in L then
Add(A,x);
fi;
od;

#Print("Length(A)= ",Length(A),"\n");
A:=Difference(A,S!.simplicesLst[3]);
#Print("Length(A)= ",Length(A),"\n");

L:=Difference(K!.simplicesLst[3],A);
L:=Concatenation(L,K!.simplicesLst[2],K!.simplicesLst[1]);
L:=MaximalSimplicesToSimplicialComplex(L);
Y:=SimplicialComplexToRegularCWComplex(L);

return FundamentalGroup(Y);
end);
#####################################################
#####################################################




