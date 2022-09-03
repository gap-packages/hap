
###############################################################
###############################################################
InstallGlobalFunction(EilenbergMacLaneSimplicialFreeAbelianGroup,
function(G,n,m)
local K, L, VectorToWord, M, D, i, j,  dimension, boundary, degeneracy,
MFUN, DFUN;

K:=EilenbergMacLaneSimplicialGroup(G,n,m);;

L:=List([0..m],i->K!.groupsList(i));
L:=List(L,AbelianInvariants);
L:=List(L,Length);

###################################
VectorToWord:=function(v)
local w,i;
w:=[];
for i in [1..Length(v)] do
if not v[i]=0 then Add(w,[i,v[i]]); fi;
od;
return w;
end;
###################################

M:=[];;
for i in [1..m] do
MFUN:=function(i);
M[i]:=[];
for j in [0..i] do
M[i][j+1]:=K!.boundariesList(i,j);
M[i][j+1]:=M[i][j+1]!.MappingGeneratorsImages[2];
M[i][j+1]:=StructuralCopy(M[i][j+1]);
M[i][j+1]:=List(M[i][j+1],Exponents);
M[i][j+1]:=List(M[i][j+1],v->VectorToWord(v));
od;
end;
#MFUN(i);
od;

D:=[];;
for i in [1..m] do  #so we are in degree i-1
DFUN:=function(i);
D[i]:=[];
if L[i]>0 then
for j in [0..i-1] do 
D[i][j+1]:=K!.degeneraciesList(i-1,j);
D[i][j+1]:=D[i][j+1]!.MappingGeneratorsImages[2];
D[i][j+1]:=StructuralCopy(D[i][j+1]);
D[i][j+1]:=List(D[i][j+1],Exponents);
D[i][j+1]:=List(D[i][j+1],v->VectorToWord(v));
od;
fi; 
end;
#DFUN(i);
od;


#####################################
dimension:=function(n);
if n<0 or n>Length(L)-1 then return 0; fi;
return L[n+1];
end;
#####################################

#####################################
boundary:=function(n,j,k); #degree n, boundary map j, generator k
if n<1 or n>Length(L)-1 then return []; fi;
if not IsBound(M[n]) then MFUN(n); fi;
return 1*M[n][j+1][k];
end;
#####################################

#####################################
degeneracy:=function(n,j,k); #degree n, degeneracy map j, generator k
if n<0 or n>Length(L)-2 then return []; fi;
if not IsBound(D[n+1]) then DFUN(n+1); fi;
return 1*D[n+1][j+1][k];
end;
#####################################


return 
Objectify(HapSimplicialFreeAbelianGroup,
                rec(
                dimension:=dimension,
                boundary:=boundary,
                degeneracy:=degeneracy,
                properties:=
                [["length",m],
                ["moorecomplexlength",n],
                ["type", "simplicialFreeAbelianGroup"],
                ["characteristic", 0 ] ]));

end);
###############################################################
###############################################################

###############################################################
###############################################################
InstallGlobalFunction(E1HomologyPage,
function(K,q,l)
local dimension, boundary, facemap, negativefacemap,properties, bases, dims, len, BoundRec, degs, x, y, i, n, newbases, newdims;

len:=Minimum(EvaluateProperty(K,"length") , l);
bases:=[];
dims:=[];
BoundRec:=[];
for n in [0..len] do
dims[n+1]:=NrCombinations([1..K!.dimension(n)],q);
bases[n+1]:=Combinations([1..K!.dimension(n)],q) ;
BoundRec[n+1]:=[];
od;

#################################################
############### WARNING:  THIS IS ONLY VALID FOR E-M COMPLEXES!!!!
newbases:=[];
newdims:=[];
newbases[1]:=bases[1];
newdims[1]:=dims[1];
for n in [1..len] do
degs:=[];
for x in bases[n] do
for i in [0..n-1] do
y:=List(x,k->K!.degeneracy(n-1,i,k)[1][1]);
y:=SSortedList(y);
Add(degs,y);
od;
od;
degs:=SSortedList(degs);
newbases[n+1]:=Difference(bases[n+1],degs);
newbases[n+1]:=SSortedList(newbases[n+1]);
newdims[n+1]:=Length(newbases[n+1]);

od;

bases:=newbases;
dims:=newdims;
################ WARNING: FINISHED
###########################################

###########################
dimension:=function(n)
local p;
if n>len or n<0 then return 0; fi;
return dims[n+1];
end;
###########################
facemap:=function(n,i,k)
local c,d,F,newF,x,y,yy,j,m,s, z,ans,pos;

c:=bases[n+1][k];
d:=[];     
for x in c do
y:=K!.boundary(n,i,x);
yy:=[];
for z in y do
s:=SignInt(z[2]);
for j in [1..AbsInt(z[2])] do
Add(yy,[z[1],s]);
od;
od;
Unbind(y);
Add(d,yy);
od;


d:=Cartesian(d);

ans:=[];

for y in d do
s:=Product(List(y,a->a[2]));
y:=List(y,a->a[1]);
yy:=y;
y:=SortedList(y);  
#RT:=RT-Runtime();  #THIS TAKES A CRAZY AMOUNT OF TIME USING POSITION()
pos:=PositionSorted(bases[n],y);
if IsBound(bases[n][pos]) then 
     if not bases[n][pos]=y then pos:=fail; fi;
else pos:=fail;
fi; 
#RT:=RT+Runtime();
if not pos=fail then 
yy:=List(yy,i->PositionSorted(y,i));
s:=s*SignPerm(PermList(yy)); 
Add(ans, [pos,s]);
fi;
od;

return ans;
end;
###########################
negativefacemap:=function(n,i,k)
local ans;
ans:=facemap(n,i,k);
ans:=List(ans,x->[x[1],-x[2]]);
return ans;
end;
###########################
boundary:=function(n,k)
local b, i;
if not IsBound(BoundRec[n+1][k]) then 
b:=[];
for i in [0..n] do
if IsEvenInt(i) then Append(b,facemap(n,i,k));
else Append(b,negativefacemap(n,i,k)); fi;
od;
BoundRec[n+1][k]:=b;
fi;
return BoundRec[n+1][k];
end;
###########################

return
Objectify(HapSparseChainComplex,
                rec(
                dimension:=dimension,
                boundary:=boundary,
bases:=bases,
                properties:=
                [["length",EvaluateProperty(K,"length")],
                ["connected",true],
                ["type", "chainComplex"],
                ["characteristic", 0] ]));

end);
###############################################################
###############################################################


###############################################################
###############################################################
InstallGlobalFunction(HomologySimplicialFreeAbelianGroup,
function(K,n)
local p, q, H, D, ans, mlen, bool;

bool:=false;

if n<0 then return []; fi;
if n=0 then return [0]; fi;

mlen:=EvaluateProperty(K,"moorecomplexlength");
ans:=[];
for  p in [0..n-1] do
   q:=n-p;
   if p<=q*mlen then
if p+1>EvaluateProperty(K,"length") then bool:=true; fi;
      D:=ContractedComplex(E1HomologyPage(K,q,p+1));
      H:= Homology(D,p) ;
      if Length(H)>0 then 
         Append(ans,H);
      fi;
   fi;
od;

#if n>=EvaluateProperty(K,"length") then
if bool then 
Print("WARNING: The computed answer is guaranteed only to be a summand of the correct answer. You'll need to use more terms of the simplicial abelian group to make sure there is no other summand.\n");
fi;

return ans;
end);
###############################################################
###############################################################

#####################################################################
InstallOtherMethod(Homology,
"integral homology of simplicial free abelian group",
[IsHapSimplicialFreeAbelianGroup,IsInt],
function(K,n) return HomologySimplicialFreeAbelianGroup(K,n);
end);
#####################################################################
