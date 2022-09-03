
#####################################################
InstallGlobalFunction(TransferCochainMap,
function(arg)
local R,H,A,B,rnk,S,CR,CS, T, TT, fn, SS,map, HhomH, trans, transTr,
LT, TietzeReduced, hom, elts, mult;

TietzeReduced:=function(R); return R; end;
#TietzeReduced:=TietzeReducedResolution;

R:=arg[1];
H:=arg[2];
S:=ResolutionFiniteSubgroup(R,H);;
trans:=S!.transversal;
trans:=List(trans,x->x^-1);
transTr:=List(trans,x->TransposedMat(x));
LT:=Length(trans);
if Length(arg)>=3 then B:=arg[3]; rnk:=Length(One(Image(B)));
A:=B;
CR:=HomToIntegralModule(R,A);
CS:=HomToIntegralModule(S,B);
trans:=List(trans,g->ImagesRepresentative(A,g));
transTr:=List(trans,x->TransposedMat(x));
fi;
if Length(arg)=2 then
CR:=HomToIntegers(R);;
CS:=HomToIntegers(S);;
fi;

if Length(arg)=2 then
##################################
fn:=function(v,n)
local w, i, x;
w:=0*[1..R!.dimension(n)];
for i in [1..S!.dimension(n)] do
x:=S!.Int2Pair(i);
w[x[1]]:= w[x[1]]+1*v[i];
od;
return w;
end;
##################################
else
##################################
fn:=function(v,n)
local w, u, i, j, k, x, a;
w:=0*[1..CR!.dimension(n)];
for i in [1..S!.dimension(n)] do
   x:=S!.Int2Pair(i);
   #a:=trans[x[2]];
   #a:=TransposedMat(a);
   a:=transTr[x[2]];
   for j in [1..rnk] do
      u:=[1..rnk]*0; u[j]:=1;
      u:=u*a;
      for k in [1..rnk] do
         w[(x[1]-1)*rnk+k]:= w[(x[1]-1)*rnk+k]+v[(i-1)*rnk+j]*u[k];
      od;
od;
od;
return w;
end;
##################################
fi;

T:=Objectify(HapCochainMap,
        rec(
           source:=CS,
           target:=CR,
           mapping:=fn,
           properties:=[ ["type","cochainMap"],
           ["characteristic", Maximum(
           EvaluateProperty(CR,"characteristic"),
           EvaluateProperty(CS,"characteristic"))],
           ]
           ));

T!.resH:=S;
return T;


elts:=S!.elts;               
mult:=elts!.mult;
if IsPseudoList(elts) then   ##A hack to allow TietzeReduced to work!!
if IsBound(elts!.mult) then  ##
Unbind(elts!.mult);          ##
fi;
fi;
SS:=TietzeReduced(S);;
elts!.mult:=mult;

HhomH:=GroupHomomorphismByFunction(H,H,x->x);

###################################
if Length(arg)=4 then
map:=EquivariantChainMap(SS,S,HhomH);;  #The "wrong direction" mapping
map:=HomToIntegralModule(map,B);;
return [T, map, SS];
fi;
###################################

map:=EquivariantChainMap(S,SS,HhomH);;
if Length(arg)>=3 then
map:=HomToIntegralModule(map,B);  #This mapping is slow
else
map:=HomToIntegers(map); 
fi;
TT:=Compose(T,map);
TT!.resH:=SS;
return TT;
end);
#####################################################


#####################################################
InstallGlobalFunction(TransferChainMap,
function(R,H)
local S,CR,CS, T, TT, fn, SS,map, HhomH, TietzeReduced;

TietzeReduced:=function(R); return R; end;
TietzeReduced:=TietzeReducedResolution;
S:=ResolutionFiniteSubgroup(R,H);;
CR:=TensorWithIntegers(R);;
CS:=TensorWithIntegers(S);;

##################################
fn:=function(v,n)
local w, i, x;

w:=0*[1..CS!.dimension(n)];
for i in [1..CS!.dimension(n)] do
x:=S!.Int2Pair(i);
w[i]:= 1*v[x[1]];
od;
return w;
end;
##################################

T:=Objectify(HapChainMap,
        rec(
           source:=CR,
           target:=CS,
           mapping:=fn,
           properties:=[ ["type","chainMap"],
           ["characteristic", Maximum(
           EvaluateProperty(CR,"characteristic"),
           EvaluateProperty(CS,"characteristic"))],
           ]
           ));
SS:=TietzeReduced(S);;
HhomH:=GroupHomomorphismByFunction(H,H,x->x);
map:=EquivariantChainMap(S,SS,HhomH);;
map:=TensorWithIntegers(map);
TT:=Compose(map,T);
TT!.resH:=SS;
return TT;
end);
#####################################################

#####################################################
InstallOtherMethod(Compose,
"Composition of chain maps",
[IsHapChainMap,IsHapChainMap],
function(F,G)
local C,D,FG;
#if not Target(G)=Source(F) then return fail; fi;
C:=Source(G);
D:=Target(F);

##################################
FG:=function(v,n);
return F!.mapping(G!.mapping(v,n),n);
end;
##################################

return  Objectify(HapChainMap,
        rec(
           source:=C,
           target:=D,
           mapping:=FG,
           properties:=[ ["type","chainMap"],
           ["characteristic", Maximum(
           EvaluateProperty(C,"characteristic"),
           EvaluateProperty(D,"characteristic"))]
           ]));

end);
#####################################################

#####################################################
InstallOtherMethod(Compose,
"Composition of cochain maps",
[IsHapCochainMap,IsHapCochainMap],
function(F,G)
local FG,C, D ;
#if not Target(G)=Source(F) then return fail; fi;  #Need to implement equality here
C:=Source(G);
D:=Target(F);
##################################
FG:=function(v,n);
return F!.mapping(G!.mapping(v,n),n);
#a:=G!.mapping(v,n);  This takes time
end;
##################################

return  Objectify(HapCochainMap,
        rec(
           source:=C,
           target:=D,
           mapping:=FG,
           properties:=[ ["type","cochainMap"],
           ["characteristic", Maximum(
           EvaluateProperty(C,"characteristic"),
           EvaluateProperty(D,"characteristic"))]
           ]));

end);
#####################################################


