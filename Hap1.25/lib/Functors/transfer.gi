

#####################################################
InstallGlobalFunction(TransferChainMap,
function(R,H)
local S,CR,CS, T, TT, fn, SS,map, HhomH;

S:=ResolutionFiniteSubgroup(R,H);;
CR:=TensorWithIntegers(R);;
CS:=TensorWithIntegers(S);;

##################################
fn:=function(v,n)
local w, i, x;

w:=0*[1..CS!.dimension(n)];
for i in [1..CS!.dimension(n)] do
x:=S!.Int2Pair(i);
#if x[2]=1 then
w[i]:= 1*v[x[1]];
#fi;
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

SS:=TietzeReducedResolution(S);;
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
if not Target(G)=Source(F) then return fail; fi;
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
local C,D,FG;
#if not Target(G)=Source(F) then return fail; fi;  #Need to implement equality here
C:=Source(G);
D:=Target(F);

##################################
FG:=function(v,n);
return F!.mapping(G!.mapping(v,n),n);
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


