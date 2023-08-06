###############################################
###############################################
InstallGlobalFunction(DirectProductOfSimplicialComplexes,
function(K,L)
local vertices, ProdSimplices, F, FK, FL,x,y;

vertices:=Cartesian(K!.vertices,L!.vertices);
FK:=MaximalSimplicesOfSimplicialComplex(K);
FL:=MaximalSimplicesOfSimplicialComplex(L);
F:=[];

########################################
ProdSimplices:=function(a,b)
local dim,la,lb,ca,cb,ea,eb,x,y,i,j,expa,expb,pair2simp,ans;

la:=Length(a);
lb:=Length(b);
dim:=la+lb-1;
ca:=Combinations([2..dim],la-1);
cb:=Combinations([2..dim],lb-1);
expa:=[];
expb:=[];

for x in ca do
ea:=List([1..x[1]-1],i->a[1]);
for i in [1..la-2] do
for j in [x[i]..x[i+1]-1] do
ea[j]:=a[i+1];
od;
for j in [x[la-1]..dim] do
ea[j]:=a[la];
od;
od;
Add(expa,ea);
od;

for x in cb do
eb:=List([1..x[1]-1],i->b[1]);
for i in [1..lb-2] do
for j in [x[i]..x[i+1]-1] do
eb[j]:=b[i+1];
od;
for j in [x[lb-1]..dim] do
eb[j]:=b[lb];
od;
od;
Add(expb,eb);
od;

#######################
pair2simp:=function(v,w);
return List([1..Length(v)], k-> Position(vertices,[v[k],w[k]]));
end;
#######################

ans:=[];
for x in expa do
for y in expb do
Add(ans,pair2simp(x,y));
od;od;

ans:=Filtered(ans,x->Length(x)=Length(SSortedList(x)));
return ans;
end;
########################################

for x in FK do
for y in FL do
Append(F,ProdSimplices(x,y));
od;od;

return SimplicialComplex(F);;
end);
######################################################
######################################################
#
#


