#######################################
#######################################
InstallGlobalFunction(SimplicialNerveOfTwoComplex,
function(K,dim)
local  A, Vertices, NrSimplices, Simplices, SimplicesLst, EnumeratedSimplex,
       bool, s, VL,x, y, d, i,j,k,l,m,n,RedundantFaces;

##########################
if not Dimension(K)<=2 then
Print("This function can only be applied to 2-dimensional simplicial complexes.\n");
return fail;
fi;
##########################


Vertices:=StructuralCopy(K!.vertices);
VL:=Length(Vertices);

SimplicesLst:=StructuralCopy(K!.simplicesLst);
Apply(SimplicesLst[3], x->SSortedList(x));
SimplicesLst[3]:=SSortedList(SimplicesLst[3]);


if dim>=3 then
SimplicesLst[4]:=[];
for s in SimplicesLst[3] do
i:=s[1];j:=s[2];k:=s[3];
for l in [k+1..VL] do
if 
SSortedList([i,j,l]) in SimplicesLst[3] and
SSortedList([i,l,k]) in SimplicesLst[3] and
SSortedList([l,j,k]) in SimplicesLst[3] 
then Add(SimplicesLst[4],SSortedList([i,j,k,l])); fi;
od;od;
SimplicesLst[4]:=SSortedList(SimplicesLst[4]);
if Length(SimplicesLst[4])=0 then dim:=2; fi;
fi;

if dim>=4 then
SimplicesLst[5]:=[];
for s in SimplicesLst[4] do
i:=s[1];j:=s[2];k:=s[3];l:=s[4];
for m in [l+1..VL] do
if
SSortedList([m,i,j]) in SimplicesLst[3] and
SSortedList([m,i,k]) in SimplicesLst[3] and
SSortedList([m,i,l]) in SimplicesLst[3] and
SSortedList([m,j,k]) in SimplicesLst[3] and
SSortedList([m,j,l]) in SimplicesLst[3] and
SSortedList([m,k,l]) in SimplicesLst[3] 
then Add(SimplicesLst[5],SSortedList([i,j,k,l,m])); fi;
od;od;
SimplicesLst[5]:=SSortedList(SimplicesLst[5]);
if Length(SimplicesLst[5])=0 then dim:=3; fi;
fi;

###########################GOT THIS FAR
if dim>=5 then
SimplicesLst[6]:=[];
for s in SimplicesLst[5] do
i:=s[1];j:=s[2];k:=s[3];l:=s[4];m:=s[5];
for n in [m+1..VL] do

if
SSortedList([n,i,j]) in SimplicesLst[3] and
SSortedList([n,i,k]) in SimplicesLst[3] and
SSortedList([n,i,l]) in SimplicesLst[3] and
SSortedList([n,i,m]) in SimplicesLst[3] and
SSortedList([n,j,k]) in SimplicesLst[3] and
SSortedList([n,j,l]) in SimplicesLst[3] and
SSortedList([n,j,m]) in SimplicesLst[3] and
SSortedList([n,k,l]) in SimplicesLst[3] and
SSortedList([n,k,m]) in SimplicesLst[3] and
SSortedList([n,l,m]) in SimplicesLst[3] 
then Add(SimplicesLst[6],SSortedList([i,j,k,l,m,n])); fi;
od;od;
SimplicesLst[6]:=SSortedList(SimplicesLst[6]);
if Length(SimplicesLst[6])=0 then dim:=4; fi;
fi;

if dim>=6 then Print("The function is only implemented up to dimension 5.\n");
return fail;
fi;

for d in [6..dim] do
SimplicesLst[d+1]:=[];
for y in Combinations(Vertices,d+1) do
bool:=true;
for x in Combinations(y,2) do
if A[x[1]][x[2]]=0 then bool:=false; break; fi;
od;
if bool then Add(SimplicesLst[d],y); fi;
od;
SimplicesLst[d+1]:=SSortedList(SimplicesLst[d+1]);
if Length(SimplicesLst[d+1])=0 then dim:=d-1; fi;
od;
##########################

##########################
Simplices:=function(d,k);
return SimplicesLst[d+1][k];
end;
#########################


##########################
NrSimplices:=function(d);
return Length(SimplicesLst[d+1]);
end;
#########################

#########################
EnumeratedSimplex:=function(v);
return PositionSet(SimplicesLst[Length(v)],v);
end;
#########################

Add(SimplicesLst,[]);

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
#####################################################################


