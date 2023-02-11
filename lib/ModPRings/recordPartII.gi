#C 2008 Graham Ellis

#####################################################################
#####################################################################
InstallGlobalFunction(ModPCohomologyRing_part_2,
function(A)
local

	PowGens,
	NiceBas,
	B,dim,one,zero,degree,bas,nicbas,SCT,gens,gens1,gens2,i,j,Compose;

######################################
PowGens:=function(A,gens,n)
local tmpbas,V,dim,dim2,powgens,lastgens,lastgens1,x,y,z;

if n=1 then return List(gens,x->[x,[x]]);fi;

dim:=0;
lastgens:=PowGens(A,gens,n-1);

if Length(lastgens)=Dimension(A) then
return lastgens; fi;

powgens:=[];
V:=SubspaceNC(A,List(powgens,x->x[1]));
lastgens1:=SSortedList(List(lastgens,x->x[1]));

for x in gens do
for y in lastgens do
z:=x*y[1]; 
if not (IsZero(z) or z in lastgens) then

  if not z in V then
  Add(powgens,[z,Concatenation([x],y[2])]); 
  tmpbas:=List(powgens,w->w[1]);
  V:=SubspaceNC(A,tmpbas);
  V!.Basis:=BasisNC(V,tmpbas);
  if not Length(Basis(V))>dim then return powgens;else dim:=Dimension(V);fi;
  fi;

fi;
od;od;

return powgens;
end;
######################################


#######################################################
NiceBas:=function(A)
local n,basA,basAcoeffs,nicbas,nicbas1Mat,nicbas1, gens,Coeff;;

n:=20;
nicbas:=[];
gens:=ModPRingGenerators(A);

while Length(nicbas)<Dimension(A) do
n:=n+1;
nicbas:=PowGens(A,gens,n);
od;

nicbas1:=BasisNC(A,List(nicbas,x->x[1]));
basA:=Basis(A);
nicbas1Mat:=MutableCopyMat(List(nicbas1,v->Coefficients(basA,v)));
basAcoeffs:=MutableCopyMat(List(basA,x->Coefficients(basA,x)));
#Coeff:=List(basA,v->Coefficients(nicbas1,v));
Coeff:=SolutionsMatDestructive(nicbas1Mat,basAcoeffs);
Apply(Coeff,v->List(v,x->IntFFE(x)));

Apply(nicbas,x->x[2]);
Apply(nicbas,x->List(x,i->Position(Basis(A),i)));

return [Coeff,nicbas];
end;
########################################################

dim:=Dimension(A);
zero:=Zero(A);
#one:=One(A);
bas:=Basis(A);
SetOne(A,bas[1]);  #I think this is OK
nicbas:=NiceBas(A);
#SCT:=StructuralCopy(StructureConstantsTable(bas));
SCT:=EmptySCTable(dim,Zero(LeftActingDomain(A)));
#gens:=ModPRingGenerators(A);
#gens1:=Filtered(gens,a->A!.degree(a)<2);
#gens1:=List(gens1,a->Position(bas,a));
#gens2:=Filtered(gens,a->A!.degree(a)>1);
#gens2:=List(gens2,a->Position(bas,a));
degree:=A!.degree;

if "gensSupport" in NamesOfComponents(A) then
##########################################
Compose:=function(i,j)
local a,b,c,d,x,y,z,coeff,sc,ans;;
a:=[];

for x in [1..dim] do
if not IsZero(nicbas[1][i][x]) then 
for y in [1..nicbas[1][i][x]] do
Add(a,nicbas[2][x]);
od;
fi;
od;

b:=[];
for x in [1..dim] do
if not IsZero(nicbas[1][j][x]) then
for y in [1..nicbas[1][j][x]] do
Add(b,nicbas[2][x]);
od;
fi;
od;

ans:=zero;
for x in a do
for y in b do
c:=Concatenation(
Filtered(x,t->t in A!.gensSupport),
Filtered(y,t->t in A!.gensSupport),
Filtered(x,t->not t in A!.gensSupport),
Filtered(y,t->not t in A!.gensSupport)
);
d:=one;
for z in c do
d:=bas[z]*d;
od;
ans:=ans+d;
od;
od;

sc:=[];
coeff:=Coefficients(bas,ans);
for x in [1..Length(coeff)] do
if not IsZero(coeff[x]) then 
Append(sc,[coeff[x],x]); fi;
od;

return sc;
end;
##########################################
else
##########################################
Compose:=function(i,j)
local a,b,c,x,y,coeff,sc;
a:=[];

for x in [1..dim] do
if not IsZero(nicbas[1][i][x]) then
for y in [1..nicbas[1][i][x]] do
Add(a,nicbas[2][x]);
od;
fi;
od;

b:=zero;
for x in a do
c:=StructuralCopy(bas[j]);
for y in Reversed(x) do
c:=bas[y]*c;
od;
b:=b+c;
od;

sc:=[];
coeff:=Coefficients(bas,b);
for x in [1..Length(coeff)] do
if not IsZero(coeff[x]) then
Append(sc,[coeff[x],x]); fi;
od;

return sc;
end;
##########################################
fi;


for i in [1..Length(bas)] do
for j in [i..Length(bas)] do
#if not ((i in gens1) or (j in gens1) or (i in gens2 and j in gens2)) then
SetEntrySCTable(SCT,i,j,Compose(i,j));
SetEntrySCTable(SCT,j,i,Compose(j,i));
#fi;
od;od;

B:=AlgebraByStructureConstants(LeftActingDomain(A),SCT);


#####################################################################
degree:=function(x)
local i;

i:=Position(GeneratorsOfAlgebra(B),x);

if i=1 then return 0; fi;

return A!.intToPair(i-1)[1];
end;
#####################################################################


B!.degree:=degree;
B!.intToPair:=A!.intToPair;
B!.intToPairModified:=A!.intToPairModified;
B!.pairToInt:=A!.pairToInt;
B!.pairToIntModified:=A!.pairToIntModified;
B!.niceBasis:=nicbas;

return B;

end);
###########################################################
###########################################################
