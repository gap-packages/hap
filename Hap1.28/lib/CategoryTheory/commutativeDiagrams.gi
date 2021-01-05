#(C) 2009 Graham Ellis

##########################################################
##########################################################
InstallGlobalFunction(HomomorphismChainToCommutativeDiagram,
function(L)
local
	Ln,Arrows,Objects,i;

Ln:=Length(L);
Arrows:=List([1..Ln],i->[]);
Objects:=[1..Ln+1];

for i in [1..Ln] do
Arrows[i]:=[L[i],i,i+1];
od;

return Objectify(HapCommutativeDiagram,
           rec(
           objects:=Objects,
	   arrows:=Arrows,
           properties:=[
           ["type","chain"]]
           ));

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(NerveOfCommutativeDiagram,
function(D)
local
	Obs,Arrs,ExitArrs,Tree,Answer,N,i;

###############################
if not IsHapCommutativeDiagram(D) then
Print("The function must be applied to a commutative diagram.\n");
fi;
###############################


Obs:=D!.objects;
Arrs:=D!.arrows;
ExitArrs:=List([1..Size(Obs)],i->[]);
#ExitArrs[i] is a list of integers j such that Arr[j] has source i.

for i in [1..Length(Arrs)] do
Add(ExitArrs[Arrs[i][2]], Arrs[i]);
od;

########################################
Tree:=function(x)
#Tree(x) will be a list of composed arrows corresponding to
#a directed maximal tree in the directed diagram with root x.
local
	T,Leaves,NewLeaves,y,f,yf;
##The tree T will grow, and at any stage Leaves is the latest
# round of elements that are to be added to T.  
 
T:=[];
Leaves:=StructuralCopy(ExitArrs[x]);

while Size(Leaves)>0 do
Append(T,Leaves);
NewLeaves:=[];
for y in Leaves do
for f in ExitArrs[y[3]] do
yf:=[y[1]*f[1],y[2],f[3]];
Add(NewLeaves,yf);
od;
od;
Leaves:=NewLeaves;
od;

return T;
end;
########################################

Answer:=[];
for i in Obs do
Append(Answer,Tree(i));
od;

return Objectify(HapCommutativeDiagram,
           rec(
           objects:=Obs,
           arrows:=Answer,
           properties:=[
           ["type","nerve", EvaluateProperty(D,"type")]]
           ));


end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(GroupHomologyOfCommutativeDiagram,
function(arg)
	local D, deg, prime, Resol, ResLst, Obs, Arrs, tmp, 
	hom,RG,G,RQ,eqmap,map,HomologyArrs, x,y;;

#####INPUT#################
D:=arg[1]; #A commutative diagram of groups
deg:=arg[2]; #The degree in which homology will be calculated.
if Length(arg)>2 then
prime:=arg[3];
else prime:=PrimePGroup(Source(D!.arrows[1][1])); fi;
if Length(arg)>3 then
Resol:=arg[4]; else
  if not prime=0 then
  Resol:=ResolutionPrimePowerGroup;
  else
  Resol:=ResolutionFiniteGroup;
  fi;
fi;
###########################

Obs:=D!.objects;
Arrs:=D!.arrows;

#####################
ResLst:=[];
tmp:=[];

for x in Arrs do
if not x[2] in tmp then
Add(tmp,x[2]);
G:=Source(x[1]);
if Order(G)=1 then
RG:=ResolutionFiniteGroup(G,deg+1);
else
RG:=Resol(G,deg+1);
fi;
ResLst[x[2]]:=RG;
fi;
if not x[3] in tmp then
Add(tmp,x[3]);
G:=Range(x[1]);
if Order(G)=1 then
RG:=ResolutionFiniteGroup(G,deg+1);
else
RG:=Resol(G,deg+1);
fi;
ResLst[x[3]]:=RG;
fi;
od;
############################

HomologyArrs:=[];

for x in Arrs do
RG:=ResLst[x[2]];
RQ:=ResLst[x[3]];
hom:=x[1];
eqmap:=EquivariantChainMap(RG,RQ,hom);
map:=TensorWithIntegersModP(eqmap,prime);
Add(HomologyArrs,[HomologyVectorSpace(map,deg), x[2], x[3]]);
od;

return Objectify(HapCommutativeDiagram,
           rec(
           objects:=StructuralCopy(Obs),
           arrows:=HomologyArrs,
           properties:=[
           ["type",EvaluateProperty(D,"type")]]
           ));

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(PersistentHomologyOfCommutativeDiagramOfPGroups,
function(D,n)
local H;

H:=GroupHomologyOfCommutativeDiagram(D,n);
H:=List(H!.arrows, x->
[Dimension(Source(x[1])),
Dimension(Image(x[1])),
Dimension(Range(x[1]))]);

return H;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(NormalSeriesToQuotientDiagram,
function(arg)
#SS is a decreasing or increasing series of groups, and RR is 
#a subseries of SS;
local
        SS,RR,S,R,H, gens,
        HOMS,N,P,G,Q,iso,HhomHmodS,HhomHmodR,RShoms,
	SHOMS,RHOMS,RSHOMS,hom,hom1,Objects,Arrows,i;

#####################
if Length(arg)=1 then
S:=NormalSeriesToQuotientHomomorphisms(arg[1]);
return HomomorphismChainToCommutativeDiagram(S);
fi;
#####################

SS:=arg[1]; RR:=arg[2];

##### ERROR ##########
if Length(SS)<Length(RR) then
Print("Error: the second series is too long.\n");
return fail;
fi; 
######################

if Size(SS[1])<=Size(SS[Length(SS)])
then S:=SS;
else
S:=Reversed(SS);
fi;
if Size(RR[1])<=Size(RR[Length(RR)])
then R:=RR;
else
R:=Reversed(RR);
fi;

##### ERROR #########
if not S[Length(S)]=R[Length(R)] then 
Print("The two series are subseries of different groups.\n");
return fail;
fi;
#####################  

for i in [Length(R)+1..Length(S)] do
R[i]:=R[Length(R)];
od;
#NWe now  hope that R is a subseries of S, and both have same length.
#######################

H:=S[Length(S)];


HhomHmodS:=[];
for i in [1..Length(S)] do
HhomHmodS[i]:=NaturalHomomorphismByNormalSubgroup(H,S[i]);
od;

HhomHmodR:=[];
for i in [1..Length(S)] do
HhomHmodR[i]:=NaturalHomomorphismByNormalSubgroup(H,R[i]);
od;


########################################
SHOMS:=[HhomHmodS[1]];
for i in [1..Length(S)-1] do
P:=Range(HhomHmodS[i]);
N:=List(GeneratorsOfGroup(S[i+1]),x->Image(HhomHmodS[i],x));
Add(N,Identity(P));
N:=Subgroup(P,N);
hom:=NaturalHomomorphismByNormalSubgroup(P,N);
gens:=GeneratorsOfGroup(Range(hom));
iso:=GroupHomomorphismByImagesNC(Range(hom),Range(HhomHmodS[i+1]),
gens, List(gens, x->
Image(HhomHmodS[i+1], PreImagesRepresentative(HhomHmodS[i],PreImagesRepresentative(hom,x)))   )   );
hom1:=GroupHomomorphismByImagesNC(P,Range(HhomHmodS[i+1]),
GeneratorsOfGroup(P),
List(GeneratorsOfGroup(P), x-> Image(iso,Image(hom,x))));
Add(SHOMS,hom1);

od;
#########################################


########################################
RHOMS:=[HhomHmodR[1]];
for i in [1..Length(R)-1] do
P:=Range(HhomHmodR[i]);
N:=List(GeneratorsOfGroup(R[i+1]),x->Image(HhomHmodR[i],x));
Add(N,Identity(P));
N:=Subgroup(P,N);
hom:=NaturalHomomorphismByNormalSubgroup(P,N);
gens:=GeneratorsOfGroup(Range(hom));
iso:=GroupHomomorphismByImagesNC(Range(hom),Range(HhomHmodR[i+1]),
gens, List(gens, x->
Image(HhomHmodR[i+1], PreImagesRepresentative(HhomHmodR[i],PreImagesRepresentative(hom,x)))   )   );
hom1:=GroupHomomorphismByImagesNC(P,Image(HhomHmodR[i+1]),
GeneratorsOfGroup(P),
List(GeneratorsOfGroup(P), x-> Image(iso,Image(hom,x))));
Add(RHOMS,hom1);

od;
#########################################


########################################
iso:=GroupHomomorphismByImagesNC(H,H, 
GeneratorsOfGroup(H), GeneratorsOfGroup(H) );
RSHOMS:=[iso];
for i in [1..Length(RHOMS)] do
P:=Range(RHOMS[i]);
Q:=Range(SHOMS[i]);
gens:=GeneratorsOfGroup(P);;
if i>1 then
hom:=GroupHomomorphismByImagesNC(P,Q,
gens, List(gens , x->
Image(SHOMS[i], Image(HhomHmodS[i-1],PreImagesRepresentative(HhomHmodR[i-1], PreImagesRepresentative(RHOMS[i],x)))) 
));
else
hom:=GroupHomomorphismByImagesNC(P,Q,
gens, List(gens , x->
Image(SHOMS[i], PreImagesRepresentative(RHOMS[i],x))
));
fi;

RSHOMS[i+1]:=hom;
od;

#########################################


for i in [1..Length(SHOMS)] do
SHOMS[i]:=[SHOMS[i],2*(i-1)+1,2*i+1];
RHOMS[i]:=[RHOMS[i],2*i,2*(i+1)];
od;

for i in [1..Length(SHOMS)+1] do
RSHOMS[i]:=[RSHOMS[i],2*i,2*(i-1)+1];
od;

Arrows:=Concatenation( [SHOMS,RHOMS, RSHOMS] );
Objects:=[1..2*(Length(S)+1)];


return Objectify(HapCommutativeDiagram,
           rec(
           objects:=Objects,
           arrows:=Arrows,
           properties:=[
           ["type","bichain"]]
           ));


end);
##########################################################
##########################################################

