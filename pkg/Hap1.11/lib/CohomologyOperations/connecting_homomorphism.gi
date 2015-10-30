
##################################################
##################################################
InstallGlobalFunction(ConnectingCohomologyHomomorphism,
#function(phi, n)
function(arg)
local phi, n, hom, R, G, A, B, K, CA, CB, CK, CAn, CAn1, CBn, CKn1,
      HCB,HCK, hcb, hck,delta, genshcb, imgenshcb,
      natK, natB, x, xtilde, FBA1, FBA, FAK, FAK1, CAnboundary; 

# n is a nonnegative integer
# phi is a surjective homomomorphism phi:A-->B of G-modules represented as G-outer groups
# K = kernel(phi)
# We will return the homomorphism hom:H^n(G,B)-->H^{n+1}(G,K) represented as G-outer groups

phi:=arg[1];
n:=arg[2];

A:=phi!.Source;
B:=phi!.Target;
K:=KernelOfGOuterGroupHomomorphism(phi);
G:=A!.ActingGroup;
if Length(arg)=3 then R:=arg[3]; else
R:=ResolutionFiniteGroup(G,n+2); fi;
CB:=HomToGModule(R,B);
CBn:=Source(CB!.boundary(n)); CBn:=CBn!.ActedGroup;
CA:=HomToGModule(R,A);
CAn:=Source(CA!.boundary(n)); CAn:=CAn!.ActedGroup;
CAn1:=Source(CA!.boundary(n+1)); CAn1:=CAn1!.ActedGroup;
CAnboundary:=CA!.boundary(n); CAnboundary:=CAnboundary!.Mapping;
CK:=HomToGModule(R,K);
CKn1:=Source(CK!.boundary(n+1)); CKn1:=CKn1!.ActedGroup;
HCB:=CohomologyModule(CB,n);
HCK:=CohomologyModule(CK,n+1);
hcb:=HCB!.ActedGroup;
hck:=HCK!.ActedGroup;

genshcb:=GeneratorsOfGroup(hcb);
#imgenshcb:=List(genshcb,x->Identity(hck));  #WRONG
imgenshcb:=[];

######################## create imgenshcb  ###########

natB:=HCB!.nat; natB:=natB!.Mapping;
natK:=HCK!.nat; natK:=natK!.Mapping;

################
FBA1:=function(i,x)
local y;

y:=Image(Projection(CBn,i),x);
y:=PreImagesRepresentative(phi!.Mapping,y);
y:=Image(Embedding(CAn,i),y);

return y;
end;

################

################
FBA:=function(x)
local i , y;

y:=One(CAn);
for i in [1..R!.dimension(n)] do
y:=y*FBA1(i,x);
od;

return y;
end;
################

################
FAK1:=function(i,x)
local y;

y:=Image(Projection(CAn1,i),x);
y:=Image(Embedding(CKn1,i),y);

return y;
end;

################

################
FAK:=function(x)
local i , y;

y:=One(CKn1);
for i in [1..R!.dimension(n+1)] do
y:=y*FAK1(i,x);
od;

return y;
end;
################



for x in genshcb do
xtilde:=PreImagesRepresentative(natB,x);
xtilde:=FBA(xtilde);
xtilde:=Image( CAnboundary, xtilde);
xtilde:=FAK(xtilde);
xtilde:=Image(natK,xtilde);
Add(imgenshcb,xtilde);
od;
######################## imgenshcb created ###########

delta:=GroupHomomorphismByImages(hcb,hck,genshcb,imgenshcb);  

 
hom:=GOuterGroupHomomorphism();;
hom!.Source:=HCB;
hom!.Target:=HCK;
hom!.Mapping:=delta;  

return hom;
end);
##################################################
##################################################




