#############################################
#############################################
HomomorphismOfDirectProduct:=function(phi,n)
local A, B, L1, L2, C, D, F, gensC, imgensC, i, x, y;

A:=Source(phi);
B:=Target(phi);
L1:=List([1..n],i->A);
L2:=List([1..n],i->B);
C:=DirectProduct(L1);
D:=DirectProduct(L2);

################
F:=function(i,x)
local xi;

xi:=Image(Projection(C,i),x);
xi:=Image(phi,xi);
xi:=Image(Embedding(D,i),xi);

return xi;
end;
################

gensC:=GeneratorsOfGroup(C);
imgensC:=[];

for x in gensC do
y:=Identity(D);
        for i in [1..n] do
        y:=y*F(i,x);
        od;
Add(imgensC,y);
od;

return GroupHomomorphismByImages(C,D,gensC,imgensC);
end;
#############################################
#############################################


##################################################
##################################################
InstallGlobalFunction(CohomologyHomomorphism,
#function(phi, n)
function(arg)
local phi,n,
      hom, R, G, A, B, C,D,HC,HD, hc, hd,nat, natd, p, fdelta, hcdelta, N, indhom,
      genshc, imgenshc, x, y, xtilde, ytilde, homdir, iso, isod;

phi:=arg[1];
n:=arg[2];

# n is a nonnegative integer
# phi is a homomomorphism phi:A-->B of G-modules represented as G-outer groups
# We will return the homomorphism hom:H^n(G,A)-->H^n(G,B) represented as G-outer groups


A:=phi!.Source;
B:=phi!.Target;
G:=A!.ActingGroup;
if Length(arg)=3 then R:=arg[3]; else
R:=ResolutionFiniteGroup(G,n+1);
fi;
C:=HomToGModule(R,A);
D:=HomToGModule(R,B);
HC:=CohomologyModule(C,n);
HD:=CohomologyModule(D,n);
hc:=HC!.ActedGroup;
hd:=HD!.ActedGroup;
nat:=HC!.nat;
nat:=nat!.Mapping;
natd:=HD!.nat;
natd:=natd!.Mapping;


p:=phi!.Mapping;
#Source(nat);
#hcdelta:=Source(nat);
# hcdelta!.ParentAttr;
#N:=hcdelta!.ParentAttr;


# indhom:hc-->hd is the induced homomorphism of groups
genshc:=GeneratorsOfGroup(hc);
#imgenshc:=List(genshc,x->Identity(hd));
imgenshc:=[];

homdir:=HomomorphismOfDirectProduct(p,R!.dimension(n));

iso:=GroupHomomorphismByImages(Source(nat)!.ParentAttr,Source(homdir),
GeneratorsOfGroup(Source(nat)!.ParentAttr),
GeneratorsOfGroup(Source(homdir))  );

isod:=GroupHomomorphismByImages(Target(homdir),Source(natd)!.ParentAttr,
GeneratorsOfGroup(Target(homdir)),
GeneratorsOfGroup(Source(natd)!.ParentAttr)  );
##################
for x in genshc do

xtilde:=PreImagesRepresentative(nat,x);
xtilde:=Image(iso,xtilde);
ytilde:=Image(homdir,xtilde);
ytilde:=Image(isod,ytilde);
y:=Image(natd,ytilde);

Add(imgenshc,y);
od;
##################

indhom:=GroupHomomorphismByImages(hc,hd,genshc,imgenshc);  


 
hom:=GOuterGroupHomomorphism();;
hom!.Source:=HC;
hom!.Target:=HD;
hom!.Mapping:=indhom;  

return hom;
end);
##################################################
##################################################



