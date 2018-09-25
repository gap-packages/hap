
############################################################
############################################################
InstallGlobalFunction(AbelianGOuterGroupToCatOneGroup,
function(N)
    local
    A,G,GC,alpha,
    AutA,GensA,GensG,GensGC,phi,
    Pro,Emb1,sC,tC;
   
    A:=ActedGroup(N);
    G:=ActingGroup(N);
    alpha:=OuterAction(N);
    AutA:=AutomorphismGroup(A);
    GensA:=GeneratorsOfGroup(A);
    GensG:=GeneratorsOfGroup(G);
    phi:=GroupHomomorphismByImages(G,AutA,GensG,List(GensG,g->GroupHomomorphismByImages(A,A,GensA,List(GensA,a->alpha(g,a)))));
    GC:=SemidirectProduct(G,phi,A);
    Pro:=Projection(GC);
    Emb1:=Embedding(GC,1);
    GensGC:=GeneratorsOfGroup(GC);
    sC:=GroupHomomorphismByImages(GC,GC,GensGC,List(GensGC,x->Image(Emb1,Image(Pro,x))));
    tC:=GroupHomomorphismByImages(GC,GC,GensGC,List(GensGC,x->Image(Emb1,Image(Pro,x))));
    return Objectify(HapCatOneGroup,
                    rec(sourceMap:=sC,
                    targetMap:=tC));
end);
############################################################
############################################################

############################################################
############################################################
InstallMethod(MooreComplex,
"Moore complex of cat-1-groups",
[IsHapCatOneGroup],

function(C)
local M,N,d,fn,ans, gensM;

if not "mooreComplex" in NamesOfComponents(C) then 
M:=Kernel(C!.sourceMap);
N:=Image(C!.sourceMap);
gensM:=GeneratorsOfGroup(M);
if Length(gensM)=0 then gensM:=[One(M)]; fi;
d:=GroupHomomorphismByImagesNC(M,N,gensM,
List(gensM,x->Image(C!.targetMap,x)));

#########
fn:=function(i);
if i=1 then return GOuterGroup(d); fi;
if i=2 then return
GOuterGroup(
GroupHomomorphismByFunction(Group(Identity(Source(d))),Source(d), x->x)
);
fi;
end;
#########
ans:=Objectify(HapGComplex, 
    rec(
    boundary:=fn,
    properties:=[["length",2]] 
     )
      );

C!.mooreComplex:=ans;
fi;

return C!.mooreComplex;
end);
############################################################
############################################################

############################################################
############################################################
InstallMethod(HomotopyModule,
"Homotopy groups of cat-1-groups",
[IsHapCatOneGroup,IsInt],

function(C,n)
local hm, pi,G,alpha,bnd,nat,PI;

if n=2 then
#bnd:=MooreComplex(C)[1];
bnd:=MooreComplex(C)!.boundary;
bnd:=bnd(1);
bnd:=bnd!.Mapping;

nat:=NaturalHomomorphismByNormalSubgroup(Range(bnd),Image(bnd));
G:=Image(nat);
pi:=Kernel(bnd);

	###############################
	alpha:=function(g,a)
        local x;
	x:=PreImagesRepresentative(nat,g);
        return x*a*x^-1;
	end;
	###############################

PI:=GOuterGroup();
SetActingGroup(PI,G);
SetActedGroup(PI,pi);
SetOuterAction(PI,alpha);
return PI;
fi;


Print("Only defined for n=2.\n");
return fail;;
end);
############################################################
############################################################

############################################################
############################################################
InstallMethod(HomotopyGroup,
"Homotopy groups of cat-1-groups",
[IsHapCatOneGroup,IsInt],

function(C,n)
local hm, nat, bnd;

bnd:=MooreComplex(C)!.boundary;
bnd:=bnd(1);
bnd:=bnd!.Mapping;
if n=2 then
return Kernel(bnd);
fi;

if n=1 then
nat:=NaturalHomomorphismByNormalSubgroup(
Range(bnd),  Image(bnd));
return Image(nat);
fi;

return Group(());;
end);
############################################################
############################################################




############################################################
############################################################
InstallGlobalFunction(HasTrivialPostnikovInvariant,
function(arg)
local
	C,R,M,G,phi,
	gensG,
	LiftOne,
	LiftTwo,
	LiftThree,
	OneCocycle,
	TwoCocycle,
	A,D,i,j,P,
	Prd,x,v,
	gensA, CoBound, CoBoundPlus;

####INPUT ARG########
C:=arg[1];
M:=MooreComplex(C)[1];
phi:=NaturalHomomorphismByNormalSubgroup(Range(M),Image(M));
G:=Image(phi);

if Length(arg)=2 then R:=arg[2];
else R:=ResolutionFiniteGroup(G,3); fi;
####ARG INPUT########

#####IS THE INVARIANT OBVIOUSLY TRIVIAL?#####
if Order(Image(M))=1 then return true; fi;
#############################################

#####FIND GENERATORS OF RANGE(M) CORRESPONDING TO R#######
gensG:=List([1..R!.dimension(1)],i->
R!.elts[
Filtered( Flat(R!.boundary(1,i)), y->not AbsInt(y)=1)[1] ]
);

LiftOne:=List(gensG,x->PreImagesRepresentative(phi,x));
##########################################################

#Here we are using the free crossed resolution 
#instead of the free resolution R.

P:=PresentationOfResolution(R);
OneCocycle:=GroupHomomorphismByImagesNC(P.freeGroup,Range(M),
GeneratorsOfGroup(P.freeGroup),LiftOne);

LiftTwo:=[];
for i in [1..R!.dimension(2)] do
Add(LiftTwo,PreImagesRepresentative(M,Image(OneCocycle,P.relators[i])));
od;


#######################
TwoCocycle:=function(w)
local ans,x,y;		#Will only give correct answer
			#when w is an identity among relations.

ans:=Identity(Source(M));
for x in w do
y:=PreImagesRepresentative(phi,R!.elts[x[2]]);
ans:=ans*y*
LiftTwo[AbsInt(x[1])]^(SignInt(x[1]))*y^-1;
od;

return ans;
end;
#######################

LiftThree:=[];
for i in [1..R!.dimension(3)] do
Add(LiftThree,TwoCocycle(R!.boundary(3,i)));
od;

#############################################################
#Now we'll compute the cokernel of
#Hom(R_2,A) --> Hom(R_3,A) where
#A is the second homotopy group of the cat-1-group C.

A:=HomotopyGroup(C,2);
D:=List([1..R!.dimension(3)],i->A);
D:=DirectProduct(D);
gensA:=MinimalGeneratingSet(A);

CoBound:=[];
for i in [1..R!.dimension(2)] do
for x in gensA do
v:=Identity(D);
for j in [1..R!.dimension(3)] do
Prd:=List(
Filtered(R!.boundary(3,j),s ->AbsInt(s[1])=i) ,
w->
PreImagesRepresentative(phi, R!.elts[w[2]])
*x^SignInt(w[1])*
PreImagesRepresentative(phi,R!.elts[w[2]])^-1
);

if Length(Prd)>0 then
v:=v*Image(Embedding(D,j),Product(Prd));
fi;
od;
Add(CoBound,v);
od;
od;

####################################################
#Now we'll add the Postnikov invariant to the coboundary image.
v:=Identity(D);
for j in [1..R!.dimension(3)] do
v:=v*Image(Embedding(D,j),LiftThree[j]);
od;
CoBoundPlus:=Concatenation(CoBound,[v]);
####################################################


return 
AbelianInvariants(D/Group(CoBound))
=
AbelianInvariants(D/Group(CoBoundPlus));

end);
############################################################
############################################################

############################################################
############################################################
InstallOtherMethod(HomotopyGroup,
"Homotopy group of simplicial group",
[IsHapSimplicialGroup,IsInt],
function(G,n)
local M, d1,d2;
M:=MooreComplex(G);

if n=1 then
d2:=M!.boundary(n);
d2:=d2!.Mapping;
return Range(d2)/Image(d2);
fi;

d1:=M!.boundary(n-1);
d1:=d1!.Mapping;
d2:=M!.boundary(n);
d2:=d2!.Mapping;
return Kernel(d1)/Image(d2);
end);
############################################################
############################################################

