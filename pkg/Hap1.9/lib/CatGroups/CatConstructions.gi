##################################################################
##################################################################
InstallGlobalFunction(AutomorphismGroupAsCatOneGroup,
function(G)
local
	A,C,
	Smap,Tmap,
	CmapG;

A:=AutomorphismGroup(G);
C:=SemidirectProduct(A,G);
Smap:=GroupHomomorphismByImagesNC(C,C,GeneratorsOfGroup(C),
List(GeneratorsOfGroup(C),x->
  Image( Embedding(C,1), Image( Projection(C),x)   )));

############################
CmapG:=function(x);
return
PreImagesRepresentative(Embedding(C,2), Image(Embedding(C,1), Image(Projection(C),x)^-1) *x);
end;
############################

Tmap:=GroupHomomorphismByImagesNC(C,C,GeneratorsOfGroup(C),
List(GeneratorsOfGroup(C),x->
  Image(Embedding(C,1),InnerAutomorphism(G,CmapG(x)))*
  Image( Embedding(C,1), Image( Projection(C),x)   )));



return Objectify(HapCatOneGroup,
rec(sourceMap:=Smap,
targetMap:=Tmap));
end);
##################################################################
##################################################################

##################################################################
##################################################################
InstallGlobalFunction(GModuleAsCatOneGroup,
function(G,alpha,M)
local
	C,Smap;
C:=SemidirectProduct(G,alpha,M);

Smap:=GroupHomomorphismByImagesNC(C,C,GeneratorsOfGroup(C),
List(GeneratorsOfGroup(C),x->
  Image( Embedding(C,1), Image( Projection(C),x)   )));

return Objectify(HapCatOneGroup,
rec(sourceMap:=Smap,
targetMap:=Smap));

end);
##################################################################
##################################################################


##################################################################
##################################################################
InstallGlobalFunction(NormalSubgroupAsCatOneGroup,
function(G,N)
local
        C,hom,Im,ImGrp,
        Smap,Tmap,
        CmapG;

Im:=List(GeneratorsOfGroup(G),x->ConjugatorAutomorphism(N,x));
ImGrp:=Group(Im);
hom:=GroupHomomorphismByImages(G,ImGrp,GeneratorsOfGroup(G),Im);

C:=SemidirectProduct(G,hom,N);
Smap:=GroupHomomorphismByImagesNC(C,C,GeneratorsOfGroup(C),
List(GeneratorsOfGroup(C),x->
  Image( Embedding(C,1), Image( Projection(C),x)   )));

############################
CmapG:=function(x);
return
PreImagesRepresentative(Embedding(C,2), Image(Embedding(C,1), Image(Projection(C),x)^-1) *x);
end;
############################

Tmap:=GroupHomomorphismByImagesNC(C,C,GeneratorsOfGroup(C),
List(GeneratorsOfGroup(C),x->
  Image(Embedding(C,1),CmapG(x))*
  Image( Embedding(C,1), Image( Projection(C),x)   )));



return Objectify(HapCatOneGroup,
rec(sourceMap:=Smap,
targetMap:=Tmap));
end);
##################################################################
##################################################################

##################################################################
##################################################################
InstallGlobalFunction(QuasiIsomorph,
function(C)
local fn,D,G,s,t,L,bool,x,  K, M, Pi2, gens, newG, news, newt, epi;

##############
fn:=function(G,H);
return Order(G)>Order(H);
end;
##############

s:=C!.sourceMap;
t:=C!.targetMap;
G:=Source(s);
M:=Kernel(s);
Pi2:=Intersection(M,Kernel(t));

L:=NormalSubgroups(M);;
L:=Filtered(L,x->IsNormal(G,x));
Sort(L,fn);;

K:=PositionProperty(L,x->Order(Intersection(x,Pi2))=1);
K:=L[K];

gens:=GeneratorsOfGroup(K);
gens:=List(gens,x->x);
for x in GeneratorsOfGroup(K) do
Add(gens,Image(t,x));
od;
if Length(gens)=0 then Add(gens, Identity(K)); fi;
K:=Group(gens);


   epi:=NaturalHomomorphismByNormalSubgroup(G,K);
   newG:=Image(epi);
   news:=GroupHomomorphismByFunction(newG,newG,x->
         Image(epi,Image(s,PreImagesRepresentative(epi,x)))    );
   newt:=GroupHomomorphismByFunction(newG,newG,x->
         Image(epi,Image(t,PreImagesRepresentative(epi,x)))    );

   D:= Objectify(HapCatOneGroup,
         rec( sourceMap:=news,
              targetMap:=newt,
         ));

return D;

end);
##################################################################
##################################################################

##################################################################
##################################################################
InstallOtherMethod(Size,
"Order of the underlying group of a cat-1-group",
[IsHapCatOneGroup],
function(C)
return Order(Source(C!.sourceMap));
end);
##################################################################
##################################################################

##################################################################
##################################################################
InstallOtherMethod(Order,
"Order of the underlying group of a cat-1-group",
[IsHapCatOneGroup],
function(C)
return Order(Source(C!.sourceMap));
end);
##################################################################
##################################################################

#############################################
#############################################
InstallGlobalFunction(SylowSubgroupOfCatOneGroup,
function(CC,p)
local

      C,P,s,t,sp,tp;

C:=XmodToHAP(CC);
s:=C!.sourceMap;
t:=C!.targetMap;
P:=SylowSubgroup(Source(s),p);

sp:=GroupHomomorphismByImages(P,P,
GeneratorsOfGroup(P),List(GeneratorsOfGroup(P),x->Image(s,x)));
tp:=GroupHomomorphismByImages(P,P,
GeneratorsOfGroup(P),List(GeneratorsOfGroup(P),x->Image(t,x)));

#I guess we don't need the following line!!!!
if not (IsSubgroup(P,Image(sp)) and  IsSubgroup(P,Image(tp))) then return fail; fi;

return Objectify(HapCatOneGroup, rec(
                     sourceMap:=sp,
                     targetMap:=tp,
                     ));

end);
############################################
############################################

