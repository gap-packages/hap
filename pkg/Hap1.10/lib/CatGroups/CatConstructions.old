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
  ImageElm( Embedding(C,1), ImageElm( Projection(C),x)   )));

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
InstallGlobalFunction(QuotientQuasiIsomorph,
function(C)
local fn,D,G,s,t,L,bool,x, LL, K, M, Pi2, gens, newG, news, newt, epi,
newGgens;


s:=C!.sourceMap;
t:=C!.targetMap;
G:=Source(s);
M:=KernelWG(s);
Pi2:=Intersection(M,KernelWG(t));

L:=NormalSubgroups(M);;
L:=Filtered(L,x->IsNormal(G,x));

L:=Filtered(L,x->Order(Intersection(x,Pi2))=1);



#############################
fn:=function(K)
local gens,x;
gens:=GeneratorsOfGroup(K);
gens:=List(gens,x->x);
for x in GeneratorsOfGroup(K) do
Add(gens,Image(t,x));
od;
if Length(gens)=0 then Add(gens, Identity(K)); fi;
K:=Group(gens);
return K;
end;
############################
L:=List(L,x->fn(x));
LL:=List(L,x->Order(x));
K:=L[Position(LL,Maximum(LL))];

   epi:=NaturalHomomorphismByNormalSubgroup(G,K);
   newG:=Image(epi);
   newGgens:=GeneratorsOfGroup(newG);
   if Length(newGgens)=0 then newGgens:=[Identity(newG)]; fi;
   news:=GroupHomomorphismByImagesNC(newG,newG,newGgens,
         List(newGgens, x->
         Image(epi,Image(s,PreImagesRepresentative(epi,x))))    );
   newt:=GroupHomomorphismByImagesNC(newG,newG,newGgens,
         List(newGgens,x->
         Image(epi,Image(t,PreImagesRepresentative(epi,x))))    );

   D:= Objectify(HapCatOneGroup,
         rec( sourceMap:=news,
              targetMap:=newt,
         ));

return D;

end);
##################################################################
##################################################################


###########################################################
###########################################################
InstallGlobalFunction(SubQuasiIsomorph,
function(C)
local n,s,t,G,L,LS,CC,c,g,Kers,K,H,news,newt,D,
SortOrder,CheckConditions,newHgens;
########################
SortOrder:=function(H,G)
        return Order(G)>Order(H);
end;
#######################
n:=Order(HomotopyGroup(C,1));
s:=C!.sourceMap;
t:=C!.targetMap;
G:=Source(s);
########################
CheckConditions:=function(H)
        return IsSubgroup(H,Image(t,H)) and IsSubgroup(H,Image(s,H)) and IsSubgroup(H,Intersection(Kernel(t),Kernel(s)));
end;
#######################
Kers:=KernelWG(s);
#L:=NormalSubgroups(G);
LS:=LatticeSubgroups(G);;
CC:=ConjugacyClassesSubgroups(LS);;
L:=[];;
for c in CC do
 for g in c do
 Add(L,g);
 od;
od;
Sort(L,SortOrder);

for K in L do
   if CheckConditions(K) then
        if Order( Image(s,K) / Image(t,Intersection(K,Kers)) )=n 
        and Order(Image(s,K) / Intersection(Image(s,K), Image(t, Kers)))=n 
        then
                        H:=K;
                        break;
        fi;
   fi;
od;
newHgens:=GeneratorsOfGroup(H);
if Length(newHgens)=0 then newHgens:=[Identity(H)]; fi;
news:=GroupHomomorphismByImagesNC(H,H,newHgens, 
  List(newHgens, x-> Image(s,x))   );
newt:=GroupHomomorphismByImagesNC(H,H,newHgens, 
  List(newHgens, x-> Image(t,x))   );

D:=Objectify(HapCatOneGroup, rec(
        sourceMap:=news,
        targetMap:=newt));
return D;
end);
#######################################################
#######################################################

#######################################################
#######################################################
InstallGlobalFunction(QuasiIsomorph,
function(C)
local
 	D;
D:=SubQuasiIsomorph(C);
D:=QuotientQuasiIsomorph(D);

while Size(D)<Size(C) do
C:=D;
D:=SubQuasiIsomorph(C);
D:=QuotientQuasiIsomorph(D);
od;

return D;

end);
#######################################################
#######################################################


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
    C,G,P,gens,s,t,sp,tp,
        Num,i,k;
    C:=XmodToHAP(CC);
    s:=C!.sourceMap;
    t:=C!.targetMap;
    G:=Source(s);
    P:=SylowSubgroup(Image(s,G),p);
        k:=1;
    Num:=Factors(Order(G));
    for i in Num do
        if i=p then
            k:=k*p;
        fi;
    od;
    while Order(P)<k do
        P:=SylowSubgroup(Normalizer(G,P),p);
    od;
        gens:=GeneratorsOfGroup(P);
    sp:=GroupHomomorphismByImages(P,P,gens,List(gens,x->Image(s,x)));
    tp:=GroupHomomorphismByImages(P,P,
    gens,List(gens,x->Image(t,x)));
    return Objectify(HapCatOneGroup, rec(
                         sourceMap:=sp,
                         targetMap:=tp,
                         ));
end);
############################################
############################################

