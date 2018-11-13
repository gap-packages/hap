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

