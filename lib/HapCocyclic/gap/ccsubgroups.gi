#############################################################################
##
##  cocyclic.gi                HAPcocyclic                    
##
##

#############################################################################
##
##
InstallMethod( TrivialSubgroup, 
    [IsCcGroup], 
    function( G )
        local ccg,one,F,B,T,type;

F:=Group(One(HapFibre(G)));
B:=Group(One(Base(G)));

        ccg := rec();
        type := NewType
                (   CollectionsFamily(ElementsFamily(G)),
                    IsGroup and
                    IsCcGroup and
                    IsComponentObjectRep and
                    IsAttributeStoringRep
                );

T:=     ObjectifyWithAttributes
        (  ccg,                type,
           ElementsFamily,     ElementsFamily(G),
           HapFibre,           F,
           Base,               B,
           OuterGroup,         OuterGroup(G),
           Cocycle,            Cocycle(G),
           ParentAttr,         G
        );

Order(T);
SetGeneratorsOfGroup(T,[CcElement( ElementsFamily(G), One(F), One(B), G )]);
return T;

    end );

#############################################################################
##
##
InstallMethod( Fibre,
    [IsCcGroup],
    function( G )
        local F, B, T, type, ccg;

F:=HapFibre(G);
B:= Group(One(Base(G)));

        ccg := rec();
        type := NewType
                (   CollectionsFamily(ElementsFamily(G)),
                    IsGroup and
                    IsCcGroup and
                    IsComponentObjectRep and
                    IsAttributeStoringRep
                );

T:=     ObjectifyWithAttributes
        (  ccg,                type,
           ElementsFamily,     ElementsFamily(G),
           HapFibre,           F,
           Base,               B,
           OuterGroup,         OuterGroup(G),
           Cocycle,            Cocycle(G),
           ParentAttr,         G
        );

Order(T);
return T;

    end );


#############################################################################
##
##
InstallMethod( NaturalHomomorphismOntoBase,
    [IsCcGroup],
    function( G )
        local B, hom;
B:=Base(G);
hom:=GroupHomomorphismByFunction(G,B,x->BaseElement(x));
SetKernelOfMultiplicativeGeneralMapping(hom,Fibre(G));
return hom;
end );

#############################################################################
##
##
InstallGlobalFunction(ResolutionFiniteCcGroup,
function(G,n)
local R,S,T,hom,gens,imgens;

R:=ResolutionGenericGroup(Base(G),n);
S:=ResolutionGenericGroup(Fibre(G),n);
hom:=NaturalHomomorphismOntoBase(G);
gens:=GeneratorsOfGroup(G);
imgens:=List(gens,x->Image(hom,x));
T:=ResolutionFiniteExtension(gens,imgens,R,n,false,S);
return T;
end);

#############################################################################
##
##
InstallGlobalFunction(ResolutionInfiniteCcGroup,
function(G,n)
local R,S,T,hom,preimage;

R:=ResolutionGenericGroup(Base(G),n);
S:=ResolutionGenericGroup(HapFibre(G),n);
S!.group:=Fibre(G);
Apply(S!.elts,x->CcElement( ElementsFamily(G), x, One(Base(G)),G ));
hom:=NaturalHomomorphismOntoBase(G);

##############
preimage:=function(x);
return CcElement(ElementsFamily(G), One(HapFibre(G)), x, G  );

end;
##############
T:=ResolutionExtension(hom,S,R,false,preimage);
return T;
end);

#############################################################################
##
##
InstallOtherMethod( IsomorphismFpGroup,
    [IsCcGroup],
    function( G )
        local R, P, F, iso1, iso2, iso, gens, imgens;
if IsFinite(G) then
iso1:=IsomorphismPermGroup(G);
iso2:=IsomorphismFpGroup(Image(iso1));
gens:=GeneratorsOfGroup(G);
imgens:=List(gens, x-> Image(iso2, Image(iso1,x)));
iso:=GroupHomomorphismByImagesNC(G,Range(iso2),gens,imgens);
return iso;
fi;

R:=ResolutionInfiniteCcGroup(G,2);
P:=PresentationOfResolution(R);
F:=P.freeGroup/P.relators;
gens:=R!.elts{P.gens};
imgens:=GeneratorsOfGroup(F);  #Am I sure that gens of F 
                               #are the "same" as those of P
iso:=GroupHomomorphismByImagesNC(G,F,gens,imgens);
#WARNING: The implementation of this isomorphism is not finished: we need to
#         implement the expression of an arbitrary element of G as a word 
#         in the generators. The contracting homotopy provides this expression.
return iso;
end );


