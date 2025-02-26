#############################################################################
##
##  cocyclic.gi                HAPcocyclic                    Robert F. Morse
##
##

#############################################################################
##
##  Constructors for Cc-groups 
##
##
InstallMethod
(  CcGroup,
   "Create a CcGroup via basic components",
   [ IsGOuterGroup, IsStandardNCocycle ],

   function( OA, SCo)

       local G,            ##  Cc group to be constructed
             type,         ##  Unique type of Cc-group
             gens,         ##  Generators of group
             elmsfam,      ##  Family of elements for the group
             f,b,          ##  Index elements
             modSCo,       ##  f+SCo with f a nonabelian cocycle  #GRAHAM
             nafn, fn,     ##  nonabelian cocycle nafn:(B,B)-->N, abelian 
             newfn;        ##  cocycle fn:(B,B)-->A  and newfn:(B,B)-->N 
                                                                  #GRAHAM

#############GRAHAM#####################
if IsAbelian(OA!.ActedGroup) then
modSCo:=SCo;
else
fn:=SCo!.Mapping;
nafn:=OA!.nonabeliancocycle;
newfn:=function(x,y); return nafn(x,y)*fn(x,y);end;
modSCo:=SCo;
modSCo!.Mapping:=newfn;
fi;
#############GRAHAM#####################

       ##  Create elements family and type of the group
       ## 
       elmsfam := NewFamily("cce", IsCcElement);

       type := NewType
               ( CollectionsFamily(elmsfam),
                 IsGroup and 
                 IsCcGroup and
                 IsComponentObjectRep and 
                 IsAttributeStoringRep
               );

       ##  Construct a Cc-group with attribute known so far
       ##
       G := rec();
       ObjectifyWithAttributes
       (  G,              type, 
          Base,           ActingGroup(OA),
          HapFibre,       ActedGroup(OA),
          OuterGroup,     OA,
          Cocycle,        modSCo,
          ElementsFamily, elmsfam
       );

       ##  Set the multiplicative identity
       ##
       SetOne( G,
               CcElement
               ( elmsfam,
                 Mapping(SCo)(One(Base(G)), One(Base(G)))^-1,
                 One(Base(G)), 
                 G
             ) ); 
  
       ##  Set generators of the group  
       ##
       gens := [];
       for f in Filtered( GeneratorsOfGroup(HapFibre(G)), 
                          x->not x = One(HapFibre(G)) ) 
       do
           Add( gens, 
                CcElement
                ( elmsfam,
                  f,
                  One(Base(G)),
                  G
                ) );
       od;

       for b in Filtered( GeneratorsOfGroup(Base(G)),
                          x->not x = One(Base(G)) ) 
       do
           Add( gens, 
                CcElement
                ( elmsfam,
                  One(HapFibre(G)),
                  b,
                  G
                ) );
       od;

       SetGeneratorsOfGroup(G, gens);
       Size(G);

       return G; 
    end
);

#############################################################################
##
##  Construct an "empty" Cc-group. Only attribute set is
##    the family for elements in this group.  
##
InstallMethod
(  CcGroup,
   "Create an empty CcGroup",
   [ ],
   function( )

       local G,            ##  Cc group to be constructed
             elmsfam,      ##  Family of elements for the group
             type;         ##  Unique type for this group


       ##  Create elements family and type of the group
       ## 
       elmsfam      := NewFamily("cce", IsCcElement);

       type := NewType
               (  CollectionsFamily(elmsfam),
                  IsGroup and
                  IsCcGroup and
                  IsComponentObjectRep and
                  IsAttributeStoringRep
               );

        ##  Create group object
        ## 
        G := rec();
        ObjectifyWithAttributes
        (  G,              type,
           ElementsFamily, elmsfam
        );

        return G;
    end );

#############################################################################
##
## IdGroup
##
InstallMethod
(  IdGroup,
   true,
   [ IsCcGroup ],
   0,

   function( Cc )
       return IdGroup( Image(IsomorphismPermGroup(Cc)) );
   end
);

#############################################################################
##
#M  AsList Method to create a list of elements of a CcGroup
##
##
InstallMethod
(  AsList, 
   true, 
   [ IsCcGroup ], 
   0,

   function( Cc )
       local lst,        ##  List of elements
             f,b;        ##  HapFibre and base index elements

       lst :=[];
       for f in HapFibre(Cc) do
           for b in Base(Cc) do
                
               Add( lst, 
                    CcElement
                    ( FamilyObj(One(Cc)),
                      f,
                      b,
                      Cc
                    ) );
           od; 
       od; 
       return lst; 
   end 
);


#############################################################################
##
#M  AsSSortedList
##
##      Returns the list of elements of a Cc-group as a set
##
InstallMethod
(  AsSSortedList, 
   true, 
   [IsCcGroup], 
   0,

   Cc -> AsSet(AsList(Cc)) 
);


#############################################################################
##
##  Size of a Cc-group
##
InstallMethod
(  Size, 
   true, 
   [ IsCcGroup ], 
   0,

   function( G )
       ##  Size of the group
       ##
       if not IsFinite(Base(G)) then
           return Size(Base(G));
       fi;
       if not IsFinite(HapFibre(G)) then
           return Size(HapFibre(G));
       fi;

       return Size(Base(G))*Size(HapFibre(G));
         
    end 
); 

#############################################################################
##
##  PrintObj and ViewObj Methods for Cc-groups.
##
InstallMethod
(  PrintObj, 
   true, 
   [ IsCcGroup ], 
   0,

   function( G )
       if HasSize(G) then
           Print("<Cc-group of Size ",Size(G),">");
       else
           Print("<Cc-group>");
       fi;
   end 
);

InstallMethod
(  ViewObj, 
   true, 
   [ IsCcGroup ], 
   SUM_FLAGS,

   function( G )
       if HasSize(G) then
           Print("<Cc-group of Size ",Size(G),">");
       else
           Print("<Cc-group>");
       fi;
   end 
);


####Added below January 2025

#####################################################
#####################################################
InstallGlobalFunction(HAP_IsomorphismCcFpGroup,
function(G)
local Q, N, IsoQ, IsoN, F, FQ, FN, FFQ, FFN, relsQ, relsN, relsG,
FFQhomF, FQhomF,FFNhomF, FFQmappingG,gensF, gensFQ, gensFFQ, gensFFN, imN,
gensFN, FNhomF, gensG, FG, r,x,y,i,xx,yy,zz;

Q:=G!.Base;
N:=G!.HapFibre;
IsoQ:=IsomorphismFpGroup(Q);
IsoN:=IsomorphismFpGroup(N);
FQ:=Range(IsoQ);
FFQ:=FreeGroupOfFpGroup(FQ);
FN:=Range(IsoN);
FFN:=FreeGroupOfFpGroup(FN);
gensFFQ:=GeneratorsOfGroup(FFQ);
gensFQ:=GeneratorsOfGroup(FQ);
gensFFN:=GeneratorsOfGroup(FFN);
gensFN:=GeneratorsOfGroup(FN);
F:=FreeGroup(Length(gensFFN)+Length(gensFFQ));
gensF:=GeneratorsOfGroup(F);
FFQhomF:=GroupHomomorphismByImagesNC(FFQ,F,gensFFQ,gensF{[Length(gensFFN)+1..Length(gensFFN)+Length(gensFFQ)]});
FFNhomF:=GroupHomomorphismByImagesNC(FFN,F,gensFFN,gensF{[1..Length(gensFFN)]});
FNhomF:=GroupHomomorphismByImagesNC(FN,F,gensFN,gensF{[1..Length(gensFN)]});
FQhomF:=GroupHomomorphismByImagesNC(FQ,F,gensFQ,gensF{[1+Length(gensFN)..Length(gensFQ)+Length(gensFN)]});
relsQ:=RelatorsOfFpGroup(FQ);
relsN:=RelatorsOfFpGroup(FN);
relsG:=List(relsN,r->Image(FFNhomF,r));

imN:=[];
for i in [1..Length(gensFQ)] do
x:=PreImagesRepresentative(IsoQ,gensFQ[i]);
x:=CcElement( FamilyObj(One(G)),One(N),x,InCcGroup(One(G)));
Add(imN, x);
od;

FFQmappingG:=GroupHomomorphismByImagesNC(FFQ,G,gensFFQ,imN   );
for r in relsQ do
x:=Image(FFQhomF,r);
y:=Image(FFQmappingG,r);
y:=y!.FibreElement;
y:=Image(IsoN,y);
y:=Image(FNhomF,y);

Add(relsG,Image(FFQhomF,r)*y^-1);
od;

for x in gensFQ do
for y in gensFN do

xx:=PreImagesRepresentative(IsoQ,x);;
xx:=CcElement( FamilyObj(One(G)),One(N),xx,InCcGroup(One(G)));
yy:=PreImagesRepresentative(IsoN,y);;
yy:=CcElement( FamilyObj(One(G)),yy,One(Q),InCcGroup(One(G)));
zz:=xx*yy*xx^-1;
zz:=zz!.FibreElement;
zz:=Image(IsoN,zz);
Add(relsG, Image(FQhomF,x)*Image(FNhomF,y)*Image(FQhomF,x)^-1*Image(FNhomF,zz)^-1);
od;od;

FG:= F/relsG;

gensG:=[];
for x in gensFN do
xx:=PreImagesRepresentative(IsoN,x);
Add(gensG, CcElement( FamilyObj(One(G)),xx,One(Q),InCcGroup(One(G))) );
od;
for x in gensFQ do
xx:=PreImagesRepresentative(IsoQ,x);
Add(gensG, CcElement( FamilyObj(One(G)),One(N),xx,InCcGroup(One(G))) );
od;

return GroupHomomorphismByImagesNC(G,FG,gensG,GeneratorsOfGroup(FG));
end);
#####################################################
#####################################################


#####################################################
#####################################################
InstallMethod(IsomorphismFpGroup,
"Isomorphism from cocyclic to fp group",
[IsCcGroup],
1000000,
function(G)
return HAP_IsomorphismCcFpGroup(G);
end);
#####################################################
#####################################################


#####################################################
#####################################################
InstallOtherMethod(IsomorphismPermGroup,
"Isomorphism from cocyclic to perm group",
[IsCcGroup],
function(G)
local isoFp, isoperm, iso, gens, ims;
isoFp:= HAP_IsomorphismCcFpGroup(G);
isoperm:=IsomorphismPermGroup(Range(isoFp));
gens:=GeneratorsOfGroup(G);
ims:=List(gens,x->Image(isoperm,Image(isoFp,x)));

return GroupHomomorphismByImages(G,Target(isoperm),gens,ims);
end);
#####################################################
#####################################################

