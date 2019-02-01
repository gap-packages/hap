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
          Fibre,          ActedGroup(OA),
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
       for f in Filtered( GeneratorsOfGroup(Fibre(G)), 
                          x->not x = One(Fibre(G)) ) 
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
                  One(Fibre(G)),
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
             f,b;        ##  Fibre and base index elements

       lst :=[];
       for f in Fibre(Cc) do
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
       if not IsFinite(Fibre(G)) then
           return Size(Fibre(G));
       fi;

       return Size(Base(G))*Size(Fibre(G));
         
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
