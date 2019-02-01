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
InstallMethod( CcGroup,
    "Create a CcGroup via basic components",
    [ IsGOuterGroup, IsStandard2Cocycle ],

    function( OA, SCo)

        local G,            ##  Cc group to be constructed
              type,         ## 
              gens,         ## 
              elmsfam,      ##
              f,b;          ##  Index elements

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
        ( G,              type, 
          Base,           ActingGroup(OA),
          Fibre,          ActedGroup(OA),
          FibreMapping,   IdentityMapping(ActedGroup(OA)),
          OuterGroup,     OA,
          Cocycle,        SCo,
          ElementsFamily, elmsfam
        );

        ##  Set the multiplicative identity
        ##
        SetOne( G,
                CcElement
                ( elmsfam,
                  Cocycle(SCo)(One(Base(G)), One(Base(G)))^-1,
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
    end);

##  Construct an "empty" Cc-group. Only attribute set is
##    the family for elements in this group.  
##
InstallMethod(CcGroup,
    "Create an empty CcGroup",
    [ ],
    function( )

        local G,            ##  Cc group to be constructed
              elmsfam,      ##
              type;

        elmsfam      := NewFamily("cce", IsCcElement);

        type := NewType
                (   CollectionsFamily(elmsfam),
                    IsGroup and
                    IsCcGroup and
                    IsComponentObjectRep and
                    IsAttributeStoringRep
                );

        G := rec();
        ObjectifyWithAttributes
        (  G,              type,
           ElementsFamily, elmsfam
        );

        return G;
    end );


InstallMethod
( CcGroup,
  "construct a new group with the same elements family",
  [ IsObject ],
  function( elmsfam )
      local G, type;
      type := NewType
              (   CollectionsFamily(elmsfam),
                  IsGroup and
                  IsCcGroup and
                  IsComponentObjectRep and
                  IsAttributeStoringRep
              );
      G := rec();
      ObjectifyWithAttributes
      (  G,              type,
         ElementsFamily, elmsfam
      );
      return G;  

  end
);
#############################################################################
##
#M  AsList Method to create a list of Element of a CcGroup
##
##
InstallMethod(AsList, true, [IsCcGroup], 0,
    function( Cc )
        local n,g, lst;

        lst :=[];
        for n in Fibre(Cc) do
            for g in Base(Cc) do
                
                Add( lst, 
                     CcElement
                     ( FamilyObj(One(Cc)),
                       n,
                       g,
                       Cc
                   ) );
            od; 
        od; 
        return lst; 
    end );

#############################################################################
##
#M  AsSSortedList
##
##      Returns the list of elements of a Cc-group as a set
##
InstallMethod(AsSSortedList, true, [IsCcGroup], 0,
    Cc -> AsSet(AsList(Cc)) );

InstallMethod(MagmaWithInversesByGenerators, true, [IsCcGroup, IsList], 0,
    function(ccg,gens)
        return "this is a joke";
    end );


#############################################################################
##
#M  Centre
##      
##      Returns the centre of a Cc-group      
##
##      Method: Take
##
InstallMethod(Centre, true, [IsCcGroup], 0,
    function(G)
        local Zgens,
              ZF, ZB, 
              B, F, fam,
              ccg, type,
              cocyc, fm,
              act;

        B     := Base(G);
        F     := Fibre(G);
        fam   := ElementsFamily(G);
        cocyc := Cocycle(Cocycle(G));

        ZB := Centre(B);
        
        Zgens := List(Cartesian(F,ZB), x->CcElement(fam, x[1], x[2], G));
        Zgens := Filtered(Zgens, x->ForAll(GeneratorsOfGroup(G), y->IsOne(Comm(x,y))));
        ZB    := Subgroup(B, List(Zgens, BaseElement));

        ##  This subgroup my be too large
        ##
        ZF    := Subgroup(F, List(Zgens, FibreElement)); 
        fm    := NaturalHomomorphismByNormalSubgroup
                 (  ZF, 
                    Subgroup
                    ( ZF, 
                      List(Cartesian(ZB,ZB), x->cocyc(x[1],x[2]))
                    )
                 );
        ZF    := Image(fm);
        ccg := rec();
        type := NewType
                (   CollectionsFamily(ElementsFamily(G)),
                    IsGroup and
                    IsCcGroup and
                    IsComponentObjectRep and
                    IsAttributeStoringRep
                );

        ObjectifyWithAttributes
        (  ccg,                type,
           ElementsFamily,     ElementsFamily(G),
           Fibre,              ZF,
           Base,               ZB,
           FibreMapping,       fm,
           OuterGroup,         OuterGroup(G),
           Cocycle,            Cocycle(G),
           ParentAttr,         G
        );                                

        Zgens := List
                 ( GeneratorsOfGroup(ZF), 
                   x->CcElement
                      ( ElementsFamily(G), 
                        x, 
                        One(ZB), 
                        ccg 
                      )
                 );
        Append( Zgens, List
                       ( GeneratorsOfGroup(ZB), 
                         x->CcElement
                            ( ElementsFamily(G), 
                              One(ZF), 
                              x, 
                              ccg 
                            )
                        )
                );
        SetOne(ccg, CcElement( ElementsFamily(G),   
                               One(ZF),
                               One(ZB),
                               ccg
                             ) );
              
        SetGeneratorsOfGroup(ccg,Zgens);
        Size(ccg);

        return ccg;

    end );

InstallMethod( TrivialSubgroup, true, [IsCcGroup], 0,
    function( ccg )
        local G,one;

        G   := CcGroup(ElementsFamily(ccg));

        one := CcElement
               (  ElementsFamily(ccg), 
                  One(Fibre(ccg)), 
                  One(Base(ccg)), 
                  ccg 
               );

        SetFibre             ( G, TrivialSubgroup(Fibre(ccg)) );
        SetBase              ( G, TrivialSubgroup(Base(ccg))  );
        SetFibreMapping      ( G, IdentityMapping(TrivialSubgroup(Fibre(ccg))) );
        SetOuterGroup        ( G, OuterGroup(ccg)             );
        SetCocycle           ( G, Cocycle(ccg)                );
        SetParentAttr        ( G, ccg                         );        
        SetGeneratorsOfGroup ( G, [one]                       );
        SetOne               ( G, one                         );
        SetSize              ( G, 1                           );

        return G;

    end );


InstallMethod( IdGroup, true, [ IsCcGroup ], 0,
    function( G )
        return IdGroup
               ( Image
                 ( IsomorphismPermGroup
                   ( GroupByMultiplicationTable
                     ( MultiplicationTable
                       ( AsList(G) )
                 ) ) )
               );   
    end );

#############################################################################
##
##  Size of a Cc-group
##
InstallMethod( Size, true, [ IsCcGroup ], 0,
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
         
    end ); 

#############################################################################
##
##  PrintObj and ViewObj Methods for Cc-groups.
##
InstallMethod( PrintObj, true, [ IsCcGroup ], 0,
    function( G )
        if HasSize(G) then
            Print("<Cc-group of Size ",Size(G),">");
        else
            Print("<Cc-group>");
        fi;
    end );

InstallMethod( ViewObj, true, [ IsCcGroup ], SUM_FLAGS,
    function( G )
        if HasSize(G) then
            Print("<Cc-group of Size ",Size(G),">");
        else
            Print("<Cc-group>");
        fi;
    end );
