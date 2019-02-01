#############################################################################
##
#W  cocycle.gi                     HAP                        Robert F. Morse
##                                                               Graham Ellis
##

#############################################################################
##
##  Construct an 'empty" Standard 2-cocycle -- the attributes must be set
##     by the user.
##
##
InstallMethod( Standard2Cocycle,
    "create an empty cocycle",
    [ ],

    function()
        local S,       ## Standard cocycle to be returned
              type;    ## The type of S

        S := rec();
        type := NewType
                (  NewFamily("scc"),
                   IsStandard2Cocycle and IsStandardNCocycle and
                   IsComponentObjectRep and IsAttributeStoringRep
                );
        ObjectifyWithAttributes
        (  S,
           type
        );

        return S;
    end
);


#############################################################################
##
##  Construct a Standard 2-cocycle from E and N where
##     N --> E --> G
##  
##
InstallMethod
(  Standard2Cocycle,
   "basic method for creating a cocycle",
   [ IsGroup, IsGroup ],

   function( E, N )
       local S,            ##  Final Cocycle object returned
             GO,           ##  G-Outer group
             type,         ##  Type object for Cocycle
             nat,          ##  Natural homomorphism E --> E/N
             cocycle;      ##  Cocyle function E/N x E/N --> Z(N)

       ##  Check to make sure N is normal
       ##   
       if not IsNormal(E,N) then
           Error("N must be a normal subgroup of E");
       fi;

       nat := NaturalHomomorphismByNormalSubgroup(E,N);

       cocycle := 
         function(g,h) 
           return 
             Representative(PreImages(nat,g))*
             Representative(PreImages(nat,h))*
             Representative(PreImages(nat,g*h))^-1; 
         end;
         
        ##  Create the associated G-outer group
        ##   
        GO  :=  GOuterGroup();
        SetActingGroup(GO,Image(nat));
        SetActedGroup(GO,N);
        SetOuterAction
        (  GO, 
           function(g,n) 
             return n^(Representative(PreImages(nat,g))^-1); #Graham changed this 
           end 
        );

        ##  Create the cocycle -- which is both 2-Cocycle and
        ##  a N-Cocycle.
        ##
        S := rec();
        type := 
          NewType
          (  NewFamily( "scc" ),
             IsStandard2Cocycle and IsStandardNCocycle and 
             IsComponentObjectRep and IsAttributeStoringRep
          );

        ObjectifyWithAttributes
        (  S,           type,
           ActingGroup, Image(nat),
           Module,      GO, 
           Mapping,     cocycle
        );

        return S;
    end);
