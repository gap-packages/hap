#############################################################################
##
##   ccelms.gi                   HAPCocyclic                  Robert F. Morse
##
##

#############################################################################
##
##   Create a CcElement
##
InstallGlobalFunction
(  CcElement,
   function
   (  Fam,          ## Elements family    
      f, b,         ## fibre and base components
      Ccg           ## Cc-group the element belongs to
   )
       local  elm,      ## Element to create 
              type;     ## type for the element
        
       elm  := rec();
       type := NewType( Fam, IsCcElement and IsCcElementRep );

       ##  Create element 
       ##
       ObjectifyWithAttributes
       (  elm,           type,
          FibreElement,  f,
          BaseElement,   b,
          Name,          "c",
          InCcGroup,     Ccg
       );

       return elm; 
   end 
);


#############################################################################
##
##   Multiplication method for cc-elements
##
InstallMethod
(  \*,
   "for cc elements",
   IsIdenticalObj,
   [IsCcElement, IsCcElement],
   0,

   function( l, r )
       local act,         ##  Outer action B --> F 
             cocyc,       ##  Cocycle (B,B) --> F
             f1, b1,      ##  fibre, base components
             f2, b2;      ##    of l and r respectively
        
       f1    := FibreElement(l); 
       b1    := BaseElement(l); 
      
       f2    := FibreElement(r);
       b2    := BaseElement(r); 

       act   := OuterAction(OuterGroup(InCcGroup(l)));
       cocyc := Mapping(Cocycle(InCcGroup(l)));     

       return CcElement
              ( FamilyObj(l), 
                #f1*act(b1^-1,f2)*cocyc(b1,b2), 
		f1*act(b1,f2)*cocyc(b1,b2),   #Modified by Graham 
                b1*b2, 
                InCcGroup(l) 
              );

   end 
);

#############################################################################
##
##   Inverse method for cc-elements
##
##   (f,b) = (f,1)(1,b)  so 
##   (f,b)^-1 = (1,b)^-1 * (f,1)^-1 = ( (b,b^-1), b^-1 ) * (f^-1 (1,1), 1) 
##
InstallMethod
(  InverseOp,
   "for ccp elements",
   true,
   [IsCcElement],
   0,

   function( g )
       local cocyc,       ##  Cocycle  BxB --> F
             f, b;        ##  fibre, base components of g
       
       cocyc := Mapping(Cocycle(InCcGroup(g)));     

       f     := CcElement
                (  FamilyObj(g), 
                   FibreElement(g)^-1 * 
                       cocyc( One(BaseElement(g)), One(BaseElement(g)) )^-1, 
                   BaseElement(One(g)), 
                   InCcGroup(g)
                );

       b     := CcElement
                (  FamilyObj(g), 
                   cocyc( BaseElement(g), BaseElement(g)^-1 )^-1, 
                   BaseElement(g)^-1, 
                   InCcGroup(g)
                );

       return b*f;

   end 
);


#############################################################################
##
##   Obtain the identity element of an element of the same family
##
InstallMethod( One,true,[IsCcElement],0,
    c -> One(InCcGroup(c)) );


#############################################################################
##
##   Determine if a Cc-element is the identity element 
##
InstallMethod( IsOne,true,[IsCcElement],0,
    c -> c = One(InCcGroup(c)) );


#############################################################################
##
##   Equality method for two CcElements.
##
##   Equality on the tuples l = (fl,bl) and r = (fr,br).
##
InstallMethod
(  \=,
   "for cc elements",
   IsIdenticalObj,
   [IsCcElement, IsCcElement],
   0,
   function( l, r )
       return ( FibreElement(l) = FibreElement(r) )  and 
              ( BaseElement(l)  = BaseElement(r)  );

   end 
);


#############################################################################
##
##   Less than method for two CcElements.
##
##   Lexicographical on the tuples l = (fl,bl) and r = (fr,br).
##
InstallMethod
(  \<,
   "for cc elements",
   IsIdenticalObj,
   [IsCcElement, IsCcElement],
   0,

   function( l, r )
       if ( FibreElement(l) < FibreElement(r) ) then return true; fi;
       if ( FibreElement(l) = FibreElement(r) ) then 
           return ( BaseElement(l) < BaseElement(r) ); 
       fi;
       return false;
   end 
);
