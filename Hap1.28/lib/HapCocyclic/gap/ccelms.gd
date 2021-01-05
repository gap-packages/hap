#############################################################################
##
#W  ccelms.gd                   HAPCocyclic                   Robert F. Morse
##
##

#############################################################################
##
##  Define the category of cc elements
##
DeclareCategory
(  "IsCcElement", 
   IsAssociativeElement and 
   IsMultiplicativeElementWithInverse 
);


#############################################################################
##
##  Define the representation of cc elements
##
DeclareRepresentation
(  "IsCcElementRep", 
   IsComponentObjectRep and IsAttributeStoringRep,
   [  "felement", "belement" ] 
);


#############################################################################
##
##  Basic attributes of cc elements
##
DeclareAttribute( "FibreElement",             IsCcElement );
DeclareAttribute( "BaseElement",              IsCcElement );
DeclareAttribute( "Name",                     IsCcElement );
DeclareAttribute( "InCcGroup",                IsCcElement );
 

#############################################################################
##
##  Constructor
##
DeclareGlobalFunction( "CcElement" );
