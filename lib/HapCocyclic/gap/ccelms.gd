#############################################################################
##
#W  ccelms.gd                   HAPCocyclic                   Robert F. Morse
##
##  $Id: ccelms.gd,v 1.3 2008-07-11 21:26:06 unialg Exp $
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


#############################################################################
##
##  History
##
##  $Log: ccelms.gd,v $
##  Revision 1.3  2008-07-11 21:26:06  unialg
##
##  Development commit for beta release. RFM
##
##  Revision 1.2  2008-07-11 21:00:34  unialg
##
##  First release version. RFM
##
##  Revision 1.1  2008-06-20 10:17:31  unialg
##  File Creation. RFM
##

