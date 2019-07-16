#############################################################################
##
#W  ccgroup.gd                   HAPcocyclic                  Robert F. Morse
##
##
##  $Id: ccgroup.gd,v 1.2 2008-07-11 21:34:33 unialg Exp $
##

#############################################################################
##
##  Declare cc groups as a type of group
##
DeclareProperty("IsCcGroup", IsGroup);


#############################################################################
##
##  Attributes of a CcGroup.
##
DeclareAttribute( "Base" ,          IsCcGroup );
DeclareAttribute( "Fibre",          IsCcGroup );
DeclareAttribute( "OuterGroup",     IsCcGroup );
DeclareAttribute( "Cocycle",        IsCcGroup );
DeclareAttribute( "ElementsFamily", IsCcGroup );

#############################################################################
##
##
DeclareOperation("CcGroup",      [IsGOuterGroup, IsStandardNCocycle] );
DeclareOperation("CcGroup",      [ ]                                 );


#############################################################################
##
#H  History
##
##  $Log: ccgroup.gd,v $
##  Revision 1.2  2008-07-11 21:34:33  unialg
##
##  Development commit for beta version. RFM
##
##  Revision 1.1  2008-07-11 21:02:00  unialg
##
##  First beta release. RFM
##
##
