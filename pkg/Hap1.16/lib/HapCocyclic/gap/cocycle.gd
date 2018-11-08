#############################################################################
##
#W  cocycle.gd                       HAP                      Robert F. Morse
##                                                               Graham Ellis
##
##
##  $Id: cocycle.gd,v 1.3 2008-07-11 21:41:08 unialg Exp $
##

#############################################################################
##
##  Declare a Standard 2-Cocycle
##
DeclareProperty("IsStandard2Cocycle", IsObject);

##############################################################################
##
##  Attributes of a Standard 2-Cocycle -- compatible with StandardNCocycle 
##
DeclareAttribute( "ActingGroup" ,        IsStandard2Cocycle );
DeclareAttribute( "Module" ,             IsStandard2Cocycle );
DeclareAttribute( "Mapping" ,            IsStandard2Cocycle );

##############################################################################
##
##  Empty Constructor 
##
DeclareOperation("Standard2Cocycle", []);

##############################################################################
##
##  Simple Constructor  
##
DeclareOperation("Standard2Cocycle", [IsGroup, IsGroup]);

##############################################################################
##
##
##  $Log: cocycle.gd,v $
##  Revision 1.3  2008-07-11 21:41:08  unialg
##
##
##  Development commit for beta release. RFM
##
##  Revision 1.2  2008-07-11 21:02:25  unialg
##
##  First beta release. RFM
##
##  Revision 1.1  2008-06-18 06:06:34  unialg
##  File Creation. RFM
##
##
##


