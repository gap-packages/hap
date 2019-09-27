#############################################################################
##
#W  cocycle.gd                       HAP                      Robert F. Morse
##                                                               Graham Ellis
##
##
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
