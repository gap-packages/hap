#############################################################################
##
#W  ccgroup.gd                   HAPcocyclic                  Robert F. Morse
##
##
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
DeclareAttribute( "HapFibre",          IsCcGroup );
DeclareAttribute( "OuterGroup",     IsCcGroup );
DeclareAttribute( "Cocycle",        IsCcGroup );
DeclareAttribute( "ElementsFamily", IsCcGroup );
DeclareAttribute( "Fibre",          IsCcGroup );
DeclareAttribute( "NaturalHomomorphismOntoBase", IsCcGroup );

#############################################################################
##
##
DeclareOperation("CcGroup",      [IsGOuterGroup, IsStandardNCocycle] );
DeclareOperation("CcGroup",      [ ]                                 );

#############################################################################
##
##
DeclareGlobalFunction("ResolutionFiniteCcGroup");
DeclareGlobalFunction("ResolutionInfiniteCcGroup");

DeclareGlobalFunction("HAP_IsomorphismCcFpGroup");
