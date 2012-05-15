#C Ellis 

#####################################################################
#####################################################################
DeclareCategory("IsHapRegularCWSpace",IsObject);

DeclareRepresentation(  "IsHapRegularCWSpaceRep",
                        IsComponentObjectRep,
                        ["nrCells",
                         "boundaries",
                         "coboundaries",
                         "vectorField"
                         ]);

HapRegularCWSpaceFamily:=NewFamily( "HapRegularCWSpaceFamily",
                                          IsHapRegularCWSpace,
                                          IsHapRegularCWSpace);

HapRegularCWSpace:=NewType(HapRegularCWSpaceFamily,
                                IsHapRegularCWSpaceRep);


InstallMethod( ViewObj,
"for HapRegularCWSpaces",
 [IsHapRegularCWSpace],
 function(R)
 Print("Regular CW-space of dimension ",EvaluateProperty(R,"dimension"));
#if IsList(R!.vectorField) then
# Print(" with discrete vector field\n");
#else 
 Print("\n");
#fi;
end);

InstallMethod( PrintObj,
"for HapRegularCWSpaces",
 [IsHapRegularCWSpace],
 function(R)
 Print("Regular CW-space of dimension ",EvaluateProperty(R,"dimension"));
#if IsList(R!.vectorField) then
# Print(" with discrete vector field\n");
#else
 Print("\n");
#fi;
end);

 
###################################################################
#####################################################################

 
