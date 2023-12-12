#(C) Graham Ellis 

#####################################################################
#####################################################################
DeclareCategory("IsHapRegularCWComplex",IsObject);

DeclareRepresentation(  "IsHapRegularCWComplexRep",
                        IsComponentObjectRep,
                        ["nrCells",
                         "boundaries",
                         "coboundaries",
                         "vectorField"
                         ]);

HapRegularCWComplexFamily:=NewFamily( "HapRegularCWComplexFamily",
                                          IsHapRegularCWComplex,
                                          IsHapRegularCWComplex);

HapRegularCWComplex:=NewType(HapRegularCWComplexFamily,
                                IsHapRegularCWComplexRep);


InstallMethod( ViewObj,
"for HapRegularCWComplexes",
 [IsHapRegularCWComplex],
 function(R)
 Print("Regular CW-complex of dimension ",EvaluateProperty(R,"dimension"));
#if IsList(R!.vectorField) then
# Print(" with discrete vector field\n");
#else 
 Print("\n");
#fi;
end);

InstallMethod( PrintObj,
"for HapRegularCWComplexes",
 [IsHapRegularCWComplex],
 function(R)
 Print("Regular CW-complex of dimension ",EvaluateProperty(R,"dimension"));
#if IsList(R!.vectorField) then
# Print(" with discrete vector field\n");
#else
 Print("\n");
#fi;
end);
 
###################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareCategory("IsHapFilteredRegularCWComplex",IsObject);

DeclareRepresentation(  "IsHapFilteredRegularCWComplexRep",
                        IsComponentObjectRep,
                        ["nrCells",
                         "boundaries",
                         "coboundaries",
                         "vectorField",
                         "filtration",
                         ]);

HapFilteredRegularCWComplexFamily:=NewFamily( "HapFilteredRegularCWComplexFamily",
                                          IsHapFilteredRegularCWComplex,
                                          IsHapFilteredRegularCWComplex);

HapFilteredRegularCWComplex:=NewType(HapFilteredRegularCWComplexFamily,
                                IsHapFilteredRegularCWComplexRep);


InstallMethod( ViewObj,
"for HapFilteredRegularCWComplexes",
 [IsHapFilteredRegularCWComplex],
 function(R)
 Print("Filtered regular CW-complex of dimension ",EvaluateProperty(R,"dimension"));
#if IsList(R!.vectorField) then
# Print(" with discrete vector field\n");
#else
 Print("\n");
#fi;
end);

InstallMethod( PrintObj,
"for HapFilteredRegularCWComplexes",
 [IsHapFilteredRegularCWComplex],
 function(R)
 Print("Filtered regular CW-complex of dimension ",EvaluateProperty(R,"dimension"));
#if IsList(R!.vectorField) then
# Print(" with discrete vector field\n");
#else
 Print("\n");
#fi;
end);

###################################################################
#####################################################################




#####################################################################
#####################################################################
DeclareCategory("IsHapRegularCWMap",IsObject);

DeclareRepresentation(  "IsHapRegularCWMapRep",
                        IsComponentObjectRep,
                        ["source",
                         "target",
                         "mapping"
                         ]);

HapRegularCWMapFamily:=NewFamily( "HapRegularCWMapFamily",
                                          IsHapRegularCWMap,
                                          IsHapRegularCWMap);

HapRegularCWMap:=NewType(HapRegularCWMapFamily,
                                IsHapRegularCWMapRep);


InstallMethod( ViewObj,
"for HapRegularCWMaps",
 [IsHapRegularCWMap],
 function(R)
 Print("Map of regular CW-complexes\n");
end);

InstallMethod( PrintObj,
"for HapRegularCWMaps",
 [IsHapRegularCWMap],
 function(R)
 Print("Map of regular CW-complexes\n");
end);


#####################################################################
#####################################################################
DeclareCategory("IsHapEquivariantCWComplex",IsObject);


DeclareRepresentation(  "IsHapEquivariantCWComplexRep",
                        IsComponentObjectRep,
                        ["dimension",
                         "boundary",
                         "stabilizer",
                         "group",
                         "elts",
                         "properties"
                         ]);

HapEquivariantCWComplexFamily:=NewFamily( "HapEquivariantCWComplexFamily",
                                          IsHapEquivariantCWComplex,
                                          IsHapEquivariantCWComplex);

HapEquivariantCWComplex:=NewType(HapEquivariantCWComplexFamily,
                                IsHapEquivariantCWComplexRep);


InstallMethod( ViewObj,
"for HapEquivariantCWComplexes",
 [IsHapEquivariantCWComplex],
 function(R)
 Print("Equivariant CW-complex of dimension ",EvaluateProperty(R,"dimension"));
#if IsList(R!.vectorField) then
# Print(" with discrete vector field\n");
#else
 Print("\n");
#fi;
end);

InstallMethod( PrintObj,
"for HapEquivariantCWComplexes",
 [IsHapEquivariantCWComplex],
 function(R)
 Print("Equivariant CW-complex of dimension ",EvaluateProperty(R,"dimension"));
#if IsList(R!.vectorField) then
# Print(" with discrete vector field\n");
#else
 Print("\n");
#fi;
end);
###################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareCategory("IsHapEquivariantChainComplex",IsObject);


DeclareRepresentation(  "IsHapEquivariantChainComplexRep",
                        IsComponentObjectRep,
                        ["dimension",
                         "boundary",
                         "group",
                         "elts",
                         "properties"
                         ]);

HapEquivariantChainComplexFamily:=NewFamily( "HapEquivariantChainComplexFamily",
                                          IsHapEquivariantChainComplex,
                                          IsHapEquivariantChainComplex);

HapEquivariantChainComplex:=NewType(HapEquivariantChainComplexFamily,
                                IsHapEquivariantChainComplexRep);


InstallMethod( ViewObj,
"for HapEquivariantChainComplexes",
 [IsHapEquivariantChainComplex],
 function(R)
 Print("Equivariant chain complex of dimension ",EvaluateProperty(R,"dimension"));
 Print("\n");
end);

InstallMethod( PrintObj,
"for HapEquivariantChainComplexes",
 [IsHapEquivariantChainComplex],
 function(R)
 Print("Equivariant chain complex of dimension ",EvaluateProperty(R,"dimension"));
 Print("\n");
end);
###################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareCategory("IsHapEquivariantNonFreeChainComplex",IsObject);


DeclareRepresentation(  "IsHapEquivariantNonFreeChainComplexRep",
                        IsComponentObjectRep,
                        ["dimension",
                         "boundary",
                         "group",
                         "elts",
                         "stabilizer",
                         "action",
                         "properties"
                         ]);

HapEquivariantNonFreeChainComplexFamily:=NewFamily( "HapEquivariantNonFreeChainComplexFamily",
                                          IsHapEquivariantNonFreeChainComplex,
                                          IsHapEquivariantNonFreeChainComplex);

HapEquivariantNonFreeChainComplex:=NewType(HapEquivariantNonFreeChainComplexFamily,
                                IsHapEquivariantNonFreeChainComplexRep);


InstallMethod( ViewObj,
"for HapEquivariantNonFreeChainComplexes",
 [IsHapEquivariantNonFreeChainComplex],
 function(R)
 Print("Equivariant non-free chain complex of dimension ",EvaluateProperty(R,"dimension"));
 Print("\n");
end);

InstallMethod( PrintObj,
"for HapEquivariantNonFreeChainComplexes",
 [IsHapEquivariantChainComplex],
 function(R)
 Print("Equivariant non-free chain complex of dimension ",EvaluateProperty(R,"dimension"));
 Print("\n");
end);
###################################################################
#####################################################################



#####################################################################
#####################################################################
DeclareCategory("IsHapSparseChainComplex",IsObject);

DeclareRepresentation( "IsHapSparseChainComplexRep",
                        IsComponentObjectRep
                        and IsHapChain
                        and IsHapComplex,
                        ["dimension",
                         "boundary",
                         "properties"]);

HapSparseChainComplexFamily:=NewFamily( "HapSparseChainComplexFamily",
                                   IsHapSparseChainComplex,
                                   IsHapSparseChainComplex);

HapSparseChainComplex:=NewType(HapSparseChainComplexFamily,IsHapSparseChainComplexRep);

InstallMethod( ViewObj,
"for HapSparseChainComplex",
[IsHapSparseChainComplex],
 function(R)
Print("Sparse chain complex of length ",
EvaluateProperty(R,"length"), " in characteristic ",
EvaluateProperty(R,"characteristic"), " . \n");
 end);

InstallMethod( PrintObj,
"for HapSparseChainComplex",
[IsHapSparseChainComplex],
function(R)
Print("Sparse chain complex of length ",
EvaluateProperty(R,"length"), " in characteristic ",
EvaluateProperty(R,"characteristic"), " . \n");


end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
DeclareCategory("IsHapSparseChainMap",IsObject);

DeclareRepresentation(  "IsHapSparseChainMapRep",
                        IsComponentObjectRep,
                        ["source",
                        "target",
                        "mapping",
                        "properties"]);

HapSparseChainMapFamily:=NewFamily( "HapSparseChainMapFamily",
                               IsHapSparseChainMap,
                               IsHapSparseChainMap);

HapSparseChainMap:=NewType(HapSparseChainMapFamily, IsHapSparseChainMapRep);


InstallMethod( ViewObj,
"for HapSparseChainMap",
 [IsHapSparseChainMap],
 function(R)
  Print("Sparse Chain Map between complexes of length ",
   Minimum(EvaluateProperty(R!.source,"length"),
    EvaluateProperty(R!.target,"length")), " . \n");

 end);

InstallMethod( PrintObj,
"for HapSparseChainMap",
 [IsHapSparseChainMap],
 function(R)
 Print("Sparse Chain Map between complexes of length ",
 Minimum(EvaluateProperty(R!.source,"length"),
 EvaluateProperty(R!.target,"length")), " . \n");
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareCategory("IsHapFilteredSparseChainComplex",IsObject);

DeclareRepresentation( "IsHapFilteredSparseChainComplexRep",
                        IsComponentObjectRep
                        and IsHapComplex and IsHapChain and
                        IsHapSparseChainComplex,
                        ["dimension",
                         "filteredDimension",
                         "boundary",
                         "properties"]);

HapFilteredSparseChainComplexFamily:=NewFamily( "HapFilteredSparseChainComplexFamily",
                                   IsHapFilteredSparseChainComplex,
                                   IsHapFilteredSparseChainComplex);

HapFilteredSparseChainComplex:=NewType(HapFilteredSparseChainComplexFamily,IsHapFilteredSparseChainComplexRep);

InstallMethod( ViewObj,
"for HapFilteredSparseChainComplex",
[IsHapFilteredSparseChainComplex],
 function(R)
Print("Filtered sparse chain complex of length ",
EvaluateProperty(R,"length"), " in characteristic ",
EvaluateProperty(R,"characteristic"), " . \n");
 end);

InstallMethod( PrintObj,
"for HapFilteredSparseChainComplex",
[IsHapFilteredSparseChainComplex],
function(R)
Print("Filtered sparse chain complex of length ",
EvaluateProperty(R,"length"), " in characteristic ",
EvaluateProperty(R,"characteristic"), " . \n");


end);
#####################################################################
#####################################################################

