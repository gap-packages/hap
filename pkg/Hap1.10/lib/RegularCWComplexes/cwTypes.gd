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

