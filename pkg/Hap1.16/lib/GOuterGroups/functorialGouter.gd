#####################################################################
#####################################################################
DeclareCategory("IsHapGComplexMap",IsObject);

DeclareRepresentation(  "IsHapGComplexMapRep",
                        IsComponentObjectRep,
                        ["source",
                        "target",
                        "mapping",
                        "properties"]);

HapGComplexMapFamily:=NewFamily( "HapGComplexMapFamily",
                               IsHapGComplexMap,
                               IsHapGComplexMap);

HapGComplexMap:=NewType(HapGComplexMapFamily, IsHapGComplexMapRep);


InstallMethod( ViewObj,
"for HapGComplexMap",
 [IsHapGComplexMap],
 function(R)
  Print("GComplex Map between complexes of length ",
   Minimum(EvaluateProperty(R!.source,"length"),
    EvaluateProperty(R!.target,"length")), " . \n");

 end);

InstallMethod( PrintObj,
"for HapGComplexMap",
 [IsHapGComplexMap],
 function(R)
 Print("GComplex Map between complexes of length ",
 Minimum(EvaluateProperty(R!.source,"length"),
 EvaluateProperty(R!.target,"length")), " . \n");
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareGlobalFunction("HomToGModule_hom");
DeclareGlobalFunction("CohomologyModule_Gmap");
DeclareGlobalFunction("CommonEndomorphisms");
DeclareGlobalFunction("CohomologyModule_AsAutModule");
