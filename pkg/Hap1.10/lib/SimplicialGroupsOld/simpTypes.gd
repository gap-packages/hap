#C Ellis & Le

#####################################################################
#####################################################################
DeclareCategory("IsHapSimplicialGroup",IsObject);

DeclareRepresentation(  "IsHapSimplicialGroupRep",
                        IsComponentObjectRep,
                        ["sourceMap",
                         "targetMap"
                         ]);

HapSimplicialGroupFamily:=NewFamily( "HapSimplicialGroupFamily",
                                          IsHapSimplicialGroup,
                                          IsHapSimplicialGroup);

HapSimplicialGroup:=NewType(HapSimplicialGroupFamily,
                                IsHapSimplicialGroupRep);


InstallMethod( ViewObj,
"for HapSimplicialGroup",
 [IsHapSimplicialGroup],
 function(R)
 Print("Simplicial group of length ",EvaluateProperty(R,"length"),  "\n");
end);

InstallMethod( PrintObj,
"for HapSimplicialGroup",
 [IsHapSimplicialGroup],
 function(R)
Print("Simplicial group of length ",EvaluateProperty(R,"length"),  "\n");
 end);
 
###################################################################
#####################################################################

DeclareCategory("IsHapSimplicialGroupMap",IsObject);

DeclareRepresentation(  "IsHapSimplicialGroupMapRep",
                        IsComponentObjectRep,
                        ["sourceMap",
                         "targetMap"
                         ]);

HapSimplicialGroupMapFamily:=NewFamily( "HapSimplicialGroupMapFamily",
                                          IsHapSimplicialGroupMap,
                                          IsHapSimplicialGroupMap);

HapSimplicialGroupMap:=NewType(HapSimplicialGroupMapFamily,
                                IsHapSimplicialGroupMapRep);


InstallMethod( ViewObj,
"for HapSimplicialGroupMap",
 [IsHapSimplicialGroupMap],
 function(R)
 Print("Simplicial group map of length ",EvaluateProperty(R,"length"),  "\n");
end);

InstallMethod( PrintObj,
"for HapSimplicialGroupMap",
 [IsHapSimplicialGroupMap],
 function(R)
Print("Simplicial group map of length ",EvaluateProperty(R,"length"),  "\n");
 end);
 
###################################################################
#####################################################################


 
