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

 
