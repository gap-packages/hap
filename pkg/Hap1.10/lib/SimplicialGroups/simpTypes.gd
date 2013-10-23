#C Ellis & Le

#1####################################################################
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
 
#2##################################################################
#####################################################################

DeclareCategory("IsHapSimplicialGroupMap",IsObject);

DeclareRepresentation(  "IsHapSimplicialGroupMapRep",
                        IsComponentObjectRep,
                        ["sourceMap",
                         "targetMap",
						 "mapping"
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
 
#3##################################################################
#####################################################################

DeclareCategory("IsHapCatOneGroupHomomorphism",IsObject);

DeclareRepresentation(  "IsHapCatOneGroupHomomorphismRep",
                        IsComponentObjectRep,
                        ["source",
                         "target",
						 "mapping"
                         ]);

HapCatOneGroupHomomorphismFamily:=NewFamily( "HapCatOneGroupHomomorphismFamily",
                                          IsHapCatOneGroupHomomorphism,
                                          IsHapCatOneGroupHomomorphism);

HapCatOneGroupHomomorphism:=NewType(HapCatOneGroupHomomorphismFamily,
                                IsHapCatOneGroupHomomorphismRep);


InstallMethod( ViewObj,
"for HapCatOneGroupHomomorphism",
 [IsHapCatOneGroupHomomorphism],
 function(R)
 Print("Homormophism of two cat one groups",  "\n");
end);

InstallMethod( PrintObj,
"for HapCatOneGroupHomomorphism",
 [IsHapCatOneGroupHomomorphism],
 function(R)
 Print("Homormophism of two cat one groups",  "\n");
 end);
 


 
