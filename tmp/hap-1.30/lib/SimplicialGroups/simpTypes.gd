#C Ellis & Le


#1##################################################################
#####################################################################

DeclareCategory("IsHapSimplicialGroup",IsObject);

DeclareRepresentation(  "IsHapSimplicialGroupRep",
                        IsComponentObjectRep,
                        ["groupsList",
                         "boundariesList",
						 "degeneraciesList", 
						 "properties"
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

DeclareCategory("IsHapSimplicialGroupMorphism",IsObject);

DeclareRepresentation(  "IsHapSimplicialGroupMorphismRep",
                        IsComponentObjectRep,
                        ["source",
                         "target",
						 "mapping", 
						 "properties"
                         ]);

HapSimplicialGroupMorphismFamily:=NewFamily( "HapSimplicialGroupMorphismFamily",
                                          IsHapSimplicialGroupMorphism,
                                          IsHapSimplicialGroupMorphism);

HapSimplicialGroupMorphism:=NewType(HapSimplicialGroupMorphismFamily,
                                IsHapSimplicialGroupMorphismRep);


InstallMethod( ViewObj,
"for HapSimplicialGroupMorphism",
 [IsHapSimplicialGroupMorphism],
 function(R)
 Print("Morphism of simplicial groups of length ",EvaluateProperty(R,"length"),  "\n");
end);

InstallMethod( PrintObj,
"for HapSimplicialGroupMorphism",
 [IsHapSimplicialGroupMorphism],
 function(R)
Print("Morphism of simplicial groups of length ",EvaluateProperty(R,"length"),  "\n");
 end);
 
#3##################################################################
#####################################################################

DeclareCategory("IsHapCatOneGroupMorphism",IsObject);

DeclareRepresentation(  "IsHapCatOneGroupMorphismRep",
                        IsComponentObjectRep,
                        ["source",
                         "target",
						 "mapping"
                         ]);

HapCatOneGroupMorphismFamily:=NewFamily( "HapCatOneGroupMorphismFamily",
                                          IsHapCatOneGroupMorphism,
                                          IsHapCatOneGroupMorphism);

HapCatOneGroupMorphism:=NewType(HapCatOneGroupMorphismFamily,
                                IsHapCatOneGroupMorphismRep);


InstallMethod( ViewObj,
"for HapCatOneGroupMorphism",
 [IsHapCatOneGroupMorphism],
 function(R)
 Print("Morphism of two cat-1-groups",  "\n");
end);

InstallMethod( PrintObj,
"for HapCatOneGroupMorphism",
 [IsHapCatOneGroupMorphism],
 function(R)
 Print("Morphism of two cat-1-groups",  "\n");
 end);
 
#4##################################################################
#####################################################################

DeclareCategory("IsHapCrossedModule",IsObject);

DeclareRepresentation(  "IsHapCrossedModuleRep",
                        IsComponentObjectRep,
                        ["map",
                         "action"
                         ]);

HapCrossedModuleFamily:=NewFamily( "HapCrossedModuleFamily",
                                          IsHapCrossedModule,
                                          IsHapCrossedModule);

HapCrossedModule:=NewType(HapCrossedModuleFamily,
                                IsHapCrossedModuleRep);


InstallMethod( ViewObj,
"for HapCrossedModule",
 [IsHapCrossedModule],
 function(R)
 Print("Crossed module with group homomorphism ",R!.map,"\n");
end);

InstallMethod( PrintObj,
"for HapCrossedModule",
 [IsHapCrossedModule],
 function(R)
 Print("Crossed module with group homomorphism ",R!.map,"\n");
end);
 
#5##################################################################
#####################################################################

DeclareCategory("IsHapCrossedModuleMorphism",IsObject);

DeclareRepresentation(  "IsHapCrossedModuleMorphismRep",
                        IsComponentObjectRep,
                        ["source",
                         "target",
						 "mapping"
                         ]);

HapCrossedModuleMorphismFamily:=NewFamily( "HapCrossedModuleMorphismFamily",
                                          IsHapCrossedModuleMorphism,
                                          IsHapCrossedModuleMorphism);

HapCrossedModuleMorphism:=NewType(HapCrossedModuleMorphismFamily,
                                IsHapCrossedModuleMorphismRep);


InstallMethod( ViewObj,
"for HapCrossedModuleMorphism",
 [IsHapCrossedModuleMorphism],
 function(R)
 Print("Morphism of two crossed modules",  "\n");
end);

InstallMethod( PrintObj,
"for HapCrossedModuleMorphism",
 [IsHapCrossedModuleMorphism],
 function(R)
 Print("Morphism of two crossed modules",  "\n");
 end);
 
 
 
 
 


#1##################################################################
#####################################################################

DeclareCategory("IsHapSimplicialFreeAbelianGroup",IsObject);

DeclareRepresentation(  "IsHapSimplicialFreeAbelianGroupRep",
                        IsComponentObjectRep,
                        ["dimension",
                         "boundary",
                         "properties"
                         ]);

HapSimplicialFreeAbelianGroupFamily:=NewFamily( "HapSimplicialFreeAbelianGroupFamily",
                                          IsHapSimplicialFreeAbelianGroup,
                                          IsHapSimplicialFreeAbelianGroup);

HapSimplicialFreeAbelianGroup:=NewType(HapSimplicialFreeAbelianGroupFamily,
                                IsHapSimplicialFreeAbelianGroupRep);


InstallMethod( ViewObj,
"for HapSimplicialFreeAbelianGroup",
 [IsHapSimplicialFreeAbelianGroup],
 function(R)
 Print("Simplicial free abelian group of length ",EvaluateProperty(R,"length"),  "\n");
end);

InstallMethod( PrintObj,
"for HapSimplicialFreeAbelianGroup",
 [IsHapSimplicialFreeAbelianGroup],
 function(R)
Print("Simplicial free abelian group of length ",EvaluateProperty(R,"length"),  "\n");
 end);
#2##################################################################
#####################################################################
 
