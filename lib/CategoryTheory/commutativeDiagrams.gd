#(C) 2009 Graham Ellis

#####################################################################
#####################################################################
DeclareCategory("IsHapCommutativeDiagram",IsObject);

DeclareRepresentation(  "IsHapCommutativeDiagramRep",
                        IsComponentObjectRep,
                        ["objects",
			 "arrows",
                         "properties"]);

HapCommutativeDiagramFamily:=NewFamily( "HapCommutativeDiagramFamily",
                                 IsHapCommutativeDiagram,
                                 IsHapCommutativeDiagram);

HapCommutativeDiagram:=NewType(HapCommutativeDiagramFamily,IsHapCommutativeDiagramRep);


InstallMethod( ViewObj,
"for HapCommutativeDiagram",
[IsHapCommutativeDiagram],
function(T)
Print("Commutative diagram with ",Size(T!.objects)," objects and ",Size(T!.arrows)," arrows.\n"     );
end);

InstallMethod( PrintObj,
"for HapCommutativeDiagram",
[IsHapCommutativeDiagram],
function(T)
Print("Commutative diagram.\n");
end);


#####################################################################
#####################################################################



