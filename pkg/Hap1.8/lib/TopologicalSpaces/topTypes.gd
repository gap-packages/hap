#####################################################################
#####################################################################
DeclareCategory("IsHapTopologicalSpace",IsObject);

DeclareRepresentation(  "IsHapTopologicalSpaceRep",
                        IsComponentObjectRep,
                        ["rawData",
                         "properties"]);

HapTopologicalSpaceFamily:=NewFamily( "HapTopologicalSpaceFamily",
                                 IsHapTopologicalSpace,
                                 IsHapTopologicalSpace);

HapTopologicalSpace:=NewType(HapTopologicalSpaceFamily,IsHapTopologicalSpaceRep);


InstallMethod( ViewObj,
"for HapTopologicalSpace",
[IsHapTopologicalSpace],
function(T)
Print("Topological space of dimension ", 
EvaluateProperty(T,"dimension"), ".\n");
end);

InstallMethod( PrintObj,
"for HapTopologicalSpace",
[IsHapTopologicalSpace],
function(T)
Print("Topological space of dimension ",
EvaluateProperty(T,"dimension"), ".\n");
end);


#####################################################################
#####################################################################

