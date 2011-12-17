#####################################################################
#####################################################################
DeclareCategory("IsHapTopologicalManifold",IsObject);

DeclareRepresentation(  "IsHapTopologicalManifoldRep",
                        IsComponentObjectRep,
                        ["BinaryList",
                         "properties"]);

HapTopologicalManifoldFamily:=NewFamily( "HapTopologicalManifoldFamily",
                                 IsHapTopologicalManifold,
                                 IsHapTopologicalManifold);

HapTopologicalManifold:=NewType(HapTopologicalManifoldFamily,IsHapTopologicalManifoldRep);


InstallMethod( ViewObj,
"for HapTopologicalManifold",
[IsHapTopologicalManifold],
function(T)
Print("Topological manifold of dimension ", 
EvaluateProperty(T,"dimension"), ".\n");
end);

InstallMethod( PrintObj,
"for HapTopologicalManifold",
[IsHapTopologicalManifold],
function(T)
Print("Topological Manifold of dimension ",
EvaluateProperty(T,"dimension"), ".\n");
end);


#####################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareCategory("IsHapTopologicalSpace",IsObject);

DeclareRepresentation(  "IsHapTopologicalSpaceRep",
                        IsComponentObjectRep,
                        ["BinaryList",
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

