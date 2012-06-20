#####################################################################
#####################################################################
DeclareCategory("IsHapPureCubicalComplex",IsObject);

DeclareRepresentation(  "IsHapPureCubicalComplexRep",
                        IsComponentObjectRep,
                        ["binaryArray",
                         "properties"]);

HapPureCubicalComplexFamily:=NewFamily( "HapPureCubicalComplexFamily",
                                 IsHapPureCubicalComplex,
                                 IsHapPureCubicalComplex);

HapPureCubicalComplex:=NewType(HapPureCubicalComplexFamily,IsHapPureCubicalComplexRep);


InstallMethod( ViewObj,
"for HapPureCubicalComplex",
[IsHapPureCubicalComplex],
function(T)
Print("Pure cubical complex of dimension ", 
EvaluateProperty(T,"dimension"), ".\n");
end);

InstallMethod( PrintObj,
"for HapPureCubicalComplex",
[IsHapPureCubicalComplex],
function(T)
Print("Pure cubical complex of dimension ",
EvaluateProperty(T,"dimension"), ".\n");
end);


#####################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareCategory("IsHapCubicalComplex",IsObject);

DeclareRepresentation(  "IsHapCubicalComplexRep",
                        IsComponentObjectRep,
                        ["binaryArray",
                         "properties"]);

HapCubicalComplexFamily:=NewFamily( "HapCubicalComplexFamily",
                                 IsHapCubicalComplex,
                                 IsHapCubicalComplex);

HapCubicalComplex:=NewType(HapCubicalComplexFamily,IsHapCubicalComplexRep);


InstallMethod( ViewObj,
"for HapCubicalComplex",
[IsHapCubicalComplex],
function(T)
if EvaluateProperty(T,"nonregular")=true then
  Print("Non-regular cubical complex of dimension ",
  EvaluateProperty(T,"dimension"));
else
  Print("Cubical complex of dimension ",
  EvaluateProperty(T,"dimension"));
fi;
if IsBound(T!.vectors) then 
  Print(" with discrete vector field.\n");
else
  Print(".\n");
fi;
end);

InstallMethod( PrintObj,
"for HapCubicalComplex",
[IsHapCubicalComplex],
function(T)
if EvaluateProperty(T,"nonregular")=true then
  Print("Non-regular cubical complex of dimension ",
  EvaluateProperty(T,"dimension"));
else
  Print("Cubical complex of dimension ",
  EvaluateProperty(T,"dimension"));
fi;
if IsBound(T!.vectors) then
  Print(" with discrete vector field.\n");
else
  Print(".\n");
fi;
end);

#####################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareCategory("IsHapSimplicialComplex",IsObject);

DeclareRepresentation(  "IsHapSimplicialComplexRep",
                        IsComponentObjectRep,
                        ["vertices",
			 "simplices",
			 "nrSimplices",
			 "enumeratedSimplex",
                         "properties"]);

HapSimplicialComplexFamily:=NewFamily( "HapSimplicialComplexFamily",
                                 IsHapSimplicialComplex,
                                 IsHapSimplicialComplex);

HapSimplicialComplex:=NewType(HapSimplicialComplexFamily,IsHapSimplicialComplexRep);


InstallMethod( ViewObj,
"for HapSimplicialComplex",
[IsHapSimplicialComplex],
function(T)
Print("Simplicial complex of dimension ",
EvaluateProperty(T,"dimension"), ".\n");
end);

InstallMethod( PrintObj,
"for HapSimplicialComplex",
[IsHapSimplicialComplex],
function(T)
Print("Simplicial complex of dimension ",
EvaluateProperty(T,"dimension"), ".\n");
end);


#####################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareCategory("IsHapSimplicialMap",IsObject);

DeclareRepresentation(  "IsHapSimplicialMapRep",
                        IsComponentObjectRep
                        and IsHapMap,
                        ["source",
                        "target",
                        "mapping",
                        "properties"]);

HapSimplicialMapFamily:=NewFamily( "HapSimplicialMapFamily",
                               IsHapSimplicialMap,
                               IsHapSimplicialMap);

HapSimplicialMap:=NewType(HapSimplicialMapFamily, IsHapSimplicialMapRep);


InstallMethod( ViewObj,
"for HapSimplicialMap",
 [IsHapSimplicialMap],
 function(R)
  Print("Simplicial Map\n  ");
 end);


InstallMethod( PrintObj,
"for HapSimplicialMap",
 [IsHapSimplicialMap],
 function(R)
  Print("Simplicial Map\n  ");
 end);

#####################################################################
#####################################################################

