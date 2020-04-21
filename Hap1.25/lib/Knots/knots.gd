#####################################################################
#####################################################################
DeclareCategory("IsHapPureCubicalLink",IsObject);

DeclareRepresentation(  "IsHapPureCubicalLinkRep",
                        IsComponentObjectRep and
                        IsHapPureCubicalComplex,
                        ["binaryArray",
                         "name",
                         "properties"]);

HapPureCubicalLinkFamily:=NewFamily( "HapPureCubicalLinkFamily",
                                 IsHapPureCubicalLink,
                                 IsHapPureCubicalLink);

HapPureCubicalLink:=NewType(HapPureCubicalLinkFamily,IsHapPureCubicalLinkRep);


InstallMethod( ViewObj,
"for HapPureCubicalLink",
[IsHapPureCubicalLink],
function(T)
if EvaluateProperty(T,"knot") then
Print(T!.name,"\n");
else
Print("Pure cubical link.\n");
fi;
end);

InstallMethod( PrintObj,
"for HapPureCubicalLink",
[IsHapPureCubicalLink],
function(T)
if EvaluateProperty(T,"knot") then
Print(T!.name,"\n");
else
Print("Pure cubical link.\n");
fi;
end);
#####################################################################
#####################################################################

