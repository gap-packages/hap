DeclareOperation("TensorWithGModule",[IsHapResolution,IsInt,IsGOuterGroup,IsGOuterGroup,IsGOuterGroup]);

DeclareOperation("HomologyModule",[IsHapGComplex,IsInt]);

#####################################################################
#####################################################################
##
InstallMethod( ViewObj,
"for HapGComplex",
[IsHapGComplex],
 function(R)
Print("G-complex of length ",
EvaluateProperty(R,"length"), " . \n");
 end);

InstallMethod( PrintObj,
"for HapGComplex",
[IsHapGComplex],
function(R)
Print("G-complex of length ",
EvaluateProperty(R,"length"), " . \n");
end);
#######################################################################
#######################################################################


