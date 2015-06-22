#(C) GrahamEllis 2005-2006

#####################################################################
InstallOtherMethod(Length,
"length of a resolution",
[IsHapResolution],
function(R) return EvaluateProperty(R,"length");
end);
#####################################################################

#####################################################################
InstallOtherMethod(Length,
"length of a non-free resolution",
[IsHapNonFreeResolution],
function(R) return EvaluateProperty(R,"length");
end);
#####################################################################


#####################################################################
InstallOtherMethod(Source,
"source of chain map",
[IsHapMap],
function(f) return f!.source;
end);
#####################################################################

#####################################################################
InstallMethod(Target,
"target of chain map",
[IsHapMap],
function(f) return f!.target;
end);
#####################################################################

#####################################################################
InstallMethod(Map,
"target of chain map",
[IsHapMap],
function(f) return f!.mapping;
end);
#####################################################################

#####################################################################
InstallOtherMethod(BoundaryMap,
"Boundary function of complex",
[IsHapComplex],
function(R) return R!.boundary;
end);
#####################################################################

#####################################################################
InstallOtherMethod(Dimension,
"rank function of a complex",
[IsHapComplex],
function(f) return f!.dimension;
end);
#####################################################################

#####################################################################
InstallMethod(GroupOfResolution,
"Group underlying a resolution",
[IsHapResolution],
function(R) return R!.group;
end);
#####################################################################

#####################################################################
InstallOtherMethod(Dimension,
"Dimension of FpG module",
[IsHapFPGModule],
function(M) return M!.dimension;
end);
#####################################################################

