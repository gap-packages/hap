#############################################################################
##
#W  goutergroup.gd                   HAP                      Robert F. Morse
##                                                               Graham Ellis
##
##


#############################################################################
##
##  Declare G-OuterGroups and their homomorphisms to be component 
##  objects. Also do StandardNCocycles.
##
DeclareProperty("IsGOuterGroup", IsComponentObjectRep);
DeclareProperty("IsGOuterGroupHomomorphism", IsComponentObjectRep);
DeclareProperty("IsStandardNCocycle", IsComponentObjectRep);



#############################################################################
##
##  Tester operaton  to check if a group homomorphism is mathematically a 
##  G-outer group homomorphisms.
##
DeclareOperation("GOuterHomomorphismTester", 
[IsGOuterGroup,IsGOuterGroup,IsGroupHomomorphism]);


#############################################################################
##
##  Basic attributes of a GOuterGroup
##
DeclareAttribute( "ActingGroup" ,       IsGOuterGroup );
DeclareAttribute( "ActedGroup" ,        IsGOuterGroup );
DeclareAttribute( "OuterAction" ,       IsGOuterGroup );

#############################################################################
##
##  Basic attributes of a GOuterGroup homomorphism
##
DeclareAttribute( "Source" ,        IsGOuterGroupHomomorphism );
DeclareAttribute( "Target" ,        IsGOuterGroupHomomorphism );
DeclareAttribute( "Mapping" ,       IsGOuterGroupHomomorphism );

#############################################################################
##
##  Basic attributes of a standardNcocycle 
##
DeclareAttribute( "CoefficientModule" ,  IsStandardNCocycle );
DeclareAttribute( "Mapping", IsStandardNCocycle);
DeclareAttribute( "Arity", IsStandardNCocycle);


#############################################################################
##
##  Empty Constructors -- each attribute must be set later using 
##  a setter function
##  
##  Example:
##     N := GOuterGroup();
##     SetActingGroup(N,G);
##     SetActedGroup(N,A);
##     SetOuterAction(N,alpha);
## 
DeclareOperation("GOuterGroup", []);
DeclareOperation("GOuterGroupHomomorphism", []);
DeclareOperation("StandardNCocycle",[]);

#############################################################################
##
##  Constructor for G-outer group from abelian group A (module) and
##  group G (assumed to act triviall on A.
##
DeclareOperation("TrivialGModuleAsGOuterGroup", [IsGroup,IsGroup]);

#############################################################################
##
##  Constructor for G-outer group from abelian group A (module) and
##  group G act via a function alpha(g,a) on A.
##
DeclareOperation("GModuleAsGOuterGroup", [IsGroup,IsGroup,IsFunction]);



#############################################################################
##
##  Constructor for G-outer group from group E (Extension) and 
##  normal subgroup A 
##
DeclareOperation("GOuterGroup", [IsGroup,IsGroup]);

#############################################################################
##
##  Constructor for G-outer group from a group E 
##
DeclareOperation("GOuterGroup", [IsGroup]);

#############################################################################
##
##  Constructor for G-outer group homomorphism from a group homomorphism  
##
DeclareOperation("GOuterGroup", [IsGroupHomomorphism]);

#############################################################################
##
##  Constructor for  G-outer group homomorphisms 
##
DeclareOperation("GOuterGroupHomomorphism", 
		[IsGOuterGroup,IsGOuterGroup,IsGroupHomomorphism]);

#############################################################################
##
##  Constructor for a standard N-cocycle 
##
DeclareOperation("StandardNCocycle",
                [IsGOuterGroup,IsFunction,IsInt]);


#############################################################################
##
##  Direct products of G-outer groups 
##
DeclareOperation("DirectProductGog",
                [IsGOuterGroup,IsGOuterGroup]);

DeclareOperation("DirectProductGog",
                [IsList]);


#############################################################################
##
##  Operation for Hom_ZG(R,n,A) and Hom_ZG(R,A) where R is a free ZG-module 
##  resolution, n is an integer and A is an abelian G-outer group. 
##
DeclareOperation("HomToGModule",
		[IsHapResolution,IsInt,IsGOuterGroup,IsGOuterGroup,IsGOuterGroup]);
DeclareOperation("HomToGModule",
                [IsHapResolution,IsGOuterGroup]);


#####################################################################
##
##  Declaration of the G-cocomplex data type.
##
DeclareCategory("IsHapGCocomplex",IsObject);
DeclareRepresentation( "IsHapGCocomplexRep",
                        IsComponentObjectRep,
                        ["boundary",
                         "properties"]);

HapGCocomplexFamily:=NewFamily( "HapGCocomplexFamily",
                                   IsHapGCocomplex,
                                   IsHapGCocomplex);

HapGCocomplex:=NewType(HapGCocomplexFamily,IsHapGCocomplexRep);


#############################################################################
##
##  Operation for returning the cohomology of a G-cochain complex as
##  a G-outer group.
DeclareOperation("CohomologyModule",
                [IsHapGCocomplex,IsInt]);

#####################################################################
#####################################################################
DeclareCategory("IsHapGComplex",IsObject);

DeclareRepresentation(  "IsHapGComplexRep",
                        IsComponentObjectRep,
                        ["boundary",
                         "properties"]);

HapGComplexFamily:=NewFamily("HapGComplexFamily",
                                IsHapGComplex,
                                IsHapGComplex);

HapGComplex:=NewType(HapGComplexFamily,IsHapGComplexRep);

InstallMethod( ViewObj,
"for HapGComplex",
 [IsHapGComplex],
 function(R)
 Print("G-Complex of length ", EvaluateProperty(R,"length"), "\n");
 end);

InstallMethod( PrintObj,
"for HapGComplex",
 [IsHapGComplex],
 function(R)
 Print("G-Complex of length ", EvaluateProperty(R,"length"), "\n");
 end);

#####################################################################
#####################################################################


DeclareOperation("GDerivedSubgroup",[IsGOuterGroup]);
DeclareOperation("LowerGCentralSeries",[IsGOuterGroup]);
DeclareGlobalFunction("AbelianGOuterGroupToCatOneGroup");
DeclareGlobalFunction("ImageOfGOuterGroupHomomorphism");
DeclareGlobalFunction("KernelOfGOuterGroupHomomorphism");
DeclareOperation("CohomologyClass",[IsGOuterGroup,IsStandardNCocycle]);
DeclareGlobalFunction("CocyclicMatrices");
DeclareGlobalFunction("IsHadamardMatrix");
DeclareGlobalFunction("HadamardGraph");
DeclareGlobalFunction("CocyclicHadamardMatrices");

