DeclareFilter("QuadraticNF");
DeclareProperty("IsQuadraticNumberField",QuadraticNF);
DeclareFilter("RingOfQuadraticIntegers");
DeclareProperty("IsRingOfQuadraticIntegers",RingOfQuadraticIntegers);
DeclareAttribute("AssociatedNumberField",IsNumberField);
DeclareAttribute("AssociatedRing",IsRing);
DeclareAttribute("NormOfIdeal",IsInt);
DeclareGlobalFunction("QuadraticNumberField");
DeclareOperation("RingOfIntegers",[IsNumberField]);
DeclareOperation("QuadraticIdeal",[IsRing and IsRingOfQuadraticIntegers, IsCyclotomic]);
DeclareFilter("IdealOfQuadraticIntegers");
DeclareProperty("IsIdealOfQuadraticIntegers",IdealOfQuadraticIntegers);
DeclareGlobalFunction("PartsOfQuadraticInteger");
DeclareGlobalFunction("SL2QuadraticIntegers");


#####################################################################
#####################################################################
DeclareCategory("IsHapQuadraticNumber",IsObject);

DeclareRepresentation( "IsHapQuadraticNumberRep",
                        IsComponentObjectRep,
                        ["rational",
                         "irrational",
                         "bianchiInteger"]);

HapQuadraticNumberFamily:=NewFamily( "HapQuadraticNumberFamily",
                                   IsHapQuadraticNumber,
                                   IsHapQuadraticNumber);

HapQuadraticNumber:=NewType(HapQuadraticNumberFamily,IsHapQuadraticNumberRep);

InstallMethod( ViewObj,
"for HapQuadraticNumber",
[IsHapQuadraticNumber],
 function(R)
if R!.irrational=0 then Print(R!.rational);
else
    if R!.rational=0 then Print(R!.irrational," Sqrt(",R!.bianchiInteger,")");
    else
    Print(R!.rational, " + ", R!.irrational, " Sqrt(",R!.bianchiInteger,")");
    fi;
fi;
 end);

InstallMethod( PrintObj,
"for HapQuadraticNumber",
[IsHapQuadraticNumber],
 function(R)
if R!.irrational=0 then Print(R!.rational);
else
    if R!.rational=0 then Print(R!.irrational," Sqrt(",R!.bianchiInteger,")");
    else
    Print(R!.rational, " + ", R!.irrational, " Sqrt(",R!.bianchiInteger,")");
    fi;
fi;

 end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareGlobalFunction( "QuadraticNumber" );
DeclareGlobalFunction( "QuadraticNumberConjugate" );


DeclareGlobalFunction("DisplayUnimodularPairs");
DeclareGlobalFunction("Display3DUnimodularPairs");
DeclareGlobalFunction("UnimodularIntersectingLine");
DeclareGlobalFunction("NeighbourhoodOfUnimodularPairs");
DeclareGlobalFunction("HAP_VertexHeights");
DeclareGlobalFunction("SwanBianchiCriterion");
DeclareGlobalFunction("QQNeighbourhoodOfUnimodularPairs");
DeclareGlobalFunction("QNeighbourhoodOfUnimodularPairs");
DeclareGlobalFunction("IsQUnimodularPair");
DeclareGlobalFunction("IsQQUnimodularPair");
DeclareGlobalFunction("HAP_BianchiFundamentalRectangle");
DeclareGlobalFunction("UnimodularPairStandardForm");
DeclareGlobalFunction("QuadraticIntegersByNorm");
DeclareGlobalFunction("HAP_UnimodularComplements");
DeclareGlobalFunction("HAP_SqrtInequality");
DeclareGlobalFunction("HAP_SqrtStrictInequality");
DeclareGlobalFunction("HAP_IsRedundantUnimodularPair");
DeclareGlobalFunction("UnimodularPairCoordinates");
DeclareGlobalFunction("UnimodularPairs");
DeclareGlobalFunction("UnimodularIntersectingPoint");
DeclareGlobalFunction("HAP_HeightOfPointOnSphere");
DeclareGlobalFunction("HAP_AreIntersectingUnimodularPairs");
DeclareGlobalFunction("HAP_AreStrictlyIntersectingUnimodularPairs");
DeclareGlobalFunction("HAP_Are3IntersectingUnimodularPairs");
DeclareGlobalFunction("HAP_PrintFloat");
DeclareGlobalFunction("IsStrictlyFundamentalUnimodularPair");
DeclareGlobalFunction("AreStrictlyFundamentalCoordinates");
DeclareGlobalFunction("SimplicialComplexOfUnimodularPairs");
DeclareGlobalFunction("UnimodularPairsReduced");
DeclareGlobalFunction("BianchiPolyhedron");
DeclareGlobalFunction("CoverOfUnimodularPairs");
DeclareGlobalFunction("IsUnimodularCollection");

####################################################################
#####################################################################
DeclareCategory("IsHapBianchiPolyhedron",IsObject);

DeclareRepresentation( "IsHapBianchiPolyhedronRep",
                        IsComponentObjectRep,
                        ["unimodularPairs",
                         "bianchiInteger",
                         "ring"]);

HapBianchiPolyhedronFamily:=NewFamily( "HapBianchiPolyhedronFamily",
                                   IsHapBianchiPolyhedron,
                                   IsHapBianchiPolyhedron);

HapBianchiPolyhedron:=NewType(HapBianchiPolyhedronFamily,IsHapBianchiPolyhedronRep);

InstallMethod( ViewObj,
"for HapBianchiPolyhedron",
[IsHapBianchiPolyhedron],
 function(P)
    Print("3-dimensional Bianchi polyhedron over OQ( Sqrt(",P!.bianchiInteger,") ) involving hemispheres of minimum squared radius ",P!.minRadius," and non-cuspidal vertices of minimum squared height ",P!.minVertexHeight, " . \n");
 end);

InstallMethod( PrintObj,
"for HapBianchiPolyhedron",
[IsHapBianchiPolyhedron],
 function(P)
    Print("3-dimensional Bianchi polyhedron over OQ( Sqrt(",P!.bianchiInteger,") )\n");
 end);

InstallMethod( Display,
"for HapBianchiPolyhedron",
[IsHapBianchiPolyhedron],
 function(P)
DisplayUnimodularPairs(P!.ring,P!.unimodularPairs,P!.cusps);
 end);

DeclareOperation("Display3D",[IsHapBianchiPolyhedron]);
DeclareOperation("Display2D",[IsHapBianchiPolyhedron]);

InstallMethod( Display2D,
"for HapBianchiPolyhedron",
[IsHapBianchiPolyhedron],
 function(P)
DisplayUnimodularPairs(P!.ring,P!.unimodularPairs,P!.cusps);
 end);

InstallMethod( Display3D,
"for HapBianchiPolyhedron",
[IsHapBianchiPolyhedron],
 function(P)
Display3DUnimodularPairs(P!.ring,P!.unimodularPairs);
 end);
#####################################################################
#####################################################################



