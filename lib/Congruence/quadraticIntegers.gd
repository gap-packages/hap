DeclareFilter("QuadraticNF");
DeclareProperty("IsQuadraticNumberField",QuadraticNF);
DeclareFilter("RingOfQuadraticIntegers");
DeclareProperty("IsRingOfQuadraticIntegers",RingOfQuadraticIntegers);
DeclareAttribute("AssociatedNumberField",IsNumberField);
DeclareAttribute("AssociatedRing",IsRing);
DeclareAttribute("NormOfIdeal",IsInt);
DeclareGlobalFunction("QuadraticNumberField");
DeclareGlobalFunction("HAPQuadraticRing");
DeclareOperation("RingOfIntegers",[IsNumberField]);
DeclareOperation("QuadraticIdeal",[IsRing and IsRingOfQuadraticIntegers, IsCyclotomic]);
DeclareFilter("IdealOfQuadraticIntegers");
DeclareProperty("IsIdealOfQuadraticIntegers",IdealOfQuadraticIntegers);
DeclareGlobalFunction("PartsOfQuadraticInteger");
DeclareGlobalFunction("SL2QuadraticIntegers");
DeclareGlobalFunction( "QuadraticNumber" );
DeclareGlobalFunction( "IsHapQuadraticInteger" );

#####################################################################
#####################################################################
DeclareCategory("IsHapQuadraticNumber",IsScalar);
cat:=CategoryCollections(IsHapQuadraticNumber);
cat:=CategoryCollections(cat);
#cat:=CategoryCollections(cat);
#cat:=CategoryCollections(cat);

DeclareRepresentation( "IsHapQuadraticNumberRep",
                        IsComponentObjectRep,
                        ["rational",
                         "irrational",
                         "bianchiInteger"]);

HapQuadraticNumberFamily:=NewFamily( "HapQuadraticNumberFamily",
                                   #IsHapQuadraticNumber,
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

Hap2x2:=CategoryCollections(CategoryCollections(IsHapQuadraticNumber));
BindGlobal("IsHap2x2matrix",Hap2x2);
Unbind(Hap2x2);
Hap2x2:=CategoryCollections(CategoryCollections(CategoryCollections(IsHapQuadraticNumber)));
Hap2x2:=Hap2x2 and IsGroup;
BindGlobal("IsHap2x2matrixGroup",Hap2x2);
Hap2x2:=Hap2x2 and IsAbelian;
BindGlobal("IsBianchiAbelianGroup",Hap2x2);
Unbind(Hap2x2);


InstallGlobalFunction( HAPQuadraticRing, function( d )
    local F, R;

    if not IsInt( d ) then
      Error( "<d> must be an  ideal" );
    fi;

    # Construct the family of element objects of our ring.
    F:= NewFamily( "HAPQuadraticRing" ,
                   IsHapQuadraticNumber );

    # Install the data.
    F!.bianchiNumber:= d;

    # Make the domain.
R:= RingWithOneByGenerators(  [ QuadraticNumber( 0,1,d) ]  );
    SetIsWholeFamily( R, true );
    SetName( R,  Concatenation("Q(Sqrt( ", String(d), ")"  ));
    R!.bianchiInteger:=d;
    SetCharacteristic(R,0);
    SetSize(R,infinity);
    # Return the ring.
    return R;
 end );


#####################################################################
#####################################################################
DeclareGlobalFunction( "QuadraticNumberConjugate" );


DeclareGlobalFunction("DisplayUnimodularPairs");
DeclareGlobalFunction("Display3DUnimodularPairs");
DeclareGlobalFunction("UnimodularIntersectingLine");
DeclareGlobalFunction("NeighbourhoodOfUnimodularPairs");
DeclareGlobalFunction("HAP_VertexHeights");
DeclareGlobalFunction("SwanBianchiCriterion");
DeclareGlobalFunction("SwanBianchiCriterion_alt");
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
DeclareGlobalFunction("UnimodularPairsReduced_NN");
DeclareGlobalFunction("BianchiPolyhedron");
DeclareGlobalFunction("CoverOfUnimodularPairs");
DeclareGlobalFunction("IsUnimodularCollection");
DeclareGlobalFunction("HAP_BianchiRegularCWComplex");
DeclareGlobalFunction("HAP_BianchiTransformations");
DeclareGlobalFunction("Hap_int");

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



