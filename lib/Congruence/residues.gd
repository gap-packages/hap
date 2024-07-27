
DeclareCategory( "IsHAPRingModIdealObj", IsScalar );
cat:= CategoryCollections( IsHAPRingModIdealObj );;
cat:= CategoryCollections( cat );;
cat:= CategoryCollections( cat );;

DeclareRepresentation("IsHAPIdealRep", IsPositionalObjectRep, [1]);

DeclareGlobalFunction( "HAPRingModIdealObj" );

InstallGlobalFunction( HAPRingModIdealObj, function( Fam, residue )
    return Objectify( NewType( Fam, IsHAPRingModIdealObj and IsHAPIdealRep ),
                      [ residue mod Fam!.modulus ] );
 end );

DeclareGlobalFunction( "HAPRingModIdeal" );
 

InstallGlobalFunction( HAPRingModIdeal, function( I )
    local F, R;
 
    if not IsIdealOfQuadraticIntegers( I ) then
      Error( "<I> must be an  ideal" );
    fi;
 
    # Construct the family of element objects of our ring.
    F:= NewFamily( "HAPRingModIdeal" ,
                   IsHAPRingModIdealObj );
 
    # Install the data.
    F!.modulus:= I;
 
    # Make the domain.
R:= RingWithOneByGenerators(  List(GeneratorsOfRing(AssociatedRing(I)) , x -> HAPRingModIdealObj( F, x) )  );
    SetIsWholeFamily( R, true );
    SetName( R,  Concatenation("ring mod ideal of norm ", String(Norm(I)))  );
    R!.associatedIdeal:=I; 
    # Return the ring.
    return R;
 end );


InstallMethod( PrintObj,
    "for element in Z/nZ (ModulusRep)",
    [ IsHAPRingModIdealObj and IsHAPIdealRep ],
    function( x )
    Print(  x![1], " mod ideal of norm ", String(Norm(FamilyObj(x)!.modulus))  );
    end );
InstallMethod( \=,
    "for two elements in Z/nZ (ModulusRep)",
    IsIdenticalObj,
    [IsHAPRingModIdealObj and IsHAPIdealRep, IsHAPRingModIdealObj and IsHAPIdealRep],
    function( x, y ) return x![1]  = y![1] ; 
end );
 
InstallMethod( \<,
    "for two elements in Z/nZ (ModulusRep)",
    IsIdenticalObj,
    [IsHAPRingModIdealObj and IsHAPIdealRep, IsHAPRingModIdealObj and IsHAPIdealRep],
    function( x, y ) return x![1]  < y![1] ;
 end );

InstallMethod( \+,
    "for two elements in Z/nZ (ModulusRep)",
    IsIdenticalObj,
    [IsHAPRingModIdealObj and IsHAPIdealRep, IsHAPRingModIdealObj and IsHAPIdealRep],
    function( x, y )
    return HAPRingModIdealObj( FamilyObj( x ), x![1] + y![1]  ) ;
    end );
 
InstallMethod( ZeroOp,
    "for element in Z/nZ (ModulusRep)",
    [ IsHAPRingModIdealObj ],
    x -> HAPRingModIdealObj( FamilyObj( x ), 0)  );
 
InstallMethod( AdditiveInverseOp,
    "for element in Z/nZ (ModulusRep)",
    [ IsHAPRingModIdealObj and IsHAPIdealRep ],
    x -> HAPRingModIdealObj( FamilyObj( x ), 
AdditiveInverse( x![1] ) ) );



InstallMethod( \+,
    "for element in Z/nZ (ModulusRep) and integer",
    [ IsHAPRingModIdealObj and IsHAPIdealRep, IsCyclotomic ],
    function( x, y )
    return HAPRingModIdealObj( FamilyObj( x ), 
    (x![1] + y)  );
    end );

InstallMethod( \^,
    "for element in Z/nZ (ModulusRep) and integer",
    [ IsHAPRingModIdealObj and IsHAPIdealRep, IsInt ],
    function( x, n )
    if n>=0 then
    return HAPRingModIdealObj( FamilyObj( x ), 
    ( (x![1])^n )  );
    else
    return HAPRingModIdealObj( FamilyObj( x ), (InverseOp(FamilyObj(x)!.modulus,x![1]))^AbsInt(n) );
    fi; 
    end );

 
InstallMethod( \+,
    "for integer and element in Z/nZ (ModulusRep)",
    [ IsCyclotomic, IsHAPRingModIdealObj and IsHAPIdealRep ],
    function( x, y )
    return HAPRingModIdealObj( FamilyObj( y ), 
(x + y![1])  );
    end );

InstallMethod( \*,
    "for two elements in Z/nZ (ModulusRep)",
    IsIdenticalObj,
    [IsHAPRingModIdealObj and IsHAPIdealRep, IsHAPRingModIdealObj and IsHAPIdealRep],
    function( x, y )
#if not IsIntegralCyclotomic(x![1]) then return fail; fi;
    return HAPRingModIdealObj( FamilyObj( x ), x![1] * y![1] );
    end );

InstallMethod( \*,
    "for integer and element in Z/nZ (ModulusRep)",
    [IsCyclotomic, IsHAPRingModIdealObj and IsHAPIdealRep],
    function( x, y )
#    if not IsIntegralCyclotomic(x) then return fail; fi;
    return HAPRingModIdealObj( FamilyObj( y ), x * y![1] );
    end );

InstallMethod( \*,
    "for integer and element in Z/nZ (ModulusRep)",
    [IsHAPRingModIdealObj and IsHAPIdealRep,IsCyclotomic],
    function( x, y )
    return HAPRingModIdealObj( FamilyObj( x ), x![1] * y );
    end );


 
InstallMethod( OneOp,
    "for element in Z/nZ (ModulusRep)",
    [ IsHAPRingModIdealObj ],
    elm -> HAPRingModIdealObj( FamilyObj( elm ), 1 ) );
 
InstallMethod( InverseOp,
    "for element in Z/nZ (ModulusRep)",
    [ IsHAPRingModIdealObj and IsHAPIdealRep ],
    function( elm )
    local residue;
    #residue:=   elm![1]^-1 mod FamilyObj( elm )!.modulus ;
    residue:= InverseOp(FamilyObj( elm )!.modulus,elm![1]);
    if residue <> fail then
      residue:= HAPRingModIdealObj( FamilyObj( elm ), residue );
    fi;
    return residue;
    end );

InstallMethod( Int,
    "for element in Z/nZ (ModulusRep)",
    [ IsHAPRingModIdealObj and IsHAPIdealRep ],
    z -> z![1] );

InstallMethod( PrintObj,
    "for full collection Z/nZ",
    [ CategoryCollections( IsHAPRingModIdealObj ) and IsWholeFamily ],
    function( R )
    Print( "(Integers mod ",
           ElementsFamily( FamilyObj(R) )!.modulus, ")" );
    end );
 
InstallMethod( \mod,
    "for `quadratic Integers', and an ideal",
    [ IsRingOfQuadraticIntegers, IsIdealOfQuadraticIntegers ],
    function( Integers, I ) return HAPRingModIdeal( I ); end );

InstallMethod( InverseOp,
    "for an ordinary matrix over a ring Z/nZ",
    [ IsMatrix and IsOrdinaryMatrix
      and CategoryCollections( CategoryCollections( IsHAPRingModIdealObj ) ) ],
    function( mat )
    local one, modulus;
    one:= One( mat[1][1] );
    if not Determinant(mat)=one then return fail; fi;
    modulus:= FamilyObj( one )!.modulus;
    mat:=  List( mat, row -> List( row, Int ) ) ;
    mat:= Determinant(mat)*InverseOp(mat);              

    if mat <> fail then
      #mat:= ( mat mod modulus ) * one;
      mat:= List(mat, row ->List(row, r -> HAPRingModIdealObj(FamilyObj( one ),r))); 
    fi;
    if not IsMatrix( mat ) then
      mat:= fail;
    fi;
    return mat;
    end );

InstallMethod( InverseSameMutability,
    "for an ordinary matrix over a ring Z/nZ",
    [ IsMatrix and IsOrdinaryMatrix
      and CategoryCollections( CategoryCollections( IsHAPRingModIdealObj ) ) ],
    function( mat );
    return InverseOp(mat);
    end );

InstallMethod( \^,
    "for two ordinary matrices over a ring Z/nZ",
    [ IsMatrix and IsOrdinaryMatrix
      and CategoryCollections( CategoryCollections( IsHAPRingModIdealObj ) ),
      IsMatrix and IsOrdinaryMatrix
      and CategoryCollections( CategoryCollections( IsHAPRingModIdealObj ) ) ],
    function( mat1 , mat2 )
    local one, modulus;
    return  InverseOp(mat2)*mat1*mat2;
    end );

InstallMethod( \^,
    "for a  matrix over a ring Z/nZ and an integer",
    [ IsMatrix and IsOrdinaryMatrix
      and CategoryCollections( CategoryCollections( IsHAPRingModIdealObj ) ),
      IsInt ],
    function( mat , n )
    local one;
    one:= One( mat[1][1] );
    mat:=mat^n;
    mat:= List( mat, row -> List( row, Int ) ) ;
    mat:= List(mat, row ->List(row, r -> HAPRingModIdealObj(FamilyObj( one ),r)));
   return mat;
    end );



InstallMethod( DefaultFieldOfMatrixGroup,
     "for a matrix group over a ring Z/nZ",
     [ IsMatrixGroup and CategoryCollections( CategoryCollections(
           CategoryCollections( IsHAPRingModIdealObj ) ) ) ],
     G -> RingWithOneByGenerators([ One( Representative( G )[1][1] ) ]));

InstallMethod( Enumerator,
    "for full collection Z/nZ",
    [ CategoryCollections( IsHAPRingModIdealObj ) and IsWholeFamily ],
    function( R )
    local F;
    F:= ElementsFamily( FamilyObj(R) );
    return List( [ 0 .. Size( R ) - 1 ], x -> HAPRingModIdealObj( F, x ) );
    end );

InstallMethod( Random,
    "for full collection Z/nZ",
    [ CategoryCollections( IsHAPRingModIdealObj ) and IsWholeFamily ],
    R -> HAPRingModIdealObj( ElementsFamily( FamilyObj(R) ),
                    Random( 0, Size( R ) - 1 ) ) );
 
InstallMethod( Size,
    "for full ring Z/nZ",
    [ CategoryCollections( IsHAPRingModIdealObj ) and IsWholeFamily ],
    #R -> ElementsFamily( FamilyObj(R) )!.modulus );
    R -> Size( RightTransversal ( R!.associatedIdeal )) );
 
InstallMethod( Units,
    "for full ring Z/nZ",
    [     CategoryCollections( IsHAPRingModIdealObj )
      and IsWholeFamily and IsRing ],
    function( R )
    local F;
    F:= ElementsFamily( FamilyObj( R ) );
    return List( PrimeResidues( Size(R) ), x -> HAPRingModIdealObj( F, x ) );
    end );

InstallMethod( Units,
    "for full ring Z/nZ",
    [     CategoryCollections( IsHAPRingModIdealObj )
      and IsWholeFamily and IsRing ],
    function( R )
    local G, gens;
 
    gens:= GeneratorsPrimeResidues( Size( R ) ).generators;
    if not IsEmpty( gens ) and gens[ 1 ] = 1 then
      gens:= gens{ [ 2 .. Length( gens ) ] };
    fi;
    gens:= Flat( gens ) * One( R );
    return GroupByGenerators( gens, One( R ) );
    end );

InstallTrueMethod( IsFinite,
    CategoryCollections( IsHAPRingModIdealObj ) and IsDomain );



InstallOtherMethod(InverseSameMutability,
    "for an element in Z/nZ (ModulusRep)",
#    IsIdenticalObj,
    [IsHAPRingModIdealObj and IsHAPIdealRep],
    function( x )
#    return HAPRingModIdealObj( FamilyObj( x ), x![1]^-1  );
    return HAPRingModIdealObj( FamilyObj( x ), InverseOp(FamilyObj(x)!.modulus,x![1])  );
    end );

