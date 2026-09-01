
#############################################################################
##
## A congruence subgroup G is viewed as a pair of groups G < GG where GG
## is some matrix group such as SL(n,R), GL(n,R), Sp(n,R) and G < GG is a 
## finite index subgroup satisfying the "congruence subgroup" property. We'll
## refer to GG as the "ambient group" and say that G is defined over the ring 
## R. The implementation allows GG to itself be some congruence subgroup.
##
## This file is copied/modified from the GAP package Congruence authored by 
## Ann Dooms, Eric Jespers, Olexandr Konovalov, Helena Verrill.

#############################################################################
##
## Constructors for congruence subgroups 
DeclareOperation("PrincipalCongruenceSubgroup", [IsMatrixGroup,IsRingElement]);
DeclareOperation("CongruenceSubgroupGamma0", [IsMatrixGroup,IsRingElement]);
DeclareOperation("CongruenceSubgroupGamma1", [IsMatrixGroup,IsRingElement]);

#############################################################################
##
## IsHAPCongruenceSubgroup( <G> )
## 
## We create the category of congruence subgroups as a subcategory of matrix 
## groups, and declare properties that are used to distinguish several
## important classes of congruence subgroups
DeclareCategory( "IsHAPCongruenceSubgroup", IsMatrixGroup );

#############################################################################
##
## IsHAPConjugatedCongruenceSubgroup( <G> )
##
## We create a category of congruence subgroups conjugated by a matrix
DeclareCategory("IsHAPConjugatedCongruenceSubgroup", IsHAPCongruenceSubgroup);


#############################################################################
##
## AmbientGroupOfCongruenceSubgroup( <G> )
##
## The ambient group GG of a congruence subgroup G is the group in which G is
## considered to be a finite index subgroup.
##
DeclareAttribute( "AmbientGroupOfCongruenceSubgroup", IsHAPCongruenceSubgroup );

#############################################################################
##
## LevelOfCongruenceSubgroup( <G> )
##
## The (arithmetic) level of a congruence subgroup G is the smallest ideal 
## N such that G contains the principal congruence subgroup of level N. We
## allow N to be an ideal in the ring R or an element in the ring (generating 
## a principal ideal in R). 
##
DeclareAttribute( "LevelOfCongruenceSubgroup", IsHAPCongruenceSubgroup );

#############################################################################
##
## IsHAPPrincipalCongruenceSubgroup( <G> )
## 
## The principal congruence subgroup of level N consists of all matrices
## (in the 2x2 case) of the form   [ 1+N    N ]
##                                 [   N  1+N ]
##
DeclareProperty( "IsHAPPrincipalCongruenceSubgroup", IsHAPCongruenceSubgroup );
InstallTrueMethod(IsHAPCongruenceSubgroup, IsHAPPrincipalCongruenceSubgroup);


#############################################################################
##
## IsHAPCongruenceSubgroupGamma0( <G> )
## 
## The congruence subgroup CongruenceSubgroupGamma0(N) consists of all matrices
## (in the 2x2 case) of the form   [   *    * ]
##                                 [   N    * ]
##
DeclareProperty( "IsHAPCongruenceSubgroupGamma0", IsHAPCongruenceSubgroup );
InstallTrueMethod(IsHAPCongruenceSubgroup, IsHAPCongruenceSubgroupGamma0);


#############################################################################
##
## IsCongruenceSubgroupGamma1( <G> )
## 
## The congruence subgroup CongruenceSubgroupGamma1(N) consists of all matrices
## of the form   [ 1+N    * ]
##               [   N  1+N ]
##
DeclareProperty( "IsHAPCongruenceSubgroupGamma1", IsHAPCongruenceSubgroup );
InstallTrueMethod(IsHAPCongruenceSubgroup, IsHAPCongruenceSubgroupGamma1);

#############################################################################
##
## IndexInAmbientGroup( <G> )
##
## The index of a congruence subgroup in SL_2(Z) will be stored as its 
## attribute. This also will allow us to install a method for Index(G,H) when
## G is SL_2(Z) and H is a congruence subgroup. You should remember that we
## are working with the SL_2(Z), because it is available in GAP, and not with
## the PSL_2(Z) since the latter is not implemented in GAP.
##
DeclareAttribute( "IndexInAmbientGroup", IsHAPCongruenceSubgroup );

############################################################################
##
## StabilizerSubgroup( <G> )
##
## It is sometimes convenient to fix a finite subgroup U in the ambient
## group GG of G and then make calculations using the coset graph of U with
## respect to  generators of GG. We can always choose U to be the identity 
## group. But for, say, the principal congruence subgroup of SL(2,U) it is
## handy to choose U equal to a group of order 6 so that the coset graph of U
## is the cubic tree.
DeclareAttribute( "StabilizerSubgroup", IsHAPCongruenceSubgroup );

############################################################################
##
## ProjectiveSpave( <G> )
##
DeclareAttribute( "ProjectiveSpace", IsHAPCongruenceSubgroup );

############################################################################
##
## AmbientTransversal( <G> )
##
## returns a right transversal for the congruence subgroup G in the ambient 
## group GG.
DeclareAttribute( "AmbientTransversal", IsHAPCongruenceSubgroup );

############################################################################
##
## AmbientPosition( <G> )
##
## returns a function cosetPos(g) that inputs an element g in the ambient 
## group GG and returns an integer corresponding to the position of the left
## coset gG in the lists of cosets of G in GG. 
DeclareAttribute( "AmbientPosition", IsHAPCongruenceSubgroup );

############################################################################
##
## AmbientRepresentation( <G> )
##
## returns a function cosetRep(g) that inputs an element g in the ambient 
## group GG and returns a canonical representative in the ambient group GG
## of the left coset gG. 
DeclareAttribute( "AmbientRepresentation", IsHAPCongruenceSubgroup );

############################################################################
##
## AmbientTree( <G> )
##
## returns a tree representing the left 'cosets' of G in the ambient group GG 
## whose construction depends on the given generating set for GG. 
DeclareAttribute( "AmbientTree", IsHAPCongruenceSubgroup );

############################################################################
##
## AmbientTreeDisplay( <G> )
##
## displays the tree associated to the congruence subgroup G
DeclareGlobalFunction("AmbientTreeDisplay");

############################################################################
##
## IntersectionWithConjugate( <G> , m )
##
## returns the intersection of a group G with its conjugate G^m
DeclareOperation( "IntersectionWithConjugate", [ IsHAPCongruenceSubgroup, IsMatrix ]);
############################################################################
##
## 
     IsHAPRightTransversalCongruenceSubgroup:=NewFilter("IsHAPRightTransversalCongruenceSubgroup");;

     HAPRightTransversalCongruenceSubgroup:=NewType(FamilyObj(rec()),
                             IsHAPRightTransversalCongruenceSubgroup
                             and IsList
                              );;

     InstallMethod( ViewObj,
     "for HAPRightTransversalCongruenceSubgroup",
     [IsHAPRightTransversalCongruenceSubgroup],
     function(R)
     Print("Transversal of ");
     ViewObj(R!.subgroup); Print(" in "); ViewObj(R!.group);Print(" . \n");
     end);

     InstallMethod( PrintObj,
     "for HAPRightTransversalCongruenceSubgroup",
     [IsHAPRightTransversalCongruenceSubgroup],
     function(R)
     Print("Transversal of ");
     ViewObj(R!.subgroup); Print(" in "); ViewObj(R!.group);Print(" . \n");
     end);

DeclareGlobalFunction("VoronoiGenerators");
DeclareGlobalFunction("VoronoiBaseMatrix");
DeclareGlobalFunction("MatrixFromWord");
DeclareGlobalFunction("VoronoiReduce");
DeclareGlobalFunction("VoronoiFactorize");
DeclareGlobalFunction("VoronoiCoordinates");
DeclareGlobalFunction("VoronoiWord");
DeclareGlobalFunction("VoronoiBuildStabilizer");

