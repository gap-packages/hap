#############################################################################
##
#W  categories.gd                   HAP                          Graham Ellis 
##
##
##  2008-06-27 
##

#############################################################################
##
##  1. Motivation.
##
##  Category Theory (as introduced by Eilenberg and Mac Lane) is a 
##  ubiquitous and well-tested framework for relating different algebraic 
##  theories. Some categories (such as Vector Spaces or ZG-Modules or 
##  Abelian groups) are abelian and thus admit all of the general 
##  constructions of homological algebra (such as chain complexs, mapping 
##  cones, spectral sequences and so forth). Other categories (such as 
##  Groups) are not abelian but do admit some important categorical 
##  properties (such as the existence of an initial object and  
##  a terminal object). In a given category some object may have special
##  properties (such as being projective); some arrows may also have 
##  special properties (such as being normal or being an epimorphism). 
##  These special properties are relative to the category (e.g. a ZG-module
##  may be projective in the category of Abelian Groups but not in the 
##  category of ZG-Modules). Different algebraic theories can be related via 
##  Functors (such as the forgetful functor from ZG-Modules to Abelian Groups,
##  or the functorial equivalence between Crossed Modules and Cat-1-Groups).
##  Functors can also be used to describe certain algebraic constructions 
##  (such as group cohomology). Functors are themselves arrows and can 
##  have properties (such as that of being an equivalence or an isomorphism). 
##  Some mathematical properties are actually properties of functors. (For 
##  example, "freeness" is a property of an adjoint pair of functors.)
##
##  Although (small) categories are interesting algebraic structures in 
##  their own right, this document is only interested in Category Theory as a
##  language for moving objects between algebraic structures and as a language 
##  for implementing general algorithms/constructions in such a way that 
##  they can be applied in a range of different algebraic situations. 
##
##  We need to distinguish between algebraic properties (such as "the group A
##  is an abelian group") and categorical properties (such as "the group A is 
##  an object in the category of abelian groups"). If we have a method 
##  which works for all permutation groups, then it also works for all 
##  abelian permutation groups. However, if we have an implementation of a 
##  construction which works for all permutation groups in the category of 
##  groups, it does not follow that it works correctly for permutation groups 
##  in the category of abelian groups. For instance, coproducts are different
##  in these two categories. 
##  
##  We are aiming at a hierarchical tree in which any higher-level function 
##  can be applied (indirectly) to all entities below it in the tree. A 
##  piece of the tree is as follows.
##
##				Homological diagrams 
##				in abelian categories
##				    |
##				    |
##				Abelian categories
##				    |
##				    |
##                -------------	Categories ------------
##	         /		 /	\	       \	
##		/		/	 \		\
##	     Groups	FpGModules	ZG-modules	MtxModules	
##	     /  |   \
##          /   |    \
##    Perms   Pc-     Fp-
##    groups  groups  groups	
##
##
##  2. Implementation.
##
##  We introduce Eilenberg and Mac Lane's notion into GAP using the 
##  following table. We choose our terminology so that it is consistent with 
##  the standard mathematical notion of "category" yet does not clash with 
##  GAP's existing use of the word "category".
##
##  #####################################################################
##  #####################################################################
##  ##  Eilenberg-Mac Lane	##  Our GAP implememntation 
##  ##  concept			##
##  #####################################################################
##  #####################################################################
##  ##				##
##  ##  Name of a category	##  Immutable global variable
##  ##				##  such as Category_Of_Groups
##  ##				##
##  #####################################################################
##  ##				##
##  ##  Properties of a		##  Property of the category name
##  ##  category	 	##  such as IsAbelian(Category_Of_Groups)
##  ##				##
##  #####################################################################
##  ##  			##
##  ##  An object X in the 	##  Component object X with attributes:
##  ##  category C.     	##  * CategoryName(X) 
##  ##				##  * Object(X) has a value such
##  ##				##    as a group or an FpGModule.
##  ##				##
##  #####################################################################
##  ##				##
##  ##  For each object X in	##  Operation IdentityArrow(X)
##  ##  C there is an identity 	##  which returns the identity 1:X-->X. 
##  ##  arrow.			##
##  ##				##
##  #####################################################################
##  ##				##
##  ##  An arrow f in the       ##  Component object f with attributes          
##  ##  category C.	        ##  * CategoryName(f)
##  ##				##  * Source(f)  a category object
##  ##				##  * Target(f)  a category object
##  ##				##  * Mapping(f) such that something like
##  ##				##    IsGroupHomomorphism() is true
##  ##				##
##  #####################################################################
##  ##				##
##  ##  Categorical property    ##  A GAP property of X. 
##  ##  of an object or arrow   ## 
##  ##  X    			##  
##  ##				##  
##  ##				##
##  #####################################################################
##  ##				##
##  ##  Functors and Natural	##  Yet to be added. Functors are functions
##  ##  Transformations		##  with attributes such as HasFunctorName(X).
##  ## 				##
##  #####################################################################
##  ##				##
##  ##  Composition of arrows	##  g*f
##  ##  f:A->B and g:B->C	##
##  ##  and, in an additive 	##  and
##  ##  category, addition of	##
##  ##  arrows g,h:B->C		##  g+h
##  ##				##
##  #####################################################################
##
##  We plan to implement general homological constructions such as 
##  "total complex of a multi-complex" using a small and portable language. 
##  This language will apply to "diagrams in a category" enabling one, for
##  example, to transform a multi-complex into its total chain complex. 
##  A "diagram" will be a function D with an indexing set S and category C; 
##  for each s in S the function will yield an arrow D(s) in C.   



#############################################################################
##
##  Declare category object and arrows to be component objects.
##
DeclareProperty("IsCategoryObject", IsComponentObjectRep);
DeclareProperty("IsCategoryArrow",  IsComponentObjectRep);
DeclareProperty("IsCategoryName",  IsComponentObjectRep);

#############################################################################
##
##  Basic attributes of a CategoryObject
##
DeclareAttribute( "CategoryName" ,    IsCategoryObject );
DeclareAttribute( "Object" ,          IsCategoryObject  );

#############################################################################
##
##  Basic attributes of a CategoryArrow
##
DeclareAttribute( "CategoryName" ,   IsCategoryArrow );
DeclareAttribute( "Source" ,         IsCategoryArrow );
DeclareAttribute( "Target" ,         IsCategoryArrow );
DeclareAttribute( "Mapping" ,        IsCategoryArrow );

#############################################################################
##
##  Basic properties and attributes of  Categories
##
DeclareProperty("IsAbelianCategory", IsString);
DeclareProperty("IsAdditiveCategory", IsString);
DeclareAttribute("Pullbacks", IsString); #changed from Propertry
DeclareAttribute("Pushouts", IsString); #changed from Propertry

#############################################################################
##
##  Constructor for converting an object G into an  
##  object in the category "string". 
##  
##  CategoryObject(G,"string")
DeclareOperation("CategoryObject", [IsObject,IsString]);

#############################################################################
##
##  Constructor for converting a morphism F into an
##  arrow in the category "string".
##
##CategoryArrow(A,B,phi,"string")
DeclareOperation("CategoryArrow", 
  [IsCategoryObject,IsCategoryObject,IsObject,IsString]);

#############################################################################
##
##  Operations for associating the identity arrow, the initial arrow and
##  the terminal arrow to a category object X.
##
##  IdentityArrow(X)
##  InitialArrow(X)
##  TerminalArrow(X)
DeclareOperation("IdentityArrow", [IsCategoryObject]);
DeclareOperation("InitialArrow", [IsCategoryObject]);
DeclareOperation("TerminalArrow", [IsCategoryObject]);



#############################################################################
##
##  More sophisticated constructor for enriching objects and morphism with 
##  the structure of category "string".
##  
##  CategoricalEnrichment(X,"string").
DeclareOperation("CategoricalEnrichment", [IsObject,IsString]);
