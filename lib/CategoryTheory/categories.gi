#############################################################################
##
#W  categories.gi                 HAP                            Graham Ellis
##                                                               
##
##  2008-06-27 

#############################################################################
##
## The category of groups. As the implementation improves the properties
## of this category will increase! 
##
BindGlobal("Category_Of_Groups", Objectify(
           NewType(NewFamily("Category_Of_Groups"),
                       IsString), []));


#############################################################################
##
##  Enrichment of an object A to a category object A in the category 
##  "NAME" .  
##
InstallMethod( CategoryObject,
    "basic method for creating a category object",
    [ IsObject, IsString ],

    function( A, Name  )
    local obj,type;

		obj:=rec();
		type:= NewType(NewFamily("category_object"),
                       IsCategoryObject and
                       IsComponentObjectRep and IsAttributeStoringRep);
	ObjectifyWithAttributes(obj,type);
        SetCategoryName(obj,Name);
	#SetCategoryStatus(obj,true);
        SetObject(obj,A);
        return obj;
    end);

#############################################################################
##
##  Method for constructing the identity arrow on an object in category  
##  "NAME" .
##
InstallMethod( IdentityArrow,
    "method for creating the identity arrow on a category object",
    [ IsCategoryObject ],

    function( A )
	local phi;

    	if not CategoryName(A)=Category_Of_Groups then
	TryNextMethod(); fi;

	phi:=GroupHomomorphismByFunction(Object(A),Object(A),x->x);
	return CategoryArrow(A,A,phi,CategoryName(A));

	end);
	
#############################################################################
##
##  Method for constructing the unique arrow from an initial object to an 
##  object A in the category "NAME" .
##
InstallMethod( InitialArrow,
    "method for creating the initial arrow to a category object",
    [ IsCategoryObject ],

    function( A )
        local initial, phi, INITIAL;

        if not CategoryName(A)=Category_Of_Groups then
        TryNextMethod(); fi;

        initial:=Group( () );
        INITIAL:=CategoricalEnrichment(initial,CategoryName(A));
        phi:=GroupHomomorphismByFunction(initial,Object(A),
                                             x->Identity(Object(A)));
        return CategoryArrow(INITIAL,A,phi,CategoryName(A));

        end);
#########################################################################
##
##  Method for constructing the unique arrow from an object in the  
##  category "NAME" to an object A in the category "NAME" .
##
InstallMethod( TerminalArrow,
    "method for creating the initial arrow to a category object",
    [ IsCategoryObject ],

    function( A )
        local terminal, phi, TERMINAL;

        if not CategoryName(A)=Category_Of_Groups then
        TryNextMethod(); fi;

        terminal:=Group( () );
        TERMINAL:=CategoricalEnrichment(terminal,CategoryName(A));
        phi:=GroupHomomorphismByFunction(Object(A), terminal,
                                             x->Identity(terminal));
        return CategoryArrow(A,TERMINAL,phi,CategoryName(A));

        end);




#############################################################################
##
##  Install method for viewing and printing category objects.
##

InstallMethod( ViewObj,
        "for category Objects",
        [IsCategoryObject],
    function(X)
    Print("Object in ", FamilyObj(CategoryName(X))!.NAME, "\n");
end);

InstallMethod( PrintObj,
        "for category Objects",
        [IsCategoryObject],
    function(X)
    Print("Object in ", FamilyObj(CategoryName(X))!.NAME, "\n");
end);


#############################################################################
##
##  Enrichment of a "morphism" Phi:A-->B to a category arrow Phi in the
##  category "NAME". Both A and B must already be appropriate category
##  objects. No checking is done.
##
InstallMethod( CategoryArrow,
    "basic method for creating a category object",
    [ IsCategoryObject, IsCategoryObject, IsObject, IsString ], 20,

    function( A, B, phi, Name  )
	local arrow,type;
 
	arrow:=rec();
                type:= NewType(NewFamily("category_object"),
                       IsCategoryArrow and
                       IsComponentObjectRep and IsAttributeStoringRep);
        ObjectifyWithAttributes(arrow,type);
        SetSource(arrow,A);
        SetTarget(arrow,B);
        SetCategoryName(arrow,Name);
        SetMapping(arrow,phi);
        return arrow;
    end);

#############################################################################
##
##  Install method for viewing and printing category arrows.
##

InstallMethod( ViewObj,
        "for category arrows",
        [IsCategoryArrow],
    function(Phi)
    Print("Arrow in ", FamilyObj(CategoryName(Phi))!.NAME, "\n");
end);

InstallMethod( PrintObj,
        "for category arrows",
        [IsCategoryArrow],
    function(Phi)
    Print("Arrow in ", FamilyObj(CategoryName(Phi))!.NAME, "\n");
end);


#############################################################################
##
## Method for using Target instead of Range
##
InstallMethod( Target,
    "Method for making Target return Range",
    [ HasRange ],

    function( X );
    return Range(X);
    end);


#############################################################################
##
##  More sophisticated method for converting objects and morphism into  
##  the  category "string".
##
##  CategoricalEnrichment(X,"string").

InstallMethod( CategoricalEnrichment,
    "method for creating a category object",
    [ IsObject, IsString ],

    function( X , Name  );
        return
        CategoryObject(X,Name);
    end);

InstallMethod( CategoricalEnrichment,
    "method for creating a category arrow",
    [ HasSource , IsString ], #HasSource is not the best

    function( X , Name  );
        return
        CategoryArrow(
		CategoryObject(Source(X),Name),
		CategoryObject(Target(X),Name),
		X,Name);
    end);

#############################################################################
##
##  Install a category composition method on "*". So f(g(x)) can be written
##  f*g(x)
##
InstallOtherMethod( \*,
    "composition of arrows in a category",
    [ IsCategoryArrow, IsCategoryArrow], 

    function( PHI, THETA  )
        local PhiOfTheta;

 	if not Source(PHI)=Target(THETA) then
        Print("Arrows are not composable. \n");
        return fail; fi;

	if not CategoryName(PHI)=Category_Of_Groups then
	TryNextMethod(); fi;

        PhiOfTheta:=GroupHomomorphismByFunction(
        Object(Source(THETA)),
        Object(Target(PHI)),
        x -> Image(Mapping(PHI),Image(Mapping(THETA),x))  );
        
	PhiOfTheta:=CategoryArrow(
	Source(THETA),Target(PHI), PhiOfTheta, CategoryName(PHI));

        return PhiOfTheta;
    end);


#############################################################################
##
##  Install a category arrow addition method on "\+". So f(x)+g(x) can be 
##  (f+g)(x)
##
InstallOtherMethod( \+,
    "addition of arrows in a category",
    [ IsCategoryArrow, IsCategoryArrow],

    function( PHI, THETA  )
        local PhiOfTheta;

        if not 
	(Source(PHI)=Source(THETA)
	and
	Target(PHI)=Target(THETA))
	then
        Print("Arrows are not composable. \n");
        return fail; fi;

        if not CategoryName(PHI)=Category_Of_Groups then
        TryNextMethod(); fi;

        PhiOfTheta:=GroupHomomorphismByFunction(
        Object(Source(THETA)),
        Object(Target(PHI)),
        x -> Image(Mapping(PHI),x)
                    *
             Image(Mapping(THETA),x)  );

        PhiOfTheta:=CategoryArrow(
        Source(THETA),Target(PHI), PhiOfTheta, CategoryName(PHI));

        return PhiOfTheta;
    end);

#############################################################################
##
##  Install a category arrow subtraction method on "\-". So f(x)-g(x) can be
##  (f-g)(x)
##
InstallOtherMethod( \-,
    "addition of arrows in a category",
    [ IsCategoryArrow, IsCategoryArrow],

    function( PHI, THETA  )
        local PhiOfTheta;

        if not
        (Source(PHI)=Source(THETA)
        and
        Target(PHI)=Target(THETA))
        then
        Print("Arrows are not composable. \n");
        return fail; fi;

        if not CategoryName(PHI)=Category_Of_Groups then
        TryNextMethod(); fi;

        PhiOfTheta:=GroupHomomorphismByFunction(
        Object(Source(THETA)),
        Object(Target(PHI)),
        x -> Image(Mapping(PHI),x)
                    *
             Image(Mapping(THETA),x)^-1  );

        PhiOfTheta:=CategoryArrow(
        Source(THETA),Target(PHI), PhiOfTheta, CategoryName(PHI));

        return PhiOfTheta;
    end);

#############################################################################
##
##  Install an equality test on category objects 
##
##
InstallOtherMethod( \=,
    "equality of objects in a category",
    [ IsCategoryObject, IsCategoryObject],

    function( X, Y  );
	
	if not CategoryName(X)=CategoryName(Y)
	then return false; fi;
	if not Object(X)=Object(Y) then 
	return false; fi;
	return true;
    end);

#############################################################################
##
##  Install an equality test on category arrows 
##  
##
InstallOtherMethod( \=,
    "equality of arrows in a category",
    [ IsCategoryArrow, IsCategoryArrow],

    function( PHI, THETA  )
        local G, g;

        if not
        (Source(PHI)=Source(THETA)
        and
        Target(PHI)=Target(THETA))
	then
        return false; fi;

	if not CategoryName(PHI)=Category_Of_Groups then
	TryNextMethod(); fi;

	G:=Object(Source(PHI));	
	for g in GeneratorsOfGroup(G) do
	if not 
	Image(Mapping(PHI),g) = Image(Mapping(THETA),g)
	then return false;
	fi;
	od;
	return true;

    end);



