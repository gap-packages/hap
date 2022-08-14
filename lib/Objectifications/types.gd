#(C) Graham Ellis, 2005-2006

IsHapChain:=NewFilter("IsHapChain");
IsHapCochain:=NewFilter("IsHapCochain");
IsHapMap:=NewFilter("IsHapMap");
IsHapComplex:=NewFilter("IsHapComplex");

#####################################################################
#####################################################################
DeclareCategory("IsHapResolution",IsObject);

DeclareRepresentation(	"IsHapResolutionRep",
			IsComponentObjectRep
			and IsHapComplex,
			["dimension",
			 "boundary",
			  "homotopy",
			  "group",
			  "elts",
			  "properties"]);

HapResolutionFamily:=NewFamily(	"HapResolutionFamily",
				IsHapResolution,
				IsHapResolution);

HapResolution:=NewType(HapResolutionFamily,IsHapResolutionRep);

InstallMethod( ViewObj,
"for HapResolution",
 [IsHapResolution],
 function(R)
 Print("Resolution of length ", EvaluateProperty(R,"length"),
 " in characteristic ", EvaluateProperty(R,"characteristic"), 
 " for "); ViewObj(R!.group); Print(" . \n");
 if R!.homotopy=fail then
 Print("No contracting homotopy available. \n");
 fi;
 if "partialHomotopy" in NamesOfComponents(R) then
 Print("A partial contracting homotopy is available. \n");
  fi;

 end);

InstallMethod( PrintObj,
"for HapResolution",
 [IsHapResolution],
  function(R)
   Print("Resolution of length ", EvaluateProperty(R,"length"),  
   " in characteristic ", EvaluateProperty(R,"characteristic"),
   " for ", R!.group," . \n");
if R!.homotopy=fail then
Print("No contracting homotopy available. \n");
fi;
if "partialHomotopy" in NamesOfComponents(R) then
 Print("A partial contracting homotopy is available. \n");
   fi;

    end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
DeclareCategory("IsHapEquivariantChainMap",IsObject);

DeclareRepresentation(  "IsHapEquivariantChainMapRep",
                        IsComponentObjectRep
			and IsHapMap,
                        ["source",                          
			 "target",
	                 "mapping",
	                 "properties"]);

HapEquivariantChainMapFamily:=NewFamily( "HapEquivariantChainMapFamily",
					  IsHapEquivariantChainMap,
					  IsHapEquivariantChainMap);

HapEquivariantChainMap:=NewType(HapEquivariantChainMapFamily,
				IsHapEquivariantChainMapRep);


InstallMethod( ViewObj,
"for HapEquivariantChainMap",
 [IsHapEquivariantChainMap],
function(R)
Print("Equivariant Chain Map between resolutions of length ",
Minimum(EvaluateProperty(R!.source,"length"),
EvaluateProperty(R!.target,"length")), " . \n");
end);

InstallMethod( PrintObj,
"for HapEquivariantChainMap",
 [IsHapEquivariantChainMap],
 function(R)
 Print("Equivariant Chain Map  \n");
 end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareCategory("IsHapChainComplex",IsObject);

DeclareRepresentation( "IsHapChainComplexRep",
                        IsComponentObjectRep 
			and IsHapChain
			and IsHapComplex,
			["dimension",
			 "boundary",
			 "properties"]);
			 
HapChainComplexFamily:=NewFamily( "HapChainComplexFamily",
				   IsHapChainComplex,
				   IsHapChainComplex);

HapChainComplex:=NewType(HapChainComplexFamily,IsHapChainComplexRep);

InstallMethod( ViewObj,
"for HapChainComplex",
[IsHapChainComplex],
 function(R)
Print("Chain complex of length ",
EvaluateProperty(R,"length"), " in characteristic ",
EvaluateProperty(R,"characteristic"), " . \n");
 end);

InstallMethod( PrintObj,
"for HapChainComplex",
[IsHapChainComplex],
function(R)
Print("Chain complex of length ",
EvaluateProperty(R,"length"), " in characteristic ",
EvaluateProperty(R,"characteristic"), " . \n");


end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareCategory("IsHapFilteredChainComplex",IsObject);

DeclareRepresentation( "IsHapFilteredChainComplexRep",
                        IsComponentObjectRep
                        and IsHapComplex and IsHapChain and
                        IsHapChainComplex,
                        ["dimension",
                         "filteredDimension",
                         "boundary",
                         "properties"]);

HapFilteredChainComplexFamily:=NewFamily( "HapFilteredChainComplexFamily",
                                   IsHapFilteredChainComplex,
                                   IsHapFilteredChainComplex);

HapFilteredChainComplex:=NewType(HapFilteredChainComplexFamily,IsHapFilteredChainComplexRep);

InstallMethod( ViewObj,
"for HapFilteredChainComplex",
[IsHapFilteredChainComplex],
 function(R)
Print("Filtered chain complex of length ",
EvaluateProperty(R,"length"), " in characteristic ",
EvaluateProperty(R,"characteristic"), " . \n");
 end);

InstallMethod( PrintObj,
"for HapFilteredChainComplex",
[IsHapFilteredChainComplex],
function(R)
Print("Filtered chain complex of length ",
EvaluateProperty(R,"length"), " in characteristic ",
EvaluateProperty(R,"characteristic"), " . \n");


end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
DeclareCategory("IsHapCochainComplex",IsObject);

DeclareRepresentation( "IsHapCochainComplexRep",
                        IsComponentObjectRep 
			and IsHapCochain
			and IsHapComplex,
			["dimension",
			"boundary",
			"properties"]); 
									 
HapCochainComplexFamily:=NewFamily( "HapCochainComplexFamily",
				     IsHapCochainComplex,
				     IsHapCochainComplex);

HapCochainComplex:=NewType(HapCochainComplexFamily,IsHapCochainComplexRep);

InstallMethod( ViewObj,
"for HapCochainComplex",
 [IsHapCochainComplex],
 function(R)
 Print("Cochain complex of length ",
 EvaluateProperty(R,"length"), " in characteristic ",
 EvaluateProperty(R,"characteristic"), " . \n");

 end);

InstallMethod( PrintObj,
"for HapCochainComplex",
[IsHapCochainComplex],
function(R)
Print("Cochain complex of length ",
EvaluateProperty(R,"length"), " in characteristic ",
EvaluateProperty(R,"characteristic") , " . \n");

end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareCategory("IsHapChainMap",IsObject);

DeclareRepresentation(  "IsHapChainMapRep",
                        IsComponentObjectRep 
			and IsHapChain
			and IsHapMap,
		   	["source",
			"target",
			"mapping",
			"properties"]);

HapChainMapFamily:=NewFamily( "HapChainMapFamily",
			       IsHapChainMap,
			       IsHapChainMap);

HapChainMap:=NewType(HapChainMapFamily, IsHapChainMapRep);


InstallMethod( ViewObj,
"for HapChainMap",
 [IsHapChainMap],
 function(R)
  Print("Chain Map between complexes of length ",
   Minimum(EvaluateProperty(R!.source,"length"),
    EvaluateProperty(R!.target,"length")), " . \n");

 end);

InstallMethod( PrintObj,
"for HapChainMap",
 [IsHapChainMap],
 function(R)
 Print("Chain Map between complexes of length ",
 Minimum(EvaluateProperty(R!.source,"length"),
 EvaluateProperty(R!.target,"length")), " . \n");
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareCategory("IsHapCochainMap",IsObject);

DeclareRepresentation(  "IsHapCochainMapRep",
                        IsComponentObjectRep 
			and IsHapCochain
			and IsHapMap,
                        ["source",
		         "target",
		         "mapping",
		         "properties"]);
			
HapCochainMapFamily:=NewFamily( "HapCochainMapFamily",
				IsHapCochainMap,
				IsHapCochainMap);  

HapCochainMap:=NewType(HapCochainMapFamily, IsHapCochainMapRep);


InstallMethod( ViewObj,
"for HapCochainMap",
 [IsHapCochainMap],
 function(R)
 Print("Cochain Map between complexes of length ",
  Minimum(EvaluateProperty(R!.source,"length"),
   EvaluateProperty(R!.target,"length")), " . \n");

end);

InstallMethod( PrintObj,
"for HapCochainMap",
 [IsHapCochainMap],
 function(R)
 Print("Cochain Map between complexes of length ",
   Minimum(EvaluateProperty(R!.source,"length"),
      EvaluateProperty(R!.target,"length")), " . \n");
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
IsHapNonFreeResolution:=NewFilter("IsHapNonFreeResolution");;
HapNonFreeResolution:=NewType(FamilyObj(rec()),
                   IsHapNonFreeResolution 
		   and IsComponentObjectRep
		   and IsHapComplex);;

InstallMethod( ViewObj,
"for HapNonFreeResolution",
[IsHapNonFreeResolution],
function(R)
Print("Non-free resolution in characteristic ", EvaluateProperty(R,"characteristic"),
 " for "); ViewObj(R!.group); Print(" . \n");
if R!.homotopy=fail then
Print("No contracting homotopy available. \n");
fi;
 end);

InstallMethod( PrintObj,
"for HapNonFreeResolution",
[IsHapNonFreeResolution],
function(R)
Print("Non-free resolution of length ", 
EvaluateProperty(R,"length"),
" in characteristic ", EvaluateProperty(R,"characteristic"),
" for ", R!.group," . \n");
if R!.homotopy=fail then
Print("No contracting homotopy available. \n");
fi;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallTrueMethod(IsHapChain,IsHapChainComplex);
InstallTrueMethod(IsHapChain,IsHapChainMap);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
IsHapFPGModule:=NewFilter("IsHapFPGModule");;
HapFPGModule:=NewType(FamilyObj(rec()),
IsHapFPGModule
and IsComponentObjectRep);;
InstallMethod( ViewObj,
"for HapFPGModule",
[IsHapFPGModule],
function(M)
Print("Module of dimension "); ViewObj(M!.dimension);
Print(" over the group ring of "); ViewObj(M!.group);
Print(" in characteristic ", (M!.characteristic), " \n");
end);
  
InstallMethod( PrintObj,
"for HapFPGModule",
[IsHapFPGModule],
function(M)
Print("Module of dimension ", M!.dimension, " over the group ring of ", M!.group, 
" in characteristic ", M!.characteristic, " \n");
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
IsHapFPGModuleHomomorphism:=NewFilter("IsHapFPGModuleHomomorphism");;
HapFPGModuleHomomorphism:=NewType(FamilyObj(rec()),
IsHapFPGModuleHomomorphism
and IsComponentObjectRep);;
InstallMethod( ViewObj,
"for HapFPGModuleHomomorphism",
[IsHapFPGModuleHomomorphism],
function(M)
Print("FpG-module homomorphism  "); 
end);

InstallMethod( PrintObj,
"for HapFPGModuleHomomorphism",
[IsHapFPGModuleHomomorphism],
function(M)
Print("FpG-module homomorphism from ", M!.source, " to ", M!.target,
 " \n");
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
DeclareCategory("IsHapCatOneGroup",IsObject);

DeclareRepresentation(  "IsHapCatOneGroupRep",
                        IsComponentObjectRep,
                        ["sourceMap",
                         "targetMap"
                         ]);

HapCatOneGroupFamily:=NewFamily( "HapCatOneGroupFamily",
                                          IsHapCatOneGroup,
                                          IsHapCatOneGroup);

HapCatOneGroup:=NewType(HapCatOneGroupFamily,
                                IsHapCatOneGroupRep);


InstallMethod( ViewObj,
"for HapCatOneGroup",
 [IsHapCatOneGroup],
function(R)
Print("Cat-1-group with underlying group ",
Source(R!.sourceMap),
 " . \n");
end);

InstallMethod( PrintObj,
"for HapCatOneGroup",
 [IsHapCatOneGroup],
 function(R)
Print("Cat-1-group with underlying group ",
Source(R!.sourceMap),
 " . \n");
 end);
#####################################################################
#####################################################################

#####################################################################
DeclareCategory("IsHapGraph",IsObject);

DeclareRepresentation(  "IsHapGraphRep",
                        IsComponentObjectRep,
                        ["incidenceMatrix",
                         ]);

HapGraphFamily:=NewFamily( "HapGraphFamily",
                                          IsHapGraph,
                                          IsHapGraph);

HapGraph:=NewType(HapGraphFamily,
                                IsHapGraphRep);


InstallMethod( ViewObj,
"for HapGraph",
 [IsHapGraph],
function(R)
Print("Graph on ", EvaluateProperty(R,"numberofvertices")," vertices.\n");
end);

InstallMethod( PrintObj,
"for HapGraph",
 [IsHapGraph],
 function(R)
Print("Graph on ", EvaluateProperty(R,"numberofvertices")," vertices.\n");
 end);
#####################################################################
#####################################################################

#####################################################################
DeclareCategory("IsHapFilteredGraph",IsObject);

DeclareRepresentation(  "IsHapFilteredGraphRep",
                        IsComponentObjectRep,
                        ["incidenceMatrix",
                         ]);

HapFilteredGraphFamily:=NewFamily( "HapFilteredGraphFamily",
                                          IsHapFilteredGraph,
                                          IsHapFilteredGraph);

HapFilteredGraph:=NewType(HapFilteredGraphFamily,
                                IsHapFilteredGraphRep);


InstallMethod( ViewObj,
"for HapFilteredGraph",
 [IsHapFilteredGraph],
function(R)
Print("Filtered graph on ", EvaluateProperty(R,"numberofvertices")," vertices.\n");
end);

InstallMethod( PrintObj,
"for HapFilteredGraph",
 [IsHapFilteredGraph],
 function(R)
Print("Filtered graph on ", EvaluateProperty(R,"numberofvertices")," vertices.\n");
 end);
#####################################################################
#####################################################################




#####################################################################
DeclareCategory("IsHapOppositeElement",IsMultiplicativeElementWithInverse);

DeclareRepresentation(  "IsHapOppositeElementRep",
                        IsComponentObjectRep,
                        IsMultiplicativeElementWithInverse,
                        ["oppositeElement",
                         ]);

HapOppositeElementFamily:=NewFamily( "HapOppositeElementFamily",
                                          IsHapOppositeElement,
                                          IsHapOppositeElement);

HapOppositeElement:=NewType(HapOppositeElementFamily,
                                IsHapOppositeElementRep);


InstallMethod( ViewObj,
"for HapOppositeElement",
 [IsHapOppositeElement],
function(R)
Print(R!.element, "_op ");
end);

InstallMethod( PrintObj,
"for HapOppositeElement",
 [IsHapOppositeElement],
 function(R)
Print(R!.element, "_op ");
 end);
#####################################################################
#####################################################################

#####################################################################
DeclareCategory("IsHapQuotientElement",IsMultiplicativeElementWithInverse);

DeclareRepresentation(  "IsHapQuotientElementRep",
                        IsComponentObjectRep and IsMultiplicativeElementWithInverse,
                        ["oppositeElement",
                         ]);

HapQuotientElementFamily:=NewFamily( "HapQuotientElementFamily",
                                          IsHapQuotientElement,
                                          IsHapQuotientElement);

HapQuotientElement:=NewType(HapQuotientElementFamily,
                                IsHapQuotientElementRep);


InstallMethod( ViewObj,
"for HapQuotientElement",
 [IsHapQuotientElement],
function(R)
Print(R!.element, "_quotient ");
end);

InstallMethod( PrintObj,
"for HapQuotientElement",
 [IsHapQuotientElement],
 function(R)
Print(R!.element, "_quotient ");
 end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
IsHapGChainComplex:=NewFilter("IsHapGChainComplex");;
HapGChainComplex:=NewType(FamilyObj(rec()),
                   IsHapGChainComplex
                   and IsComponentObjectRep
                   and IsHapComplex);;

InstallMethod( ViewObj,
"for HapGChainComplex",
[IsHapGChainComplex],
function(R)
Print("G-chain complex in characteristic ", EvaluateProperty(R,"characteristic"),
 " for "); ViewObj(R!.group); Print(" . \n");
 end);

InstallMethod( PrintObj,
"for HapGChainComplex",
[IsHapGChainComplex],
function(R)
Print("G-chain complex of length ",
EvaluateProperty(R,"length"),
" in characteristic ", EvaluateProperty(R,"characteristic"),
" for ", R!.group," . \n");
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
IsHapRightTransversalSL2ZSubgroup:=NewFilter("IsHapRightTransversalSL2ZSubgroup");;

HapRightTransversalSL2ZSubgroup:=NewType(FamilyObj(rec()),
                   IsHapRightTransversalSL2ZSubgroup
                   and IsList
                   );;

InstallMethod( ViewObj,
"for HapRightTransversalSL2ZSubgroup",
[IsHapRightTransversalSL2ZSubgroup],
function(R)
Print("Transversal of "); 
  ViewObj(R!.subgroup); Print(" in "); ViewObj(R!.group);Print(" . \n");
 end);

InstallMethod( PrintObj,
"for HapTransversalSL2ZSubgroup",
[IsHapRightTransversalSL2ZSubgroup],
function(R)
Print("Transversal of ");
  ViewObj(R!.subgroup); Print(" in "); ViewObj(R!.group);Print(" . \n");
 end);

#####################################################################
#####################################################################


#####################################################################
#####################################################################
IsHapSL2Subgroup:=NewFilter("IsHapSL2Subgroup");;
IsHapSL2ZSubgroup:=NewFilter("IsHapSL2ZSubgroup");;

HapSL2ZSubgroup:=NewType(FamilyObj(rec()),
                   IsHapSL2Subgroup
                   and IsHapSL2ZSubgroup
                   and IsGroup
                   );;

InstallMethod( ViewObj,
"for HapSL2ZSubgroup",
[IsHapSL2ZSubgroup and IsGroup],
function(G)
Print(G!.name,"(",G!.level,") ");
 end);

InstallMethod( PrintObj,
"for HapSL2ZSubgroup",
[IsHapSL2ZSubgroup and IsGroup],
function(G)
Print(G!.name,"(",G!.level,") ");
 end);

#####################################################################
#####################################################################

#####################################################################
#####################################################################
IsHapSL2ConjugatedSubgroup:=NewFilter("IsHapSL2ConjugatedSubgroup");;
IsHapSL2ZConjugatedSubgroup:=NewFilter("IsHapSL2ZConjugatedSubgroup");;

HapSL2ZConjugatedSubgroup:=NewType(FamilyObj(rec()),
                   IsHapSL2ConjugatedSubgroup
                   and IsHapSL2ZConjugatedSubgroup
                   and IsGroup
                   );;

InstallMethod( ViewObj,
"for HapSL2ZConjugatedSubgroup",
[IsHapSL2ZConjugatedSubgroup and IsGroup],
function(G)
Print(G!.name,"(",G!.level,") ");
 end);

InstallMethod( PrintObj,
"for HapSL2ZConjugatedSubgroup",
[IsHapSL2ZConjugatedSubgroup and IsGroup],
function(G)
Print(G!.name,"(",G!.level,") ");
 end);

#####################################################################
#####################################################################

#####################################################################
#####################################################################
IsHapSL2OSubgroup:=NewFilter("IsHapSL2OSubgroup");;

HapSL2ZSubgroup:=NewType(FamilyObj(rec()),
                   IsHapSL2Subgroup
                   and IsHapSL2OSubgroup
                   and IsGroup
                   );;

InstallMethod( ViewObj,
"for HapSL2OSubgroup",
[IsHapSL2OSubgroup and IsGroup],
function(G)
Print(G!.name,"(",G!.level,") ");
 end);

InstallMethod( PrintObj,
"for HapSL2OSubgroup",
[IsHapSL2OSubgroup and IsGroup],
function(G)
Print(G!.name,"(",G!.level,") ");
 end);

#####################################################################
#####################################################################



