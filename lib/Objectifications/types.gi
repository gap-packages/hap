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
Print("Non-free resolution of length ", EvaluateProperty(R,"length"),
" in characteristic ", EvaluateProperty(R,"characteristic"),
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
