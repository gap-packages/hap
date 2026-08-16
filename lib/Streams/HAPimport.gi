#(C) Graham Ellis 2008

HAPTEMPORARYFUNCTION:=0;
MakeReadOnlyGlobal("HAPTEMPORARYFUNCTION");

#####################################################################
#####################################################################
InstallMethod(HAPRead,
   "method for reading in a HAp resolution from a text file",
   [IsString],

function(file)
	local 	R,
		Dimension,
		ChangeSign,
		Boundary,
		Homotopy;

##################### READ FILE TO R ################################
MakeReadWriteGlobal("HAPTEMPORARYFUNCTION");
Read(file);
if IsFunction(HAPTEMPORARYFUNCTION) then
R:=HAPTEMPORARYFUNCTION();
else
Print("YES\n");
R:=HAPTEMPORARYFUNCTION;
fi;
HAPTEMPORARYFUNCTION:=0;
MakeReadOnlyGlobal("HAPTEMPORARYFUNCTION");
##################### FILE READ TO R ################################

if IsList(R) then return R; fi;


if FiltersType(R.type)=FiltersType(HapResolution) then 
###########################HapResolution#############################
Dimension:=function(i)
if i<0 then return 0; fi;
if i=0 then return 1; fi;
return R.ranks[i+1];
end;
#####################################################################

#####################################################################
ChangeSign:=function(j,b)
if j>0 then return b; else
return NegateWord(b); fi;
end;
#####################################################################

#####################################################################
Boundary:=function(i,j)
if i=0 then return []; else
return ChangeSign(j,R.boundaries[i][AbsoluteValue(j)]); fi;
end;
#####################################################################

#####################################################################
Homotopy:=function(n,p)        
#if i <0 then return fail; fi;   
if p[1]>0 then 
return R.contracting_homotopy[n+1][p[1]][p[2]];
else
return NegateWord(R.contracting_homotopy[n+1][-p[1]][p[2]]);
fi;
end;
#####################################################################


return Objectify(HapResolution,
                rec(
                dimension:=Dimension,
                boundary:=Boundary,
                homotopy:=Homotopy,
                elts:=R.elements,
                group:=Group(R.elements),
                properties:=R.special_properties
                    ));
###############End HapResolution#####################################
fi;

if FiltersType(R.type)=FiltersType(HapSparseChainComplex) then
###########################HapSparseChainComplex#####################
Dimension:=function(i)
if i<0 or i+1>Length(R.ranks) then return 0; fi;
return R.ranks[i+1];
end;
#####################################################################


#####################################################################
Boundary:=function(i,j)
if i=0 or i>Length(R.boundaries) then return []; else
return R.boundaries[i][j]; fi;
end;
#####################################################################

return Objectify(HapSparseChainComplex,
                rec(
                dimension:=Dimension,
                boundary:=Boundary,
                properties:=R.special_properties
                    ));
###############End HapSparseChainComplex#############################
fi;


end);
#####################################################################
#####################################################################
