######################################
#	NumberSmallQuasiCatOneGroups
#	SmallQuasiCatOneGroup
#	IdQuasiCatOneGroup
#	NumberSmallQuasiCrossedModules
#	SmallQuasiCrossedModule
#	IdQuasiCrossedModule
#######################################


#############################################################################
#0
#F	NumberSmallQuasiCatOneGroups
##
InstallGlobalFunction(NumberSmallQuasiCatOneGroups, function(n)
	
	if (n > QUASICATONEGROUP_DATA_SIZE) or (n in QUASICATONEGROUP_DATA_NOT) then
		Print("This function only apply for order < ",QUASICATONEGROUP_DATA_SIZE+1);
		Print(" and not in ",QUASICATONEGROUP_DATA_NOT,"\n");
		return fail;
	fi;
	return Length(SMALLQUASICATONEGROUP_DATA[n]);
end);
##
#############################################################################

#############################################################################
#0
#F	SmallQuasiCatOneGroup
##
InstallGlobalFunction(SmallQuasiCatOneGroup, function(n,k)
local m,x;

	if (n > QUASICATONEGROUP_DATA_SIZE) or (n in QUASICATONEGROUP_DATA_NOT) then
		Print("This function only apply for order < ",QUASICATONEGROUP_DATA_SIZE+1);
		Print(" and not in ",QUASICATONEGROUP_DATA_NOT,"\n");
		return fail;
	fi;
	m:=Length(SMALLQUASICATONEGROUP_DATA[n]);
	if k> m then
		Print("This function only apply for order < ",QUASICATONEGROUP_DATA_SIZE+1);
		Print(" and not in ",QUASICATONEGROUP_DATA_NOT,"\n");
		return fail;
	fi;
	x:=SMALLQUASICATONEGROUP_DATA[n][k];
    return SmallCatOneGroup(n,x[1],x[2]);
end);
##
#############################################################################

#############################################################################
#0
#F	IdQuasiCatOneGroup
##
InstallGlobalFunction(IdQuasiCatOneGroup, function(C)
local D,x;

	D:=QuasiIsomorph(C);
	x:=IdCatOneGroup(D);
	if (x[1] > QUASICATONEGROUP_DATA_SIZE) or (x[1] in QUASICATONEGROUP_DATA_NOT) then
		Print("This function only apply for order < ",QUASICATONEGROUP_DATA_SIZE+1);
		Print(" and not in ",QUASICATONEGROUP_DATA_NOT,"\n");
		return fail;
	fi;
	return IDQUASICATONEGROUP_DATA[x[1]][x[2]][x[3]];
end);
##
#############################################################################

#############################################################################
#0
#F	NumberSmallQuasiCrossedModules
##
InstallGlobalFunction(NumberSmallQuasiCrossedModules, function(n)
	
	if (n > QUASICATONEGROUP_DATA_SIZE) or (n in QUASICATONEGROUP_DATA_NOT) then
		Print("This function only apply for order < ",QUASICATONEGROUP_DATA_SIZE+1);
		Print(" and not in ",QUASICATONEGROUP_DATA_NOT,"\n");
		return fail;
	fi;
	return Length(SMALLQUASICATONEGROUP_DATA[n]);
end);
##
#############################################################################

#############################################################################
#0
#F	SmallQuasiCrossedModule
##
InstallGlobalFunction(SmallQuasiCrossedModule, function(n,k)
local m,x;

	if (n > QUASICATONEGROUP_DATA_SIZE) or (n in QUASICATONEGROUP_DATA_NOT) then
		Print("This function only apply for order < ",QUASICATONEGROUP_DATA_SIZE+1);
		Print(" and not in ",QUASICATONEGROUP_DATA_NOT,"\n");
		return fail;
	fi;
	m:=Length(SMALLQUASICATONEGROUP_DATA[n]);
	if k> m then
		Print("There are only ",m," quasi-isomorphism classes of order ",n,"\n");
		return fail;
	fi;
	x:=SMALLQUASICATONEGROUP_DATA[n][k];
    return CrossedModuleByCatOneGroup(SmallCatOneGroup(n,x[1],x[2]));
end);
##
#############################################################################

#############################################################################
#0
#F	IdQuasiCrossedModule
##
InstallGlobalFunction(IdQuasiCrossedModule, function(X)
local C,x;

	C:=QuasiIsomorph(CatOneGroupByCrossedModule(X));
	x:=IdCatOneGroup(C);
	if (x[1] > QUASICATONEGROUP_DATA_SIZE) or (x[1] in QUASICATONEGROUP_DATA_NOT) then
		Print("This function only apply for order < ",QUASICATONEGROUP_DATA_SIZE+1);
		Print(" and not in ",QUASICATONEGROUP_DATA_NOT,"\n");
		return fail;
	fi;
	return IDQUASICATONEGROUP_DATA[x[1]][x[2]][x[3]];
end);
##
#############################################################################














