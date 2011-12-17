#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionAbelianGroup,
function(arg)
local  
	ResolutionAbelianInvariants,
	ResolutionAbGroup;

#####################################################################
#####################################################################
ResolutionAbelianInvariants:=function(coeffs,n)
local	head,tail,R;

if Length(coeffs)=0 then return 
ResolutionFiniteGroup(CyclicGroup(1),n); fi;

if Length(coeffs)=1 then
	if coeffs[1]=0 then 
	return ResolutionAsphericalPresentation(FreeGroup(1),[],n); 
	else
	return ResolutionFiniteGroup(CyclicGroup(coeffs[1]),n);
	fi;
fi;

head:=[coeffs[1]];
tail:=List([2..Length(coeffs)],i->coeffs[i]);

return ResolutionDirectProduct(ResolutionAbelianInvariants(head,n),
	ResolutionAbelianInvariants(tail,n));

end;
#####################################################################
#####################################################################

#####################################################################
#####################################################################
ResolutionAbGroup:=function(G,n)
local gens, C, head, tail, R, hom ;

gens:=TorsionGeneratorsAbelianGroup(G);

if Length(gens)=0 then return
ResolutionFiniteGroup([Identity(G)],n); fi;

if Length(gens)=1 then
C:=Group(gens[1]);
        if IsFreeGroup(C) then
        return ResolutionAsphericalPresentation(FreeGroup(1),[],n);
		#Need to fix this as it will cause problems!!
        else
return ResolutionFiniteGroup(C,n);
        fi;
fi;

head:=Group([gens[1]]);
tail:=Group(List([2..Length(gens)],i->gens[i]));

R:=ResolutionDirectProduct(ResolutionAbGroup(head,n),
		        ResolutionAbGroup(tail,n));
hom:=GroupHomomorphismByFunction(R.group,G,x->
Image(Projection(R.group,1),x)*
Image(Projection(R.group,2),x));
R.elts:=List(R.elts,x->Image(hom,x));
R.group:=G;
return R;

end;
#####################################################################
#####################################################################

if IsList(arg[1]) and IsInt(arg[2]) then
return ResolutionAbelianInvariants(arg[1],arg[2]); fi;

if IsGroup(arg[1]) and IsInt(arg[2]) then
return ResolutionAbGroup(arg[1],arg[2]); fi;

Print("The first argument must be a list of nonnegative integers or an abelian group. The second argument must be a positive integer. \n");
return fail;
end);
