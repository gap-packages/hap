#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionAbelianGroup,
function(arg)
local  
	ResolutionFreeAbelianGroup,
	ResolutionAbelianInvariants,
	ResolutionAbGroup;

#####################################################################
#####################################################################
ResolutionFreeAbelianGroup:=function(d,n)
local 
	gens,
	F, FrAb, G,
	FhomFrAb, GhomFrAb,
	FirstProj,SecondProj,
	R,T,k;

if d=0 then return ResolutionFiniteGroup(Group([()]),n); fi;

F:=FreeGroup(1);
R:=ResolutionAsphericalPresentation(F,[],n);

FrAb:=CoxeterDiagramFpArtinGroup(List([1..d],i->[i]));
FrAb:=FrAb[1]/FrAb[2];
gens:=GeneratorsOfGroup(FrAb);
FhomFrAb:=[];
FirstProj:=[];
SecondProj:=[];

for k in [1..d] do
FhomFrAb[k]:=GroupHomomorphismByImagesNC(F,FrAb,[F.1],[gens[k]]);
od;

T:=R;
for k in [2..d] do
T:=ResolutionDirectProduct(R,T);
FirstProj[k]:=T.firstProjection;
SecondProj[k]:=T.secondProjection;
od;

	#############################################################
	GhomFrAb:=function(k,g)
	local g1, g2;
	
	if k=1 then return Image(FhomFrAb[d],g); fi;
	
	g1:=Image(FhomFrAb[d-k+1],Image(FirstProj[k],g));
	g2:=Image(SecondProj[k],g);
	g2:=GhomFrAb(k-1,g2);
	return g1*g2;
	end;
	#############################################################

T.elts:=List(T.elts,x->GhomFrAb(d,x));
T.group:=FrAb;

return T;
end;
#####################################################################
#####################################################################


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
return ResolutionFiniteGroup(C,n);
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
	if Sum(arg[1])=0 then return 
	ResolutionFreeAbelianGroup(Length(arg[1]),arg[2]);
	else
	return ResolutionAbelianInvariants(arg[1],arg[2]);
	fi;
fi;

if IsGroup(arg[1]) and IsInt(arg[2]) then
if IsFinite(arg[1]) then
return ResolutionAbGroup(arg[1],arg[2]); fi; fi;

Print("The first argument must be a list of nonnegative integers or a finite abelian group. The second argument must be a positive integer. \n");
return fail;
end);
