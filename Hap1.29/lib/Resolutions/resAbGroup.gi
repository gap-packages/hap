#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionAbelianGroup_alt,
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

if HAPconstant>49 then
T:=ResolutionDirectProduct(R,T);
else
T:=ResolutionFiniteDirectProduct(R,T);
fi;
FirstProj[k]:=T!.firstProjection;
SecondProj[k]:=T!.secondProjection;
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

T!.elts:=List(T!.elts,x->GhomFrAb(d,x));
T!.group:=FrAb;

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

#if Minimum(coeffs)=0 then
if HAPconstant>49 then 
return ResolutionDirectProduct(ResolutionAbelianInvariants(head,n),
	ResolutionAbelianInvariants(tail,n));
else
return ResolutionFiniteDirectProduct(ResolutionAbelianInvariants(head,n),
        ResolutionAbelianInvariants(tail,n));
fi;

end;
#####################################################################
#####################################################################

#####################################################################
#####################################################################
ResolutionAbGroup:=function(G,n)
local gens, C, head, tail, R, hom, OriginalElts,OriginalGroup ;

if Order(G)=1 then
   if not IsPcpGroup(G) then
   R:=ResolutionFiniteGroup(Group(Identity(G)),n);
   return R;
   fi;
#The following code gets around the bug that Elements(G) does not
#work for a pcp group G
R:=ResolutionFiniteGroup(Group(()),n);
R!.group:=G;
R!.elts:=[Identity(G)];
return R;
fi;

gens:=TorsionGeneratorsAbelianGroup(G); 

#if Length(gens)=0 then return
#ResolutionFiniteGroup([Identity(G)],n); fi;

if Length(gens)=1 then
C:=Group(gens[1]);
return ResolutionFiniteGroup(C,n);
fi;

head:=Group([gens[1]]);
tail:=Group(List([2..Length(gens)],i->gens[i]));

#if Minimum(AbelianInvariants(G))=0 then
if HAPconstant>49 then 
R:=ResolutionDirectProduct(ResolutionAbGroup(head,n),
		        ResolutionAbGroup(tail,n));
else
R:=ResolutionFiniteDirectProduct(ResolutionAbGroup(head,n),
                        ResolutionAbGroup(tail,n));
fi;

OriginalElts:=R!.elts;
if IsMutable(R!.elts) then
Append(R!.elts,Elements(R!.group));
fi;
OriginalGroup:=R!.group;

hom:=GroupHomomorphismByFunction(OriginalGroup,G,x->
Image(Projection(OriginalGroup,1),x)*
Image(Projection(OriginalGroup,2),x));

R!.elts:=List(R!.elts,x->Image(hom,x));
R!.group:=G;
return R;

end;
#####################################################################
#####################################################################

if IsList(arg[1]) and IsInt(arg[2]) then
	if Sum(arg[1])=0 and Length(arg[1])>1 then return 
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
