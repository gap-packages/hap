#############################################
InstallMethod(MooreComplex,
"Moore complex of simplicial group",
[IsHapSimplicialGroup],
function(N)
local 
	ListG,ListB,
	ListGroups,Boundaries,BoundariesList,
	n,i,j,
	G,H,
	phi;
	
	Boundaries:=[];
	ListGroups:=[];
	
	if not IsHapSimplicialGroup(N) then
		Print("This function must be applied to a simplicial group.\n");
		return fail;
	fi;
	ListG:=N!.groupsList;
	ListB:=N!.boundariesList;
	n:=EvaluateProperty(N,"length");
	H:=ListG(0);
	Add(ListGroups,H);
	for i in [1..n] do
		G:=ListG(i);
		for j in [1..i] do
			G:=Intersection(G,Kernel(ListB(i,j)));
		od;
		Add(ListGroups,G);
		phi:=GroupHomomorphismByImages(G,H,GeneratorsOfGroup(G),List(GeneratorsOfGroup(G),g->Image(ListB(i,0),g)));
		Add(Boundaries,phi);
		H:=G;
	od;	
Boundaries:=List(Boundaries,h->GOuterGroup(h));
#################	
BoundariesList:=function(i)
	return Boundaries[i];
end;
######################
return Objectify(HapGComplex,
       rec(
           boundary:=BoundariesList,
		   properties:=[["length",n]] 
         ));
end);
###################################
###################################


###################################
###################################
InstallMethod(Homology,
"Homology of a G Complex",
[IsHapGComplex,IsInt],
function(C,n) local phi1,phi2;

if n>0 then
	##############
	phi1:=C!.boundary(n);
	phi2:=C!.boundary(n+1);
return 
	AbelianInvariants(
	Kernel( phi1!.Mapping )
	/
	Image(  phi2!.Mapping )
	);
	##############
fi;

if n=0 then
##############
	phi2:=C!.boundary(n+1);
return
	AbelianInvariants(
	Range( phi2!.Mapping )
	/
	Image(  phi2!.Mapping )
	);
	##############
fi;

end);

############################
