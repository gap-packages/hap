###########################################################
##	HomotopyGroup
##	Size
##	Order
##	CrossedModuleByAutomorphismGroup
##	CrossedModuleByNormalSubgroup
##	CrossedModuleByCatOneGroup
##	HomotopyCrossedModule
##	NumberSmallCrossedModules
##	SmallCrossedModule
##	IdCrossedModule
##	IsomorphismCrossedModules
############################################################

#############################################################################
#0
#O	HomotopyGroup
##	Input:	A crossed module X and n=1,2
##	Output:	The nth homotopy groups of X 
##
InstallOtherMethod(HomotopyGroup, "Homology group of crossed modules",
[IsHapCrossedModule, IsInt], function(X,n)
local d;

	d:=X!.map;
	if n = 1 then
		return Range(d)/Image(d);
	fi;
	if n =2 then
		return Kernel(d);
	fi;
	Print("Only apply for n=1,2");
	return fail;
end);
##
#################### end of HomotopyGroup ###################################

#############################################################################
#0
#O	Size
##	Input:	A crossed module X
##	Output:	The size of X.
##
InstallOtherMethod(Size, "Size of crossed modules",
[IsHapCrossedModule], function(X)
	return Size(Source(X!.map))*Size(Range(X!.map));
end);
##
#################### end of Size ############################################

#############################################################################
#0
#O	Order
##	Input:	A crossed module X
##	Output:	The order of X.
##
InstallOtherMethod(Order, "Order of crossed modules",
[IsHapCrossedModule], function(X)
	return Order(Source(X!.map))*Order(Range(X!.map));	
end);
##
#################### end of Order ###########################################

#############################################################################
#0
#F	CrossedModuleByAutomorphismGroup
##	Input:	A group G
##	Output: The crossed module d:G->Aut(G)
##
InstallGlobalFunction(CrossedModuleByAutomorphismGroup, function(G)
local	AutG,GensG,d,act;
	
	AutG:=AutomorphismGroup(G);
	GensG:=GeneratorsOfGroup(G);
	d:=GroupHomomorphismByImages(G,AutG,GensG,List(GensG,
			g->GroupHomomorphismByImages(G,G,GensG,List(GensG,x->g^(-1)*x*g))));
	act:=function(f,g)
		return Image(f,g);
	end;
	return Objectify(HapCrossedModule,rec(
						map:=d,
						action:=act
					));
end);
##
################### end of CrossedModuleByAutomorphismGroup #################

#############################################################################
#0
#F	CrossedModuleByNormalSubgroup
##	Input:	A group G with normal subgroup N
##	Output: The inclusion crossed module i:N->G
##
InstallGlobalFunction(CrossedModuleByNormalSubgroup, function(G,N)
local	d,act;
	
	if not IsNormal(G,N) then
		Print("Only apply for a normal subgroup of group");
		return fail;
	fi;
	d:=GroupHomomorphismByFunction(N,G,x->x);
	act:=function(g,h)
		return g^(-1)*h*g;
	end;
	return Objectify(HapCrossedModule,rec(
						map:=d,
						action:=act
					));
end);
##
################### end of CrossedModuleByNormalSubgroup ####################

#############################################################################
#0
#F	CrossedModuleByCatOneGroup
##	Input:	Either a cat-1-group, or a morphism of cat-1-groups, or 
##			a sequence of morphisms of cat-1-groups
##	Output:	The image of the input under the functor 
##				(cat-1-groups)->(crossed modules)
##
InstallGlobalFunction(CrossedModuleByCatOneGroup,
function(X)
local
	Cat1ToCM_Obj,
	Cat1ToCM_Morpre,
	Cat1ToCM_Mor,
	Cat1ToCM_Seq;
	
    ######################################################################
	#1
	#F	Cat1ToCM_Obj
	##	Input:	A cat-1-group C
	##	Output:	The crossed module corresponds to C
	##
	Cat1ToCM_Obj:=function(C)
	local s, t,M,P,GensM,d,act;
	
		s:=C!.sourceMap;
		t:=C!.targetMap;
		M:=Kernel(s);
		P:=Image(s);
		GensM:=GeneratorsOfGroup(M);
		d:=GroupHomomorphismByImages(M,P,GensM,List(GensM,m->Image(t,m)));
		act:=function(p,m)
			return p^(-1)*m*p;
		end;
		return Objectify(HapCrossedModule,
						rec(map:=d,
						action:=act)
						);
	end;
	##
    ############### end of Cat1ToCM_Obj ##################################
	
	######################################################################
	#1
	#F	Cat1ToCM_Morpre
	##
	Cat1ToCM_Morpre:=function(XC,XD,f)  
	local 
		phiC,MC,PC,
		GensM,GensP,mapM,mapP,Map,
		phiD,MD,PD;
	
		phiC:=XC!.map;
		phiD:=XD!.map;
		MC:=Source(phiC);
		PC:=Range(phiC);
		MD:=Source(phiD);
		PD:=Range(phiD);
		GensM:=GeneratorsOfGroup(MC);
		GensP:=GeneratorsOfGroup(PC);
		mapM:=GroupHomomorphismByImages(MC,MD,GensM,List(GensM,m->Image(f,m)));
		mapP:=GroupHomomorphismByImages(PC,PD,GensP,List(GensP,p->Image(f,p)));
		
		#############################################################
		#2
		Map:=function(n)
			if n=1 then
				return mapM;
			fi;
			if n=2 then
				return mapP;
			fi;
			Print("Only apply for n =1,2");
			return fail;
		end;
		##
		#############################################################
		
		return Objectify(HapCrossedModuleMorphism,
						rec(source:=XC,
							target:=XD,
							mapping:=Map
							));
	end;
	##
	################ end of Cat1ToCM_Morpre ##############################
	
	######################################################################
	#1
	#F	Cat1ToCM_Mor
	##	Input:	A morphism fC of cat-1-groups 
	##	Output:	The morphism of crossed modules corresponds to fC
	##
	Cat1ToCM_Mor:=function(fC)
	local	XC,XD,f;
			
		XC:=Cat1ToCM_Obj(fC!.source);
		XD:=Cat1ToCM_Obj(fC!.tatget);
		f:=fC!.mapping;
		return Cat1ToCM_Morpre(XC,XD,f);	
	end;
	##
	############### end of Cat1ToCM_Mor ##################################
	
	######################################################################
	#1
	#F	Cat1ToCM_Seq
	##	Input:	A sequence LfC of morphisms of cat-1-groups 
	##	Output:	The sequence of morphisms of crossed modules corresponds to LfC
	##	
	Cat1ToCM_Seq:=function(LfC)
	local n,i,XC,Res;

		n:=Length(LfC);
		XC:=[];
		for i in [1..n] do
			XC[i]:=Cat1ToCM_Obj(LfC[i]!.source);
		od;
		XC[n+1]:=Cat1ToCM_Obj(LfC[i]!.target);
		Res:=[];
		for i in [1..n] do
			Res[i]:=Cat1ToCM_Morpre(XC[i],XC[i+1],LfC[i]!.mapping);
		od;
		return Res;
	end;
	##
    ############### end of Cat1ToCM_Seq ##################################
	
	if IsHapCatOneGroup(X) then
		return Cat1ToCM_Obj(X);
	fi;
	if IsHapCatOneGroupMorphism(X) then
		return Cat1ToCM_Mor(X);
	fi;
	if IsList(X) then
		return Cat1ToCM_Seq(X);
	fi;
end);
##
################### end of CrossedModuleByCatOneGroup #######################

#############################################################################
#0
##	HomotopyCrossedModule
##	Input: A crossed module X
##	Output: The crossed module 0:A->G
##
InstallGlobalFunction(HomotopyCrossedModule, function(X)
local	phi,act,P,A,nat,G,Gens,alpha;

	phi:=X!.map;
	act:=X!.action;
	P:=Range(phi);
	A:=Kernel(phi);
	nat:=NaturalHomomorphismByNormalSubgroup(P,Image(phi));
	G:=Range(nat);
	#################################################################
	#1
	alpha:=function(g,a)
	local x;
		
		x:=PreImagesRepresentative(nat,g);
		return act(x,a);
	end;
	##
	#################################################################
	Gens:=GeneratorsOfGroup(A);
	return Objectify(HapCrossedModule,rec(
				map:=GroupHomomorphismByImages(A,G,Gens,List(Gens,x->One(G))),
				action:=alpha
			));
end);
##
################### end of HomotopyCrossedModule ############################

#############################################################################
#0
#F	NumberSmallCrossedModules
##	Input: An integer n
##	Output: The number of non-isomorphic crossed modules of order n
##
InstallGlobalFunction(NumberSmallCrossedModules, function(n)

	if n >CATONEGROUP_DATA_SIZE then
		Print("This function only apply for order less than or equal to ",
				CATONEGROUP_DATA_SIZE,".\n");
		return fail;
	fi;
	return Sum(CATONEGROUP_DATA[n],x->Length(x));
end);
##
################### NumberSmallCrossedModules ###############################

#############################################################################
#0
#F	SmallCrossedModule
##	Input: An integers n,k.
##	Output: The kth crossed module of order n.
##
InstallGlobalFunction(SmallCrossedModule, function(n,k)
local sum,m,i;

	if n >CATONEGROUP_DATA_SIZE then
		Print("This function only apply for order less than or equal to ",
				CATONEGROUP_DATA_SIZE,".\n");
		return fail;
	fi;
	if k>NumberSmallCrossedModules(n) then
		Print("There are only ",NumberSmallCrossedModules(n),
				" crossed modules of order ",n,"\n");
		return fail;
	fi;
	sum:=0;
	m:=0;
	while sum<k do
		i:=k-sum;
		m:=m+1;
		sum:=sum+Length(CATONEGROUP_DATA[n][m]);
	od;
	return CrossedModuleByCatOneGroup(SmallCatOneGroup(n,m,i));
end);
##
################### SmallCrossedModule ######################################

#############################################################################
#0
#F	IdCrossedModule
##	Input: A crossed module X.
##	Output: The position of X.
##
InstallGlobalFunction(IdCrossedModule, function(X)
local T;

	if Order(X) > CATONEGROUP_DATA_SIZE then
		Print("This function only apply for crossed module of order less than or equal to",
				CATONEGROUP_DATA_SIZE,".\n");
		return fail;
	fi;
	T:=IdCatOneGroup(CatOneGroupByCrossedModule(X));
	return [T[1],Sum(List([1..T[2]-1],m->Length(CATONEGROUP_DATA[T[1]][m])))+T[3]];
end);
##
################### IdCrossedModule #########################################

#############################################################################
#0
#F	IsomorphismCrossedModules
##	Input: Two crossed modules XC, XD.
##	Output: A morphism between XC and XD if XC and XD are isomorphic.
##			return false for other else
##
InstallGlobalFunction(IsomorphismCrossedModules, function(XC,XD)
local 
	C,D,Iso,f,GC,GD,MC,MD,
	proD,emb1C,emb2C,emb1D,emb2D,
	Gens,Imgs,m,x,px,Map;
	
	C:=CatOneGroupByCrossedModule(XC);
	D:=CatOneGroupByCrossedModule(XD);
	Iso:=IsomorphismCatOneGroups(C,D);
	if Iso=fail then
		return fail;
	fi;
	f:=Iso!.mapping;
	GC:=Source(f);
	GD:=Range(f);
	
	MC:=Source(XC!.map);
	MD:=Source(XD!.map);
	
	proD:=Projection(GD);
	emb1C:=Embedding(GC,1);
	emb2C:=Embedding(GC,2);
	emb1D:=Embedding(GD,1);
	emb2D:=Embedding(GD,2);

	Gens:=GeneratorsOfGroup(MC);
	Imgs:=[];
	for m in Gens do
		x:=Image(f,Image(emb2C,m));
		px:=Image(emb1D,Image(proD,x));
		Add(Imgs,PreImagesRepresentative(emb2D,px^(-1)*x));
	od;
	
	######################################################################
	##
	Map:=function(n)
		if n=1 then
			return GroupHomomorphismByImages(MC,MD,Gens,Imgs);
		fi;
		if n=2 then
			return emb1C*f*proD;
		fi;
		Print("Only apply for n =1,2");
		return fail;
	end;
	##
	######################################################################
	return Objectify(HapCrossedModuleMorphism,
					rec(source:=XC,
						target:=XD,
						mapping:=Map
					));
end);
##
################### IsomorphismCrossedModules ###############################








