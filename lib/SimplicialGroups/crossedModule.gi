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
##	HomotopyCrossedModule
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
#F	CatOneGroupByCrossedModule
##	Input:	A crossed module or a morphism of crossed modules or 
##			a sequence of morphisms of crossed modules
##	Output:	The image of the input under the functor 
##				(crossed modules)->(cat-1-groups).
##
InstallGlobalFunction(CatOneGroupByCrossedModule, function(X)
local
	CMToCat1_Obj,
	CMToCat1_Morpre,
	CMToCat1_Mor,
	CMToCat1_Seq;
	
	######################################################################
	#1
	##	CMToCat1_Obj
	##	Input:	A crossed module X
	##	Output:	The cat-1-group corresponds to X
	##
	CMToCat1_Obj:=function(XC)
	local 
		d,act,p,m,M,P,AutM,GensM,GensP,alpha,
		G,GensG,pro,emb1,emb2,Elts,Eltt,g,pg,s,t;
	
		d:=XC!.map;
		act:=XC!.action;
		M:=Source(d);
		P:=Range(d);

		AutM:=AutomorphismGroup(M);
		GensM:=GeneratorsOfGroup(M);
		GensP:=GeneratorsOfGroup(P);
		alpha:=GroupHomomorphismByImages(P,AutM,GensP,List(GensP,p->
				GroupHomomorphismByImages(M,M,GensM,List(GensM,m->act(p,m)))));
		G:=SemidirectProduct(P,alpha,M);
		GensG:=GeneratorsOfGroup(G);
		pro:=Projection(G);
		emb1:=Embedding(G,1);
		emb2:=Embedding(G,2);
		Elts:=[];
		Eltt:=[];
		for g in GensG do
			p:=Image(pro,g);
			pg:=Image(emb1,p);
			Add(Elts,pg);
			m:=PreImagesRepresentative(emb2,pg^(-1)*g);
			Add(Eltt,Image(emb1,p*Image(d,m)));
		od;
		s:=GroupHomomorphismByImages(G,G,GensG,Elts);
		t:=GroupHomomorphismByImages(G,G,GensG,Eltt);
		return Objectify(HapCatOneGroup,
					rec(sourceMap:=s,
					targetMap:=t
					));

	end;
	##
	############### end of CMToCat1_Obj ##################################
	
	######################################################################
	#1
	##	CMToCat1_Morpre
	##	
	CMToCat1_Morpre:=function(CC,CD,map)
	local 
		GC,GensGC,proC,emb1C,emb2C,g,
		GD,fM,fP,p,m,emb1D,emb2D,ImGensGC;
	
		GC:=Source(CC!.sourceMap);
		GensGC:=GeneratorsOfGroup(GC);
		GD:=Source(CD!.sourceMap);
		fM:=map(1);
		fP:=map(2);
		
		proC:=Projection(GC);
		emb1C:=Embedding(GC,1);
		emb2C:=Embedding(GC,2);
		emb1D:=Embedding(GD,1);
		emb2D:=Embedding(GD,2);
		ImGensGC:=[];
		for g in GensGC do
			p:=Image(proC,g);
			m:=PreImagesRepresentative(emb2C,(Image(emb1C,p))^(-1)*g);
			Add(ImGensGC,Image(emb1D,Image(fP,p))*Image(emb2D,Image(fM,m)));
		od;
		return Objectify(HapCatOneGroupMorphism,
					   rec(
							source:= CC,
							target:= CD,
							mapping:= GroupHomomorphismByImages(GC,
									GD,GensGC,ImGensGC)
						  ));
	end;
	##
	############### end of CMToCat1_Morpre ###############################
	
	######################################################################
	#1
	##	CMToCat1_Mor
	##	Input:	An morphism of crossed modules fX
	##	Output:	The morphism of cat-1-groups corresponds to fX
	##
	CMToCat1_Mor:=function(fX)
	local CC,CD,map;
	
		CC:=CMToCat1_Obj(fX!.source);
		CD:=CMToCat1_Obj(fX!.target);
		map:=fX!.mapping;
		return CMToCat1_Morpre(CC,CD,map);
	end;
	##
	############### end of CMToCat1_Mor ##################################
	
	######################################################################
	#1
	##	CMToCat1_Seq
	##	Input:	An sequence of morphisms of crossed modules LfX
	##	Output:	The sequence of morphisms of cat-1-groups corresponds to fX
	##	
	CMToCat1_Seq:=function(LfX)
		local n,i,GC,Res;

		n:=Length(LfX);
		GC:=[];
		for i in [1..n] do
			GC[i]:=CMToCat1_Obj(LfX[i]!.source);
		od;
		GC[n+1]:=CMToCat1_Obj(LfX[i]!.target);
		Res:=[];
		for i in [1..n] do
			Res[i]:=CMToCat1_Morpre(GC[i],GC[i+1],LfX[i]!.mapping);
		od;
		return Res;
	end;
	##
	############### end of CMToCat1_Seq ##################################

	if IsHapCrossedModule(X) then
		return CMToCat1_Obj(X);
	fi;
	if IsHapCrossedModuleMorphism(X) then
		return CMToCat1_Mor(X);
	fi;
	if IsList(X) then
		return CMToCat1_Seq(X);
	fi;
end);
##
#################### end of CatOneGroupByCrossedModule ######################

#############################################################################
#0
#F	CrossedModuleByCatOneGroup
##	Input:	A cat-1-group or a morphism of cat-1-groups or 
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
	##	Input:	A morphism of cat-1-groups fC
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
	##	Input:	A sequence of morphisms of cat-1-groups LfC
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






