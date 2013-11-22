#############################################################################
#0
#F	HomotopyLowerCentralSeriesOfCrossedModule
##	Input:	A crossed module X
##	Output:	The homotopy lower central series of X
##
InstallGlobalFunction(HomotopyLowerCentralSeriesOfCrossedModule, function(X)
local 
	del,act,M,P,GensM,ImgGensM,A,G,nat,Gs,Ps,
	nOne,XOne,i,phi,MorphismOne,
	Gens,As,nTwo,a,g,natMs,Ms,XTwo,GenMs,
	PreImGenMs,MorphismTwo,
	ActOne,MapOne,MapTwo;
	
	del:=X!.map;
	act:=X!.action;
	M:=Source(del);
	P:=Range(del);
	GensM:=GeneratorsOfGroup(M);
	ImgGensM:=List(GensM,m->Image(del,m));
	nat:=NaturalHomomorphismByNormalSubgroup(P,Image(del));
	A:=Kernel(del);
	G:=Range(nat);
	Gs:=[G];
	Ps:=[P];
	nOne:=1;
	while not IsTrivial(Gs[nOne]) do
		nOne:=nOne+1;
		Gs[nOne]:=CommutatorSubgroup(Gs[nOne-1],G);
		Ps[nOne]:=PreImage(nat,Gs[nOne]);
	od;
	MorphismOne:=[];
	if nOne>1 then	
		Ps:=Reversed(Ps);
		XOne:=[];
		for i in [1..nOne-1] do
			phi:=GroupHomomorphismByImages(M,Ps[i],GensM,ImgGensM);
			XOne[i]:=Objectify(HapCrossedModule,
							rec(map:=phi,
								action:=act
								));
		od;
		XOne[nOne]:=X;
		#########################################################
		#1
		MapOne:= function(i)
			return function(n)
				local Gens;
				if n = 1 then
					return IdentityMapping(M);
				fi;
				if n =2 then
					Gens:=GeneratorsOfGroup(Ps[i]);
					return GroupHomomorphismByImages(Ps[i],
							Ps[i+1],Gens,Gens);
				fi;
			end;
		end;
		##
		#########################################################
		
		for i in [1..nOne-1] do
			MorphismOne[i]:=Objectify(HapCrossedModuleMorphism,
						rec(source:=XOne[i],
							target:=XOne[i+1],
							mapping:=MapOne(i)
						));
		od;
	fi;  ## end of nOne>1

	G:=List(G,g->PreImagesRepresentative(nat,g));
	As:=[A];
	nTwo:=1;
	while not IsTrivial(As[nTwo]) do
		Gens:=[];
		for a in As[nTwo] do
			for g in G do
				Add(Gens,a*act(g,a^(-1)));
			od;
		od;
		nTwo:=nTwo+1;
		As[nTwo]:=Group(Gens);
	od;
	MorphismTwo:=[];
	if nTwo>1 then
		As:=Reversed(As);
		natMs:=[IdentityMapping(M)];
		Ms:=[M];
		XTwo:=[X];
		GenMs:=[GensM];
		PreImGenMs:=[GensM];
		
		#########################################################
		#1
		ActOne:=function(i)
			return	function(p,mA)
				return Image(natMs[i],act(p,
						PreImagesRepresentative(natMs[i],mA)));
			end;
		end;
		##
		#########################################################
		
		for i in [2..nTwo] do
			natMs[i]:=NaturalHomomorphismByNormalSubgroup(M,As[i]);
			Ms[i]:=Range(natMs[i]);
			GenMs[i]:=GeneratorsOfGroup(Ms[i]);
			PreImGenMs[i]:=List(GenMs[i],
					m->PreImagesRepresentative(natMs[i],m));
			phi:=GroupHomomorphismByImages(Ms[i],P,GenMs[i],
					List(PreImGenMs[i],m->Image(del,m)));
			XTwo[i]:=Objectify(HapCrossedModule,
							rec(map:=phi,
								action:=ActOne(i)
								));
		od;
		
		#########################################################
		#1
		MapTwo:=function(i)
			return function(n)
				if n = 1 then
					return GroupHomomorphismByImages(Ms[i],Ms[i+1],
						GenMs[i],List(PreImGenMs[i],
							m->Image(natMs[i+1],m)));
				fi;
				if n =2 then
					return IdentityMapping(P);
				fi;
			end;
		end;
		##
		#########################################################
		
		for i in [1..nTwo-1] do
			MorphismTwo[i]:=Objectify(HapCrossedModuleMorphism,
							rec(source:=XTwo[i],
								target:=XTwo[i+1],
								mapping:=MapTwo(i)
								));
		od;
	fi; 	##end of nTwo>1
	return Concatenation(MorphismOne,MorphismTwo);
end);
##
################### end of HomotopyLowerCentralSeriesOfCrossedModule ########


#############################################################################
#0
#F	HomotopyLowerCentralSeries
##	Input:	A crossed module or a cat-1-group
##	Output:	The homotopy lower central series of the input
##
InstallGlobalFunction(HomotopyLowerCentralSeries, function(X)
local XC,LXC;

	if IsHapCrossedModule(X) then
		return HomotopyLowerCentralSeriesOfCrossedModule(X);
	fi;	
	if IsHapCatOneGroup(X) then
		XC:=CrossedModuleByCatOneGroup(X);
		LXC:=HomotopyLowerCentralSeriesOfCrossedModule(XC);
		return CatOneGroupByCrossedModule(LXC);
	fi;
end);
##
################### end of HomotopyLowerCentralSeries #######################