###########################################################
##	CatOneGroupByCrossedModule
##	XmodToHAP
##	CatOneGroupsByGroup
##  HomotopyCatOneGroup
##	NumberSmallCatOneGroups
##	SmallCatOneGroup
##	IdCatOneGroup
############################################################


#############################################################################
#0
#F	CatOneGroupByCrossedModule
##	Input:	Either a crossed module, or a morphism of crossed modules, or 
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
	##	Input:	A morphism fX of crossed modules
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
	##	Input:	A sequence LfX of morphisms of crossed modules LfX
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
#F	XmodToHAP
##	Input: 	A cat-1-group from the Xmod package 
##	Output:	The cat-1-group with data type from the HAP package 
##
InstallGlobalFunction(XmodToHAP,function (X)
local G,Gens,s,t,e,se,te;

    if IsHapCatOneGroup(X)  then
        return X;
    fi;
    if IsComponentObjectRep(X)  then
        if "TailMap" in NamesOfComponents(X) and 
           "HeadMap" in NamesOfComponents(X) then
		        G:=X!.Source;
			Gens:=GeneratorsOfGroup(G);
			s:=X!.TailMap;
			t:=X!.HeadMap;
			e := X!.RangeEmbedding; # Added 21/7/2015 as 
            		se := s*e;              # suggested by Chris
            		te := t*e;              # Wensley
            return Objectify(HapCatOneGroup,rec(
			sourceMap:=GroupHomomorphismByImagesNC(G,G,Gens,
					List(Gens, g->Image(se,g))),
			targetMap:=GroupHomomorphismByImagesNC(G,G,Gens,
					List(Gens, g->Image(te,g)))
					));      
        fi;
    fi;
    Print("Argument is not a cat-1-group from the Xmod package.\n");
    return fail;
end);
##
#################### end of XmodToHAP #######################################

#############################################################################
#0
#F	CatOneGroupsByGroup
##	Input:	A group G
##	Output: The list of all non-isomorphic cat-1-groups with underlying group G
##
InstallGlobalFunction(CatOneGroupsByGroup, function(G)
local 
	nk,n,k,Lst,S,p,x,tmp,i,Imgs,s,h,hinv,Gens,C,ResCats,
	ClassifyPairsByOrbit,CreatePairsByAbelianGroup,
	CreatePairsByNonAbelianGroup,CreatePairsByGroup;
	
	######################################################################
	#1
	#F	ClassifyPairsByOrbit
	##	Input: A list L of pairs of group homomorphisms [s,t]:G->G
	##				and the automorphism group of G.
	##	Output: A list of non-isomorphic pairs of L
	##
	ClassifyPairsByOrbit:=function(A,Lst)
	local 
		CL,TmpCL,Lx,Res,
		ActToMap,ActToPair,
		RefineClassesUnderGroup,
		processDuplicates;
		
		if Length(Lst)<=1 then
			return Lst;
		fi;
		
		#############################################################
		#2
		ActToMap:=function(s,f)
		  return InverseGeneralMapping(f)*s*f;
		end;
		##
		#############################################################
		#2
		ActToPair:=function(p,f)
			local h;
			h:=InverseGeneralMapping(f);
			return [h*p[1]*f,h*p[2]*f];
		end;
		##
		#############################################################
		
		#############################################################
		#2
		#F	RefineClassesUnderGroup
		##
		RefineClassesUnderGroup:=function(A,Lst,Indx,attr,actAttr)
		local 
			ValAttr,i,NC,LC,Sel,Orb,T,Dict,
			S,Gens,op,qs,g,h,img,p,x,cnt;

			ValAttr:=[];
			for i in [1..Length(Indx)] do
				ValAttr[i]:=attr(Lst[Indx[i]]);
			od;
			NC:=[];
			Gens:=SmallGeneratingSet(A);
			Sel:=[1..Length(Indx)];
			while Length(Sel)>0 do
			  # orbit algorithm on attributes, regular transversal
				LC:=[Indx[Sel[1]]];
				Orb:=[ValAttr[Sel[1]]];
				Unbind(Sel[1]);
				T:=[One(A)];
				Dict:=NewDictionary(Orb[1],true);
				AddDictionary(Dict,Orb[1],1);
				S:=TrivialSubgroup(A);
				op:=1;
				qs:=Size(A);
				while op<=Length(Orb) and Size(S)<qs do
					for g in Gens do
						img:=actAttr(Orb[op],g);
						p:=LookupDictionary(Dict,img);
						if p=fail then
							Add(Orb,img);
							AddDictionary(Dict,img,Length(Orb));
							Add(T,T[op]*g);
							qs:=Size(A)/Length(Orb);
						elif Size(S)<=qs/2 then
							x:=T[op]*g/T[p];
							S:=ClosureSubgroup(S,x);
						fi;
					od;
					op:=op+1;
				od;

				# which other values are in the orbit
				for i in [2..Length(Sel)] do
					p:=LookupDictionary(Dict,ValAttr[Sel[i]]);
					if p<>fail then
						
						x:=Indx[Sel[i]];
						AddSet(LC,x);
						Unbind(Sel[i]);
						if p>1 then # not identity
							h:=T[p]^-1; 
							Lst[x]:=ActToPair(Lst[x],h); 
						fi;
					fi;
				od;
				Add(NC,[S,LC]); 
				Sel:=Set(Sel); 
			od;
			return NC;
		end;
		##
		########## end of RefineClassesUnderGroup ###################

		#############################################################
		#2
		#P 	processDuplicates   #remove duplicates
		##
		processDuplicates:=function()
		local Lx,i,p,Sel;
		  
			TmpCL:=[];
			for Lx in CL do
				Sel:=[];
				for i in Lx[2] do
					p:=First(Sel,x->Lst[i]=Lst[x]);
					if p=fail then
						Add(Sel,i);
					fi;
				od;
				Add(TmpCL,[Lx[1],Sel]);
			od;
			CL:=TmpCL;
		end;
		##
		########## end of processDuplicates #########################

		############### Images of first component ###################
		CL:=RefineClassesUnderGroup(A,Lst,[1..Length(Lst)],x->Image(x[1]),
					function(s,a) return Image(a,s);end);
		processDuplicates();
		
		Res:=[];
		
		############### Kernels of first component ##################
		TmpCL:=[];
		for Lx in CL do
			if Length(Lx[2])= 1 then
				Add(Res,Lst[Lx[2][1]]);
			else
				Append(TmpCL,RefineClassesUnderGroup(Lx[1],Lst,Lx[2],
					x->Kernel(x[1]),function(s,a) return Image(a,s);end));
			fi;
		od;
		CL:=TmpCL;
		processDuplicates();

		############### Kernels of second component #################
		TmpCL:=[];
		for Lx in CL do
			if Length(Lx[2])= 1 then
				Add(Res,Lst[Lx[2][1]]);
			else
				Append(TmpCL,RefineClassesUnderGroup(Lx[1],Lst,Lx[2],
					x->Kernel(x[2]),function(s,a) return Image(a,s);end));
			fi;
		od;
		CL:=TmpCL;
		processDuplicates();

		############### First component #############################
		TmpCL:=[];
		for Lx in CL do
			if Length(Lx[2])= 1 then
				Add(Res,Lst[Lx[2][1]]);
			else
				Append(TmpCL,RefineClassesUnderGroup(Lx[1],Lst,Lx[2],
					x->x[1],ActToMap));
			fi;
		od;
		CL:=TmpCL;
		processDuplicates();

		############### Second component ############################
		TmpCL:=[];
		for Lx in CL do
			if Length(Lx[2])= 1 then
				Add(Res,Lst[Lx[2][1]]);
			else
				Append(TmpCL,RefineClassesUnderGroup(Lx[1],Lst,Lx[2],
					x->x[2],ActToMap));
			fi;
		od;
		CL:=TmpCL;
		processDuplicates();
		
		if IsEmpty(CL) then
			return Res;
		fi;

		if ForAny(CL,x->Length(x[2])>1) then
			Error("Uniqueness failure");
		fi;
		return Concatenation(Res,List(CL,x->Lst[x[2][1]]));
	end;
	##
	################ end of ClassifyPairsByOrbit #########################

	######################################################################
	#1
	#F	CreatePairsByAbelianGroup
	##	Input: 	An abelian group G
	##	Output: A list of non-isomorphic pairs [s,t]:G->G such that (G,s,t) 
	##			is a cat-1-group
	##
	CreatePairsByAbelianGroup:=function(G)
	local
		AbIn,NumX,SizeX,nX,GroupX,GensX,
		i,LSX,e,Gens,sum,
		GensLK,LK,sLK,LComX,tmp,GensK,Imgs,xK,xG,
		M,S,f,finv,nLK,Aut,n,LfK,K,CompK,N,GensN,t,
		ResPairs,FCombination;
		
		AbIn:=AbelianInvariants(G);
		
		if IsEmpty(AbIn) then
			return [IdentityMapping(G),IdentityMapping(G)];
		fi;

		NumX:=Set(AbIn);
		SizeX:=List(NumX,x->Length(Filtered(AbIn,i->i=x)));
		
		#################################################################
		##
		FCombination:=function(n)
		local 
			T,ST,tmp,Res,x;
			if n=1 then
				return List([0..SizeX[n]],x->[x]);
			fi;
			if n>1 then
				Res:=[];
				T:=FCombination(n-1);
				for x in [0..SizeX[n]] do
					ST:=StructuralCopy(T);
					for tmp in ST do
						Add(tmp,x);
					od;
					Append(Res,ST);
				od;
				return Res;
			fi;
		end;
		##
		#################################################################
		nX:=Length(NumX);
		GroupX:=[];
		GensX:=[];
		for i in [1..nX] do
			GroupX[i]:=CyclicGroup(NumX[i]);
			GensX[i]:=First(GroupX[i],g->Order(g)=NumX[i]);
		od;

		LSX:=[];
		for i in [1..nX] do
			Append(LSX,List([1..SizeX[i]],m->GroupX[i]));
		od;
		S:=DirectProduct(LSX);	
		e:=One(S);
		Gens:=[];
		sum:=0;
		for i in [1..nX] do
			Append(Gens,List([1..SizeX[i]],m->Image(Embedding(S,sum+m),GensX[i])));
			sum:=sum+SizeX[i];
		od;
		GensLK:=[];
		LK:=[];
		sLK:=[];

		LComX:=FCombination(Length(NumX));
		for tmp in LComX do
			GensK:=[];
			Imgs:=[];
			sum:=0;
			for i in [1..Length(tmp)] do
				xK:=List([1..tmp[i]],m->Gens[sum+m]);
				xG:=Concatenation(xK,List([tmp[i]+1..SizeX[i]],m->e));
				Append(GensK,xK);
				Append(Imgs,xG);
				sum:=sum+SizeX[i];
			od;
			Add(GensLK,GensK);
			if IsEmpty(GensK) then
				Add(LK,Group(e));
			else
				Add(LK,Group(GensK));
			fi;
			Add(sLK,GroupHomomorphismByImages(S,S,Gens,Imgs));
		od;

		f:=IsomorphismGroups(S,G);
		finv:=InverseGeneralMapping(f);
		LK:=List(LK,K->Image(f,K));
		sLK:=List(sLK,s->finv*s*f);
		GensLK:=List(LK,K->GeneratorsOfGroup(K));
		nLK:=Length(LK);

		Aut:=AutomorphismGroup(G);
		e:=One(G);
		ResPairs:=[];
		for n in [1..nLK] do
			LfK:=[];
			K:=LK[n];
			#CompK:=Complementclasses(G,K); #CHANGED 25/11/2018
                        CompK:=ComplementClassesRepresentatives(G,K);
			for N in CompK do
				GensN:=SmallGeneratingSet(N);
				Gens:=Concatenation(GensLK[n],GensN);
				Imgs:=Concatenation(GensLK[n],List(GensN,g->e));
				t:=GroupHomomorphismByImages(G,G,Gens,Imgs);
				Add(LfK,[sLK[n],t]);
			od;
			Append(ResPairs,ClassifyPairsByOrbit(Aut,LfK));
		od;
		return ResPairs;
	end;
	##	
	############### end of CreatePairsByAbelianGroup #####################
	
	######################################################################
	#1
	#F	CreatePairsByNonAbelianGroup
	##	Input: 	An non-abelian group G
	##	Output: A list of non-isomorphic pairs [s,t]:G->G such that (G,s,t) 
	##			is a cat-1-group
	##
	CreatePairsByNonAbelianGroup:=function(G)
	local	
		Aut,GensG,LN,nLN,IdLN,SizeLN,SetSizeLN,CLN,
		i,j,k,nCLN,LS,IdLS,nLS,
		n,L,nL,LfN,dem,N,nat,GoN,SizeGoN,IdGoN,K,GensK,f,h,
		x,y,s,Li,SKer,NotSKer,a,b,
		T,PairSKer,PairNotSKer,ResPairs;
		
		Aut:=AutomorphismGroup(G);
		GensG:=GeneratorsOfGroup(G);
		LN:=NormalSubgroups(G);
		nLN:=Length(LN);
		IdLN:=List(LN,x->IdGroup(x));
		SizeLN:=List(LN,x->Size(x));
		SetSizeLN:=Set(SizeLN);
		CLN:=List([1..Length(SetSizeLN)],x->[]);
		for i in [1..nLN] do
			Add(CLN[Position(SetSizeLN,SizeLN[i])],i);
		od;
		nCLN:=Length(CLN);
		LS:=LatticeSubgroups(G)!.conjugacyClassesSubgroups;
		if not IsMutable(LS) then
			LS:= ShallowCopy(LS);
		fi;
		LS:=List(LS,x->x[1]);
		IdLS:=List(LS,x->IdGroup(x));
		nLS:=Length(LS);
		PairSKer:=[];
		PairNotSKer:=[];
		for n in [1..nCLN] do
			L:=CLN[n];
			nL:=Length(L);
			LfN:=List([1..nLN],x->[]);
			for i in L do
				N:=LN[i];
				nat:=NaturalHomomorphismByNormalSubgroup(G,N);
				GoN:=Range(nat);
				SizeGoN:=Size(GoN);
				IdGoN:=IdGroup(GoN);
				for k in [1..nLS] do
					if IdLS[k]= IdGoN then
						K:=LS[k];
						if Size(Image(nat,K))=SizeGoN then
							GensK:=GeneratorsOfGroup(K);
							f:=GroupHomomorphismByImages(K,GoN,GensK,
									List(GensK,g->Image(nat,g)));
							h:=InverseGeneralMapping(f);
							s:=GroupHomomorphismByImages(G,G,GensG,
									List(GensG,g->Image(h,Image(nat,g))));
							Add(LfN[i],[k,s]);
						fi;
					fi;
				od;
			od;
			
			SKer:=[];
			for a in [1..nL] do
				i:=L[a];
				Li:=[];
				if  Size(CommutatorSubgroup(LN[i],LN[i]))=1 then
					Li:=List(LfN[i],x->[x[2],x[2]]);
				fi;
				Append(SKer,ClassifyPairsByOrbit(Aut,Li));
			od;
			Append(PairSKer,ClassifyPairsByOrbit(Aut,SKer));
			
			NotSKer:=[];
			for a in [2..nL] do
				i:=L[a];
				Li:=[];
				for b in [1..a-1] do
					j:=L[b];
					if  Size(CommutatorSubgroup(LN[i],LN[j]))=1 then
						for x in LfN[i] do
						for y in LfN[j] do
							if x[1]=y[1] then
								Add(Li,[x[2],y[2]]);
							fi;
						od;od;	
					fi;
				od;
				Append(NotSKer,ClassifyPairsByOrbit(Aut,Li));
			od;
			NotSKer:=ClassifyPairsByOrbit(Aut,NotSKer);
			T:=ShallowCopy(NotSKer);
			for x in NotSKer do
				Add(T,[x[2],x[1]]);
			od;
			Append(PairNotSKer,ClassifyPairsByOrbit(Aut,T));	
		od;
		ResPairs:=Concatenation(PairSKer,PairNotSKer);
		return ResPairs;
	end;
	##
	############### end CreatePairsByNonAbelianGroup #####################

	######################################################################
	#1
	CreatePairsByGroup:=function(G)
		if IsAbelian(G) then
			return CreatePairsByAbelianGroup(G);
		else
			return CreatePairsByNonAbelianGroup(G);
		fi;
	end;
	##
	######################################################################
	
    n:=Size(G);
	if n<=CATONEGROUP_DATA_SIZE then
		nk:=IdGroup(G);
		n:=nk[1];
		k:=nk[2];
		Gens:=GeneratorsOfGroup(G);
		S:=SmallGroup(n,k);
		h:=IsomorphismGroups(G,S);
		hinv:=InverseGeneralMapping(h);
		Lst:=[];
		if nk in CATONEGROUP_DATA_PERM then 
			for x in CATONEGROUP_DATA[n][k] do
				tmp:=[];
				for i in [1..2] do
					s:=GroupHomomorphismByImages(S,S,x[i][1],x[i][2]);
					Imgs:=List(Gens,g->Image(hinv,Image(s,(Image(h,g)))));
					tmp[i]:=GroupHomomorphismByImages(G,G,Gens,Imgs);
				od;
				Add(Lst,tmp);
			od;
		else
			p:=Pcgs(S);
			for x in CATONEGROUP_DATA[n][k] do
				tmp:=[];
				for i in [1..2] do
					Imgs:=List(x[i],m->PcElementByExponents(p,m));
					s:=GroupHomomorphismByImages(S,S,p,Imgs);
					Imgs:=List(Gens,g->Image(hinv,Image(s,(Image(h,g)))));
					tmp[i]:=GroupHomomorphismByImages(G,G,Gens,Imgs);		
				od;
				Add(Lst,tmp);
			od;
		fi;
	else 
		Lst:=CreatePairsByGroup(G);
	fi;
	
	ResCats:=[];
	for x in Lst do
		C:=Objectify(HapCatOneGroup,
					rec(sourceMap:=x[1],
					targetMap:=x[2]
					));
		Add(ResCats,C);
	od;
	return ResCats;
end);
##
################### end of CatOneGroupsByGroup ##############################

#############################################################################
#0
#F	HomotopyCatOneGroup
##	Input:	
##	Output:	
##
InstallGlobalFunction(HomotopyCatOneGroup, function(C)
	if Order(HomotopyGroup(C,1))*Order(HomotopyGroup(C,2))=Order(C) then
		return C;
	fi;
	return	CatOneGroupByCrossedModule(HomotopyCrossedModule(
			CrossedModuleByCatOneGroup(C)));
end);
##
################### end of NumberSmallCatOneGroups ##########################

#############################################################################
#0
#F	NumberSmallCatOneGroups
##	Input:	Integers n,k
##	Output:	The number of non-isomorphic cat-1-groups of SmallGroup(n,k)
##
InstallGlobalFunction(NumberSmallCatOneGroups, function(arg)
local n,k;
	
	n:=arg[1];
	if n >CATONEGROUP_DATA_SIZE then
			Print("This function only apply for order less than or equal to ",
					CATONEGROUP_DATA_SIZE,".\n");
			return fail;
	fi;
	if Length(arg)=1 then
		return Sum(List([1..NumberSmallGroups(n)],k->Length(CATONEGROUP_DATA[n][k])));
	fi;
	if Length(arg)=2 then
		k:=arg[2];
		if k>NumberSmallGroups(n) then
			Print("There are only ",NumberSmallGroups(n)," groups of order ",n,"\n");
			return fail;
		fi;
		return Length(CATONEGROUP_DATA[n][k]);
	fi;
	
    return fail;
end);
##
################### end of NumberSmallCatOneGroups ##########################

#############################################################################
#0
#F	SmallCatOneGroup
##	Input: 	Integer numbers n,k,i. 
##	Output: The ith cat-1-group of the list of all non-isomorphic cat-1-groups 
##			of underlying group SmallGroup(n,k).
##
InstallGlobalFunction(SmallCatOneGroup, function(arg)
local	S,m,sum,x,s,t,p,
		n,k,i;
		
	n:=arg[1];
	k:=arg[2];
	if n > CATONEGROUP_DATA_SIZE then
		Print("This function only apply for order less than or equal to ",
				CATONEGROUP_DATA_SIZE,".\n");
		return fail;
	fi;
	if Length(arg)=2 then
		m:=NumberSmallCatOneGroups(n);
		if k>m then
			Print("There are only ",m," cat-1-groups of order ",n,"\n");
			return fail;
		fi;
		sum:=0;
		m:=0;
		while sum<k do
			i:=k-sum;
			m:=m+1;
			sum:=sum+Length(CATONEGROUP_DATA[n][m]);
		od;
		k:=m;
	fi;
	if Length(arg)=3 then
		i:=arg[3];
		m:=Length(CATONEGROUP_DATA[n][k]);
		if i>m then
			Print("There are only ",m," cat-1-groups of SmallGroup(",n,",",k,")\n");
			return fail;
		fi;
	fi;
	S:=SmallGroup(n,k);
	x:=CATONEGROUP_DATA[n][k][i];
	if [n,k] in CATONEGROUP_DATA_PERM then 
		s:=GroupHomomorphismByImages(S,S,x[1][1],x[1][2]);
		t:=GroupHomomorphismByImages(S,S,x[2][1],x[2][2]);
	else
		p:=Pcgs(S);
		s:=GroupHomomorphismByImages(S,S,p,List(x[1],
				m->PcElementByExponents(p,m)));
		t:=GroupHomomorphismByImages(S,S,p,List(x[2],
				m->PcElementByExponents(p,m)));
	fi;
	return Objectify(HapCatOneGroup,
						rec(sourceMap:=s,
						targetMap:=t
					));
end);
##
################### end of SmallCatOneGroup #################################

#############################################################################
#0
#F	IdCatOneGroup
##
InstallGlobalFunction(IdCatOneGroup, function(C)

local 	
	ActToMap,ActToPair,ActToSubgroup,FindOrbit,processOrbit,
	s,t,G,nk,n,k,S,f,Lst,x,tmp,i,p,Imgs,A,xC,Ln,CLn,attr,actAttr,M;

	######################################################################
	#1
	ActToMap:=function(s,f)
	  return InverseGeneralMapping(f)*s*f;
	end;
	##
	######################################################################
	#1
	ActToPair:=function(p,f)
		local h;
		h:=InverseGeneralMapping(f);
		return [h*p[1]*f,h*p[2]*f];
	end;
	##
	######################################################################
	#1
	ActToSubgroup:=function(K,f) 
		return Image(f,K); 
	end;
	##
	#######################################################################
	
	######################################################################
	#1
	#F	FindOrbit
	##
	FindOrbit:=function(A,attr,actAttr)
	local 
		Gens,Orb,T,Dict,S,op,qs,g,p,img,h;

		Gens:=SmallGeneratingSet(A);
        Orb:=[attr(xC)];
		T:=[One(A)];
		Dict:=NewDictionary(Orb[1],true);
		AddDictionary(Dict,Orb[1],1);
		S:=TrivialSubgroup(A);
		op:=1;
		qs:=Size(A);
		while op<=Length(Orb) and Size(S)<qs do
			for g in Gens do
				img:=actAttr(Orb[op],g);
				p:=LookupDictionary(Dict,img);
				if p=fail then
					Add(Orb,img);
					AddDictionary(Dict,img,Length(Orb));
					Add(T,T[op]*g);
					qs:=Size(A)/Length(Orb);
				elif Size(S)<=qs/2 then # otherwise stabilizer cant grow
					h:=T[op]*g/T[p];
					S:=ClosureSubgroup(S,h);
				fi;
			od;
			op:=op+1;
		od;
		return [Dict,S,T];
	end;
	##
	######################################################################
	##
	processOrbit:=function()
		M:=FindOrbit(A,attr,actAttr);
		for i in Ln do
			p:=LookupDictionary(M[1],attr(Lst[i]));
			if p<>fail then
				Add(CLn,i);
				Lst[i]:=ActToPair(Lst[i],M[3][p]^-1);
			fi;
		od;
	end;
	##
	######################################################################
	s:= C!.sourceMap;
    t:= C!.targetMap;
    G:=Source(s);
	n:=Size(G);
	if n>CATONEGROUP_DATA_SIZE then
		Print("This function only apply for cat-1-groups of order less than ",
			CATONEGROUP_DATA_SIZE+1,"\n");
		return fail;
	fi;
	nk:=IdGroup(G);
	k:=nk[2];
	S:=SmallGroup(n,k);
	f:=IsomorphismGroups(G,S);
	Lst:=[];
	if nk in CATONEGROUP_DATA_PERM then 
		for x in CATONEGROUP_DATA[n][k] do
			tmp:=[];
			for i in [1..2] do
				tmp[i]:=GroupHomomorphismByImages(S,S,x[i][1],x[i][2]);
			od;
			Add(Lst,tmp);
		od;
	else
		p:=Pcgs(S);
		for x in CATONEGROUP_DATA[n][k] do
			tmp:=[];
			for i in [1..2] do
				Imgs:=List(x[i],m->PcElementByExponents(p,m));
				tmp[i]:=GroupHomomorphismByImages(S,S,p,Imgs);
			od;
			Add(Lst,tmp);
		od;
	fi;
	
	A:=AutomorphismGroup(S);
	xC:=ActToPair([s,t],f);
	
	############## Image of first component #########################
	A:=AutomorphismGroup(S);
	Ln:=[1..Length(Lst)];
	CLn:=[];
	attr:=s->Image(s[1]);
	actAttr:=ActToSubgroup;
	processOrbit();
	if Length(CLn) =1 then
		return [n,k,CLn[1]];
	fi;
	A:=M[2];
	Ln:=CLn;
	CLn:=[];

	############## Kernel of first component ########################
	attr:=s->Kernel(s[1]);
	actAttr:=ActToSubgroup;
	processOrbit();
	if Length(CLn) =1 then
		return [n,k,CLn[1]];
	fi;
	A:=M[2];
	Ln:=CLn;
	CLn:=[];
	
	############## Image of second component ########################
	attr:=s->Image(s[2]);
	actAttr:=ActToSubgroup;
	processOrbit();
	if Length(CLn) =1 then
		return [n,k,CLn[1]];
	fi;
	A:=M[2];
	Ln:=CLn;
	CLn:=[];
	
	############## Kernel of second component #######################
	attr:=s->Kernel(s[1]);
	actAttr:=ActToSubgroup;
	processOrbit();
	if Length(CLn) =1 then
		return [n,k,CLn[1]];
	fi;
	A:=M[2];
	Ln:=CLn;
	CLn:=[];
	
	############## First component ##################################
	attr:=s->s[1];
	actAttr:=ActToMap;
	processOrbit();
	if Length(CLn) =1 then
		return [n,k,CLn[1]];
	fi;
	A:=M[2];
	Ln:=CLn;
	CLn:=[];
	
	############## Second component #################################
	attr:=s->s[2];
	actAttr:=ActToMap;
	processOrbit();
	if Length(CLn) =1 then
		return [n,k,CLn[1]];
	fi;
	return fail;
end);
##
##################### end of IdCatOneGroup ##################################


#############################################################################
#0
#F	IsomorphismCatOneGroups	
##	Input: Two cat-1-groups C and D
##	Output: A morphism between C and D if C and D are isomorphic.
##			return false for other else
##
InstallGlobalFunction(IsomorphismCatOneGroups, function(C,D)
local 
	sC,tC,G,sD,tD,GD,
	xC,xD,A,attr,actAttr,M,p,f,h,
	ActToMap,ActToPair,ActToSubgroup,FindOrbit,processOrbit,
	Map,map;
	
	######################################################################
	#1
	ActToMap:=function(s,f)
	  return InverseGeneralMapping(f)*s*f;
	end;
	##
	######################################################################
	#1
	ActToPair:=function(p,f)
		local h;
		h:=InverseGeneralMapping(f);
		return [h*p[1]*f,h*p[2]*f];
	end;
	##
	######################################################################
	#1
	ActToSubgroup:=function(K,f) 
		return Image(f,K); 
	end;
	##
	######################################################################
	#1
	#F	FindOrbit
	##
	FindOrbit:=function(A,attr,actAttr)
	local 
		Gens,Orb,T,Dict,S,op,qs,g,p,img,h;

		Gens:=SmallGeneratingSet(A);
        Orb:=[attr(xC)];
		T:=[One(A)];
		Dict:=NewDictionary(Orb[1],true);
		AddDictionary(Dict,Orb[1],1);
		S:=TrivialSubgroup(A);
		op:=1;
		qs:=Size(A);
		while op<=Length(Orb) and Size(S)<qs do
			for g in Gens do
				img:=actAttr(Orb[op],g);
				p:=LookupDictionary(Dict,img);
				if p=fail then
					Add(Orb,img);
					AddDictionary(Dict,img,Length(Orb));
					Add(T,T[op]*g);
					qs:=Size(A)/Length(Orb);
				elif Size(S)<=qs/2 then # otherwise stabilizer cant grow
					h:=T[op]*g/T[p];
					S:=ClosureSubgroup(S,h);
				fi;
			od;
			op:=op+1;
		od;
		return [Dict,S,T];
	end;
	##
	######################################################################
	##
	processOrbit:=function()
	
		M:=FindOrbit(A,attr,actAttr);
		p:=LookupDictionary(M[1],attr(xD));
		if p=fail then
			return fail;
		fi;
		A:=M[2];
		if p<>1 then
			h:=M[3][p]^-1;
			f:=f*h;
			xD:=ActToPair(xD,h);
		fi;
	end;
	##
	######################################################################
	
	######################################################################
	sC:=C!.sourceMap;
	tC:=C!.targetMap;
	G:=Source(sC);
	sD:=D!.sourceMap;
	tD:=D!.targetMap;
	GD:=Source(sD);
	
	f:=IsomorphismGroups(GD,G);
	if f = fail then
		return fail;
	fi;
	if IsomorphismGroups(HomotopyGroup(C,1),HomotopyGroup(D,1))=fail then
		return fail;
	fi;
	if IsomorphismGroups(HomotopyGroup(C,2),HomotopyGroup(D,2))=fail then
		return fail;
	fi;
	if IsomorphismGroups(Kernel(tC),Kernel(tD))=fail then
		return fail;
	fi;
	if IsomorphismGroups(Kernel(sC),Kernel(sD))=fail then
		return fail;
	fi;
	if IsomorphismGroups(Image(tC),Image(tD))=fail then
		return fail;
	fi;
	if IsomorphismGroups(Image(sC),Image(sD))=fail then
		return fail;
	fi;
	######################################################################
	##
	Map:=function()
		xC:=[sC,tC]; 
		xD:=ActToPair([sD,tD],f);
		if xC=xD then
			return InverseGeneralMapping(f);
		fi;
		
		############## Image of first component #####################
		A:=AutomorphismGroup(G);
		attr:=s->Image(s[1]);
		actAttr:=ActToSubgroup;
		processOrbit();
		
		############## Kernel of first component #####################
		attr:=s->Kernel(s[1]);
		actAttr:=ActToSubgroup;
		processOrbit();
		if xC=xD then
			return InverseGeneralMapping(f);
		fi;
		
		############## Image of second component #####################
		attr:=s->Image(s[2]);
		actAttr:=ActToSubgroup;
		processOrbit();
		if xC=xD then
			return InverseGeneralMapping(f);
		fi;	
		############## Kernel of second component #####################
		attr:=s->Kernel(s[1]);
		actAttr:=ActToSubgroup;
		processOrbit();
		if xC=xD then
			return InverseGeneralMapping(f);
		fi;	
		############## First component #####################
		attr:=s->s[1];
		actAttr:=ActToMap;
		processOrbit();
		if xC=xD then
			return InverseGeneralMapping(f);
		fi;	
		############## Second component #####################
		attr:=s->s[2];
		actAttr:=ActToMap;
		processOrbit();
		if xC=xD then
			return InverseGeneralMapping(f);
		fi;
		return fail;
	end;
	##
	########################################################################
	map:=Map();
	if map=fail then 
		return fail;
	fi;
	return Objectify(HapCatOneGroupMorphism,
		   rec(
				source:= C,
				target:= D,
				mapping:= map
			  ));
end);
##
##################### IsomorphismCatOneGroups ###############################







	
