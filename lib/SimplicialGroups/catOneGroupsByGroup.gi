#############################################################################
#0
#F	CatOneGroupsByGroup
##	Input:	A group G
##	Output: The list of all non-isomorphic cat-1-groups with underlying group G
##
InstallGlobalFunction(CatOneGroupsByGroup, function(G)
local 
	nk,n,k,Lst,S,p,x,tmp,i,Imgs,s,h,hinv,Gens,C,Res,
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
	##	Output: A list of non-isomorphic pairs [s,t]:G->G such that G,s,t) 
	##			is a cat-1-groups
	##
	CreatePairsByAbelianGroup:=function(G)
	local
		GensG,LN,IdLN,n,LfN,m,i,j,N,phi,GoN,SizeGoN,IdGoN,
		k,K,GensK,f,h,s,x,y,AutG,Lst,ResCats,T,RL,
		SetIdLN,ClLfN,p,nCl;
		
		GensG:=GeneratorsOfGroup(G);
		AutG:=AutomorphismGroup(G);
		LN:=NormalSubgroups(G);
		IdLN:=List(LN,x->IdGroup(x));
		n:=Length(LN);;
		LfN:=List([1..n],x->[]);;
		for i in [1..n] do
			N:=LN[i];
			phi:=NaturalHomomorphismByNormalSubgroup(G,N);
			GoN:=Range(phi);
			SizeGoN:=Size(GoN);
			IdGoN:=IdGroup(GoN);
			for k in [1..n] do
				if IdLN[k]= IdGoN then
					K:=LN[k];
					if Size(Image(phi,K))=SizeGoN then
						GensK:=GeneratorsOfGroup(K);
						f:=GroupHomomorphismByImages(K,GoN,GensK,
								List(GensK,g->Image(phi,g)));
						h:=InverseGeneralMapping(f);
						s:=GroupHomomorphismByImages(G,G,GensG,
								List(GensG,g->Image(h,Image(phi,g))));
						Add(LfN[i],[k,s]);
					fi;
				fi;
			od;
		od;
		
		SetIdLN:=Set(IdLN);
		ClLfN:=List([1..Length(SetIdLN)],x->[]);;
		for i in [1..n] do
			if not IsEmpty(LfN[i]) then
				Add(ClLfN[Position(SetIdLN,IdLN[i])],LfN[i]);
			fi;
		od;
		p:=Position(ClLfN,[]);
		while p<> fail do
			Remove(ClLfN,p);
			p:=Position(ClLfN,[]);
		od;
		nCl:=Length(ClLfN);
		ResCats:=List(ClLfN,x->[x[1][1][2],x[1][1][2]]);;
		for m in [1..nCl] do
			Lst:=[];
			k:=Length(ClLfN[m]);
			for i in [2..k] do
			for j in [1..i-1] do
				for x in ClLfN[m][i] do
				for y in ClLfN[m][j] do
					if x[1]=y[1] then
						Add(Lst,[x[2],y[2]]);
					fi;
				od;od;
			od;od;
			if not IsEmpty(Lst) then
				RL:=ClassifyPairsByOrbit(AutG,Lst);
				T:=ShallowCopy(RL);
				for x in RL do
					Add(T,[x[2],x[1]]);
				od;
				Append(ResCats,ClassifyPairsByOrbit(AutG,T));
			fi;
		od;
		return ResCats;
	end;
	##	
	############### end of CreatePairsByAbelianGroup #####################
	
	######################################################################
	#1
	#F	CreatePairsByNonAbelianGroup
	##	Input: 	An non-abelian group G
	##	Output: A list of non-isomorphic pairs [s,t]:G->G such that (G,s,t) 
	##			is a cat-1-groups
	##
	CreatePairsByNonAbelianGroup:=function(G)
	local 
		AutG,GensG,LN,IdLN,SizeLN,NumST,N,i,j,ij,nLN,
		LS,nLS,SizeLS,IdLS,
		LfN,phi,n,m,k,GoN,SizeGoN,IdGoN,K,GensK,f,h,
		x,y,s,SKer,NotSKer,
		Lst,ResCats,T,RL;

		AutG:=AutomorphismGroup(G);	
		GensG:=GeneratorsOfGroup(G);
		
		LN:=NormalSubgroups(G);;
		IdLN:=List(LN,x->IdGroup(x));
		SizeLN:=List(LN,x->Size(x));
		NumST:=[];
		nLN:=Length(LN);
		for i in [1..nLN] do
			for j in [i..nLN] do
				if SizeLN[i]=SizeLN[j] then
					if  Size(CommutatorSubgroup(LN[i],LN[j]))=1 then
						Add(NumST,[i,j]);
					fi;
				fi;
			od;
		od;
		LS:=LatticeSubgroups(G)!.conjugacyClassesSubgroups;
		if not IsMutable(LS) then
			LS:= ShallowCopy(LS);
		fi;
		LS:=List(LS,x->x[1]);
		nLS:=Length(LS);
		SizeLS:=List(LS,x->Size(x));
		IdLS:=List(LS,x->IdGroup(x));
		
		LfN:=List([1..nLN],x->[]);
		for i in [1..nLN] do
			N:=LN[i];
			phi:=NaturalHomomorphismByNormalSubgroup(G,N);
			GoN:=Range(phi);
			SizeGoN:=Size(GoN);
			IdGoN:=IdGroup(GoN);
			for k in [1..nLS] do
				if IdLS[k]= IdGoN then
					K:=LS[k];
					if Size(Image(phi,K))=SizeGoN then
						GensK:=GeneratorsOfGroup(K);
						f:=GroupHomomorphismByImages(K,GoN,GensK,
								List(GensK,g->Image(phi,g)));
						h:=InverseGeneralMapping(f);
						s:=GroupHomomorphismByImages(G,G,GensG,
								List(GensG,g->Image(h,Image(phi,g))));
						Add(LfN[i],[k,s]);
					fi;
				fi;
			od;
		od;
		
		SKer:=List([1..Size(G)],x->[]);
		NotSKer:=List([1..Size(G)],x->[]);;
		for ij in NumST do
			i:=ij[1];
			j:=ij[2];
			if i=j then
				for x in LfN[i] do 
					Add(SKer[SizeLN[i]],[x[2],x[2]]);
				od;
			else
				for x in LfN[i] do
				for y in LfN[j] do
					if x[1]=y[1] then
						Add(NotSKer[SizeLN[i]],[x[2],y[2]]);
					fi;
				od;od;
			fi;
		od;
		
		ResCats:=[];
		for Lst in SKer do
			if Length(Lst)>0 then
				Append(ResCats,ClassifyPairsByOrbit(AutG,Lst));
			fi;
		od;
		for Lst in NotSKer do
			if Length(Lst)>0 then
				RL:=ClassifyPairsByOrbit(AutG,Lst);
				T:=ShallowCopy(RL);
				for x in RL do
					Add(T,[x[2],x[1]]);
				od;
				Append(ResCats,ClassifyPairsByOrbit(AutG,T));
			fi;
		od;
		return ResCats;
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
	Res:=[];
	for x in Lst do
		C:=Objectify(HapCatOneGroup,
					rec(sourceMap:=x[1],
					targetMap:=x[2]
					));
		Add(Res,C);
	od;
	return Res;
end);
##
################### end of CatOneGroupsByGroup ##############################


#############################################################################
#0
#F	NumberSmallCatOneGroups
##	Input:	Integers n,k
##	Output:	The number of non-isomorphic cat-1-groups of SmallGroup(n,k)
##
InstallGlobalFunction(NumberSmallCatOneGroups, function(n,k)

	if n >CATONEGROUP_DATA_SIZE then
		Print("This function only apply for order less than or equal ",
				CATONEGROUP_DATA_SIZE,".\n");
		return fail;
	fi;
	
	if k>NumberSmallGroups(n) then
		Print("There are only ",NumberSmallGroups(n)," groups of order ",n,"\n");
		return fail;
	fi;
	return Length(CATONEGROUP_DATA[n][k]);
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
InstallGlobalFunction(SmallCatOneGroup, function(n,k,i)
local	S,m,x,s,t,p;
	
	if n > CATONEGROUP_DATA_SIZE then
		Print("This function only apply for order less than or equal ",
				CATONEGROUP_DATA_SIZE,".\n");
		return fail;
	fi;
	m:=NumberSmallCatOneGroups(n,k);
	if i>m then
		Print("There are only ",m," cat-1-groups of SmallGroup(",n,",",k,")\n");
		return fail;
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


	
