###############################################################################
#0
#F	NerveOfCatOneGroup
##	Input:	A cat-1-group or a morphism of cat-1-group or a sequence of 
##			morphisms of cat-1-groups
##	Output:	The image of the input under the functor
##				Nerve: (cat-1-groups)->(simplicial groups)
##
InstallGlobalFunction(NerveOfCatOneGroup,
function(X,n)
local	
	NerveOfCatOneGroup_Morpre, 
	C,NerveOfCatOneGroup_Obj,	
	NerveOfCatOneGroup_Mor,
	NerveOfCatOneGroup_Seq;
	  
	######################################################################
	#1	
	#F	NerveOfCatOneGroup_Obj
	##	Input: Cat-1-group C and an integer number 
	##	Output: The nerve of C
	##
	NerveOfCatOneGroup_Obj:=function(C,number)
	local 
		LGs,LBs, LDs, 
		s,t,e,
		N,M,AutM,phi,
		g,pg,m,ConjTmp,Tmp,TmpBs,TmpDs,
		Gens,GensToLists,ImgGens,
		LTmpB,LTmpD,
		i,j,k,n,len,
		EmbOnes,EmbTwos,Pros,
		ListToOne,BoundariesOfToList,DegeneraciesOfToList,
		GroupsList,BoundariesList,DegeneraciesList;

		if not IsHapCatOneGroup(C) then
			Print("This function must be applied to a cat-1-group.\n");
			return fail;
		fi;		
		s:=C!.sourceMap;
		t:=C!.targetMap;
		N:=Image(s);
		M:=Kernel(s);
		AutM:=AutomorphismGroup(M);
		e:=One(M);
		LGs:=[];
		LBs:=[];
		LDs:=[];
		EmbOnes:=[];
		EmbTwos:=[];
		Pros:=[];
		Gens:=[];
		GensToLists:=[];
		
		########## Compute the list of group G_i for i=1..n ##############
		LGs[1]:=s!.Source;
		Gens[1]:=GeneratorsOfGroup(LGs[1]);
		GensToLists[1]:=List(Gens[1],g->[g]);
		for n in [2..number] do
			ConjTmp:=[];
			len:=Length(Gens[n-1]);	
			for i in [1..len] do
				m:=GensToLists[n-1][i][1];
				for j in [2..n-1] do
					m:=m*GensToLists[n-1][i][j];
				od;
				Add(ConjTmp,ConjugatorAutomorphismNC(M,Image(t,m)));
			od;
			phi:=GroupHomomorphismByImagesNC(LGs[n-1],AutM,Gens[n-1],ConjTmp);  
			LGs[n]:=SemidirectProduct(LGs[n-1],phi,M);
			EmbOnes[n]:=Embedding(LGs[n],1);
			EmbTwos[n]:=Embedding(LGs[n],2);
			Pros[n]:=Projection(LGs[n]);	
			Gens[n]:=GeneratorsOfGroup(LGs[n]);
			len:=Length(Gens[n]);
			GensToLists[n]:=List([1..len],x->[]);
			for i in [1..len] do
				g:=Gens[n][i];
				Tmp:=[];
				for j in [1..n-1] do
					pg:=Image(Pros[n-j+1],g);
					m:=PreImagesRepresentative(EmbTwos[n-j+1],
						(Image(EmbOnes[n-j+1],pg))^(-1)*g);	
					Tmp[n-j+1]:=m;
					g:=pg;
				od;
				Tmp[1]:=g;
				GensToLists[n][i]:=Tmp;
			od;	
		od;
		#############################################################
		#2
		#F	BoundariesOfToList
		##	Input: List m:=[g_1,m_2,...,m_n]
		##	Output: List of the image of d_i(m) with i:=0..n 
		##
		BoundariesOfToList:=function(Lm,n)
		local i,j,TmpB,LB;
		
			if n=2 then 
				LB:=[[Image(t,Lm[1])*Lm[2]],[Lm[1]*Lm[2]],[Lm[1]]];	
			fi;

			if n>2 then
				LB:=[];
				
				########## Compute d_0 #########################
				TmpB:=[Image(t,Lm[1])*Lm[2]];
				for i in [2..n-1] do
					TmpB[i]:=Lm[i+1];
				od;
				Add(LB,TmpB);
				
				########## Compute d_1-->d_{n-1} ###############
				for i in [2..n] do
					TmpB:=[];
					for j in [1..i-2] do
						TmpB[j]:=Lm[j];
					od;
					TmpB[i-1]:= Lm[i-1]*Lm[i];
					for j in [i..n-1] do
						TmpB[j]:=Lm[j+1];
					od;
					Add(LB,TmpB);
				od;
				
				######### Compute d_n ##########################
				TmpB:=[];
				for i in [1..n-1] do
					TmpB[i]:=Lm[i];
				od;
				Add(LB,TmpB);
			fi;
			return LB;
		end;
		##
		########## end of BoundariesOfToList ########################
		
		#############################################################
		#2
		#F	DegeneraciesOfToList		
		##	Input: List m:=[g_1,m_2,...,m_n]
		##	Output: List of the image of s_i(m) with i:=0..n 
		##
		DegeneraciesOfToList:=function(Lm,n)
		local i,j,TmpD,LD,g;
		
			g:=Lm[1]; 
			if n=1 then 
				LD:=[[Image(s,g),Image(s,g^(-1))*g],[g,e]];	
			fi;
			if n>1 then
				LD:=[];
				
				########## Compute s_0 #########################
				TmpD:=[Image(s,g),Image(s,g^(-1))*g];
				for i in [3..n+1] do
					TmpD[i]:=Lm[i-1];
				od;
				Add(LD,TmpD);
				
				########## Compute s_1 -> s_n ##################
				for i in [2..n+1] do
					TmpD:=[];
					for j in [1..i-1] do
						TmpD[j]:=Lm[j];
					od;
					TmpD[i]:=e;
					for j in [i+1..n+1] do
						TmpD[j]:=Lm[j-1];
					od;
					Add(LD,TmpD);
				od;
			fi;
			return LD;
		end;
		##
		########## end of DegeneraciesOfToList ######################
		
		#############################################################
		#2   
		#F	ListToOne
	    ##	Input: List [g_1,m_2,m_3,...,m_n]
        ##	Output: The semi-product g_1 x| m_2 x| m_3 x| ... x| m_n
		##
		ListToOne:=function(Lm,n) 
		local i,m;
		
			if n=1 then
				m:=Lm[1];
			fi;
			if n>1 then
				m:=Lm[1];
				for i in [2..n] do    
					m:=Image(EmbOnes[i],m)*Image(EmbTwos[i],Lm[i]);
				od;	
			fi;
			return m;
		end;
		##
	    ########## end of ListToOne #################################
		
		############### Compute boundary maps #######################
		LBs:=[[t,s]];
		for n in [2..number] do
			len:=Length(Gens[n]);
			Tmp:=[];
			TmpBs:=[];
			for i in [1..len] do
				Tmp[i]:=BoundariesOfToList(GensToLists[n][i],n);
				TmpBs[i]:=List(Tmp[i],Lm->ListToOne(Lm,n-1));
			od;
			ImgGens:=[];
			for k in [1..n+1] do
				ImgGens[k]:=List([1..len],i->TmpBs[i][k]);
			od;
			LTmpB:=[];
			for k in [1..n+1] do
				LTmpB[k]:=GroupHomomorphismByImagesNC(LGs[n],LGs[n-1],
						Gens[n],ImgGens[k]);
			od;
			LBs[n]:=LTmpB;
		od;
		
		############### Compute degeneracy maps #####################
		for n in [1..number-1] do
			len:=Length(Gens[n]);
			Tmp:=[];
			TmpDs:=[];
			for i in [1..len] do
				Tmp[i]:=DegeneraciesOfToList(GensToLists[n][i],n);
				TmpDs[i]:=List(Tmp[i],Lm->ListToOne(Lm,n+1));
			od;
			ImgGens:=[];
			for k in [1..n+1] do
				ImgGens[k]:=List([1..len],i->TmpDs[i][k]);
			od;
			LTmpD:=[];
			for k in [1..n+1] do
				LTmpD[k]:=GroupHomomorphismByImagesNC(LGs[n],LGs[n+1],
						Gens[n],ImgGens[k]);
			od;
			LDs[n]:=LTmpD;	
		od;

		#############################################################
		#2
		GroupsList:=function(n)
			if n=0 then
				return N;
			fi;
			return LGs[n];
		end;
		##
		#############################################################
		
		#############################################################
		#2
		BoundariesList:=function(n,k)
			return LBs[n][k+1];
		end;
		##
		#############################################################
		
		#############################################################
		#2
		DegeneraciesList:=function(n,k)
			if n=0 and k = 0 then
				return GroupHomomorphismByFunction(N,LGs[1],
						function(x) return x; end);
			fi;
			return LDs[n][k+1];
		end;
		##
		#############################################################
		
		return Objectify(HapSimplicialGroup,
			   rec(
					groupsList:=GroupsList,
					boundariesList:=BoundariesList,
					degeneraciesList:=DegeneraciesList,
					properties:=[["length",number]]
				  ));
	end;
	##
	############### end of NerveOfCatOneGroup_Obj ########################
	
	######################################################################
	#1
	#F	NerveOfCatOneGroup_Morpre
	##	Input: Nerve of G, nerve of H, map f:G-->H
	##	Output: simplicial map between nerve of G and nerve of H
	##
	NerveOfCatOneGroup_Morpre:=function(NG,NH,f,number)
	local 
		GLs,GEmbOnes,GEmbTwos,GPros,HLs,HEmbOnes,HEmbTwos,HPros,
		Gens,GensToLists,
		i,j,n,m,g,pg,len,
		Tmp,ImgGens,Maps,
		HListToOne,Mapping;
		
		GLs:=[];
		GEmbOnes:=[];
		GEmbTwos:=[];
		GPros:=[];
		HLs:=[];
		HEmbOnes:=[];
		HEmbTwos:=[];
		HPros:=[];
		Gens:=[];
		GensToLists:=[];
		for n in [2..number] do
			GLs[n]:=NG!.groupsList(n);
			GEmbOnes[n]:=Embedding(GLs[n],1);
			GEmbTwos[n]:=Embedding(GLs[n],2);
			GPros[n]:=Projection(GLs[n]);
			HLs[n]:=NH!.groupsList(n);
			HEmbOnes[n]:=Embedding(HLs[n],1);
			HEmbTwos[n]:=Embedding(HLs[n],2);
			HPros[n]:=Projection(HLs[n]);	
			Gens[n]:=GeneratorsOfGroup(GLs[n]);
			len:=Length(Gens[n]);
			GensToLists[n]:=List([1..len],x->[]);
			for i in [1..len] do
				g:=Gens[n][i];
				Tmp:=[];
				for j in [1..n-1] do
					pg:=Image(GPros[n-j+1],g);
					m:=PreImagesRepresentative(GEmbTwos[n-j+1],(Image(
							GEmbOnes[n-j+1],pg))^(-1)*g);	
					Tmp[n-j+1]:=m;
					g:=pg;
				od;
				Tmp[1]:=g;
				GensToLists[n][i]:=Tmp;
			od;
		od;
		
		#############################################################
		#2   
		#F	HListToOne
	    ##	Input: List [h_1,m_2,m_3,...,m_n]
        ##	Output: The semi-product h_1 x| m_2 x| m_3 x| ... x| m_n
		##
		HListToOne:=function(Lm,n) 
		local i,m;
		
			m:=Lm[1];
			for i in [2..n] do    
				m:=Image(HEmbOnes[i],m)*Image(HEmbTwos[i],Lm[i]);
			od;	
			return m;
		end;
		##
		########## end of HListToOne ################################
		
		Maps:=[];
		for n in [2..number] do
			len:=Length(Gens[n]);
			ImgGens:=[];
			for i in [1..len] do
				ImgGens[i]:=HListToOne(List(GensToLists[n][i],m->Image(f,m)),n);
			od;	
			Maps[n]:=GroupHomomorphismByImages(GLs[n],HLs[n],Gens[n],ImgGens);
		od;

		#############################################################
		#2
		Mapping:=function(n)
			if n=0 then
				return GroupHomomorphismByFunction(NG!.groupsList(0),
						NH!.groupsList(0),function(x) return Image(f,x); end);
			fi;
			if n=1 then
				return f;
			fi;
			return Maps[n];
		end;
		##
		#############################################################
		
		return Objectify(HapSimplicialGroupMorphism,
		   rec(
				source:=NG,
				target:=NH,
				mapping:=Mapping,
				properties:=[["length",number]]
			  ));
	end;
	##
	############### end of NerveOfCatOneGroup_Morpre #####################
	
	######################################################################
    #1
    #F	NerveOfCatOneGroup_Mor
    ##	Input: Morphism of cat-1-groups	
    ##	Output: Simplicial map of their nerves	
	##
	NerveOfCatOneGroup_Mor:=function(Cf,n)
	local NG,NH,f;
		
		NG:=NerveOfCatOneGroup_Obj(Cf!.source,n);  
		NH:=NerveOfCatOneGroup_Obj(Cf!.target,n);  
		f:=Cf!.mapping;
		return NerveOfCatOneGroup_Morpre(NG,NH,f,n);
	end;	
	##
	############### end of NerveOfCatOneGroup_Mor ########################
	
	######################################################################
	#1
	#F	NerveOfCatOneGroup_Seq
	##	Input: Sequence of morphisms of cat1groups
	##	Output: Sequence of simplicial maps
	NerveOfCatOneGroup_Seq:=function(Lf,n)   
	local len,i,NC,Res;

		len:=Length(Lf);
		NC:=[];
		for i in [1..len] do
			NC[i]:=NerveOfCatOneGroup_Obj(Lf[i]!.source,n);
		od;
		NC[len+1]:=NerveOfCatOneGroup_Obj(Lf[len]!.target,n);
		Res:=[];
		for i in [1..len] do
			Res[i]:=NerveOfCatOneGroup_Morpre(NC[i],NC[i+1],Lf[i]!.mapping,n);
		od;
		return Res;
	end;
	##
	############### end of NerveOfCatOneGroup_Seq ########################

	if IsHapCatOneGroup(X) then
		return NerveOfCatOneGroup_Obj(X,n);
	fi;

	if IsHapCatOneGroupMorphism(X) then
		return NerveOfCatOneGroup_Mor(X,n);
	fi;
	if IsList(X) then
		if IsEmpty(X) then
			return [];
		fi;
		return NerveOfCatOneGroup_Seq(X,n);
	fi;
 end);
##
#################### end of NerveOfCatOneGroup ##############################






