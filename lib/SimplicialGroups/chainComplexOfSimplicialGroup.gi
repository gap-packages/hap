#############################################################################
#0
#F	ChainComplexOfSimplicialGroup
##	Input:	A simplicial group or a morphism of simplicial groups or 
##				a sequence of morphisms of simplicial groups
##	Output:	The image of the input under the functor 
##				C: (simplicial groups)->(chain complexes)
##
InstallGlobalFunction(ChainComplexOfSimplicialGroup, function(X)
local 
	AddElement,
	ChainComplexOf_Objpre,
	ChainComplexOf_Obj,
	ChainComplexOf_Morpre,
	ChainComplexOf_Mor,
	ChainComplexOf_Seq;
	
	######################################################################
	#1
	#F	AddElement
	##	Input:	A list L=[[m_1,h_1,g_11,...,g_1n],...,[m_k,h_k,g_k1,...,g_kn]]
	##				and an element x=[m',h',g'_1,...,g'_n]
	##	Output:	Add the element x into the list L
	##
	AddElement:=function(L,x)
	local	sx,nx,nL,flag,i,j;

		sx:=StructuralCopy(x);
		nx:=Length(sx);	
		nL:=Length(L);
		for i in [1..nL] do
			flag:=0;
			for j in [2..nx] do
				if L[i][j]<>sx[j] then
				   flag:=1;
				   break;
				fi;
			od;
			if flag=0 then
				L[i][1]:=L[i][1]+sx[1];
				if L[i][1]=0 then
					Remove(L,i);
				fi;
				return;
			fi;
		od;
		Add(L,sx);
	end;
	##
	#################### end of AddElement  ##############################
		  
	######################################################################
	#1
	#F	ChainComplexOf_Objpre
	##
	ChainComplexOf_Objpre:=function(N,Maps,Bar,Hap)
	local 
		HapBoundary,Phi,Psi,Equiv,HapDimension,BarMap,
		Dim,BoundChain,ImageBasis,
		i,j,k,n,tmp,t,
		M,m,ii,jj,
		SearchPosition,Dimension,BelowDim,Boundary,
		d0,dm;
		
		#############################################################
		#2
		HapBoundary:= function(i,j,k)  
			return Hap[i+1]!.boundary(j,k);
		end;
		##
		#############################################################
		#2
		Psi:=function(i,j,w)  
			return Bar[i+1]!.psi(j,w);
		end;
		##
		#############################################################
		#2
		Phi:=function(i,j,w)  
			return Bar[i+1]!.phi(j,w);
		end;#
		##
		#############################################################
		#2
		Equiv:=function(i,j,w) 
			return Bar[i+1]!.equiv(j,w);
		end;
		##
		#############################################################
		#2
		HapDimension:=function(i,j) 
			return Hap[i+1]!.dimension(j);
		end;
		##
		#############################################################
		
		#############################################################
		#2
		#F	BarMap
		##	Input: A position (i,j) and a word w 
		##	Output: The image of w under del_n: B_{i}^j->B_{i-1}^j
		##
		BarMap:=function(i,j,w) 
		local Rew,sign,ii,jj,x,d,tmp;
		
			if j mod 2 = 0 then
				sign:=1;
			else
				sign:=-1;
			fi;
			Rew:=[];
			for ii in [0..i] do
				d:=Maps(i,ii);
				for x in w do
					tmp:=[sign*x[1]];
					for jj in [1..j] do
						Add(tmp,Image(d,x[jj+1]));
					od;
					AddElement(Rew,tmp);
				od;
				sign:=-sign;
			od;
			return Rew;		
		end;
		##
		########## end of BarMap ####################################

        ########## Compute the dimendion of K_n #####################
		Dim:=[];
		for i in [0..N] do
			k:=0;
			for j in [0..i] do;
				k:=k+HapDimension(j,i-j);
			od;
			Dim[i+1]:=k;
		od;
		
		#############################################################
		#2
		Dimension:=function(n)
			return Dim[n+1];
		end;
		##
		#############################################################
		
		########### Compute the sum of dimensions under the position(i,j)
		BelowDim:=List([0..N],n->[]);  ###i+1,j+1
		for i in [0..N] do
			for j in [0..N-i] do
				n:=i+j;
				tmp:=0;
				for k in [0..j-1] do
					tmp:=tmp+HapDimension(n-k,k);
				od;
				BelowDim[i+1][j+1]:=tmp;
			od;	
		od;

		############################################################
		#2
		#F	SearchPosition
		##	Input:	Numbers n and t 
		##	Output: The position (i,j) and the order of e_t 
		SearchPosition:=function(n,t) 
		local count,j,k; 

			count:=0;
			for j in [0..n] do
				k:=t-count;
				count:=count+HapDimension(n-j,j);
				if t <= count then
					return [n-j,j,k];
					break;
				fi;
			od;
		end;
		##		
		########## end of SearchPosition ############################
		
		#############################################################
		#2
		#F	d0
		##	Input:	A position (i,j) and a basis element e_k
		##	Output:	The image of e_k under the map d_0
		##
		d0:=function(i,j,k) 
		local n,t,beg,Bound,Rew;
		
			n:=i+j;
			if n=0 then
				return [0];
			fi;
			Rew:=List([1..Dimension(n-1)],x->0);   
			if j=0 then 
				return Rew;
			fi;
			beg:=BelowDim[i+1][j];		##below(i,j-1)
			Bound:=HapBoundary(i,j,k);
			for t in [1..HapDimension(i,j-1)] do
				Rew[beg+t]:=Bound[t];
			od;
		return Rew;
		end;
		##
		########## end of d0 ########################################
		
		########## Compute d_m(e_k) at the postion (i,j) ############
		
		ImageBasis:=[];				
		for i in [1..N] do 
			ImageBasis[i]:=[];              
			for j in [0..N-i]do
				ImageBasis[i][j+1]:=[];
				for m in [1..i] do
					ImageBasis[i][j+1][m]:=[];
				od;
			od;
		od;
		
		######## Compute ImageBasis[i][j+1][1][k] ################
		for i in [1..N] do
			for j in [0..N-i] do
				for k in [1..HapDimension(i,j)] do
					ImageBasis[i][j+1][1][k]:=BarMap(i,j,Psi(i,j,[[1,k]]));
				od;
			od;
		od;	
		
		######## Compute for m>1 ###############################
		for i in [2..N]do
			for j in [0..N-i]do
				for m in [2..i] do
					for k in [1..HapDimension(i,j)] do
						tmp:=StructuralCopy(ImageBasis[i][j+1][m-1][k]);
						ImageBasis[i][j+1][m][k]:=BarMap(i-m+1,j+m-1,
								Equiv(i-m+1,j+m-2,tmp));		
					od;
				od;
			od;
		od;	
		
		#############################################################
		#2
		#F	dm
		##	Input:	A position (i,j) and a basis element e_k
		##	Output:	The image of e_k under the map d_m
		##
		dm:=function(i,j,m,k) 	
		local n,t,beg,Rew,Phiw;
			
			n:=i+j;
			if n=0 then
				return [0];
			fi;
			Rew:=List([1..Dimension(n-1)],x->0);   
			if m>i then 
				return Rew;
			fi;
			Phiw:= Phi(i-m,j+(m-1),ImageBasis[i][j+1][m][k]);
			beg:=BelowDim[(i-m)+1][j+(m-1)+1];
			for t in [1..HapDimension(i-m,j+(m-1))] do
				Rew[beg+t]:=Phiw[t];
			od;
			return Rew;
		end;
		##		
		########## end of dm ########################################
		
		BoundChain:=List([0..N],x->[]);  
		BoundChain[1][1]:=[0];  
		for n in [1..N] do
			for t in [1..Dimension(n)] do
				M:=SearchPosition(n,t);
				i:=M[1];
				j:=M[2];
				k:=M[3];
				tmp:=d0(i,j,k);
				for m in [1..i] do
					tmp:=tmp+dm(i,j,m,k);  
				od;
			BoundChain[n+1][t]:=tmp;
			od;
		od;
		#############################################################
		#2
		Boundary:=function(n,k)
			return BoundChain[n+1][k];
		end;
		##
		#############################################################
						
		return Objectify( HapChainComplex, rec(
							boundary:=Boundary,
							dimension:=Dimension,
							properties:= [ [ "length",N],
								[ "type", "chainComplex" ], 
								[ "characteristic",0 ] ] 
						) );
	end;
	##
	############### end of ChainComplexOf_Objpre ##########
	
	######################################################################
	#1
	#F	ChainComplexOf_Obj
	##	Input:	A simplicial group G
	##	Output: A chain complex of G
	##
	ChainComplexOf_Obj:=function(G)
	local	N,Maps,Grps,Bar,Hap,Res;
			
		N:=EvaluateProperty(G,"length");
		Maps:=G!.boundariesList;
		Grps:=G!.groupsList;
		Res:=List([0..N],i->ResolutionGenericGroup(Grps(i),N-i));
		Bar:=List([0..N],i->BarComplexEquivalence(Res[i+1]));
		Hap:=List([0..N],i->TensorWithIntegers(Res[i+1]));
		return ChainComplexOf_Objpre(N,Maps,Bar,Hap);
	end;
	##
	############### end of ChainComplexOf_Obj ############################
	
	######################################################################
	#1
	#F	ChainComplexOf_Morpre
	##
	ChainComplexOf_Morpre:=function(N,map,MapsH,BarH,HapH,MapsG,BarG,HapG)
	local 
		PsiH,PhiH,EquivH,HapDimensionH,DimensionH,BarMapH,DimH,
		SearchPosH,ZeroVectorH,ImageBasisH,
		PsiG,PhiG,EquivG,HapDimensionG,DimensionG,BarMapG,DimG,
		SearchPosG,ZeroVectorG,ImgG,Imgs,
		i,j,k,m,n,t,tmp,Tmp,itmp,w,LM,FF,LowDimG,Mapping,
		BarMapHG;
		
		############### All functions for H #########################
		#############################################################
		#2
		PsiH:=function(i,j,w)  
			return BarH[i+1]!.psi(j,w);
		end;
		##
		#############################################################
		#2
		PhiH:=function(i,j,w)  
			return BarH[i+1]!.phi(j,w);
		end;
		##
		#############################################################
		#2
		EquivH:=function(i,j,w) 
			return BarH[i+1]!.equiv(j,w);
		end;
		##
		#############################################################
		#2
		HapDimensionH:=function(i,j) 
			return HapH[i+1]!.dimension(j);
		end;	
		##
		#############################################################

		#############################################################
		#2
		#F	BarMapH
		##	Input: A position (i,j) and a word w 
		##	Output: The image of w under del_n: BH_{i}^j->BH_{i-1}^j
		##
		BarMapH:=function(i,j,w) 
		local Rew,sign,ii,jj,x,d,tmp;
		
			if j mod 2 = 0 then
				sign:=1;
			else
				sign:=-1;
			fi;
			Rew:=[];
			for ii in [0..i] do
				d:=MapsH(i,ii);
				for x in w do
					tmp:=[sign*x[1]];
					for jj in [1..j] do
						Add(tmp,Image(d,x[jj+1]));
					od;
					AddElement(Rew,tmp);
				od;
				sign:=-sign;
			od;
			return Rew;		
		end;
		##
		########## end of BarMapH ###################################
	
		########## Compute the dimendion of KH_n ####################
		DimH:=[];
		for i in [0..N] do
			k:=0;
			for j in [0..i] do;
				k:=k+HapDimensionH(j,i-j);
			od;
			DimH[i+1]:=k;
		od;
		#############################################################
		#2
		DimensionH:=function(n)
			return DimH[n+1];
		end;
		##
		#############################################################
		
		#############################################################
		#2
		#F	SearchPosH
		##	Input:	Numbers n and t 
		##	Output: The position (i,j) and the order of e_t 
		SearchPosH:=function(n,t)  
		local count,j,k; 
		
			if t>DimensionH(n) then
				return fail;
			fi;
			count:=0;
			for j in [0..n] do
				k:=t-count;
				count:=count+HapDimensionH(n-j,j);
				if t <= count then
					return [n-j,j,k];
					break;
				fi;
			od;
		end;
		##		
		########## end of SearchPosition ############################
		
		############### All functions for G #########################
		#############################################################
		#2
		PsiG:=function(i,j,w)  
			return BarG[i+1]!.psi(j,w);
		end;
		##
		#############################################################
		#2
		PhiG:=function(i,j,w)  
			return BarG[i+1]!.phi(j,w);
		end;
		##
		#############################################################
		#2
		EquivG:=function(i,j,w) 
			return BarG[i+1]!.equiv(j,w);
		end;
		##
		#############################################################
		#2
		HapDimensionG:=function(i,j) 
			return HapG[i+1]!.dimension(j);
		end;	
		##
		#############################################################

		#############################################################
		#2
		#F	BarMapG
		##	Input: A position (i,j) and a word w 
		##	Output: The image of w under del_n: BG_{i}^j->BG_{i-1}^j
		##
		BarMapG:=function(i,j,w) 
		local Rew,sign,ii,jj,x,d,tmp;
		
			if j mod 2 = 0 then
				sign:=1;
			else
				sign:=-1;
			fi;
			Rew:=[];
			for ii in [0..i] do
				d:=MapsG(i,ii);
				for x in w do
					tmp:=[sign*x[1]];
					for jj in [1..j] do
						Add(tmp,Image(d,x[jj+1]));
					od;
					AddElement(Rew,tmp);
				od;
				sign:=-sign;
			od;
			return Rew;		
		end;
		##
		########## end of BarMapG ###################################
		
		########## Compute the dimendion of KG_n ####################
		DimG:=[];
		for i in [0..N] do
			k:=0;
			for j in [0..i] do;
				k:=k+HapDimensionG(j,i-j);
			od;
			DimG[i+1]:=k;
		od;
		#############################################################
		#2
		DimensionG:=function(n)
			return DimG[n+1];
		end;
		##
		#############################################################
		
		#############################################################
		#2
		#F	SearchPosG
		##	Input:	Numbers n and t 
		##	Output: The position (i,j) and the order of e_t 
		SearchPosG:=function(n,t)  
		local count,j,k; 
		
			if t>DimensionG(n) then
				return fail;
			fi;
			count:=0;
			for j in [0..n] do
				k:=t-count;
				count:=count+HapDimensionG(n-j,j);
				if t <= count then
					return [n-j,j,k];
					break;
				fi;
			od;
		end;
		##		
		########## end of SearchPosition ############################
		
		#############################################################
		#2
		#F	BarMapHG
		##	Input:	A position (i,j) and a word w
		##	Output:	The image of w under the map_i: BH->BG
		##
		BarMapHG:=function(i,j,w)  
		local x,f,jj,tmp,Rew;

			f:=map(i);
			Rew:=[];
			for x in w do
				tmp:=[x[1]];
				for jj in [2..j+1] do
				   Add(tmp,Image(f,x[jj]));
				od;
				AddElement(Rew,tmp);
			od;
			return Rew;
		end;	
		##
		########## end of BarMapHG ##################################
		
		########## Compute d_m(e_k) at the postion (i,j) in BH ######
		ImageBasisH:=[];				
		for i in [0..N] do 
			ImageBasisH[i+1]:=[];              
			for j in [0..N-i]do
				ImageBasisH[i+1][j+1]:=[];
				for m in [0..i] do
					ImageBasisH[i+1][j+1][m+1]:=[];
				od;
			od;
		od;
		
		######### Compute for m =1 ##################################
		for i in [0..N] do
			for j in [0..N-i] do
				for m in [0..i] do
					for k in [1..HapDimensionH(i,j)] do
						ImageBasisH[i+1][j+1][m+1][k]:=PsiH(i,j,[[1,k]]);
					od;
				od;
			od;
		od;
		
		######### Compute for m >1 ##################################
		for i in [1..N]do
			for j in [0..N-i]do
				for m in [1..i] do
					for k in [1..HapDimensionH(i,j)] do
						tmp:=StructuralCopy(ImageBasisH[i+1][j+1][m][k]);
						ImageBasisH[i+1][j+1][m+1][k]:=EquivH(i-m,j+(m-1),
								(BarMapH(i-(m-1),j+(m-1),tmp)));		
					od;
				od;
			od;
		od;	
		
		########## Compute d_m(e_k) at the postion (i,j) from BH ->BG
		ImgG:=[];				
		for i in [0..N] do 
			ImgG[i+1]:=[];              
			for j in [0..N-i]do
				ImgG[i+1][j+1]:=[];
				for m in [0..i] do
					ImgG[i+1][j+1][m+1]:=[];
				od;
			od;
		od;
		
		for i in [0..N]do
			for j in [0..N-i]do
				for m in [0..i] do
					for k in [1..HapDimensionH(i,j)] do	
					ImgG[i+1][j+1][m+1][k]:=BarMapHG(i-m,j+m,
							ImageBasisH[i+1][j+1][m+1][k]);
					od;
				od;
			od;
		od;		

		Imgs:=[];				
		for i in [0..N] do 
			Imgs[i+1]:=[];              
			for j in [0..N-i]do
				Imgs[i+1][j+1]:=[];
				for m in [0..i] do
					Imgs[i+1][j+1][m+1]:=[];
				od;
			od;
		od;
			
		for i in [0..N]do
			for j in [0..N-i]do
				for k in [1..HapDimensionH(i,j)] do
					Tmp:=List([0..i],m->[]);
					for m in [0..i] do   ##Create Tmp(m,0)
						Tmp[m+1][0+1]:=ImgG[i+1][j+1][m+1][k];
					od;
					for m in [0..i-1] do
						for n in [1..i-m] do
							w:=StructuralCopy(Tmp[m+1][n]); 
							Tmp[m+1][n+1]:=BarMapG(i-(m+(n-1)),j+(m+n),
									EquivG(i-(m+(n-1)),j+(m+(n-1)),w));
						od;
					od;
					for m in [0..i] do
						itmp:=[];
						for n in [0..m] do
							Append(itmp,Tmp[n+1][m-n+1]);
						od;
						Imgs[i+1][j+1][m+1][k]:=PhiG(i-m,j+m,itmp);
					od;
				od;
			od;
		od;		
		
		########### Compute the sum of dimensions under position(i,j)
		LowDimG:=List([0..N],n->[]);
		for i in [0..N] do
			for j in [0..N-i] do
				n:=i+j;
				tmp:=0;
				for k in [0..j-1] do
					tmp:=tmp+HapDimensionG(n-k,k);
				od;
				LowDimG[i+1][j+1]:=tmp;
			od;	
		od;		
		####################################
		FF:=List([0..N],n->[]);
		for n in [0..N] do
			for t in [1..DimensionH(n)] do
				LM:=SearchPosH(n,t);
				i:=LM[1];
				j:=LM[2];
				k:=LM[3];		
				Tmp:=List([1..LowDimG[i+1][j+1]],m->0);
				for m in [0..i] do
					Append(Tmp,Imgs[i+1][j+1][m+1][k]);
				od;
				FF[n+1][t]:=Tmp;			
			od;
		od;

		#############################################################
		#2
		#F	Mapping
		##
		Mapping:=function(v,n)
		local Rew,len,k;
			Rew:=List([1..DimensionG(n)],x->0);
			len:=Length(v);
			for k in [1..len] do;
				if v[k] <> 0 then
					Rew:=Rew+v[k]*FF[n+1][k];
				fi;
			od;	
			return Rew;
		end;
		##
		########## end of Mapping ###################################
		
	return Mapping;
	end;
	##
	############## end of ChainComplexOf_Morpre ##########################
	
	######################################################################
	#1
	#F	ChainComplexOf_Mor
	##	Input:	A simplicial map Sf:G->G'
	##	Output:	A chain complex of Sf
	##
	ChainComplexOf_Mor:=function(Sf)
	local 
		map,N,
		H,GrpsH,MapsH,RH,HapH,BarH,
		G,GrpsG,MapsG,RG,HapG,BarG;
		
		H:=Sf!.source;
		G:=Sf!.target;
		map:=Sf!.mapping;
		N:=EvaluateProperty(Sf,"length");

		GrpsH:=H!.groupsList;
		MapsH:=H!.boundariesList;
		RH:=List([0..N],i->ResolutionGenericGroup(GrpsH(i),(N+1)-i));
		HapH:=List([0..N],i->TensorWithIntegers(RH[i+1])); 
		BarH:=List([0..N],i->BarComplexEquivalence(RH[i+1]));
		
		GrpsG:=G!.groupsList;
		MapsG:=G!.boundariesList;	
		RG:=List([0..N],i->ResolutionGenericGroup(GrpsG(i),(N+1)-i));
		HapG:=List([0..N],i->TensorWithIntegers(RG[i+1])); 
		BarG:=List([0..N],i->BarComplexEquivalence(RG[i+1]));
		
		return Objectify( HapChainMap, rec(
					source:=ChainComplexOf_Objpre(N,MapsH,BarH,HapH),
					target:=ChainComplexOf_Objpre(N,MapsG,BarG,HapG),
					mapping:=ChainComplexOf_Morpre(N,map,MapsH,BarH,
								HapH,MapsG,BarG,HapG),
					properties:= [ [ "type", "chainMap" ],
							[ "characteristic",0 ] ] 
					));
	end;
	##
	############### end of ChainComplexOf_Mor ############################
	
	######################################################################
	#1
	#F	ChainComplexOf_Seq
	##	Input:	A squence of simplicial maps L
	##	Output:	A chain complex of L
	##
	ChainComplexOf_Seq:=function(L)
	local  
		nL,k,
		LSG,G,N,Grps,Maps,R,Hap,Bar,KG,RewL,map;

		nL:=Length(L);
		LSG:=[];
		for k in [1..nL] do;
			LSG[k]:=L[k]!.source;
		od;
		LSG[nL+1]:=L[nL]!.target;
		Hap:=[];
		Bar:=[];
		Maps:=[];
		KG:=[];
		for k in [1..nL+1] do
			G:=LSG[k];
			N:=EvaluateProperty(G,"length");
			Grps:=G!.groupsList;
			Maps[k]:=G!.boundariesList;	
			R:=List([0..N],i->ResolutionGenericGroup(Grps(i),(N+1)-i));
			Hap[k]:=List([0..N],i->TensorWithIntegers(R[i+1])); 
			Bar[k]:=List([0..N],i->BarComplexEquivalence(R[i+1]));
			KG[k]:=ChainComplexOf_Objpre(N,Maps[k],Bar[k],Hap[k]);
		od;
		
		RewL:=[];
		for k in [1..nL] do
			N:=EvaluateProperty(L[k],"length");
			map:=L[k]!.mapping;
			RewL[k]:=Objectify( HapChainMap, rec(
						source := KG[k],
						target := KG[k+1],
						mapping := ChainComplexOf_Morpre(N,map,Maps[k],
								Bar[k],Hap[k],Maps[k+1],Bar[k+1],Hap[k+1]),
						properties:= [ [ "type", "chainMap" ],
							[ "characteristic",0 ] ] 
					) );
		od;
		return RewL;
	end;
	##
	############## end of ChainComplexOf_Seq #############################

	if IsHapSimplicialGroup(X) then
		return ChainComplexOf_Obj(X);
	fi;

	if IsHapSimplicialGroupMorphism(X) then
		return ChainComplexOf_Mor(X);
	fi;

	if IsList(X) then
		return ChainComplexOf_Seq(X);
	fi;
end);
##
################### end of ChainComplexOfSimplicialGroup ####################

##########################################################
##########################################################
InstallOtherMethod(ChainComplex,
"for simplicial groups",
[IsHapSimplicialGroup],
function(G)
return 
ChainComplexOfSimplicialGroup(G);
end);
##########################################################
##########################################################

