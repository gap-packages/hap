#############################################################################
#0
#F	EilenbergMacLaneSimplicialGroup
##	Input:	An abelian group A or a homomorphism of abelian groups f
##			and a positive integer n, and a positive integer dim 
##	Output:	The first 1+dim terms of K(A,n) or K(f,n) under functor 
##				K(-,n): (abelian groups)->(simplicial abelian groups)
##	
InstallGlobalFunction(EilenbergMacLaneSimplicialGroup, function(X,n,dim)
local
	EilenbergMacLane_Obj,
	EilenbergMacLane_Map,
        GroupHomomorphismByImages;
GroupHomomorphismByImages:=GroupHomomorphismByImagesNC;
	######################################################################
	#1
	#F	EilenbergMacLane_Obj
	##	Input:	
	EilenbergMacLane_Obj:= function(A,NN,nK)
	local 
		nn,zero,i,j,n,k,pos,tmp,x,y,
		GensA,ZeroGrp,CoF,CoD,MN,
		ListGensG,T,N,G,H,GensG,nG,Pro,Emb,nNumG,nNumH,
		NumberFace,NumberDegen,ImgGens,
		ZeroToZero,AToZero,ZeroToA,IdA,
		Surjection,NumOfSur,
		Faces,Degens,ListGroups,
		GroupsList,FaciesList,DegeneraciesList,
		AllSurjections,CompositeOfMaps,CoDegeneracies,CoFaces;
		
		nn:=NN-1;
		
		#################################################################
		#2
		#F	AllSurjections
		##	Input:	Two numbers m, n and m>=n
		##	Output:	All surjections from [m] to [n]
		AllSurjections:=function(m,n)
		local
			CreatCouple,Res,y,M;
			
			if m<n then
				return fail;
			fi;
			if m=0 then
				return [[[0,0]]];
			fi;
			########################################################
			#3
			#F
			##
			CreatCouple:=function(m,n)
			local 
				LA,LB,M,x,i;
				if n=0 then
					M:=[];
					for i in [0..m] do
						Add(M,[i,0]);
					od;
					return [M];
				fi;
				if m=n-1 then
					M:=[];
					for i in [0..m] do
						Add(M,[i,i]);
					od;
					return [M];
				fi;
				LA:=CreatCouple(m-1,n-1);
					for x in LA do
						Add(x,[m,n-1]);
					od;
				
				LB:=CreatCouple(m-1,n);
					for x in LB do
						Add(x,[m,n]);
					od;
				
				return Concatenation(LA,LB);		
			end;
			##
			########################################################
			Res:=CreatCouple(m-1,n);
			for y in Res do
				Add(y,[m,n]);
			od;
			return Res;
		end;

		#################################################################
		#2
		#F	CoFaces
		##	Input:	A integer number n>=0
		##	Output:	All cofaces: [n]->[n+1]
		##
		CoFaces:=function(n)
		local	i,k,M,Res;
		
			Res:=[];
			for i in [0..n+1] do
				M:=[];
				for k in [0..i-1] do
					Add(M,[k,k]);
				od;
				for k in [i..n] do
					Add(M,[k,k+1]);
				od;
				Add(Res,M);
			od;
			return Res;
		end;
		##
		########## end of CoFaces ###################################
		
		#############################################################
		#2
		#F	CoDegeneracies
		##	Input:	A integer number n>=0
		##	Output:	All codegeneracies:[n]-->[n-1]
		##
		CoDegeneracies:=function(n)
		local	i,k,M,Res;
		
			Res:=[];
			for i in [0..n-1] do
				M:=[];
				for k in [0..i-1] do
					Add(M,[k,k]);
				od;
				Add(M,[i,i]);
				for k in [i+1..n] do
					Add(M,[k,k-1]);
				od;
				Add(Res,M);
			od;
			return Res;
		end;
		##
		########## end of CoDegeneracies ############################
		
		#############################################################
		#2
		#F	CompositeOfMaps
		##	Input:	Two maps m]->[n] and [n]->>[k]
		##	Output:	The map [m]->>[k] if it exists or 0 for other else
		##
		CompositeOfMaps:=function(M,N)
		local Res,k,m,Temp,i,x,y;
			k:=nn;
			m:=Length(M)-1;
			x:=M[m+1][2];
			y:=N[x+1][2];
			if y<>k then
				return 0;
			fi;
			Res:=[];
			Temp:=[];
			for i in [0..m] do
				x:=M[i+1][2];
				y:=N[x+1][2];
				Add(Res,[i,y]);
				Add(Temp,y);
			od;
			if Length(Set(Temp))<k+1 then
				return 0;
			fi;
			return Res;
		end;
		##
		########## end of CompositeOfMaps ###########################
			
		Surjection:=[];
		NumOfSur:=[];
		for i in [nn..nK] do
			Surjection[i+1]:=AllSurjections(i,nn); 	 ##[i+1]
			NumOfSur[i+1]:=Length(Surjection[i+1]);      ##[i+1]
		od;
		
		zero:=Identity(A);
		ZeroGrp:=Group(zero);
		ListGroups:=[];
		for i in [0..nn-1] do
			ListGroups[i+1]:=ZeroGrp;               
		od;
		ListGroups[nn+1]:=A;
		for i in [nn+1..nK] do
			ListGroups[i+1]:=DirectProduct(List([1..NumOfSur[i+1]],x->A));
		od;
		
		GensA:=GeneratorsOfGroup(A);
		ZeroToZero:=GroupHomomorphismByImages(ZeroGrp,ZeroGrp,[],[]);
		AToZero:=GroupHomomorphismByImages(A,ZeroGrp,GensA,List(GensA,x->zero));
		ZeroToA:=GroupHomomorphismByImages(ZeroGrp,A,[],[]);
		IdA:=GroupHomomorphismByImages(A,A,GensA,GensA);
		
		########## Compute the faces map: K_n-->K_{n-1}##############
	
		Faces:=List([1..nK],i->[]); 
		
		########### Compute the face maps d_k^i with k<n ############
		for i in [1..nn-1] do
			for j in [0..i] do
				Faces[i][j+1]:=ZeroToZero;
			od;
		od;
		
		########### Compute the face maps d_nn^i ####################
		if nn>0 then
			for j in [0..nn] do
				Faces[nn][j+1]:=AToZero;
			od;
		fi;
		
		########### Compute the face map d_n^i ######################
		NumberFace:=[];
		for n in [1..nK] do
			NumberFace[n]:=[];
			for i in [0..n] do
				NumberFace[n][i+1]:=[];
			od;
		od;
		
		for n in [nn+1..nK] do
			CoF:=CoFaces(n-1);
			MN:=Surjection[n];
			for i in [0..n] do
				for k in [1..NumOfSur[n+1]] do
					N:=Surjection[n+1][k];
					T:=CompositeOfMaps(CoF[i+1],N);
					if T=0 then
						NumberFace[n][i+1][k]:=0;
					else
						NumberFace[n][i+1][k]:=Position(MN,T);  
					fi;
				od;
			od;
		od;
		
		for n in [nn+1..nK] do
			G:=ListGroups[n+1];
			H:=ListGroups[n];
			nNumG:=NumOfSur[n+1];
			nNumH:=NumOfSur[n];
			GensG:=GeneratorsOfGroup(G);
			nG:=Length(GensG);
			Pro:=List([1..nNumG],k->Projection(G,k));
			ListGensG:=[];
			for i in [1..nG] do
				x:=GensG[i];
				tmp:=List([1..nNumG],k->Image(Pro[k],x));
				Add(ListGensG,tmp);
			od;
	
			if n=nn+1 then      ## at position nn, there is only A
				Emb:=[IdA];
			else
				Emb:=List([1..nNumH],k->Embedding(H,k));
			fi;
			
			for i in [0..n] do
				ImgGens:=[];
				for j in [1..nG] do	
					x:=Identity(H);
					for k in [1..nNumG] do
						pos:=NumberFace[n][i+1][k];
						if pos<>0 then
							x:=x*Image(Emb[pos],ListGensG[j][k]);
						fi;
					od;
					ImgGens[j]:=x;
				od; 
				Faces[n][i+1]:=GroupHomomorphismByImages(G,H,GensG,ImgGens);
			od;
		od;				
		
		########### Compute the degeneracy map s_n^i ################
		Degens:=List([0..nK-1],i->[]); 		

		########### Compute the degeneracy maps s_k^i with k<n-1 ####
		for i in [0..nn-2] do
			for j in [0..i] do
				Degens[i+1][j+1]:=ZeroToZero;
			od;
		od;

		########## Compute the degeneracy maps at s_{n-1}^i #########
		if nn>1 then
			for j in [0..nn-1] do
				Degens[nn][j+1]:=ZeroToA;
			od;
		fi;
		
		NumberDegen:=[];  ####[n+1][i+1][k]
		for n in [0..nK-1] do
			NumberDegen[n+1]:=[];       
			for i in [0..n] do
				NumberDegen[n+1][i+1]:=[];
			od;
		od;
		
		for n in [nn..nK-1] do
			CoD:=CoDegeneracies(n+1);
			MN:=Surjection[n+2];
			for i in [0..n] do
				for k in [1..NumOfSur[n+1]] do
					N:=Surjection[n+1][k];
					T:=CompositeOfMaps(CoD[i+1],N);
					if T=0 then
						NumberDegen[n+1][i+1][k]:=0;
					else
						NumberDegen[n+1][i+1][k]:=Position(MN,T);  ##i+1
					fi;
				od;
			od;
		od;
	   ###########################################
		for n in [nn..nK-1] do
			G:=ListGroups[n+1];
			H:=ListGroups[n+2];
			nNumG:=NumOfSur[n+1];
			nNumH:=NumOfSur[n+2];
			GensG:=GeneratorsOfGroup(G);
			nG:=Length(GensG);
			
			if n=nn then   
					Pro:=[IdA];
			else
				Pro:=List([1..nNumG],k->Projection(G,k));
			fi;
			ListGensG:=[];
			for i in [1..nG] do
				x:=GensG[i];
				tmp:=List([1..nNumG],k->Image(Pro[k],x));
				Add(ListGensG,tmp);
				
			od;
			Emb:=List([1..nNumH],k->Embedding(H,k));
			
			for i in [0..n] do
				ImgGens:=[];
				for j in [1..nG] do	
					x:=Identity(H);
					for k in [1..nNumG] do
						pos:=NumberDegen[n+1][i+1][k];
						if pos<>0 then
							x:=x*Image(Emb[pos],ListGensG[j][k]);
						fi;
					od;
					ImgGens[j]:=x;
				od; 
				Degens[n+1][i+1]:=GroupHomomorphismByImages(G,H,GensG,ImgGens);
			od;
		od;						
			
		#######################################################
		GroupsList:=function(i)
			return ListGroups[i+1];
		end;
		#######################################################
		FaciesList:=function(i,j)
			return Faces[i][j+1];
		end;
		#######################################################
		DegeneraciesList:=function(i,j)
			return Degens[i+1][j+1];
		end;
		###########################################################
		return Objectify(HapSimplicialGroup,
			   rec(
					groupsList:=GroupsList,
					boundariesList:=FaciesList,
					degeneraciesList:=DegeneraciesList,
					properties:=[["length",nK]]
				  ));
	end;		
	##
	############### end of EilenbergMacLane_Obj ##########################
	
	######################################################################
	#1
	#F	EilenbergMacLane_Map
	##	Input:	A homomorphism of abelian groups f:A->B and n, nK
	##	Output:	The morphism fK:K(A,n)->K(B,n)
	##
	EilenbergMacLane_Map:=function(f,n,nK)
	local 
		A,B,KA,KB,
		Maps,Mapping,GrpKA,GrpKB,
		Gens,Pro,Emb,ImgGens,
		i,j,k,t,nGens,
		h,g;
			
	A:=f!.Source;
	B:=f!.Range;
	KA:=EilenbergMacLane_Obj(A,n,nK);
	KB:=EilenbergMacLane_Obj(B,n,nK);
	Maps:=[];
	for i in [0..n-2] do
		Maps[i+1]:=GroupHomomorphismByImages(Group(Identity(A)),Group(Identity(B)),[],[]);
	od;
	Maps[n]:=f;  	##n-1
	for i in [n..nK] do
		GrpKA:=KA!.groupsList(i);
		GrpKB:=KB!.groupsList(i);
		Gens:=GeneratorsOfGroup(GrpKA);
		Pro:=[];
		Emb:=[];
		k:=Length(GrpKA!.DirectProductInfo!.groups);	
		for j in [1..k] do
			Pro[j]:=Projection(GrpKA,j);
			Emb[j]:=Embedding(GrpKB,j);
		od;
		ImgGens:=[];
		nGens:=Length(Gens);
		for j in [1..nGens] do
			h:=Gens[j];
			g:=Identity(GrpKB);
			for t in [1..k] do
				g:=g*Image(Emb[t],Image(f,Image(Pro[t],h)));
			od;
			ImgGens[j]:=g;
		od;
		Maps[i+1]:=GroupHomomorphismByImages(GrpKA,GrpKB,Gens,ImgGens);
	od;
	###################
	Mapping:=function(i)
		return Maps[i+1];
	end;
	###################
	return Objectify(HapSimplicialGroupMorphism,
		   rec(
				source:=KA,
				target:=KB,
				mapping:=Mapping,
				properties:=[["length",nK]]
			  ));
	end;		
	##
	############### end of EilenbergMacLane_Map ##########################
	
	if IsGroup(X) then
		return EilenbergMacLane_Obj(X,n,dim);
	fi;

	if IsGroupHomomorphism(X) then
		return EilenbergMacLane_Map(X,n,dim);
	fi;	
end);
##
################### end of EilenbergMacLaneSimplicialGroup ##################
