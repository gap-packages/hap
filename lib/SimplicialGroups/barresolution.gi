## Input : w =[[m1,h1,g11,g12,g13,..,g1n],...[mk,hk,gk1,...gk]]
## Output: image of w under d_n
#############################################
InstallGlobalFunction(AddWord,
function(L,iw)
	local siw,nsiw,nL,flag,i,j;
	siw:=StructuralCopy(iw);
	nsiw:=Length(siw);	
	nL:=Length(L);
	for i in [1..nL] do
	    flag:=0;
	    for j in [2..nsiw] do
			if L[i][j]<>siw[j] then
			   flag:=1;
			   break;
			fi;
		od;
		if flag=0 then
			L[i][1]:=L[i][1]+siw[1];
			if L[i][1]=0 then
				Remove(L,i);
			fi;
			return;
		fi;
	od;
	Add(L,siw);
end);
######################################################################
#######################################################################
InstallGlobalFunction(BarResolutionBoundary,    ### w:=[[m1,h1,g1,g2,...]... ]
function(n,w)
	local i,j,
	tmp,iw,Rew;
	Rew:=[];
	for iw in w do 
	#### Creat 0####	 
		tmp:=[iw[1],iw[2]*iw[3]];
		for j in [2..n] do
			Add(tmp,iw[j+2]);
		od;
		AddWord(Rew,tmp);
	#### Creat 1 --> n-1#####
		for i in [1..n-1] do
		   tmp:=[(-1)^i*iw[1],iw[2]];
		   for j in [1..i-1]	do
			   Add(tmp,iw[j+2]);
		   od;
		   Add(tmp,iw[i+2]*iw[i+3]);
		   for j in [i+2..n]do
				Add(tmp,iw[j+2]);
		   od;
		   AddWord(Rew,tmp); 
		od;
	#### Creat n #########
		tmp:=[(-1)^n*iw[1],iw[2]];
		for j in [1..n-1] do
			Add(tmp,iw[j+2]);
		od;
		AddWord(Rew,tmp);
	od;
return Rew;
end);

#######################################################################
#######################################################################
InstallGlobalFunction(BarResolutionEquivalence,
	function(R)
	local lenE,Elts,e,lenR,
	hoto,dim,bound,
	SearchPos,HotoR,BarResolutionHomotopy,
	PsiBase,BoundBase,i,j,TmpPsi,tmp1,tmp2,sign,base,g,Jtmp,
	Phi,Psi,Equiv;

	Elts:=R!.elts;
	lenE:=Length(Elts);
	e:=Identity(R!.group);
	hoto:=R!.homotopy;
	dim:=R!.dimension;
	bound:=R!.boundary;
	lenR:=EvaluateProperty(R,"length");
	
####################################################
####################################################
SearchPos:=function(g)
   local m,k;
   m:=Length(Elts);
   for k in [1..m] do
		if Elts[k]=g then
		    return k;
		fi;	
   od;
   Add(Elts,g);           #These two lines added by Graham
   return Length(Elts);   #
end;
######################################################
######################################################
HotoR:=function(n,w)  ### w:=[[m1,e1,pos1],...[mk,ek,posk]]
	local
		Rew, iw, Hw,iHw,m;
	Rew:=[];
	for iw in w do
	    m:=iw[1];
		Hw:=hoto(n,[iw[2],iw[3]]);
		for iHw in Hw do
			if iHw[1]>0 then 
				AddWord(Rew,[m,iHw[1],iHw[2]]);
			else
				AddWord(Rew,[-m,-iHw[1],iHw[2]]);
			fi;
		od;
	od;
	return Rew;
end;

###########################################################
############################################################
BarResolutionHomotopy:=function(n,w) #######Input w =[[m1,h1,g11,g12,g13,..,g1n],...[mk,hk,gk1,...gk]]
	local k,iw,Rew,tmp;
	Rew:=[];
	for iw in w do
		tmp:=[iw[1],e,iw[2]];
		for k in [1..n] do
			Add(tmp,iw[k+2]);
		od;
		AddWord(Rew,tmp);
	od;	
	return Rew;
end;
#############################################################
#############################################################
PsiBase:=List([0..lenR],x->[]);
	PsiBase[1][1]:=[[1,e]];  ##################[0+1][1]
	for i in [1..lenR] do
		for j in [1..dim(i)] do
			TmpPsi:=[];
			BoundBase:=bound(i,j); 			###ex:[[2,3],[-3,5]]
			for tmp1 in BoundBase do
				if tmp1[1]<0 then 
					sign:=-1;
					base:=-tmp1[1];
				else
					sign:=1;
					base:=tmp1[1];	
				fi;
				g:=Elts[tmp1[2]];
				Jtmp:=StructuralCopy(PsiBase[i][base]); ##i = (i-1) +1
				for tmp2 in Jtmp do
					tmp2[1]:=sign*tmp2[1];
					tmp2[2]:=g*tmp2[2];	
				od;
				Append(TmpPsi,Jtmp);
			od;
			PsiBase[i+1][j]:=BarResolutionHomotopy(i-1,TmpPsi);
		od;	
	od;	
####################################################################
####################################################################
Psi:= function(n,w)        ## w:=[[m1,e1,pos1],...,[mk,ek,posk]]	ex: [[-2,1,5],...[5,4,27]]
	local Rew,m,h,iw,u,Psiiw;	
	Rew:=[];
	for iw in w do
		m:=iw[1];
		h:=Elts[iw[3]];
		Psiiw:=StructuralCopy(PsiBase[n+1][iw[2]]);
		for u in Psiiw do
			u[1]:=m*u[1];
			u[2]:=h*u[2];
			AddWord(Rew,u);
		od;
	od;
	return Rew;
end;

######################################################################	
Phi:=function(n,w)  ###### Input w =[[m1,h1,g11,g12,g13,..,g1n],...[mk,hk,gk1,...gk]]
	local       ### Output:  [[m,order,postion],...]
	   iw,Rew,Reiw,h,u,cw;
	cw:=StructuralCopy(w);
	Rew:=[];	
	if n=0 then  		#### w:=[[m1,h1],...[mk,hk]]
	   for iw in cw do
			AddWord(Rew,[iw[1],1,SearchPos(iw[2])]);
	   od;
	   return Rew;
	fi;
	for iw in cw do  #### w:=[[m1,h1,g11,...],...[mk,hk,gk1,...]]
		h:=iw[2];
		iw[2]:=e;
		Reiw:=HotoR(n-1,Phi(n-1,BarResolutionBoundary(n,[iw])));
		for u in Reiw do
			u[3]:=SearchPos(h*Elts[u[3]]);
			AddWord(Rew,u);
		od;
	od; 
	return Rew;
end;
#########################################################################
Equiv:=function(n,w)         ### Input w =[[m1,h1,g11,g12,g13,..,g1n],...[mk,hk,gk1,...gk]
    local cw,h,iw,PsiPhiiw,HBiw,HLiw,tmp,Reiw,Rew,u,L;
    cw:=StructuralCopy(w);
    if n = 0 then 
		return [];
	fi;
	
	Rew:=[];
	for iw in cw do
		h:=iw[2];
		iw[2]:=e;
		HBiw:=Equiv(n-1,BarResolutionBoundary(n,[iw]));
		AddWord(HBiw,iw);
		for tmp in HBiw do
			tmp[1]:=-tmp[1];
		od;
		PsiPhiiw:= Psi(n,Phi(n,[iw]));
		HLiw:= Concatenation(PsiPhiiw,HBiw);
		Reiw:=BarResolutionHomotopy(n,HLiw);
		for u in Reiw do
			u[2]:=h*u[2];
			AddWord(Rew,u);
		od;	
	od;
    return Rew;
end;
######################################################################
######################################################################
return rec(
            phi:=Phi,
			psi:=Psi,
			equiv:=Equiv
          );
end);	
				
#############################################
#############################################

