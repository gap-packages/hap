##################################################################################
InstallGlobalFunction(ChainComplexOfSimplicialGroup,
	function(G)
	local 
	n,LG,LB,R,B,Rz,
	boundR,phi,psi,equiv,dim,boundCB,
	Dimen,
	i,j,k,m,psiw,ImageBase,t,count,
	SearchOrder,Dimension,ZeroVector,Boundary,
	M,ReB,ii,jj,len,u,Bound,
	d0,dd,
	Sign;
#################
##################

	n:=EvaluateProperty(G,"length");
	LG:=G!.groupsList;
	LB:=G!.boundariesList;
	######
	#R:=List([0..n],i->ResolutionFiniteGroup(LG(i),n-i));
	R:=List([0..n],i->ResolutionGenericGroup(LG(i),n-i));
	B:=List([0..n],i->BarComplexEquivalence(R[i+1]));
    Rz:=List([0..n],i->TensorWithIntegers(R[i+1]));
	boundR:= function(i,j,k)  ###input position [i,j] and order element k
		return Rz[i+1]!.boundary(j,k);
	end;
	###########################
	psi:=function(i,j,w)  ###input position [i,j] and word w
		return B[i+1]!.psi(j,w);
	end;
	############################
	phi:=function(i,j,w)  ###input position [i,j] and word w
		return B[i+1]!.phi(j,w);
	end;
	#########################
	equiv:=function(i,j,w) ###input position [i,j] and word w
		return B[i+1]!.equiv(j,w);
	end;
	###########################
	dim:=function(i,j) ###input position [i,j] 
		return Rz[i+1]!.dimension(j);
	end;	
#################################################################################	
boundCB:=function(n,w)  ##input position [n] and word w=[[m,g1..gn]..] of columB{n} ouput: w of columB{n-1}
	local Rew,sign,i,k,iw,d,len,tmp;
	Rew:=[];
	sign:=1;
	for i in [0..n] do
		d:=LB(n,i);
		for iw in w do;
			len:=Length(iw);
			tmp:=[sign*iw[1]];
			for k in [2..len] do
				Add(tmp,Image(d,iw[k]));
			od;
			AddWord(Rew,tmp);
		od;
		sign:=-sign;
	od;
	return Rew;		
end;

##################################################################################
	Dimen:=[];
	for i in [0..n] do
		count:=0;
		for j in [0..i] do;
			count:=count+dim(j,i-j);
		od;
		Dimen[i+1]:=count;
	od;
######################################################
Dimension:=function(n)
	return Dimen[n+1];
end;
##########################################################
ZeroVector:=function(n)
	return List([1..Dimension(n)],x->0);
end;
############################################################
SearchOrder:=function(n,k)  ### n: element Chaincomplex of Totalcomplex, k is order element basis, output: [i,j]
	local count,i,ord; 
	if k>Dimension(n) then
		return fail;
	fi;
	count:=0;
	for i in [0..n] do
		ord:=k-count;
		count:=count+dim(n-i,i);
		if k <= count then
			return [n-i,i,ord];
			break;
		fi;
	od;
end;		
##################################################################	
d0:=function(i,j,k)   ##input R(i,j) k is order of elemmet of R(i,j)
	local n,ii,jj,beg,Bound,Rew;
	n:=i+j;
	if n=0 then
		return [0];
	fi;
	Rew:=ZeroVector(n-1);   
	if j=0 then 
		return Rew;
	fi;
	beg:=0;
	for jj in [0..j-2] do
		beg:=beg+dim(n-1-jj,jj);
	od;
	Bound:=boundR(i,j,k);
	for ii in [1..dim(i,j-1)] do
		Rew[beg+ii]:=Bound[ii];
	od;
return Rew;
end;
###################################################################
ImageBase:=[];				##[i][j+1][m][k]   R_n->B_n--B-->
for i in [1..n] do 
	ImageBase[i]:=[];              
	for j in [0..n-i]do
		ImageBase[i][j+1]:=[];
		for m in [1..i] do
			ImageBase[i][j+1][m]:=[];
		od;
	od;
od;
#####Create ImageBase[i][j+1][1][k]
for i in [1..n] do
	for j in [0..n-i] do
		for k in [1..dim(i,j)] do
			psiw:=psi(i,j,[[1,k]]);	
			ImageBase[i][j+1][1][k]:=boundCB(i,psiw);
		od;
	od;
od;	
#####Creat m>1###################
for i in [2..n]do
	for j in [0..n-i]do
		for m in [2..i] do
			for k in [1..dim(i,j)] do
				psiw:=StructuralCopy(ImageBase[i][j+1][m-1][k]);
				ImageBase[i][j+1][m][k]:=boundCB(i-m+1,equiv(i-m+1,j+m-2,psiw));		
			od;
		od;
	od;
od;	

###################################################################
Sign:=function(j,m)  #the sign depend on j and m. The function return sign of d_m at (i,j) 
	local x;
		x:=j mod 2;
		m:=m+x;
		x:=Int(m/2);
		if x mod 2=1 then 
			return 1;
		else
			return -1;
		fi;
end;
####################################################
####################################################
dd:=function(i,j,m,k) 	## input d_m input R(i,j) k is order of elemmet of R(i,j) m>0
	local n,ii,jj,beg,Rew,sign,phiw;
	n:=i+j;
	if n=0 then
		return [0];
	fi;
	Rew:=ZeroVector(n-1);   
	if m>i then 
		return Rew;
	fi;
	sign:=Sign(j,m);
	phiw:= phi(i-m,j+(m-1),ImageBase[i][j+1][m][k]);
	###############
	beg:=0;
	for jj in [0..j+(m-2)] do
		beg:=beg+dim(n-1-jj,jj);
	od;
	for ii in [1..dim(i-m,j+(m-1))] do
		Rew[beg+ii]:=sign*phiw[ii];
	od;
	#if IsZero(Rew) then
	#	Print(i,j,m,k,"\n");
	#fi;
return Rew;
end;		
#######################################################################
Bound:=List([0..n],x->[]);  
Bound[1][1]:=[0];  ###(d(0,1)=0
for m in [1..n] do
	len:=Dimension(m-1);
	for t in [1..Dimension(m)] do
		M:=SearchOrder(m,t);
		i:=M[1];
		j:=M[2];
		k:=M[3];
		ReB:=d0(i,j,k);
		for ii in [1..i] do
			u:=dd(i,j,ii,k);  #compute d_ii
			for jj in [1..len] do
				ReB[jj]:=ReB[jj]+u[jj];
			od;
		od;
	Bound[m+1][t]:=ReB;
	od;
od;
##################
Boundary:=function(n,k)
	return Bound[n+1][k];
end;
				
#########################################################################
return Objectify( HapChainComplex, rec(
				boundary:=Boundary,
				dimension:=Dimension,
              properties:= [ [ "length",n],
                    [ "type", "chainComplex" ], 
                    [ "characteristic",0 ] ] ) );
end);						
