##################################################################################
InstallGlobalFunction(ChainComplexOfSimplicialGroup,
	function(G)
	local 
		N,Grps,Maps,R,Bar,Hap,tmp,
		BoundHap,Phi,Psi,Equiv,Dim,MapBar,
		LDim,LowDim,
		i,j,k,n,Psiw,ImageBase,t,
		SearchPos,Dimension,Boundary,
		M,ReD,m,ii,jj,DM,BoundChain,
		d0,dm;
#################
##################
	N:=EvaluateProperty(G,"length");
	Grps:=G!.groupsList;
	Maps:=G!.boundariesList;
	######
	R:=List([0..N],i->ResolutionGenericGroup(Grps(i),N-i));
	Bar:=List([0..N],i->BarComplexEquivalence(R[i+1]));
    Hap:=List([0..N],i->TensorWithIntegers(R[i+1]));
	BoundHap:= function(i,j,k)  ###input position [i,j] and order element k
		return Hap[i+1]!.boundary(j,k);
	end;
	###########################
	Psi:=function(i,j,w)  ###input position [i,j] and word w
		return Bar[i+1]!.psi(j,w);
	end;
	############################
	Phi:=function(i,j,w)  ###input position [i,j] and word w
		return Bar[i+1]!.phi(j,w);
	end;
	#########################
	Equiv:=function(i,j,w) ###input position [i,j] and word w
		return Bar[i+1]!.equiv(j,w);
	end;
	###########################
	Dim:=function(i,j) ###input position [i,j] 
		return Hap[i+1]!.dimension(j);
	end;	
#################################################################################	
MapBar:=function(i,j,w)  ##input position [n] and word w=[[m,h1..hn]..] of columB{n} ouput: w of columB{n-1}
	local Rew,sign,k,t,iw,d,n,tmp;
	if j mod 2 = 0 then
		sign:=1;
	else
		sign:=-1;
	fi;
	Rew:=[];
	n:=j+1;
	for k in [0..i] do
		d:=Maps(i,k);
		for iw in w do;
			tmp:=[sign*iw[1]];
			for t in [2..n] do
				Add(tmp,Image(d,iw[t]));
			od;
			AddWord(Rew,tmp);
		od;
		sign:=-sign;
	od;
	return Rew;		
end;

##################################################################################
	LDim:=[];
	for i in [0..N] do
		k:=0;
		for j in [0..i] do;
			k:=k+Dim(j,i-j);
		od;
		LDim[i+1]:=k;
	od;
	################  lower (i,j);
	LowDim:=List([0..N],n->[]);  ###i+1,j+1
	for i in [0..N] do
		for j in [0..N-i] do
			n:=i+j;
			tmp:=0;
			for k in [0..j-1] do
				tmp:=tmp+Dim(n-k,k);
			od;
			LowDim[i+1][j+1]:=tmp;
		od;	
	od;
######################################################
Dimension:=function(n)
	return LDim[n+1];
end;
############################################################
SearchPos:=function(n,t)  ### n: element Chaincomplex of Totalcomplex, k is order element basis, output: [i,j]
	local count,j,k; 
	if t>Dimension(n) then
		return fail;
	fi;
	count:=0;
	for j in [0..n] do
		k:=t-count;
		count:=count+Dim(n-j,j);
		if t <= count then
			return [n-j,j,k];
			break;
		fi;
	od;
end;		
##################################################################	
d0:=function(i,j,k)   ##input R(i,j) k is order of elemmet of R(i,j)
	local n,t,beg,Bound,Rew;
	n:=i+j;
	if n=0 then
		return [0];
	fi;
	Rew:=List([1..Dimension(n-1)],x->0);   
	if j=0 then 
		return Rew;
	fi;
	beg:=LowDim[i+1][j];
	Bound:=BoundHap(i,j,k);
	for t in [1..Dim(i,j-1)] do
		Rew[beg+t]:=Bound[t];
	od;
return Rew;
end;
###################################################################
ImageBase:=[];				##[i][j+1][m][k]   R_n->B_n--B-->
for i in [1..N] do 
	ImageBase[i]:=[];              
	for j in [0..N-i]do
		ImageBase[i][j+1]:=[];
		for m in [1..i] do
			ImageBase[i][j+1][m]:=[];
		od;
	od;
od;
#####Create ImageBase[i][j+1][1][k]
for i in [1..N] do
	for j in [0..N-i] do
		for k in [1..Dim(i,j)] do
			ImageBase[i][j+1][1][k]:=MapBar(i,j,Psi(i,j,[[1,k]]));
		od;
	od;
od;	
#####Creat m>1###################
for i in [2..N]do
	for j in [0..N-i]do
		for m in [2..i] do
			for k in [1..Dim(i,j)] do
				tmp:=StructuralCopy(ImageBase[i][j+1][m-1][k]);
				ImageBase[i][j+1][m][k]:=MapBar(i-m+1,j+m-1,Equiv(i-m+1,j+m-2,tmp));		
			od;
		od;
	od;
od;	
####################################################
####################################################
dm:=function(i,j,m,k) 	## input d_m input R(i,j) k is order of elemmet of R(i,j) m>0
	local n,t,beg,Rew,Phiw;
	n:=i+j;
	if n=0 then
		return [0];
	fi;
	Rew:=List([1..Dimension(n-1)],x->0);   
	if m>i then 
		return Rew;
	fi;
	Phiw:= Phi(i-m,j+(m-1),ImageBase[i][j+1][m][k]);
	beg:=LowDim[(i-m)+1][j+(m-1)+1];
	for t in [1..Dim(i-m,j+(m-1))] do
		Rew[beg+t]:=Phiw[t];
	od;
return Rew;
end;		
#######################################################################
BoundChain:=List([0..N],x->[]);  
BoundChain[1][1]:=[0];  ###(d(0,1)=0  n+1,k
for n in [1..N] do
	for t in [1..Dimension(n)] do
		M:=SearchPos(n,t);
		i:=M[1];
		j:=M[2];
		k:=M[3];
		ReD:=d0(i,j,k);
		for ii in [1..i] do
			DM:=dm(i,j,ii,k);  #compute d_ii
			ReD:=ReD+DM;
		od;
	BoundChain[n+1][t]:=ReD;
	od;
od;
##################
Boundary:=function(n,k)
	return BoundChain[n+1][k];
end;
				
#########################################################################
return Objectify( HapChainComplex, rec(
				boundary:=Boundary,
				dimension:=Dimension,
              properties:= [ [ "length",N],
                    [ "type", "chainComplex" ], 
                    [ "characteristic",0 ] ] ) );
end);						
