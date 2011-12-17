##################################################################################
InstallGlobalFunction(ChainComplexOfSimplicialGroup,
	function(G)
	local n,LG,LB,R,B,Rz,
	i,j,k,m,psiw,Imagebase,t,ti,tj,count,
	boundR,phi,psi,equiv,dim,boundCB,SearchOrder,
	Dimen,Dimension,
	M,ReB,ii,jj,len,u,Bound,Boundary,
	d0,dd,
	Zerovector;

	n:=EvaluateProperty(G,"length");
	LG:=G!.groupsList;
	LB:=G!.boundariesList;
	######
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
	local Rew,sign,i,iw,d,len,temp;
	Rew:=[];
	sign:=1;
	for i in [0..n] do
		d:=LB(n,i);
		for iw in w do;
			len:=Length(iw);
			temp:=[sign*iw[1]];
			for i in [2..len] do
				Add(temp,Image(d,iw[i]));
			od;
			Add(Rew,temp);
		od;
		sign:=-sign;
	od;
	return ReduceAlg(Rew);		
end;

##################################################################################
Dimen:=[];
	for k in [0..n] do
		count:=0;
		for j in [0..k] do;
			count:=count+dim(j,k-j);
		od;
		Dimen[k+1]:=count;
	od;

Dimension:=function(n)
	return Dimen[n+1];
end;

##########################################################
Zerovector:=[];
for i in [0..n] do
	Zerovector[i+1]:=List([1..Dimension(i)],x->0);
od;
########################################################
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
	local n,jj,ii,beg,bound,Rew;
	n:=i+j;
	if n=0 then
		return [0];
	fi;
	Rew:=StructuralCopy(Zerovector[n]);   ### n-1
	if j=0 then 
		return Rew;
	fi;
	beg:=0;
	for jj in [0..j-2] do
		beg:=beg+dim(n-1-jj,jj);
	od;
	bound:=boundR(i,j,k);
	for ii in [1..dim(i,j-1)] do
		Rew[beg+ii]:=bound[ii];
	od;
return Rew;
end;
###################################################################
Imagebase:=[];				#[j][j+1][m][k]   R_n->B_n--B-->
for i in [1..n] do 
	Imagebase[i]:=[];              
	for j in [0..n-i]do
		Imagebase[i][j+1]:=[];
		for m in [1..i] do
			Imagebase[i][j+1][m]:=[];
		od;
	od;
od;
#Creat Imagebase[i][j+1][1][k]
for i in [1..n] do
	for j in [0..n-i] do
		for k in [1..dim(i,j)] do
		psiw:=psi(i,j,[[1,k]]);	
		Imagebase[i][j+1][1][k]:=boundCB(i,psiw);
		od;
	od;
od;	
##Creat m>1
for i in [2..n]do
	for j in [0..n-i]do
		for m in [2..i] do
			for k in [1..dim(i,j)] do
				psiw:=StructuralCopy(Imagebase[i][j+1][m-1][k]);
				ti:=i-m+1;
				tj:=j+m-2;
				psiw:=boundCB(ti,equiv(ti,tj,psiw));	
				Imagebase[i][j+1][m][k]:=StructuralCopy(psiw);
			od;
		od;
	od;
od;	
##############################
dd:=function(m,i,j,k) 	## input d_m input R(i,j) k is order of elemmet of R(i,j)
	local n,ii,jj,beg,Rew,sign,phiw;
	n:=i+j;
	if k > dim(i,j) then
		return fail;
	fi;
	if n=0 then
		return [0];
	fi;
	Rew:=StructuralCopy(Zerovector[n]);   ##n-1
	if m>i then 
		return Rew;
	fi;
	if j mod 2 =0 then 
		sign:=1;
	else
		sign:=-1;
	fi;
	phiw:= phi(i-m,j+m-1,Imagebase[i][j+1][m][k]);
	beg:=0;
	###############
	for jj in [0..j+(m-2)] do
		beg:=beg+dim(n-1-jj,jj);
	od;
	for ii in [1..dim(i-m,j+m-1)] do
		Rew[beg+ii]:=sign*phiw[ii];
	od;
return Rew;
end;		
#######################################################################
Bound:=List([0..n],x->[]);  
for m in [0..n] do
	for t in [1..Dimension(m)] do
		M:=SearchOrder(m,t);
		i:=M[1];
		j:=M[2];
		k:=M[3];
		ReB:=d0(i,j,k);
		len:=Length(ReB);
		for jj in [1..i] do
			u:=dd(jj,i,j,k);
			for ii in [1..len] do
				ReB[ii]:=ReB[ii]+u[ii];
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
