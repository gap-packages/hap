#############################################################################
#0
#F	BarComplexEquivalence
##	Input:	A HAP resolution
##	Output: An equivariant chain homotopy between the bar and the HAP complex
##  
InstallGlobalFunction(BarComplexEquivalence,function(R)
local 
	e,dim,
	BarResEqui,Phi,Psi,Equiv,
	CPhi,CPsi,CEquiv;	

	e:=Identity(R!.group);
	dim:=R!.dimension;
	BarResEqui:=BarResolutionEquivalence(R);
	Phi:=BarResEqui!.phi;
	Psi:=BarResEqui!.psi;
	Equiv:=BarResEqui!.equiv;

	######################################################################
	#1
	#F	CPsi
	##	Input: A word w =[[m1,e_1],...[m_k,e_k]] with k:=R!.dimension(n)
	##	Output: The image of w under the map cpsi: cR_n->cB_n
	##	
	CPsi:=function(n,w)  
	local Rew,x,cw;
	
		cw:=StructuralCopy(w);  
		for x in cw do
			Add(x,1);
		od;
		Rew:=Psi(n,cw);
		for x in Rew do
			Remove(x,2);
		od;
		return Rew;
	end;
	##	
	############### end of CPsi ##########################################	
	
	######################################################################
	#1
	#F	CPhi
	##	Input: A word w =[[m_1,g_11,..,g_1n],...[m_k,g_k1,...,g_kn]]
	##	Output: The image of w under the map cphi: cB_n->cR_n
	##	
	CPhi:=function(n,w) 
	local Zw,x,tmp,PhiZw,i,Rew;
		
		Zw:=[];  
		for x in w do  
			tmp:=[x[1],e]; 
			for i in [2..n+1] do
				Add(tmp,x[i]);
			od;
			Add(Zw,tmp);
		od;
		PhiZw:=Phi(n,Zw);
		Rew:= List([ 1..dim(n)],x->0);
		for tmp in PhiZw do
			i:=tmp[2];
			Rew[i]:=Rew[i]+tmp[1];
		od;
	return Rew;
	end;
	##	
	############### end of CPhi ##########################################
	
	######################################################################
	#1
	#F	CEquiv
	##	Input: A word w =[[m_1,g_11,...,g_1n],...,[m_k,g_k1,...,g_kn]]
	##	Output: The image of w under the homotopy map cH_n: cB_n->cB_{n+1}
	##	
	CEquiv:=function(n,w) 
	local Zw,x,i,tmp,Rew;
	
		Zw:=[];
		for x in w do  
			tmp:=[x[1],e]; 
			for i in [2..n+1] do
				Add(tmp,x[i]);
			od;
			Add(Zw,tmp);
		od;
		Rew:=Equiv(n,Zw);
		for tmp in Rew do
			Remove(tmp,2);
		od;
		return Rew;
	end;
	##
	############### end of CEquiv ########################################

	return rec(
				phi:=CPhi,
				psi:=CPsi,
				equiv:=CEquiv
			);
end);
##				
################### end of BarComplexEquivalence ############################
	