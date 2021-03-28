#############################################################################
#0
#F	BarResolutionEquivalence
##	Input:	A HAP resolution
##	Output:	An equivariant chain homotopy between the bar & the HAP resolution
##  
InstallGlobalFunction(BarResolutionEquivalence,function(R)
local 
	nElts,Elts,e,nR,n,k,
	SearchPosition,AddElement,HapHomotopy,BarHomotopy,BarBoundary,
	PsiBasis,BoundBasis,TmpPsi,tmp1,tmp2,sign,base,g,Stmp,
	Phi,Psi,Equiv;

	Elts:=R!.elts;
	nElts:=Length(Elts);
	e:=Identity(R!.group);
	nR:=EvaluateProperty(R,"length");
	
	######################################################################
	#1
	##	SearchPosition
	##	Input:	An element g of G
	##	Output:	The position of g in Elts
	##
	SearchPosition:=function(g)
	local 	n,i,pos;
	        pos:=Position(Elts,g);
                if not pos=fail then return pos;
                else Add(Elts,g); return Length(Elts);
                fi;
#Adjusted 27/03/2021
    
		#n:=Length(Elts);
		#for i in [1..n] do
		#	if Elts[i]=g then
		#		return i;
		#	fi;	
		#od;
		#Add(Elts,g);        	
		#return n+1;   	
	end;
	##
	############### end of SearchPosition ################################
	
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
	#F	BarBoundary
	##	Input:	A word w =[[m_1,h_1,g_11,...,g_1n],...,[m_k,h_k,g_k1,...,g_kn]]
	##				and number n
	##	Output:	The image of w under the boudary d_n:B_n->B_{n-1} 
	##
	BarBoundary:=function(n,w)
	local i,j,tmp,x,Rew;

		Rew:=[];
		for x in w do 
		
			############### Compute 0 #####################
			tmp:=[x[1],x[2]*x[3]];
			for j in [2..n] do
				Add(tmp,x[j+2]);
			od;
			AddElement(Rew,tmp);
			
			############### Compute 1 -> n-1 ##############
			for i in [1..n-1] do
				tmp:=[(-1)^i*x[1],x[2]];
				for j in [1..i-1]	do
					Add(tmp,x[j+2]);
				od;
				Add(tmp,x[i+2]*x[i+3]);
				for j in [i+2..n]do
					Add(tmp,x[j+2]);
				od;
				AddElement(Rew,tmp); 
			od;
			
			############### Compute n #####################
			tmp:=[(-1)^n*x[1],x[2]];
			for j in [1..n-1] do
				Add(tmp,x[j+2]);
			od;
			AddElement(Rew,tmp);
		od;
	return Rew;
	end;
	##
	################### end of BarBoundary ###############################
	
	######################################################################
	#1
	#F	HapHomotopy	
	##	Input: A word w:=[[m_1,e_1,pos_1],...,[m_k,e_k,pos_k]] and number n
	##	Output: The image of w under the boundary h_n: R_n->R_{n+1}
	##
	HapHomotopy:=function(n,w)  
	local	Rew, x, Hw,iHw,m;
	
		Rew:=[];
		for x in w do
			m:=x[1];
			Hw:=R!.homotopy(n,[x[2],x[3]]);
			for iHw in Hw do
				if iHw[1]>0 then 
					AddElement(Rew,[m,iHw[1],iHw[2]]);
				else
					AddElement(Rew,[-m,-iHw[1],iHw[2]]);
				fi;
			od;
		od;
		return Rew;
	end;
	##
	############### end of HapHomotopy ###################################
	
	######################################################################
	#1
	#F	BarHomotopy	
	##	Input:	A word w =[[m_1,h_1,g_11,...,g_1n],...,[m_k,h_k,g_k1,...,g_kn]]
	##				and number n
	##	Output: The image of w under the boundary h_n: B_n->B_{n+1}
	##
	BarHomotopy:=function(n,w) 
	local i,x,Rew,tmp;
	
		Rew:=[];
		for x in w do
			tmp:=[x[1],e,x[2]];
			for i in [1..n] do
				Add(tmp,x[i+2]);
			od;
			AddElement(Rew,tmp);
		od;	
		return Rew;
	end;
	##
	############### end of HapHomotopy ###################################
		
	########## Compute the image of basis of R under the map psi:R->B ####
	PsiBasis:=List([0..nR],x->[]);
	PsiBasis[1][1]:=[[1,e]];  					##[0+1][1]
	for n in [1..nR] do
		for k in [1..R!.dimension(n)] do
			TmpPsi:=[];
			BoundBasis:=R!.boundary(n,k); 		##ex:[[2,3],[-3,5]]
			for tmp1 in BoundBasis do
				if tmp1[1]<0 then 
					sign:=-1;
					base:=-tmp1[1];
				else
					sign:=1;
					base:=tmp1[1];	
				fi;
				g:=Elts[tmp1[2]];
				Stmp:=StructuralCopy(PsiBasis[n][base]); ##n = (n-1) +1
				for tmp2 in Stmp do
					tmp2[1]:=sign*tmp2[1];
					tmp2[2]:=g*tmp2[2];	
				od;
				Append(TmpPsi,Stmp);
			od;
			PsiBasis[n+1][k]:=BarHomotopy(n-1,TmpPsi);
		od;	
	od;	
	
	######################################################################
	#1
	#F	Psi
	##	Input: A word w:=[[m1,e1,pos1],...,[mk,ek,posk]] and number n
	##	Output: The image of w under the map psi_n: R_n->B_n
	##
	Psi:= function(n,w)        
	local Rew,m,h,x,u,Psix;	
	
		Rew:=[];
		for x in w do
			m:=x[1];
			h:=Elts[x[3]];
			Psix:=StructuralCopy(PsiBasis[n+1][x[2]]);
			for u in Psix do
				u[1]:=m*u[1];
				u[2]:=h*u[2];
				AddElement(Rew,u);
			od;
		od;
		return Rew;
	end;
	##
	################ end of Psi ##########################################
	
	######################################################################
	#1
	#F	Phi
	##	Input: A word w =[[m1,h1,g11,g12,g13,..,g1n],...,[mk,hk,gk1,...gkn]]
	##	Output: The image of w under the map phi_n: B_n->R_n
	##
	Phi:=function(n,w)  
	local	x,Rew,Rex,h,u,cw;
	
		cw:=StructuralCopy(w);
		Rew:=[];	
		if n=0 then  		
		   for x in cw do
				AddElement(Rew,[x[1],1,SearchPosition(x[2])]);
		   od;
		   return Rew;
		fi;
		for x in cw do  
			h:=x[2];
			x[2]:=e;
			Rex:=HapHomotopy(n-1,Phi(n-1,
					BarBoundary(n,[x])));
			for u in Rex do
				u[3]:=SearchPosition(h*Elts[u[3]]);
				AddElement(Rew,u);
			od;
		od; 
		return Rew;
	end;
	##
	################ end of Phi ##########################################
	
	######################################################################
	#1
	#F	Equiv
	##	Input: A word w =[[m1,h1,g11,g12,g13,..,g1n],...,[mk,hk,gk1,...gkn]]
	##	Output: The image of w under the homotopy map H_n: B_n->B_{n+1}
	##
	Equiv:=function(n,w) 
    local 
		cw,h,x,
		PsiPhix,HBx,HLx,
		tmp,Rex,Rew,u;
	
		cw:=StructuralCopy(w);
		if n = 0 then 
			return [];
		fi;	
		Rew:=[];
		for x in cw do
			h:=x[2];
			x[2]:=e;
			HBx:=Equiv(n-1,BarBoundary(n,[x]));
			AddElement(HBx,x);
			for tmp in HBx do
				tmp[1]:=-tmp[1];
			od;
			PsiPhix:= Psi(n,Phi(n,[x]));
			HLx:= Concatenation(PsiPhix,HBx);
			Rex:=BarHomotopy(n,HLx);
			for u in Rex do
				u[2]:=h*u[2];
				AddElement(Rew,u);
			od;	
		od;
		return Rew;
	end;
	##
	################ end of Equiv ########################################

	return rec(
				phi:=Phi,
				psi:=Psi,
				equiv:=Equiv,
                                barBoundary:=BarBoundary,
                                barHomotopy:=BarHomotopy
			);
end);
##				
################### end of BarResolutionEquivalence #########################


