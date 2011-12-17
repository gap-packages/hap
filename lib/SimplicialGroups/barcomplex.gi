

## Input : w =[[m1,g11,g12,g13,..,g1n],...[mk,gk1,..,gk]]
## Output: image of w under d_n
InstallGlobalFunction(BarComplexBoundary,
function(w)
	local n,i,j,sign,
	temp,iw,Dw;
	Dw:=[];
	for iw in w do 
		n:=Length(iw)-1;
	## Creat 0	 
		temp:=[iw[1]];
		for j in [2..n] do
			Add(temp,iw[j+1]);
		od;
		Add(Dw,temp);
	## Creat 1 --> n 
		sign:=-1;
		for i in [1..n-1] do
		   temp:=[sign*iw[1]];
		   for j in [1..i-1]	do
			   Add(temp,iw[j+1]);
		   od;
		   Add(temp,iw[i+1]*iw[i+2]);
		   for j in [i+2..n]do
				Add(temp,iw[j+1]);
		   od;
		   Add(Dw,temp);
		   sign:=-sign;  
		od;
	## Creat n+1
		temp:=[sign*iw[1]];
		for j in [1..n-1] do
			Add(temp,iw[j+1]);
		od;
		Add(Dw,temp);
	od;
return Dw;
end);

#################################################################################
#################################################################################
InstallGlobalFunction(BarComplexEquivalence,
function(R)

local 
	e,dim,Elts,
	BHe,Phi,Psi,Equiv,
	CPhi,CPsi,CEquiv;
	Elts:=R!.elts;
	e:=Identity(R!.group);
	dim:=R!.dimension;
	BHe:=BarResolutionEquivalence(R);
	Phi:=BHe!.phi;
	Psi:=BHe!.psi;
	Equiv:=BHe!.equiv;
#################################################	
CPhi:=function(n,w)  ###w:=[[m1,g11..,g1n],[mk,gk1,..gkn]]
	local ew,iw,temp,phiw,i,Zw;
	ew:=[];
	for iw in w do  
		temp:=[iw[1],e]; #### add indentity element
		for i in [2..n+1] do
			Add(temp,iw[i]);
		od;
		Add(ew,temp);
	od;
    phiw:=Phi(n,ew);
    Zw:= 0 * [ 1 ..dim(n)];
	for temp in phiw do
		i:=temp[2];
		Zw[i]:=Zw[i]+temp[1];
    od;
return Zw;
end;	
#######################################################
CPsi:=function(n,w)  ###w:=[[m1,e1],..,[m1,ek]] with k:=dim(n); ex w:=[[-2,1],[4,2]...]
local 
    psiw,temp;
	for temp in w do
		Add(temp,1);
	od;
	psiw:=Psi(n,w);
	for temp in psiw do
		Remove(temp,2);
	od;
return psiw;
end;	
########################################################		
CEquiv:=function(n,w) ######w:=[[m1,g11..,g1n],[mk,gk1,..gkn]
	local ew,iw,i,temp,equivw;
	ew:=[];
	for iw in w do  
		temp:=[iw[1],e]; #### add indentity element
		for i in [2..n+1] do
			Add(temp,iw[i]);
		od;
		Add(ew,temp);
	od;
    equivw:=Equiv(n,ew);
	for temp in equivw do
		Remove(temp,2);
	od;
	return equivw;
end;
#############################################################
return rec(
            phi:=CPhi,
			psi:=CPsi,
			equiv:=CEquiv
          );
end);
	