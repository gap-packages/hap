InstallGlobalFunction(SL2ZResolution,
function(m,n)
local l,p,k,
      C,R,T,RH,RK,RGamma,H,K,Gamma,D,G,F,RF;
l:=Factors(m);
p:=l[1];
k:=m/p;
C:=SL2ZTree(0,0);
R:=ResolutionGTree(C,n);
if m=1 then
return R;
else
RH:=SL2ZResolution(k,n);
H:=RH!.group;
#SetName(H,"H");
## Create resolution for K
RK:=ConjugatedResolution(RH,[[1,0],[0,p]]);
RK!.group:=ConjugateSL2ZGroup(H,[[1,0],[0,p]]);
#RK!.group:=SL2Z(p);
## Create resolution for Gamma
Gamma:=CongruenceSubgroup(k,p);
RGamma:=ResolutionFiniteSubgroup(RH,Gamma);
SetName(Gamma,"Gamma");
## Create tree of groups
D:=[RH,RK,RGamma];
G:=SL2Z(1/m); 
F:=TreeOfResolutionsToSL2Zcomplex(D,G); #Compute a non-free complex for SL(2,Z[1/3])
RF:=ResolutionGTree(F,n);
return RF;
fi;
end);