InstallGlobalFunction(SL2ZResolution_alt,
function(arg)
local l,p,k,
      m,n,tietze,C,R,T,RH,RK,RGamma,H,K,Gamma,D,G,F,RF;
m:=arg[1];
n:=arg[2];
if IsBound(arg[3]) then tietze:=arg[3]; else tietze:=0; fi;
l:=Factors(m);
#p:=l[1];
p:=l[Length(l)];
k:=m/p;
#C:=SL2ZTree(0,0);
#R:=ResolutionGTree(C,n);
R:=ResolutionSL2Z(1,n);
if m=1 then
return R;
else
RH:=SL2ZResolution_alt(k,n);
H:=RH!.group;
#SetName(H,"H");
## Create resolution for K
RK:=ConjugatedResolution(RH,[[1,0],[0,p]]);
RK!.group:=ConjugateSL2ZGroup(H,[[1,0],[0,p]]);
#RK!.group:=SL2Z(p);
## Create resolution for Gamma
Gamma:=CongruenceSubgroup(k,p);
RGamma:=ResolutionFiniteSubgroup(RH,Gamma);
if tietze>0 then
RGamma:=HAPTietzeReduction_Inf(RGamma,tietze);
fi;
SetName(Gamma,"Gamma");
## Create tree of groups
D:=[RH,RK,RGamma];
G:=SL2Z(1/m); 
F:=TreeOfResolutionsToSL2Zcomplex(D,G); #Compute a non-free complex for SL(2,Z[1/3])
if tietze >0 then
RF:=FreeGResolution(F,n);
else
RF:=ResolutionGTree(F,n);
fi;
return RF;
fi;
end);
