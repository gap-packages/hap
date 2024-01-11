
#######################################################################
#0
#F  SL2ZResolution
##  Input: A pair of positive integers (m,n) 
##         
##  Output: The first n+1 terms of a free ZG-resolution  
##          where G is SL2Z(1/m)
##

InstallGlobalFunction(SL2ZResolution,
function(m,n)
local l,p,k,
    C,R,T,RH,RK,RGamma,H,K,Gamma,D,G,F,RF;
    
	l:=Factors(m);
    p:=l[Length(l)];
    k:=m/p;
    R:=ResolutionSL2Z(1,n);
    if m=1 then
        return R;
    else
    RH:=SL2ZResolution(k,n);
    H:=RH!.group;
	
## Create resolution for K
    RK:=ConjugatedResolution(RH,[[1,0],[0,p]]);
    RK!.group:=ConjugateSL2ZGroup(H,[[1,0],[0,p]]);

## Create resolution for Gamma
    Gamma:=CongruenceSubgroup(k,p);
    RGamma:=ResolutionFiniteSubgroup(RH,Gamma);
    SetName(Gamma,"Gamma");
	
## Create tree of groups
    D:=[RH,RK,RGamma];
    G:=SL2Z(1/m); 

## Compute a non-free complex for SL(2,Z[1/p])
    F:=TreeOfResolutionsToSL2Zcomplex(D,G); 
    RF:=ResolutionGTree(F,n);

    return RF;
    fi;
end);

################### end of SL2ZResolution ############################
