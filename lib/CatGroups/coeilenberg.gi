
###############################################################
###############################################################
InstallGlobalFunction(CohomologySimplicialFreeAbelianGroup,
function(K,n)
local p, q, H, D, ans, mlen, bool;

bool:=false;

if n<0 then return []; fi;
if n=0 then return [0]; fi;

mlen:=EvaluateProperty(K,"moorecomplexlength");
ans:=[];
for  p in [0..n-1] do
   q:=n-p;
   if p<=q*mlen then
if p+1>EvaluateProperty(K,"length") then bool:=true; fi;
      D:=ContractedComplex(E1HomologyPage(K,q,p+1));
      D:=SparseChainComplexToChainComplex(D);
      D:=HomToIntegers(D);
      H:= Cohomology(D,p) ;
      if Length(H)>0 then 
         Append(ans,H);
      fi;
   fi;
od;

#if n>=EvaluateProperty(K,"length") then
if bool then 
Print("WARNING: The computed answer is guaranteed only to be a summand of the correct answer. You'll need to use more terms of the simplicial abelian group to make sure there is no other summand.\n");
fi;

return ans;
end);
###############################################################
###############################################################

#####################################################################
InstallOtherMethod(Cohomology,
"integral homology of simplicial free abelian group",
[IsHapSimplicialFreeAbelianGroup,IsInt],
function(K,n) return CohomologySimplicialFreeAbelianGroup(K,n);
end);
#####################################################################
