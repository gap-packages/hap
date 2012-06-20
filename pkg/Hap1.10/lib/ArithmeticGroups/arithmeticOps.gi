
#######################################################
############################
# InstallMethod( \in,
#                "for SL(2,Z)_p",
#               [ IsMatrix,  IsHAPRationalSpecialLinearGroup],
# function ( g, G )
# local P,n,d,facs,m,p,H,K;
# if IsBound(G!.coprimes) then
# m:=G!.coprimes[1];
# p:=G!.coprimes[2];
# P:=[[1,0],[0,p]];
# return P^-1*g*P in SL2Z(1/m);
# fi;
# ##############
# if IsBound(G!.levels) then
# m:=G!.levels[1];
# p:=G!.levels[2];
# H:=SL2Z(1/m);
# K:=ConjugateSL2ZGroup(H,[[1,0],[0,p]]);
# return (g in H and g in K);
# fi;
# end );
# #################################


InstallOtherMethod(RightTransversal,
"Right Transversal for SL2Z(1/m) and its congruence subgroup",
[IsHAPRationalSpecialLinearGroup,IsMatrixGroup],
function(G,H)
local T,trans,NameG,i,j,m,n,p,S,d,enum;
NameG:=Name(G);
i:=Position(NameG,'/');
j:=Position(NameG,']');
m:=Int(NameG{[i+1..j-1]});
n:=H!.levels[1];
p:=H!.levels[2];
if not n=m then return "two groups are not in the same level";fi;
T:=[[1,0],[1,1]];
S:=[[0,-1],[1,0]];

enum:=EnumeratorByFunctions(CollectionsFamily( FamilyObj( G ) ),rec(
      ElementNumber:=function(enum,n)
	  if n<p+1 then return CanonicalRightCountableCosetElement(H,T^n);
	  else return CanonicalRightCountableCosetElement(H,S);
	  fi;
      end,
      NumberElement:=function(enum,elm)
	  local i;
	  if elm*S^-1 in H then return p+1;fi;
	  for i in [1..p] do
	      if elm*T^-i in H then
		   return i;
	      fi;
	  od;
	  return fail;
      end,
      Length:=function(enum)
	  return p+1;
      end,
      PrintObj:=function(enum)
	  Print("RightTransversal of ",Name(H)," in ",Name(G));
      end,
      group:=G,
      subgroup:=H,
));
SetIsSSortedList( enum, true );
return enum;
end);
##########################
InstallOtherMethod(PositionCanonical,
"Position Canonical connects to Right Transversal for SL2Z(1/m) and its congruence subgroup",
[IsFinite,IsDuplicateFree],
function(trans,g)
if IsHAPRationalSpecialLinearGroup(trans!.group) then
return Position(trans,CanonicalRightCountableCosetElement(trans!.subgroup,g));
fi;
end);
##########################
InstallMethod( \in,
               "for CongruenceSubgroupGamma0(p)",
              [ IsMatrix,  IsCongruenceSubgroupGamma0 ],
function ( g, G )
local p;
p:=G!.LevelOfCongruenceSubgroup;
if g in SL(2,Integers) then
if Determinant(g)=1 then
if IsInt(g[2][1]/p) then return true;fi;
fi; 
fi;
############################
return false; 
end );

###################################################################
