
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
###################################################################
