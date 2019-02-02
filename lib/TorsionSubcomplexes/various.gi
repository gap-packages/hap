DeclareGlobalFunction("CellComplexBoundaryCheck");

InstallGlobalFunction(CellComplexBoundaryCheck,
function(C)
local N,i,j,Elts,pos,CLeftCosetElt,Mult,bdr,w,k, Boundary;

i:=0;
while C!.dimension(i)>0 do
    i:=i+1;
od;
N:=i-1;

Elts:=C!.elts;

pos:=function(g)
local posit;

    posit:=Position(Elts,g);
    if posit=fail then
        Add(Elts,g);
        return Length(Elts);
    else
        return posit;
    fi;
end;

CLeftCosetElt:=function(i,j,g)

return pos(CanonicalRightCountableCosetElement
                        (C!.stabilizer(AbsInt(i),j),Elts[g]^-1)^-1);
end;

Mult:=function(L,k,g)
local x,w,t,h,y,vv,LL;
    vv:=[];
    LL:=ShallowCopy(L);
    for x in [1..Length(LL)] do
        w:=Elts[g]*Elts[LL[x][2]];
        t:=CLeftCosetElt(k,AbsInt(LL[x][1]),pos(w));
        Add(vv,[C!.action(k,AbsInt(LL[x][1]),t)*C!.action(k,AbsInt(LL[x][1]),pos(w))*LL[x][1],t]);
    od;
    return vv;
end;

################################
for i in [2..N+1] do
for j in [1..C!.dimension(i)] do
bdr:=StructuralCopy(C!.boundary(i,j));
w:=[];
Boundary:=[];
for k in bdr do
  w:=Mult(C!.boundary(i-1,AbsInt(k[1])),i-2,k[2]);
  if k[1]<0 then w:=NegateWord(w);fi;
  Append(Boundary,w);
od;
Boundary:=AlgebraicReduction(Boundary);
Print([i,j],Boundary);
if not IsEmpty(Boundary) then return false;fi;
od;
od;
return true;
end);


#####################################################################
DeclareProperty("IsGammaSubgroupInSL3Z",IsHAPRationalSpecialLinearGroup);

InstallGlobalFunction( GammaSubgroupInSL3Z,
    function(n)
    local groupname, C, G;
    groupname:=Concatenation("Gamma_",String(n),"inSL3Z");
    C:=ContractibleGcomplex(groupname);
    G:=C!.group;
    G!.number:=n;
    SetName(G,groupname);
    SetIsHAPRationalMatrixGroup(G,true);
    SetIsHAPRationalSpecialLinearGroup(G,true);
    SetIsGammaSubgroupInSL3Z(G,true);
return G;
############################

end);
####################################################################
#####################################################################


######################
InstallMethod( \in,
               "for GammaSubgroupInSL3Z",
              [ IsMatrix,  IsGammaSubgroupInSL3Z ],
function ( g, G )
local groupname, p, n;
    groupname:=Name(G);
    n:=G!.number;
    if n=1 then

        if (g[2][1] mod 2 =0)and(g[3][1] mod 2 =0) then

            return true;
        else return false;
        fi;
    elif n=2 then
      if (g[2][1] mod 2 =0)and(g[3][1] mod 2 =0) and (g[3][2] mod 2 =0)
          then return true;
      else return false;
      fi;
    fi;
end );
######################
