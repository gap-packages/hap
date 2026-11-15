
###########################################################
###########################################################
InstallGlobalFunction(BuildLocalHomotopy,
function(K,N)
local Y, homotopy, PseudoBoundary, Boundary, HomotopyRec,
FinalHomotopy, pos,C,D,i,k,n,b,x,s, dim;

Y:=GComplexToFiniteRegularCWRegion(K,N);
Print("Contracting homotopy is valid in characteristic 2 only at present.\n");
#At present this function only works mod 2. We need to adjust
##the orientation on Y so that it agrees with the orientation on K.
if Length(CriticalCells(Y))<>1 then Print("Error\n"); return fail; fi;
C:=ChainComplexWithChainHomotopy(Y);
D:=ChainComplexOfRegularCWComplex(Y);

dim:=-1;;
while K!.dimension(dim+1)>0 do
   dim:=dim+1;
od;

if false then ##############################################
##Each cell in Y is represented as g.e with g a preferred 
##group element and e is the preferred representative cell 
##in the orbit. We modify the boundary of K so that only preferred
##g occur in the expression for the boundary.
PseudoBoundary:=[];
for n in [0..0] do
   PseudoBoundary[n+1]:=
   List([1..K!.dimension(n)],k->[1,0]);
od;
for n in [1..dim] do
   PseudoBoundary[n+1]:=[];
   for k in [1..K!.dimension(n)] do
      pos:=PositionSet(Y!.G_Cells[n+1], Y!.StandardForm(n,[k,One(K!.group)]));
      b:=[];
      for i in [1..-1+Length(Y!.boundaries[n+1][pos])] do
         x:= Y!.boundaries[n+1][pos][i+1];
         x:= Y!.G_Cells[n][x];
         x:=Y!.StandardForm(n-1,x);
         s:= Y!.orientation[n+1][pos][i];
         x:=Y!.StandardForm(n-1,x);
         x[1]:=s*x[1];
         Add(K!.elts,x[2]);
         x[2]:=Length(K!.elts);
         Add(b,x);
      od;
      PseudoBoundary[n+1][k]:=b;
   od;
od;

#############################
Boundary:=function(n,k);
if k>0 then
return PseudoBoundary[n+1][k];
else
return NegateWord(PseudoBoundary[n+1][-k]);
fi;
end;
#############################
K!.boundary:=Boundary;
##
##
fi; ##############################################


############################################
############################################
homotopy:=function(n,x)
local id, sf, pos, v, h, H, z, i, j, g, sgn;
#x=[k,g] is an integer pair representing g.e^n_k
#For the moment will give a very clumsy implementation. 
#The plan is to refine all the "searching" and to use sparse
#matrix notation to speed things up. 

g:=K!.elts[x[2]];
sf:=Y!.StandardForm(n,1*[x[1],g]);
#sgn:=sf[2]^-1*g;
#sgn:=K!.action(n,x[1],Position(K!.elts,sgn)); 
sf:=PositionSet(Y!.G_Cells[n+1],sf); #Y!.G_Cells is a strictly sorted list
v:=[1..Y!.nrCells(n)]*0;
v[sf]:=SignInt(x[1]);
h:=C!.chainHomotopy(v,n);
H:=[];
for i in [1..Length(h)] do
if h[i]<>0 then
  z:=ShallowCopy(Y!.G_Cells[n+2][i]);
  Add(K!.elts,z[2]);
  z[2]:=Length(K!.elts);
    if h[i]<0 then
      z[1]:=-z[1];
      h[i]:=-h[i];
    fi;
      for j in [1..h[i]] do
        Add(H,z);
      od;
fi;
od;

return H;

end;
############################################
############################################

HomotopyRec:=List([0..dim+1],i->[]);
for n in [0..dim] do
for k in [1..K!.dimension(n)] do
HomotopyRec[n+1][k]:=[];
od;
od;

############################################
############################################
FinalHomotopy:=function(n,w)
local absw;

absw:=[AbsInt(w[1]),w[2]];
if not IsBound(HomotopyRec[n+1][absw[1]][w[2]]) then
   HomotopyRec[n+1][absw[1]][w[2]]:=homotopy(n,absw);
fi;

if w[1]>1 then 
   return HomotopyRec[n+1][w[1]][w[2]];
else
   return NegateWord(HomotopyRec[n+1][absw[1]][w[2]]);
fi;

end;
############################################
############################################

K!.homotopy:=FinalHomotopy;
K!.StandardForm:=Y!.StandardForm;
K!.space:=Y;
return true;
end);
###########################################################
###########################################################


