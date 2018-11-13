
##########################################################
##########################################################
InstallGlobalFunction(UniversalCover,
function(X)
local G, dim, Elts, FreeElts, Boundary,PseudoBoundary, gamma, epi;

if not IsHapRegularCWComplex(X) then
Print("The function applies only to regular CW complexes.\n");
return fail;
fi;
if Length(PiZero(X)[1])>1 then
Print("The function applies only to path connected regular CW complexes.\n");
return fail;
fi;

dim:=EvaluateProperty(X,"dimension");
G:=FundamentalGroupOfRegularCWComplex(X,"nosimplify");
epi:=EpimorphismFromFreeGroup(G);
gamma:=G!.edgeToWord;
Elts:=[One(G)];
FreeElts:=[PreImagesRepresentative(epi,One(G))];
PseudoBoundary:=List([1..dim],i->[]);


##################################
Boundary:=function(n,k)
local vts,g, fg, ng, B, bb, bbb, xx, ii, j,kk, pos, bool, indx, bnd;

if n=0 then return []; fi;
if IsBound(PseudoBoundary[n][k]) then
return PseudoBoundary[n][k]; fi;

if n=1 then
#########################
vts:=X!.boundaries[2][k]; 
g:=gamma(k); 
fg:=PreImagesRepresentative(epi,g);
ng:=Position(FreeElts,fg);
if ng=fail then
   Add(Elts,g);
   Add(FreeElts,fg);
   ng:=Length(Elts);
fi;
PseudoBoundary[n][k]:= [[vts[2],1], [vts[3],ng]];
return PseudoBoundary[n][k];
#########################
fi;

#########################
B:=1*X!.boundaries[n+1][k];
B:=1*B{[2..Length(B)]};
bnd:=List(B,b->[b,1]);
indx:=1*B{[2..Length(B)]};

while Length(indx)>0 do

bool:=false;
for ii in indx do
  bb:=1*Boundary(n-1,ii);
  bbb:=1*List(bb,x->x[1]);
  for j in Difference(B,indx) do
  xx:=Boundary(n-1,j);
     for kk in [1..Length(xx)] do
     pos:= Position(bbb, xx[kk][1]);
     if not pos=fail then
       #g:=(Elts[xx[kk][2]]*Elts[bb[pos][2]]^-1); I THINK THIS WAS THE SLIP!
g:=(Elts[bnd[Position(B,j)][2]]*Elts[xx[kk][2]]*Elts[bb[pos][2]]^-1);
       fg:=PreImagesRepresentative(epi,g);
       ng:=Position(FreeElts,fg);
          if ng=fail then
              Add(Elts,g);
              Add(FreeElts,fg);
              ng:=Length(Elts);
          fi;
       bnd[Position(B,ii)][2]:=ng; 
       indx:=Filtered(indx,a->not a=ii); bool:=true;
       break;
     fi;
     if bool then break; fi;
     od;
  if bool then break; fi;
  od;
  if bool then break; fi;
od;

od;
PseudoBoundary[n][k]:= bnd;
#########################

return PseudoBoundary[n][k];
end;
##################################


return Objectify(HapEquivariantCWComplex,
            rec(
            dimension:=X!.nrCells,
            boundary:=Boundary,
            elts:=Elts,
            group:=G,
            stabilizer:=Group(One(G)),
            properties:=
            [["dimension",dim], 
            ]  ));


end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(EquivariantCWComplexToRegularCWComplex,
function(U,H)
local boundaries, bnd, gbnd, ind, trans,pair2int, dimU, n, k,g;

if not IsHapEquivariantCWComplex(U) then
    Print("This function applies only to G-CW-complexes.\n");
    return fail;
fi;

if not IsSubgroup(U!.group,H) then
    Print("The provided group is not a subgroup of the fundamental group of the G-CW-complex.\n");
    return fail;
fi;

##Apply(U!.elts,x->x^-1); #So now we have a right action!
trans:=RightCosets(U!.group,H);
ind:=Length(trans);
if not ind < infinity then
    Print("The provided subgroup is not of finite index.\n");
    return fail;
fi;
dimU:=EvaluateProperty(U,"dimension");


###############################
pair2int:=function(e,gH);
return (e-1)*ind + gH;
end;
###############################

boundaries:=[];
boundaries[1]:=List([1..ind*U!.dimension(0)],i->[1,0]);
for n in [1..dimU] do
boundaries[n+1]:=[];
for k in [1..U!.dimension(n)] do
bnd:=U!.boundary(n,k);

for g in trans do
gbnd:=List(bnd,x->[x[1],
Position(trans, g*U!.elts[x[2]])]);

gbnd:=List(gbnd,x->pair2int(AbsInt(x[1]),x[2]));
gbnd:=Concatenation([Length(gbnd)],gbnd);
Add(boundaries[n+1],gbnd);
od;
od;
od;
boundaries[dimU+2]:=[];

return  HAPRegularCWComplex(boundaries);

end);
##########################################################
##########################################################

