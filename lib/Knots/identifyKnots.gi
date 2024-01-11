
#################################################
#################################################
InstallGlobalFunction(HAP_KnotGroupInv,
function(arg)
local GG, N,G, L, LL, M, i;

GG:=arg[1];
N:=7;
if Length(arg)>1 then
N:=arg[2];
fi;

G:=SimplifiedFpGroup(GG);
L:=LowIndexSubgroupsFpGroup(G,N);
LL:=List(L,g->Index(G,g));
M:=List([1..N],i->Filtered([1..Length(LL)],a->LL[a]=i)  );
M:=List(M,m->L{m});
M:=List(M,m->SortedList(List(m,g->AbelianInvariants(g))));
return M;
end);
#################################################
#################################################


#################################################
#################################################
InstallGlobalFunction(IdentifyKnot,
function(gc)
local  N, L, x, y, str, W, Inv;

if not IsBoundGlobal("HAP_knot_census") then
ReadPackage("HAP", "lib/Knots/census.txt");
fi;

if IsList(gc) then
if Length(SSortedList( List(Flat(gc[1]),i->AbsInt(i)) ))>11 then
Print("Currently this function is implemented for knots presented with fewer than 12 crossings. \n");
fi;
W:=WirtingerGroup(gc);
fi;
if IsHapPureCubicalComplex(gc) then
if EvaluateProperty(gc,"knot")=true then
W:=WirtingerGroup(GaussCodeOfPureCubicalKnot(gc));;
else
W:=KnotGroup(gc);  
fi;
fi;

for N in [3..7] do
   Inv:=function(W);
   return HAP_KnotGroupInv(W,N);
   end;

   L:=Filtered(ValueGlobal("HAP_knot_census"),x->
      x[2]{[1..Minimum(Length(x[2]),N)]} = Inv(W){[1..Minimum(Length(x[2]),N)]});

   if Length(L)=0 then Print("Identification has failed.\n");
   return fail; fi;

   if Length(L)=1 then break; fi;
od;

L:=List(L,x->x[1]);
str:=0;
for x in L do
for y in x do
if str=0 then
Print("PrimeKnot("); str:=1;
else
Print(" + PrimeKnot(");
fi;
Print(y[1]);
Print(",");
Print(y[2]);
Print(")");
od;
od;

Print("    modulo reflections of components. \n");
return L;
end);
#################################################
#################################################
