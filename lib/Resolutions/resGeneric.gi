#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(ResolutionGenericGroup,
function(arg)

    local G,N, bool, prime, L, M, GG, A, R, fn, epi, F, FhomG, AhomG, Tz,i;

#This function will try to return an efficient resolution, 
#and any resolution returned will have a contracting homotopy.

G:=arg[1];
N:=arg[2];
if Length(arg)>2 then bool:=arg[3]; fi;
if Length(arg)>3 then prime :=arg[4]; fi;

####FINITE################
if Order(G)<infinity then 

#Tz:=TietzeReducedResolution;
Tz:=function(R) return R; end;

   if IsAbelian(G) then
   return ResolutionAbelianGroup_alt(G,N);    ###  CHECK THIS
   if IsPcpGroup(G) then return ResolutionAbelianPcpGroup(G,N);
   else return ResolutionAbelianGroup(G,N); fi;
   fi;

   if Order(G)<64 then
   return Tz(ResolutionFiniteGroup(G,N));
   fi;

   if IsPGroup(G) then
   #return ResolutionNormalSeries(LowerCentralSeries(G),N);
   return ResolutionNormalSeries(BigStepLCS(G,8),N);  #OPTIMIZE!!!
   fi;

L:=LowerCentralSeries(G);
if Maximum(List(L,Order)) <1000 then
return  ResolutionNormalSeries(LowerCentralSeries(G),N);
fi;

   if IsNilpotent(G) then
   return ResolutionNilpotentGroup(G,N);
   fi;

   if IsSolvable(G) then
   return ResolutionSubnormalSeries(CompositionSeries(G),N);
   fi;

   #And is all else fails on a finite group:
   return Tz(ResolutionFiniteGroup(G,N));

fi;
####FINITE DONE###########

####PCP GROUPS############
if IsPcpGroup(G) then

   if IsAbelian(G) then
   return ResolutionAbelianPcpGroup(G,N);
   fi;

   if IsNilpotent(G) then 
   return ResolutionNilpotentGroup(G,N);
   fi;

   if IsAlmostCrystallographic(G) then
   return ResolutionAlmostCrystalGroup(G,N);
   fi;
   
fi;
####PCP DONE##############

####MATRIX GROUP##########

if IsMatrixGroup(G) then

   ###ABELIAN#####
   if IsAbelian(G) then
     ######
     fn:=function(g);
     if Order(g)=infinity then return 0;
     else return Order(g);
     fi;
     end;
     ######
   L:=List(GeneratorsOfGroup(G),x->x);
   A:=AbelianPcpGroup(Length(L),List(L,fn));
   R:=ResolutionAbelianPcpGroup(A,N);
   epi:=EpimorphismFromFreeGroup(A);
   F:=Source(epi);
   FhomG:=GroupHomomorphismByImagesNC(F,G,
                         GeneratorsOfGroup(F),
                         GeneratorsOfGroup(G));
   AhomG:=GroupHomomorphismByFunction(A,G,x->
   Image(FhomG,PreImagesRepresentative(epi,x)));

   R!.group:=G;
   R!.elts:=List(R!.elts,x->Image(AhomG,x));
#if IsBound(G!.gname) then
#if G!.gname="sl2inf" then
#for i in [-100..100] do
#Add(R!.elts,R!.elts[3]^i);
#Add(R!.elts,-R!.elts[3]^i);
#od;
#fi;fi;
   return R;
   fi;
   ####ABELIAN DONE#########

fi;

####MATRIX DONE###########

return fail;
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(ResolutionArithmeticGroup,
function(string, N)
local C, R;

C:=ContractibleGcomplex(string);
R:=FreeGResolution(C,N);
return R;

end);
#####################################################################
#####################################################################

