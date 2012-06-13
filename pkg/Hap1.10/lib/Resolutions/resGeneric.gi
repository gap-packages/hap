#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(ResolutionGenericGroup,
function(G,N)

    local L, M, GG, A, R, fn, epi, F, FhomG, AhomG, Tz;

#This function will try to return an efficient resolution, 
#and any resolution returned will have a contracting homotopy.

####FINITE################
if Order(G)<infinity then 

Tz:=TietzeReducedResolution;
   if IsAbelian(G) then
   return ResolutionAbelianGroup(G,N);
   fi;

   if Order(G)<32 then
   return Tz(ResolutionFiniteGroup(G,N));
   fi;

   if IsPGroup(G) and Order(G)<256 then
   return ResolutionNormalSeries(LowerCentralSeries(G),N);
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

