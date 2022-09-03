InstallGlobalFunction(CuspidalCohomologyHomomorphism,
function(arg) 
local H,N,A,P,R,S,B,RB,SB,T,ShomT,f,eqmap,t,fn,mapping,eqmap2,cmap,cmap1,cmap2, k,c;

H:=arg[1];
N:=arg[2];
N:=Maximum(N,1);
if Length(arg)>2 then k:=arg[3]-2; else k:=0; fi;
if k<0 then
Print("The forms must be of weight k>=0.\n");
return fail;
fi;

P:=ContractibleGcomplex("bsSL2Z");;
R:=FreeGResolution(P,N+1);;
S:=ResolutionFiniteSubgroup(R,H);;

B:=ContractibleGcomplex("bsSL2Ztmp");;B!.group:=P!.group;;
RB:=FreeGResolution(B,N+1);;
SB:=ResolutionFiniteSubgroup(RB,H);;

f:=GroupHomomorphismByFunction(SB!.group,S!.group,x->x);;
eqmap:=EquivariantChainMap(SB,S,f);;

#t:=Length(H!.tree);  #Not true!!!!!
t:=SB!.dimension(0);
   fn:=function(x);
   if AbsInt(x[1])<=t then
      return [SignInt(x[1])*(AbsInt(x[1])+t),x[2]];
   else 
      return [SignInt(x[1])*(3*t+AbsInt(x[1])),x[2]];
   fi;
   end;

##############################
mapping:=function(w,n);
return List(w,x->fn(x));
end;
##############################
eqmap!.mapping:=mapping;

T:=ResolutionSL2Z_alt(N+1);
T:=ResolutionFiniteSubgroup(T,H);
T:=TietzeReducedResolution(T);
f:=GroupHomomorphismByFunction(S!.group,T!.group,x->x);;
eqmap2:=EquivariantChainMap(S,T,f);


A:=HomogeneousPolynomials(R!.group,k);
cmap1:=HomToIntegralModule(eqmap,A);
cmap2:=HomToIntegralModule(eqmap2,A);
cmap:=Compose(cmap1,cmap2);

c:=Cohomology(cmap,N);
return c;
end);
