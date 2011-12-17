#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionAbelianPcpGroup,
function(arg)
local  
	ResolutionAbGroup;

#####################################################################
#####################################################################
ResolutionAbGroup:=function(G,n)
local 
	gens,gensInf, gensFin, C, head, tail, 
	R, hom,x,tmp,FreeElts,gensG ,
	PcpGT;

gensFin:= GeneratorsOfGroup(TorsionSubgroup(G));
PcpGT:=Pcp(G,TorsionSubgroup(G));
gensInf:=List([1..Length(PcpGT)],i->PcpGT[i]);


if Length(gensFin)>0 then
gensFin:=TorsionGeneratorsAbelianGroup(Group(gensFin));
else
gensFin:=[];
fi;

if Length(gensInf)>0 then
gensInf:=ReduceGenerators(gensInf,Group(gensInf));
#gensInf:=GeneratorsOfGroup(Group(gensInf));
fi;

gens:=Concatenation(gensInf,gensFin);
gensG:=gens;

if Length(gens)=1 then

if IsFinite(Group(gens)) then return ResolutionFiniteGroup(gens,n);
else
tmp:=ResolutionAbelianGroup([0],n);
FreeElts:=tmp.elts;
	tmp.appendToElts:=function(x)
	local a,i,j;
	a:=MappedWord(FreeElts[2],[FreeElts[2]], gensG);
	for i in [0..10000] do
	j:=false;
	if a^i=x then j:=i; break; fi;
	if a^-i=x then j:=-i; break; fi;
	od;
	Append(FreeElts,[FreeElts[2]^j]);
	return(List(FreeElts, b-> 
	MappedWord(FreeElts[2],[FreeElts[2]], gensG)));
	end;

tmp.elts:=List(tmp.elts,x->MappedWord(x,GeneratorsOfGroup(tmp.group),
                                        gensG));


tmp.group:=G;
return tmp;
fi;
fi;

if Length(gens)=0 then return
ResolutionFiniteGroup([Identity(G)],n); fi;

head:=Subgroup(G,[gens[1]]);
tail:=Subgroup(G,List([2..Length(gens)],i->gens[i]));

R:=ResolutionDirectProduct(ResolutionAbGroup(head,n),
		        ResolutionAbGroup(tail,n),"internal");
hom:=GroupHomomorphismByFunction(R.group,G,x->
Image(R.firstProjection,x)*
Image(R.secondProjection,x));
R.elts:=List(R.elts,x->Image(hom,x));
R.group:=G;
return R;

end;
#####################################################################
#####################################################################

if IsPcpGroup(arg[1]) and IsAbelian(arg[1]) and IsInt(arg[2]) then
return ResolutionAbGroup(arg[1],arg[2]);  fi;

if (not IsPcpGroup(arg[1])) and IsAbelian(arg[1]) and IsInt(arg[2]) then
return ResolutionAbelianGroup(arg[1],arg[2]);  fi;


Print("The first argument must be an abelian Pcp group. The second argument must be a positive integer. \n");
return fail;
end);
