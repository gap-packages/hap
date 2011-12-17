#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionAbelianPcpGroup,
function(arg)
local  
	ResolutionAbGroup;

#####################################################################
#####################################################################
ResolutionAbGroup:=function(G,n)
local gens,gensInf, gensFin, C, head, tail, R, hom,x,tmp ;

gensInf:=[];
gensFin:=[];

for x in GeneratorsOfGroup(G) do
if Order(x)=infinity then Append(gensInf,[x]); fi;
if Order(x)<infinity  then Append(gensFin,[x]); fi;
od;

if Length(gensFin)>0 then
gensFin:=TorsionGeneratorsAbelianGroup(Group(gensFin));
else
gensFin:=[];
fi;

if Length(gensInf)>0 then
gensInf:=ReduceGenerators(gensInf,Group(gensInf));
fi;

gens:=Concatenation(gensInf,gensFin);


if Length(gens)=1 then

if IsFinite(G) then return ResolutionFiniteGroup(gens,n);
else
tmp:=ResolutionAbelianGroup([0],n);

tmp.elts:=List(tmp.elts,x->MappedWord(x,GeneratorsOfGroup(tmp.group),
					GeneratorsOfGroup(G)));
tmp.group:=G;
return tmp;
fi;
fi;

if Length(gens)=0 then return
ResolutionFiniteGroup([Identity(G)],n); fi;

head:=Group([gens[1]]);
tail:=Group(List([2..Length(gens)],i->gens[i]));

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
