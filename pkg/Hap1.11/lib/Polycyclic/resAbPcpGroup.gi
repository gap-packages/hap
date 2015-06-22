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
	gens, C, head, tail, 
	R, hom,x,tmp,FreeElts,
	PcpG,
	gens1,hom1,OriginalAppend;

#if IsFinite(G) then return ResolutionAbelianGroup(G,n); fi;
###The last command is a real cheat!!

PcpG:=Pcp(G,"snf");
gens:=List([1..Length(PcpG)],i->PcpG[i]);


if Length(gens)=1 then

if IsFinite(Group(gens)) then return ResolutionFiniteGroup(gens,n);
#if Order(Group(gens))<1000 then return ResolutionFiniteGroup(gens,n);
else



tmp:=ResolutionAbelianGroup([0],n);
FreeElts:=tmp!.elts;
	tmp!.appendToElts:=function(x)
	local a,i,j;
	a:=gens[1];		######################HERE
	for i in [0..10000] do
	j:=false;
	if a^i=x then j:=i; break; fi;
	if a^-i=x then j:=-i; break; fi;
	od;
	Append(FreeElts,[FreeElts[2]^j]);
	Append(tmp!.elts,
	[MappedWord(x,GeneratorsOfGroup(tmp!.group),
	                                        gens)]);
	end;

tmp!.elts:=List(tmp!.elts,x->MappedWord(x,GeneratorsOfGroup(tmp!.group),
                                        gens));


tmp!.group:=G;
return tmp;
fi;
fi;

if Length(gens)=0 then return
ResolutionFiniteGroup([Identity(G)],n); fi;

head:=Subgroup(G,[gens[1]]);
tail:=Subgroup(G,
List([2..Length(gens)],i->gens[i]));



R:=ResolutionDirectProduct(ResolutionAbGroup(head,n),
		        ResolutionAbGroup(tail,n),"internal");

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
