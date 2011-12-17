##########################################################
#DeclareCategory("IsPseudoList",IsObject);

#DeclareRepresentation(  "IsPseudoListRep",
#                        IsComponentObjectRep,
#                        ["elts",
#                         "pos" ]);

PseudoListFamily:=NewFamily( "PseudoListFamily",
                                IsPseudoList,
                                IsPseudoList);

PseudoList:=NewType(PseudoListFamily,IsPseudoListRep);

InstallMethod( ViewObj,
"for PseudoList",
 [IsPseudoList],
 function(R)
Print( R!.lst);
end);

InstallMethod( PrintObj,
"for PseudoList",
 [IsPseudoList],
 function(R)
Print(R!.lst);
end);
#############################################################

#############################################################
InstallGlobalFunction(ListToPseudoList,
function(L)
local Lst,Sslst,Pos;

if IsPseudoList(L) then return L; fi;  #New line ( Nov 2007)

Sslst:=SSortedList(L);

Lst:=L;

Pos:=List(Sslst,x->Position(Lst,x));

return Objectify(PseudoList,
	rec(
	lst:=Lst,
	pos:=Pos,
	sslst:=Sslst));

end);
#############################################################

#############################################################
InstallOtherMethod(ELM_LIST,
"for PseudoLists",
[IsPseudoList,IsInt],
function(PL,i)
return PL!.lst[i];
end);
#############################################################

#############################################################
InstallOtherMethod(Position,
"for PseudoLists",
[IsPseudoList,IsObject,IsInt],
function(PL,g,i)
local p;
p:=Position(PL!.sslst,g,i);
if IsInt(p) then
return PL!.pos[p];
else
return p;
fi;
end);
#############################################################

#############################################################
InstallOtherMethod(Length,
"for PseudoLists",
[IsPseudoList],
function(PL)
return Length(PL!.lst);
end);
#############################################################

#############################################################
InstallOtherMethod(IN,
"for PseudoLists",
[IsObject,IsPseudoList],
function(x,PL)
return x in PL!.sslst;
end);
#############################################################

#############################################################
InstallOtherMethod(Add,
"for PseudoLists",
[IsPseudoList,IsObject],
function(L,x)
local p,shft;
if not x in L then
Add(L!.lst,x);
AddSet(L!.sslst,x);
p:=Position(L!.sslst,x);
shft:=L!.pos{[p..Length(L!.pos)]};
L!.pos:=L!.pos{[1..p-1]};
Add(L!.pos,Length(L!.lst));
Append(L!.pos,shft);
fi;
end);
#############################################################

#############################################################
InstallOtherMethod(Append,
"for PseudoLists",
[IsPseudoList,IsObject],
function(L,K)
local x;
for x in K do
Add(L,x);
od;
end);
#############################################################

#############################################################
InstallOtherMethod(Iterator,
"for PseudoLists",
[IsPseudoList],
function(L)
return Iterator(L!.lst);
end);
#############################################################

#############################################################
InstallOtherMethod(ListOp,
"for PseudoLists",
[IsPseudoList,IsObject],
function(L,f)
return ListOp(L!.lst,f);
end);
#############################################################

