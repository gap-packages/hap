##########################################################
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
InstallGlobalFunction(LazyList,
function(arg)
local fn, name, length, opts, L, len, pos, x;

fn:=arg[1];
if Length(arg)=2 then opts:=arg[2];
else opts:=[];
fi;

name:="lazy list";
length:=infinity;
for x in opts do
if x[1]="name" then name:=x[2]; fi;
if x[1]="length" then length:=x[2]; fi;
od;

##################
len:=function()
return length;
end;
##################

##################
pos:=function(t)
local i, bool;
i:=0;
bool:=true;
while bool and i< length do
i:=i+1;
bool:= not (fn(i)=t) ;
od;
if not bool then return i;
else return fail; fi;
end;
##################


L:=     Objectify(PseudoList,
        rec(
        lst:=[name],
        pos:=fail,
        sslst:=fail,
        eltfun:=fn,
        lnthfun:=len,
        posfun:=pos));
SetIsPseudoListWithFunction(L,true);

return L;
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
InstallOtherMethod(ELM_LIST,
"for PseudoLists with functions",
[IsPseudoListWithFunction,IsInt],
function(PL,i)
return PL!.eltfun(i);
end);
#############################################################


#############################################################
InstallOtherMethod(Position,
"for PseudoLists",
[IsPseudoList,IsObject,IsInt],
function(PL,g,i)
local p;
#if i=0 and IsBound(PL!.posfun) then return PL!.posfun(g); fi;
p:=Position(PL!.sslst,g,i);
if IsInt(p) then
return PL!.pos[p];
else
return p;
fi;
end);
#############################################################
 


#############################################################
InstallOtherMethod(Position,
"for PseudoLists",
[IsPseudoListWithFunction and IsPseudoList,IsObject,IsInt],
function(PL,g,i)
local j;
if i=0 and IsBound(PL!.posfun) then return PL!.posfun(g); fi;
j:=i;   #CHANGED 30July 2019
while j<PL!.lnthfun() do
j:=j+1;
if PL!.eltfun(j)=g then return j; fi;
od;
return fail;
end);
#############################################################

#############################################################
InstallOtherMethod(Position,
"for PseudoLists",
[IsPseudoListWithFunction and IsPseudoList,IsObject], 1000,
function(PL,g)
local p;
return PL!.posfun(g);
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
InstallOtherMethod(Length,
"for PseudoLists",
[IsPseudoListWithFunction],
function(PL) 
return  PL!.lnthfun();
end);
#############################################################

#############################################################
InstallOtherMethod(IsEmpty,
"for PseudoLists",
[IsPseudoList],
function(PL)
return 0=Length(PL!.lst);
end);
#############################################################

#############################################################
InstallOtherMethod(IsEmpty,
"for PseudoLists",
[IsPseudoListWithFunction],
function(PL)
return  0=PL!.lnthfun();
end);
#############################################################

#############################################################
InstallOtherMethod(Size,
"for PseudoLists",
[IsPseudoList],
function(PL)
return Length(PL!.lst);
end);
#############################################################

#############################################################
InstallOtherMethod(Size,
"for PseudoLists",
[IsPseudoListWithFunction],
function(PL)
return  PL!.lnthfun();
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
InstallOtherMethod(IN,
"for PseudoLists",
[IsObject,IsPseudoListWithFunction],
function(x,PL)
return not fail = PL!.posfun(x);
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
InstallOtherMethod(Add,
"for PseudoLists",
[IsPseudoListWithFunction,IsObject],
function(L,x)
L!.addfun(x);
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

#############################################################
InstallOtherMethod(ListOp,
"for PseudoLists",
[IsPseudoListWithFunction,IsObject],
function(L,f)
return List(List([1..Length(L)],i->L!.eltfun(i)) , f);
end);
#############################################################


