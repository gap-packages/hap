#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(Negate,
function(p);
return [-p[1],p[2]];
end);
#####################################################################

#####################################################################
InstallGlobalFunction(NegateWord,
function(b);
return List(b, x->Negate(x));
end);
#####################################################################

#####################################################################
InstallGlobalFunction(AlgebraicReduction,
function(arg)      
local w,p,v,x,j,k,pos;

w:=arg[1];
if Length(arg)=2 then p:=arg[2]; else p:=0; fi;

if p= 0 or p>2 then
	v:=[];
	for x in w do
	k:=Position(v,Negate(x));
	if (k=fail) then Append(v,[x]); else
	v[k]:=0;
	fi;
	od;

	return Filtered(v,y->(not y=0));
fi;

if p= 2 then
	v:=[];
	for x in w do
	if x in v then RemoveSet(v,x); else
	if Negate(x) in v then RemoveSet(v,Negate(x));
	else AddSet(v,x); fi; fi;
	od;

	return v;
fi;

end);
#####################################################################

#####################################################################
InstallGlobalFunction(AddWords,
function(v,w)
local x,u,AddLetter, w2;

if Length(w)=0 then return v; fi;
if Length(v)=0 then return w; fi;

w2:=ShallowCopy(w);

   ##################################################################
   AddLetter:=function(x,u)
   local r;

   r:=Position(w2,Negate(x));
   if r=fail then Append(u,[x]);
   else u[r]:=0; w2[r]:=0;fi; 
   end;
   ##################################################################

u:=ShallowCopy(w);
for x in v do
AddLetter(x,u);
od;

u:=Filtered(u,i->(not i=0));
return u;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(AddWordsModP,
function(v,w,p)
local  x, w2,  AddLetter;

if Length(v)=0 then return w; fi;
if Length(w)=0 then return v; fi;

	########################### if p=2 #############################
if p = 2 then
w2:=SSortedList(w);

        #############################################################
        AddLetter:=function(x,u)
        local nx,r ;

        nx:=Negate(x);
        if nx in w2 then RemoveSet(w2,nx);
        else if x in w2 then RemoveSet(w2,x);
        else AddSet(w2,nx);fi;fi;

        end;
        #############################################################

for x in v do
AddLetter(x,w2);
od;

return w2;
fi;
	########################## fi p=2 ###########################

return AddWords(v,w);
end);
#####################################################################

#####################################################################
InstallGlobalFunction(PrintZGword,
function(w,Elts)
local word, basis, coeffs, i, x, PrintGroupElt;

word:=AlgebraicReduction(w);
basis:=SSortedList(List(word, x->AbsoluteValue(x[1])));
coeffs:=[];

for i in basis do
coeffs[i]:=[];
od;

for x in word do
Append(coeffs[AbsoluteValue(x[1])],[SignInt(x[1])*x[2]]);
od;

	#############################################################
	PrintGroupElt:=function(i);
	if i>0 then Print("+ ", Elts[i], " ");
	else
	Print("- ", Elts[-i], " ");
	fi;
	end;
	#############################################################

for i in basis do
Print("( ");
	for x in coeffs[i] do
	PrintGroupElt(x);
	od;

if not i=Maximum(basis) then 
Print(")E",i," +  "); 
else
Print(")E",i,"\n");
fi;
od;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(MultiplyWord,
function(n,w)
local v, u, i;
v:=[];

if n>0 then u:=w; else u:=NegateWord(w); fi;

for i in [1..AbsoluteValue(n)] do
Append(v,u);
od;

return v;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(WordModP,
function(w,p)
local u, v, y;

v:=Collected(w);
v:=List(v,x->[x[1],x[2] mod p]);
u:=[];
for y in v do
u:=AddWords(u, MultiplyWord(y[2],[y[1]]));
od;

return u;
end);
#####################################################################
