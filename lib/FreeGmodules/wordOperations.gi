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
return List(b, x->[-x[1],x[2]]);
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
	k:=Position(v,[-x[1],x[2]]);
	if (k=fail) then Add(v,x); else
	v[k]:=0;
	fi;
	od;

	return Filtered(v,y->(not y=0));
fi;

if p= 2 then #Words mod 2 are automatically simplified when added.
        v:=List(w,x->[AbsInt(x[1]),x[2]]);
	v:=Collected(v);
        Apply(v,x->[x[1],x[2] mod 2]);
        Apply(v, x->MultiplyWord(x[2],[x[1]]));

return Concatenation(v);
fi;

end);
#####################################################################


#####################################################################
InstallGlobalFunction(AddFreeWords,
function(arg)
local x,u,v,w,r, w2;

if Length(arg[2])<Length(arg[1]) then
v:=arg[1];w:=arg[2];
else
w:=arg[1];v:=arg[2];fi;


if Length(w)=0 then return v; fi;
if Length(v)=0 then return w; fi;

w2:=ShallowCopy(w);

u:=ShallowCopy(w);

for x in v do
   r:=Position(w2,[-x[1],x[2]]);
   if r=fail then Add(u,x);
   else u[r]:=0; w2[r]:=0;fi;
od;

u:=Filtered(u,i->(not i=0));

return u;
end);
#####################################################################


#####################################################################
InstallGlobalFunction(AddFreeWordsModP,
function(v,w,p)
local x, SM,vw,i,j,ab;

if Length(v)=0 then return w; fi;
if Length(w)=0 then return v; fi;


	########################### if p=2 #############################
if p = 2 then

SM:=[];

for x in v do
if not IsBound(SM[x[2]]) then SM[x[2]]:=[]; fi;
ab:=AbsInt(x[1]);
if not IsBound(SM[x[2]][ab]) then
	SM[x[2]][ab]:=[ab,x[2]];
else
Unbind(SM[x[2]][ab]);
fi;
od;

for x in w do
if not IsBound(SM[x[2]]) then SM[x[2]]:=[]; fi;
ab:=AbsInt(x[1]);
if not IsBound(SM[x[2]][ab]) then
        SM[x[2]][ab]:=[ab,x[2]];
else
Unbind(SM[x[2]][ab]);
fi;
od;


vw:=[];
for x in SM do
for j in x do
Add(vw,j);
od;
od;


return vw;

fi;
	########################## fi p=2 ###########################

return AddFreeWords(v,w);
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

if n>0 then u:=w; else u:=List(w,x->[-x[1],x[2]]); fi;

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
u:=AddFreeWords(u, MultiplyWord(y[2],[y[1]]));
od;

return u;
end);
#####################################################################
