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
local w,p,v,vv,vvv,vvvv,x,j,k,pos,pos2;

w:=arg[1];
if Length(arg)=2 then p:=arg[2]; else p:=0; fi;

if p=0  then

vv:=SSortedList(List(w,x->[AbsInt(x[1]),x[2]]));   #ADDED AUGUST 2022
if Length(vv) = Length(w) then return SortedList(w); fi;

#######################################
         v:=Collected(w);
         vv:=List(v,x->x[1]);
         vvv:=List(vv,x->[AbsInt(x[1]),x[2]]);
         vvv:=SSortedList(vvv);
         vvvv:=[];
         for x in vvv do
             pos:=Position(vv,x);
             pos2:=Position(vv,[-x[1],x[2]]);
             k:=0;
             if not pos=fail then k:=k+v[pos][2]; fi;
             if not pos2=fail then k:=k-v[pos2][2]; fi;
             if not k=0 then Append(vvvv,MultiplyWord(k,[x])); fi;
od;
return vvvv;
#######################################
fi;

        v:=Collected(w);
        Apply(v,x->[x[1],x[2] mod p]);
        Apply(v, x->MultiplyWord(x[2],[x[1]]));
        v:=Collected(Concatenation(v));
        Apply(v,x->[x[1],x[2] mod p]);
        Apply(v, x->MultiplyWord(x[2],[x[1]]));
        return Concatenation(v);

end);
#####################################################################


#####################################################################
InstallGlobalFunction(AlgebraicReduction_alt,
function(arg)      
local ww,y,s,w,p,v,x,j,k,pos;

w:=arg[1];
if Length(arg)=2 then p:=arg[2]; else p:=0; fi;

if p= 0  then

        s:=SSortedList(w);
        Apply(s,x->[AbsInt(x[1]),x[2]]);
        s:=SSortedList(s);
        v:=[1..Length(s)]*0; 
        for x in w do
        y:=[AbsInt(x[1]),x[2]];
        pos:=PositionSorted(s,y);
        v[pos]:=v[pos]+SignInt(x[1]); 
        od;
        ww:=[];
        for j in [1..Length(v)] do
        if v[j]>0 then 
        for k in [1..v[j]] do
        Add(ww,s[j]);
        od;
        else
        for k in [1..-v[j]] do
        Add(ww,[-s[j][1],s[j][2]]);
        od;
        fi;

        od;

        return ww;
fi;


	v:=Collected(w);
        Apply(v,x->[x[1],x[2] mod p]);
        Apply(v, x->MultiplyWord(x[2],[x[1]]));
        v:=Collected(Concatenation(v));
        Apply(v,x->[x[1],x[2] mod p]);
        Apply(v, x->MultiplyWord(x[2],[x[1]]));

	return Concatenation(v);


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


if Length(w)=0 then  return v; fi;
if Length(v)=0 then  return w; fi;

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
InstallGlobalFunction(AppendFreeWord,
function(v,w)
local x,u,r, w2;

if Length(w)=0 then return false; fi;

for x in w do
   r:=Position(v,[-x[1],x[2]]);
   if r=fail then Add(v,x);
   else v[r]:=0;fi;
od;

v:=Filtered(v,x->not x=0);

return true;
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
local v, u, i,x;
v:=[];

u:=1*w; 
if n<=0 then for x in u do x[1]:=-x[1]; od; fi;

for i in [1..AbsoluteValue(n)] do
Append(v,u);
od;

return 1*v;
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

#####################################################################
InstallGlobalFunction(OppositeGroup,
function(w)
local u, v, y;

return Objectify(HapOppositeElement,
                    rec(
                    element:=w) );
end);
#####################################################################

#####################################################################
InstallGlobalFunction(QuotientGroup,
function(w,D)
local u, v, y;

return Objectify(HapQuotientElement,
                    rec(
                    element:=w*D) );
end);
#####################################################################


#####################################################################
InstallOtherMethod( \*,
    "composition in opposite group",
    [ IsHapOppositeElement, IsHapOppositeElement],

function(x,y) local w;
w:=y!.element*x!.element;
return Objectify(HapOppositeElement,
                    rec(
                    element:=w) );
end);
#####################################################################

#####################################################################
InstallOtherMethod( \*,
    "composition in quotient group",
    [ IsHapQuotientElement, IsHapQuotientElement],

function(x,y) local w;
w:=x!.element[1]*y!.element;
return Objectify(HapQuotientElement,
                    rec(
                    element:=w) );
end);
#####################################################################


#####################################################################
InstallOtherMethod( \^,
    "inverse in opposite group",
    [ IsHapOppositeElement, IsInt],

function(x,n) local w;
w:=(x!.element)^n;
return Objectify(HapOppositeElement,
                    rec(
                    element:=w) );
end);
#####################################################################

#####################################################################
InstallOtherMethod( \^,
    "inverse in quotient group",
    [ IsHapQuotientElement, IsInt],

function(x,n) local w;
w:=(x!.element[1])^n;
w:=(w*x!.element[1]^-1)*x!.element;
return Objectify(HapQuotientElement,
                    rec(
                    element:=w) );
end);
#####################################################################


#####################################################################
InstallOtherMethod( \^,
    "inverse in opposite group",
    [ IsHapOppositeElement, IsHapOppositeElement],

function(x,y) local w;
w:=(x!.element)^y!.element;
return Objectify(HapOppositeElement,
                    rec(
                    element:=w) );
end);
#####################################################################

#####################################################################
InstallOtherMethod( \^,
    "inverse in quotient group",
    [ IsHapQuotientElement, IsHapQuotientElement],

function(x,y) local w,D;
w:=x!.element[1]^y!.element[1];
D:=x!.element[1]^-1*x!.element;
w:=w*D;
return Objectify(HapQuotientElement,
                    rec(
                    element:=w) );
end);
#####################################################################


#####################################################################
InstallOtherMethod( One,
    "identity in opposite group",
    [ IsHapOppositeElement],

function(x) local w;
w:=One(x!.element);
return Objectify(HapOppositeElement,
                    rec(
                    element:=w) );
end);
#####################################################################

#####################################################################
InstallOtherMethod( One,
    "identity in quotient group",
    [ IsHapQuotientElement],

function(x) local w;
w:=x!.element[1]^-1;
w:=w*x!.element;
return Objectify(HapQuotientElement,
                    rec(
                    element:=w) );
end);
#####################################################################


#####################################################################
InstallOtherMethod( InverseMutable,
    "inverse in opposite group",
    [ IsHapOppositeElement],

function(x) local w;
w:=InverseMutable(x!.element);
return Objectify(HapOppositeElement,
                    rec(
                    element:=w) );
end);
#####################################################################

#####################################################################
InstallOtherMethod( InverseMutable,
    "inverse in quotient group",
    [ IsHapQuotientElement],

function(x) local w;
w:=List(x!.element,InverseMutable);
return Objectify(HapOppositeElement,
                    rec(
                    element:=w) );
end);
#####################################################################


#####################################################################
InstallOtherMethod( \=,
    "equality in opposite group",
    [ IsHapOppositeElement, IsHapOppositeElement],

function(x,y) ;
return  x!.element = y!.element;
end);
#####################################################################

#####################################################################
InstallOtherMethod( \=,
    "equality in quotient group",
    [ IsHapQuotientElement, IsHapQuotientElement],

function(x,y) ;
return  SSortedList(x!.element) = SSortedList(y!.element);
end);
#####################################################################


#####################################################################
InstallOtherMethod( \<,
    "equality in opposite group",
    [ IsHapOppositeElement, IsHapOppositeElement],

function(x,y) ;
return  x!.element < y!.element;
end);
#####################################################################

#####################################################################
InstallOtherMethod( \<,
    "equality in opposite group",
    [ IsHapQuotientElement, IsHapQuotientElement],

function(x,y) ;
return  SSortedList(x!.element) < SSortedList(y!.element);
end);
#####################################################################



###############################################################
###############################################################
InstallGlobalFunction(ResolutionBoundaryOfWord,
function(R,n,W)
local x, DW, Boundary, Dimension,Elts,pos, ans;

Dimension:=R!.dimension;
Boundary:=R!.boundary;
Elts:=R!.elts;
DW:=[];

for x in W do
ans:=Boundary(n,x[1]);
ans:=List(ans, a->[a[1],Elts[a[2]]]);
ans:=List(ans, a->[a[1],Elts[x[2]]*a[2]]);
Append(DW,ans);
od;

DW:= AlgebraicReduction(DW);
for x in DW do
if not x[2] in Elts then Add(Elts,x[2]);fi;
od;
DW:=List(DW,x->[x[1],Position(Elts,x[2])]);
return DW;

end);
###############################################################
###############################################################



