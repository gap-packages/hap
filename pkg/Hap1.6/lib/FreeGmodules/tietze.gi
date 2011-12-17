#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(TietzeReduction,
function(S,c)
local
	ElementaryReduction,
	ElementaryReductionPosNeg,
	b;

#####################################################################
ElementaryReduction:=function(b,c)
local d;
d:=AddFreeWords(b,c);
if Length(d) < Length(c) then return ElementaryReduction(b,d);
else return c; fi;
end;
#####################################################################

#####################################################################
ElementaryReductionPosNeg:=function(b,c)
local d;
d:=ElementaryReduction(b,c);
if Length(d)<Length(c) then return d;
else d:=ElementaryReduction(Negate(b),c);
   if Length(d) < Length(c) then return d;
   else return c;
   fi;
fi;
end;
#####################################################################

for b in S do
c:=ElementaryReductionPosNeg(b,c);
od;

return c;
end);
#####################################################################

