#(C) Graham Ellis 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(Syzygy,
function(R,g)
local
	Homotopy, Syz,  i, Elts, Mult, n,
	R203, R1023, R1203, R102;

n:=Length(g);

#####################################################################
Elts:=function(x);
return Position(R!.elts,x);
end;
#####################################################################

#####################################################################
Mult:=function(x,w);	#x is in the group G.
return List(w,y->[y[1],Position(R!.elts,x*R!.elts[y[2]])]);
end;
#####################################################################

#####################################################################
Homotopy:=function(i,w)
local ans, x;
ans:=[];

for x in w do
ans:=AddFreeWords(ans,R!.homotopy(i,x));
od;
return ans;
end;
#####################################################################



if n<1 or n>3 then 
Print("ERROR: Syzygy() is so far only implemented for 1-, 2- and 3-syzygies. \n");
return fail;
fi;

##################### IF N=1 ########################################
if n=1 then 
Syz:=[[1,Elts(g[1])]];
Syz:=AddFreeWords(Syz,[[-1,Elts(g[1]*g[1]^-1)]]);
return Homotopy(0,Syz);
fi;
#####################################################################

##################### IF N=2 ########################################
if n=2 then
Syz:= Syzygy(R,[g[1]]);
Syz:=AddFreeWords(Syz,Mult(g[1],Syzygy(R,[g[2]])));
Syz:=AddFreeWords(Syz,NegateWord(Syzygy(R,[g[1]*g[2]]))) ;
return Homotopy(1,Syz);
fi;
##################### FI N=2 ########################################

##################### IF N=3 ########################################
Syz:=Mult(g[1], Syzygy(R,[g[2],g[3]]));
Syz:=AddFreeWords(Syz,Syzygy(R,[g[1],g[2]*g[3]]));
Syz:=AddFreeWords(Syz,NegateWord(Syzygy(R,[g[1]*g[2],g[3]])));
Syz:=AddFreeWords(Syz,NegateWord(Syzygy(R,[g[1],g[2]])));

return Homotopy(2,Syz);
#################### FI N=3 #########################################


end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(StandardCocycle,
function(arg)
local 
	 R,f,n,q,Standard;
R:=arg[1];
f:=arg[2];
n:=arg[3];
if Length(arg)>3 then q:=arg[4]; else q:=0; fi;


#####################################################################
Standard:=function(arg)
local S,v,x,g,h,k,lst;

g:=arg[1]; lst:=[g];
if Length(arg)>1 then h:=arg[2]; lst:=[g,h]; fi;
if Length(arg)>2 then k:=arg[3]; lst:=[g,h,k]; fi;

S:=Syzygy(R,lst);
Apply(S,x->x[1]);

v:=List([1..R!.dimension(n)],x->0);

for x in S do
v[AbsoluteValue(x)]:=v[AbsoluteValue(x)] + SignInt(x);
od;

if q=0 then
return v*f;
else
return v*f mod q; fi;
end;
#####################################################################

return Standard;
end);
#####################################################################
#####################################################################
