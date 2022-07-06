
######################################################
######################################################
InstallGlobalFunction(PolymakeFaceLattice,
function(arg)
local polygon, toggle, F, FF, FA, D, FL, i, j, cnt, lst,bool;

#WARNING: This fixed is used in HAPcryst. Don't fiddle with toggle
polygon:=arg[1];
if Length(arg)=2 then toggle:=false;
else toggle:=true;
fi;
bool:=false;

FF:=Polymake(polygon,"FACES");
FA:=Polymake(polygon,"ADJACENCY");
D:=Polymake(polygon,"F_VECTOR");
if  not Size(Filtered(FF,x->Size(x)=1))=D[1] then D:=Reversed(D); fi;
if not FF[1]=[] then FF:=Reversed(FF); FA:=Reversed(FA); bool:=true; fi;
if toggle then F:=FF; else F:=FA;  fi;

FL:=[];;
FL[1]:=[F[1]];
cnt:=1;
for i in [1..Length(D)] do
FL[i+1]:=[];
lst:=[cnt+1..cnt+D[i]]; 
if bool then lst:=Reversed(lst); fi;
for j in lst do
Add(FL[i+1],F[j]);
od;

cnt:=cnt+D[i];
od;

if FL[1]=[[]] then FL:=Reversed(FL); else Add(FL,[[]]); fi;
if not toggle then Remove(FL,1); fi;
return FL;
end);
############################################################
############################################################
