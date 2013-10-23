####################################################
####################################################
InstallGlobalFunction(QuillenComplex,
function(G,p)
local
	SubsCl, Subs, fn, cl, MaxSubs, bool, MaxSimps, s,t,x,m,mm,tmp;


SubsCl:=ConjugacyClassesSubgroups(LatticeByCyclicExtension(G, 
      IsElementaryAbelian, true));;

SubsCl:=Filtered(SubsCl,cl->
        IsElementaryAbelianPGroup(ClassElementLattice(cl,1),p));
SubsCl:=Filtered(SubsCl,cl->Order(ClassElementLattice(cl,1))>1);
Subs:=[];
for cl in SubsCl do
for x in [1..Size(cl)] do
Add(Subs,ClassElementLattice(cl,x));
od;
od;

Unbind(SubsCl);

###################
fn:=function(A,B);
return Order(A)>=Order(B);
end;
###################

Sort(Subs,fn);

MaxSimps:=[];
for s in Subs do
bool:=true;
for t in MaxSimps do
if IsSubgroup(t,s) then bool:=false; break; fi;
od;
if bool then Add(MaxSimps,s); fi;
od;

Unbind(Subs);
MaxSimps:=List(MaxSimps,s->[s]);

bool:=true;

while bool do
bool:=false;

for x in [1..Length(MaxSimps)] do
if IsBound(MaxSimps[x]) then

m:=MaxSimps[x];
if Order(m[Length(m)])>p then 
for t in MaximalSubgroups(m[Length(m)]) do
mm:=Concatenation(m,[t]);
Add(MaxSimps,mm); if Order(t)>p then bool:=true; fi;
Unbind(MaxSimps[x]);
od;
fi;

fi;
od;
od;

tmp:=MaxSimps;
MaxSimps:=[];
for m in tmp do
Add(MaxSimps,m);
od;
Unbind(tmp);
return MaximalSimplicesToSimplicialComplex(MaxSimps);

end);
####################################################
####################################################
