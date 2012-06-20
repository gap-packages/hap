####################################################
####################################################
InstallGlobalFunction(QuillenComplex,
function(G,p)
local
	S,Subs,ElAb,cl,MaxSimps,i,dim,arrows,simp,IsComp,s,t,x,P,P1;


S:=LatticeSubgroups(G);
S:=ConjugacyClassesSubgroups(S);
Subs:=[];

#####Create list of all elementary abelian subgroups#####
for cl in S do
ElAb:=[];
for i in [1..Size(cl)] do
Add(ElAb,ClassElementLattice(cl,i));
od;
ElAb:=Filtered(ElAb,x->IsElementaryAbelian(x) and p in PrimePowersInt(Order(x)) );
Append(Subs,ElAb);
od;
#####Created#############################################

MaxSimps:=[];
dim:=Prank(SylowSubgroup(G,p));
arrows:=[];

######################Loops##############################
for i in [1..dim-1] do
arrows[i]:=[]; 
P:=p^i; P1:=p*P;
for s in Filtered(Subs,x->Size(x)=P) do
for t in Filtered(Subs,x->Size(x)=P1) do
if IsSubgroup(t,s) then Add(arrows[i],[s,t]); fi;
od;
od;
od;
#####################Loops done###########################

#########################
IsComp:=function(x)
local i,bool;
bool:=true;

for i in [1..Length(x)-1] do
if not IsSubgroup(x[i+1][1],x[i][2]) then bool:=false; break; fi;
od;

return bool;
end;
#########################

#########Final loop#########################
for x in Cartesian(List([1..Length(arrows)],i->arrows[i])) do
if IsComp(x) then 
simp:=List(x,a->a[1]);
Add(simp,x[Length(x)][2]);
Add(MaxSimps,simp); fi;
od;
############################################

return MaximalSimplicesToSimplicialComplex(MaxSimps);

end);
####################################################
####################################################
