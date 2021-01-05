#The method HAP_EquivalenceClasses is taken from the ResClasses package by Stefann Kohl

###########################################################
###########################################################
InstallOtherMethod( HAP_EquivalenceClasses,
                    "for a list and a relation or a class invariant (RCWA)",
                    ReturnTrue, [ IsList, IsFunction ], 0,

  function ( list, relation )

    local  classes, invs, longestfirst, byinvs, elm, pos, inserted, count;

    if IsEmpty(list) then return []; fi;

    longestfirst := function(c1,c2) return Length(c1) > Length(c2); end;
    byinvs := function(c1,c2) return relation(c1[1]) < relation(c2[1]); end;

      classes := [[list[1]]]; count := 0;
      for elm in list{[2..Length(list)]} do
        inserted := false; count := count + 1;
        for pos in [1..Length(classes)] do
          if relation(elm,classes[pos][1]) then
            Add(classes[pos],elm);
            inserted := true;
            break;
          fi;
        od;
        if   not inserted
        then classes := Concatenation(classes,[[elm]]); fi;
        if   count mod 100 = 0 # rough performance heuristics ...
        then Sort(classes,longestfirst); fi;
      od;
      Sort(classes,longestfirst);
 return classes;
  end );
###########################################################
###########################################################

#############################################
#############################################
GroupIsomorphismRepresentatives:=function(arg)
local L,C,D,inv,G,d,x,AreIsomorphic,bool;

L:=arg[1];
if Length(arg)>1 then bool:=true; else bool:=false; fi;

# Here L is a list of groups. The function returns a set of isomorphism
# class representatives for L. If bool=true then a possibly redundant list
# will be returned.

###########
inv:=function(G)
if Order(G)<1024 and not Order(G)=512 then return IdGroup(G); fi;
return G;
end;
###########

D:=Classify(L,inv);

if bool then return List(D,c->c[1]);; fi;  #A possibly redundant 
                                           #list will be returned.

#############
AreIsomorphic:=function(H)
if IsomorphismGroups(G,H)=fail then return false;
else return true; fi;
end;
#############

C:=[];

for d in D do
x:=HAP_EquivalenceClasses(d,AreIsomorphic);
Append(C,x);
od;

return C;


end;
#############################################
#############################################

#############################################
#############################################
StemGroups:=function(G)
# Inputs a group G and returns a list of stem groups,
# one in each isoclinism class.
local S, ZS, DS, M, L, stems, lems, stems1, pos, pos1, K, SS, Aut,
      Inn,fn, x,n,m;

if Order(Epicentre(G))>1 then
Print("This function can only be applied to central factor groups.\n");
return fail; fi;

S:=SchurCover(G);
ZS:=Center(S);
DS:=DerivedSubgroup(S);
M:=Intersection(ZS,DS);
L:=SubgroupsSolvableGroup(M);

lems:=[];
stems:=[];
for K in L do
SS:=S/K;
if IdGroup(SS/Center(SS))=IdGroup(G) then
Add(stems,SS); Add(lems,K); fi;
od;

stems:=GroupIsomorphismRepresentatives(stems,"with possible redundancies");
stems:=IsoclinismClasses(stems);
fn:=function(x,y); return Order(x)<Order(y); end;
for x in stems do Sort(x,fn);  od;
stems:=List(stems,x->x[1]);
return stems;
end;
#############################################
#############################################
