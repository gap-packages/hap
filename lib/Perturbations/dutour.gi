HAP_GCOMPLEX_LIST:=0;

################################################
################################################
InstallGlobalFunction(ContractibleGcomplex,
function(groupname)
local
	C, 
        G, StabilizerGroups, Stabilizer,
        lnth,
        dims,Dimension,
        Boundary,
        boundaryList,
        Elts,
	Rot,Stab,
        RotSubGroups,Action, ActionRecord,
        lst,n,k,s,BI,SGN,tmp, LstEl , bool;


groupname:=Filtered(groupname,x->not(x='(' or x=')' or x=',' or x='[' or x=']'));
groupname:=Concatenation("lib/Perturbations/Gcomplexes/",groupname);
bool:=ReadPackage("HAP",groupname);

if bool = false then 
Print("G-complexes are implemented for the following groups only:  \n\n");

lst:="  SL(3,Z) , PGL(3,Z[i]) , PGL(3,Eisenstein_Integers) , \n PSL(4,Z) , PSL(4,Z)_b , PSL(4,Z)_c , PSL(4,Z)_d \n"; 

Print(lst,"\n", "(where subscripts _b etc denote alternative complexes for a given group).\n");
return fail;
fi;

C:=StructuralCopy(HAP_GCOMPLEX_LIST);
lnth:=Length(C)-1;

dims:=List([1..lnth+1],n->Length(C[n]));

###################
Dimension:=function(n);
if n>lnth then return 0; fi;
return dims[n+1];
end;
###################

Elts:=[];
StabilizerGroups:=[];
RotSubGroups:=[];
boundaryList:=[];

#######
for n in [1..lnth+1] do
boundaryList[n]:=[];
StabilizerGroups[n]:=[];
RotSubGroups[n]:=[];
for k in [1..Dimension(n-1)] do
tmp:=C[n][k].BoundaryImage;
Stab:=C[n][k].TheMatrixStab;;
Stab:=List(GeneratorsOfGroup(Stab),w->TransposedMat(w));
C[n][k].TheMatrixStab:=Group(Stab);
Rot:=C[n][k].TheRotSubgroup;;
Rot:=List(GeneratorsOfGroup(Rot),w->TransposedMat(w));
C[n][k].TheRotSubgroup:=Group(Rot);;
Append(Elts,Elements(C[n][k].TheMatrixStab));
Add(StabilizerGroups[n],C[n][k].TheMatrixStab);
Add(RotSubGroups[n],C[n][k].TheRotSubgroup);
od;
od;
####

Elts:=SSortedList(Elts);

#######
for n in [1..lnth+1] do
boundaryList[n]:=[];
for k in [1..Dimension(n-1)] do
tmp:=C[n][k].BoundaryImage;
BI:=tmp.ListIFace;
SGN:=tmp.ListSign;
LstEl:=List(tmp.ListElt,w->TransposedMat(w));
Append(Elts,LstEl);
for s in [1..Length(BI)] do
BI[s]:=[SGN[s]*BI[s],Position(Elts,LstEl[s])];
od;
Add(boundaryList[n],BI);
od;
od;
####

G:=Group(Elts);


####################
Boundary:=function(n,k);
return boundaryList[n+1][k];
end;
####################

####################
Stabilizer:=function(n,k);
return StabilizerGroups[n+1][k];
end;
####################

####################
Action:=function(n,k,g)
local r,u,H;

H:=StabilizerGroups[n+1][AbsInt(k)];
r:=CanonicalRightCosetElement(H,Elts[g]^-1)^-1;
u:=Elts[g]^-1*r;
if u in RotSubGroups[n+1][AbsInt(k)] then return 1;
else return -1; fi;

end;
####################

return Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=fail,
            elts:=Elts,
            group:=G,
            stabilizer:=Stabilizer,
            action:=Action,
            properties:=
            [["length",Maximum(1000,lnth)],
             ["characteristic",0],
             ["type","resolution"],
             ["reduced",true]]  ));

end);
################################################
################################################



