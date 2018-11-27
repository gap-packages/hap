
InstallGlobalFunction(QuotientByTorsionSubcomplex,
function(C,p)
local pTorsionSubcomplex, celldata, quotient, torsionCells,
      i, j, k, t, l, currentCell, BoundaryImage, dim0, ReIndexCell,
      lnth, dims,n, Elts, Boundary, StabilizerGroups, RotSubGroups, s,
      boundaryList, BI, SGN, LstEl, tmp, G, Stabilizer, Action, Dimension;

    pTorsionSubcomplex:=GetTorsionSubcomplex(C,p);
    celldata:=StructuralCopy(pTorsionSubcomplex!.originalData);
    torsionCells:=[];
    for i in [1..Length(celldata)] do
          torsionCells[i]:=List(pTorsionSubcomplex!.torsionPosition[i],a->a[2]);
          Sort(torsionCells[i],function(v,w) return v>w;end);
    od;
#    Print(torsionCells,"\n");
    quotient:=StructuralCopy(celldata);
    i:=Length(quotient);
    while i > 0 do
        if Length(quotient[i])=0 then
            Remove(quotient,i);
            i:=i-1;
        else break;
        fi;
    od;
#    Print(quotient,"\n");
    if not Length(torsionCells[1])=0 then

    ## ** Delete torsion vertices ** ##
          for i in torsionCells[1] do
#              Print(torsionCells[1],"\n");
              Remove(quotient[1],i);
#              Print(quotient[1],"\n");
          od;
    ## ** Add a new vertex which is the collapsing point ** ##
          Add(quotient[1],rec(TheMatrixStab :=C!.group,
                      TheRotSubgroup:=C!.group,
                      BoundaryImage :=rec(
                            ListIFace:=[],
                            ListSign:=[],
                            ListElt:=[])
                      )
          );
#          Print(quotient[1],"\n");
          dim0:=Length(quotient[1]);
    ## ** Search for those edges that having torsion vertex and replace
    ## ** them with the collapsing vertex
          for j in torsionCells[2] do
              Remove(quotient[2],j);
          od;
          for i in [1..Length(quotient[2])] do
                currentCell:=quotient[2][i];
                BoundaryImage:=currentCell.BoundaryImage;
                for j in [1..Length(BoundaryImage.ListIFace)] do
                    if BoundaryImage.ListIFace[j] in torsionCells[1] then
                        BoundaryImage.ListIFace[j]:=dim0;
                        BoundaryImage.ListElt[j]:=One(C!.group);
                    fi;
                od;
                currentCell.BoundaryImage:=BoundaryImage;
                quotient[2][i]:=currentCell;
          od;
          for i in [3..Length(quotient)] do
    ## ** Remove the torsion i-cells
              for j in torsionCells[i] do
                  Remove(quotient[i],j);
              od;
    ## ** Search the boundary of an i-cell for the torsion (i-1)-cells
    ## ** and delete them
              for j in [1..Length(quotient[i])] do
                  currentCell:=quotient[i][j];
                  BoundaryImage:=currentCell.BoundaryImage;
                  l:=Length(BoundaryImage.ListIFace);
                  for k in [1..l] do
                      t:=l+1-k;
                      if BoundaryImage.ListIFace[t] in torsionCells[i-1] then
                          Remove(BoundaryImage.ListIFace,t);
                          Remove(BoundaryImage.ListSign,t);
                          Remove(BoundaryImage.ListElt,t);
                      fi;
                  od;
                  currentCell.BoundaryImage:=BoundaryImage;
                  quotient[i][j]:=currentCell;
              od;
          od;
    fi;

    ## ** A subfunction which reindexes cells after deleting torsion cells
    ReIndexCell:=function(k,currentCell)
    local i,j, BoundaryImage, x;
        BoundaryImage:=currentCell.BoundaryImage;
        for j in [1..Length(BoundaryImage.ListIFace)] do
            for i in torsionCells[k] do
                if BoundaryImage.ListIFace[j]>i then
                    BoundaryImage.ListIFace[j]:=BoundaryImage.ListIFace[j]-1;
                fi;
            od;
        od;
        currentCell.BoundaryImage:=BoundaryImage;
        return currentCell;
    end;
    ##

    for i in [2..Length(quotient)] do
        for j in [1..Length(quotient[i])] do
            quotient[i][j]:=ReIndexCell(i-1,quotient[i][j]);
        od;
    od;

    ############Taken from dutour.gi implemented by Graham Ellis#########
    lnth:=Length(quotient)-1;

    dims:=List([1..lnth+1],n->Length(quotient[n]));

    ###################
    Dimension:=function(n);
    if n>lnth then return 0; fi;
    return dims[n+1];
    end;
    ###################

    Elts:=[Identity(quotient[1][1].TheMatrixStab)];
    StabilizerGroups:=[];
    RotSubGroups:=[];
    boundaryList:=[];


    #######
    for n in [1..lnth+1] do
    boundaryList[n]:=[];
    StabilizerGroups[n]:=[];
    RotSubGroups[n]:=[];
      for k in [1..Dimension(n-1)] do
      #if not name in InfGrps then
#      Append(Elts,Elements(quotient[n][k].TheMatrixStab));
      #fi;
      Add(StabilizerGroups[n],quotient[n][k].TheMatrixStab);
      Add(RotSubGroups[n],quotient[n][k].TheRotSubgroup);
      od;
    od;
    ####

    Elts:=SSortedList(Elts);

    #######
    for n in [1..lnth+1] do
    boundaryList[n]:=[];
    for k in [1..Dimension(n-1)] do
    tmp:=quotient[n][k].BoundaryImage;
    BI:=tmp.ListIFace;
    SGN:=tmp.ListSign;
    LstEl:=StructuralCopy(tmp.ListElt);
    Append(Elts,Difference(LstEl,Elts));
    for s in [1..Length(BI)] do
    BI[s]:=[SGN[s]*BI[s],Position(Elts,LstEl[s])];
    od;
    Add(boundaryList[n],BI);
    od;
    od;
    ####

    G:=C!.group;

    ####################
    Boundary:=function(n,k);
    if k>0 then
    return boundaryList[n+1][k];
    else
    return NegateWord(boundaryList[n+1][-k]);
    fi;
    end;
    ####################

    ####################
    Stabilizer:=function(n,k);
    return StabilizerGroups[n+1][k];
    end;
    ####################


    ####################
    Action:=function(n,k,g)
    return 1;
    end;


    ##################################################################
return Objectify(HapNonFreeResolution,
    rec(
    quotientCells:=quotient,
    stabilizer:=Stabilizer,
    elts:=Elts,
    group:=G,
    length:=lnth,
    boundary:=Boundary,
    homotopy:=fail,
    action:=Action,
    dimension:=Dimension,
    properties:=
    [["length",Maximum(1000,lnth)],
    ["characteristic",0],
    ["type","resolution"]]  ));
end);


################### end of ControlledSubdivision ############################




InstallGlobalFunction( "ExtractTorsionSubcomplex", function(C, p)
#####################################################################
## Here, p is the prime for which to take the torsion subcomplex.  ##
## We extract the cells the stabilizer of which contains p-torsion.##
#####################################################################
local vcd, stabilizerCardinalities, celldata, data, torsionPosition,
torsionCells, numberOfTorsionCells, n, j, returnedData, warned,
groupname, admissibilityCheck, x, i, b, tmpCell, cell, boundary, groupName,
lnth, dims, Elts, Boundary, StabilizerGroups, RotSubGroups, s,k,
boundaryList, BI, SGN, LstEl, tmp, G, Stabilizer, Action, Dimension;

    admissibilityCheck := function(celldata)
    #########################################################
    ## A cell complex is admissible in the sense of Brown, ##
    ## if each cell stabilizer fixes its cell pointwise.   ##
    ## Additionally,				       ##
    ## we gather the cardinalities of the stabilizers.     ##
    #########################################################
    local stabilizerCardinalities, G, card, n, j, R, vcd, warned;
       warned := false;
       stabilizerCardinalities := [];
       vcd := Length(celldata)-1;

       for n in [0..vcd] do
	    stabilizerCardinalities[n+1] := [];
	    for j in [1..Length(celldata[n+1])] do
	       G :=   celldata[n+1][j]!.TheMatrixStab;
	       if IsFinite(G) then
	          card := Order(G);
	          stabilizerCardinalities[n+1][j] := card;
	          ## *** Now we have to compare              *** ##
	          ## *** with the order of "TheRotSubgroup"  *** ##
	          R := celldata[n+1][j]!.TheRotSubgroup;
	          if card > Order(R) and warned = false then
		    Print("****Warning: cell complex not admissible ",
			    "in the sense of Brown!****\n",
		    " Torsion subcomplex reduction requires cell subdivision.\n");
		    warned := true;
	          fi;
	       fi;
	    od;
       od;
       return [stabilizerCardinalities, warned];
   end;

   # Case 1: the input is a group name
   if IsString(C) then
   groupName:=C;
   groupname := Filtered( C, function ( x )
            return not (x = '(' or x = ')' or x = ',' or x = '[' or x = ']');
   end );
   Read(Concatenation( 	DirectoriesPackageLibrary("HAP")[1]![1],
			"Perturbations/Gcomplexes/",groupname));
   celldata := StructuralCopy(HAP_GCOMPLEX_LIST);

   # Case 2: the input is a variable
   else
       groupName:=fail;
       celldata:=[];
       i:=0;
       while C!.dimension(i) > 0 do
           cell:=[];
           for j in [1..C!.dimension(i)] do
               if not i=0 then
               boundary:=C!.boundary(i,j);
               Add(cell,rec(TheMatrixStab :=C!.stabilizer(i,j),
                           TheRotSubgroup:=C!.stabilizer(i,j),
                           BoundaryImage :=rec(
                                 ListIFace:=List(boundary,w->AbsInt(w[1])),
                                 ListSign:=List(boundary,w->SignInt(w[1])),
                                 ListElt:=List(boundary,w->C!.elts[w[2]])
                           )
                       )
               );
               else
               Add(cell,rec(TheMatrixStab :=C!.stabilizer(i,j),
                           TheRotSubgroup:=C!.stabilizer(i,j),
                           BoundaryImage :=rec(
                                 ListIFace:=[],
                                 ListSign:=[],
                                 ListElt:=[]
                           )
                       )
               );
               fi;
           od;
           Add(celldata,cell);
           i:=i+1;
       od;
   fi;
   vcd := Length(celldata) -1;
#   Print("Extracting the ",p,"-torsion subcomplex of the ",
#		vcd,"-dimensional ",groupName,"-cell complex ... \n");
   returnedData := admissibilityCheck(celldata);
   stabilizerCardinalities := returnedData[1];
   warned := returnedData[2];
   torsionCells := [];
   numberOfTorsionCells := [];
   for n in [0..vcd] do
	torsionCells[n+1] := [];
	numberOfTorsionCells[n+1] := 0;
	for j in [1..Length(celldata[n+1])] do
	   ## Check if the stabilizer contains p-torsion ##
	   if stabilizerCardinalities[n+1][j] mod p = 0 then
#		Print("Extracted ",n,"-cell numero ",j,
#			" of stabilizer cardinality ",
#			stabilizerCardinalities[n+1][j],".\n");
		numberOfTorsionCells[n+1]
			:= numberOfTorsionCells[n+1]+1;
	        torsionCells[n+1][numberOfTorsionCells[n+1]]
			:=[n, j];
	   fi;
	od;
   od;
#   return
#     [torsionCells, numberOfTorsionCells, celldata, stabilizerCardinalities, warned];
  data:=[];
  for i in [1..Length(torsionCells)] do
      data[i]:=[];
      for x in torsionCells[i] do
          Add(data[i],celldata[i][x[2]]);
      od;
  od;
for j in [2..Size(data)] do
  for i in [1..Size(data[j])] do
      tmpCell:=StructuralCopy(data[j][i]!.BoundaryImage);
      b:=List(tmpCell!.ListIFace,w->Position(torsionCells[j-1], [j-2,w]));
      tmpCell!.ListIFace:=b;
      data[j][i]!.BoundaryImage:=tmpCell;
  od;
od;
  torsionPosition:=torsionCells;
  torsionCells:=[];
  for i in [1..Size(data)] do
     torsionCells[i]:=[];
     for j in [1..Size(data[i])] do
         torsionCells[i][j]:=[i-1,j];
     od;
  od;
  ############################
  ############Taken from dutour.gi implemented by Graham Ellis#########
  lnth:=Length(data)-1;

  dims:=List([1..lnth+1],n->Length(data[n]));

  ###################
  Dimension:=function(n);
  if n>lnth then return 0; fi;
  return dims[n+1];
  end;
  ###################

  Elts:=[Identity(data[1][1].TheMatrixStab)];
  StabilizerGroups:=[];
  RotSubGroups:=[];
  boundaryList:=[];


  #######
  for n in [1..lnth+1] do
  boundaryList[n]:=[];
  StabilizerGroups[n]:=[];
  RotSubGroups[n]:=[];
    for k in [1..Dimension(n-1)] do
    #if not name in InfGrps then
    Append(Elts,Elements(data[n][k].TheMatrixStab));
    #fi;
    Add(StabilizerGroups[n],data[n][k].TheMatrixStab);
    Add(RotSubGroups[n],data[n][k].TheRotSubgroup);
    od;
  od;
  ####

  Elts:=SSortedList(Elts);

  #######
  for n in [1..lnth+1] do
  boundaryList[n]:=[];
  for k in [1..Dimension(n-1)] do
  tmp:=data[n][k].BoundaryImage;
  BI:=tmp.ListIFace;
  SGN:=tmp.ListSign;
  LstEl:=StructuralCopy(tmp.ListElt);
  Append(Elts,Difference(LstEl,Elts));
  for s in [1..Length(BI)] do
  BI[s]:=[SGN[s]*BI[s],Position(Elts,LstEl[s])];
  od;
  Add(boundaryList[n],BI);
  od;
  od;
  ####

  G:=C!.group;

  ####################
  Boundary:=function(n,k);
  if k>0 then
  return boundaryList[n+1][k];
  else
  return NegateWord(boundaryList[n+1][-k]);
  fi;
  end;
  ####################

  ####################
  Stabilizer:=function(n,k);
  return StabilizerGroups[n+1][k];
  end;
  ####################


  ####################
  Action:=function(n,k,g)
  return 1;
  end;



  #############################
  return Objectify(HapNonFreeResolution,
            rec(
            stabilizer:=Stabilizer,
            elts:=Elts,
            group:=G,
            length:=lnth,
            boundary:=Boundary,
            homotopy:=fail,
            action:=Action,
            dimension:=Dimension,
            torsion:=p,
            groupname:=groupName,
            torsionCells:=torsionCells,
            torsionPosition:=torsionPosition,
            originalData:=celldata,
            celldata:= data,
            numberOfTorsionCells:= numberOfTorsionCells,
            stabilizerCardinalities:= stabilizerCardinalities,
            warned:= warned,
            properties:=
            [["length",Maximum(1000,lnth)],
            ["characteristic",0],
            ["type","resolution"]]  ));

end);

###########################################################################
DeclareGlobalFunction("PrintTorsionSubcomplex");
InstallGlobalFunction("PrintTorsionSubcomplex",
function(R)
local N, reducedTorsionCells,J,n,j,G,B,b, ReduceModP, p, celldata, Stab,freq,
      t, StabList, StabListStr;

reducedTorsionCells:=R!.torsionCells;
celldata:=R!.celldata;
p:=R!.torsion;
#########
ReduceModP := function( G, p)
local K, h, H, Q, P;

	Q := G;
	K := ConjugacyClassesSubgroups(G);
	for h in K do
          if Size(h)=1 then
	       H := Representative(h);
	           if Size(H) > 1 then
	     	       if Gcd(Size(H),p) = 1 then
            	       # We pass to the quotient Q in the extension
		       # 1 -> H -> G -> Q -> 1,
		       # because H has trivial mod p cohomology.
		       # Our loop ends up using the normal subgroup H of maximal size,
		       # because the list K is ordered from small to large.
		           Q := Image( NaturalHomomorphismByNormalSubgroup( G,H ));
	              fi;
	           fi;
	    fi;
	od;
       P:=Q;
       if IsPNormal(Q,p) then P:=Normalizer(Q,Center(SylowSubgroup(Q,p)));fi;
  return P;
end;
#########

freq:=[];
StabList:=[];
for N in [1..Size(reducedTorsionCells)] do
  freq[N]:=[];
  for J in reducedTorsionCells[N] do
	n := J[1];
	j := J[2];
	G := celldata[n+1][j]!.TheMatrixStab;;

	Stab:=IdSmallGroup(ReduceModP(G,p));
  t:=Position(StabList,Stab);
  if t=fail then
      Add(StabList,Stab);
      freq[N][Length(StabList)]:=1;
  else
      if IsBound(freq[N][t]) then
          freq[N][t]:=freq[N][t]+1;
      else
          freq[N][t]:=1;
      fi;
  fi;
  od;

od;
StabListStr:=StructuralCopy(StabList);
for J in [1..Length(StabList)] do
    StabList[J]:=StructureDescription(SmallGroup(StabList[J]));
od;
for N in [1..Size(reducedTorsionCells)] do
    for J in [1..Length(StabList)] do
        if not IsBound(freq[N][J]) then
            freq[N][J]:=0;
        fi;
    od;

od;



return [StabListStr,StabList,freq];
end);



#################################
DeclareGlobalFunction("ConvertTorsionComplexToGcomplex");
InstallGlobalFunction("ConvertTorsionComplexToGcomplex",
function(C)
local data, Elts,n,lnth,dims,Dimension, StabilizerGroups,
      RotSubGroups,boundaryList,k, tmp, BI, SGN, LstEl,
      s, G, Boundary, Stabilizer, Action;

data:=StructuralCopy(C!.celldata);

lnth:=Length(data)-1;

dims:=List([1..lnth+1],n->Length(data[n]));

###################
Dimension:=function(n);
if n>lnth then return 0; fi;
return dims[n+1];
end;
###################

Elts:=[Identity(data[1][1].TheMatrixStab)];
StabilizerGroups:=[];
RotSubGroups:=[];
boundaryList:=[];


#######
for n in [1..lnth+1] do
boundaryList[n]:=[];
StabilizerGroups[n]:=[];
RotSubGroups[n]:=[];
  for k in [1..Dimension(n-1)] do
  #if not name in InfGrps then
  Append(Elts,Elements(data[n][k].TheMatrixStab));
  #fi;
  Add(StabilizerGroups[n],data[n][k].TheMatrixStab);
  Add(RotSubGroups[n],data[n][k].TheRotSubgroup);
  od;
od;
####

Elts:=SSortedList(Elts);

#######
for n in [1..lnth+1] do
boundaryList[n]:=[];
for k in [1..Dimension(n-1)] do
tmp:=data[n][k].BoundaryImage;
# Print(tmp,"\n");
BI:=tmp.ListIFace;
SGN:=tmp.ListSign;
LstEl:=StructuralCopy(tmp.ListElt);
Append(Elts,Difference(LstEl,Elts));
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
if k>0 then
return boundaryList[n+1][k];
else
return NegateWord(boundaryList[n+1][-k]);
fi;
end;
####################

####################
Stabilizer:=function(n,k);
return StabilizerGroups[n+1][k];
end;
####################


####################
Action:=function(n,k,g)
return 1;
end;



#############################
return Objectify(HapNonFreeResolution,
          rec(
          stabilizer:=Stabilizer,
          elts:=Elts,
          group:=G,
          length:=lnth,
          boundary:=Boundary,
          homotopy:=fail,
          action:=Action,
          dimension:=Dimension,
          torsion:=C!.torsion,
          groupname:=C!.groupname,
          torsionCells:=C!.torsionCells,
#          torsionPosition:=torsionPosition,
#          originalData:=celldata,
          celldata:= data,
#          numberOfTorsionCells:= numberOfTorsionCells,
#          stabilizerCardinalities:= stabilizerCardinalities,
#          warned:= warned,
          properties:=
          [["length",Maximum(1000,lnth)],
          ["characteristic",0],
          ["type","resolution"]]  ));


end);
