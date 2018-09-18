
#######################################################################
##
##  RigidFacetsSubdivision
##  Input: A non-rigid cell complex C, with finite cell stabilizers and such
##  that we could equip it with a metric such that each of the cells of C is
##  convex, the restriction of the metric to each cell is CAT(0)
##  and the cell stabilizers act by CAT(0) isometries.
##  An optional second argument allows to carry out the subdivision
##  only up to a given upper bound on the cell dimension.
##
##  Output: A rigidification of C.
##
####################################################################
DeclareGlobalFunction("ListsOfCellsToRegularCWComplex");
InstallGlobalFunction(ListsOfCellsToRegularCWComplex,
function(dim,boundaries,coboundaries,orientation)
local NrCells,Properties;

Properties:=[["dimension",dim]];

NrCells:=function(n);
if n>dim then return 0; fi;
return Length(Filtered(boundaries[n+1],x->not x[1]=0));
end;

return Objectify(HapRegularCWComplex,
       rec(
           nrCells:=NrCells,
           boundaries:=boundaries,
           coboundaries:=coboundaries,
           orientation:=orientation,
           vectorField:=fail,
           inverseVectorField:=fail,
           criticalCells:=fail,
           properties:=Properties));

end);
##################################################################
InstallGlobalFunction(RigidFacetsSubdivision,
function(arg)
local W, StabRec, i, j, N, M, x, bdry, s1, s2, p, k, w, t,
      DimRec, BoundaryRec, id, dims, NotRigid, NewCell,
      Cell, Elts, Boundary, Dimension, CLeftCosetElt, LCoset,
      pos, Stab, Mult,DimTemp, BoundaryTemp, Partition, ConnectedCheck,
      Stabilizer, Action, IsRigidCell, ReplaceCell, SubdividingCell,
      Orbit, C, NaiveEulerChar, CellHomology, exportFile, IsContractible,
      CosetRec, IsSimplyConnected, BoundaryMult, BoundaryMultRec, jj, jjj;

    C:=arg[1];
    i:=0;
    while C!.dimension(i)>0 do
        i:=i+1;
    od;
    N:=i-1; # Length of the chain complex
    M:=N;


if Length(arg)=2 then M:=Minimum(N,arg[2]);
fi;


    Elts:=C!.elts;
    StabRec:=[];
    DimRec:=[];
    CosetRec:=[];

    ##################################################################
    # If g in Elts return the position of g in the list,
    # otherwise, add g to Elts and return the position.
    pos:=function(g)
    local posit;

        posit:=Position(Elts,g);
        if posit=fail then
            Add(Elts,g);
            return Length(Elts);
        else
            return posit;
        fi;
    end;
    ##################################################################
    id:=pos(One(C!.group));
    ##################################################################
    # return the stabilizer of g*e,
    #
    Stab:=function(e,g)
    return ConjugateGroup(StabRec[e[1]+1][AbsInt(e[2])],Elts[g]^-1);
    end;
    ##################################################################
    # returns  a  "canonical"  representative  of  the  right  coset
    # Elts[g]*Stab[i+1][j]

    ######Add a record for left coset
    for i in [0..N] do
      CosetRec[i+1]:=[];
    od;

    CLeftCosetElt:=function(i,j,g)
      local p;
      if IsBound(CosetRec[AbsInt(i)+1][j]) then
        if IsBound(CosetRec[AbsInt(i)+1][j][g]) then
          return CosetRec[AbsInt(i)+1][j][g];
        else
          p:=pos(CanonicalRightCountableCosetElement
                                  (StabRec[AbsInt(i)+1][j],Elts[g]^-1)^-1);
          CosetRec[AbsInt(i)+1][j][g]:=p;
          return p;
        fi;
      else
        CosetRec[i+1][j]:=[];
        p:=pos(CanonicalRightCountableCosetElement
                                (StabRec[AbsInt(i)+1][j],Elts[g]^-1)^-1);
        CosetRec[AbsInt(i)+1][j][g]:=p;
        return p;
      fi;

    end;
    ##################################################################
    ##
    ##  Input:  A list L, degree k, position g of an element
    ##  Output: Product of g and L.
    ##
    Mult:=function(L,k,g)
    local x,w,t,h,y,vv,LL;
        vv:=[];
        LL:=ShallowCopy(L);
        for x in [1..Length(LL)] do
            w:=Elts[g]*Elts[LL[x][2]];
            t:=CLeftCosetElt(k,AbsInt(LL[x][1]),pos(w));
            Add(vv,[LL[x][1],t]);
        od;
        return vv;
    end;
    ###################################################################
    # Store essential data: stabilizers, boundaries, dimensions

    Orbit:=[];
    for i in [1..N+1] do
        Orbit[i]:=[];
    od;

    for i in [0..N] do
        StabRec[i+1]:=[];
        DimRec[i+1]:=C!.dimension(i);
        for j in [1..C!.dimension(i)] do
            StabRec[i+1][j]:=C!.stabilizer(i,j);
        od;
    od;

    BoundaryRec:=[];
    for i in [1..N] do
        BoundaryRec[i]:=[];
        for j in [1..DimRec[i+1]] do
            bdry:=C!.boundary(i,j);
            BoundaryRec[i][j]:=[];
            for x in bdry do
                s1:=C!.action(i-1,AbsInt(x[1]),x[2]);
                p:=pos(CanonicalRightCountableCosetElement
                            (C!.stabilizer(i-1,AbsInt(x[1])),Elts[x[2]]^-1)^-1);
                s2:=C!.action(i-1,AbsInt(x[1]),p);

                Add(BoundaryRec[i][j],[s1*s2*x[1],p]);
            od;

        od;
    od;

    #############
    BoundaryMultRec:=[];
    for jj in [1..N] do
      BoundaryMultRec[jj]:=[];
      for jjj in [1..Length(BoundaryRec[jj])] do
        BoundaryMultRec[jj][jjj]:=[];
      od;
    od;
    #############
    BoundaryMult:=function(i,j,g)
      local w;
      if IsBound(BoundaryMultRec[i]) then
        if IsBound(BoundaryMultRec[i][j]) then
          if IsBound(BoundaryMultRec[i][j][g]) then
            return BoundaryMultRec[i][j][g];
          else
            BoundaryMultRec[i][j][g]:=Mult(BoundaryRec[i][j],i-1,g);
            return BoundaryMultRec[i][j][g];
          fi;
        else
          BoundaryMultRec[i][j]:=[];
          BoundaryMultRec[i][j][g]:=Mult(BoundaryRec[i][j],i-1,g);
          return BoundaryMultRec[i][j][g];
        fi;
      else
        BoundaryMultRec[i]:=[];
        BoundaryMultRec[i][j]:=[];
        BoundaryMultRec[i][j][g]:=Mult(BoundaryRec[i][j],g);
        return BoundaryMultRec[i][j][g];
      fi;
    end;
    ##################################################################


    ##################################################################
    # Check whether the cell is rigid or not

    IsRigidCell:=function(k,m)
    local bdry, intst, L;
        bdry:=BoundaryRec[k][m];
        L:=List(bdry,w->Elements(ConjugateGroup(StabRec[k][AbsInt(w[1])],Elts[w[2]]^-1)));
        intst:=Intersection(L);
        if not Elements(StabRec[k+1][m])=Elements(intst) then
            return false;
        else return true;
        fi;

    end;
    ##################################################################


    ##################################################################
    # Subdividing (k+1)-cell number i (denoted by sigma) with the
    # Rigid Facets Subdivision Algorithm :
    SubdividingCell:=function(k,i)
    local bdry, w, x, d, y, a, b, T, j, j1,j2, z, t, c, Flag, s, s1, ConnectToCenter,
          SearchComponent,temp, l1, bdry1, bdry2, Orb, IsSameOrbit, OrbFlag, OrbElm,
          ConnectedCheck,tau,Ctau,F,exportBoundary,exportStabilizer,
          Boundaries,Stabilizers,Tcells, convexTcells,p;

	# Record the barycenter of sigma as a vertex, with stabilizer Gamma_sigma :
        DimRec[1]:=DimRec[1]+1;
        Add(StabRec[1],StabRec[k+1][i]);

        bdry:=ShallowCopy(BoundaryRec[k][i]);


        #############################################################
        # In each orbit {g_{jt} sigma_j with j >= 2,
	# choose one cell such that its union with the cells in T is connected,
	# and such that the naive Euler characteristic of this union remains 1.
  SearchComponent:=function(q)
  local x,z,j,w,TT, tau, boundaryOFtau, S, Ps, ps,countPrint, bool;
      x:=[];
      w:=[];
      Ps:=[[1,1]];

      j:=1;
      while not Length(T)=Length(Orb) do


          while j<=Length(Orb) do
              if OrbFlag[j]=0 then

                  while j1<=Length(Orb[j]) do
                      if Flag[j][j1]=0 then
                          TT:=Union(T,[Orb[j][j1]]);



                          if IsContractible(k-1,TT) then

                            if k>=2 then
                              S:=[];
                              for tau in TT do

                                   boundaryOFtau:=BoundaryMult(k-1,AbsInt(tau[1]),tau[2]);
                                if  tau[1]<0 then
                                    boundaryOFtau:=NegateWord(boundaryOFtau);
                                fi;
                                Append(S,boundaryOFtau);

                              od;
                              S:=AlgebraicReduction(S);

                              if (Length(T)=Length(Orb)-1) and k>=4 then bool:=true;
                              else bool:=false;
                              fi;
                              if IsSimplyConnected(k-2,S,bool) then

                                  Add(T,Orb[j][j1]);
                                  Flag[j][j1]:=1;
                                  Add(Ps,[j,j1]);
                                  OrbFlag[j]:=1;
                                  j:=1;
                                  break;
                              fi;
                            else
                              Add(T,Orb[j][j1]);
                              Flag[j][j1]:=1;
                              Add(Ps,[j,j1]);
                              OrbFlag[j]:=1;
                              j:=1;
                              break;

                            fi;

                          fi;

                      fi;

                      j1:=j1+1;
                    od;
                  if j1<=Length(Orb[j]) then
                      break;

                  fi;
              fi;
              j:=j+1;
              j1:=1;
          od;
          if j>Length(Orb) then
              Remove(T);

              ps:=Ps[Length(Ps)];
              OrbFlag[ps[1]]:=0;
              Flag[ps[1]][ps[2]]:=0;
              Remove(Ps);
              j:=ps[1];
              j1:=ps[2]+1;
              while j1>Length(Orb[j]) do
                ps:=Ps[Length(Ps)];
                OrbFlag[ps[1]]:=0;
                Flag[ps[1]][ps[2]]:=0;
                Remove(Ps);
                Remove(T);
                j:=ps[1];
                j1:=ps[2]+1;

              od;
              if Length(T)=0 then
                  Print("\n Could not find any fundamental domain for the ",i,"th of the",k,"-cells! \n");
                  Print(1/0);
              fi;
          fi;

      od;

      return T;
  end;
        #############################################################

    #################################################################
    # Connect the collection of (k-1)-cells T to the barycenter of the cell sigma.
    # The elements of T, as well as sigma, are in the format [k,i,gamma]:
    # dimension k, obtained by sending the ith orbit representative to its image under the element gamma (specified by its number in Gamma via C!.elts).
    # Input: e := [k-1, T].
    ConnectToCenter:=function(e)
    local bdry1, x, stab, bdrye, w, stablst, redbdry, tau, boundaryOFtau,
          LCoset, AddCell;


    ##################################################################
    # Add a k-cell with stabilizer stab and boundary bdry
    # to the cell complex
    AddCell:=function(m,stab,bdry,e)
    local a,j,g,s,w,u,v;

    if m<k then
        w:=e[1];

        for j in [1..DimTemp[m+1]] do
            u:=BoundaryTemp[m+1][j];
	    # The cells u and w are in the format [orbit number, gamma],
	    # such that gamma (pointing to C!.elts) sends the orbit representative to the cell.
	    # We make use of the common vertex, namely the barycenter of sigma,
	    # to test if there exists gamma in Gamma_sigma with gamma u = w.
            s:=DimRec[m+1]-DimTemp[m+1]+j;
            if AbsInt(w[1])=AbsInt(u[1]) then

            	# v:=StabRec[m][AbsInt(u[1]);
              for v in StabRec[m][AbsInt(u[1])] do
            	   g:=Elts[w[2]]*v*Elts[u[2]]^-1;
            	    if g in StabRec[k+1][i] then
                    return [s, CLeftCosetElt(m,s,pos(g))];
                  fi;
              od;

            fi;
        od;

        DimTemp[m+1]:=DimTemp[m+1]+1;

        Add(BoundaryTemp[m+1],e[1]);

	#Update the dimension :
        DimRec[m+1]:=DimRec[m+1]+1;

	# Record the stabilizer :
        Add(StabRec[m+1],stab);\

	# Record the boundary of F :
        Add(BoundaryRec[m],bdry);

    else
        DimRec[m+1]:=DimRec[m+1]+1;
        Add(StabRec[m+1],stab);
         Add(BoundaryRec[m],bdry);
    fi;
    return [DimRec[m+1],CLeftCosetElt(m,DimRec[m+1],id)];
    end;
    ##################################################################

        if e[1]=0 # k-1 = 0 : e[2] = T is a collection of vertices.
	then
            p:=Position(Tcells[1],e[2][1]);
            if not p=fail then return convexTcells[1][p];
            else
            bdry1:=[[-SignInt(e[2][1][1])*DimRec[1],CLeftCosetElt(0,DimRec[1],id)],[e[2][1][1],e[2][1][2]]];

            stab:=Intersection(Stab([0,e[2][1][1]],e[2][1][2]),StabRec[k+1][i]);
            w:=AddCell(1,stab,bdry1,e[2]);
            Add(Tcells[1],e[2][1]);
            Add(convexTcells[1],w);
            return w;
            fi;
        fi;

        stablst:=[];
        bdry1:=[];
        Append(bdry1,e[2]);
        bdrye:=[];

	## Enumerate the set S := boundary(geometric realization of T) =: bdrye,
	## by letting tau run through the cells of T = e[2],
	## and afterwards cancelling cells in S which appear in pairs with opposite orientation.
        for tau in e[2] do

            boundaryOFtau:=BoundaryMult(e[1],AbsInt(tau[1]),tau[2]);
            if  tau[1]<0 then boundaryOFtau:=NegateWord(boundaryOFtau);fi;
            Append(bdrye,boundaryOFtau);
            Add(stablst,Stab([e[1],tau[1]],tau[2]));
        od;
        Add(stablst,StabRec[k+1][i]);
	## Cancel cells in S which appear in pairs with opposite orientation :
        bdrye:=AlgebraicReduction(bdrye);
	## Construction of S = bdrye is now completed.

	## For every x in S, take the convex envelope w := e(x) of x and the barycenter of sigma :

        for x in bdrye do
          p:=Position(Tcells[e[1]],[AbsInt(x[1]),x[2]]);
          if not p=fail then
            w:=convexTcells[e[1]][p];
          else
            w:=ConnectToCenter([e[1]-1,[[AbsInt(x[1]),x[2]]]]);
            Add(Tcells[e[1]],[AbsInt(x[1]),x[2]]);
            Add(convexTcells[e[1]],w);
          fi;
            Add(bdry1,[-SignInt(x[1])*w[1],w[2]]);
	    ## Here, the orientation of w := e(x) has been recorded in -SignInt(x[1]).
        od;

	# Calculate the stabilizer :
        stab:=Intersection(stablst);
        if not Length(e[2])=1 then

            w:=AddCell(e[1]+1,stab,bdry1,e[2]);

            return w;
        else
            return AddCell(e[1]+1,stab,bdry1,e[2]);
        fi;
    end;
    ##################################################################
    ## End of function ConnectToCenter, end of Algorithm 3.
    ##################################################################


    IsSameOrbit:=function(e,f)
    local s1,s;
        if AbsInt(e[1])=AbsInt(f[1]) then
            s:=StabRec[k+1][i];
            s1:=StabRec[k][AbsInt(e[1])];
            for x in s do
                if Elts[f[2]]^-1*x*Elts[e[2]] in s1 then
                    return SignInt(e[1])*SignInt(f[1])*pos(x);
                fi;
            od;
            return false;
         else
            return false;
         fi;
    end;

##### 2017 - LUXEMBOURG


##### Check if the fundamental domain formed by sub[1] is connected ####

ConnectedCheck:=function(F,n)
	local b, lst, count, flag, i, d, cnr, m;

	###

	lst:=[];
        flag:=[];
	b:=[];
	for i in [1..Length(F)] do
		flag[i]:=0;
    m:=BoundaryMult(n,AbsInt(F[i][1]),F[i][2]);
		b[i]:=Set(List([1..Length(m)],k->[AbsInt(m[k][1]),m[k][2]]));
	od;
	count:=0;
	lst:=Union(lst,b[1]);
	flag[1]:=1;
	count:=1;
	d:=1;
	while count<Length(F) do
		for cnr in [2..Length(F)] do
			if flag[cnr]=0 then
				if not IsEmpty(Intersection(lst,b[cnr])) then
					flag[cnr]:=1;
					count:=count+1;
					lst:=Union(lst,b[cnr]);
					break;
				fi;
			fi;
		od;
		d:=d+1;

		if not count=d then return false; fi;
	od;
	return true;
end;


############################ END CHECK##########################

    ############################################################################
    # Sort the (m-1)-faces of sigma into orbits under the action of Gamma_sigma


        Orb:=[];
        Orb[1]:=[];
        OrbElm:=[];
        OrbElm[1]:=[];
	## Start with the cells of boundary(sigma) = bdry[1] :
        Add(Orb[1],bdry[1]);
        Add(OrbElm[1],id);

        for j in [2..Length(bdry)] do
            t:=0;
            for j1 in [1..Length(Orb)] do
                s:=IsSameOrbit(bdry[j],Orb[j1][1]);
                if not (s=false) then
                    Add(Orb[j1],bdry[j]);
                    Add(OrbElm[j1],s);
                    break;
                fi;
                t:=t+1;
            od;
            if t=Length(Orb) then Add(Orb,[bdry[j]]);Add(OrbElm,[id]);fi;
        od;


     ##################################################################
     # "Divide the boundary into rigid parts in the big cell [k,i]"
     #
     # Choose a fundamental domain T for the boundary of sigma
     # under the action of Gamma_sigma,
     # and tessellate the boundary of sigma with copies of T.
     ##################################################################

        Flag:=[];
        for j in [1..Length(Orb)] do
            Flag[j]:=[];
            for j1 in [1..Length(Orb[j])] do
                Flag[j][j1]:=0;
            od;
        od;

        for j in [1..Length(Orb[1])] do
            Flag[1][j]:=1;
        od;

        OrbFlag:=[];
        for j1 in [1..Length(Orb)] do
               OrbFlag[j1]:=0;
        od;
        OrbFlag[1]:=1;
if (k<=N) then
        T :=[Orb[1][1]];
        # In each orbit {g_{jt} sigma_j with j >= 2,
	# choose one cell such that its union with the cells in T is connected,
	# and such that the naive Euler characteristic of this union remains 1.



        T := SearchComponent(T[1]);


        if not Length(T)=Length(Orb) then
            Print("    Can not find a contractible fundamental domain! \n");
            Print(1/0);

        fi;


     ###################################################################
        w:=[];
        BoundaryTemp:=[];
        DimTemp:=[];
        for j in [0..(k)] do
            BoundaryTemp[j+1]:=[];
            DimTemp[j+1]:=0;
        od;

	# Use Algorithm 3 to construct a fundamental domain F for sigma


        # Create a list of already-taken-envelop
        Tcells:=[];
        for j in [1..k] do
          Tcells[j]:=[];
        od;
        convexTcells:=[];
        for j in [1..k] do
          convexTcells[j]:=[];
        od;

        d:=ConnectToCenter([k-1,T]);

        for j in [1..Length(OrbElm[1])] do
            Add(w,[SignInt(OrbElm[1][j])*d[1],pos(Elts[AbsInt(OrbElm[1][j])]*Elts[d[2]])]);
        od;
        return w;
      else

        T :=List(Orb,i->i[1]);
     ###################################################################
        w:=[];
        BoundaryTemp:=[];
        DimTemp:=[];
        for j in [0..(k)] do
            BoundaryTemp[j+1]:=[];
            DimTemp[j+1]:=0;
        od;

  # Use Algorithm 3 to construct a fundamental domain F for sigma
        F:=[];
        Tcells:=[];
        for j in [1..k] do
          Tcells[j]:=[];
        od;
        convexTcells:=[];
        for j in [1..k] do
          convexTcells[j]:=[];
        od;
        for tau in T do
             Ctau:=ConnectToCenter([k-1,[tau]]);
             Add(F,Ctau);
        od;

        for j in [1..Length(OrbElm[1])] do
            Append(w,List(F,d->[SignInt(OrbElm[1][j])*d[1],pos(Elts[AbsInt(OrbElm[1][j])]*Elts[d[2]])]));
        od;

        return w;
      fi;
    end;
    ##################################################################
    # Replacing a cell by its subdivision
    ReplaceCell:=function(k,m)
    local i, j, p, w, x, bdry, y, ww;
        w:=ShallowCopy(SubdividingCell(k,m));
        if k=N then
            Partition:=StructuralCopy(w);
            for i in [1..Length(Partition)] do
                Partition[i][1]:=AbsInt(Partition[i][1]-SignInt(Partition[i][1]));
            od;
        else Partition:=fail;
        fi;

        if k<=M and k<N then
        for i in [1..DimRec[k+2]] do
            bdry:=ShallowCopy(BoundaryRec[k+1][i]);
            p:=PositionsProperty(bdry,w->AbsInt(w[1])=m);
            for j in p do
                x:=bdry[j];
                ww:=ShallowCopy(w);
                if x[1]<0 then ww:=NegateWord(ww);fi;
                ww:=Mult(ww,k,x[2]);
                Append(bdry,ww);
            od;
            y:=bdry{p};
            bdry:=Set(bdry);
            SubtractSet(bdry,y);
            BoundaryRec[k+1][i]:=bdry;
        od;
        fi;
        BoundaryRec[k][m]:="del";
        StabRec[k+1][m]:="del";
    end;
    ##################################################################
    # Calculate the naive Euler characteristic for the list of k-cells L
    NaiveEulerChar:=function(k,L)
    local Cells, nrCells,x,w,j,s,id, Mat,t,M,b,p,h,r;

    Cells:=[];
    # id:=pos(One(C!.group));
    for j in [1..k+1] do
        Cells[j]:=[];

    od;

    j:=k;
    Append(Cells[k+1],L);


    while j>0 do
        for s in [1..Length(Cells[j+1])] do
            x:=Cells[j+1][s];
            w:=StructuralCopy(BoundaryRec[j][AbsInt(x[1])]);
            w:=Mult(w,j-1,x[2]);
            w:=List(w,a->[AbsInt(a[1]),a[2]]);
            Cells[j]:=Union(Cells[j],w);
        od;
        j:=j-1;
    od;
    nrCells:=List([1..k+1],a->Length(Cells[a]));
    w:=0;
    for j in [1..Length(nrCells)] do
        w:=w+(-1)^(j-1)*nrCells[j];
    od;

     return [w,nrCells];


    end;
    ################
    IsContractible:=function(k,L)
      local Boundaries,Coboundaries, Cells, x, w, b, j, s, t, Y, crit, l, flag,
            count, d, S, nrCells, eulerChar;



      # Construct the list of cells and the corresponding boundary and
      # coboundary of those cells
      Cells:=[];

      for j in [1..k+1] do
          Cells[j]:=[];

      od;
      j:=k;
      Append(Cells[k+1],L);
      Boundaries:=[];

      while j>0 do
          w:=[];
          for s in [1..Length(Cells[j+1])] do
              x:=Cells[j+1][s];
              w[s]:=BoundaryMult(j,AbsInt(x[1]),x[2]);
              w[s]:=List(w[s],a->[AbsInt(a[1]),a[2]]);
              Cells[j]:=Union(Cells[j],w[s]);
          od;
          if j=k then
              flag:=[];
              l:=Length(Cells[j+1]);
              for s in [1..l] do
                flag[s]:=0;
              od;
              flag[1]:=1;
              count:=1;
              S:=w[1];
              d:=1;
              while count<l do
                for s in [2..l] do
                  if flag[s]=0 then
                    if not IsEmpty(Intersection(S,w[s])) then
                      flag[s]:=1;
                      count:=count+1;
                      S:=Union(S,w[s]);
                      break;
                    fi;
                  fi;
                od;
                d:=d+1;
                if not count=d then return false;fi;
              od;
          fi;

          Boundaries[j+1]:=[];
          for s in [1..Length(Cells[j+1])] do
                  b:=List(w[s],a->Position(Cells[j],a));
                  Boundaries[j+1][s]:=Concatenation([Length(b)],b);
          od;

          j:=j-1;
      od;
      nrCells:=List([1..k+1],i->Length(Cells[i]));
      # Naive Euler Char of L=Union of T and \tau
      eulerChar:=Sum(List([1..k+1],i->(-1)^(i+1)*nrCells[i]));
      if not eulerChar=1 then return false;fi;

      Boundaries[1]:=List(Cells[1],a->[1,0]); # denote the boundary of 0-cell
      Boundaries[k+2]:=[];
      ######BOUNDARIES END########
      ######COBOUNDARIES BEGIN####
      Coboundaries:=[];; #Coboundaries[n+1] contains the info on n-cells.
      for j in [0..k] do
#        t:=j+3;
        Coboundaries[j+1]:=List(Boundaries[j+1],i->[0]);
        for s in [1..Length(Boundaries[j+2])] do
          b:=Boundaries[j+2][s];
          t:=Length(b);
          for i in b{[2..t]} do
            Coboundaries[j+1][i][1]:=Coboundaries[j+1][i][1]+1;
            Add(Coboundaries[j+1][i],s);
          od;
        od;

      od;
      Coboundaries[k+1]:=List(Boundaries[k+1],a->[0]);
      #####COBOUNDARIES END######
      Y:=ListsOfCellsToRegularCWComplex(k,Boundaries,Coboundaries,fail);
      crit:=CriticalCellsOfRegularCWComplex(Y);
      if (Length(crit)=1 and crit[1][1]=0) then return true;
      else return false;
      fi;

    end;
    ################
    IsSimplyConnected:=function(k,L,bool)
      local Boundaries,Coboundaries, Cells, x, w, b, j, s, t, Y, crit, l, flag,
            count, d, S, nrCells, eulerChar, P, Orientation;
      # Construct the list of cells and the corresponding boundary and
      # coboundary of those cells
      Cells:=[];

      for j in [1..k+1] do
          Cells[j]:=[];

      od;
      j:=k;
      Append(Cells[k+1],L);
      Boundaries:=[];
      Orientation:=[];

      while j>0 do
          w:=[];
          Orientation[j+1]:=[];
          for s in [1..Length(Cells[j+1])] do
              x:=Cells[j+1][s];
              w[s]:=BoundaryMult(j,AbsInt(x[1]),x[2]);
              Orientation[j+1][s]:=List(w[s],a->SignInt(a[1]));
              w[s]:=List(w[s],a->[AbsInt(a[1]),a[2]]);
              Cells[j]:=Union(Cells[j],w[s]);
          od;

          Boundaries[j+1]:=[];

          for s in [1..Length(Cells[j+1])] do
                  b:=List(w[s],a->Position(Cells[j],a));

                  Boundaries[j+1][s]:=Concatenation([Length(b)],b);
          od;


          j:=j-1;
      od;
      nrCells:=List([1..k+1],i->Length(Cells[i]));
      # Naive Euler Char of the boundary S of L
      eulerChar:=Sum(List([1..k+1],i->(-1)^(i+1)*nrCells[i]));
      if not eulerChar=1+(-1)^k then return false;fi;


      if bool=true then # bool=true if dimension of cells in S >=2 and
                        # Length(T)=Length(Orb)-1 so that the function is
                        # seeking the last cell in T


      Boundaries[1]:=List(Cells[1],a->[1,0]); # denote the boundary of 0-cell
      Orientation[1]:=List(Cells[1],a->[1]);
      Boundaries[k+2]:=[];
      ######BOUNDARIES END########
      ######COBOUNDARIES BEGIN####
      Coboundaries:=[];; #Coboundaries[n+1] contains the info on n-cells.
      for j in [0..k] do

        Coboundaries[j+1]:=List(Boundaries[j+1],i->[0]);
        for s in [1..Length(Boundaries[j+2])] do
          b:=Boundaries[j+2][s];
          t:=Length(b);
          for i in b{[2..t]} do
            Coboundaries[j+1][i][1]:=Coboundaries[j+1][i][1]+1;
            Add(Coboundaries[j+1][i],s);
          od;
        od;

      od;
      Coboundaries[k+1]:=List(Boundaries[k+1],a->[0]);
      #####COBOUNDARIES END######
      Y:=ListsOfCellsToRegularCWComplex(k,Boundaries,Coboundaries,Orientation);

      P:=FundamentalGroupOfRegularCWComplex(Y);

      if (IsTrivial(P)) then return true;
      else return false;
      fi;
    fi;
    return true;

    end;
    ###############
    CellHomology:=function(k,L)
    local Cells, nrCells,x,w,j,s,id, Mat,t,M,b,p,h,r;

    Cells:=[];

    for j in [1..k+1] do
        Cells[j]:=[];

    od;

    j:=k;
    Append(Cells[k+1],L);


    while j>0 do
        for s in [1..Length(Cells[j+1])] do
            x:=Cells[j+1][s];

            w[s]:=BoundaryMult(j,AbsInt(x[1]),x[2]);
            w:=List(w,a->[AbsInt(a[1]),a[2]]);
            Cells[j]:=Union(Cells[j],w);
        od;
        j:=j-1;
    od;
    nrCells:=List([1..k+1],a->Length(Cells[a]));

  Mat:=[];
  for j in [2..k+1] do
    M:=[];
    for s in [1..Length(Cells[j])] do
       M[s]:=[];
       for t in [1..Length(Cells[j-1])] do
          M[s][t]:=0;
       od;
    od;

    for s in [1..Length(Cells[j])] do
        x:=Cells[j][s];

        b:=BoundaryMult(j-1,AbsInt(x[1]),x[2]);
        for t in [1..Length(b)] do
            p:=Position(Cells[j-1],[AbsInt(b[t][1]),b[t][2]]);
            M[s][p]:=SignInt(b[t][1]);
        od;
    od;
    Mat[j]:=SmithNormalFormIntegerMat(TransposedMat(M));


  od;
    w:=[];
    # Create Mat[1] for d_0
    M:=[];
    for s in [1..Length(Cells[1])] do
        M[s]:=0;
    od;
    Mat[1]:=TransposedMat([M]);

    # Create Mat[k+1] do d_(k+1)
    M:=[];
    for s in [1..Length(Cells[k+1])] do
        M[s]:=0;
    od;
    Mat[k+2]:=[M];

    for j in [1..k+1] do
        h:=[];
        r:=Length(Cells[j])-RankMat(Mat[j])-RankMat(Mat[j+1]);
        for s in [1..r] do
	    h[s]:=0;
        od;
        for s in [1..RankMat(Mat[j+1])] do
            if not Mat[j+1][s][s]=1 then
                Add(h,Mat[j+1][s][s]);
            fi;
        od;
        Add(w,h);

    od;



     return w;


    end;



    ##################################################################

    ##################################################################
    # Main part: subdividing the fundamental domain
    NotRigid:=[];
    dims:=ShallowCopy(DimRec);
    i:=1;

    while i<=M do

        j:=1;
        while j<=dims[i+1] do
            if not IsRigidCell(i,j) then
                Add(NotRigid,[i,j]);
            fi;
            j:=j+1;
        od;

        i:=i+1;

    od;
    for x in NotRigid do
#       Apply Rigid Facets Subdivision to the cell x :
        ReplaceCell(x[1],x[2]);
    od;


    #Delete cells which are already replaced by its subdivision

    t:=1;
    for w in [1..Length(NotRigid)] do
        k:=NotRigid[w][1];
        j:=NotRigid[w][2];
        if k<N then
        for i in [1..DimRec[k+2]] do
            bdry:=BoundaryRec[k+1][i];

            if not IsString(bdry) then
               for x in bdry do
                    if AbsInt(x[1])>j then
                        x[1]:=x[1]-SignInt(x[1]);
                    fi;
                od;
             fi;
             BoundaryRec[k+1][i]:=bdry;

         od;
         fi;
         dims[k+1]:=dims[k+1]-1;
         DimRec[k+1]:=DimRec[k+1]-1;
         Remove(BoundaryRec[k],j);
         Remove(StabRec[k+1],j);
         if IsBound(NotRigid[w+1]) and NotRigid[w+1][1]=NotRigid[w][1] then
             NotRigid[w+1][2]:=NotRigid[w+1][2]-t;
             t:=t+1;
         else
             t:=1;
         fi;

    od;

    ##################################################################
    Boundary:=function(k,m)
        if m<0 then return NegateWord(BoundaryRec[k][AbsInt(m)]);fi;
        return BoundaryRec[k][m];
    end;

    Stabilizer:=function(k,m)
        return StabRec[k+1][m];
    end;

    Dimension:=function(k)
        if k>N then return 0;fi;
        return DimRec[k+1];
    end;

    Action:=function(k,i,j)
        return 1;
    end;
    ##################################################################
return Objectify(HapNonFreeResolution,
    rec(
    dimension:=Dimension,
    Partition:=Partition,
    boundary:=Boundary,
    homotopy:=fail,
    elts:=Elts,
    group:=C!.group,
    stabilizer:=Stabilizer,
    action:=Action,
    subdividing:=SubdividingCell,
    replacecell:=ReplaceCell,
    isrigid:=IsRigidCell,
    properties:=
    [["length",Maximum(1000,N)],
    ["characteristic",0],
    ["type","resolution"]]  ));
end);


################### end of ControlledSubdivision ############################


################## HybridSubdivision #######################################
## Hybrid Subdivision combine the goods of both RFS and VSS. It will Use
## RFS to subdivide the cells in dimension N-1
## and continue using VSS to subdivide top dimensional cells.

DeclareGlobalFunction("HybridSubdivision");
InstallGlobalFunction(HybridSubdivision,
function(arg)
local W, StabRec, i, j, N, M, x, bdry, s1, s2, p, k, w, t,
      DimRec, BoundaryRec, id, dims, NotRigid, NewCell,
      Cell, Elts, Boundary, Dimension, CLeftCosetElt, LCoset,
      pos, Stab, Mult,DimTemp, BoundaryTemp, Partition, ConnectedCheck,
      Stabilizer, Action, IsRigidCell, ReplaceCell, SubdividingCell,
      Orbit, C, NaiveEulerChar, CellHomology, exportFile, IsContractible,
      CosetRec, IsSimplyConnected, BoundaryMult, BoundaryMultRec, jj, jjj;

    C:=arg[1];
    i:=0;
    while C!.dimension(i)>0 do
        i:=i+1;
    od;
    N:=i-1; # Length of the chain complex
    M:=N;

if Length(arg)=2 then M:=Minimum(N,arg[2]);
fi;


    Elts:=C!.elts;
    StabRec:=[];
    DimRec:=[];
    CosetRec:=[];

    ##################################################################
    # If g in Elts return the position of g in the list,
    # otherwise, add g to Elts and return the position.
    pos:=function(g)
    local posit;

        posit:=Position(Elts,g);
        if posit=fail then
            Add(Elts,g);
            return Length(Elts);
        else
            return posit;
        fi;
    end;
    ##################################################################
    id:=pos(One(C!.group));
    ##################################################################
    # return the stabilizer of g*e,
    #
    Stab:=function(e,g)
    return ConjugateGroup(StabRec[e[1]+1][AbsInt(e[2])],Elts[g]^-1);
    end;
    ##################################################################
    # returns  a  "canonical"  representative  of  the  right  coset
    # Elts[g]*Stab[i+1][j]

    ######Add a record for left coset
    for i in [0..N] do
      CosetRec[i+1]:=[];
    od;

    CLeftCosetElt:=function(i,j,g)
      local p;
      if IsBound(CosetRec[AbsInt(i)+1][j]) then
        if IsBound(CosetRec[AbsInt(i)+1][j][g]) then
          return CosetRec[AbsInt(i)+1][j][g];
        else
          p:=pos(CanonicalRightCountableCosetElement
                                  (StabRec[AbsInt(i)+1][j],Elts[g]^-1)^-1);
          CosetRec[AbsInt(i)+1][j][g]:=p;
          return p;
        fi;
      else
        CosetRec[i+1][j]:=[];
        p:=pos(CanonicalRightCountableCosetElement
                                (StabRec[AbsInt(i)+1][j],Elts[g]^-1)^-1);
        CosetRec[AbsInt(i)+1][j][g]:=p;
        return p;
      fi;

    end;
    ##################################################################
    ##
    ##  Input:  A list L, degree k, position g of an element
    ##  Output: Product of g and L.
    ##
    Mult:=function(L,k,g)
    local x,w,t,h,y,vv,LL;
        vv:=[];
        LL:=ShallowCopy(L);
        for x in [1..Length(LL)] do
            w:=Elts[g]*Elts[LL[x][2]];
            t:=CLeftCosetElt(k,AbsInt(LL[x][1]),pos(w));
            Add(vv,[LL[x][1],t]);
        od;
        return vv;
    end;
    ###################################################################
    # Store essential data: stabilizers, boundaries, dimensions

    Orbit:=[];
    for i in [1..N+1] do
        Orbit[i]:=[];
    od;

    for i in [0..N] do
        StabRec[i+1]:=[];
        DimRec[i+1]:=C!.dimension(i);
        for j in [1..C!.dimension(i)] do
            StabRec[i+1][j]:=C!.stabilizer(i,j);
        od;
    od;

    BoundaryRec:=[];
    for i in [1..N] do
        BoundaryRec[i]:=[];
        for j in [1..DimRec[i+1]] do
            bdry:=C!.boundary(i,j);
            BoundaryRec[i][j]:=[];
            for x in bdry do
                s1:=C!.action(i-1,AbsInt(x[1]),x[2]);
                p:=pos(CanonicalRightCountableCosetElement
                            (C!.stabilizer(i-1,AbsInt(x[1])),Elts[x[2]]^-1)^-1);
                s2:=C!.action(i-1,AbsInt(x[1]),p);

                Add(BoundaryRec[i][j],[s1*s2*x[1],p]);
            od;

        od;
    od;

    #############
    BoundaryMultRec:=[];
    for jj in [1..N] do
      BoundaryMultRec[jj]:=[];
      for jjj in [1..Length(BoundaryRec[jj])] do
        BoundaryMultRec[jj][jjj]:=[];
      od;
    od;
    #############
    BoundaryMult:=function(i,j,g)
      local w;
      if IsBound(BoundaryMultRec[i]) then
        if IsBound(BoundaryMultRec[i][j]) then
          if IsBound(BoundaryMultRec[i][j][g]) then
            return BoundaryMultRec[i][j][g];
          else
            BoundaryMultRec[i][j][g]:=Mult(BoundaryRec[i][j],i-1,g);
            return BoundaryMultRec[i][j][g];
          fi;
        else
          BoundaryMultRec[i][j]:=[];
          BoundaryMultRec[i][j][g]:=Mult(BoundaryRec[i][j],i-1,g);
          return BoundaryMultRec[i][j][g];
        fi;
      else
        BoundaryMultRec[i]:=[];
        BoundaryMultRec[i][j]:=[];
        BoundaryMultRec[i][j][g]:=Mult(BoundaryRec[i][j],g);
        return BoundaryMultRec[i][j][g];
      fi;
    end;
    ##################################################################


    ##################################################################
    # Check whether the cell is rigid or not

    IsRigidCell:=function(k,m)
    local bdry, intst, L;
        bdry:=BoundaryRec[k][m];
        L:=List(bdry,w->Elements(ConjugateGroup(StabRec[k][AbsInt(w[1])],Elts[w[2]]^-1)));
        intst:=Intersection(L);
        if not Elements(StabRec[k+1][m])=Elements(intst) then
            return false;
        else return true;
        fi;

    end;
    ##################################################################


    ##################################################################
    # Subdividing (k+1)-cell number i (denoted by sigma) with the
    # Rigid Facets Subdivision Algorithm :
    SubdividingCell:=function(k,i)
    local bdry, w, x, d, y, a, b, T, j, j1,j2, z, t, c, Flag, s, s1, ConnectToCenter,
          SearchComponent,temp, l1, bdry1, bdry2, Orb, IsSameOrbit, OrbFlag, OrbElm,
          ConnectedCheck,tau,Ctau,F,exportBoundary,exportStabilizer,
          Boundaries,Stabilizers,Tcells, convexTcells,p;

	# Record the barycenter of sigma as a vertex, with stabilizer Gamma_sigma :
        DimRec[1]:=DimRec[1]+1;
        Add(StabRec[1],StabRec[k+1][i]);

        bdry:=ShallowCopy(BoundaryRec[k][i]);


        #############################################################
        # In each orbit {g_{jt} sigma_j with j >= 2,
	# choose one cell such that its union with the cells in T is connected,
	# and such that the naive Euler characteristic of this union remains 1.
  SearchComponent:=function(q)
  local x,z,j,w,TT, tau, boundaryOFtau, S, Ps, ps,countPrint, bool;
      x:=[];
      w:=[];
      Ps:=[[1,1]];

      j:=1;
      while not Length(T)=Length(Orb) do


          while j<=Length(Orb) do
              if OrbFlag[j]=0 then

                  while j1<=Length(Orb[j]) do
                      if Flag[j][j1]=0 then
                          TT:=Union(T,[Orb[j][j1]]);



                          if IsContractible(k-1,TT) then

                            if k>=2 then
                              S:=[];
                              for tau in TT do

                                   boundaryOFtau:=BoundaryMult(k-1,AbsInt(tau[1]),tau[2]);
                                if  tau[1]<0 then
                                    boundaryOFtau:=NegateWord(boundaryOFtau);
                                fi;
                                Append(S,boundaryOFtau);

                              od;
                              S:=AlgebraicReduction(S);
                              if (Length(T)=Length(Orb)-1) and k>=4 then bool:=true;
                              else bool:=false;
                              fi;
                              if IsSimplyConnected(k-2,S,bool) then

                                  Add(T,Orb[j][j1]);
                                  Flag[j][j1]:=1;
                                  Add(Ps,[j,j1]);
                                  OrbFlag[j]:=1;
                                  j:=1;
                                  break;
                              fi;
                            else
                              Add(T,Orb[j][j1]);
                              Flag[j][j1]:=1;
                              Add(Ps,[j,j1]);
                              OrbFlag[j]:=1;
                              j:=1;
                              break;

                            fi;

                          fi;

                      fi;

                      j1:=j1+1;
                    od;
                  if j1<=Length(Orb[j]) then
                      break;

                  fi;
              fi;
              j:=j+1;
              j1:=1;
          od;
          if j>Length(Orb) then
              Remove(T);

              ps:=Ps[Length(Ps)];
              OrbFlag[ps[1]]:=0;
              Flag[ps[1]][ps[2]]:=0;
              Remove(Ps);
              j:=ps[1];
              j1:=ps[2]+1;
              while j1>Length(Orb[j]) do
                ps:=Ps[Length(Ps)];
                OrbFlag[ps[1]]:=0;
                Flag[ps[1]][ps[2]]:=0;
                Remove(Ps);
                Remove(T);
                j:=ps[1];
                j1:=ps[2]+1;

              od;
              if Length(T)=0 then
                  Print("\n Could not find any fundamental domain for the ",i,"th of the",k,"-cells! \n");
                  Print(1/0);
              fi;
          fi;

      od;

      return T;
  end;
        #############################################################

    #################################################################
    # Connect the collection of (k-1)-cells T to the barycenter of the cell sigma.
    # The elements of T, as well as sigma, are in the format [k,i,gamma]:
    # dimension k, obtained by sending the ith orbit representative to its image under the element gamma (specified by its number in Gamma via C!.elts).
    # Input: e := [k-1, T].
    ConnectToCenter:=function(e)
    local bdry1, x, stab, bdrye, w, stablst, redbdry, tau, boundaryOFtau,
          LCoset, AddCell;


    ##################################################################
    # Add a k-cell with stabilizer stab and boundary bdry
    # to the cell complex
    AddCell:=function(m,stab,bdry,e)
    local a,j,g,s,w,u,v;

    if m<k then
        w:=e[1];

        for j in [1..DimTemp[m+1]] do
            u:=BoundaryTemp[m+1][j];
	    # The cells u and w are in the format [orbit number, gamma],
	    # such that gamma (pointing to C!.elts) sends the orbit representative to the cell.
	    # We make use of the common vertex, namely the barycenter of sigma,
	    # to test if there exists gamma in Gamma_sigma with gamma u = w.
            s:=DimRec[m+1]-DimTemp[m+1]+j;
            if AbsInt(w[1])=AbsInt(u[1]) then


              for v in StabRec[m][AbsInt(u[1])] do
            	   g:=Elts[w[2]]*v*Elts[u[2]]^-1;
            	    if g in StabRec[k+1][i] then
                    return [s, CLeftCosetElt(m,s,pos(g))];
                  fi;
              od;

            fi;
        od;

        DimTemp[m+1]:=DimTemp[m+1]+1;

        Add(BoundaryTemp[m+1],e[1]);

	#Update the dimension :
        DimRec[m+1]:=DimRec[m+1]+1;

	# Record the stabilizer :
        Add(StabRec[m+1],stab);\

	# Record the boundary of F :
        Add(BoundaryRec[m],bdry);

    else
        DimRec[m+1]:=DimRec[m+1]+1;
        Add(StabRec[m+1],stab);
         Add(BoundaryRec[m],bdry);
    fi;
    return [DimRec[m+1],CLeftCosetElt(m,DimRec[m+1],id)];
    end;
    ##################################################################

        if e[1]=0 # k-1 = 0 : e[2] = T is a collection of vertices.
	then
            p:=Position(Tcells[1],e[2][1]);
            if not p=fail then return convexTcells[1][p];
            else
            bdry1:=[[-SignInt(e[2][1][1])*DimRec[1],CLeftCosetElt(0,DimRec[1],id)],[e[2][1][1],e[2][1][2]]];

            stab:=Intersection(Stab([0,e[2][1][1]],e[2][1][2]),StabRec[k+1][i]);
            w:=AddCell(1,stab,bdry1,e[2]);
            Add(Tcells[1],e[2][1]);
            Add(convexTcells[1],w);
            return w;
            fi;
        fi;

        stablst:=[];
        bdry1:=[];
        Append(bdry1,e[2]);
        bdrye:=[];

	## Enumerate the set S := boundary(geometric realization of T) =: bdrye,
	## by letting tau run through the cells of T = e[2],
	## and afterwards cancelling cells in S which appear in pairs with opposite orientation.
        for tau in e[2] do
            # boundaryOFtau:=Mult(BoundaryRec[e[1]][AbsInt(tau[1])],e[1]-1,tau[2]);
            boundaryOFtau:=BoundaryMult(e[1],AbsInt(tau[1]),tau[2]);
            if  tau[1]<0 then boundaryOFtau:=NegateWord(boundaryOFtau);fi;
            Append(bdrye,boundaryOFtau);
            Add(stablst,Stab([e[1],tau[1]],tau[2]));
        od;
        Add(stablst,StabRec[k+1][i]);
	## Cancel cells in S which appear in pairs with opposite orientation :
        bdrye:=AlgebraicReduction(bdrye);
	## Construction of S = bdrye is now completed.

	## For every x in S, take the convex envelope w := e(x) of x and the barycenter of sigma :

        for x in bdrye do
          p:=Position(Tcells[e[1]],[AbsInt(x[1]),x[2]]);
          if not p=fail then
            w:=convexTcells[e[1]][p];
          else
            w:=ConnectToCenter([e[1]-1,[[AbsInt(x[1]),x[2]]]]);
            Add(Tcells[e[1]],[AbsInt(x[1]),x[2]]);
            Add(convexTcells[e[1]],w);
          fi;
            Add(bdry1,[-SignInt(x[1])*w[1],w[2]]);
	    ## Here, the orientation of w := e(x) has been recorded in -SignInt(x[1]).
        od;


	# Calculate the stabilizer :
        stab:=Intersection(stablst);
        if not Length(e[2])=1 then

            w:=AddCell(e[1]+1,stab,bdry1,e[2]);

            return w;
        else
            return AddCell(e[1]+1,stab,bdry1,e[2]);
        fi;
    end;
    ##################################################################
    ## End of function ConnectToCenter, end of Algorithm 3.
    ##################################################################


    IsSameOrbit:=function(e,f)
    local s1,s;
        if AbsInt(e[1])=AbsInt(f[1]) then
            s:=StabRec[k+1][i];
            s1:=StabRec[k][AbsInt(e[1])];
            for x in s do
                if Elts[f[2]]^-1*x*Elts[e[2]] in s1 then
                    return SignInt(e[1])*SignInt(f[1])*pos(x);
                fi;
            od;
            return false;
         else
            return false;
         fi;
    end;

##### 2017 - LUXEMBOURG


##### Check if the fundamental domain formed by sub[1] is connected ####

ConnectedCheck:=function(F,n)
	local b, lst, count, flag, i, d, cnr, m;

	###

	lst:=[];
        flag:=[];
	b:=[];
	for i in [1..Length(F)] do
		flag[i]:=0;

    m:=BoundaryMult(n,AbsInt(F[i][1]),F[i][2]);
		b[i]:=Set(List([1..Length(m)],k->[AbsInt(m[k][1]),m[k][2]]));
	od;
	count:=0;
	lst:=Union(lst,b[1]);
	flag[1]:=1;
	count:=1;
	d:=1;
	while count<Length(F) do
		for cnr in [2..Length(F)] do
			if flag[cnr]=0 then
				if not IsEmpty(Intersection(lst,b[cnr])) then
					flag[cnr]:=1;
					count:=count+1;
					lst:=Union(lst,b[cnr]);
					break;
				fi;
			fi;
		od;
		d:=d+1;

		if not count=d then return false; fi;
	od;
	return true;
end;


############################ END CHECK##########################

    ############################################################################
    # Sort the (m-1)-faces of sigma into orbits under the action of Gamma_sigma
        Orb:=[];
        Orb[1]:=[];
        OrbElm:=[];
        OrbElm[1]:=[];
	## Start with the cells of boundary(sigma) = bdry[1] :
        Add(Orb[1],bdry[1]);
        Add(OrbElm[1],id);

        for j in [2..Length(bdry)] do
            t:=0;
            for j1 in [1..Length(Orb)] do
                s:=IsSameOrbit(bdry[j],Orb[j1][1]);
                if not (s=false) then
                    Add(Orb[j1],bdry[j]);
                    Add(OrbElm[j1],s);
                    break;
                fi;
                t:=t+1;
            od;
            if t=Length(Orb) then Add(Orb,[bdry[j]]);Add(OrbElm,[id]);fi;
        od;

     ##################################################################
     # "Divide the boundary into rigid parts in the big cell [k,i]"
     #
     # Choose a fundamental domain T for the boundary of sigma
     # under the action of Gamma_sigma,
     # and tessellate the boundary of sigma with copies of T.
     ##################################################################

        Flag:=[];
        for j in [1..Length(Orb)] do
            Flag[j]:=[];
            for j1 in [1..Length(Orb[j])] do
                Flag[j][j1]:=0;
            od;
        od;

        for j in [1..Length(Orb[1])] do
            Flag[1][j]:=1;
        od;

        OrbFlag:=[];
        for j1 in [1..Length(Orb)] do
               OrbFlag[j1]:=0;
        od;
        OrbFlag[1]:=1;
if (k<=N-1) then

        T :=[Orb[1][1]];
        # In each orbit {g_{jt} sigma_j with j >= 2,
	# choose one cell such that its union with the cells in T is connected,
	# and such that the naive Euler characteristic of this union remains 1.

        T := SearchComponent(T[1]);


        if not Length(T)=Length(Orb) then
            Print("    Can not find a contractible fundamental domain! \n");
            Print(1/0);
        else


        fi;


     ###################################################################
        w:=[];
        BoundaryTemp:=[];
        DimTemp:=[];
        for j in [0..(k)] do
            BoundaryTemp[j+1]:=[];
            DimTemp[j+1]:=0;
        od;

	# Use Algorithm 3 to construct a fundamental domain F for sigma

        # Create a list of already-taken-envelop
        Tcells:=[];
        for j in [1..k] do
          Tcells[j]:=[];
        od;
        convexTcells:=[];
        for j in [1..k] do
          convexTcells[j]:=[];
        od;

        d:=ConnectToCenter([k-1,T]);

        for j in [1..Length(OrbElm[1])] do
            Add(w,[SignInt(OrbElm[1][j])*d[1],pos(Elts[AbsInt(OrbElm[1][j])]*Elts[d[2]])]);
        od;
        return w;
      else

        T :=List(Orb,i->i[1]);
     ###################################################################
        w:=[];
        BoundaryTemp:=[];
        DimTemp:=[];
        for j in [0..(k)] do
            BoundaryTemp[j+1]:=[];
            DimTemp[j+1]:=0;
        od;

  # Use Algorithm 3 to construct a fundamental domain F for sigma
        F:=[];
        Tcells:=[];
        for j in [1..k] do
          Tcells[j]:=[];
        od;
        convexTcells:=[];
        for j in [1..k] do
          convexTcells[j]:=[];
        od;
        for tau in T do
             Ctau:=ConnectToCenter([k-1,[tau]]);
             Add(F,Ctau);
        od;

        for j in [1..Length(OrbElm[1])] do
            Append(w,List(F,d->[SignInt(OrbElm[1][j])*d[1],pos(Elts[AbsInt(OrbElm[1][j])]*Elts[d[2]])]));
        od;

        return w;
      fi;
    end;
    ##################################################################
    # Replacing a cell by its subdivision
    ReplaceCell:=function(k,m)
    local i, j, p, w, x, bdry, y, ww;
        w:=ShallowCopy(SubdividingCell(k,m));
        if k=N then
            Partition:=StructuralCopy(w);
            for i in [1..Length(Partition)] do
                Partition[i][1]:=AbsInt(Partition[i][1]-SignInt(Partition[i][1]));
            od;
        else Partition:=fail;
        fi;

        if k<=M and k<N then
        for i in [1..DimRec[k+2]] do
            bdry:=ShallowCopy(BoundaryRec[k+1][i]);
            p:=PositionsProperty(bdry,w->AbsInt(w[1])=m);
            for j in p do
                x:=bdry[j];
                ww:=ShallowCopy(w);
                if x[1]<0 then ww:=NegateWord(ww);fi;
                ww:=Mult(ww,k,x[2]);
                Append(bdry,ww);
            od;
            y:=bdry{p};
            bdry:=Set(bdry);
            SubtractSet(bdry,y);
            BoundaryRec[k+1][i]:=bdry;
        od;
        fi;
        BoundaryRec[k][m]:="del";
        StabRec[k+1][m]:="del";
    end;
    ##################################################################
    # Calculate the naive Euler characteristic for the list of k-cells L
    NaiveEulerChar:=function(k,L)
    local Cells, nrCells,x,w,j,s,id, Mat,t,M,b,p,h,r;

    Cells:=[];

    for j in [1..k+1] do
        Cells[j]:=[];

    od;

    j:=k;
    Append(Cells[k+1],L);


    while j>0 do
        for s in [1..Length(Cells[j+1])] do
            x:=Cells[j+1][s];
            w:=StructuralCopy(BoundaryRec[j][AbsInt(x[1])]);
            w:=Mult(w,j-1,x[2]);
            w:=List(w,a->[AbsInt(a[1]),a[2]]);
            Cells[j]:=Union(Cells[j],w);
        od;
        j:=j-1;
    od;
    nrCells:=List([1..k+1],a->Length(Cells[a]));
    w:=0;
    for j in [1..Length(nrCells)] do
        w:=w+(-1)^(j-1)*nrCells[j];
    od;



     return [w,nrCells];


    end;
    ################
    IsContractible:=function(k,L)
      local Boundaries,Coboundaries, Cells, x, w, b, j, s, t, Y, crit, l, flag,
            count, d, S, nrCells, eulerChar;



      # Construct the list of cells and the corresponding boundary and
      # coboundary of those cells
      Cells:=[];

      for j in [1..k+1] do
          Cells[j]:=[];

      od;
      j:=k;
      Append(Cells[k+1],L);
      Boundaries:=[];

      while j>0 do
          w:=[];
          for s in [1..Length(Cells[j+1])] do
              x:=Cells[j+1][s];

              w[s]:=BoundaryMult(j,AbsInt(x[1]),x[2]);
              w[s]:=List(w[s],a->[AbsInt(a[1]),a[2]]);
              Cells[j]:=Union(Cells[j],w[s]);
          od;
          if j=k then
              flag:=[];
              l:=Length(Cells[j+1]);
              for s in [1..l] do
                flag[s]:=0;
              od;
              flag[1]:=1;
              count:=1;
              S:=w[1];
              d:=1;
              while count<l do
                for s in [2..l] do
                  if flag[s]=0 then
                    if not IsEmpty(Intersection(S,w[s])) then
                      flag[s]:=1;
                      count:=count+1;
                      S:=Union(S,w[s]);
                      break;
                    fi;
                  fi;
                od;
                d:=d+1;
                if not count=d then return false;fi;
              od;
          fi;

          Boundaries[j+1]:=[];
          for s in [1..Length(Cells[j+1])] do
                  b:=List(w[s],a->Position(Cells[j],a));
                  Boundaries[j+1][s]:=Concatenation([Length(b)],b);
          od;


          j:=j-1;
      od;
      nrCells:=List([1..k+1],i->Length(Cells[i]));
      # Naive Euler Char of L=Union of T and \tau
      eulerChar:=Sum(List([1..k+1],i->(-1)^(i+1)*nrCells[i]));
      if not eulerChar=1 then return false;fi;
      #
      # # Naive Euler Char of the boundary S of L

      Boundaries[1]:=List(Cells[1],a->[1,0]); # denote the boundary of 0-cell
      Boundaries[k+2]:=[];
      ######BOUNDARIES END########
      ######COBOUNDARIES BEGIN####
      Coboundaries:=[];; #Coboundaries[n+1] contains the info on n-cells.
      for j in [0..k] do

        Coboundaries[j+1]:=List(Boundaries[j+1],i->[0]);
        for s in [1..Length(Boundaries[j+2])] do
          b:=Boundaries[j+2][s];
          t:=Length(b);
          for i in b{[2..t]} do
            Coboundaries[j+1][i][1]:=Coboundaries[j+1][i][1]+1;
            Add(Coboundaries[j+1][i],s);
          od;
        od;

      od;
      Coboundaries[k+1]:=List(Boundaries[k+1],a->[0]);
      #####COBOUNDARIES END######
      Y:=ListsOfCellsToRegularCWComplex(k,Boundaries,Coboundaries,fail);
      crit:=CriticalCellsOfRegularCWComplex(Y);

      if (Length(crit)=1 and crit[1][1]=0) then return true;
      else return false;
      fi;

    end;
    ################
    IsSimplyConnected:=function(k,L,bool)
      local Boundaries,Coboundaries, Cells, x, w, b, j, s, t, Y, crit, l, flag,
            count, d, S, nrCells, eulerChar, P, Orientation;



      # Construct the list of cells and the corresponding boundary and
      # coboundary of those cells
      Cells:=[];

      for j in [1..k+1] do
          Cells[j]:=[];

      od;
      j:=k;
      Append(Cells[k+1],L);
      Boundaries:=[];
      Orientation:=[];

      while j>0 do
          w:=[];
          Orientation[j+1]:=[];
          for s in [1..Length(Cells[j+1])] do
              x:=Cells[j+1][s];

              w[s]:=BoundaryMult(j,AbsInt(x[1]),x[2]);
              Orientation[j+1][s]:=List(w[s],a->SignInt(a[1]));
              w[s]:=List(w[s],a->[AbsInt(a[1]),a[2]]);
              Cells[j]:=Union(Cells[j],w[s]);
          od;

          Boundaries[j+1]:=[];

          for s in [1..Length(Cells[j+1])] do
                  b:=List(w[s],a->Position(Cells[j],a));

                  Boundaries[j+1][s]:=Concatenation([Length(b)],b);
          od;


          j:=j-1;
      od;
      nrCells:=List([1..k+1],i->Length(Cells[i]));

      eulerChar:=Sum(List([1..k+1],i->(-1)^(i+1)*nrCells[i]));
      if not eulerChar=1+(-1)^k then return false;fi;


      if bool=true then


      Boundaries[1]:=List(Cells[1],a->[1,0]); # denote the boundary of 0-cell
      Orientation[1]:=List(Cells[1],a->[1]);
      Boundaries[k+2]:=[];
      ######BOUNDARIES END########
      ######COBOUNDARIES BEGIN####
      Coboundaries:=[];; #Coboundaries[n+1] contains the info on n-cells.
      for j in [0..k] do

        Coboundaries[j+1]:=List(Boundaries[j+1],i->[0]);
        for s in [1..Length(Boundaries[j+2])] do
          b:=Boundaries[j+2][s];
          t:=Length(b);
          for i in b{[2..t]} do
            Coboundaries[j+1][i][1]:=Coboundaries[j+1][i][1]+1;
            Add(Coboundaries[j+1][i],s);
          od;
        od;

      od;
      Coboundaries[k+1]:=List(Boundaries[k+1],a->[0]);
      #####COBOUNDARIES END######
      Y:=ListsOfCellsToRegularCWComplex(k,Boundaries,Coboundaries,Orientation);

      P:=FundamentalGroupOfRegularCWComplex(Y);

      if (IsTrivial(P)) then return true;
      else return false;
      fi;
      fi;
      return true;
    end;
    ###############
    CellHomology:=function(k,L)
    local Cells, nrCells,x,w,j,s,id, Mat,t,M,b,p,h,r;

    Cells:=[];

    for j in [1..k+1] do
        Cells[j]:=[];

    od;

    j:=k;
    Append(Cells[k+1],L);


    while j>0 do
        for s in [1..Length(Cells[j+1])] do
            x:=Cells[j+1][s];

            w[s]:=BoundaryMult(j,AbsInt(x[1]),x[2]);
            w:=List(w,a->[AbsInt(a[1]),a[2]]);
            Cells[j]:=Union(Cells[j],w);
        od;
        j:=j-1;
    od;
    nrCells:=List([1..k+1],a->Length(Cells[a]));

  Mat:=[];
  for j in [2..k+1] do
    M:=[];
    for s in [1..Length(Cells[j])] do
       M[s]:=[];
       for t in [1..Length(Cells[j-1])] do
          M[s][t]:=0;
       od;
    od;

    for s in [1..Length(Cells[j])] do
        x:=Cells[j][s];

        b:=BoundaryMult(j-1,AbsInt(x[1]),x[2]);
        for t in [1..Length(b)] do
            p:=Position(Cells[j-1],[AbsInt(b[t][1]),b[t][2]]);
            M[s][p]:=SignInt(b[t][1]);
        od;
    od;
    Mat[j]:=SmithNormalFormIntegerMat(TransposedMat(M));


  od;
    w:=[];
    # Create Mat[1] for d_0
    M:=[];
    for s in [1..Length(Cells[1])] do
        M[s]:=0;
    od;
    Mat[1]:=TransposedMat([M]);

    # Create Mat[k+1] do d_(k+1)
    M:=[];
    for s in [1..Length(Cells[k+1])] do
        M[s]:=0;
    od;
    Mat[k+2]:=[M];

    for j in [1..k+1] do
        h:=[];
        r:=Length(Cells[j])-RankMat(Mat[j])-RankMat(Mat[j+1]);
        for s in [1..r] do
	    h[s]:=0;
        od;
        for s in [1..RankMat(Mat[j+1])] do
            if not Mat[j+1][s][s]=1 then
                Add(h,Mat[j+1][s][s]);
            fi;
        od;
        Add(w,h);

    od;



     return w;


    end;



    ##################################################################

    ##################################################################
    # Main part: subdividing the fundamental domain
    NotRigid:=[];
    dims:=ShallowCopy(DimRec);
    i:=1;

    while i<=M do

        j:=1;
        while j<=dims[i+1] do
            if not IsRigidCell(i,j) then
                Add(NotRigid,[i,j]);
            fi;
            j:=j+1;
        od;

        i:=i+1;

    od;
    for x in NotRigid do
#       Apply Rigid Facets Subdivision to the cell x :
        ReplaceCell(x[1],x[2]);
    od;


    #Delete cells which are already replaced by its subdivision

    t:=1;
    for w in [1..Length(NotRigid)] do
        k:=NotRigid[w][1];
        j:=NotRigid[w][2];
        if k<N then
        for i in [1..DimRec[k+2]] do
            bdry:=BoundaryRec[k+1][i];

            if not IsString(bdry) then
               for x in bdry do
                    if AbsInt(x[1])>j then
                        x[1]:=x[1]-SignInt(x[1]);
                    fi;
                od;
             fi;
             BoundaryRec[k+1][i]:=bdry;

         od;
         fi;
         dims[k+1]:=dims[k+1]-1;
         DimRec[k+1]:=DimRec[k+1]-1;
         Remove(BoundaryRec[k],j);
         Remove(StabRec[k+1],j);
         if IsBound(NotRigid[w+1]) and NotRigid[w+1][1]=NotRigid[w][1] then
             NotRigid[w+1][2]:=NotRigid[w+1][2]-t;
             t:=t+1;
         else
             t:=1;
         fi;

    od;

    ##################################################################
    Boundary:=function(k,m)
        if m<0 then return NegateWord(BoundaryRec[k][AbsInt(m)]);fi;
        return BoundaryRec[k][m];
    end;

    Stabilizer:=function(k,m)
        return StabRec[k+1][m];
    end;

    Dimension:=function(k)
        if k>N then return 0;fi;
        return DimRec[k+1];
    end;

    Action:=function(k,i,j)
        return 1;
    end;
    ##################################################################
return Objectify(HapNonFreeResolution,
    rec(
    dimension:=Dimension,
    Partition:=Partition,
    boundary:=Boundary,
    homotopy:=fail,
    elts:=Elts,
    group:=C!.group,
    stabilizer:=Stabilizer,
    action:=Action,
    subdividing:=SubdividingCell,
    replacecell:=ReplaceCell,
    isrigid:=IsRigidCell,
    properties:=
    [["length",Maximum(1000,N)],
    ["characteristic",0],
    ["type","resolution"]]  ));
end);



################### end of ControlledSubdivision ############################

################## End of Hybrid Subdivision ##############################
DeclareGlobalFunction("VirtuallySimplicialSubdivision");
InstallGlobalFunction(VirtuallySimplicialSubdivision,
function(arg)
local W, StabRec, i, j, N, M, x, bdry, s1, s2, p, k, w, t,
      DimRec, BoundaryRec, id, dims, NotRigid, NewCell,
      Cell, Elts, Boundary, Dimension, CLeftCosetElt, LCoset,
      pos, Stab, Mult,DimTemp, BoundaryTemp, Partition, ConnectedCheck,
      Stabilizer, Action, IsRigidCell, ReplaceCell, SubdividingCell,
      Orbit, C, NaiveEulerChar, CellHomology, exportFile, IsContractible,
      CosetRec;

    C:=arg[1];
    i:=0;
    while C!.dimension(i)>0 do
        i:=i+1;
    od;
    N:=i-1; # Length of the chain complex
    M:=N;




if Length(arg)=2 then M:=Minimum(N,arg[2]);
fi;


    Elts:=C!.elts;
    StabRec:=[];
    DimRec:=[];
    CosetRec:=[];

    ##################################################################
    # If g in Elts return the position of g in the list,
    # otherwise, add g to Elts and return the position.
    pos:=function(g)
    local posit;

        posit:=Position(Elts,g);
        if posit=fail then
            Add(Elts,g);
            return Length(Elts);
        else
            return posit;
        fi;
    end;
    ##################################################################
    id:=pos(One(C!.group));
    ##################################################################
    # return the stabilizer of g*e,
    #
    Stab:=function(e,g)
    return ConjugateGroup(StabRec[e[1]+1][AbsInt(e[2])],Elts[g]^-1);
    end;
    ##################################################################
    # returns  a  "canonical"  representative  of  the  right  coset
    # Elts[g]*Stab[i+1][j]

    ######Add a record for left coset
    for i in [0..N] do
      CosetRec[i+1]:=[];
    od;

    CLeftCosetElt:=function(i,j,g)
      local p;
      if IsBound(CosetRec[AbsInt(i)+1][j]) then
        if IsBound(CosetRec[AbsInt(i)+1][j][g]) then
          return CosetRec[AbsInt(i)+1][j][g];
        else
          p:=pos(CanonicalRightCountableCosetElement
                                  (StabRec[AbsInt(i)+1][j],Elts[g]^-1)^-1);
          CosetRec[AbsInt(i)+1][j][g]:=p;
          return p;
        fi;
      else
        CosetRec[i+1][j]:=[];
        p:=pos(CanonicalRightCountableCosetElement
                                (StabRec[AbsInt(i)+1][j],Elts[g]^-1)^-1);
        CosetRec[AbsInt(i)+1][j][g]:=p;
        return p;
      fi;

    end;
    ##################################################################
    ##
    ##  Input:  A list L, degree k, position g of an element
    ##  Output: Product of g and L.
    ##
    Mult:=function(L,k,g)
    local x,w,t,h,y,vv,LL;
        vv:=[];
        LL:=ShallowCopy(L);
        for x in [1..Length(LL)] do
            w:=Elts[g]*Elts[LL[x][2]];
            t:=CLeftCosetElt(k,AbsInt(LL[x][1]),pos(w));
            Add(vv,[LL[x][1],t]);
        od;
        return vv;
    end;
    ###################################################################
    # Store essential data: stabilizers, boundaries, dimensions

    Orbit:=[];
    for i in [1..N+1] do
        Orbit[i]:=[];
    od;

    for i in [0..N] do
        StabRec[i+1]:=[];
        DimRec[i+1]:=C!.dimension(i);
        for j in [1..C!.dimension(i)] do
            StabRec[i+1][j]:=C!.stabilizer(i,j);
        od;
    od;

    BoundaryRec:=[];
    for i in [1..N] do
        BoundaryRec[i]:=[];
        for j in [1..DimRec[i+1]] do
            bdry:=C!.boundary(i,j);
            BoundaryRec[i][j]:=[];
            for x in bdry do
                s1:=C!.action(i-1,AbsInt(x[1]),x[2]);
                p:=pos(CanonicalRightCountableCosetElement
                            (C!.stabilizer(i-1,AbsInt(x[1])),Elts[x[2]]^-1)^-1);
                s2:=C!.action(i-1,AbsInt(x[1]),p);

                Add(BoundaryRec[i][j],[s1*s2*x[1],p]);
            od;

        od;
    od;
    ##################################################################


    ##################################################################
    # Check whether the cell is rigid or not

    IsRigidCell:=function(k,m)
    local bdry, intst, L;
        bdry:=BoundaryRec[k][m];
        L:=List(bdry,w->Elements(ConjugateGroup(StabRec[k][AbsInt(w[1])],Elts[w[2]]^-1)));
        intst:=Intersection(L);
        if not Elements(StabRec[k+1][m])=Elements(intst) then
            return false;
        else return true;
        fi;

    end;
    ##################################################################


    ##################################################################
    # Subdividing (k+1)-cell number i (denoted by sigma) with the
    # Rigid Facets Subdivision Algorithm :
    SubdividingCell:=function(k,i)
    local bdry, w, x, d, y, a, b, T, j, j1,j2, z, t, c, Flag, s, s1, ConnectToCenter,
          SearchComponent,temp, l1, bdry1, bdry2, Orb, IsSameOrbit, OrbFlag, OrbElm,
          ConnectedCheck,tau,Ctau,F,exportBoundary,exportStabilizer,
          Boundaries,Stabilizers,Tcells, convexTcells,p;

	# Record the barycenter of sigma as a vertex, with stabilizer Gamma_sigma :
        DimRec[1]:=DimRec[1]+1;
        Add(StabRec[1],StabRec[k+1][i]);

        bdry:=ShallowCopy(BoundaryRec[k][i]);


        #############################################################
        # In each orbit {g_{jt} sigma_j with j >= 2,
	# choose one cell such that its union with the cells in T is connected,
	# and such that the naive Euler characteristic of this union remains 1.
  SearchComponent:=function(q)
  local x,z,j,w,TT, tau, boundaryOFtau, S, Ps, ps,countPrint;
      x:=[];
      w:=[];
      Ps:=[[1,1]];

      j:=1;
      while not Length(T)=Length(Orb) do
          while j<=Length(Orb) do
              if OrbFlag[j]=0 then

                  while j1<=Length(Orb[j]) do
                      if Flag[j][j1]=0 then
                          TT:=Union(T,[Orb[j][j1]]);
                          S:=[];
                          for tau in TT do
                              boundaryOFtau:=Mult(BoundaryRec[k-1]
                                 [AbsInt(tau[1])],k-2,tau[2]);
                              if  tau[1]<0 then
                                  boundaryOFtau:=NegateWord(boundaryOFtau);
                              fi;
                              Append(S,boundaryOFtau);
                              S:=AlgebraicReduction(S);

                          od;
                          if k<2 or (NaiveEulerChar(k-2,S)[1]=1+(-1)^(k-2)) then

                          if IsContractible(k-1,TT) then

                                  Add(T,Orb[j][j1]);
                                  Flag[j][j1]:=1;
                                  Add(Ps,[j,j1]);
                                  OrbFlag[j]:=1;
                                  j:=1;
                                  break;


                          fi;
                          fi;


                      fi;

                      j1:=j1+1;
                    od;
                  if j1<=Length(Orb[j]) then
                      break;

                  fi;
              fi;
              j:=j+1;
              j1:=1;
          od;
          if j>Length(Orb) then
              Remove(T);

              ps:=Ps[Length(Ps)];
              OrbFlag[ps[1]]:=0;
              Flag[ps[1]][ps[2]]:=0;
              Remove(Ps);
              j:=ps[1];
              j1:=ps[2]+1;
              while j1>Length(Orb[j]) do
                ps:=Ps[Length(Ps)];
                OrbFlag[ps[1]]:=0;
                Flag[ps[1]][ps[2]]:=0;
                Remove(Ps);
                Remove(T);
                j:=ps[1];
                j1:=ps[2]+1;

              od;
              if Length(T)=0 then
                  Print("\n Could not find any fundamental domain for the ",i,"th of the",k,"-cells! \n");

                  Print(1/0);
              fi;
          fi;

      od;

      return T;
  end;
        #############################################################

    #################################################################
    # Connect the collection of (k-1)-cells T to the barycenter of the cell sigma.
    # The elements of T, as well as sigma, are in the format [k,i,gamma]:
    # dimension k, obtained by sending the ith orbit representative to its image under the element gamma (specified by its number in Gamma via C!.elts).
    # Input: e := [k-1, T].
    ConnectToCenter:=function(e)
    local bdry1, x, stab, bdrye, w, stablst, redbdry, tau, boundaryOFtau,
          LCoset, AddCell;


    ##################################################################
    # Add a k-cell with stabilizer stab and boundary bdry
    # to the cell complex
    AddCell:=function(m,stab,bdry,e)
    local a,j,g,s,w,u,v;

    if m<k then
        w:=e[1];

        for j in [1..DimTemp[m+1]] do
            u:=BoundaryTemp[m+1][j];
	    # The cells u and w are in the format [orbit number, gamma],
	    # such that gamma (pointing to C!.elts) sends the orbit representative to the cell.
	    # We make use of the common vertex, namely the barycenter of sigma,
	    # to test if there exists gamma in Gamma_sigma with gamma u = w.
            s:=DimRec[m+1]-DimTemp[m+1]+j;
            if AbsInt(w[1])=AbsInt(u[1]) then

            	# v:=StabRec[m][AbsInt(u[1]);
              for v in StabRec[m][AbsInt(u[1])] do
            	   g:=Elts[w[2]]*v*Elts[u[2]]^-1;
            	    if g in StabRec[k+1][i] then
                    return [s, CLeftCosetElt(m,s,pos(g))];
                  fi;
              od;

            fi;
        od;

        DimTemp[m+1]:=DimTemp[m+1]+1;

        Add(BoundaryTemp[m+1],e[1]);

	#Update the dimension :
        DimRec[m+1]:=DimRec[m+1]+1;

	# Record the stabilizer :
        Add(StabRec[m+1],stab);\

	# Record the boundary of F :
        Add(BoundaryRec[m],bdry);

    else
        DimRec[m+1]:=DimRec[m+1]+1;
        Add(StabRec[m+1],stab);
         Add(BoundaryRec[m],bdry);
    fi;
    return [DimRec[m+1],CLeftCosetElt(m,DimRec[m+1],id)];
    end;
    ##################################################################

        if e[1]=0 # k-1 = 0 : e[2] = T is a collection of vertices.
	then
            p:=Position(Tcells[1],e[2][1]);
            if not p=fail then return convexTcells[1][p];
            else
            bdry1:=[[-SignInt(e[2][1][1])*DimRec[1],CLeftCosetElt(0,DimRec[1],id)],[e[2][1][1],e[2][1][2]]];

            stab:=Intersection(Stab([0,e[2][1][1]],e[2][1][2]),StabRec[k+1][i]);
            w:=AddCell(1,stab,bdry1,e[2]);
            Add(Tcells[1],e[2][1]);
            Add(convexTcells[1],w);
            return w;
            fi;
        fi;

        stablst:=[];
        bdry1:=[];
        Append(bdry1,e[2]);
        bdrye:=[];

	## Enumerate the set S := boundary(geometric realization of T) =: bdrye,
	## by letting tau run through the cells of T = e[2],
	## and afterwards cancelling cells in S which appear in pairs with opposite orientation.
        for tau in e[2] do
            boundaryOFtau:=Mult(BoundaryRec[e[1]][AbsInt(tau[1])],e[1]-1,tau[2]);
            if  tau[1]<0 then boundaryOFtau:=NegateWord(boundaryOFtau);fi;
            Append(bdrye,boundaryOFtau);
            Add(stablst,Stab([e[1],tau[1]],tau[2]));
        od;
        Add(stablst,StabRec[k+1][i]);
	## Cancel cells in S which appear in pairs with opposite orientation :
        bdrye:=AlgebraicReduction(bdrye);
	## Construction of S = bdrye is now completed.

	## For every x in S, take the convex envelope w := e(x) of x and the barycenter of sigma :

        for x in bdrye do
          p:=Position(Tcells[e[1]],[AbsInt(x[1]),x[2]]);
          if not p=fail then
            w:=convexTcells[e[1]][p];
          else
            w:=ConnectToCenter([e[1]-1,[[AbsInt(x[1]),x[2]]]]);
            Add(Tcells[e[1]],[AbsInt(x[1]),x[2]]);
            Add(convexTcells[e[1]],w);
          fi;
            Add(bdry1,[-SignInt(x[1])*w[1],w[2]]);
	    ## Here, the orientation of w := e(x) has been recorded in -SignInt(x[1]).
        od;
	# Calculate the stabilizer :
        stab:=Intersection(stablst);
        if not Length(e[2])=1 then

            w:=AddCell(e[1]+1,stab,bdry1,e[2]);

            return w;
        else
            return AddCell(e[1]+1,stab,bdry1,e[2]);
        fi;
    end;
    ##################################################################
    ## End of function ConnectToCenter, end of Algorithm 3.
    ##################################################################


    IsSameOrbit:=function(e,f)
    local s1,s;
        if AbsInt(e[1])=AbsInt(f[1]) then
            s:=StabRec[k+1][i];
            s1:=StabRec[k][AbsInt(e[1])];
            for x in s do
                if Elts[f[2]]^-1*x*Elts[e[2]] in s1 then
                    return SignInt(e[1])*SignInt(f[1])*pos(x);
                fi;
            od;
            return false;
         else
            return false;
         fi;
    end;

##### 2017 - LUXEMBOURG


##### Check if the fundamental domain formed by sub[1] is connected ####

ConnectedCheck:=function(F,n)
	local b, lst, count, flag, i, d, cnr, m;

	###

	lst:=[];
        flag:=[];
	b:=[];
	for i in [1..Length(F)] do
		flag[i]:=0;
		m:=Mult(BoundaryRec[n][AbsInt(F[i][1])],n-1,F[i][2]);
		b[i]:=Set(List([1..Length(m)],k->[AbsInt(m[k][1]),m[k][2]]));
	od;
	count:=0;
	lst:=Union(lst,b[1]);
	flag[1]:=1;
	count:=1;
	d:=1;
	while count<Length(F) do
		for cnr in [2..Length(F)] do
			if flag[cnr]=0 then
				if not IsEmpty(Intersection(lst,b[cnr])) then
					flag[cnr]:=1;
					count:=count+1;
					lst:=Union(lst,b[cnr]);
					break;
				fi;
			fi;
		od;
		d:=d+1;
		if not count=d then return false; fi;
	od;
	return true;
end;


############################ END CHECK##########################

    ############################################################################
    # Sort the (m-1)-faces of sigma into orbits under the action of Gamma_sigma
        Orb:=[];
        Orb[1]:=[];
        OrbElm:=[];
        OrbElm[1]:=[];
	## Start with the cells of boundary(sigma) = bdry[1] :
        Add(Orb[1],bdry[1]);
        Add(OrbElm[1],id);

        for j in [2..Length(bdry)] do
            t:=0;
            for j1 in [1..Length(Orb)] do
                s:=IsSameOrbit(bdry[j],Orb[j1][1]);
                if not (s=false) then
                    Add(Orb[j1],bdry[j]);
                    Add(OrbElm[j1],s);
                    break;
                fi;
                t:=t+1;
            od;
            if t=Length(Orb) then Add(Orb,[bdry[j]]);Add(OrbElm,[id]);fi;
        od;


     ##################################################################
     # "Divide the boundary into rigid parts in the big cell [k,i]"
     #
     # Choose a fundamental domain T for the boundary of sigma
     # under the action of Gamma_sigma,
     # and tessellate the boundary of sigma with copies of T.
     ##################################################################

        Flag:=[];
        for j in [1..Length(Orb)] do
            Flag[j]:=[];
            for j1 in [1..Length(Orb[j])] do
                Flag[j][j1]:=0;
            od;
        od;

        for j in [1..Length(Orb[1])] do
            Flag[1][j]:=1;
        od;

        OrbFlag:=[];
        for j1 in [1..Length(Orb)] do
               OrbFlag[j1]:=0;
        od;
        OrbFlag[1]:=1;

        T :=List(Orb,i->i[1]);
     ###################################################################
        w:=[];
        BoundaryTemp:=[];
        DimTemp:=[];
        for j in [0..(k)] do
            BoundaryTemp[j+1]:=[];
            DimTemp[j+1]:=0;
        od;

  # Use Algorithm 3 to construct a fundamental domain F for sigma
        F:=[];
        Tcells:=[];
        for j in [1..k] do
          Tcells[j]:=[];
        od;
        convexTcells:=[];
        for j in [1..k] do
          convexTcells[j]:=[];
        od;
        for tau in T do
             Ctau:=ConnectToCenter([k-1,[tau]]);
             Add(F,Ctau);
        od;

        for j in [1..Length(OrbElm[1])] do
            Append(w,List(F,d->[SignInt(OrbElm[1][j])*d[1],pos(Elts[AbsInt(OrbElm[1][j])]*Elts[d[2]])]));
        od;

        return w;

    end;
    ##################################################################
    # Replacing a cell by its subdivision
    ReplaceCell:=function(k,m)
    local i, j, p, w, x, bdry, y, ww;
        w:=ShallowCopy(SubdividingCell(k,m));
        if k=N then
            Partition:=StructuralCopy(w);
            for i in [1..Length(Partition)] do
                Partition[i][1]:=AbsInt(Partition[i][1]-SignInt(Partition[i][1]));
            od;
        else Partition:=fail;
        fi;

        if k<=M and k<N then
        for i in [1..DimRec[k+2]] do
            bdry:=ShallowCopy(BoundaryRec[k+1][i]);
            p:=PositionsProperty(bdry,w->AbsInt(w[1])=m);
            for j in p do
                x:=bdry[j];
                ww:=ShallowCopy(w);
                if x[1]<0 then ww:=NegateWord(ww);fi;
                ww:=Mult(ww,k,x[2]);
                Append(bdry,ww);
            od;
            y:=bdry{p};
            bdry:=Set(bdry);
            SubtractSet(bdry,y);
            BoundaryRec[k+1][i]:=bdry;
        od;
        fi;
        BoundaryRec[k][m]:="del";
        StabRec[k+1][m]:="del";
    end;
    ##################################################################
    # Calculate the naive Euler characteristic for the list of k-cells L
    NaiveEulerChar:=function(k,L)
    local Cells, nrCells,x,w,j,s,id, Mat,t,M,b,p,h,r;

    Cells:=[];
    # id:=pos(One(C!.group));
    for j in [1..k+1] do
        Cells[j]:=[];

    od;

    j:=k;
    Append(Cells[k+1],L);


    while j>0 do
        for s in [1..Length(Cells[j+1])] do
            x:=Cells[j+1][s];
            w:=StructuralCopy(BoundaryRec[j][AbsInt(x[1])]);
            w:=Mult(w,j-1,x[2]);
            w:=List(w,a->[AbsInt(a[1]),a[2]]);
            Cells[j]:=Union(Cells[j],w);
        od;
        j:=j-1;
    od;
    nrCells:=List([1..k+1],a->Length(Cells[a]));
    w:=0;
    for j in [1..Length(nrCells)] do
        w:=w+(-1)^(j-1)*nrCells[j];
    od;

     return [w,nrCells];


    end;
    ################
    IsContractible:=function(k,L)
      local Boundaries,Coboundaries, Cells, x, w, b, j, s, t, Y, crit, l, flag,
            count, d, S, nrCells, eulerChar;



      # Construct the list of cells and the corresponding boundary and
      # coboundary of those cells
      Cells:=[];

      for j in [1..k+1] do
          Cells[j]:=[];

      od;
      j:=k;
      Append(Cells[k+1],L);
      Boundaries:=[];

      while j>0 do
          w:=[];
          for s in [1..Length(Cells[j+1])] do
              x:=Cells[j+1][s];
              w[s]:=StructuralCopy(BoundaryRec[j][AbsInt(x[1])]);
              w[s]:=Mult(w[s],j-1,x[2]);
              w[s]:=List(w[s],a->[AbsInt(a[1]),a[2]]);
              Cells[j]:=Union(Cells[j],w[s]);
          od;
          if j=k then
              flag:=[];
              l:=Length(Cells[j+1]);
              for s in [1..l] do
                flag[s]:=0;
              od;
              flag[1]:=1;
              count:=1;
              S:=w[1];
              d:=1;
              while count<l do
                for s in [2..l] do
                  if flag[s]=0 then
                    if not IsEmpty(Intersection(S,w[s])) then
                      flag[s]:=1;
                      count:=count+1;
                      S:=Union(S,w[s]);
                      break;
                    fi;
                  fi;
                od;
                d:=d+1;
                if not count=d then return false;fi;
              od;
          fi;

          Boundaries[j+1]:=[];
          for s in [1..Length(Cells[j+1])] do
                  b:=List(w[s],a->Position(Cells[j],a));
                  Boundaries[j+1][s]:=Concatenation([Length(b)],b);
          od;


          j:=j-1;
      od;
      nrCells:=List([1..k+1],i->Length(Cells[i]));
      # Naive Euler Char of L=Union of T and \tau
      eulerChar:=Sum(List([1..k+1],i->(-1)^(i+1)*nrCells[i]));
      if not eulerChar=1 then return false;fi;

      # # Naive Euler Char of the boundary S of L

      Boundaries[1]:=List(Cells[1],a->[1,0]); # denote the boundary of 0-cell
      Boundaries[k+2]:=[];
      ######BOUNDARIES END########
      ######COBOUNDARIES BEGIN####
      Coboundaries:=[];; #Coboundaries[n+1] contains the info on n-cells.
      for j in [0..k] do

        Coboundaries[j+1]:=List(Boundaries[j+1],i->[0]);
        for s in [1..Length(Boundaries[j+2])] do
          b:=Boundaries[j+2][s];
          t:=Length(b);
          for i in b{[2..t]} do
            Coboundaries[j+1][i][1]:=Coboundaries[j+1][i][1]+1;
            Add(Coboundaries[j+1][i],s);
          od;
        od;

      od;
      Coboundaries[k+1]:=List(Boundaries[k+1],a->[0]);
      #####COBOUNDARIES END######
      Y:=ListsOfCellsToRegularCWComplex(Boundaries,Coboundaries,k);
      crit:=CriticalCellsOfRegularCWComplex(Y);

      if (Length(crit)=1 and crit[1][1]=0) then return true;
      else return false;
      fi;

    end;
    ###############
    CellHomology:=function(k,L)
    local Cells, nrCells,x,w,j,s,id, Mat,t,M,b,p,h,r;

    Cells:=[];

    for j in [1..k+1] do
        Cells[j]:=[];

    od;

    j:=k;
    Append(Cells[k+1],L);


    while j>0 do
        for s in [1..Length(Cells[j+1])] do
            x:=Cells[j+1][s];
            w:=StructuralCopy(BoundaryRec[j][AbsInt(x[1])]);
            w:=Mult(w,j-1,x[2]);
            w:=List(w,a->[AbsInt(a[1]),a[2]]);
            Cells[j]:=Union(Cells[j],w);
        od;
        j:=j-1;
    od;
    nrCells:=List([1..k+1],a->Length(Cells[a]));

  Mat:=[];
  for j in [2..k+1] do
    M:=[];
    for s in [1..Length(Cells[j])] do
       M[s]:=[];
       for t in [1..Length(Cells[j-1])] do
          M[s][t]:=0;
       od;
    od;

    for s in [1..Length(Cells[j])] do
        x:=Cells[j][s];
        b:=StructuralCopy(BoundaryRec[j-1][AbsInt(x[1])]);
        b:=Mult(b,j-2,x[2]);
        for t in [1..Length(b)] do
            p:=Position(Cells[j-1],[AbsInt(b[t][1]),b[t][2]]);
            M[s][p]:=SignInt(b[t][1]);
        od;
    od;
    Mat[j]:=SmithNormalFormIntegerMat(TransposedMat(M));


  od;
    w:=[];
    # Create Mat[1] for d_0
    M:=[];
    for s in [1..Length(Cells[1])] do
        M[s]:=0;
    od;
    Mat[1]:=TransposedMat([M]);

    # Create Mat[k+1] do d_(k+1)
    M:=[];
    for s in [1..Length(Cells[k+1])] do
        M[s]:=0;
    od;
    Mat[k+2]:=[M];

    for j in [1..k+1] do
        h:=[];
        r:=Length(Cells[j])-RankMat(Mat[j])-RankMat(Mat[j+1]);
        for s in [1..r] do
	    h[s]:=0;
        od;
        for s in [1..RankMat(Mat[j+1])] do
            if not Mat[j+1][s][s]=1 then
                Add(h,Mat[j+1][s][s]);
            fi;
        od;
        Add(w,h);

    od;



     return w;


    end;



    ##################################################################

    ##################################################################
    # Main part: subdividing the fundamental domain
    NotRigid:=[];
    dims:=ShallowCopy(DimRec);
    i:=1;
    while i<=M do


        j:=1;
        while j<=dims[i+1] do
            if not IsRigidCell(i,j) then
                Add(NotRigid,[i,j]);
            fi;
            j:=j+1;
        od;

        i:=i+1;

    od;
    for x in NotRigid do
#       Apply Rigid Facets Subdivision to the cell x :
        ReplaceCell(x[1],x[2]);
    od;


    #Delete cells which are already replaced by its subdivision

    t:=1;
    for w in [1..Length(NotRigid)] do
        k:=NotRigid[w][1];
        j:=NotRigid[w][2];
        if k<N then
        for i in [1..DimRec[k+2]] do
            bdry:=BoundaryRec[k+1][i];

            if not IsString(bdry) then
               for x in bdry do
                    if AbsInt(x[1])>j then
                        x[1]:=x[1]-SignInt(x[1]);
                    fi;
                od;
             fi;
             BoundaryRec[k+1][i]:=bdry;

         od;
         fi;
         dims[k+1]:=dims[k+1]-1;
         DimRec[k+1]:=DimRec[k+1]-1;
         Remove(BoundaryRec[k],j);
         Remove(StabRec[k+1],j);
         if IsBound(NotRigid[w+1]) and NotRigid[w+1][1]=NotRigid[w][1] then
             NotRigid[w+1][2]:=NotRigid[w+1][2]-t;
             t:=t+1;
         else
             t:=1;
         fi;

    od;

    ##################################################################
    Boundary:=function(k,m)
        if m<0 then return NegateWord(BoundaryRec[k][AbsInt(m)]);fi;
        return BoundaryRec[k][m];
    end;

    Stabilizer:=function(k,m)
        return StabRec[k+1][m];
    end;

    Dimension:=function(k)
        if k>N then return 0;fi;
        return DimRec[k+1];
    end;

    Action:=function(k,i,j)
        return 1;
    end;
    ##################################################################
return Objectify(HapNonFreeResolution,
    rec(
    dimension:=Dimension,
    Partition:=Partition,
    boundary:=Boundary,
    homotopy:=fail,
    elts:=Elts,
    group:=C!.group,
    stabilizer:=Stabilizer,
    action:=Action,
    subdividing:=SubdividingCell,
    replacecell:=ReplaceCell,
    isrigid:=IsRigidCell,
    properties:=
    [["length",Maximum(1000,N)],
    ["characteristic",0],
    ["type","resolution"]]  ));
end);
