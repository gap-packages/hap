
##########################################################################
#0
#F  SL2ZTree
##
##  This function will compute a G-equivalent CW-space for 
##  group G=SL2Z[1/(m*p)] and a tree X associated with 
##  amalgamated product SL2Z[1/m]*SL2z[1/m] over
##  congruence subgroup gamma 0 of level p.
##  Input:  A pair of positive integers (m,p)
##         
##  Output: A G-equivarient CW-space where G is the free  
##          product of two copiesof SL(1/m) amalgamated by 
##          its congruence group of level p  
##
InstallGlobalFunction(SL2ZTree,
function(m,p)
local t1,t2,
    Elts,H,K,G,Gamma,
    Id,ID,Idcoset,BoundaryList,
    Boundary,Dimension,Action,Stabilizer,Homotopy,StabGrps,
    pos,RemoveLoops,HtpyRec;
    
    Elts:=[[[1,0],[0,1]]];
    
    ## In this case G=SL2Z
    if p=0 then
        H:=Group([[0,-1],[1,0]]);
        K:=Group([[0,-1],[1,1]]);
        G:=SL(2,Integers);
        Gamma:=Group([[-1,0],[0,-1]]);
        Append(Elts,Elements(H));
        Append(Elts,Elements(K));
        Append(Elts,Elements(Gamma));
        Elts:=SSortedList(Elts);

    else 
        if m=1 then 
            H:=SL(2,Integers);
            K:=SL2Z(p);
            Gamma:=CongruenceSubgroupGamma0(p);
        else 
            H:=SL2Z(1/m);
            K:=ConjugateSL2ZGroup(H,[[1,0],[0,p]]);
            Gamma:=CongruenceSubgroup(m,p);
        fi;
        G:=SL2Z(1/(m*p));
        ID:=Group(One(G));
        Append(Elts,GeneratorsOfGroup(H));
        Append(Elts,GeneratorsOfGroup(K));
        Append(Elts,GeneratorsOfGroup(Gamma));
        Elts:=SSortedList(Elts);
    fi;
    SetName(Gamma,"Gamma");
    Id:=Position(Elts,[[1,0],[0,1]]);
    
    ###################################################################     
    #1
    #F  pos
    ##  Find the position of an element in a list 
    ##  Input:  A list Elts and a matrix g
    ##  Output: If g in Elts return the position of g in the list,
    ##          otherwise, add g to Elts and return the position.
    ##
    pos:=function(Elts,g)
    local posit;

        posit:=Position(Elts,g);
        if posit=fail then 
            Add(Elts,g);  
            return Length(Elts);
        else  
            return posit;
        fi;
    end;
    ###################################################################
    
    BoundaryList:=[];
    t1:=pos(Elts,CanonicalRightCountableCosetElement(H,Elts[Id]^-1)^-1); 
    t2:=pos(Elts,CanonicalRightCountableCosetElement(K,Elts[Id]^-1)^-1);
    Append(BoundaryList,[[[1,t1],[-2,t2]]]);
    
    ###################################################################     
    #1
    #F  Boundary
    ##  
    ##  This function presents the boundary map d_n: C_n -> C_{n-1}
    ##
    ##  Input:  a pair of integers (n,k) where n is the dimension and 
    ##          k is the position of the generator.
    ##  Output: a list of words [g,f] 
    ##
    Boundary:=function(n,k)
    local w;
        if not n=1 then 
            return [];
        fi;
        w:=BoundaryList[AbsInt(k)];
        if k>0 then 
            return w;
        else 
            return NegateWord(w);
        fi;
    end;
    
    ###################################################################     
    #1
    #F  Dimension
    ##  
    ##  This function computes the G-rank of ZG-module C_n
    ##
    ##  Input:  a positive integer n 
    ##          
    ##  Output: the G-rank of ZG-module C_n 
    ##
    Dimension:=function(n)
        if not n in [0,1] then return 0;fi;
        if n=0 then return 2;fi;
        if n=1 then return 1;fi;
    end;
    ###################################################################

    ###################################################################     
    #1
    #F  Action
    ##
    ##  Input:  a triple (n,k,l) of integers
    ##          
    ##  Output: 1 or -1  
    ##
    Action:=function(n,k,l);
        return 1;
    end;
    ###################################################################
    
    StabGrps:=[];
    Add(StabGrps,[H,K]);
    Add(StabGrps,[Gamma]);
    
    ###################################################################     
    #1
    #F  Stabilizer
    ##
    ##  Input:  a pair of integers (n,k)
    ##          
    ##  Output: the kith stabiliser subgroup in dimension n   
    ##
    Stabilizer:=function(n,k);
        return StabGrps[n+1][k];
    end;
    ###################################################################

    Idcoset:=pos(Elts,
            CanonicalRightCountableCosetElement(Gamma,Elts[Id]^-1)^-1);
    
    ###################################################################     
    #1
    #F  RemoveLoops
    ##
    ##  This function will remove loops in a path.
    ##
    ##  Input:  a list of integers d 
    ##          
    ##  Output: a list of integers with no loop.   
    ##
    RemoveLoops:=function(d)
    local i,h,j,l;
        l:=StructuralCopy(d);
        h:=[[1,0],[0,1]];
        i:=1;
        while i<Length(d) do
            h:=h*d[i];
            if h in H or h in K then
                for j in [1..i-1] do
                    Remove(l,1);
                od;
                l[1]:=h;
            fi;
            i:=i+1;
        od;
        return l;
    end;
    ###################################################################

    ## Create a record for the homotopy 
    HtpyRec:=[];
    HtpyRec[1]:=[];
    HtpyRec[2]:=[];
    
    ###################################################################     
    #1
    #F  Homotopy
    ##
    ##  This function presents the homotopy map h_n: C_n -> C_{n+1}
    ##
    ##  Input:  a positive integer n and a word w. 
    ##          
    ##  Output: a list of words.   
    ##
    Homotopy:=function(n,w)
    local d,path,i,h,k,g,pk,r,t;

    if not n=0 then 
        return [];
    fi;
    k:=w[1];
    g:=w[2];
    pk:=AbsInt(k);
    if not IsBound(HtpyRec[pk][g]) then
        d:=SL2ZmElementsDecomposition(Elts[g],p);
        r:=[];
        Add(r,d[1]);
        for i in [2..Length(d)] do
            if d[i]*d[i-1] in H then
                r[Length(r)]:=r[Length(r)]*d[i];
            else
                Add(r,d[i]);
            fi;
        od;
        d:=r;
        r:=[];
        Add(r,d[1]);
        for i in [2..Length(d)] do
            if d[i]*d[i-1] in K then
                r[Length(r)]:=r[Length(r)]*d[i];
            else
                Add(r,d[i]);
            fi;
        od;
        d:=StructuralCopy(r);

        ## kill all the loops in the path d

        d:=RemoveLoops(d);
        if (d[1] in K) and (not d[1] in Gamma) then
            r:=[[[1,0],[0,1]]];
            Append(r,d);
            d:=StructuralCopy(r);
        fi;
        h:=[[1,0],[0,1]];
        path:=[];
        if d[Length(d)] in Gamma and Length(d)>1 then 
            Remove(d,Length(d));
        fi;

        if pk=1 then 
            if d[Length(d)] in H then 
                Remove(d,Length(d));fi;
            for i in [1..Length(d)] do
                h:=h*d[i];
                t:=CanonicalRightCountableCosetElement(Gamma,h^-1)^-1;
                Add(path,[(-1)^(i),pos(Elts,t)]);
            od;
        else
            if Elts[g] in Gamma then 
                Add(path,[-1,Idcoset]);
            fi;
            if d[Length(d)] in K then 
                Remove(d,Length(d));
            fi;
            for i in [1..Length(d)] do
                h:=h*d[i];
                t:=CanonicalRightCountableCosetElement(Gamma,h^-1)^-1;
                Add(path,[(-1)^(i),pos(Elts,t)]);
            od;
        fi;
        HtpyRec[pk][g]:=path;
    fi;
    if k>0 then
        return HtpyRec[pk][g];
    else 
        return NegateWord(HtpyRec[pk][g]);
    fi;
    end;
 
   #######################################################################

return Objectify(HapNonFreeResolution,
    rec(
    dimension:=Dimension,
    boundary:=Boundary,
    homotopy:=Homotopy,
    elts:=Elts,
    group:=G,
    stabilizer:=Stabilizer,
    action:=Action,
    properties:=
    [["length",100],
    ["characteristic",0],
    ["type","resolution"]]  ));
end);
################### end of SL2ZTree ######################################

##########################################################################
#0
#F  TreeOfResolutionsToSL2Zcomplex
##  Input: A list of resolutions D and an arithmetic group G=SL2Z[1/m] 
##         
##  Output: A G-equivarient CW-space 
##             
##
InstallGlobalFunction(TreeOfResolutionsToSL2Zcomplex,
function(D,G)
local RH,RK,RGamma,
    H,K,Gamma,
    NameH,
    i,j,m,p,
    C,Resolutions,NamesOfGroups;
    
    RH:=D[1];
    RK:=D[2];
    RGamma:=D[3];
    H:=RH!.group;
    K:=RK!.group;
    Gamma:=RGamma!.group;
    NameH:=H!.Name;
    if H=SL(2,Integers) then 
        m:=1;
        p:=Gamma!.LevelOfCongruenceSubgroup;
    else
        i:=Position(NameH,'/');
        j:=Position(NameH,']');
        m:=Int(NameH{[i+1..j-1]});
        p:=Gamma!.levels[2];
    fi;
    NamesOfGroups:=[Name(H),Name(K),Name(Gamma)];
    Resolutions:=[RH,RK,RGamma];
    C:=SL2ZTree(m,p);
    C!.resolutions:=[Resolutions,NamesOfGroups];
    return C;
end);
################### end of TreeOfResolutionsToSL2Zcomplex ################

