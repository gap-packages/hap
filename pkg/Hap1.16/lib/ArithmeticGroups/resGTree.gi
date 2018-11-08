
##########################################################################
#0
#F  ResolutionGTree
##  
##  This function use C.T.C Wall’s method to build a free ZG-resolution
##  from a G-equivariant CW-complex given resolutions for all stabilisers
##
##  Input: A non free ZG-resolution with resolutions for all stabilizers, 
##         a positive integer n 
##  Output: The first n+1 terms of a free ZG-reoslution of Z. 
##             
##

InstallGlobalFunction(ResolutionGTree,
function(arg)
local 
	R,n,G,i,L,
	StabRes,StabGrps,Triple2Pair,Quad2One,GrpsRes,Pair2Quad,Quad2Pair,
	Dimension,Hmap,Gmult,Gmultrec,Mult,AlgRed,CorrectList,Boundary,
	Homotopy,FinalHomotopy,
	StRes,hmap,Action,PseudoBoundary,PseudoHomotopy, ZeroDimensionHmap, 
        ZeroDimensionHtpy,
	HmapRec,p,q,r,s,HtpyRec,k,g, ZeroDimensionHmapRec, Pair2QuadRec,
        QuadToPairRec,DimensionRec;
    
    R:=arg[1];
    n:=arg[2];
    G:=R!.group;

    ###################################################################     
    #1
    #F  Action
    ##
    ##  Input:  a triple of integers (p,r,q)
    ##  Output: 1 or -1 
    ##

    Action:=function(p,r,g)
    if not IsBound(R!.action) then 
        return 1;
    else 
        return R!.action(p,r,g);
    fi;
    end;

    ###################################################################
    AlgRed:=AlgebraicReduction;

    Gmultrec:=[];

    ###################################################################     
    #1
    #F  Gmult
    ##
    ##  The product of 2 elements in Elts 
    ##
    ##  Input:  a pair of positive integers (i,j)
    ##  Output: the position of the product in the list Elts  
    ##
    Gmult:=function(i,j)
    local posit,g;

    if not IsBound(Gmultrec[i]) then
        Gmultrec[i]:=[];
    fi;
    if not IsBound(Gmultrec[i][j]) then
        g:=R!.elts[i]*R!.elts[j];
        posit:=Position(R!.elts,g);
        if posit=fail then
            Add(R!.elts,g);
            posit:= Length(R!.elts);
        fi;
        Gmultrec[i][j]:=posit;
    fi;

    return Gmultrec[i][j];
    end;
    ################################################################### 

    ###################################################################     
    #1
    #F  Mult(g,w)
    ##
    ##  Multiply gth-element with a list of words w
    ##
    ##  Input:  an integer g and a list of words w 
    ##  Output: the product of q and w   
    ##    

    Mult:=function(g,w) # Multiply gth-element with a word
    local 
        l,x;

    l:=StructuralCopy(w);
    if R!.elts[g]=[] then 
        return [];
    fi;
    Apply(l,y->[y[1],Gmult(g,y[2])]);
    return l;
    end;
    ################################################################### 

    ###################################################################     
    #1
    #F  GrpsRes
    ##
    ##  Return a free resolution for a given group 
    ##
    ##  Input:  A group G and a positive integer n 
    ##  Output: The first n+1 term of a free ZG-resolution of Z 
    ##  
    GrpsRes:=function(G,n) # Resolutions of Group
    local 
        iso,Q,res,x;
    if IsBound(R!.resolutions) and HasName(G) then 
        x:=Position(R!.resolutions[2], Name(G)); 
        if not x=fail then 
            return R!.resolutions[1][x]; 
        fi;
    fi;
    iso:=RegularActionHomomorphism(G);
    Q:=Image(iso);
    res:=ResolutionFiniteGroup(Q,n);
    res!.group:=G;
    res!.elts:=List(res!.elts,x->PreImagesRepresentative(iso,x));
    return res;
    end;
    ###################################################################

    #Create list of stabilizer subgroups and their resolutions
    StabGrps:= List([0..Length(R)],n->
                   List([1..R!.dimension(n)], k->R!.stabilizer(n,k))); 
    StabRes:=[];
    for L in StabGrps do
        Add(StabRes,List(L,g->ExtendScalars(GrpsRes(g,n),G,R!.elts))  ); 
    od;
    ###################################################################
    CorrectList:=function(list)
    local 
        l,i;
    if list=[] then return [];fi;
    l:=StructuralCopy(list[1]);
    for i in [2..Length(list)] do
        Append(l,StructuralCopy(list[i]));
    od;
    return l;
    end;
    ###################################################################

 
    ###################################################################     
    #1
    #F  Quad2One
    ##
    ##  return nth-generator of F_(p,q) from (r,s)th-generator of 
    ##  stabilizer 
    ##
    ##  Input:  A quadruple of integers (p,q,r,s) 
    ##  Output: n-th generator of F_{p,q} 
    ## 
    Quad2One:=function(p,q,r,s)
    local
        n,d,i,j;
    n:=0;
    i:=SignInt(s);
    s:=AbsInt(s);
    d:=List([1..R!.dimension(p)],x->StabRes[p+1][x]!.dimension(q));
    for j in [1..r-1] do
        n:=n+d[j];
    od;
    n:=n+s;
    if q=0 and n>R!.dimension(p) then 
        n:=R!.dimension(p);
    fi;
    return i*n;
    end;
    ###################################################################

    ###################################################################     
    #1
    #F  Triple2Pair
    ##
    ##  Input:  A triple of integers (p,q,n) 
    ##  Output: a pair of integers (r,s) 
    ## 
    Triple2Pair:=function(p,q,n)
    local 
    r,s,d,i;
    r:=0;
    d:=List([1..R!.dimension(p)],x->StabRes[p+1][x]!.dimension(q));
    i:=SignInt(n);
    n:=AbsInt(n);
    while n>0 do
        r:=r+1;
        s:=n;
        n:=n-d[r];
    od;
    return [r,i*s];
    end;
    ###################################################################

    # Create a record for horizontal map: d_1 
    HmapRec:=[];
    for p in [1..2] do
        HmapRec[p]:=[];
        for q in [1..n+1] do
            HmapRec[p][q]:=[];
            for r in [1..R!.dimension(p-1)] do
                HmapRec[p][q][r]:=[];
            od;
        od;
    od;
    ZeroDimensionHmapRec:=[];


    ###################################################################     
    #1
    #F  ZeroDimensionHmap
    ##
    ##  Input:  An integer k 
    ##  Output: The map d_1 at degree 0 
    ## 
    ZeroDimensionHmap:=function(k)
    local i,j,pk;
    pk:=AbsInt(k);
    if not IsBound(ZeroDimensionHmapRec[pk]) then
        j:=0;
        for i in [1..pk-1] do
            j:=j+StabRes[1][i]!.dimension(0);
        od;
        j:=j+1;
        ZeroDimensionHmapRec[pk]:=j;
    fi;
    if k>0 then 
        return ZeroDimensionHmapRec[pk];
    else return -ZeroDimensionHmapRec[pk];fi;
    end;
    ################################################################### 

    ###################################################################     
    #1
    #F  ZeroDimensionHmap
    ##
    ##  Horiziontal map d_1:A(p,q)->A(p-1,q), acts on the
    ##  (r,s) th-generator of A(p,q)
    ##
    ##  Input:  An integer k 
    ##  Output: The map d_1 at degree 0 
    ## 
    Hmap:=function(p,q,r,s)     
    local 
        i,l,d0,m,bdr,ps,d1d0,w;
    ps:=AbsInt(s);
    if p<>1 then 
        return [];
    else
        if not IsBound(HmapRec[p+1][q+1][r][ps]) then
            if q=0 then bdr:=StructuralCopy(R!.boundary(1,1));
                Apply(bdr,w->[ZeroDimensionHmap(w[1]),w[2]]);
                HmapRec[p+1][q+1][r][ps]:=bdr;
            else
                l:=[];m:=[];
                d0:=StructuralCopy(List(StabRes[p+1][r]!.boundary(q,ps),
                                      x->[Action(p,r,x[2])*x[1],x[2]]));
                for w in d0 do
                    Append(m,Mult(w[2],Hmap(p,q-1,r,w[1])));
                od;
                Apply(m,x->[Triple2Pair(p-1,q-1,x[1]),x[2]]);
                for w in m do
                    Append(l,List(StabRes[p][w[1][1]]!.homotopy
                                      (q-1,[w[1][2],w[2]]),y->
                                  [Quad2One(p-1,q,w[1][1],y[1]),y[2]]));
                od;
                HmapRec[p+1][q+1][r][ps]:=AlgRed(l);
            fi;
        fi;
    fi;
    if SignInt(s)=1 then 
        return HmapRec[p+1][q+1][r][ps];
    else 
        return NegateWord(HmapRec[p+1][q+1][r][ps]);
    fi;
    end;
    ###################################################################

    Pair2QuadRec:=[];

    ###################################################################     
    #1
    #F  Pair2Quad
    ##
    ##  Input:  A pair of integers (k,n) 
    ##  Output: A triple of integers  
    ## 
    Pair2Quad:=function(k,nn)
    local 
        x,n,nnn, p,q,r,s,i,temp,j1,j2;
    i:=SignInt(nn);
    n:=AbsInt(nn);
    nnn:=n;

    if not IsBound(Pair2QuadRec[k+1]) then 
        Pair2QuadRec[k+1]:=[]; 
    fi;

    if not IsBound(Pair2QuadRec[k+1][n]) then 

        temp:=0;
        for j1 in [0..k] do
            for j2 in [1..R!.dimension(j1)] do
                temp:=temp+StabRes[j1+1][j2]!.dimension(k-j1);
            od;
        od;

        p:=-1;
        while n>0 do;
            p:=p+1;
            r:=0;
            while (n>0 and r<R!.dimension(p)) do
                r:=r+1;
                s:=n;
                n:=n-StabRes[p+1][r]!.dimension(k-p);
            od;
        od;
        q:=k-p;
        Pair2QuadRec[k+1][nnn]:=[p,q,r,s];

    fi;

    x:=Pair2QuadRec[k+1][nnn];
    return [x[1],x[2],x[3],i*x[4]];

    end;
    ###################################################################

    QuadToPairRec:=[];

    ###################################################################     
    #1
    #F  Pair2Quad
    ##
    ##  Input:  A pair of integers (k,n) 
    ##  Output: A triple of integers  
    ## 
    Quad2Pair:=function(p,q,r,s)
    local
        k,n,i,j,p1,q1,r1,s1;

    p1:=p+1;q1:=q+1;r1:=r+1;s1:=AbsInt(s)+1;
    if not IsBound(QuadToPairRec[p1]) then 
        QuadToPairRec[p1]:=[]; 
    fi;
    if not IsBound(QuadToPairRec[p1][q1]) then 
        QuadToPairRec[p1][q1]:=[]; 
    fi;
    if not IsBound(QuadToPairRec[p1][q1][r1]) then 
        QuadToPairRec[p1][q1][r1]:=[]; 
    fi;
    if not IsBound(QuadToPairRec[p1][q1][r1][s1]) then 
        k:=p+q;
        n:=0;
        for i in [0..p-1] do
            for j in [1..R!.dimension(i)] do
                n:=n+StabRes[i+1][j]!.dimension(k-i);
            od;
        od;	
        for i in [1..r-1] do
            n:=n+StabRes[p+1][i]!.dimension(k-p);
        od;
        n:=n+AbsInt(s);
        QuadToPairRec[p+1][q+1][r+1][AbsInt(s)+1]:=[k,n];

    fi;
    k:=QuadToPairRec[p1][q1][r1][s1];
    return [k[1],SignInt(s)*k[2]];
    end;
    ###################################################################

    # Create an empty list for the pseudo boundary 
    PseudoBoundary:=[];
    for k in [1..n+1] do
        PseudoBoundary[k]:=[];
    od;

    ###################################################################     
    #1
    #F  Boundary
    ##
    ##  Input:  A pair of integers (k,n) 
    ##  Output: the boundary d(k,n)  
    ## 
    Boundary:=function(k,n)
    local
        d,l,p,q,r,s,w,pn;
    pn:=AbsInt(n);
    if not IsBound(PseudoBoundary[k+1][pn]) then
        w:=Pair2Quad(k,pn);
        p:=w[1];q:=w[2];r:=w[3];s:=w[4];
        d:=[];
        if q<>0 then 
	    l:=StructuralCopy(List(StabRes[p+1][r]!.boundary(q,s),
              x->[Quad2Pair(p,q-1,r,Action(p,r,x[2])*x[1])[2],x[2]]));
            Append(d,StructuralCopy(l));
        fi;
        if IsEvenInt(q) then 
            Append(d,StructuralCopy(Hmap(p,q,r,s)));
        else
            Append(d,StructuralCopy(NegateWord(Hmap(p,q,r,s))));
        fi;
        PseudoBoundary[k+1][pn]:=AlgRed(d);	
    fi;
    if SignInt(n)=1 then 
        return PseudoBoundary[k+1][pn];
    else 
        return NegateWord(PseudoBoundary[k+1][pn]);
    fi;
    end;
    ###################################################################

    DimensionRec:=[]; 

    ###################################################################     
    #1
    #F  Boundary
    ##
    ##  Input:  A pair of integers (k,n) 
    ##  Output: A triple of integers  
    ## 
    Dimension:=function(n)
    local
        dim,p,i;
 
    if not IsBound(DimensionRec[n+1]) then
        dim:=0;
        for p in [0..n] do
            for i in [1..R!.dimension(p)] do
                dim:=dim+StabRes[p+1][i]!.dimension(n-p);
            od;
        od;
        DimensionRec[n+1]:= dim;
    fi;

    return DimensionRec[n+1];
    end;
    ###################################################################

    # Create a record for the homotopy map 
    HtpyRec:=[];
    for k in [1..n] do
        HtpyRec[k]:=[];
        for s in [1..Dimension(k-1)] do
            HtpyRec[k][s]:=[];
        od;
    od;

    ###################################################################     
    #1
    #F  ZeroDimensionHtpy
    ##
    ##  Input:  An integer k 
    ##  Output: the homotopy h(0,[1,k]) of the chain complex  
    ## 
    ZeroDimensionHtpy:=function(k)
    local i,j,r;
    i:=0;
    while k>0 do
        i:=i+1;
        k:=k-StabRes[1][i]!.dimension(0);
        r:=i;
    od;
    return r;
    end;
    ###################################################################

    ###################################################################     
    #1
    #F  Homotopy
    ##
    ##  Input:  An integer n and a word w=[f,g] 
    ##  Output: the homotopy h(n,w)   
    ##
    Homotopy:=function(n,w)
    local 
        t,g,h0,h11,e,h,dh,
        p,q,r,s,v,m,pt,ppt,
        h1,d1h1,x,k,y,ps;
    t:=w[1];
    g:=w[2];
    e:=[];
    h:=[];
    dh:=[];
    pt:=AbsInt(t);
    v:=Pair2Quad(n,pt);#Print(v);
    p:=v[1];q:=v[2];r:=v[3];s:=v[4];
    if not IsBound(HtpyRec[n+1][pt][g]) then
        if n=0 then
            ppt:=ZeroDimensionHtpy(pt);
            h1:=StructuralCopy(R!.homotopy(n,[ppt,g]));
            d1h1:=StructuralCopy(AlgRed(CorrectList(List(h1,x->
                                   Mult(x[2],Hmap(p+1,q,1,x[1]))))));
	        for x in d1h1 do
	            k:=Pair2Quad(n,x[1]);
	            y:=StructuralCopy(StabRes[k[1]+1][k[3]]!.
                                          homotopy(q,[k[4],x[2]]));
	            Apply(y,w->[Quad2Pair(k[1],k[2]+1,k[3],w[1])[2],w[2]]);
	            Append(e,y);
	        od;
	        h0:=StructuralCopy(StabRes[p+1][r]!.homotopy(0,[s,g]));
	        Apply(h0,w->[Quad2Pair(p,q+1,r,w[1])[2],w[2]]);
	        h11:=List(h1,x->[Quad2Pair(p+1,q,Triple2Pair(p+1,q,x[1])[1],
                               Triple2Pair(p+1,q,x[1])[2])[2],x[2]]);
	        Append(h,NegateWord(e));
	        Append(h,h0);
	        Append(h,h11);
	        HtpyRec[n+1][pt][g]:=AlgRed(h);
        else
	        if p=0 then 
	        h0:=StructuralCopy(StabRes[p+1][r]!.homotopy(q,[s,g]));
	        Apply(h0,w->[Quad2Pair(p,q+1,r,w[1])[2],w[2]]);
	        Append(h,h0);
	        else
	            ps:=Action(1,1,g)*s;
	            m:=StructuralCopy(StabRes[p+1][r]!.homotopy(q,[ps,g]));
	            Apply(m,x->[Action(1,1,x[2])*x[1],x[2]]);
	            h0:=List(m,x->[Quad2Pair(p,q+1,r,x[1])[2],x[2]]);
	            Append(h,h0);
	            dh:=AlgRed(CorrectList(List(m,x->
                                 Mult(x[2],Hmap(p,q+1,1,x[1])))));
	            for x in dh do
	                k:=Pair2Quad(n,x[1]);
	                y:=StructuralCopy(StabRes[k[1]+1][k[3]]!.
                                   homotopy(k[2],[k[4],x[2]]));
	                Apply(y,w->[Quad2Pair(k[1],k[2]+1,k[3],w[1])[2],w[2]]);
	                Append(e,y);
	            od;
	            if IsEvenInt(q) then 
	                Append(h,e);
	            else 
				    Append(h,NegateWord(e));
                fi;
	        fi;
            HtpyRec[n+1][pt][g]:=AlgRed(h);
        fi;
    fi;
    if SignInt(t)=1 then 
        return HtpyRec[n+1][pt][g];
    else 
        return NegateWord(HtpyRec[n+1][pt][g]);
    fi;
    end;

    ###################################################################

    ###################################################################     
    #1
    #F  FinalHomotopy
    ##
    ##  Input:  An integer n and a word w=[f,g] 
    ##  Output: the homotopy h(n,w)   
    ##
    FinalHomotopy:=function(n,g)
    if R!.homotopy=fail then
        return fail;
    else 
        return Homotopy(n,g);
    fi;
    end;


    ####ADDED MAY 2012################################################
    StRes:=function(n,k)
        return StabRes[n+1][k];
    end;
    
######################################################################
return Objectify(HapResolution,
                rec(
                dimension:=Dimension,
                boundary:=Boundary,
                homotopy:=FinalHomotopy,
                elts:=R!.elts,
                group:=R!.group,
                stabres:=StRes,
                properties:=
                   [["length",n],
                    ["initial_inclusion",true],
                    ["type","resolution"],
                    ["characteristic",EvaluateProperty(
                                  R,"characteristic")]  ]	));
end);

##
################### end of ResolutionGTree ###########################

