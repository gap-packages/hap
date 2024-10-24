

#Note:    The first part of this function written by Graham Ellis
#         The "homotpy" part is written by Bui Anh Tuan 

##########################################################################
#0
#F  sl2zngens.gi
##  Input: A G-equivalent CW-space with resolutions for all stabilisers 
##         a positive integer n. 
##  Output: The first n terms of a free ZG-resolution of Z. If the input 
##          together with a contracting homotopy then the output together 
##          with a contracting homotopy.
##
InstallGlobalFunction(FreeZGResolution,
function(arg)
local 
    P,N,prime,
    Dimension, DimensionRecord, DimRecs, FiltDimRecs, BinGp,
    Boundary, BoundaryP, Pair2Quad, Pair2QuadRec,
    Quad2Pair,Quad2PairRec, HtpyGen, HtpyWord,
    StabGrps, StabResls, ResolutionFG,
    Action, AlgRed,  EltsG, G, Mult, MultRecord,
    DelGen, DelWord, DelGenRec,PseudoBoundary,FinalBoundary,
    FilteredLength, FilteredDimension, FilteredDimensionRecord,
    L,i,k,n,q,r,s,t,
    InducedHtpyGen, InducedHtpyWord, DelListSum, #Added Feb 2014 BUI A.T.
    Homotopy, HomotopyGen, NegateListWord, VertHtpy, 
    InducedHtpyList, IndHtpyRec;    

    SetInfoLevel(InfoWarning,0);

    P:=arg[1];
    N:=arg[2];
    if Length(arg)>2 then 
        prime:=Gcd(arg[3],EvaluateProperty(P,"characteristic"));
    else 
        prime:=EvaluateProperty(P,"characteristic"); 
    fi;

    N:=Minimum(EvaluateProperty(P,"length"),N);
    G:=P!.group;
    EltsG:=P!.elts;
    BoundaryP:=P!.boundary;

    BinGp:=ContractibleGcomplex("SL(2,O-2)");
    BinGp:=BinGp!.stabilizer(0,4);;
    BinGp:=Image(RegularActionHomomorphism(BinGp));

    ###################################################################     
    #1
    #F  ResolutionFG
    ##
    ##  Return resolutions for stabiliser groups if they are given
    ##  otherwise, find a new one.
    ##
    ##  Input:  group G and a positive integer n
    ##  Output: the first n+1 terms of a free ZG-resolution of Z  
    ##
    ResolutionFG:=function(G,n)
    local 
        x, tmp, iso,iso1,iso2,iso3,res,Q, fn;

        ##Added Jan 2012
        if IsBound(P!.resolutions) and HasName(G) then
            x:=Position(P!.resolutions[2], Name(G));
            if not x=fail then 
                return P!.resolutions[1][x]; 
            fi;
        fi;

        if Order(G)=infinity and IsAbelian(G) then
            #This will only be correct if G is abelian of "rank" equal
            #to the number of generators GAP has for G 

            res:=ResolutionGenericGroup(G,n);
            return res;
        fi;
        iso:=RegularActionHomomorphism(G);
        Q:=Image(iso);

        if IdGroup(Image(iso))=[24,3] then 
            iso1:=IsomorphismGroups(Q,BinGp); 
            res:=ResolutionFiniteGroup(BinGp,n);
            res!.group:=G;
            res!.elts:=List(res!.elts,x->
            PreImagesRepresentative(iso,PreImagesRepresentative(iso1,x)));
            return res;
        fi;

        res:=ResolutionFiniteGroup(Q,n); 
        res!.group:=G;
        res!.elts:=List(res!.elts,x->PreImagesRepresentative(iso,x));
        return res;
    end;    
    ############### end of ResolutionFG ###############################


    if prime>0 then 

        ###############################################################     
        #1
        #F  AlgRed
        ##
        ##  Algebraic reduction for list of words mod p
        ##
        ##  Input:  a list of words
        ##  Output: reduced list of words  
        ##
        AlgRed:= function(ww)
        local 
            w,x,v,pos,u;

            w:=StructuralCopy(ww);

            v:=Collected(w);
            for x in v do
                if x[1][1]<0 then 
                    x[1][1]:=-x[1][1]; x[2]:=-x[2] mod prime; 
                fi;
                if x[1][2]<0 then 
                    x[1][2]:=-x[1][2]; x[2]:=-x[2] mod prime; 
                fi;
                x[2]:=x[2] mod prime;
            od;

            u:=[];
            for x in v do
                Append(u,MultiplyWord(x[2],[x[1]]));
            od;

            v:=Collected(u);
            for x in v do
                x[2]:=x[2] mod prime;
            od;

            u:=[];
            for x in v do
                Append(u,MultiplyWord(x[2],[x[1]]));
            od;

            return u;
        end;
        ######################## end of AlgRed ########################

    else


        ###############################################################    
        #1
        #F  AlgRed
        ##
        ##  Algebraic reduction for list of words 
        ##
        ##  Input:  a list of words
        ##  Output: reduced list of words  
        ##
        AlgRed:= function(ww)
        local x,i,v,k,u,w;
            w:=ww;#w:=StructuralCopy(ww);

            for x in w do
                if x[2]<0 then 
                    x[1]:=-x[1];x[2]:=-x[2];
                fi;
            od;
            v:=Filtered(w,x->x[1]>0);
            for x in w do
                if x[1]<0 then
                    k:=Position(v,[-x[1],x[2],x[3]]);  
                    if (k=fail) then 
                        Add(v,x); 
                    else
                        Unbind(v[k]);
                    fi;
                fi;
            od;
            v:=Filtered(v,x->IsBound(x));

            return v;
        end;
        ###############################################################
    fi;



    ###################################################################
    if IsBound(P!.action) and not prime=2 then 
        Action:=P!.action;
    else 
        Action:=function(k,j,g) return 1; end;
    fi;
    ###################################################################

    MultRecord:=[];
    
    ###################################################################
    
    ###################################################################     
    #1
    #F  Mult
    ##
    ##  The product of 2 elements in EltsG 
    ##
    ##  Input:  a pair of positive integers (g,h)
    ##  Output: the position of the product in the list EltsG  
    ##
    Mult:=function(g,h)
    local pos;
        if not IsBound(MultRecord[g]) then 
            MultRecord[g]:=[];
        fi;
        if not IsBound(MultRecord[g][h]) then
            pos:= Position(EltsG,EltsG[g]*EltsG[h]);
            if pos=fail then 
                Add(EltsG,EltsG[g]*EltsG[h]);  
                MultRecord[g][h]:= Length(EltsG);
            else 
                MultRecord[g][h]:= pos; 
            fi;
        fi;
    return MultRecord[g][h];
    end;
    #######################end of Mult#################################

    # Create a list of stabiliser subgroups and their resolutions

    StabGrps:= List([0..Length(P)],n->
                List([1..P!.dimension(n)], k->P!.stabilizer(n,k))); 

    StabResls:=[];
    i:=N;
    if prime=0 then

        for L in StabGrps do
            Add(StabResls,List(L,    
                g->ExtendScalars(ResolutionFG(g,i),G,EltsG))
                );
            i:=Maximum(0,AbsInt(i-1));
        od;
    else
        for L in StabGrps do
            Add(StabResls,List(L,
            g->ExtendScalars(ResolutionFiniteGroup(g,i,false,prime),
            G,EltsG)));
            i:=Maximum(0,AbsInt(i-1));
        od;
    fi;

    DimRecs:=List([0..N],i->[]);

    ###################################################################
    
    ###################################################################     
    #1
    #F  Dimension
    ##
    ##  Find the ZG-rank of R_n 
    ##
    ##  Input:  a non-negative integer n
    ##  Output: ZG-rank of R_n  
    ##
    Dimension:=function(k)
    local 
        dim,i,R;
        dim:=0;
        for i in [0..k] do
            DimRecs[k+1][i+1]:=[];
            for R in StabResls[i+1] do
                dim:=dim+R!.dimension(k-i);
                Add(DimRecs[k+1][i+1],dim);
            od;
        od;
    return dim;
    end;
    ###################################################################
    
    DimensionRecord:=List([0..N],Dimension);

    Dimension:=function(k);
        return DimensionRecord[k+1];
    end;
    ###################################################################

    ###################################################################

    # Create a record for the function Quad2Pair
    
    Quad2PairRec:=[];
    for q in [0..N] do
        Quad2PairRec[q+1]:=[];
        for r in [1..Length(StabGrps[q+1])] do
            Quad2PairRec[q+1][r]:=[];
            for s in [0..N-q] do
                Quad2PairRec[q+1][r][s+1]:=[];
            od;
        od;
    od;
    ###################################################################


    ###################################################################     
    #1
    #F  Pair2Quad
    ##
    ##  The n-th generator in degree k of our final resolution is
    ##  actually the t-th generator in degree s of the resolution
    ##  of the r-th stabilizer group of the q-th chain module of the
    ##  non-free resolution. We need the function f(k,n)=[q,r,s,t]
    ##
    ##  Input:  a pair of integers (k,n)
    ##  Output: [q,r,s,t]  
    ##
    Pair2Quad:=function(k,n)
    local 
        qq,q,r,s,t;

        for qq in [0..N] do
            if n <= DimRecs[k+1][qq+1][Length(DimRecs[k+1][qq+1])] then 
                q:=qq; break; 
            fi;
        od;

        r:=PositionProperty(DimRecs[k+1][q+1],x->(n<=x));
        s:=k-q;
        if r-1>0 then
            t:=n-DimRecs[k+1][q+1][r-1];
        else
            if q>=1 then
                t:=n-DimRecs[k+1][q][Length(  DimRecs[k+1][q] )];;
            else 
                t:=n;
            fi;
        fi;
        Quad2PairRec[q+1][r][s+1][t]:=[k,n];
    return [q,r,s,t];
    end;
    ####################end of Pair2Quad ##############################

    # Create a record for the function Pair2Quad 

    Pair2QuadRec:=List([1..N+1],i->[]);
    for k in [0..N] do
        for n in [1..Dimension(k)] do
            Pair2QuadRec[k+1][n]:=Pair2Quad(k,n);
        od;
    od;
    ###################################################################

    ###################################################################     
    #1
    #F  Pair2Quad
    ##
    ##  Input:  a pair of integers (k,n)
    ##  Output: [q,r,s,t]  
    ##
    Pair2Quad:=function(k,n)
    local 
        a;
        if n>0 then
            return StructuralCopy(Pair2QuadRec[k+1][n]);
        else
            a:=StructuralCopy(Pair2QuadRec[k+1][-n]);
            a[4]:=-a[4];
            return a;
        fi;
    end;
    #################end of Pair2Quad##################################
    
    ###################################################################     
    #1
    #F  Quad2Pair 
    ##
    ##  Input:  [q,r,s,t]
    ##  Output: [k,n]  
    ##
    Quad2Pair:=function(q,r,s,t)
    local 
        a,pr,pt;
        if r>0 then 
            pr:=r;pt:=t;
        else
            pr:=-r;pt:=-t;
        fi;

        if pt>0 then
            return StructuralCopy(Quad2PairRec[q+1][pr][s+1][pt]);
        else
            a:=StructuralCopy(Quad2PairRec[q+1][pr][s+1][-pt]);
            a[2]:=-a[2];
            return a;
        fi;
    end;
    ######################end of Quad2pair ############################
    
    ###################################################################     
    #1
    #F  HtpyGen
    ##
    ##  This applies the "vertical homotopy" to the free group 
    ##  generator [r,t,g] in "dimension" [q,s].
    ##  The output is an "r-word" in "dimension" [q,s+1]. 
    ##
    ##  Input:  a free generato [r,t,g] in demension [q,s]
    ##  Output: an "r-word" in "dimension" [q,s+1]  
    ##
    HtpyGen:=function(q,s,r,t,g)
    local 
        y,pr,pt;

        if r>0 then 
            pr:=r;pt:=t;
        else
            pr:=-r;pt:=-t;
        fi;

        y:=StructuralCopy(StabResls[q+1][pr]!.homotopy(s,[pt,g]));
        Apply(y,x->[pr,x[1],x[2]]);
    return y;
    end;
    #######################end of HtpyGen##############################

    ###################################################################     
    #1
    #F  HtpyWord
    ##
    ##  This applies the "vertical homotopy" to the r-word w in 
    ##  “dimension” [q,s]. The output is an r-word in "dimension" [q,s+1].
    ##
    ##  Input: r-word w in “dimension” [q,s]
    ##  Output: an r-word in "dimension" [q,s+1] 
    ##
    HtpyWord:=function(q,s,w)
    local 
        h,z,x,y;
    #This applies the "vertical homotopy" to the r-word w in "dimension" 
    #[q,s]. The output is an r-word in "dimension" [q,s+1].

        h:=[];
        for y in w do
            x:=[Action(q,y[1],y[3])*y[1],y[2],y[3]];
            z:=HtpyGen(q,s,x[1],x[2],x[3]);
            z:=List(z,a->[Action(q,a[1],a[3])*a[1],a[2],a[3]]);
            Append(h,z);
        od;

    return AlgRed(h);
    end;
    #####################end of HtpyWord ##############################

    # Create a record for function DelGen

    DelGenRec:=[];
    for k in [1..N+1] do
        DelGenRec[k]:=[];
        for q in [1..N+1] do
            DelGenRec[k][q]:=[];
            for s in [1..N+1] do
                DelGenRec[k][q][s]:=[];
                for r in [1..P!.dimension(q-1)] do
                    DelGenRec[k][q][s][r]:=[];
                od;
            od;
        od;
    od;
    ###################################################################

    ###################################################################     
    #1
    #F  DelGen
    ##
    ##  For k=0,1,2 ... this is the equivariant homomorphism
    ##  Del_k:A_{q,s} ---> A_{q-k,s+k-1} applied to a free
    ##  r-generator [r,t] in dimension [q,s]. 
    ##
    ##  Input:  a positive integer k and 
    ##          a free r-generator [r,t] in dimension [q,s]
    ##  Output: a word in dimension [q-k,s+k-1].  
    ##
    DelGen:=function(k,q,s,r,t)
    local 
        y,pr,pt,i;

        if r>0 then 
            pr:=r;pt:=t;
        else
            pr:=-r;pt:=-t;
        fi;

        if IsBound(DelGenRec[k+1][q+1][s+1][pr][AbsInt(pt)]) then

            if pt>0 then return 
                DelGenRec[k+1][q+1][s+1][pr][pt];
            else
                return List(DelGenRec[k+1][q+1][s+1][pr][-pt], 
                                        a->[a[1],-a[2],a[3]]);
            fi;

        fi;

        if k=0 then
            if s=0 then 
                return [];
            else
                y:=List(StabResls[q+1][pr]!.boundary(s,pt),
                        x->[Action(q,r,x[2])*x[1],x[2]]);
                if pt>0 then
                    DelGenRec[k+1][q+1][s+1][pr][pt]:= 
                            AlgRed(List(y,x->[pr,x[1],x[2]]));
                    return DelGenRec[k+1][q+1][s+1][pr][pt];
                else 
                    DelGenRec[k+1][q+1][s+1][pr][-pt]:=
                            AlgRed(List(y,x->[pr,-x[1],x[2]]));
                    return AlgRed(List(y,x->[pr,x[1],x[2]]));   
                fi;
            fi;
        fi;

        if k=1 then
            if s=0 then  
                if q=0 then 
                    return [];
                fi;
                y:=BoundaryP(q,pr);
                if pt>0 then 
                    DelGenRec[k+1][q+1][s+1][pr][pt]:= 
                                AlgRed(List(y,x->[x[1],1,x[2]]));
                    return DelGenRec[k+1][q+1][s+1][pr][pt];
                else
                    DelGenRec[k+1][q+1][s+1][pr][-pt]:= 
                                AlgRed(List(y,x->[x[1],1,x[2]]));
                    return List(y,x->[x[1],-1,x[2]]);
                fi;

            else

                if pt>0 then
                    DelGenRec[k+1][q+1][s+1][pr][pt]:=
                    AlgRed(HtpyWord(q-1,s-1,DelWord(1,q,s-1,
                                        DelGen(0,q,s,pr,-pt)))) ;
                    return DelGenRec[k+1][q+1][s+1][pr][pt];
                else
                    DelGenRec[k+1][q+1][s+1][pr][-pt]:=
                    AlgRed(HtpyWord(q-1,s-1,DelWord(1,q,s-1,
                                        DelGen(0,q,s,pr,pt)))) ;

                    return 
                        List(DelGenRec[k+1][q+1][s+1][pr][-pt], 
                                        a->[a[1],-a[2],a[3]]);
                fi;

            fi;
        fi;

        y:=[];
        for i in [1..k] do
            Append(y,
            HtpyWord(q-k,s+k-2,DelWord(i,q-k+i,s+k-i-1,
                                        DelGen(k-i,q,s,pr,-pt)))
            );
        od;
        y:=AlgRed(y);

        if pt>0 then
            DelGenRec[k+1][q+1][s+1][pr][pt]:=y;
        else
            DelGenRec[k+1][q+1][s+1][pr][-pt]:=List(y,a->
                                         [a[1],-a[2],a[3]]);
        fi;

    return y;
    end;
    #########################end of DelGen  ###########################

    ###################################################################     
    #1
    #F  DelWord
    ##
    ##  The product of 2 elements in Elts 
    ##
    ##  Input:  a pair of positive integers (i,j)
    ##  Output: the position of the product in the list Elts  
    ##
    DelWord:=function(k,q,s,w)
    local 
        y,x;
    #For k=0,1,2 ... this is the equivariant homomorphism
    #Del_k:A_{q,s} ---> A_{q-k,s+k-1} applied to an r-word [[r,t,g],...]
    #in dimension [q,s].

        y:=[];
        for x in w do
            Append(y,List(DelGen(k,q,s,x[1],x[2]),
                                    a->[a[1],a[2],Mult(x[3],a[3])]));
        od;

        return y;  #Added Jan 2013. Speeds up the calculation 
                   ## in some(!!) examples.
        return AlgRed(y);

    end;
    ###################################################################

    ###################################################################     
    #1
    #F  Boundary
    ##
    ##  Boundary map d_k: C_k -> C_{k-1} 
    ##
    ##  Input:  a pair of positive integers (k,n)
    ##  Output: The boundary d_k(f_n)  
    ##
    Boundary:=function(k,n)
    local 
        q,s,r,t,x,y,z,i;
        y:=Pair2Quad(k,n); 
        q:=y[1];s:=y[3];r:=y[2];t:=y[4];

        y:=[];

        for i in [0..k] do
            if q>=i then
                z:=DelGen(i,q,s,r,t);
                Append(y,List(z,x->
                        [Quad2Pair(q-i,x[1],s+i-1,x[2])[2],x[3]])  );
            else 
                break;
            fi;
        od;

    return AlgebraicReduction(y);
    end;
    ###################################################################

    PseudoBoundary:=[];
    for n in [1..N+1] do
        PseudoBoundary[n]:=[];
    od;

    ###################################################################     
    #1
    #F  FinalBoundary
    ##
    ##  Boundary map d_k: C_k -> C_{k-1}  
    ##
    ##  Input:  a pair of positive integers (k,n)
    ##  Output: the boundary d_k(f_n)  
    ##
    FinalBoundary:=function(n,k)
    local  
        pk;
        pk:=AbsInt(k);
        if  not IsBound(PseudoBoundary[n+1][pk]) then
            PseudoBoundary[n+1][pk]:= Boundary(n,pk);
        fi;
        if k>0 then  
            return PseudoBoundary[n+1][k];
        else  
            return NegateWord(PseudoBoundary[n+1][pk]); 
        fi;
    end;
    ###################################################################



    ################spectral sequence requirements#####################    


    FiltDimRecs:=[];
    for k in [0..N] do
        FiltDimRecs[k+1]:=[];
        for i in [1..Dimension(k)] do
            FiltDimRecs[k+1][i]:=Pair2Quad(k,i)[1];
        od;
    od;

    FilteredLength:=Maximum(Flat(FiltDimRecs));

    ###################################################################
    FilteredDimension:=function(r,k);

        return Length(Filtered(FiltDimRecs[k+1],x->x<=r));

    end;
    ###################################################################



    ########### BUI ANH TUAN - FEB 2014 ###############################



    ###################################################################     
    #1
    #F  InducedHtpyGen
    ##
    ##  This constructs the "induced homotopy" h1 from
    ##  the given homotopy of the non-free complex
    ##   h1: A_qs ->A_{q+1}s
    ##  This applies the "induced homotopy" to the free group generator
    ##   [r,t,g] in "dimension" [q,s].
    ##  The output is an "r-word" in "dimension" [q+1,s].
    ##
    ##  Input:  five integers
    ##  Output: an r-word in dimension [q+1,s].  
    ##
    InducedHtpyGen:=function(q,s,r,t,g)
    local 
        y,pr,pt,w,v;

        if r>0 then 
            pr:=r;pt:=t;
        else
            pr:=-r;pt:=-t;
        fi;

        y:=StructuralCopy(P!.homotopy(q,[pr,g]));
        if pt>0 then 
            Apply(y,x->[x[1],1,x[2]]);
        else 
            Apply(y,x->[x[1],-1,x[2]]);
        fi;
    return y;
    end;
    ###################################################################
    
    ###################################################################     
    #1
    #F  InducedHtpyWord
    ##
    ##  This applies the "induced homotopy" to the r-word w in  
    ##  "dimension" [q,s]. The output is an r-word in 
    ##  "dimension" [q+1,s].
    ##
    ##  Input:  dimension q,s and a word w
    ##  Output: an r-word in dimension [q+1,s]. 
    ##
    InducedHtpyWord:=function(q,s,w)
    local 
        h,z,x,y;


        h:=[];
        for y in w do
            x:=[y[1],y[2],y[3]];
            z:=StructuralCopy(InducedHtpyGen(q,s,x[1],x[2],x[3]));
            z:=List(z,a->[a[1],a[2],a[3]]);
            Append(h,z);
        od;

    return AlgRed(h);
    end;
    ###################################################################

    ###################################################################     
    #1
    #F  InducedHtpyList
    ##
    ##  This applies the Horizontal Homotopy to a list of words
    ##  For each word, this applies the "induced homotopy" to  
    ##  the r-word w in "dimension" [q,s]. 
    ##  The output is an r-word in "dimension" [q+1,s].
    ##
    ##  Input:  an r-word w in dimension [q,s]
    ##  Output: an r-word in dimension [q+1,s].  
    ##    
    InducedHtpyList:=function(w)
    local 
        h,z,x,y,v,b;
       
        h:=[];
        for y in w do 
            z:=StructuralCopy(InducedHtpyGen(y[1],y[2],y[3],y[4],y[5]));
            z:=List(z,a->[y[1]+1,y[2],a[1],a[2],a[3]]);
            Append(h,z);
        od;
        return h;
    end;
    ################################################################### 

    ############# d+=d1+d2+..+dq of a list of words####################
    
    ###################################################################     
    #1
    #F  DelListSum
    ##
    ##  Sum of DelWord_k where k from 1 to q
    ##
    ##  Input:  a word w
    ##  Output: the differential map d(w)  
    ##

    DelListSum:=function(w)
    local 
        y,d,x,h,k;
        h:=[];
        for x in w do
            y:=[];
            for k in [1..x[1]] do
                d:=StructuralCopy(DelGen(k,x[1],x[2],x[3],x[4]));
                Apply(d,v->[x[1]-k,x[2]+k-1,v[1],v[2],Mult(x[5],v[3])]);        
                Append(y,d);
            od;
            Append(h,y);
        od;

    return h;
    end;
    ################################################################### 
    
    ###################################################################     
    #1
    #F  VertHtpy
    ##
    ##  Applies to a list of [q,s,r,t,g]: could be in different A_qs
    ##  return a list of elements of the form [q,s,r,t,g]
    ##
    ##  Input:  [q,s,r,t,g]
    ##  Output: a list of elements of the form [q,s,r,t,g]. 
    ##
    VertHtpy:=function(w)

    local 
        h,x,y,v;
        h:=[];;
        for x in w do
            v:=[x[1],x[2],Action(x[1],x[3],x[5])*x[3],x[4],x[5]];
            y:=StructuralCopy(HtpyGen(v[1],v[2],v[3],v[4],v[5]));
            Apply(y,a->[x[1],x[2]+1,Action(x[1],a[1],a[3])*a[1],a[2],a[3]]);
            Append(h,y);
        od;    
    return h;
    end;
    ###################################################################  
    
    ###################################################################     
    #1
    #F  NegateListWord 
    ##
    ##  Input:  a list of elements of the form [q,s,r,t,g]
    ##  Output: the negated word.  
    ##

    NegateListWord:=function(w)
        Apply(w,x->[x[1],x[2],-x[3],x[4],x[5]]);    
    return w;
    end;
    ###################################################################
    
    ###################################################################     
    #1
    #F  HomotopyGen
    ##
    ##  The product of 2 elements in Elts 
    ##
    ##  Input:  [q,s,r,t,g]
    ##  Output: the homotopy h(q,s,r,t,g).  
    ##
    HomotopyGen:=function(arg)
    local 
        f,g,q,s,r,t,x,e,v,y,
        h0, h1, h0dh1, e3, h2, h;

        q:=arg[1];
        s:=arg[2];
        r:=arg[3];    
        t:=arg[4];
        g:=arg[5];

        if arg=[] then 
            return [];
        else
            if s = 0 then 
                h1:=StructuralCopy(InducedHtpyList([[q,s,r,t,g]]));
                h0dh1:=StructuralCopy(VertHtpy(DelListSum(h1)));
                v:=StructuralCopy(DelListSum(h0dh1));
                e3:=[];         # e3=h(d+)h0(d+)d1
                for x in v do
                    Append(e3,HomotopyGen(x[1],x[2],x[3],x[4],x[5]));
                od;
                e:=[];          # e=h1-h0(d+)h1+h(d+)h0(d+)h1
                Append(e,h1);
                Append(e,NegateListWord(h0dh1));
                Append(e,e3);
            elif s>0 then
            # s>0 then e=0
                e:=[];
            fi;
            y:=StructuralCopy([q,s,r,t,g]);

            h0:=VertHtpy([y]);

            v:=DelListSum(h0);
            h2:=[];         # h2=h(d+)h0
            for x in v do
                Append(h2,HomotopyGen(x[1],x[2],x[3],x[4],x[5]));
            od;
            h:=[];         # h=h0-d(d+)h0 + e
            Append(h,h0);
            Append(h,NegateListWord(h2));
            Append(h,e);
            return h;
        fi;
    end;

    
    ###################################################################     
    #1
    #F  Homotopy
    ## 
    ##
    ##  Input:  degree k an a word w
    ##  Output: the homotopy h_k(w).  
    ##

    Homotopy:=function(k,w)
    local 
        f,g,q,s,r,t,v,e,h;

    ### h([])=[] 
        if w=[] then 
            return [];
        fi;
        f:=w[1];
        g:=w[2];
        v:=Pair2Quad(k,f);
        q:=v[1];
        r:=v[2];
        s:=v[3];
        t:=v[4];

        h:=HomotopyGen(q,s,r,t,g);
        Apply(h,x->[Quad2Pair(x[1],x[3],x[2],x[4])[2],x[5]]);
    return AlgebraicReduction(h);
    end;
    ###################################################################

SetInfoLevel(InfoWarning,1);

#######################################################################
return Objectify(HapResolution,
                rec(
                inputresl:=P,
                verthtpy:=VertHtpy,
                htpy:=HtpyWord,
                delword:=DelWord,
                dimension:=Dimension,
                filteredDimension:=FilteredDimension,
                boundary:=FinalBoundary,
        inducedhomotopy:=InducedHtpyGen,
        stabrels:=StabResls,
                homotopy:=Homotopy,
                elts:=P!.elts,
                group:=P!.group,
                pseudoBoundary:=PseudoBoundary,
                properties:=
                   [["length",N],
                    ["filtration_length",FilteredLength],
                    ["initial_inclusion",true],
                    ["reduced",true],
                    ["type","resolution"],
                    ["characteristic",prime]  ]));


end);

################### end of FreeZGResolution ############################

