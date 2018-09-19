##################################################
InstallGlobalFunction(LieExteriorSquare,
function(L)
local 
      BasisL, lenBL, SCTL, lisabelian, K, SCTC, u1, u2, u3, u4, l1, l2,
      index1, index2, C, t1, t2, t3, t4, t5, t6, i, j, k, m, p,
      vectorsI, vectorsII,w1,w2,
      zr1, mr1, zr2, mr2, derayeh1, derayeh2,
      lenv1, lenv2, lenv3, lenv4, lenv5, lenv6, 
      e, q, w, v1, v2, lzr1, lzr2, lmr1, lmr2,
      q1, q2, tt1, tt2, ww, h, LTL, LVL, t, v7, I, II,
      bL, bLTL, g, MLTL, v, vv1, BI, I1, I2, I3,
      LenBI, BLTL, BBLTL, liltl, llzr1,
      MLVL,BLVL,BBLVL,II1,II2,II3,BII,LenBII,lilvl,bLV,
      LL,f,g1,Pair;


lisabelian:=0;
if IsLieAbelian(L) then lisabelian:=1; fi; 
K:=L!.LeftActingDomain;
BasisL:=Basis(L);
lenBL:=Length(BasisL);
SCTL:=StructureConstantsTable(BasisL); 
SCTC:=EmptySCTable(lenBL^2,0*One(K),"antisymmetric");

for u1 in [1..lenBL] do 
 for u2 in [1..lenBL] do
  for u3 in [1..lenBL] do
   for u4 in [1..lenBL] do
     l1:=Length(SCTL[u1][u2][1]);
     l2:=Length(SCTL[u3][u4][1]);
     if l1<>0 then if l2<>0 then 
                         index1:=(u1-1)*lenBL+u2;
                         index2:=(u3-1)*lenBL+u4;
                         derayeh1:=[];
 
                         for t1 in [1..l1] do
                          for t2 in [1..l2] do
                            i:=SCTL[u1][u2][1][t1];
                            j:=SCTL[u3][u4][1][t2];
                            m:=(i-1)*lenBL+j;
                            p:=SCTL[u1][u2][2][t1]*SCTL[u3][u4][2][t2];
                            Add(derayeh1,p*One(K));
                            Add(derayeh1,m);
                          od;
                         od;  
        
                         SetEntrySCTable(SCTC,index1,index2,derayeh1);
     fi; fi;
   od;
  od;
 od;
od;

C:=AlgebraByStructureConstants(K,SCTC);

vectorsI:=[];

for u1 in [1..lenBL] do 
 for u2 in [1..lenBL] do
  for u3 in [1..lenBL] do
    zr1:=[];
    mr1:=[];
    zr2:=[];
    mr2:=[];
    lenv1:=Length(SCTL[u1][u2][1]);
    if lenv1<>0 then
      for t1 in [1..lenv1] do
        i:=SCTL[u1][u2][1][t1];
        j:=u3;
        m:=(i-1)*lenBL+j;  
        p:=SCTL[u1][u2][2][t1]*1;
        Add(zr1,p*One(K)); 
        Add(mr1,m);
      od;
    fi;
    lenv2:=Length(SCTL[u2][u3][1]);
    if lenv2<>0 then
      for t2 in [1..lenv2] do
        i:=u1;
        j:=SCTL[u2][u3][1][t2];
        m:=(i-1)*lenBL+j;  
        p:=1*SCTL[u2][u3][2][t2];
        e:=0;
        if Length(zr1)<>0 then 
          for q in [1..Length(zr1)] do
            if mr1[q]=m then zr1[q]:=zr1[q]-p*One(K); e:=1; fi;
          od;      
          if e=0 then Add(zr1,-1*p*One(K)); Add(mr1,m); fi;
        fi;
        if Length(zr1)=0 then  Add(zr1,-1*p*One(K)); Add(mr1,m); fi;
      od;
    fi;
    lenv3:=Length(SCTL[u1][u3][1]);
    if lenv3<>0 then
      for t3 in [1..lenv3] do
        i:=u2;
        j:=SCTL[u1][u3][1][t3];
        m:=(i-1)*lenBL+j;  
        p:=1*SCTL[u1][u3][2][t3];
        e:=0; 
        if Length(zr1)<>0 then 
          for q in [1..Length(zr1)] do
            if mr1[q]=m then zr1[q]:=zr1[q]+p*One(K); e:=1; fi;
          od;
          if e=0 then Add(zr1,p*One(K)); Add(mr1,m); fi;
        fi;
        if Length(zr1)=0 then Add(zr1,p*One(K)); Add(mr1,m); fi;
      od;
    fi;
    lenv4:=Length(SCTL[u2][u3][1]);
    if lenv4<>0 then
      for t4 in [1..lenv2] do
        i:=u1;
        j:=SCTL[u2][u3][1][t4];
        m:=(i-1)*lenBL+j;  
        p:=1*SCTL[u2][u3][2][t4];
        Add(zr2,p*One(K)); 
        Add(mr2,m);
      od;
    fi;
    lenv5:=Length(SCTL[u3][u1][1]);
    if lenv5<>0 then
      for t5 in [1..lenv5] do
        i:=SCTL[u3][u1][1][t5];
        j:=u2; 
        m:=(i-1)*lenBL+j;  
        p:=SCTL[u3][u1][2][t5]*1;
        e:=0;
        if Length(zr2)<>0 then 
          for q in [1..Length(zr2)] do
            if mr2[q]=m then zr2[q]:=zr2[q]-p*One(K); e:=1; fi;
          od;      
          if e=0 then Add(zr2,-1*p*One(K)); Add(mr2,m); fi;
        fi;
        if Length(zr2)=0 then  Add(zr2,-1*p*One(K)); Add(mr2,m); fi;
      od;
    fi;
    lenv6:=Length(SCTL[u2][u1][1]);
    if lenv6<>0 then
      for t6 in [1..lenv6] do
        i:=SCTL[u2][u1][1][t6];
        j:=u3;
        m:=(i-1)*lenBL+j;  
        p:=SCTL[u2][u1][2][t6]*1;
        e:=0; 
        if Length(zr2)<>0 then 
          for q in [1..Length(zr2)] do
            if mr2[q]=m then zr2[q]:=zr2[q]+p*One(K); e:=1; fi;
          od;
          if e=0 then Add(zr2,p*One(K)); Add(mr2,m); fi;
        fi;
        if Length(zr2)=0 then Add(zr2,p*One(K)); Add(mr2,m); fi;
      od;
    fi;

    lzr1:=Length(zr1);
    lmr1:=Length(mr1);
    v1:=0;
    if lzr1<>0 then 
      v1:=zr1[1]*(Elements(Basis(C))[lenBL^2-mr1[1]+1]);
      for q1 in [2..lzr1] do
        if zr1[q1]<>0 then  
          v1:=v1+zr1[q1]*(Elements(Basis(C))[lenBL^2-mr1[q1]+1]);
        fi;
      od;
    fi;
    tt1:=0;
    for ww in [1..Length(Basis(C))] do
      if v1=0*Elements(Basis(C))[ww] then tt1:=1; fi;
    od; 
    w:=0; 
    for h in [1..Length(vectorsI)] do
      if v1=vectorsI[h] then w:=1; fi;
      if -1*v1=vectorsI[h] then w:=1; fi;
    od;
    if v1<>0 then if w=0 then if tt1=0 then Add(vectorsI,v1); fi; fi; fi;
    lzr2:=Length(zr2);
    lmr2:=Length(mr2);
    v2:=0;
    if lzr2<>0 then 
      v2:=zr2[1]*(Elements(Basis(C))[lenBL^2-mr2[1]+1]);
      for q2 in [2..lzr2] do
        if zr2[q2]<>0 then  
          v2:=v2+zr2[q2]*(Elements(Basis(C))[lenBL^2-mr2[q2]+1]);
        fi;
      od;
    fi;
    tt2:=0;
    for ww in [1..Length(Basis(C))] do
      if v2=0*Elements(Basis(C))[ww] then tt2:=1; fi;
    od; 
    w:=0; 
    for h in [1..Length(vectorsI)] do
      if v2=vectorsI[h] then w:=1; fi;
      if -1*v2=vectorsI[h] then w:=1; fi;
    od;
    if v2<>0 then if w=0 then if tt2=0 then Add(vectorsI,v2); fi; fi; fi;

  od;
 od;
od;

I:=Ideal(C,vectorsI);
LTL:=C/I;
#g1:=NaturalHomomorphismByIdeal(C,I);
#LTL:=Image(g1);
vectorsII:=[];

for t in [1..Length(vectorsI)] do
  Add(vectorsII,vectorsI[t]);
od;

for i in [1..lenBL] do 
  w:=(i-1)*lenBL+i;
  v7:=Elements(Basis(C))[lenBL^2-w+1];
  Add(vectorsII,v7);
od;

for i in [1..lenBL] do 
 for j in [1..lenBL] do 
   w1:=(i-1)*lenBL+j;
   w2:=(j-1)*lenBL+i;
   v7:=Elements(Basis(C))[lenBL^2-w1+1]+Elements(Basis(C))[lenBL^2-w2+1];
   Add(vectorsII,v7);
 od;
od;

II:=Ideal(C,vectorsII);
#LVL:=C/II;
g1:=NaturalHomomorphismByIdeal(C,II);
LVL:=Image(g1);
MLTL:=[];

for i in [1..Length(Basis(C))] do
  v:=Elements(Basis(C))[lenBL^2-i+1];
  if not (v in I) then Add(MLTL,v); fi;
od;

BI:=Basis(I);
LenBI:=Length(BI);
BLTL:=[];
BBLTL:=[];

for j in [1..LenBI] do
  Add(BBLTL,Elements(BI)[j]);
od;

for i in [1..Length(MLTL)] do
  Add(BLTL,MLTL[i]);
  Add(BBLTL,MLTL[i]);
  I1:=VectorSpace(K,BBLTL);
  if Dimension(I1)=LenBI+1 then LenBI:=LenBI+1;  
    else 
      Remove(BLTL);
      Remove(BBLTL);
  fi;
  if Length(BLTL)=Dimension(C) then break; fi;
od;

I2:=VectorSpace(K,BLTL);
I3:=VectorSpace(K,Basis(LTL));
liltl:=[];

for i in [1..Length(BLTL)] do
 for j in [1..Length(Basis(C))] do
   if BLTL[i]=Elements(Basis(C))[lenBL^2-j+1] then 
     Add(liltl,j);         
     break;
   fi;
 od;
od; 

bL:=[];

for k in [1..Length(liltl)] do
  u1:=Int(liltl[k]/lenBL)+1;  
  u2:= liltl[k] mod lenBL;
  if u2=0 then u2:=lenBL; u1:=u1-1;  fi;
  u3:=SCTL[u1][u2];
  llzr1:=Length(SCTL[u1][u2][1]);
  vv1:=0*Elements(BasisL)[1];
  if llzr1<>0 then 
    vv1:=SCTL[u1][u2][2][1]*Elements(BasisL)[lenBL-SCTL[u1][u2][1][1]+1];
    for i in [2..llzr1] do
      vv1:=vv1+SCTL[u1][u2][2][i]*Elements(BasisL)[lenBL-SCTL[u1][u2][1][i]+1];
    od;
  fi;
  Add(bL,vv1);
od;

bLTL:=Basis(LTL);

g:= AlgebraHomomorphismByImages( LTL, L, bLTL , bL );;

MLVL:=[];

for i in [1..Length(Basis(C))] do
  v:=Elements(Basis(C))[lenBL^2-i+1];
  if not (v in II) then Add(MLVL,v); fi;
od;

BII:=Basis(II);
LenBII:=Length(BII);
BLVL:=[];
BBLVL:=[];

for j in [1..LenBII] do
  Add(BBLVL,Elements(BII)[j]);
od;

for i in [1..Length(MLVL)] do
  Add(BLVL,MLVL[i]);
  Add(BBLVL,MLVL[i]);
  II1:=VectorSpace(K,BBLVL);
  if Dimension(II1)=LenBII+1 then LenBII:=LenBII+1;  
    else 
      Remove(BLVL);
      Remove(BBLVL);
  fi;
  if Length(BLVL)=Dimension(C) then break; fi;
od;

II2:=VectorSpace(K,BLVL);
II3:=VectorSpace(K,Basis(LVL));
lilvl:=[];

for i in [1..Length(BLVL)] do
 for j in [1..Length(Basis(C))] do
   if BLVL[i]=Elements(Basis(C))[lenBL^2-j+1] then 
     Add(lilvl,j);         
     break;
   fi;
 od;
od; 

bLV:=[];

for k in [1..Length(lilvl)] do
  u1:=Int(lilvl[k]/lenBL)+1;  
  u2:= lilvl[k] mod lenBL;
  if u2=0 then u2:=lenBL; u1:=u1-1;  fi;
  u3:=SCTL[u1][u2];
  llzr1:=Length(SCTL[u1][u2][1]);
  vv1:=0*Elements(BasisL)[1];
  if llzr1<>0 then 
    vv1:=SCTL[u1][u2][2][1]*Elements(BasisL)[lenBL-SCTL[u1][u2][1][1]+1];
    for i in [2..llzr1] do
      vv1:=vv1+SCTL[u1][u2][2][i]*Elements(BasisL)[lenBL-SCTL[u1][u2][1][i]+1];
    od;
  fi;
  Add(bLV,vv1);
od;

BLVL:=Basis(LVL);
LL:=LieDerivedSubalgebra(L);
#####  f: hom(LwL ---> [L,L])  ##################
f:= AlgebraHomomorphismByImages( LVL, LL, BLVL , bLV );
#g1:= NaturalHomomorphismByIdeal( C, II );

################################
Pair:=function(x,y)

  local  
    z1,z2,z3,z4,l,v1,v,k,i,j;

z1:=Coefficients( Basis(L), x );
z2:=Coefficients( Basis(L), y );
l:=Length(z1);
v:=0*Random(Images(g1,Basis(C)[1]));
for i in [1..l] do
 for j in [1..l] do
   z3:=z1[i]*z2[j];
   if z3<>0 then
     k:=(i-1)*lenBL+j;
     v1:=Random(Images(g1,Basis(C)[k]));
     v:=v+z3*v1;    
   fi;  
 od;
od;
return v;

end;
################################

return rec(homomorphism:=f, pairing:=Pair );
end);
