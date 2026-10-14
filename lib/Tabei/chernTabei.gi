#############################################################################
##
##  chernTabei.gi                HAP library                  Moe Tabei
##
##  Second Chern classes of complex representations of a finite group, for
##  representations all of whose irreducible constituents have degree at most
##  two and are monomial (always the case for nilpotent groups), computed as
##  Euler classes of representation spheres.
##
##  Contributed to HAP by Moe Tabei <tabei@ryun.jp>, 2026.  MIT licence.
##  Development version, tests and calibration:
##      https://github.com/yui9696/hap-crosschecks/tree/main/chern
##  The functions are declared in decsTabei.gd.
##
#############################################################################
##
##  THE METHOD
##
##  Let W be a two-dimensional complex representation of a finite group H.
##  Then c_2(W) is the top Chern class of W, hence the Euler class of the
##  underlying oriented real 4-plane bundle, and that class in
##
##      H^4(H,Z)  =  Ext^4_{ZH}(Z,Z)
##
##  is the class of the exact sequence of cellular chains
##
##      0 -> Z -> C_3(S) -> C_2(S) -> C_1(S) -> C_0(S) -> Z -> 0
##
##  of an H-CW structure on the unit sphere S = S(W_R) = S^3.  If W is
##  MONOMIAL, with entries in mu_m wr Sigma_2, we may take
##
##      S^3  =  S^1 * S^1
##
##  the join of the two coordinate circles, each given the m-gon structure,
##  on which H acts cellularly.  The class is computed by lifting the identity
##  of Z to a chain map from a HAP free resolution R of H into C_*(S): the
##  composite
##
##      R_4 --d--> R_3 --f_3--> C_3(S) --> ker(C_3 -> C_2) = Z.[S^3]
##
##  is a 4-cocycle representing the Euler class up to a universal sign.
##
##  Every irreducible representation of degree two of a nilpotent group is
##  monomial, induced from a linear character of a subgroup of index two.
##
##  For a sum W = (+) W_j with rank(W_j) <= 2 the Whitney formula gives
##
##      c_2(W) = sum_j c_2(W_j)  +  sum_{j<k} c_1(W_j) c_1(W_k),
##
##  and c_1(W_j) c_1(W_k) = c_2(det W_j (+) det W_k) is again the Euler class
##  of a two-dimensional (diagonal) representation, so the whole computation
##  is made of Euler classes.
##
##  NOTHING here uses Bocksteins or cup products, which makes this route
##  independent of the usual one.
##
##  The identification of the Euler class with the class of the cellular
##  chain complex is classical for FREE actions (Swan's periodic
##  resolutions).  The actions here are cellular but not free; that case is
##  Proposition 2.3 of A. Guclukan and E. Yalcin, "The Euler class of a
##  subset complex", Quart. J. Math. 61 (2010), 43-68; the join structure is
##  not a strict H-CW structure, and the README in the repository above
##  explains why the proposition still applies.  The functions below are
##  also calibrated against an independent route (calib.g in the same
##  repository).
##
#############################################################################

##  A universal sign, fixed once and for all by the calibration in calib.g
##  (repository above).
##  It is read-only; a user who adopts the opposite orientation convention
##  should change the value here.
BindGlobal("HAPCHERN_EulerSign", 1);

##  Caches, keyed by the order m of the roots of unity and by the resolution.
BindGlobal("HAPCHERN_Cache", rec(join := [], cr4 := []));

#############################################################################
##
#F  HAPCHERN_RootExponent( <x>, <m> )
##
##  For x in mu_m, the a in [0..m-1] with x = E(m)^a.
##
BindGlobal("HAPCHERN_RootExponent", function(x, m)
  local d;
  d := DescriptionOfRootOfUnity(x);          # x = E(d[1])^d[2]
  if m mod d[1] <> 0 then
    Error("HAPCHERN_RootExponent: value is not an m-th root of unity");
  fi;
  return (d[2] * (m / d[1])) mod m;
end);

#############################################################################
##
#F  HAPCHERN_CharacterValue( <chi>, <g> )
##
BindGlobal("HAPCHERN_CharacterValue", function(chi, g)
  local cc;
  cc := ConjugacyClasses(UnderlyingGroup(chi));
  return chi[PositionProperty(cc, c -> g in c)];
end);

#############################################################################
##
#F  JoinComplexOfCyclicGroup( <m> )
##
##  The cellular chain complex of S^1 * S^1 = S^3, each circle carrying the
##  m-gon structure.  Cells of a circle are coded [d,k] with d = -1 (the empty
##  cell), 0 (the vertex at E(m)^k) or 1 (the edge from E(m)^k to E(m)^(k+1)).
##  A cell of the join is a pair [x,y], not both empty, of degree
##  deg(x) + deg(y) + 1.
##
##  Returns rec( m, cells, pos, bnd ) where bnd[n+1] is the matrix of
##  d : C_n -> C_{n-1} in the row-vector convention, with degree -1 the
##  augmentation to Z.
##
##  The construction is verified on the spot: d d = 0, the augmentation is
##  onto with every vertex mapping to 1, the reduced homology vanishes and is
##  torsion free in degrees 0,1,2, and H_3 has rank one.
##
InstallGlobalFunction(JoinComplexOfCyclicGroup, function(m)
  local circ, cells, n, x, y, deg, pos, bnd, i, c, terms, t, s, mat, res, dcirc;

  if IsBound(HAPCHERN_Cache.join[m]) then return HAPCHERN_Cache.join[m]; fi;

  circ := Concatenation([[-1,0]],
                        List([0..m-1], k -> [0,k]),
                        List([0..m-1], k -> [1,k]));
  cells := List([0..3], n -> []);
  for x in circ do
    for y in circ do
      deg := x[1] + y[1] + 1;
      if deg >= 0 then Add(cells[deg+1], [x,y]); fi;
    od;
  od;
  pos := List([0..3], n -> (c -> Position(cells[n+1], c)));

  ##  augmented boundary of a cell of one circle, as a list of [coeff, cell]
  dcirc := function(x)
    if x[1] = -1 then return []; fi;
    if x[1] = 0 then return [[1, [-1,0]]]; fi;
    return [[1, [0, (x[2]+1) mod m]], [-1, [0, x[2]]]];
  end;

  bnd := [];
  for n in [0..3] do
    mat := [];
    for c in cells[n+1] do
      x := c[1]; y := c[2];
      if n = 0 then res := [0]; else res := List(cells[n], z -> 0); fi;
      terms := [];
      for t in dcirc(x) do Add(terms, [t[1], [t[2], y]]); od;
      for t in dcirc(y) do Add(terms, [(-1)^(x[1]+1) * t[1], [x, t[2]]]); od;
      for t in terms do
        if n = 0 then
          res[1] := res[1] + t[1];
        else
          i := pos[n](t[2]);
          if i = fail then Error("JoinComplexOfCyclicGroup: boundary cell not found"); fi;
          res[i] := res[i] + t[1];
        fi;
      od;
      Add(mat, res);
    od;
    Add(bnd, mat);
  od;

  ##  self-checks
  for n in [1..3] do
    if not IsZero(bnd[n+1] * bnd[n]) then
      Error("JoinComplexOfCyclicGroup: d o d <> 0 in degree ", n);
    fi;
  od;
  if ForAny(bnd[1], r -> r <> [1]) then
    Error("JoinComplexOfCyclicGroup: augmentation is not the sum of coefficients");
  fi;
  for n in [0..2] do
    s := SmithNormalFormIntegerMat(bnd[n+2]);
    if Length(cells[n+1]) - RankMat(bnd[n+1]) <> RankMat(bnd[n+2]) then
      Error("JoinComplexOfCyclicGroup: not exact in degree ", n);
    fi;
    if ForAny(DiagonalOfMat(s), d -> not d in [0,1,-1]) then
      Error("JoinComplexOfCyclicGroup: torsion in degree ", n);
    fi;
  od;
  if Length(cells[4]) - RankMat(bnd[4]) <> 1 then
    Error("JoinComplexOfCyclicGroup: H_3 is not of rank one");
  fi;

  HAPCHERN_Cache.join[m] := rec(m := m, cells := cells, pos := pos, bnd := bnd);
  return HAPCHERN_Cache.join[m];
end);

#############################################################################
##
#F  HAPCHERN_JoinAction( <J>, <p>, <a> )
##
##  The action on C_*(S^1 * S^1) of a monomial element with permutation part
##  p and exponents a, as signed permutation matrices A[n+1] acting on the
##  right: g.v = v * A[n+1].
##
##  The sign of a cell is not written down a priori: it is the unique sign for
##  which g.(dc) = d(gc), and the function fails if neither +1 nor -1 works.
##
BindGlobal("HAPCHERN_JoinAction", function(J, p, a)
  local m, A, n, c, x, y, xs, ys, img, sgn, i, mat, cimg, k, v, row;
  m := J.m;
  cimg := function(x, j)
    if x[1] = -1 then return x; fi;
    return [x[1], (x[2] + a[j]) mod m];
  end;
  A := [];
  for n in [0..3] do
    mat := NullMat(Length(J.cells[n+1]), Length(J.cells[n+1]));
    for i in [1..Length(J.cells[n+1])] do
      c := J.cells[n+1][i]; x := c[1]; y := c[2];
      xs := cimg(x, 1); ys := cimg(y, 2);
      if p[1] = 1 then img := [xs, ys]; else img := [ys, xs]; fi;
      k := J.pos[n+1](img);
      if k = fail then Error("HAPCHERN_JoinAction: image cell not found"); fi;
      if n = 0 then
        sgn := 1;
      else
        v := J.bnd[n+1][i] * A[n];
        row := J.bnd[n+1][k];
        if v = row then sgn := 1;
        elif v = -row then sgn := -1;
        else Error("HAPCHERN_JoinAction: g.(dc) is not +- d(gc)"); fi;
      fi;
      mat[i][k] := sgn;
    od;
    Add(A, mat);
  od;
  for n in [1..3] do
    if A[n+1] * J.bnd[n+1] <> J.bnd[n+1] * A[n] then
      Error("HAPCHERN_JoinAction: not a chain map in degree ", n);
    fi;
  od;
  return A;
end);

#############################################################################
##
#F  HAPCHERN_FundamentalCycle( <J> )
##
##  The generator of H_3 = ker(d_3), normalised to have coefficient +1 on the
##  cell edge(1,0) (x) edge(2,0).
##
BindGlobal("HAPCHERN_FundamentalCycle", function(J)
  local N, z, i;
  N := NullspaceIntMat(J.bnd[4]);
  if Length(N) <> 1 then
    Error("HAPCHERN_FundamentalCycle: ker d_3 is not of rank one");
  fi;
  z := N[1];
  i := J.pos[4]([[1,0],[1,0]]);
  if AbsInt(z[i]) <> 1 then
    Error("HAPCHERN_FundamentalCycle: not primitive at the reference cell");
  fi;
  return z * z[i];
end);

#############################################################################
##
#F  MonomialModelOfMatrixFunction( <H>, <m>, <matfun> [, <checkallpairs>] )
##
##  A monomial model of a two-dimensional representation of H given by a
##  function matfun : H -> 2x2 cyclotomic matrices, acting on column vectors.
##  Entries must be m-th roots of unity or zero, with exactly one non-zero
##  entry in each column.
##
##  By default the homomorphism property is checked on all pairs of
##  generators and on the relations implied by the group multiplication of a
##  generating set; pass true as a fourth argument to check all |H|^2 pairs
##  instead (slow, and only worth it once per new construction).
##
InstallGlobalFunction(MonomialModelOfMatrixFunction, function(arg)
  local H, m, matfun, full, elts, mats, perm, expo, i, j, M, p, a, k, gens,
        gi, prod;
  H := arg[1]; m := arg[2]; matfun := arg[3];
  if Length(arg) >= 4 then full := arg[4]; else full := false; fi;

  elts := Elements(H);
  mats := List(elts, matfun);

  perm := []; expo := [];
  for i in [1..Length(elts)] do
    M := mats[i];
    p := []; a := [];
    for j in [1..2] do
      k := PositionProperty([1..2], r -> M[r][j] <> 0);
      if k = fail or Number([1..2], r -> M[r][j] <> 0) <> 1 then
        Error("MonomialModelOfMatrixFunction: the matrix of ", elts[i],
              " is not monomial");
      fi;
      p[j] := k;
      a[j] := HAPCHERN_RootExponent(M[k][j], m);
    od;
    if p[1] = p[2] then
      Error("MonomialModelOfMatrixFunction: the matrix of ", elts[i],
            " is not invertible");
    fi;
    Add(perm, p); Add(expo, a);
  od;

  if full then
    for i in [1..Length(elts)] do
      for j in [1..Length(elts)] do
        if mats[i] * mats[j] <> mats[Position(elts, elts[i]*elts[j])] then
          Error("MonomialModelOfMatrixFunction: not a homomorphism");
        fi;
      od;
    od;
  else
    ##  It suffices to check mats[g*x] = mats[g]*mats[x] for g in a generating
    ##  set and x in H.  Induction on the length of a positive word for h in
    ##  the generators (every element of a finite group is such a word) then
    ##  gives mats[h*x] = mats[h]*mats[x] for all h, the base case mats[1] = I
    ##  following from mats[g] = mats[g]*mats[1] and the invertibility of
    ##  mats[g], which the monomial check above has already established.
    gens := GeneratorsOfGroup(H);
    for gi in gens do
      i := Position(elts, gi);
      for j in [1..Length(elts)] do
        prod := Position(elts, gi * elts[j]);
        if mats[i] * mats[j] <> mats[prod] then
          Error("MonomialModelOfMatrixFunction: not a homomorphism");
        fi;
      od;
    od;
  fi;

  return rec(H := H, m := m, elts := elts, mat := mats, perm := perm,
             expo := expo);
end);

#############################################################################
##
#F  MonomialModelOfLinearCharacterPair( <H>, <mu1>, <mu2> )
##
##  The diagonal model of the two-dimensional representation mu1 (+) mu2.
##
InstallGlobalFunction(MonomialModelOfLinearCharacterPair, function(H, mu1, mu2)
  return MonomialModelOfMatrixFunction(H, Exponent(H),
    h -> [[HAPCHERN_CharacterValue(mu1,h), 0],
          [0, HAPCHERN_CharacterValue(mu2,h)]]);
end);

#############################################################################
##
#F  HAPCHERN_IndexTwoSubgroups( <H> )
##
##  The subgroups of index two of H.  These are the preimages of the index
##  two subgroups of H/[H,H], so they are found in the abelianisation rather
##  than by a search through the subgroup lattice.
##
BindGlobal("HAPCHERN_IndexTwoSubgroups", function(H)
  local f, A, subs;
  if IsOddInt(Size(H)) then return []; fi;
  f := NaturalHomomorphismByNormalSubgroup(H, DerivedSubgroup(H));
  A := Image(f);
  subs := Filtered(MaximalSubgroupClassReps(A), S -> Index(A, S) = 2);
  return List(subs, S -> PreImage(f, S));
end);

#############################################################################
##
#F  MonomialModelOfCharacter( <H>, <chi> )
##
##  A monomial model of an irreducible character chi of H of degree two,
##  realised as Ind_K^H(lambda) for a subgroup K of index two and a linear
##  character lambda of K.  The trace of the model is checked against chi on
##  every element of H.
##
InstallGlobalFunction(MonomialModelOfCharacter, function(H, chi)
  local m, K, lam, g, matfun, model, h, tr;
  if Degree(chi) <> 2 then
    Error("MonomialModelOfCharacter: the character is not of degree two");
  fi;
  m := Exponent(H);
  for K in HAPCHERN_IndexTwoSubgroups(H) do
    for lam in Filtered(Irr(K), c -> Degree(c) = 1) do
      if InducedClassFunction(lam, H) = chi then
        g := First(Elements(H), x -> not x in K);
        matfun := function(h)
          if h in K then
            return [[HAPCHERN_CharacterValue(lam, h), 0],
                    [0, HAPCHERN_CharacterValue(lam, g^-1*h*g)]];
          else
            return [[0, HAPCHERN_CharacterValue(lam, h*g)],
                    [HAPCHERN_CharacterValue(lam, g^-1*h), 0]];
          fi;
        end;
        model := MonomialModelOfMatrixFunction(H, m, matfun);
        for h in model.elts do
          tr := TraceMat(model.mat[Position(model.elts, h)]);
          if tr <> HAPCHERN_CharacterValue(chi, h) then
            Error("MonomialModelOfCharacter: trace mismatch at ", h);
          fi;
        od;
        model.K := K; model.lambda := lam; model.g := g;
        return model;
      fi;
    od;
  od;
  Error("MonomialModelOfCharacter: chi is not induced from a subgroup of ",
        "index two, so it is not monomial in the sense used here");
end);

#############################################################################
##
#F  EulerCocycleOfMonomialRepresentation( <R>, <model> )
##
##  The 4-cocycle on the free generators of R_4 representing the Euler class,
##  up to the universal sign HAPCHERN_EulerSign.  R must be a free
##  ZH-resolution with R!.group = model.H, of length at least four.
##
##  The computation carries its own checks: the signed permutation matrices
##  form a left action, that action fixes the fundamental cycle (so the
##  representation is orientation preserving, as it must be for a complex
##  representation), each lift exists, and the image of d(e_{4,i}) is an
##  integer multiple of the fundamental cycle.
##
InstallGlobalFunction(EulerCocycleOfMonomialRepresentation, function(R, model)
  local H, elts, act, A, i, g, z, f, n, rhs, w, t, gi, x, coc, k, idx;

  H := model.H;
  if R!.group <> H then
    Error("EulerCocycleOfMonomialRepresentation: R is a resolution of a ",
          "different group");
  fi;
  if R!.dimension(4) = fail then
    Error("EulerCocycleOfMonomialRepresentation: R is too short; length at ",
          "least four is needed");
  fi;
  elts := model.elts;

  act := [];
  for i in [1..Length(R!.elts)] do
    g := R!.elts[i];
    gi := Position(elts, g);
    if gi = fail then
      Error("EulerCocycleOfMonomialRepresentation: R!.elts contains an ",
            "element outside the group");
    fi;
    act[i] := HAPCHERN_JoinAction(JoinComplexOfCyclicGroup(model.m),
                                  model.perm[gi], model.expo[gi]);
  od;

  ##  left action: (gh).v = g.(h.v) = (v*A_h)*A_g
  for g in GeneratorsOfGroup(H) do
    i := Position(R!.elts, g);
    for k in [1..Length(R!.elts)] do
      idx := Position(R!.elts, R!.elts[i]*R!.elts[k]);
      if idx = fail then continue; fi;
      if act[idx][4] <> act[k][4]*act[i][4] then
        Error("EulerCocycleOfMonomialRepresentation: not a left action");
      fi;
    od;
  od;

  z := HAPCHERN_FundamentalCycle(JoinComplexOfCyclicGroup(model.m));
  for i in [1..Length(R!.elts)] do
    if z * act[i][4] <> z then
      Error("EulerCocycleOfMonomialRepresentation: the action does not ",
            "preserve the orientation");
    fi;
  od;

  ##  f[n+1][i] = image in C_n of the i-th free generator of R_n
  f := [];
  if R!.dimension(0) <> 1 then
    Error("EulerCocycleOfMonomialRepresentation: R_0 is not of rank one");
  fi;
  f[1] := [ List(JoinComplexOfCyclicGroup(model.m).cells[1], c -> 0) ];
  f[1][1][1] := 1;
  for n in [1..3] do
    f[n+1] := [];
    for i in [1..R!.dimension(n)] do
      w := R!.boundary(n, i);
      rhs := List(JoinComplexOfCyclicGroup(model.m).cells[n], c -> 0);
      for t in w do
        rhs := rhs + SignInt(t[1]) * (f[n][AbsInt(t[1])] * act[t[2]][n]);
      od;
      x := SolutionIntMat(JoinComplexOfCyclicGroup(model.m).bnd[n+1], rhs);
      if x = fail then
        Error("EulerCocycleOfMonomialRepresentation: the lift failed in ",
              "degree ", n);
      fi;
      Add(f[n+1], x);
    od;
  od;

  coc := [];
  for i in [1..R!.dimension(4)] do
    w := R!.boundary(4, i);
    rhs := List(JoinComplexOfCyclicGroup(model.m).cells[4], c -> 0);
    for t in w do
      rhs := rhs + SignInt(t[1]) * (f[4][AbsInt(t[1])] * act[t[2]][4]);
    od;
    if not IsZero(rhs * JoinComplexOfCyclicGroup(model.m).bnd[4]) then
      Error("EulerCocycleOfMonomialRepresentation: the image is not a cycle");
    fi;
    k := rhs[JoinComplexOfCyclicGroup(model.m).pos[4]([[1,0],[1,0]])];
    if rhs <> k * z then
      Error("EulerCocycleOfMonomialRepresentation: the image is not a ",
            "multiple of the fundamental cycle");
    fi;
    Add(coc, k);
  od;
  return coc;
end);

#############################################################################
##
#F  HAPCHERN_CR4( <R> )
##
##  CR_CocyclesAndCoboundaries(R, 4, true), cached on the resolution.  R must
##  have length at least five for this.
##
BindGlobal("HAPCHERN_CR4", function(R)
  if not IsBound(R!.HAPCHERN_cr4) then
    R!.HAPCHERN_cr4 := CR_CocyclesAndCoboundaries(R, 4, true);
  fi;
  return R!.HAPCHERN_cr4;
end);

#############################################################################
##
#F  EulerClassOfMonomialRepresentation( <R>, <model> )
##
##  The Euler class as a vector of coordinates in H^4(H,Z), in the basis and
##  torsion coefficients of CR_CocyclesAndCoboundaries(R, 4, true).
##
InstallGlobalFunction(EulerClassOfMonomialRepresentation, function(R, model)
  local CR4, coc, tor, cls;
  CR4 := HAPCHERN_CR4(R);
  coc := EulerCocycleOfMonomialRepresentation(R, model);
  tor := CR4.torsionCoefficients;
  cls := CR4.cocycleToClass(coc);
  return List([1..Length(tor)],
              k -> (HAPCHERN_EulerSign * cls[k]) mod tor[k]);
end);

#############################################################################
##
#F  SecondChernClassOfCharacter( <R>, <chi> )
##
##  c_2 of the representation with character chi, as a vector of coordinates
##  in H^4(G,Z), where G = R!.group.  Every irreducible constituent of chi
##  must have degree at most two, and those of degree two must be monomial
##  (see MonomialModelOfCharacter); otherwise an error is raised.
##
##  The Whitney formula is applied to the decomposition of chi into
##  irreducible constituents with multiplicity:
##
##      c_2(W) = sum_j c_2(W_j) + sum_{j<k} c_1(W_j) c_1(W_k),
##
##  the summands of the first kind being Euler classes of the degree two
##  constituents, and those of the second kind Euler classes of the diagonal
##  representations det(W_j) (+) det(W_k).  Constituents of degree one
##  contribute nothing of the first kind.
##
InstallGlobalFunction(SecondChernClassOfCharacter, function(R, chi)
  local G, CR4, tor, cons, pieces, c, mult, i, j, ans, dets;

  G := R!.group;
  CR4 := HAPCHERN_CR4(R);
  tor := CR4.torsionCoefficients;

  pieces := [];
  for c in ConstituentsOfCharacter(chi) do
    if Degree(c) > 2 then
      Error("SecondChernClassOfCharacter: a constituent has degree ",
            Degree(c), "; this method covers degrees one and two only");
    fi;
    mult := ScalarProduct(chi, c);
    for i in [1..mult] do Add(pieces, c); od;
  od;

  ans := List(tor, t -> 0);
  for c in pieces do
    if Degree(c) = 2 then
      ans := ans + EulerClassOfMonomialRepresentation(R,
                       MonomialModelOfCharacter(G, c));
    fi;
  od;

  dets := List(pieces, DeterminantOfCharacter);
  for i in [1..Length(pieces)] do
    for j in [i+1..Length(pieces)] do
      ans := ans + EulerClassOfMonomialRepresentation(R,
                       MonomialModelOfLinearCharacterPair(G, dets[i], dets[j]));
    od;
  od;

  return List([1..Length(tor)], k -> ans[k] mod tor[k]);
end);
