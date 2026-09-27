
########################################################
########################################################
InstallMethod(FiniteProjectiveLine,
"Finite projective line for the ring Z/nZ",
[IsInt],
function(n)
return HAP_FiniteProjectiveLineIntegers(n);
end);
########################################################
########################################################

########################################################
########################################################
InstallMethod(FiniteProjectivePlane,
"Finite projective plane for the ring Z/nZ",
[IsInt],
function(n)
return HAP_FiniteProjectivePlaneIntegers(n);
end);
########################################################
########################################################
InstallMethod(FiniteProjectiveLine_alt,
"Finite projective line for the ring Z/nZ",
[IsInt],
function(n)
return HAP_FiniteProjectiveLineIntegers_alt(n);
end);

InstallMethod(FiniteProjectivePlane_alt,
"Finite projective line for the ring Z/nZ",
[IsInt],
function(n)
return HAP_FiniteProjectivePlaneIntegers_alt(n);
end);

InstallMethod(FiniteProjectiveSpace3_alt,
"3 dimensional finite projective space for the ring Z/nZ",
[IsInt],
function(n)
return HAP_FiniteProjectiveSpace3Integers_alt(n);
end);

InstallMethod(FiniteProjectiveSpace3,
"3 dimensional finite projective space for the ring Z/nZ",
[IsInt],
function(n)
return HAP_FiniteProjectiveSpace3Integers(n);
end);


InstallGlobalFunction(HAP_FiniteProjectiveLineIntegers_alt,
function(n)
    local UnitEls, x, y, i, c, d, u, UnitsAction, Representatives, 
          RepOf, r, v, w, m, min;

    UnitEls := Units(Integers mod n);

    UnitsAction := function(c, u)
      local uu;
      uu :=Int (u);
      return List(c, x -> (uu * x) mod n);
    end;

    Representatives := [];

    RepOf := [];
    for x in [1..n] do
      RepOf[x] := [];
    od;

    for x in [0..n-1] do
        for y in [0..n-1] do
          if Gcd(x,y,n) = 1 then
            v := [x,y];
            if not IsBound(RepOf[x+1][y+1]) then
              m := Orbit(UnitEls,v,UnitsAction);
              min := Minimum(m);
              AddSet(Representatives,min);
              for w in m + 1 do
                RepOf[w[1]][w[2]] := min;
              od;
            fi;
          fi;
        od;
    od;

    return rec(
        Reps := Set(Representatives),
        RepOf:= RepOf
    );
end);

InstallGlobalFunction(HAP_FiniteProjectivePlaneIntegers_alt,
function(n)
    local UnitEls, x, y, z, i, c, d, u, UnitsAction, Representatives, 
          RepOf, r, v, w, m, min;

    UnitEls := Units(Integers mod n);

    UnitsAction := function(c, u)
      local uu;
      uu :=Int (u);
      return List(c, x -> (uu * x) mod n);
    end;

    Representatives := [];

    RepOf := [];
    for x in [1..n] do
      RepOf[x] := [];
      for y in [1..n] do
        RepOf[x][y] := [];
      od;
    od;

    for x in [0..n-1] do
      for y in [0..n-1] do
        for z in [0..n-1] do
          if Gcd(x,y,z,n) = 1 then
            v := [x,y,z];
            if not IsBound(RepOf[x+1][y+1][z+1]) then
              m := Orbit(UnitEls,v,UnitsAction);
              min := Minimum(m);
              AddSet(Representatives,min);
              for w in m + 1 do
                RepOf[w[1]][w[2]][w[3]] := min;
              od;
            fi;
          fi;
        od;
      od;
    od;

    return rec(
        Reps := Set(Representatives),
        RepOf:= RepOf
    );
end);

InstallGlobalFunction(HAP_FiniteProjectiveLineIntegers,
function(n)
    local Rep, i, factors, Znd, d, z, p, l, toFill, t;

    Rep := [[0,1]];

    for i in [1..n] do
        Add(Rep,[1,i-1]);
    od;
    
    factors := List(DivisorsInt(n));
    Remove(factors,1);
    Remove(factors);

    for p in factors do
        d := n/p;

        Znd := [0..d-1];
        toFill := [];

        for z in Znd do
            if Gcd(p,z,d) = 1 then
                Add(toFill, z);
            fi;
        od;

        t := 1;
        while not IsEmpty(toFill) do
            if Gcd(p,t) = 1 then
                if (t mod d) in toFill then
                    Add(Rep,[p,t]);
                    Remove(toFill, Position(toFill,t mod d));
                fi;
            fi;
            t := t+1;
        od;
    od;

    return Rep;
end);

InstallGlobalFunction(HAP_FiniteProjectivePlaneIntegers,
function(n)
    local divs, Rep, d, q, g, q_g, m, u, U_m, dd, qq, z, uu, zz;

    divs := List(DivisorsInt(n));
    Remove(divs);
    divs := Concatenation([n],divs);
    Rep := [];

    for d in divs do
      q := n/d;
      for g in divs do
        q_g := n/g;
        m := Gcd(q,q_g);
        U_m := Filtered([0..m-1], z -> Gcd(z,m) = 1);

        for u in U_m do
          uu := u;          
          while Gcd(uu, q_g) <> 1 do
            uu := uu + m;
          od;

          dd := Gcd(d,g);
          qq := n/dd;

          for z in [0..qq-1] do
            if Gcd(z,dd,qq) = 1 then
              zz := z;
              while Gcd(zz, dd) <> 1 do
                zz := zz + qq;
              od;

              Add(Rep, [d mod n,g*uu mod n,zz]);
            fi;
          od;
        od;
      od;
    od;

    return rec(Reps := Set(List(Rep, Immutable)));;
end);

InstallGlobalFunction(HAP_FiniteProjectiveSpace3Integers,
function(n)
    local divs, Rep, d_x, q_x, d_y, q_y, q_xy, d_z, q_z, Q, y_0, y_00, Y_0, z_0, z_00, Z_0, QQ, t, tt;

    divs := List(DivisorsInt(n));
    Remove(divs);
    divs := Concatenation([n],divs);
    Rep := [];

    for d_x in divs do
      q_x := n/d_x;
      for d_y in divs do
        q_y := n/d_y;
        q_xy := Gcd(q_x,q_y);
        Y_0 := Filtered([0..q_xy-1], e -> Gcd(e,q_xy) = 1);
        Q := Lcm(q_x,q_y);

        for y_0 in Y_0 do
          y_00 := y_0;          
          while Gcd(y_00, q_y) <> 1 do
            y_00 := y_00 + q_xy;
          od;

          for d_z in divs do
            q_z := n/d_z;
            QQ := Lcm(Q,q_z);

            Z_0 := Filtered([0..Gcd(Q,q_z)-1], e -> Gcd(e,Gcd(Q,q_z)) = 1);
            for z_0 in Z_0 do
              z_00 := z_0;          
              while Gcd(z_00, q_z) <> 1 do
                z_00 := z_00 + Gcd(Q,q_z);
              od;
              for t in [0..QQ-1] do
                if Gcd(t,QQ,n/QQ) = 1 then
                  tt := t;
                  while Gcd(tt, n/QQ) <> 1 do
                    tt := tt + QQ;
                  od;
                  Add(Rep, [d_x mod n,d_y*y_00 mod n,d_z*z_00 mod n, tt]);
                fi;
              od;
            od;
          od;
        od;
      od;
    od;

    return rec(Reps := Set(List(Rep, Immutable)));;
end);


InstallGlobalFunction(HAP_FiniteProjectiveSpace3Integers_alt,
function(n)
    local UnitEls, x, y, z, t, i, c, d, u, UnitsAction, Representatives, 
          RepOf, r, v, w, m, min;

    UnitEls := Units(Integers mod n);

    UnitsAction := function(c, u)
      local uu;
      uu :=Int (u);
      return List(c, x -> (uu * x) mod n);
    end;

    Representatives := [];

    RepOf := [];
    for x in [1..n] do
      RepOf[x] := [];
      for y in [1..n] do
        RepOf[x][y] := [];
        for z in [1..n] do
          RepOf[x][y][z] := [];
        od;
      od;
    od;

    for x in [0..n-1] do
      for y in [0..n-1] do
        for z in [0..n-1] do
          for t in [0..n-1] do
            if Gcd(x,y,z,t,n) = 1 then
              v := [x,y,z,t];
              if not IsBound(RepOf[x+1][y+1][z+1][t+1]) then
                m := Orbit(UnitEls,v,UnitsAction);
                min := Minimum(m);
                AddSet(Representatives,min);
                for w in m + 1 do
                  RepOf[w[1]][w[2]][w[3]][w[4]] := min;
                od;
              fi;
            fi;
          od;
        od;
      od;
    od;

    return rec(
        Reps := Set(Representatives),
        RepOf:= RepOf
    );
end);