DeclareFilter("QuadraticNF");
DeclareProperty("IsQuadraticNumberField",QuadraticNF);
DeclareFilter("RingOfQuadraticIntegers");
DeclareProperty("IsRingOfQuadraticIntegers",RingOfQuadraticIntegers);
DeclareAttribute("AssociatedNumberField",IsNumberField);
DeclareAttribute("AssociatedRing",IsRing);
DeclareAttribute("NormOfPrincipalIdealGenerator",IsInt);
DeclareGlobalFunction("QuadraticNumberField");
DeclareOperation("RingOfIntegers",[IsNumberField]);
DeclareOperation("PrincipalIdeal",[IsRing and IsRingOfQuadraticIntegers, IsCyclotomic]);
DeclareFilter("PrincipalIdealOfRing");
DeclareProperty("IsPrincipalIdeal",PrincipalIdealOfRing);




