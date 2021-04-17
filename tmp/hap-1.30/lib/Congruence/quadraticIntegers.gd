DeclareFilter("QuadraticNF");
DeclareProperty("IsQuadraticNumberField",QuadraticNF);
DeclareFilter("RingOfQuadraticIntegers");
DeclareProperty("IsRingOfQuadraticIntegers",RingOfQuadraticIntegers);
DeclareAttribute("AssociatedNumberField",IsNumberField);
DeclareAttribute("AssociatedRing",IsRing);
DeclareAttribute("NormOfIdeal",IsInt);
DeclareGlobalFunction("QuadraticNumberField");
DeclareOperation("RingOfIntegers",[IsNumberField]);
DeclareOperation("QuadraticIdeal",[IsRing and IsRingOfQuadraticIntegers, IsCyclotomic]);
DeclareFilter("IdealOfQuadraticIntegers");
DeclareProperty("IsIdealOfQuadraticIntegers",IdealOfQuadraticIntegers);
DeclareGlobalFunction("PartsOfQuadraticInteger");
DeclareGlobalFunction("SL2QuadraticIntegers");



