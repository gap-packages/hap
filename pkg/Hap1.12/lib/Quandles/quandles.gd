##
##
##			Cédric FRAGNAUD - 2017
##
##
##

########################################## QUANDLES ##########################################

DeclareCategory("IsHapConjQuandElt",IsMultiplicativeElement);
DeclareRepresentation(  "IsHapConjQuandEltRep", IsComponentObjectRep, IsMultiplicativeElementWithInverse);
HapConjQuandEltFamily:=NewFamily( "HapConjQuandEltFamily", IsHapConjQuandElt, IsHapConjQuandElt);
HapConjQuandElt:=NewType(HapConjQuandEltFamily, IsHapConjQuandEltRep);

########

DeclareGlobalFunction("Cedric_ConjugateQuandleElement");
DeclareGlobalFunction("ConjugationQuandle");
DeclareGlobalFunction("FirstQuandleAxiomIsSatisfied");
DeclareGlobalFunction("SecondQuandleAxiomIsSatisfied");
DeclareGlobalFunction("ThirdQuandleAxiomIsSatisfied");
DeclareGlobalFunction("IsQuandle");
DeclareGlobalFunction("Cedric_CheckThirdAxiomRow");
DeclareGlobalFunction("Cedric_Permute");
DeclareGlobalFunction("Cedric_Quandle1");
DeclareGlobalFunction("Cedric_Quandle2");
DeclareGlobalFunction("Cedric_Quandle3");
DeclareGlobalFunction("Cedric_Quandle4");
DeclareGlobalFunction("Cedric_Quandle5");
DeclareGlobalFunction("Cedric_Quandle6");
Quandles:=NewOperation("Quandles",[IsInt]);
Quandle:=NewOperation("Quandle",[IsInt,IsInt]);
DeclareGlobalFunction("IdQuandle");
IsLatin:=NewOperation("IsLatin",[IsMagma]);
IsConnected:=NewOperation("IsConnected",[IsMagma]);
ConnectedQuandles:=NewOperation("ConnectedQuandles",[IsInt]);
ConnectedQuandle:=NewOperation("ConnectedQuandle",[IsInt,IsInt]);
DeclareGlobalFunction("IdConnectedQuandle");
DeclareGlobalFunction("IsQuandleEnvelope");
DeclareGlobalFunction("QuandleQuandleEnvelope");
RightMultiplicationGroupOfQuandleAsPerm:=NewOperation("RightMultiplicationGroupOfQuandleAsPerm",[IsMagma]);
RightMultiplicationGroupOfQuandle:=NewOperation("RightMultiplicationGroupOfQuandle",[IsMagma]);
DeclareGlobalFunction("Cedric_FromAutGeReToAutQe");
AutomorphismGroupQuandleAsPerm:=NewOperation("AutomorphismGroupQuandleAsPerm",[IsMagma]);
AutomorphismGroupQuandle:=NewOperation("AutomorphismGroupQuandle",[IsMagma]);

DeclareGlobalFunction("TupleOrbitReps_perm");
DeclareGlobalFunction("TupleOrbitReps");
DeclareGlobalFunction("NumberOfHomomorphisms_connected");
#DeclareGlobalFunction("FundamentalGroupFromQuandle");
DeclareGlobalFunction("AdjointGroupOfQuandle");
DeclareGlobalFunction("DerivedGroupOfQuandle");


########################################## KNOTS ##########################################

DeclareGlobalFunction("PresentationKnotQuandle");
DeclareGlobalFunction("PD2GC");
DeclareGlobalFunction("PlanarDiagramKnot");
DeclareGlobalFunction("GaussCodeKnot");
DeclareGlobalFunction("PresentationKnotQuandleKnot");
DeclareGlobalFunction("Cedric_IsHomomorphism");
NumberOfHomomorphisms:=NewOperation("NumberOfHomomorphisms",[IsRecord,IsMagma]);
PartitionedNumberOfHomomorphisms:=NewOperation("PartitionedNumberOfHomomorphisms",[IsRecord,IsMagma]);
DeclareGlobalFunction("KnotInvariantCedric");
