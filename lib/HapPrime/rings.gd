#############################################################################
##
##  HAPPRIME - rings.gd
##  Functions, Operations and Methods for dealing with cohomology rings
##  Paul Smith
##
##  Copyright (C)  2007-2008
##  Paul Smith
##  National University of Ireland Galway
##
##  This file is part of HAPprime. 
## 
##  HAPprime is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
## 
##  HAPprime is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
## 
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <https://www.gnu.org/licenses/>.
## 
##
#############################################################################

DeclareAttribute("PresentationOfGradedStructureConstantAlgebra", IsAlgebra);

DeclareOperation("Mod2CohomologyRingPresentation", [IsGroup]);
DeclareOperation("Mod2CohomologyRingPresentation", [IsGroup, IsPosInt]);
DeclareOperation("Mod2CohomologyRingPresentation", [IsHapResolution]);
DeclareOperation("Mod2CohomologyRingPresentation", [IsAlgebra]);

DeclareOperation("LHSSpectralSequenceLastSheet", [IsGroup, IsGroup]);
DeclareOperation("LHSSpectralSequence", [IsGroup, IsGroup, IsPosInt]);

DeclareAttribute("ModPRingGeneratorDegrees", IsAlgebra);
DeclareAttribute("ModPRingNiceBasis", IsAlgebra);
DeclareAttribute("ModPRingNiceBasisAsPolynomials", IsAlgebra);
DeclareAttribute("ModPRingBasisAsPolynomials", IsAlgebra);

DeclareGlobalFunction("HAPPRIME_LHSSpectralSequence");
DeclareGlobalFunction("HAPPRIME_CohomologyRingWithoutResolution");
DeclareGlobalFunction("HAPPRIME_Polynomial2Algebra");
DeclareGlobalFunction("HAPPRIME_Algebra2Polynomial");
