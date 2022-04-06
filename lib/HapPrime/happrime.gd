#############################################################################
##
##  HAPPRIME - happrime.gd
##  General Functions, Operations and Methods
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

#####################################################################
##  <#GAPDoc Label="InfoHAPprime_manMisc">
##  <ManSection>
##    <InfoClass Name="InfoHAPprime"/>
##    <Description>
##  The <C>InfoHAPprime</C> info class is used throughout the &HAPprime; package.
##  Use <C>SetInfoLevel(InfoHAPprime, </C><A>level</A><C>)</C> to change the
##  amount of information displayed about the progress of the computation
##  (see <Ref Oper="SetInfoLevel" BookName="ref"/> in the &GAP; reference 
##  manual). The different distinct levels are:
##  <List>
##    <Item><C>0</C> print nothing (this is the default)</Item>
##    <Item><C>1</C> print some information, mainly progress information 
##    during computations that may take some time</Item>
##    <Item><C>2</C> print more detailed information, incluing details of 
##    internal calculations</Item>
##  </List>
##    </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
DeclareInfoClass("InfoHAPprime");

DeclareGlobalFunction("MakeHAPprimeDoc");
DeclareGlobalFunction("HAPPRIME_VersionWithSVN");
