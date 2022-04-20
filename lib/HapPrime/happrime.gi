#############################################################################
##
##  HAPPRIME - happrime.gi
##  General Functions, Operations and Methods
##  Paul Smith
##
##  Copyright (C) 2007-2008
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
##  <#GAPDoc Label="MakeHAPprimeDoc_manMisc">
##  <ManSection>
##  <Func Name="MakeHAPprimeDoc" Arg="[manual-name]"/>
##
##  <Returns>
##    nothing
##  </Returns>
##  <Description>
##  The two manuals supplied with &HAPprime; - this user guide and the datatypes
##  reference manual - are written using the &GAPDoc; package and are
##  available in PDF, HTML and text format. It should not be necessary to rebuild 
##  these files, but should you wish to do so then this can be done using the 
##  function  <Ref Func="MakeHAPprimeDoc"/>.
##  <P/>
##  The optional argument <A>manual-name</A> is a string specifying which 
##  manuals to build. It may be one of the following
##  <List>
##    <Item><C>"all"</C> builds both manuals. This is the default</Item> 
##    <Item><C>"userguide"</C> builds just the user guide</Item> 
##    <Item><C>"datatypes"</C> builds just the datatypes reference manual</Item> 
##    <Item><C>"internal"</C> builds both manuals, including the otherwise 
##      undocumented internal functions</Item> 
##    <Item><C>"testexamples"</C> builds neither manual, but tests all of the
##      examples using <Ref Func="TestManualExamples" BookName="GAPDoc"/></Item>
##  </List>
##  <P/>
##  As well as building the manuals, this function at the same time builds &GAP; 
##  test files <Ref Sect="Test Files" BookName="ref"/> which means that 
##  all of the testable examples in the manuals are added to the &HAPprime; test 
##  routines described in Section <Ref Sect="LoadingAndTesting"/>.
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallGlobalFunction(MakeHAPprimeDoc,
  function(arg)
    local narg, build, makeinternal, docdir, userguidedir, datatypesdir, 
    sourcefiles, preprocessfiles, sysdir, tstdir, examples, stream, 
    testexamples, infolevel, datatypesGAP4stones, userguideGAP4stones;

    userguideGAP4stones := 4550000000;
    datatypesGAP4stones := 2000000000;
    
    # Parse the input
    if Length(arg) = 0 then
      arg := ["all"];
    fi;

    if Length(arg) > 1  then
      Error("MakeHAPprimeDoc must have zero or one arguments.");
    fi;

    if not IsString(arg[1]) then
      Error("the argument to MakeHAPprimeDoc must be a string.");
    fi;

    makeinternal := false;
    testexamples := false;
    if LowercaseString(arg[1]) = "all" then
      build := [true, true];
    elif LowercaseString(arg[1]) = "userguide" then
      build := [true, false];
    elif LowercaseString(arg[1]) = "datatypes" then
      build := [false, true];
    elif LowercaseString(arg[1]) = "internal" then
      build := [true, true];
      makeinternal := true;
    elif LowercaseString(arg[1]) = "testexamples" then
      build := [true, true];
      makeinternal := true;
      testexamples := true;
    else
      Error("the argument to MakeHAPprimeDoc must one of \"all\", \"userguide\", \"datatypes\", \"internal\" or \"testexamples\".");
    fi;

    # First merge the include files
    sysdir := DirectoriesSystemPrograms();
    docdir := DirectoriesPackageLibrary("HAPprime", "doc")[1];
    tstdir := DirectoriesPackageLibrary("HAPprime", "tst")[1];
    userguidedir := DirectoriesPackageLibrary("HAPprime", "doc/userguide")[1];
    datatypesdir := DirectoriesPackageLibrary("HAPprime", "doc/datatypes")[1];
    sourcefiles := [ 
      "../../lib/happrime.gd", 
      "../../lib/happrime.gi", 
      "../../lib/basemat.gd",
      "../../lib/basemat.gi",
      "../../lib/fpgmoduleint.gd", 
      "../../lib/fpgmoduleint.gi", 
      "../../lib/fpgmoduleG.gd", 
      "../../lib/fpgmoduleG.gi", 
      "../../lib/fpgmodulehomG.gd", 
      "../../lib/fpgmodulehomG.gi", 
      "../../lib/resolutions.gd", 
      "../../lib/resolutions.gi", 
      "../../lib/poincare.gd", 
      "../../lib/poincare.gi", 
      "../../lib/gradedalgebra.gd", 
      "../../lib/gradedalgebra.gi", 
      "../../lib/polynomials.gd", 
      "../../lib/polynomials.gi", 
      "../../lib/derivation.gd", 
      "../../lib/derivation.gi", 
      "../../lib/ringhomomorphism.gd",
      "../../lib/ringhomomorphism.gi",
      "../../lib/rings.gd",
      "../../lib/rings.gi",
      "../../lib/singular.gd",
      "../../lib/singular.gi",
      "../../lib/groups.gd",
      "../../lib/groups.gi",
      "../../lib/test.gd",
      "../../lib/test.gi"];

    # Now build the requested manuals 
    if build[1] then
      preprocessfiles := [
        "intro", 
        "functions"];
      Perform(preprocessfiles, 
        i->Process(userguidedir, Filename(sysdir, "sh"),
          InputTextUser(), OutputTextUser(), 
          ["-c", Concatenation("../includesourcedoc.sh ", i, ".xml.in ", i, ".xml")]));
      if not testexamples then
        MakeGAPDocDoc(userguidedir, "userguide.xml", sourcefiles, 
          "HAPprime", "../../../../../gap4r4");;
        examples := ManualExamples(userguidedir, "userguide.xml", sourcefiles, "Single"); 
        stream := OutputTextFile(Filename(tstdir, "userguideexamples.tst"), false);
        SetPrintFormattingStatus(stream, false);
        if not IsEmpty(Flat(examples)) then
          PrintTo(stream, "gap> START_TEST(\"HAPprime version ", 
            InstalledPackageVersion("HAPprime"), " userguide examples\");\n");
          PrintTo(stream, "gap> HAPprimeInfoLevel := InfoLevel(InfoHAPprime);;\n");
          PrintTo(stream, "gap> SetInfoLevel(InfoHAPprime, 0);;\n");
          PrintTo(stream, "gap> #\n");
          PrintTo(stream, Flat(examples));
          PrintTo(stream, "gap> #\n");
          PrintTo(stream, "gap> SetInfoLevel(InfoHAPprime, HAPprimeInfoLevel);\n");
          PrintTo(stream, 
            "gap> STOP_TEST(\"HAPprime/tst/userguideexamples.tst\", ",
              String(userguideGAP4stones), ");\n");
        fi;
        CloseStream(stream);
      else
        infolevel := InfoLevel(InfoHAPprime);;
        SetInfoLevel(InfoHAPprime, 0);
        TestManualExamples(userguidedir, "userguide.xml", sourcefiles);  
        SetInfoLevel(InfoHAPprime, infolevel);
      fi;
    fi;
    
    if build[2] then
      preprocessfiles := [
        "fpgmodule", 
        "fpgmodulehom", 
        "resolution", 
        "poincare", 
        "derivation", 
        "ringhomomorphism", 
        "gradedalgebrapresentation",
        "gradedalgebra",
        "GAPfunctions"];
      Perform(preprocessfiles, 
        i->Process(datatypesdir, Filename(sysdir, "sh"),
          InputTextUser(), OutputTextUser(), 
          ["-c", Concatenation("../includesourcedoc.sh ", i, ".xml.in ", i, ".xml")]));
      if makeinternal or testexamples then
        Process(datatypesdir, Filename(sysdir, "sh"), 
        InputTextUser(), OutputTextUser(), 
        [ "-c", "../includesourcedoc.sh internal.xml.in internal.xml"]);
        Process(datatypesdir, Filename(sysdir, "cp"), 
        InputTextUser(), OutputTextUser(), 
        ["introinternal.xml.in", "introinternal.xml"]);
      else
        Process(datatypesdir, Filename(sysdir, "cp"), 
        InputTextUser(), OutputTextUser(), 
        ["internal.xml.none", "internal.xml"]);
        Process(datatypesdir, Filename(sysdir, "cp"), 
        InputTextUser(), OutputTextUser(), 
        ["introinternal.xml.none", "introinternal.xml"]);
      fi;
      if not testexamples then
        # Now make the gapdoc
        MakeGAPDocDoc(datatypesdir, "datatypes.xml", sourcefiles, 
          "HAPprime Datatypes", "../../../../../gap4r4");;
        # Now build the examples
        examples := ManualExamples(datatypesdir, "datatypes.xml", sourcefiles, "Single"); 
        stream := OutputTextFile(Filename(tstdir, "datatypesexamples.tst"), false);
        SetPrintFormattingStatus(stream, false);
        if not IsEmpty(Flat(examples)) then
          PrintTo(stream, "gap> START_TEST(\"HAPprime version ", 
            InstalledPackageVersion("HAPprime"), " datatypes reference manual examples\");\n");
          PrintTo(stream, "gap> HAPprimeInfoLevel := InfoLevel(InfoHAPprime);;\n");
          PrintTo(stream, "gap> SetInfoLevel(InfoHAPprime, 0);;\n");
          PrintTo(stream, "gap> #\n");
          PrintTo(stream, Flat(examples));
          PrintTo(stream, "gap> #\n");
          PrintTo(stream, "gap> SetInfoLevel(InfoHAPprime, HAPprimeInfoLevel);\n");
          PrintTo(stream, 
            "gap> STOP_TEST(\"HAPprime/tst/datatypesexamples.tst\", ",
              String(datatypesGAP4stones), ");\n");
        fi;
        CloseStream(stream);
      else
        infolevel := InfoLevel(InfoHAPprime);;
        SetInfoLevel(InfoHAPprime, 0);
        TestManualExamples(datatypesdir, "datatypes.xml", sourcefiles); 
        SetInfoLevel(InfoHAPprime, infolevel);
      fi;
    fi;

  end);

#####################################################################
##  <#GAPDoc Label="HAPPRIME_VersionWithSVN_manTestInt">
##  <ManSection>
##  <Func Name="HAPPRIME_VersionWithSVN" Arg=""/>
##
##  <Returns>
##    String
##  </Returns>
##  <Description>
##  Returns a string giving a current version number for the HAPprime 
##  installation assuming that it is checked out from a subversion repository.
##  This fetches the current version number from <C>PackageInfo.g</C> and
##  appends the return of the <C>svnversion</C> program to this, returning
##  the resulting composite string.
##  <Log><![CDATA[
##  gap> HAPPRIME_VersionWithSVN();
##  "0.2.1.302:319M"
##  ]]></Log>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
#####################################################################
InstallGlobalFunction(HAPPRIME_VersionWithSVN,
  function()
  local svnVersionString, version;
  version := InstalledPackageVersion("HAPprime");
  svnVersionString := "";
  Process(
    DirectoriesPackageLibrary("HAPprime")[1], 
    Filename(DirectoriesSystemPrograms(), "svnversion"),
    InputTextUser(), 
    OutputTextString(svnVersionString, false), 
    ["-n"]);
  return Concatenation(version, ".", svnVersionString);
end);
#####################################################################
  





## To do for a new release:
##  - update .Version, .Date and .ArchiveURL in PackageInfo.g
##  - update VERSION file to match
##  - update date in userguide.xml and datatypes.xml file to match
##  - update CHANGES file
##  - check the README file
##  - update the version number in tst/happrime.tst 
##  - check test ReadTest("/home/pas/GAP/pkg/HAPPrime/tst/happrime.tst");
##  - tweak GAP4stones number in tst/happrime.tst to be between 300000-400000
##  - check the manual examples using MakeHAPprimeDoc("testexamples");
##  - check test routine ReadPackage("HAPprime", "tst/testall.g");
##  - tweak GAP4stones numbers for examples in lib/happrime.gi
##  - make manual MakeHAPprimeDoc();
##  - run ./make_tarball
##  - update the website with updatewebman.sh and updatewebpkg.sh
##  - update download.html on local website 
##  - update website by running cd ~/public_html/CHA/HAPprime/;./copytolarmor.sh
##  - tag the release on subversion

