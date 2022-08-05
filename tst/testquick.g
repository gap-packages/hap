LoadPackage( "hap" );

dirs := [
    DirectoriesPackageLibrary( "hap", "tst/testall1" ),
    DirectoriesPackageLibrary( "hap", "tst/testall2" ),
    DirectoriesPackageLibrary( "hap", "tst/testextra" ),
    DirectoriesPackageLibrary( "hap", "tst/testextra2" ),
];

#if StartsWith(GAPInfo.Architecture,"x86_64-pc-linux-gnu-default64") then
if GAPInfo.BytesPerVariable >= 8 then
  Add(dirs, DirectoriesPackageLibrary( "hap", "tst/testall3" ));

# if (GAPInfo.Version{[1..4]}="4.11" or GAPInfo.Version{[1..4]}="4.12") then
  if CompareVersionNumbers(ReplacedString(GAPInfo.Version, "dev", ""), "4.11") then
    Add(dirs, DirectoriesPackageLibrary( "hap", "tst/testallV11" ));
  fi;
fi;

TestDirectory(dirs,
rec(exitGAP     := true,
    testOptions := rec(compareFunction := "uptowhitespace") ) );

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
