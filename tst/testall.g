LoadPackage( "hap" );

TestDirectory(DirectoriesPackageLibrary( "hap", "tst/testextra" ),
  rec(exitGAP     := false,
      testOptions := rec(compareFunction := "uptowhitespace") ) );

TestDirectory(DirectoriesPackageLibrary( "hap", "tst/testall" ),
  rec(exitGAP     := true,
      testOptions := rec(compareFunction := "uptowhitespace") ) );

if not GAPInfo.Architecture="x86_64-pc-linux-gnu-default32-kv8" then
TestDirectory(DirectoriesPackageLibrary( "hap", "tst/testall2" ),
  rec(exitGAP     := true,
      testOptions := rec(compareFunction := "uptowhitespace") ) );
fi;



FORCE_QUIT_GAP(1); # if we ever get here, there was an error
