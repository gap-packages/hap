LoadPackage( "hap" );
TestDirectory(DirectoriesPackageLibrary( "hap", "tst/testall" ),
  rec(exitGAP     := true,
      testOptions := rec(compareFunction := "uptowhitespace") ) );

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
