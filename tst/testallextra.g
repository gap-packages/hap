LoadPackage( "hap" );

dirs := Concatenation(
  DirectoriesPackageLibrary( "hap", "tst/testextra" ),
  DirectoriesPackageLibrary( "hap", "tst/testall" ),
  DirectoriesPackageLibrary( "hap", "tst/testall2" ),
  DirectoriesPackageLibrary( "hap", "tst/testextraextra" )
);

TestDirectory(dirs,
  rec(exitGAP     := true,
      testOptions := rec(compareFunction := "uptowhitespace") ) );

FORCE_QUIT_GAP(1); # if we ever get here, there was an error
