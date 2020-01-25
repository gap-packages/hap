#############################################################################
##
##  PackageInfo.g                HAP Package                Graham Ellis 
##
#############################################################################

SetPackageInfo( rec(

  PackageName := "HAP",
  Subtitle  := "Homological Algebra Programming",
  Version := "1.25",
  Date    := "25/01/2020",

  SourceRepository := rec(
      Type := "git",
      URL := Concatenation( "https://github.com/gap-packages/", LowercaseString(~.PackageName) ),
  ),
  IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
  PackageWWWHome  := Concatenation( "https://gap-packages.github.io/", LowercaseString(~.PackageName) ),
  README_URL      := Concatenation( ~.PackageWWWHome, "/README.md" ),
  PackageInfoURL  := Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),
  ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                   "/releases/download/v", ~.Version,
                                   "/", LowercaseString(~.PackageName), "-", ~.Version ),
  ArchiveFormats := ".tar.gz",

  Persons := [ 
    rec( 
      LastName      := "Ellis",
      FirstNames    := "Graham",
      IsAuthor      := true,
      IsMaintainer  := true,
      Email         := "graham.ellis@nuigalway.ie",
      WWWHome       := "http://hamilton.nuigalway.ie",
      PostalAddress := Concatenation( [
                         "Graham Ellis\n",
                         "Mathematics Department\n",
                         "NUI Galway\n",
                         "Galway\n",
                         "Ireland" ] ),
      Place         := "Galway",
      Institution   := "National University of Ireland, Galway"
    )
  ],  

  Status  := "accepted",
  CommunicatedBy 
          := "Derek Holt (Warwick)",
  AcceptDate 
          := "03/2006",

  AbstractHTML := 
    "This package provides some functions for group cohomology. ",
                  
  PackageDoc := rec(
    BookName  := "HAP",
    ArchiveURLSubset := ["doc", "www"],
    HTMLStart := "www/index.html",
    PDFFile   := "doc/manual.pdf",
    SixFile   := "doc/manual.six",
    LongTitle := "Homological Algebra Programming Package",
    Autoload := true 
  ),


  Dependencies := rec(
    GAP := ">= 4.5.6",
    NeededOtherPackages := [
                             [ "polycyclic", ">=1.1" ],
                             [ "crystcat",   ">=1.1" ],
                             [ "fga",        ">=1.1" ],
                             [ "aclib",      ">=1.1" ],
                             [ "nq",         ">=1.1" ],
                            #[ "gapdoc",     ">=0.0" ]
                           ],
    SuggestedOtherPackages := [
                              #[ "polycyclic", ">=1.1" ],
                              #[ "aclib",      ">=1.1" ],
                               [ "gapdoc",      ">=0.0" ],
			      #[ "nq",         ">=1.1" ],
                               [ "nql",         ">=0.0" ],
                               [ "homology",    ">=0.0"   ], 
			       [ "edim",      ">=1.2.2" ],
			       [ "singular", ">=06.07.23" ],
                               [ "congruence", ">=0,0" ],
                               [ "HAPcryst", ">0.0" ],
                               [ "xmod", ">0.0" ],
			      ],

    ExternalConditions := [["Some optional functions require Polymake software",
    "http://polymake.org/doku.php"],
    ["Some optional functions require Graphviz software",
    "http://www.graphviz.org/"],
     ["One optional function requires the Simplicial Homology GAP package",
         "http://www.cis.udel.edu/~dumas"]
    ]
  ),

AvailabilityTest := ReturnTrue,

BannerString     := Concatenation( "Loading HAP ",
                            String( ~.Version ), " ...\n" ),

TestFile := "tst/testall.g",

Keywords := [ "homology", "cohomology", "resolution", "homotopy group", 
"module of identities", "CW complex", "simplicial complex", "cubical complex", "permutahedral complex", "knots", "nonabelian tensor", "nonabelian exterior", "covering space" ]

));

