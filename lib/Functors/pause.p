# Source - https://stackoverflow.com/a
# Posted by DVK
# Retrieved 2025-11-22, License - CC BY-SA 3.0

#print "Enter a word to look up: ";
my $userword = <STDIN>; # I moved chomp to a new line to make it more readable
chomp $userword; # Get rid of newline character at the end
exit 0 if ($userword eq ""); # If empty string, exit.

