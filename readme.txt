The original readme.txt that went with the data on snps.psych.umn.edu
is below.  There also is a PDF document in this directory that was
prepared by Matt McGue.  The file phenotypes_2011-02-09.tab.txt was a
little different from the file McGue_Miller_20110209.tab.txt described
in the text below.  I changed some variable names:

FAMID became FamID
SEX became Sex

I added DadID, MomID and Zygo and I dropped SC.  I dropped the P3
(EEG) variable.

I also dropped all of the subjects who were not European American and
then I dropped the variables ETHNIC and WHITE because they were no
longer useful.

Matt's PDF document explains the variables.  I guess the file was once
called GEDI.DAT.


original readme.txt:

The data file that should be used for analysis is McGue_Miller_20110209.tab.txt.
The major phenotypes in this file should be the same as in the
20101028 data files created by Di Samek, but a few things have been
added, such as ethnicity and EIGENSTRAT principal component data.  The
file McGue_Miller_20110209_no_GWAS.tab.txt contains data for subjects
who have no GWAS data, but very few of them have usable data anyway
(e.g., no sex or age or anything else).

The file McGue_Documentation_20110202.pdf explains the phenotypes, but
the file was changed after it was written.  You should rely on it for
information about the variables, which haven't changed.  The file
described in the documentation had fixed column width and the current
form of the file is strictly tab-delimited.  The earlier data file
also included the subjects without GWAS data who have been moved to
the "no_GWAS" file, so the counts of numbers of subjects in the
documentation reflect the contents of the two current data files
combined.

Note that a directory called 20110203_McGue was deleted, when the data
were revised, but only Saonli and her group saw those files.  They
should delete their copies and work with these.  Specifically, the
file McGue_update_to_20101028.tab.txt should not be used.

Added 2/15/2011: I deleted missing value codes to make a file where
all missing data are simply missing (consecutive tabs) instead of
using codes.  That new file is here:

McGue_Miller_20110209_missing_empty.tab.txt

I also loaded the data into R and created the .Rhistory and .RData
files.

--Mike Miller (user mbmiller)
email: mbmiller+m@gmail.com
