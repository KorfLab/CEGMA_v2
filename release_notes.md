## CEGMA RELEASE NOTES ##

## [v2.5](http://korflab.ucdavis.edu/Datasets/cegma/cegma_v2.5.tar.gz) - 2014-05-19

+ Fixed `qw()` bug in foreach loops of two required Perl modules. Perl 5.13 and later no 
longer allows you use `qw()` to form lists in for/foreach loops without enclosing function
in another set of parentheses.
+ Fixed very rare bug that could have ignored a possible core gene if only 1 of the 6 core 
proteins for any KOG matched a genomic region. This occurs in the TBLASTN step which looks
to produce a list of 5 potential candidate regions for each KOG. Regions that are only
supported by a KOG protein from a single species are highly unlikely to correspond to real
core genes, but *theoretically* it was possible that this bug could cause 1 or 2 core genes
to be missed.
+ Changed use of `-max_num_description` option for tblastn to instead use `-max_target_seqs`.
+ Fixed typo in documentation which mentioned `-hmm_directory` option instead of `-hmm_profiles`.
+ Use of `-v|--verbose` option now gives much more useful progress information about
how far CEGMA has progressed through each stage. 
+ More detailed output is produced for some steps when using `-v|--verbose` mode.
+ Have moved some functions that are used by multiple scripts to a new `Cegma.pm` module.
+ Tried improving (and harmonizing) various error and warning messages.
+ Added more graceful ways of ending the script prematurely if no core genes are found at 
any stage of analysis (rather than just dying with a cryptic error message).
+ Removed duplication of exon details in GFF output. Lines with GFF feature description of
'Exon' are no longer printed (just 'First', 'Internal', and 'Terminal').
+ Added more default information to completeness.pl output, and refined use of `-v|--verbose` 
option.
+ Updated README file to include a new section about citing CEGMA and added link to the 
new CEGMA FAQ
+ Moved code to a [GitHub repository](https://github.com/KorfLab/CEGMA_v2)
+ Included release notes file as part of code distribution

## [v2.4](http://korflab.ucdavis.edu/Datasets/cegma/cegma_v2.4.010312.tar.gz)

+ Fixed use of -threads option
+ Fixed bug which happened when genome contained no CEGs at all
+ Fixed error that could occur when processing some unusual genewise predictions
+ Fixed bug in calculating completeness statistics for the subset of 248 highly conserved CEGs
+ Made some error messages a little more meaningful
+ Made hmmsearch produce more useful verbose output

## v2.3.1

+ Added --threads details to built-in documentation
+ Used 'config no_ignore_case' options for Getopt::Long to avoid case sensitivity problems with -t/-T options
+ Fixed bug where if no protein alignments were made by genewise, we would not run geneid but still try to run awk on a nonexistent geneid output file. 
+ Made some error messages a little more useful when script dies
+ parsewise script contained a bug that would manifest if a genewise prediction did not start with a  loop state. The code to fix this scenario would not be applied if the last genewise prediction did not start with a  loop.
+ Fixed ortholog statistics generated by final step of pipeline
+ Added key to final table


## [v2.3](http://korflab.ucdavis.edu/Datasets/cegma/cegma_v2.3.190711.tar.gz) 
+ Supports NCBI BLAST+ only (support for WU-BLAST has been removed)
+ New --threads option to let you utliize multiple cores when running tblastn and hmmsearch
+ Now works with latest versions of HMMER and geneid
+ Misc. small improvements to code
