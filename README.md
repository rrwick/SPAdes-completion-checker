# This tool is deprecated!

I have since written [Unicycler](https://github.com/rrwick/Unicycler) which does a good job of circularising genomes. I now use Unicycler instead of hybridSPAdes. If you still need to do this sort of thing, I'd recommend checking out [Circlator](https://sanger-pathogens.github.io/circlator/).

If you're still interested, the original README follows below:


# SPAdes Completion Checker

This tool is intended to aid the user in getting completed bacterial genome sequences from a [hybrid](http://bioinformatics.oxfordjournals.org/content/early/2015/12/16/bioinformatics.btv688) [SPAdes](http://bioinf.spbau.ru/spades) assembly.

After running SPAdes with both short and long read data, the user may aim to find one circular path through the assembly graph for each replicon (chromosome or plasmid) in the sample.  While SPAdes alone often produces long contigs by tracing long paths through the graph, these may not be complete (missing graph segments) or they may overlap with themselves (duplicated graph segments).  Through manual analysis of the assembly graph using  [Bandage](https://github.com/rrwick/Bandage), the user may be able to propose a correct circular path for each replicon.

One important piece of information in such work is the read depth of the graph segments.  A segment which occurs in the chromosomal path once is expected to have a depth close to the chromosome's median depth.  A segment which occurs in a plasmid twice is expected to have a depth of twice the plasmid's median depth.  A segment which occurs once in the chromosome and one in a plasmid is expected to have a depth which is the sum of the chromosome's median depth and the plasmid's median depth.

This tool can be used to evaluate proposed graph paths by comparing the actual read depth of graph segments to the expected read depth given the proposed paths.  It will highlight to the user graph segments which may occur too few or too many times in the proposed paths.  The user may then be able to use this information to improve the paths.

This tool does not automatically fix graph paths or analyse connections or ordering of segments.  It simply reports which graph segments possibly have an incorrect copy number in paths.  Since read depth typically varies over an assembly graph (due to biases in library preparation and sequencing), this tool will likely report false positives - the user must evaluate each reported node.  This tool will be of greater use when the sequencing reads are more evenly distributed across the genome (i.e. constant read depth per replicon).


## Input

This tool requires the following arguments as input:
* A SPAdes assembly graph file.
* A file listing proposed graph paths, one per line.
* The largest k-mer size used in assembly (so the tool can know how much graph segments overlap).
* Two values specifying the range of acceptable ratios of actual to expected read depth.
* (Optionally) a read depth cutoff, segments below which are ignored (default = 1.0).  This is useful for cases where the graph contains low-depth contamination.

Example command:

`spades_completion_checker.py assembly_graph.fastg paths.txt 107 0.8 1.2 path_check`


## Output

The tool produces three output files:
* A table showing the median read depth (by base) of each of the input paths.  This can be used to estimate the abundance of each replicon, e.g. plasmid copy number.
* A table showing anomalous graph segments.
  * For each anomalous segment, the following information is given:
    * Its actual read depth, expected read depth (based on the given paths) and the ratio of these numbers.
    * For each path given, the number of times the segment actually occurs in the path and the number of times it is expected to occur, based on depth.  The expected copy number is only given if a segment occurs only in one path - i.e. segments which occur in multiple paths are not assigned expected copy numbers.
  * Segments are included in this table if both of the following are true:
    * their actual copy number in the paths is not equal to their expected copy number (rounded to the nearest integer).
    * their ratio of actual to expected read depth is outside the acceptable range.
* A CSV file designed for import into Bandage to add informative labels to the graph.

Additionally, the tool will output a list of anomalous graph segment numbers (the same segments in the saved table) to stdout.


## Tuning

By using a larger acceptable ratio range, the tool will be less stringent and report fewer segments:

`spades_completion_checker.py assembly_graph.fastg paths.txt 107 0.5 1.5 path_check`

By using a smaller acceptable ratio range, the tool will be more stringent and report more segments:

`spades_completion_checker.py assembly_graph.fastg paths.txt 107 0.9 1.1 path_check`


## License

GNU General Public License, version 3
