# AffyMap

AffyMap is able to index and align Affymetrix microarray datasets produced by the [Axiom Analysis Suite](https://www.thermofisher.com/us/en/home/technical-resources/software-downloads/axiom-analysis-suite.html "Axiom Analysis Suite").

## basic usage

- Input microarray data files, genotypes and genotype labels are defined in a config file. Config files can be created manually or automatically using [AffyConfig.pl](https://github.com/brianabernathy/AffyMap/blob/main/bin/AffyConfig.pl "AffyConfig.pl").

- The reference config file and indexing options are provided to [AffyIndex.pl](https://github.com/brianabernathy/AffyMap/blob/main/bin/AffyIndex.pl "AffyIndex.pl") to generate a reference index.

- The reference index file, query config file and mapping options are provided to [AffyMap.pl](https://github.com/brianabernathy/AffyMap/blob/main/bin/AffyMap.pl "AffyMap.pl") to compare genotype data. For each query genotype, a set of the most similar reference genotypes are listed.

- Detailed command options are available by running the command without options or the -h or --help option. [Command options](https://github.com/brianabernathy/AffyMap#commandoptions "command options") are also listed below.

## test data

coming soon...

## command options


