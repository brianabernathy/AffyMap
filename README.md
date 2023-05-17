# AffyMap

AffyMap is able to index and align Affymetrix microarray datasets produced by the [Axiom Analysis Suite](https://www.thermofisher.com/us/en/home/technical-resources/software-downloads/axiom-analysis-suite.html "Axiom Analysis Suite").

## basic usage

- Input microarray data files, genotypes and genotype labels are defined in a config file. Config files can be created manually or automatically using [AffyConfig.pl](https://github.com/brianabernathy/AffyMap/blob/main/bin/AffyConfig.pl "AffyConfig.pl").

- The reference config file and indexing options are provided to [AffyIndex.pl](https://github.com/brianabernathy/AffyMap/blob/main/bin/AffyIndex.pl "AffyIndex.pl") to generate a reference index.

- The reference index file, query config file, [scoring matrix](https://github.com/brianabernathy/AffyMap/blob/main/misc/iupac.scores.csv "scoring matrix") and mapping options are provided to [AffyMap.pl](https://github.com/brianabernathy/AffyMap/blob/main/bin/AffyMap.pl "AffyMap.pl") to compare genotype data. For each query genotype, a set of the most similar reference genotypes are listed.

- Detailed command options are available by running the command without options or the -h or --help option. [Command options](https://github.com/brianabernathy/AffyMap#commandoptions "command options") are also listed below.

## test data

coming soon...

## command options

### AffyConfig.pl

Usage:

    AffyConfig.pl -p probes.txt [options] > probes.config.txt

Options:

     -p --probe   probe file (required)

     -g --gt      genotype pattern
                    used to include only header fields matching the pattern
                    ex: 'CEL_call_code'

     -s --sample  sample attributes file

     -w --well    append well row/col to sample names
                    useful to make sample names unique

     --skip       skip samples without sample attribute records

     -h --help    display help menu

### AffyIndex.pl 

Usage:

    AffyIndex.pl -c probe.config.txt [options] > ref.index.txt

Options:

  input:

     -c --config        probe config file (required)
                          1 genotype per line
                          format: file<tab>file_GT<tab>GT_label
                          GT labels are optional

                          example lines below
                          /my/affy/probes1.txt      geno1.CEL_call_code     geno1
                          /my/affy/probes1.txt      geno2.CEL_call_code     geno2
                          ...
                          /my/affy/probesX.txt      genoX.CEL_call_code     genoX

  filtering:

    By default, no filtering is performed. However, it is recommended to enable all options below for optimal results.

     --monomorphic      filter monomorphic probes

     --maf_filter       filter probes with >= maf_filter value
                          recommended: 0.05

     --miss_probe_pct   filter probes with >= miss_probe_pct missing data
                          recommended: 20

     --miss_gt_pct      filter genotypes with >= miss_gt_pct missing data
                          recommended: 20

     --retain_probes    file containing probes to retain, regardless of 
                          filtering status
                          1 probe per line

     --retain_gts       file containing genotypes to retain, regardless
                          of filtering status
                          1 genotype per line

  column labels:

     --maf_label        minor allele frequency column label
                          default: MinorAlleleFrequency

     --allele_a_label   allele A column label
                          used to convert numeric or AB input format
                          default: Allele_A

     --allele_b_label   allele B column label
                          used to convert numeric or AB input format
                          default: Allele_B

  output:

     --invalid_call     invalid call value
                          default: 'na'

     -v --verbose       enable verbose logging to stderr
                          list all filtered probes and GTs, etc...

  help:

     -h --help          display help menu

### AffyMap.pl

Usage:

    AffyMap.pl -r ref.index.txt -q query.probe.config.txt --iupac iupac.scores [options] > query.ref.map.txt

Options:

     -r --ref           reference index file (required)
                          created by AffyIndex.pl

     -q --query         query probe config file (required)
                          created by AffyConfigGen.pl

     --iupac            iupac scores file (required)
                          csv format, see provided file for details

     -a --alns          number of best scoring alignments to display
                          default: 3

     --allele_a_label   query allele A column label
                          used to convert numeric or AB input format
                          default: Allele_A

     --allele_b_label   query allele B column label
                          used to convert numeric or AB input format
                          default: Allele_B

     -h --help          display help menu
