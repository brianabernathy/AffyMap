#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(dirname(abs_path($0))) . '/lib';

use Affy::Utils qw(parse_probes);
use Getopt::Long;
use Pod::Usage;

my $ref_index_file;
my $query_config_file;
my $iupac_scores_file;
my $max_alns = 3;
my $query_allele_a_label = 'Allele_A';
my $query_allele_b_label = 'Allele_B';
my $help;

my %query_config = ();
my %query_probes = ();
my %ref_probes = ();
my %iupac_scores = ();
my %max_iupac_scores = ();
my %max_gt_scores = ();
my %query_scores = ();

parse_args();
parse_probe_config();

foreach my $query_probe_file (sort keys %query_config) {
	my $maf_index = parse_probes($query_probe_file, \%query_config, \%query_probes, undef, undef, $query_allele_a_label, $query_allele_b_label);
}

parse_probes($ref_index_file, undef, \%ref_probes);
parse_iupac_scores();
map_probes();
print_results();

exit(0);


sub print_results {
	print(STDOUT '#query');

	foreach my $index (1..$max_alns) {
		print(STDOUT "\tref${index}\tref${index}_id%");
	}

	print(STDOUT "\n");

	foreach my $q_gt (sort keys %query_scores) {
		my %gt_pcts = ();

		print(STDOUT "$q_gt");

		foreach my $r_gt (keys %{$query_scores{$q_gt}}) {
			my $score = $query_scores{$q_gt}{$r_gt};
			my $max_gt_score = $max_gt_scores{$r_gt};
			my $pct = 0;

			if ($max_gt_score > 0) {
				$pct = $score / $max_gt_score * 100;
			}

			$gt_pcts{$pct}{$r_gt} = $score;
		}

		
		my $count = 0;
		my $print_rec = 1;

		foreach my $pct (reverse sort { $a <=> $b } keys %gt_pcts) {
			foreach my $r_gt (sort keys %{$gt_pcts{$pct}}) {
				if ($print_rec == 1) {
					print(STDOUT "\t$r_gt\t", sprintf("%.2f", $pct), '%');
				}

				$count++;

				if ($count >= $max_alns) {
					$print_rec = 0;
				}
			}
		}

		print(STDOUT "\n");
	}

	return(0);
}


sub map_probes {
	my $shared_probe_count = 0;

	foreach my $probe_id (keys %ref_probes) {
		if (! exists($query_probes{$probe_id})) {
			next();
		}

		my $first_q_gt = 1;

		foreach my $q_gt (keys %{$query_probes{$probe_id}}) {
			my $q_call = $query_probes{$probe_id}{$q_gt};

			foreach my $r_gt (keys %{$ref_probes{$probe_id}}) {
				my $r_call = $ref_probes{$probe_id}{$r_gt};
				my $score = 0;

				if (exists($iupac_scores{$q_call}{$r_call})) {
					$score = $iupac_scores{$q_call}{$r_call};
				}

				$query_scores{$q_gt}{$r_gt} += $score;

				if ($first_q_gt == 1) {
					if (exists($max_iupac_scores{$r_call})) {
						$max_gt_scores{$r_gt} += $max_iupac_scores{$r_call};
					}
				}
			}

			$first_q_gt = 0;
		}

		$shared_probe_count++;
	}

	print(STDERR "shared probe count: $shared_probe_count\n");

	return(0);
}


sub parse_iupac_scores {
	my %col_labels = ();
	my $first_row = 1;

	open(IUPAC, '<', $iupac_scores_file) or error("can't read $iupac_scores_file: $!");

	while (my $line = <IUPAC>) {
		chomp($line);

		my @cols = split(',', $line);

		if ($first_row == 1) {
			foreach my $index (1..$#cols) {
				$col_labels{$index} = $cols[$index];
			}

			$first_row = 0;

			next();
		}

		my $row_label = $cols[0];

		foreach my $index (1..$#cols) {
			my $val = $cols[$index];
			my $col_label = $col_labels{$index};

			$iupac_scores{$row_label}{$col_label} = $val;
			$iupac_scores{$col_label}{$row_label} = $val;

			if (! defined($max_iupac_scores{$row_label}) || $val > $max_iupac_scores{$row_label}) {
				$max_iupac_scores{$row_label} = $val;
			}
		}
	}

	close(IUPAC);

	return(0);
}


sub parse_probe_config {
	my %labels = ();

	open(CONFIG, '<', $query_config_file) or error("can't read $query_config_file: $!");

	while (my $line = <CONFIG>) {
		if ($line =~ /^#/) {
			next();
		}

		chomp($line);

		my ($file, $gt, $label) = split(/\t/, $line);

		if (! defined($file) || ! defined($gt)) {
			error("query config record does not contain a valid file and GT: $line");
		}

		if (! defined($label)) {
			$label = $gt;
		}

		$query_config{$file}{$gt} = $label;
		$labels{$label}++;

		if ($labels{$label} > 1) {
			error("duplicate probe config label found: $label\trecord: $line");
		}
	}

	close(CONFIG);

	return(0);
}


sub error {
	my $msg = shift();

	if (defined($msg)) {
		print(STDERR "error: $msg\n");
	}

	exit(0);
}


sub arg_error {
	my $msg = shift();

	if (defined($msg)) {
		print(STDERR "error: $msg\n\n");
	}

	print("use option -h or --help to display help menu\n");

	exit(0);
}


sub parse_args {
	if ($#ARGV == -1) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	GetOptions ('r|ref=s' => \$ref_index_file,
				'q|query=s' => \$query_config_file,
				'iupac=s' => \$iupac_scores_file,
				'a|alns=i' => \$max_alns,
				'allele_a_label=s' => \$query_allele_a_label,
				'allele_b_label=s' => \$query_allele_b_label,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	if (! defined($ref_index_file)) {
		arg_error('reference index file required');
	}

	if (! defined($query_config_file)) {
		arg_error('query config file required');
	}

	if (! defined($iupac_scores_file)) {
		arg_error('IUPAC scores file required');
	}

	if ($max_alns < 1) {
		arg_error('number of alignments to display must be >= 1');
	}

	return(0);
}


__END__

=head1 NAME

AffyMap.pl

=head1 SYNOPSIS

AffyMap.pl -r ref.index.txt -q query.probe.config.txt --iupac iupac.scores [options] > query.ref.map.txt

=head1 DESCRIPTION

AffyMap.pl aligns query probes to a reference probe index, reporting the best alignments

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -r --ref           reference index file (required)
                      created by AffyIndex.pl

 -q --query         query probe config file (required)
                      created by AffyConfigGen.pl

                      1 genotype per line
                      format: file<tab>file_GT<tab>GT_label
                      GT labels are optional

                      example lines below
                      /my/affy/probes1.txt  geno1.CEL_call_code geno1
                      /my/affy/probes1.txt  geno2.CEL_call_code geno2
                      ...
                      /my/affy/probesX.txt  genoX.CEL_call_code genoX

 --iupac            iupac scores file (required)
                      csv format, see provided file for details

 -a --alns          number of most similar alignments to display
                      default: 3

 --allele_a_label   query allele A column label
                      used to convert numeric or AB input format
                      default: Allele_A

 --allele_b_label   query allele B column label
                      used to convert numeric or AB input format
                      default: Allele_B

 -h --help          display help menu

=cut
