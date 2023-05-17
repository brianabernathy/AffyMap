#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(dirname(abs_path($0))) . '/lib';

use Affy::Utils qw(parse_probes);
use Getopt::Long;
use Pod::Usage;

my $probe_config_file;
my $monomorphic;
my $maf_filter;
my $max_miss_probe_pct;
my $max_miss_gt_pct;
my $retain_probes_file;
my $retain_gts_file;
my $maf_label = 'MinorAlleleFrequency';
my $allele_a_label = 'Allele_A';
my $allele_b_label = 'Allele_B';
my $invalid_call = 'na';
my $verbose;
my $help;

my %config = ();
my %probes = ();
my %maf = ();
my %gts = ();
my %retain_probes = ();
my %retain_gts = ();

parse_args();
parse_probe_config();

if (defined($retain_probes_file)) {
	parse_retain_probes();
}

if (defined($retain_gts_file)) {
	parse_retain_gts();
}

foreach my $probe_file (sort keys %config) {
	my $maf_index = parse_probes($probe_file, \%config, \%probes, \%maf, $maf_label, $allele_a_label, $allele_b_label);

	if (! defined($maf_index)) {
		error("MAF column \"$maf_label\" not found");
	}
}

set_invalid_calls();

if (defined($monomorphic) || defined($maf_filter) || defined($max_miss_probe_pct) || defined($max_miss_gt_pct)) {
	filter_probes_gts();
}

print_probes();

exit(0);


sub print_probes {
	my $print_header = 1;

	foreach my $probe_id (sort keys %probes) {
		if ($print_header == 1) {
			print(STDOUT join("\t", 'probeset_id', sort keys %{$probes{$probe_id}}), "\n");

			$print_header = 0;
		}

		print(STDOUT "$probe_id");

		foreach my $gt (sort keys %{$probes{$probe_id}}) {
			print(STDOUT "\t$probes{$probe_id}{$gt}");
		}

		print(STDOUT "\n");
	}

	return(0);
}


sub filter_probes_gts {
	my %monomorphic_probes = ();
	my %maf_filter_probes = ();
	my %miss_data_probes = ();
	my %miss_data_gts = ();
	my %gt_calls = ();
	my $probe_count = scalar keys %probes;
	my $gt_count = scalar keys %gts;

	foreach my $probe_id (sort keys %probes) {
		my %probe_calls = ();

		foreach my $gt (keys %{$probes{$probe_id}}) {
			my $call = $probes{$probe_id}{$gt};

			$probe_calls{$call}++;
			$gt_calls{$gt}{$call}++;
		}

		# check monomorphic
		if (defined($monomorphic)) {
			my %explicit_calls = ();

			foreach my $call (keys %probe_calls) {
				if ($call ne 'N' && $call ne $invalid_call) {
					$explicit_calls{$call}++;
				}
			}

			if (scalar keys %explicit_calls == 1) {
				$monomorphic_probes{$probe_id}++;

				if (defined($verbose)) {
					print(STDERR "monomorphic probe: $probe_id\n");
				}
			}
		}

		# check maf
		if (defined($maf_filter)) {
			foreach my $maf (@{$maf{$probe_id}}) {
				if ($maf >= $maf_filter) {
					$maf_filter_probes{$probe_id}++;

					if (defined($verbose)) {
						print(STDERR "MAF filter probe: $probe_id\tMAF(s): ", join(',', @{$maf{$probe_id}}), "\n");
					}
				}
			}
		}

		# check missing probe data
		if (defined($max_miss_probe_pct)) {
			my $miss_probe_count = 0;

			foreach my $call (keys %probe_calls) {
				if ($call eq 'N' || $call eq $invalid_call) {
					$miss_probe_count++;
				}
			}

			my $miss_probe_pct = 0;

			if ($gt_count > 0) {
				$miss_probe_pct = $miss_probe_count / $gt_count * 100;
			}

			if ($miss_probe_pct >= $max_miss_probe_pct) {
				$miss_data_probes{$probe_id}++;

				if (defined($verbose)) {
					print(STDERR "missing data probe: $probe_id ($miss_probe_pct%)\n");
				}
			}
		}
	}

	# check missing gt data
	if (defined($max_miss_gt_pct)) {
		foreach my $gt (keys %gts) {
			my $miss_gt_count = 0;

			foreach my $probe_id (keys %probes) {
				my $call = $probes{$probe_id}{$gt};

				if ($call eq 'N' || $call eq $invalid_call) {
					$miss_gt_count++;
				}
			}

			my $miss_gt_pct = 0;

			if ($probe_count > 0) {
				$miss_gt_pct = $miss_gt_count / $probe_count * 100;
			}

			if ($miss_gt_pct >= $max_miss_gt_pct) {
				$miss_data_gts{$gt}++;

				if (defined($verbose)) {
					print(STDERR "missing data gt: $gt ($miss_gt_pct%)\n");
				}
			}
		}
	}


	my $filtered_probe_count = 0;
	my %filtered_gts = ();

	foreach my $probe_id (keys %probes) {
		if ((exists($monomorphic_probes{$probe_id}) || exists($maf_filter_probes{$probe_id}) || exists($miss_data_probes{$probe_id})) && ! exists($retain_probes{$probe_id})) {
			delete($probes{$probe_id});

			$filtered_probe_count++;
		}

		else {
			foreach my $gt (keys %miss_data_gts) {
				if (! exists($retain_gts{$gt})) {
					delete($probes{$probe_id}{$gt});
					delete($gts{$gt});

					$filtered_gts{$gt}++;
				}
			}
		}
	}

	print(STDERR "monomorphic probes:    ", scalar keys %monomorphic_probes, "\n");
	print(STDERR "MAF filter probes:     ", scalar keys %maf_filter_probes, "\n");
	print(STDERR "missing data probes:   ", scalar keys %miss_data_probes, "\n");
	print(STDERR "total filtered probes: $filtered_probe_count\n");
	print(STDERR "total retained probes: ", scalar keys %probes, "\n\n"); 
	print(STDERR "missing data GTs:      ", scalar keys %miss_data_gts, "\n");
	print(STDERR "total filtered GTs:    ", scalar keys %filtered_gts, "\n");
	print(STDERR "total retained GTs:    ", scalar keys %gts, "\n");

	return(0);
}


sub set_invalid_calls {
	foreach my $probe_id (keys %probes) {
		foreach my $gt (keys %{$probes{$probe_id}}) {
			$gts{$gt}++;

			if (! defined($probes{$probe_id}{$gt})) {
				$probes{$probe_id}{$gt} = $invalid_call;
			}
		}
	}

	return(0);
}


sub parse_probe_config {
	my %labels = ();

	open(CONFIG, '<', $probe_config_file) or error("can't read $probe_config_file: $!");

	while (my $line = <CONFIG>) {
		if ($line =~ /^#/) {
			next();
		}

		chomp($line);

		my ($file, $gt, $label) = split(/\t/, $line);

		if (! defined($file) || ! defined($gt)) {
			error("probe config record does not contain a valid file and GT: $line");
		}

		if (! defined($label)) {
			$label = $gt;
		}

		$config{$file}{$gt} = $label;
		$labels{$label}++;

		if ($labels{$label} > 1) {
			error("duplicate probe config label found: $label\trecord: $line");
		}
	}

	close(CONFIG);

	return(0);
}


sub parse_retain_probes {
	open(PROBES, '<', $retain_probes_file) or error("can't read $retain_probes_file: $!");

	while (my $line = <PROBES>) {
		chomp($line);

		$retain_probes{$line}++;
	}

	close(PROBES);
}


sub parse_retain_gts {
	open(GTS, '<', $retain_gts_file) or error("can't read $retain_gts_file: $!");

	while (my $line = <GTS>) {
		chomp($line);

		my $gt = $line;

		foreach my $file (sort keys %config) {
			if (exists($config{$file}{$gt})) {
				$retain_gts{$config{$file}{$gt}}++;
			}

			else {
				$retain_gts{$gt}++;
			}
		}
	}

	close(GTS);

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

	GetOptions ('c|config=s' => \$probe_config_file,
				'monomorphic' => \$monomorphic,
				'maf_filter=f' => \$maf_filter,
				'miss_probe_pct=f' => \$max_miss_probe_pct,
				'miss_gt_pct=f' => \$max_miss_gt_pct,
				'retain_probes=s' => \$retain_probes_file,
				'retain_gts=s' => \$retain_gts_file,
				'maf_label=s' => \$maf_label,
				'allele_a_label=s' => \$allele_a_label,
				'allele_b_label=s' => \$allele_b_label,
				'invalid_call=s' => \$invalid_call,
				'v|verbose' => \$verbose,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	if (! defined($probe_config_file)) {
		arg_error('probe file required');
	}

	if (defined($maf_filter) && $maf_filter <= 0) {
		arg_error('maf_filter value must be > 0');
	}

	if (defined($max_miss_probe_pct) && $max_miss_probe_pct <= 0) {
		arg_error('miss_probe_pct value must be > 0');
	}

	if (defined($max_miss_gt_pct) && $max_miss_gt_pct <= 0) {
		arg_error('miss_gt_pct value must be > 0');
	}

	$invalid_call =~ s/\t/_/g;

	return(0);
}


__END__

=head1 NAME

AffyIndex.pl

=head1 SYNOPSIS

AffyIndex.pl -c probe.config.txt [options] > ref.index.txt

=head1 DESCRIPTION

AffyIndex.pl produces an IUPAC formatted index suitable for alignment with AffyMap.pl

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

=head2 input

 -c --config        probe config file (required)
                      1 genotype per line
                      format: file<tab>file_GT<tab>GT_label
                      GT labels are optional

                      example lines below
                      /my/affy/probes1.txt	geno1.CEL_call_code	geno1
                      /my/affy/probes1.txt	geno2.CEL_call_code	geno2
                      ...
                      /my/affy/probesX.txt	genoX.CEL_call_code	genoX

=head2 filtering

By default, no filtering is performed. However, it is recommended to
enable all options below for optimal results.

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

=head2 column labels

 --maf_label        minor allele frequency column label
                      default: MinorAlleleFrequency

 --allele_a_label   allele A column label
                      used to convert numeric or AB input format
                      default: Allele_A

 --allele_b_label   allele B column label
                      used to convert numeric or AB input format
                      default: Allele_B

=head2 output

 --invalid_call     invalid call value
                      default: 'na'

 -v --verbose       enable verbose logging to stderr
                      list all filtered probes and GTs, etc...

=head2 help

 -h --help          display help menu

=cut
