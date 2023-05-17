#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

my $probe_file;
my $gt_pattern;
my $sample_file;
my $include_well = 0;
my $skip_no_att = 0;
my $help;

my %labels = ();
my $cel_pattern = '.CEL';

parse_args();

if (defined($sample_file)) {
	parse_samples();
}

parse_probes();

exit(0);

sub parse_probes {
	open(PROBES, '<', $probe_file) or error("can't read $probe_file: $!");

	while (my $line = <PROBES>) {
		if ($line =~ /^#/) {
			next();
		}

		chomp($line);

		my @cols = split(/\t/, $line);

		# parse header
		if ($cols[0] eq 'probeset_id') {
			foreach my $index (1..$#cols) {
				my $val = $cols[$index];
				my $print_gt = 1;

				if (defined($gt_pattern)) {
					if ($val !~ /$gt_pattern/) {
						$print_gt = 0;
					}
				}

				if ($print_gt == 1) {
					my $header = $val;
					my $label;

					$header =~ s/$cel_pattern.*$//;

					if (exists($labels{$header})) {
						$label = $labels{$header};
					}

					elsif ($skip_no_att == 1) {
						next();
					}

					print(STDOUT "$probe_file\t$val");

					if (defined($label)) {
						print(STDOUT "\t$label");
					}

					print(STDOUT "\n");
				}
			}

			last();
		}
	}

	return(0);
}


sub parse_samples {
	open(SAMPLES, '<', $sample_file) or error("can't read $sample_file: $!");

	while (my $line = <SAMPLES>) {
		if ($line =~ /^Sample/) {
			next();
		}

		chomp($line);

		my @cols = split(/\t/, $line);
		my $header = $cols[0];
		my $well = $cols[2];
		my $sample = $cols[3];
		my $label = $sample;

		if ($include_well == 1) {
			$label .= ".${well}";
		}

		$header =~ s/$cel_pattern.*$//;

		$labels{$header} = $label;
	}

	close(SAMPLES);

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

	GetOptions ('p|probe=s' => \$probe_file,
				'g|gt=s' => \$gt_pattern,
				's|sample=s' => \$sample_file,
				'w|well' => \$include_well,
				'skip' => \$skip_no_att,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -output => '-', -verbose => 1);
	}

	if (! defined($probe_file)) {
		arg_error('probe file required');
	}

	return(0);
}


__END__

=head1 NAME

AffyConfig.pl

=head1 SYNOPSIS

AffyConfig.pl -p probes.txt [options] > probes.config.txt

=head1 DESCRIPTION

AffyConfig.pl produces a probe config file for use with AffyIndex.pl

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

=head2 input

 -p --probe   probe file (required)

 -g --gt      genotype pattern
                used to include only header fields matching the pattern
                ex: 'CEL_call_code'

 -s --sample  sample attributes file

 -w --well    append well row/col to sample names
                useful to make sample names unique

 --skip       skip samples without sample attribute records

 -h --help    display help menu

=cut
