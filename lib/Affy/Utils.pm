package Affy::Utils;

use strict;
use warnings;

use Exporter qw(import);

our @EXPORT_OK = qw(parse_probes);

my %nuc_to_iupac = ();

$nuc_to_iupac{'A/A'} = 'A';
$nuc_to_iupac{'C/C'} = 'C';
$nuc_to_iupac{'G/G'} = 'G';
$nuc_to_iupac{'T/T'} = 'T';
$nuc_to_iupac{'A/C'} = 'M';
$nuc_to_iupac{'C/A'} = 'M';
$nuc_to_iupac{'A/G'} = 'R';
$nuc_to_iupac{'G/A'} = 'R';
$nuc_to_iupac{'A/T'} = 'W';
$nuc_to_iupac{'T/A'} = 'W';
$nuc_to_iupac{'C/G'} = 'S';
$nuc_to_iupac{'G/C'} = 'S';
$nuc_to_iupac{'T/C'} = 'Y';
$nuc_to_iupac{'C/T'} = 'Y';
$nuc_to_iupac{'T/G'} = 'K';
$nuc_to_iupac{'G/T'} = 'K';


sub parse_probes {
	my $probe_file = shift();
	my $config_ref = shift();
	my $probes_ref = shift();
	my $maf_ref = shift();
	my $maf_label = shift();
	my $allele_a_label = shift();
	my $allele_b_label = shift();

	my %col_labels = ();
	my $format;
	my $maf_index;
	my $allele_a_index;
	my $allele_b_index;

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

				if (defined($config_ref)) {
					if (exists($$config_ref{$probe_file}{$val})) {
						$col_labels{$index} = $$config_ref{$probe_file}{$val};
					}
				}

				else {
					$col_labels{$index} = $val;
				}


				if (defined($maf_label) && $val eq $maf_label) {
					$maf_index = $index;
				}

				if (defined($allele_a_label) && $val eq $allele_a_label) {
					$allele_a_index = $index;
				}

				if (defined($allele_b_label) && $val eq $allele_b_label) {
					$allele_b_index = $index;
				}
			}

			next();
		}

		# determine format
		if (! defined($format)) {
			my $call;
			my $col_index;

			foreach my $index (sort keys %col_labels) {
				$call = $cols[$index];
				$col_index = $index;

				$format = record_format($call);

				last();
			}

			if (! defined($format)) {
				error("file format cannot be determined, column: $col_labels{$col_index}, call: \"$call\"");
			}

			if ($format eq 'num' || $format eq 'ab') {
				if (! defined($allele_a_index)) {
					if (! defined($allele_a_label)) {
						error("allele A label is not defined");
					}

					else {
						error("allele A column \"$allele_a_label\" not found");
					}
				}

				elsif (! defined($allele_b_index)) {
					if (! defined($allele_b_label)) {
						error("allele B label is not defined");
					}

					else {
						error("allele B column \"$allele_b_label\" not found");
					}
				}
			}
		}

		# parse probes
		my $probe_id = $cols[0];

		foreach my $index (keys %col_labels) {
			my $call = $cols[$index];
			my $iupac_call;

			if ($format eq 'iupac') {
				$iupac_call = $call;
			}

			else {
				$iupac_call = convert_iupac($call, $cols[$allele_a_index], $cols[$allele_b_index]);
			}

			$$probes_ref{$probe_id}{$col_labels{$index}} = $iupac_call;
		}

		if (defined($maf_index)) {
			push(@{$$maf_ref{$probe_id}}, $cols[$maf_index]);
		}
	}

	close(PROBES);

	return($maf_index);
}


sub record_format {
	my $call = shift();

	if ($call =~ /[ACGT]\/[ACGT]/ || $call eq '---') {
		return('nuc');
	}

	if ($call eq '-1' || $call eq '0' || $call eq '1' || $call eq '2') {
		return('num');
	}

	if ($call eq 'NoCall' || $call eq 'AA' || $call eq 'AB' || $call eq 'BB') {
		return('ab');
	}

	if ($call =~ /^[ACGTMRWSYKN-]$/) {
		return('iupac');
	}

	return(undef);
}


sub convert_iupac {
	my $call = shift();
	my $allele_a = shift();
	my $allele_b = shift();

	if ($call eq '---' || $call eq '-1' || $call eq 'NoCall') {
		return('-');
	}

	if (exists($nuc_to_iupac{$call})) {
		return($nuc_to_iupac{$call});
	}

	if (! defined($allele_a) || ! defined($allele_b)) {
		return(undef);
	}

	if ($call eq '0' || $call eq 'AA') {
		return($nuc_to_iupac{"${allele_a}/${allele_a}"});
	}

	if ($call eq '1' || $call eq 'AB') {
		return($nuc_to_iupac{"${allele_a}/${allele_b}"});
	}

	if ($call eq '2' || $call eq 'BB') {
		return($nuc_to_iupac{"${allele_b}/${allele_b}"});
	}

	return(undef);
}


sub error {
	my $msg = shift();

	if (defined($msg)) {
		print(STDERR "error: $msg\n");
	}

	exit(0);
}


1;
