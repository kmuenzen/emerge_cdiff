#!/usr/bin/perl -w

use strict;

##### User Inputs #####
my $path = "base_dir/";
my $results_file_list = "gwas_results_list";
my $variant_quality_results_file = "quality_file";
my $variant_freq_file_list = "frequency_file";
#######################

my %chr_pos_variant_rsq_quality;
my %chr_pos_variant_frequency;
my %results_files;
&get_quality;
&get_frequency;
&get_results_list;
&format_results;
&write_shells;

sub get_frequency
{
	open (FILELIST, $variant_freq_file_list) or die "1 Can't open file\n";
	while (<FILELIST>)
	{
		my $file = $_;
		chomp $file;
		open (FILE, $file) or die "2 Can't open file\n";
		while (<FILE>)
		{
			my $line = $_;
			chomp $line;
			my @contents = split(/\s+/, $line);

			if (exists $contents[5] && $line !~ m/^CHROM/)
			{
				my $chr = $contents[0];
				my $pos = $contents[1];
				my $allele_freq = $contents[5];
				my @allele_freq_contents = split(/:/, $allele_freq);
				my $freq = $allele_freq_contents[1];
				my $minor_allele = $allele_freq_contents[0];

				my $major_allele_freq = $contents[4];
				my @major_allele_freq_contents = split(/:/, $major_allele_freq);
				my $major_allele = $major_allele_freq_contents[0];

				if ($freq >= 0.05)
				{
					$chr_pos_variant_frequency{$chr}{$pos}{'f'} = $freq;
					$chr_pos_variant_frequency{$chr}{$pos}{'mi'} = $minor_allele;
					$chr_pos_variant_frequency{$chr}{$pos}{'ma'} = $major_allele;
				}
			}
		}
		close FILE;
	}
	print('Frequency function finished\n');
	close FILELIST;
}

sub get_quality
{
	open (FILE, $variant_quality_results_file) or die "3 Can't open file\n";
	while (<FILE>)
	{
		my $line = $_;
		chomp $line;
		my @contents = split(/\s+/, $line);
		my $chr_pos = $contents[0];
		my $minor_allele = $contents[1];
		$minor_allele =~ s/MA\://;
		my $genotyped_count = $contents[3];
		my $imputed_count = $contents[5];
		my $rsq = $contents[8];
		my @chr_pos_contents = split(/:/, $chr_pos);
		my $chr = $chr_pos_contents[0];
		my $pos = $chr_pos_contents[1];
		if ($rsq >= 0.3)
		{
			$chr_pos_variant_rsq_quality{$chr}{$pos}{'r'} = $rsq;
			$chr_pos_variant_rsq_quality{$chr}{$pos}{'m'} = $minor_allele;
			$chr_pos_variant_rsq_quality{$chr}{$pos}{'g'} = $genotyped_count;
			$chr_pos_variant_rsq_quality{$chr}{$pos}{'i'} = $imputed_count;
		}
	}
	print('Quality function finished\n');
	close FILE;
}

sub get_results_list
{
	open (FILE, $results_file_list) or die "4 Can't open file\n";
	while (<FILE>)
	{
		my $line = $_;
		chomp $line;
		my @contents = split(/\./, $line);
		my $chr = $contents[0];
		my $start_pos = $contents[1];
		$chr =~ s/results_max_adjust_GDA_final\/chr//gi;
		$results_files{$chr}{$start_pos} = $line;
	}
	print('Results list function finished\n');
	close FILE;
}

sub format_results
{
	open (OUTFILE, ">$results_file_list.formatted.variant_rsq_mean_quality_gteq_0.3_freq_gteq_0.05") or die "5 Can't open file\n";
	open (OUTFILE1, ">$results_file_list.formatted.variant_rsq_mean_quality_gteq_0.3_freq_gteq_0.05_p_lteq_5x10-8") or die "6 Can't open file\n";
	print OUTFILE1 "chromosome\tposition\tp\tOR\trsq\tfreq\tminor_allele\tmajor_allele\tgenotyped_batches\timputed_batches\n";
	print OUTFILE "chromosome\tp\n";
	foreach my $chr (sort {$a <=> $b} keys %results_files)
	{
		foreach my $start_pos (sort {$a <=> $b} keys %{$results_files{$chr}})
		{
			my $file1 = $results_files{$chr}{$start_pos};
			my $file = $file1 =~ s/P1/CASE_CONTROL/r;
			print "working on $file\n";
			open (FILE, $file) or die "7 Can't open file\n";
			while (<FILE>)
			{
				my $line = $_;
				chomp $line;
				my @contents = split(/\s+/, $line);
				my $chr = $contents[1];
				my $pos = $contents[3];
				my $odds_ratio = $contents[7];
				my $p = $contents[9];
				if ($p eq "NA" || $p eq "P")
				{}
				elsif (exists $chr_pos_variant_rsq_quality{$chr}{$pos} && exists $chr_pos_variant_frequency{$chr}{$pos})
				{
					if ($p <= 0.00000005)
					{
						my $rsq = $chr_pos_variant_rsq_quality{$chr}{$pos}{'r'};
						my $genotyped_count =$chr_pos_variant_rsq_quality{$chr}{$pos}{'g'};
						my $imputed_count = $chr_pos_variant_rsq_quality{$chr}{$pos}{'i'};
						my $freq = $chr_pos_variant_frequency{$chr}{$pos}{'f'};
						my $minor_allele = $chr_pos_variant_frequency{$chr}{$pos}{'mi'};
						my $major_allele = $chr_pos_variant_frequency{$chr}{$pos}{'ma'};
						print OUTFILE1 "$chr\t$pos\t$p\t$odds_ratio\t$rsq\t$freq\t$minor_allele\t$major_allele\t$genotyped_count\t$imputed_count\n";
					}
					print OUTFILE "$chr\t$p\n";
				}
			}
			close FILE;
		}
	}
	print('Results formatting function finished\n');
	close OUTFILE;
}

sub write_shells
{
	my $shell = $path . "/run_$results_file_list.rsq_0.3_var_qual_freq_gteq_0.05.graph.sh";
	print "$shell\n";

		open (SHELL, ">$shell") or die "Can't open $shell\n";
            print SHELL "#!/bin/bash\n";
            print SHELL "#\$ -S /bin/bash\n";
            print SHELL "#\$ -N graphGWAS\n";
            print SHELL "#\$ -l mfree=25G\n";
            print SHELL "#\$ -q cluster.q\n";
            print SHELL "#\$ -cwd\n";
            print SHELL "#\$ -R y\n";
            print SHELL "#\$ -o $path/$results_file_list.rsq_0.3_var_qual_freq_gteq_0.05.graph.o\n";
            print SHELL "#\$ -e $path/$results_file_list.rsq_0.3_var_qual_freq_gteq_0.05.graph.e\n";
            print SHELL "cd $path\n";
            print SHELL "date\n";
            print SHELL "module load R/3.6.1\n";
            print SHELL "module load hdf5/1.8.13\n";
            my $rscript = $path . "/run_$results_file_list.rsq_0.3_var_qual_freq_gteq_0.05.graph.R";
            print "$rscript\n";
            print SHELL "R --no-save < $rscript\n";
            print SHELL "\n";
            print SHELL "date\n";
            print SHELL "\n";
		close SHELL;

		open (RSCRIPT, ">$rscript") or die "Can't open RSCRIPT\n";
	        print RSCRIPT "\n";
	        print RSCRIPT "library(\"GWASTools\")\n";
	        print RSCRIPT "\n";
	        print RSCRIPT "results <- read.table(\"$results_file_list.formatted.variant_rsq_mean_quality_gteq_0.3_freq_gteq_0.05\",sep=\"\\t\", header=T)\n";
	        print RSCRIPT "\n";
		    print RSCRIPT "png(file=\"$results_file_list.variant_rsq_mean_quality_gteq_0.3_freq_gteq_0.05.manhattan.png\", width = 1200, height = 450, type = \"cairo\")\n";
	        print RSCRIPT "manhattanPlot(results\$p, results\$chromosome, trunc.lines=TRUE, signif=5e-8)\n";
		    print RSCRIPT "dev.off()\n";
	        print RSCRIPT "\n";
		    print RSCRIPT "png(file=\"$results_file_list.variant_rsq_mean_quality_gteq_0.3_freq_gteq_0.05.qq.png\", width = 400, height = 450, type = \"cairo\")\n";
	        print RSCRIPT "qqPlot(results\$p, truncate=FALSE, ci=TRUE)\n";
		    print RSCRIPT "dev.off()\n";
	        print RSCRIPT "\n";
	        print RSCRIPT "chisq <- qchisq(1 - results\$p, 1)\n";
	        print RSCRIPT "lambda_gc <- median(chisq)/qchisq(0.5, 1)\n";
	        print RSCRIPT "lambda_gc\n";
	        print RSCRIPT "\n";
	        print RSCRIPT "library(\"GenABEL\")\n";
	        print RSCRIPT "genabel_lambda_gc <- estlambda(Pvalues, method=\"median\")\n";
	        print RSCRIPT "genabel_lambda_gc\n";
	        print RSCRIPT "\n";
	        print RSCRIPT "\n";
		close RSCRIPT;
		print('Script writing function finished\n');
}

