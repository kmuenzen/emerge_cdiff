#!/usr/bin/perl -w

use strict;

##### User Input #####
my $path ="base_dir/";
######################

my %vcf_files;
&get_vcf_files_list;
my %phenotype_files;
&get_phenotype_files_list_ancestry;
&write_shells;

sub get_vcf_files_list
{
    ##### User Input #####
	open (FILE, "vcf_file_list") or die "1 Can't open file\n";
	######################

	while (<FILE>)
	{
		my $line = $_;
		chomp $line;
		$vcf_files{$line} = 0;
	}
	close FILE;
}

sub get_phenotype_files_list_ancestry
{
    ##### User Input #####
	open (FILE, "ancestry_phenotype_files") or die "2 Can't open file\n";
	######################

	while (<FILE>)
	{
		my $line = $_;
		chomp $line;
		$phenotype_files{$line} = 0;
	}
	close FILE;
}

sub write_shells
{
	open (QSUB, ">qsub_run_assoc_max_adjust_GDA_final.sh") or die "Can't open qsub.sh\n";
        print QSUB "#!/bin/bash\n";
        print QSUB "#\$ -S /bin/bash\n";
        print QSUB "#\$ -N masterqsub\n";
        print QSUB "#\$ -l mfree=1G\n";
        print QSUB "#\$ -q cluster.q\n";
        print QSUB "#\$ -cwd\n";
        print QSUB "#\$ -R y\n";
        print QSUB "#\$ -o shell_logs/qsub_run_plink_assoc_max_GDA_final.o\n";
        print QSUB "#\$ -e shell_logs/qsub_run_plink_assoc_max_GDA_final.e\n";

	my $count = 0;
	foreach my $vcf_file (sort keys %vcf_files)
	{
		my $bfile = $vcf_file;
		$bfile =~ s/\.vcf\.gz//gi;
		my $bfile_long = $bfile;
		$bfile_long = "${bfile_long}";
		foreach my $pheno (sort keys %phenotype_files)
		{
			$count++;
			my $shell = $path . "/shells_max_adjust_GDA_final/run_assoc.$count.sh";
			print QSUB "sleep 1 && qsub $shell\n";
			open (SHELL, ">$shell") or die "Can't open $shell\n";
			print SHELL "#!/bin/bash\n";
			print SHELL "#\$ -S /bin/bash\n";
			print SHELL "#\$ -N pGWAS\n";
			print SHELL "#\$ -P jarvik_genomehogger200\n";
			print SHELL "#\$ -l mfree=5G\n";
			print SHELL "#\$ -q cluster.q\n";
			print SHELL "#\$ -cwd\n";
			print SHELL "#\$ -R y\n";
			print SHELL "#\$ -o $path/shell_logs_max_adjust_GDA_final/assoc.$count.o\n";
			print SHELL "#\$ -e $path/shell_logs_max_adjust_GDA_final/assoc.$count.e\n";
			print SHELL "cd $path\n";
            print SHELL "module load plink/1.90b2m\n";
			print SHELL "date\n";
		    print SHELL "\n";
		    print SHELL "#plink --vcf $vcf_file --make-bed --out $path$bfile\n";
		    print SHELL "\n";
		    print SHELL "\n";
			print SHELL "date\n";
			print SHELL "plink \\\n";
			print SHELL "--bfile $bfile_long \\\n";
			print SHELL "--pheno $path$pheno --1 \\\n";
			print SHELL "--all-pheno \\\n";
			print SHELL "--allow-no-sex \\\n";
		    print SHELL "--logistic \\\n";
			print SHELL " hide-covar \\\n";
			print SHELL "--covar $path/cdiff.covar  \\\n";
			print SHELL "--covar-number 1-3,14-15,17-20,23-26 \\\n";
			my $result = "results_max_adjust_GDA_final/$bfile.$pheno.result";
			$result =~ s/nobackup\///;
			print SHELL "--out $result\n";
		        print SHELL "date\n";
		        print SHELL "\n";
		        print SHELL "\n";
		        print SHELL "\n";
		        print SHELL "\n";
		        print SHELL "\n";
			close SHELL;
		}
	}
	close QSUB;
}

