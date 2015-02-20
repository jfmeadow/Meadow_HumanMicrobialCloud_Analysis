if ($#ARGV == 7) {
	$file = $ARGV[0];
	$trim1 = $ARGV[1];
	$trim2 = $ARGV[2];
	$trim3 = $ARGV[3];
	$trim4 = $ARGV[4];
	$trim5 = $ARGV[5];
	$trim6 = $ARGV[6];
	$primer = $ARGV[7];
	
} else {
	die "usage: phased_read_trimmer.pl <inputfile.fastq> <variable trim seq 1> <variable trim seq 2> <variable trim seq 3> <variable trim seq 4> <variable trim seq 5> <variable trim seq 6> <invariant primer sequence> > <outputfile.fastq>  please list trim sequences in order from shortest to longest\n";	
}

$trim1_length = length($trim1);
$trim2_length = length($trim2);
$trim3_length = length($trim3);
$trim4_length = length($trim4);
$trim5_length = length($trim5);
$trim6_length = length($trim6);
$primer_length = length($primer);

#print "$trim1_length\n";
#print "$trim2_length\n";
#print "$trim3_length\n";
#print "$trim4_length\n";
#print "$primer_length\n";


$short_trim2 = substr($trim2,0,$trim1_length);
$short_trim3 = substr($trim3,0,$trim1_length);
$short_trim4 = substr($trim4,0,$trim1_length);
$short_trim5 = substr($trim5,0,$trim1_length);
$short_trim6 = substr($trim6,0,$trim1_length);

#print "$short_trim2\n";
#print "$short_trim3\n";
#print "$short_trim4\n";

open(FILE, "<$file")
	or die;

	while (<FILE>) {

	$ID_line_1 = $_;
		$seq_line = <FILE>;
		$ID_line_2 = <FILE>;
		$Quality_line = <FILE>;
		
		chomp $seq_line;
		chomp $Quality_line;

		$trim_seq = substr($seq_line,0,$trim1_length);
		
		if ($trim_seq eq $trim1)  {
			$correct_count++;
			
			$read_length1 = length($seq_line) - ($trim1_length + $primer_length);
			$read1 = substr($seq_line,($trim1_length + $primer_length),$read_length1);

		$Quality_line_trimmed1 = substr($Quality_line,($trim1_length + $primer_length),$read_length1);
			
			print $ID_line_1;
			print "$read1\n";
			print $ID_line_2;
			print "$Quality_line_trimmed1\n";
	
		} 
		
		elsif ($trim_seq eq $short_trim2)	{
				
			$read_length2 = length($seq_line) - ($trim2_length + $primer_length);
			$read2 = substr($seq_line,($trim2_length + $primer_length),$read_length2);
			
			$Quality_line_trimmed2 = substr($Quality_line,($trim2_length + $primer_length),$read_length2); 
			
			print $ID_line_1;
			print "$read2\n";
			print $ID_line_2;
			print "$Quality_line_trimmed2\n";
		}
	
		elsif ($trim_seq eq $short_trim3)	{
				
			$read_length3 = length($seq_line) - ($trim3_length + $primer_length);
			$read3 = substr($seq_line,($trim3_length + $primer_length),$read_length3);
			
			$Quality_line_trimmed3 = substr($Quality_line,($trim3_length + $primer_length),$read_length3); 
			
			print $ID_line_1;
			print "$read3\n";
			print $ID_line_2;
			print "$Quality_line_trimmed3\n";
		}
	
		elsif ($trim_seq eq $short_trim4)	{
	
			$read_length4 = length($seq_line) - ($trim4_length + $primer_length);
			$read4 = substr($seq_line,($trim4_length + $primer_length),$read_length4);
			
			$Quality_line_trimmed4 = substr($Quality_line,($trim4_length + $primer_length),$read_length4); 
		
			print $ID_line_1;
			print "$read4\n";
			print $ID_line_2;
			print "$Quality_line_trimmed4\n";
		}
	
		elsif ($trim_seq eq $short_trim5)	{
	
			$read_length5 = length($seq_line) - ($trim5_length + $primer_length);
			$read5 = substr($seq_line,($trim5_length + $primer_length),$read_length5);
			
			$Quality_line_trimmed5 = substr($Quality_line,($trim5_length + $primer_length),$read_length5); 
		
			print $ID_line_1;
			print "$read5\n";
			print $ID_line_2;
			print "$Quality_line_trimmed5\n";
		}
	
		elsif ($trim_seq eq $short_trim6)	{
	
			$read_length6 = length($seq_line) - ($trim6_length + $primer_length);
			$read6 = substr($seq_line,($trim6_length + $primer_length),$read_length6);
			
			$Quality_line_trimmed6 = substr($Quality_line,($trim6_length + $primer_length),$read_length6); 
		
			print $ID_line_1;
			print "$read6\n";
			print $ID_line_2;
			print "$Quality_line_trimmed6\n";
		}
		
		else	{
			
			print $ID_line_1;
			print "$seq_line\n";
			print $ID_line_2;
			print "$Quality_line\n";

		}


	}
	close FILE;

