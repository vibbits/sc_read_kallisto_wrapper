#!/usr/bin/perl
# This is a simplified version of the script, for those who have no
# access to the BCL file. It takes as input the fastQ files, which must
# have names of type test_sample_S1_L001_R1_001.fastq.gz

if ($ARGV[0] =~ /^--?h(elp)?$/) {
  HELP:
  print "commandline is :\n\n";
  print "sc_read_kallisto_simplewrapper.pl --dir ??? --index ??? --transcriptome ??? --chemistry ??? --expcells ??? --distance ??? \n\n";
  print "arguments :\n";
  print "  --dir : the folder with the fastQ files\n";
  print "  --index : a kallisto index\n";
  print "  --transcriptome :  the fastA file with the transcriptome used to make the kallisto index\n";
  print "  --chemistry : the chemistry of the library (SC3Pv1 or SC3Pv2, default SC3Pv2)\n";
  print "  --expcells : expected number of cells in a sample (default : 3000)\n";
  print "  --distance : minimum distance between cell barcodes (default : 5)\n";
  print
  exit;
}

no warnings; # needed to avoid warning about goto into a construct
use Getopt::Long;
use Archive::Extract;
use File::Find;
use File::Copy;
use File::Path qw(make_path remove_tree);

# These must be adapted appropriately
$kallisto = 'XXX/kallisto_linux-v0.43.1';
$libdir = 'XXX'; # the location of the Python scripts
$python = '/usr/bin/python';
$Nthreads = 8;

GetOptions(\%options,
  "dir=s",
  "index=s",
  "transcriptome=s",
  "chemistry=s",
  "expcells=i",
  "distance=i"
);
if (not exists $options{dir}) { goto 'HELP' };
if (not exists $options{index}) { goto 'HELP' };
if (exists $options{transcriptome}) { $classes2transcripts = 1 }
if (not exists $options{chemistry}) { $options{chemistry} = 'SC3Pv2' }
if (not exists $options{expcells}) { $options{expcells} = 3000 }
if (not exists $options{distance}) { $options{distance} = 5 }

# set barcode length according to chemistry
if ($options{chemistry} =~ /SC3Pv1/i) {
  $chemistry = 'SC3Pv1' ; $cell_barcode_length = 14;
} elsif ($options{chemistry} =~ /SC3Pv2/i) {
  $chemistry = 'SC3Pv2' ; $cell_barcode_length = 16;
} else {
  print "\n$options{chemistry} is not a known chemistry.\n";
  goto 'HELP';
}

# make list of sample names
opendir DIR, $options{dir} or die "cannot open folder $options{dir}\n";
while ($file = readdir DIR) {
  if ($file !~ /^\./) { # avoid hidden files and the dirs . and ..
    if ($file =~ /(.+)_S\d+_L00\d_([IR]\d)_001.fastq.gz/) {
      $sample_names{$1} = 1;
    } else {
      die "$file is not a proper file name. You need names of type :\n  test_sample_S1_L001_R1_001.fastq.gz\n  where\n    test_sample is the sample name\n    S1 is S followed by a sample number\n    R1 is R or I followed by a number\n";
    }
  }
}
closedir DIR;
@sample_names = keys %sample_names;

# for each sample split the data per cell relying on the cell barcodes,
# then map the data to a transcriptome using kallisto.
$error = 0;
foreach $sample_name (@sample_names) {
  open OUT, ">config.json";
    print OUT "{\n";
    print OUT "    \"NUM_THREADS\": $Nthreads,\n";
    print OUT "    \"EXP_CELLS\": $options{expcells},\n";
    print OUT "    \"SOURCE_DIR\": \"$libdir\",\n";
    print OUT "    \"BASE_DIR\": \"$options{dir}/\",\n";
    print OUT "    \"sample_idx\": \"$sample_name\",\n";
    print OUT "    \"SAVE_DIR\": \"CELL_BARCODES/\",\n";
    print OUT "    \"dmin\": $options{distance},\n";
    print OUT "    \"BARCODE_LENGTH\": $cell_barcode_length,\n";
    print OUT "    \"OUTPUT_DIR\": \"FASTQ_SPLIT_PER_CELL/\"\n";
    print OUT "}";
  close OUT;
  make_path('CELL_BARCODES', 'FASTQ_SPLIT_PER_CELL', "${sample_name}_kallisto");
  if ($chemistry eq 'SC3Pv1') {
    system "$python $libdir/get_cell_barcodes_chem1.py config.json";
    system "$python $libdir/error_correct_and_split_chem1.py config.json";
  } elsif ($chemistry eq 'SC3Pv2') {
    system "$python $libdir/get_cell_barcodes_chem2.py config.json";
    system "$python $libdir/error_correct_and_split_chem2.py config.json";
  }
  move('CELL_BARCODES/umi_barcodes.png',"${sample_name}_umi_barcodes.png");
  system "$kallisto/kallisto pseudo -i $options{index} -o ${sample_name}_kallisto --umi -b FASTQ_SPLIT_PER_CELL/umi_read_list.txt -t $Nthreads 2>> stdout.txt";
    # we redirect STDERR of kallisto to STDOUT to avoid error icon
  remove_tree('config.json', 'CELL_BARCODES', 'FASTQ_SPLIT_PER_CELL');
  if (not -e "${sample_name}_kallisto/matrix.cells") { $error = 1 }
  if (not -e "${sample_name}_kallisto/matrix.ec") { $error = 1 }
  if (not -e "${sample_name}_kallisto/matrix.tsv") { $error = 1 }
}
if ($error) {
  die "Something went wrong with the kallisto mapping. Are you sure you provided a correct index ?\n";
}

# if transcriptome is provided, parse it and extract information about
# genes and transcripts for the sake of preparing input for Seurat
if ($classes2transcripts) {
  open IN, $options{transcriptome}; # perform simple test to check if file is OK
  $line = <IN>;
  # print $line; # for testing
  if ($line !~ /^>/) {
    warn "The transcriptome file does not look OK.\nThe first line should start with a >.\n";
    $classes2transcripts = 0;
  }
  close IN;
  open IN, $options{transcriptome};
  $transcript_ids[0] = 'BOGUS'; # we start effectively at number 1
  while (<IN>) {
    if (/^>([^ \n]+)/) {
      push @transcript_ids, $1;
    }
  }
  close IN;
}

# make from kallisto output input suited for Seurat
foreach $sample_name (@sample_names) {
  make_path("${sample_name}_4Seurat");
  open IN, "${sample_name}_kallisto/matrix.cells";
  open OUT, ">${sample_name}_4Seurat/barcodes.tsv";
  $Ncells = 0;
  while (<IN>) {
    print OUT;
    $Ncells++;
  }
  close IN; close OUT;
  open IN, "${sample_name}_kallisto/matrix.ec";
  open OUT, ">${sample_name}_4Seurat/genes.tsv";
  $Nclasses = 0;
  if ($classes2transcripts) {
    while (<IN>) {
      @fields = split;
      $fields[0]++;
      @items = split /,/, $fields[1];
      for($i=0;$i<=$#items;$i++) {
        $items[$i]++;
        $items[$i] = $transcript_ids[$items[$i]];
      }
      $items = join ',', @items;
      print OUT "$fields[0]\t$items\n";
      $Nclasses++;
    }
  } else {
    while (<IN>) {
      @fields = split;
      $fields[0]++;
      @items = split /,/, $fields[1];
      for($i=0;$i<=$#items;$i++) {
        $items[$i]++;
      }
      $items = join ',', @items;
      print OUT "$fields[0]\t$items\n";
      $Nclasses++;
    }
  }
  close IN; close OUT;
  $Nmappings = 0;
  open IN, "${sample_name}_kallisto/matrix.tsv";
  while (<IN>) {
    $Nmappings++;
  }
  close IN;
  open IN, "${sample_name}_kallisto/matrix.tsv";
  open OUT, ">${sample_name}_4Seurat/matrix.mtx";
  print OUT "%%MatrixMarket matrix coordinate integer general\n%\n";
  print OUT "$Nclasses\t$Ncells\t$Nmappings\n";
  while (<IN>) {
    @fields = split;
    $fields[0]++; $fields[1]++;
    print OUT "$fields[0]\t$fields[1]\t$fields[2]\n";
  }
  close IN; close OUT;
}
