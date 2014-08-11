#!/usr/bin/perl
=doc
a [dot] llave [at] irri.org 28APR2014-1445

A basic VCF to GFF3 conversion tool. Includes AF and DP (allele frequency and
read depth respectively) from INFO column of VCF.

A spin-off from /home/applications/vaast_converter.pl
which was referred from http://gmod.827538.n3.nabble.com/vcf-to-gff3-td3434214.html

Copyright (C) 2014  Abraham Darius S. Llave under the C4 Rice Project, International Rice Research Institute

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
=cut

use strict;
use warnings;
use Getopt::Long;
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
# specify your source
my $sourceName = "";
# used for identifying the feature's source as can be gleamed from the feature ID/name
my $featureIDStart = "";
# this would be part of the feature ID/name to identify the dataset.
my $datasetTag = "";
my $countSNP = 0;
my $countINDEL = 0;
my $debugMode = 0;
my $noUnderscore = 0;
my $separator = "_";
my $usage = "

Synopsis:

vcf_to_gff3.pl [options] file1.vcf [file2.vcf...]



Options:
 --help          Optional. Display this help.
 --path          Optional. Default directory to write results. Defaults to current/`pwd`.
 --orgsource    Optional. A string describing the organization/institute that is source of the features, for column 2 of GFF file. Defaults to a dot.
 --orgsource2    Optional. A string describing the organization/institute that is source of the features, put in the beginning of the name. Default to 'loc'.
 --string_tag    Required. A string of minlength 1 to describe the features.
 --nounderscore  Optional. Do not separate the --orgsource from the rest of the characters in the features' name/ID. Default false. 
 
Description:

This script will convert VCF files to GFF3.
Take note, that outputted info are limited to location, variance, allele frequency and read depth.
By default structure of feature ID/name in quasi-regex: (orgsource|'loc')(undesrcore)?<(I|S)><string_tag><order of feature among S/I>
";


my ($help, $build, $append, $path);

my $opt_success = GetOptions('help'    => \$help,
                             'nounderscore' => \$noUnderscore,
                             'orgsource=s' => \$sourceName,                             
                             'orgsource2=s' => \$featureIDStart,
					         'path=s'  => \$path,
                             'string_tag=s' => \$datasetTag
);

die $usage if $help || ! $opt_success;

$path ||= './';
$featureIDStart ||= "loc";
$sourceName ||= ".";
if ($noUnderscore) { $separator = "" };

handle_message('FATAL', 'directory_does_not_exist', $path)
    unless -e $path;
handle_message('FATAL', 'path_is_not_a_directory', $path)
    unless -d $path;
handle_message('FATAL', 'directory_is_not_writable', $path)
    unless -w $path;
handle_message('FATAL', 'string_tag should be at least 1 char', $datasetTag)
    unless ( length($datasetTag) > 0 );
handle_message('NOTICE', 'string_tag ', $datasetTag);
my @vcf_files = @ARGV;

if (! @vcf_files) {
  print $usage;
  handle_message('FATAL', 'no_vcf_files_to_convert');
}

for my $vcf_file (@vcf_files) {
  if ( $vcf_file !~ /\./ ){
    handle_message('FATAL', 'VCF file should have an extension name (i.e. dot)');
  }

  ( my $vcf_file_withoutExt = $vcf_file ) =~ s/\.\w+$//;
  my $gff_file = $vcf_file_withoutExt . ".gff";
  handle_message('FATAL', 'file_does_not_exist', $vcf_file) unless -e $vcf_file;
  open(my $IN, '<', $vcf_file) or
    handle_message('FATAL', 'cant_open_file_for_reading', $vcf_file);
  open(my $OUThandle, ">", $gff_file) or
    handle_message('FATAL', 'cant_open_file_for_writing', $gff_file);
   print $OUThandle "##gff-version 3\n";
   print $OUThandle "# Generated " . localtime() . ' via https://github.com/abrahamdsl/vcf_to_gff3' . "\n";
   

=debug
  # a [dot] llave [at] irri.org THIS MIGHT NOT BE NEEDED IN OUR CASE
  my $header;
  while (my $line = <$IN>) {
    if ($line =~ /^#CHROM/) {
      $header = $line;
      last;
    }
  }

  handle_message('FATAL', 'no_header_line_found',
		 'File was read, but no header line was found')
    unless $header;
  my @ids = split /\s/, $header;
  splice(@ids, 0, 9);
=cut
  
  handle_message( 'NOTICE', 'Processing file', $vcf_file );
  handle_message( 'NOTICE', 'Writing to file', $gff_file );
  while (my $line = <$IN>) {    
    if( $line =~ /^#/ or $line =~ /^super/ ){ next; }
    my @data = split /\t/, $line;
    #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE1
    my %record;

    @record{qw(chrom pos id ref alt qual filter info format)} =
      splice(@data, 0, 9);

    my @feature_ids = split /;/, $record{id}; # more often than not this is just '.'
    my $feature_id = shift @feature_ids;
    my $seqid      = $record{chrom};
    my $source;
    my $type;
    my $start      = $record{pos};
    my $end        = $start;       
    my $score      = $record{qual};

    # In VCF, variants are always in forward representation.
    # Citation: http://www.1000genomes.org/faq/what-strand-are-variants-your-vcf-file
    my $strand     = '+';                     
    my $phase      = '.';                     # not needed to specify.
    my $note;
    my $featureID;
    my $reference_seq = uc $record{ref};
 
    my @infos = split( /;/, $record{info} );
    # Which details in the INFO field do you want to include?
    my @wanted = qw(AF DP);
    my $wantedCount = scalar( @wanted );
    my $count = 0;
    my %includeMisc;


    my @variant_seqs = split /,/, $record{alt};
    map {$_ = uc $_} @variant_seqs;
    $featureID = $featureIDStart;

      foreach( @infos ){
        my @temp = split( /=/, $_ );
        if( grep { $_ eq $temp[0] } @wanted ){
          # are there multiple variances?
          my @diffCases = split( /,/, $temp[1] );
          foreach( @diffCases ){
            # see if it has decimal and worthy to be included (i.e !A.XY WHERE XY are not both zero and A is natural num)          
            my @digits =  split ( /\./, $_ );
            my $value;
            if( scalar( @digits ) > 1 ){
              if( $digits[1] > 0 ){                
                $value = $_;   # decimal is worthy, include
                $value =~ s/0{1,}$//; # remove trailing zeroes
              }else{
                $value = $digits[0]; # no decimal, just the whole num / left part of the num in decimal
              }
            }else{
              $value = $_;
            }
            push( @{ $includeMisc{ $temp[0] } },  $value );
            $count++;
          }
       }
       if( $count >= $wantedCount ){
           last;
         }
      }    
=info    
    I opted to add "_I" or "_S" to indicate INDEL or SNP, you might want to modify these parts
    accordingly.
=cut    
    if( length($reference_seq) > 1 or grep { length $_ > 1 } @variant_seqs ){
      $type = "indel";  # It seems to me, Chado and GBrowse is not okay if this is in CAPS!
	      $end = ($start + length($reference_seq)) - 1;
      $countINDEL += 1;       
      $featureID .= ( $separator . "I" . $datasetTag . $countINDEL );
     }else{
       $type = "SNP";
       $countSNP += 1;
       $featureID .= ( $separator . "S" . $datasetTag . $countSNP );
     }

     foreach( @variant_seqs ){   
      $note = "$reference_seq>$_,";
      foreach my $key( sort keys %includeMisc ){
        $note .= " $key=" . shift(@{ $includeMisc{$key} });
      }
      print $OUThandle "$seqid\t$sourceName\t$type\t$start\t$end\t$score\t$strand\t$phase\tID=$featureID;Name=$featureID;Note=$note\n";
     }
  }
  handle_message( 'NOTICE', 'Finished', 'Written to $gff_file .' );
}

sub handle_message {
  my $message;
  my ($level, $code, $sentMsg) = @_;
  $sentMsg ||= "No message.";
  chomp $sentMsg;
  $message .= ( localtime() . ' ' . $sentMsg . "\n" );
  if ($level eq 'FATAL') {
    die join ' : ', ($level, $code, $message);
  }else {
    print STDERR join ' : ', ($level, $code, $message);
  }
}
