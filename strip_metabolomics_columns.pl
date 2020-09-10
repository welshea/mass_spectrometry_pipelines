#!/usr/bin/perl -w

# 2020-09-10:  add check for Excel corruption of values >= 1E7
# 2020-09-08:  rename Spikein Flag header to Potential Spikein Flag
# 2020-08-26:  default to deleting "one-hit wonders" (single pre gap-fill peak)
# 2020-08-26:  attempt to deal with more pos/neg sample name typos
# 2020-08-26:  add Number of non-zero peaks, Percent pre-gap peaks,
#              Percent non-zero peaks columns
# 2020-07-24:  add Non-Spikein Identified Flag column to output
#              begin adding support for El-MAVEN
# 2020-07-13:  add additional command line options for output file names
#              more header conforming (strip .cdf extension)
#              change "assigned" to "identified"
#              no longer output assigned/identified file
#
# 2020-05-19: fix more header conforming
#
# 2020-03-20: detect and insert missing row ID column
#             case insensitive search now for .mzXML, due to case issues
#             strip .mzXML from sample names, as they aren't always there

use Scalar::Util qw(looks_like_number);
use POSIX;

sub is_number
{
    # use what Perl thinks is a number first
    # this is purely for speed, since the more complicated REGEX below should
    #  correctly handle all numeric cases
    if (looks_like_number($_[0]))
    {
        # Perl treats infinities as numbers, Excel does not
        if ($_[0] =~ /^[-+]*inf/)
        {
            return 0;
        }
        
        return 1;
    }

    # optional + or - sign at beginning
    # then require either:
    #  a number followed by optional comma stuff, then optional decimal stuff
    #  mandatory decimal, followed by optional digits
    # then optional exponent stuff
    #
    # Perl cannot handle American comma separators within long numbers.
    # Excel does, so we have to check for it.
    # Excel doesn't handle European dot separators, at least not when it is
    #  set to the US locale (my test environment).  I am going to leave this
    #  unsupported for now.
    #
    if ($_[0] =~ /^([-+]?)([0-9]+(,[0-9]{3,})*\.?[0-9]*|\.[0-9]*)([Ee]([-+]?[0-9]+))?$/)
    {
        # current REGEX can treat '.' as a number, check for that
        if ($_[0] eq '.')
        {
            return 0;
        }
        
        return 1;
    }
    
    return 0;
}

sub reformat_sci
{
    my $field = $_[0];
    my $temp;

    # n.nnEn, nn.nEn, etc.
    #
    # If it looks like scientific notation, Excel will automatically
    #  format it as scientific notation.  However, if the magnitude is
    #  >= 1E7, it will also automatically display it set to only 2 digits
    #  to the right of the decimal (it is still fine internally).  If the
    #  file is re-exported to text, truncation to 3 significant digits
    #  will occur!!!
    #
    # Reformat the number to (mostly) avoid this behavior.
    #
    # Unfortunately, the >= 11 significant digits behavior still
    #  triggers, so it still truncates to 10 digits when re-exporting
    #  General format.  10 digits is still better than 3....
    #
    # The re-export truncation behavior can only be more fully avoided by
    #  manually setting the format to Numeric and specifying a large
    #  number of digits after the decimal place for numbers with
    #  fractions, or 0 digits after the decimal for whole numbers.
    #
    # Ugh.
    #
    # There is no fully fixing this brain-damagedness automatically,
    #  I can only decrease the automatic truncation of significant
    #  digits from 3 to 10 digits :(  Any precision beyond 10 digits
    #  *WILL* be lost on re-export if the format is set to General.
    #
    # NOTE -- We truncate to 16 significant digits by going through
    #         a standard IEEE double precision intermediate.
    #         However, Excel imports numbers as double precision
    #         anyways, so we aren't losing any precision that Excel
    #         wouldn't already be discarding.
    #
    if (is_number($field))
    {
          # strip commas
          $temp = $field;
          $temp =~ s/\,//g;
          
          if (abs($temp) >= 1 &&
              $temp =~ /^([-+]?[0-9]*\.*[0-9]*)[Ee]([-+]?[0-9]+)$/)
          {
#              $number   = $1 + 0;
#              $exponent = $2 + 0;
              
              $temp /= 1;

              # replace original with new scientific notation format
              $field = $temp;
          }
    }
    
    return $field;
}


sub is_heavy_labeled
{
    my $string = $_[0];
    
    if ($string =~ /\b13C[0-9]+\b/) { return 1; }
    if ($string =~ /\bD[0-9]+\b/)   { return 1; }
    
    return 0;
}


$keep_single_pregap_flag = 0;

$syntax_error_flag    = 0;
$num_files            = 0;

for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        if ($field =~ /^--keep-single-pregap$/)
        {
            $keep_single_pregap_flag = 1;
        }
        elsif ($field =~ /^--discard-single-pregap$/)
        {
            $keep_single_pregap_flag = 0;
        }
        else
        {
            printf STDERR "ABORT -- unknown option %s\n", $field;
            $syntax_error_flag = 1;
        }
    }
    else
    {
        if ($num_files == 0)
        {
            $filename = $field;
            $num_files++;
        }
        elsif ($num_files == 1)
        {
            $output_unidentified_filename = $field;
            $num_files++;
        }
        elsif ($num_files == 2)
        {
            $output_spikeins_filename = $field;
            $num_files++;
        }
    }
}


if (!defined($filename) || $syntax_error_flag)
{
    printf STDERR "Usage: strip_mzmine_columns.pl [options] mzmine_tab_delimited.txt [unidentified.txt spikeins.txt]\n";
    printf STDERR "\n";
    printf STDERR "Options:\n";
    printf STDERR "  options which use the MZmine \"row number of detected peaks\" column:\n";
    printf STDERR "    --discard-single-pregap    discard pre gap-filled single-hit rows (default)\n";
    printf STDERR "    --keep-single-pregap       keep pre gap-filled single-hit rows\n";

    exit(1);
}


if (!defined($output_unidentified_filename))
{
    $output_unidentified_filename = 'metabolomics_auto_unidentified.txt';
}
if (!defined($output_spikeins_filename))
{
    $output_spikeins_filename     = 'metabolomics_auto_spikeins.txt';
}


open INFILE, "$filename" or die "can't open $filename\n";


# header line
$line = <INFILE>;
$line =~ s/[\r\n]+//g;
$line =~ s/\"//g;
@array = split /\t/, $line;
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    $array[$i] =~ s/\s+/ /g;

    # clean up sample names
    $array[$i] =~ s/^_+//;
    $array[$i] =~ s/_+$//;
    
    $field = $array[$i];
    $header_col_hash{$field} = $i;
    $header_col_array[$i] = $field;
}


# Excel gets very upset if first field is "ID", change it
if ($header_col_array[0] =~ /^id$/i)
{
    $header_col_array[0] = 'Index';
    $header_col_hash{'Index'} = 0;
}


# conform sample columns
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];

    # strip mzXML from sample names
    $field =~ s/\.mzX?ML( Peak \S+)$/$1/i;
    $field =~ s/\.cdf( Peak \S+)$/$1/i;

    if ($field =~ / Peak (\S+)$/i)
    {

        $second_word = $1;
        $second_word = lc $second_word;
        $field =~ s/ Peak (\S+)$/ Peak $second_word/i;
    }

    $header_col_array[$col] = $field;
}


# scan for peak height and peak area
$peak_height_flag = 0;
$peak_area_flag   = 0;

for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];

    if ($field =~ / Peak height$/i)
    {
        $peak_height_flag = 1;
        next;
    }

    if ($field =~ / Peak area$/i)
    {
        $peak_area_flag = 1;
        next;
    }
}


# flag columns to remove, as they clutter up the file
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];
    
    if ($field =~ /\.mzX?ML[^.]+$/i ||
        $field =~ /\.cdf[^.]+$/i ||
        $field =~ / Peak \S+$/ ||
        $field =~ /row identity/ ||
        $field =~ /row comment/)
    {
        # always keep peak height
        if ($field =~ / Peak height$/i)
        {
            next;
        }
        
        # only keep peak area if no peak height
        if ($peak_height_flag == 0 &&
            $field =~ / Peak area$/i)
        {
            next;
        }
        
        # only keep row identity (all IDs)
        if ($field eq 'row identity (all IDs)')
        {
            next;
        }
 
        $col_to_remove_hash{$col} = 1;
    }
    
    # mzMine exports a blank column at the end of every row,
    # which can majorly screw up other software
    # just remove ALL blank column headers
    if ($field eq '')
    {
        $col_to_remove_hash{$col} = 1;
    }
}


$first_abundance_col = 9E99;
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];

    if ($field =~ / Peak height$/i ||
        $field =~ / Peak area$/i)
    {
        if (!defined($col_to_remove_hash{$col}))
        {
            $sample_col_hash{$col} = 1;
        
            if ($col < $first_abundance_col)
            {
                $first_abundance_col = $col;
            }
        }
    }
}

# none found, check for pos/neg
if ($first_abundance_col == 9E99)
{
  for ($col = 0; $col < @header_col_array; $col++)
  {
    $field = $header_col_array[$col];

    if ($field =~ /(^|[^A-Za-z0-9]+)(pos|neg)([^A-Za-z0-9]+|$)/ ||
        $field =~ /(pos|neg)$/)
    {
        if (!defined($col_to_remove_hash{$col}))
        {
            $sample_col_hash{$col} = 1;

            if ($col < $first_abundance_col)
            {
                $first_abundance_col = $col;
            }
        }
    }
  }
}


@sample_col_array = sort {$a<=>$b} keys %sample_col_hash;
$n_total = @sample_col_array;


# scan for number of detected peaks prior to gap filling
$has_pregap_flag = 0;
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];
    
    if ($field =~ /row number of detected peaks/i)
    {
        $has_pregap_flag = 1;
        
        $pregap_col = $col;
    }
}


# we're going to use the compound names to identify spikeins
# and signal/background peaks
$rowid_col = $header_col_hash{'row ID'};
$name_col  = $header_col_hash{'row identity (all IDs)'};

# maybe it is an El-MAVEN file
if (!defined($rowid_col))
{
    $rowid_col = $header_col_hash{'groupId'};
}
if (!defined($name_col))
{
    $name_col = $header_col_hash{'compound'};
}


#if (!defined($rowid_col))
#{
#    printf STDERR "ABORT -- can't find 'row ID' column\n";
#    exit(1);
#}
if (!defined($name_col))
{
    printf STDERR "ABORT -- can't find 'row identity (all IDs)' column\n";
    exit(1);
}


# print missing rowid
$print_flag = 0;
if (!defined($rowid_col))
{
    print "row ID";
    $print_flag = 1;
}

# print header line
for ($col = 0; $col < @header_col_array; $col++)
{
    if (defined($col_to_remove_hash{$col}))
    {
        next;
    }
    
    if ($print_flag == 1)
    {
        print "\t";
    }
    else
    {
        $print_flag = 1;
    }

    if ($col == $first_abundance_col)
    {
        print "Number of non-zero peaks\t";
        
        if ($has_pregap_flag)
        {
            print "Percent pre-gap peaks\t";
        }
        
        print "Percent non-zero peaks\t";
        print "Potential Spikein Flag\t";
        print "Identified Flag\t";
        print "Non-Spikein Identified Flag\t";
    }

    print $header_col_array[$col];
}
print "\n";



open OUTFILE_UNIDENTIFIED, ">$output_unidentified_filename" or die "can't open file $output_unidentified_filename for writing\n";
open OUTFILE_SPIKEINS,     ">$output_spikeins_filename"     or die "can't open file $output_spikeins_filename for writing\n";


$row_count = 0;
$num_excel_large = 0;	# number of potential Excel-corruptible values
$sum_excel_large = 0;	# sum of number of non-zero last 5 digits

while(defined($line=<INFILE>))
{
    $row_count++;

    $line =~ s/[\r\n]+//g;
    $line =~ s/\"//g;

    @array = split /\t/, $line;

    $print_flag = 0;

    # clean up fields
    for ($col = 0; $col < @array; $col++)
    {
        $array[$col] =~ s/^\s+//;
        $array[$col] =~ s/\s+$//;
        $array[$col] =~ s/\s+/ /g;
        
        $array[$col] = reformat_sci($array[$col]);
    }
    

    if (defined($rowid_col))
    {
        $rowid = $array[$rowid_col];
    }
    else
    {
        $rowid = $row_count;
    }

    
    # flag spikeins, store identified/unidentified
    $name = $array[$name_col];
    $spikein_flag = 0;
    if (is_heavy_labeled($name))
    {
        $spikein_flag = 1;
        
        print OUTFILE_SPIKEINS "$rowid\n";
    }
    
    # consider it identified if it has two letters in a row
    # this should result in treating purely chemical formulas as unidentified
    $identified_flag = 0;
    if ($name =~ /[A-Za-z][A-Za-z]/ &&
        !($name =~ /\d m\/z adduct of \d/) &&
        !($name =~ /Complex of [0-9.]+ and \d+/))
    {
        $identified_flag = 1;
    }
    else
    {
        print OUTFILE_UNIDENTIFIED "$rowid\n";
    }
    
    $nsid_flag = 0;
    if ($spikein_flag == 0 && $identified_flag == 1)
    {
        $nsid_flag = 1;
    }


    # print missing rowid
    if (!defined($rowid_col))
    {
        print $rowid;
        $print_flag = 1;
    }
    
    # count number of non-zero samples
    $n = 0;
    foreach $col (@sample_col_array)
    {
        $value = $array[$col];
        
        if (is_number($value) && $value != 0)
        {
            $n++;
        }
    }
    
    $n_percent = 0;
    if ($n_total)
    {
        $n_percent = 100 * $n / $n_total;
    }

    $n_pregap = 0;
    $n_percent_pregap = 0;    
    if ($has_pregap_flag && $n_total)
    {
        $n_pregap = $array[$pregap_col];
        
        if (!is_number($n_pregap))
        {
            $n_pregap = 0;
        }
    
        $n_percent_pregap = 100 * $n_pregap / $n_total;
    }
    
    # skip "one-hit wonders"
    if ($has_pregap_flag && $keep_single_pregap_flag == 0)
    {
        if ($n_pregap < 2)
        {
            next;
        }
    }

    for ($col = 0; $col < @array; $col++)
    {
        $field = $array[$col];
    
        if (defined($col_to_remove_hash{$col}))
        {
            next;
        }
    
        if ($print_flag == 1)
        {
            print "\t";
        }
        else
        {
            $print_flag = 1;
        }
        
        if ($col == $first_abundance_col)
        {
            print "$n\t";
            
            if ($has_pregap_flag)
            {
                print "$n_percent_pregap\t";
            }
            
            print "$n_percent\t";
            print "$spikein_flag\t";
            print "$identified_flag\t";
            print "$nsid_flag\t";
        }
        
        print $field;
    }
    print "\n";


    # scan samples for number of significant digits
    for ($col = $first_abundance_col; $col < @header_col_array; $col++)
    {
        if (defined($col_to_remove_hash{$col}))
        {
            next;
        }

        $value = $array[$col];
        
        # largeish numbers succeptible to Excel corruption
        if (is_number($value) && $value != 0 &&
            $value >= 1E7)
        {
            if ($value =~ /([0-9]{5})$/)
            {
                $last_five_digits = $1;

                $num_excel_large += 1;

                @temp_array = ($last_five_digits =~ m/[1-9]/g);
                $count = @temp_array;
                
                $sum_excel_large += $count;
            }
        }
    }
}

$avg_excel_large = 9E99;
if ($num_excel_large)
{
    $avg_excel_large = $sum_excel_large / $num_excel_large;
}

if ($avg_excel_large < 3.1)
{
    printf STDERR "WARNING -- !!! Excel-corrupted values >= 1E7 detected !!!\n";
}
