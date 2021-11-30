#!/usr/bin/perl -w


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


%guess_ion_charge_hash = ();    # crude guess at charge, purely based on +/-

# known ion charges
$known_ion_charge_hash{'+H'}     = 'pos';
$known_ion_charge_hash{'+NH4'}   = 'pos';
$known_ion_charge_hash{'+Na'}    = 'pos';
$known_ion_charge_hash{'+H-H2O'} = 'pos';
$known_ion_charge_hash{'-H'}     = 'neg';
$known_ion_charge_hash{'+HCOO'}  = 'neg';


$num_files = 0;
$filename_pos = 'lipidomics_split_pos.txt';
$filename_neg = 'lipidomics_split_neg.txt';
for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        printf STDERR "ABORT -- unknown option %s\n", $field;
        $syntax_error_flag = 1;
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
            $filename_pos = $field;
            $num_files++;
        }
        elsif ($num_files == 2)
        {
            $filename_neg = $field;
            $num_files++;
        }
    }
}


if (!defined($filename) || $syntax_error_flag)
{
    printf STDERR "Usage: split_lipidomics_into_pos_neg.pl tab_delimited.txt [outfile_pos outfile_neg]\n";
    exit(1);
}


open INFILE,      "$filename"      or die "can't open $filename\n";
open OUTFILE_POS, ">$filename_pos" or die "can't open output $filename_pos\n";
open OUTFILE_NEG, ">$filename_neg" or die "can't open output $filename_neg\n";


# skip down to first line that has anything on it
# lipidomics data has this issues sometimes
while($line=<INFILE>)
{
    # skip comment lines
    if ($line =~ /^#/)
    {
        next;
    }

    # this line isn't purely whitespace, assume it is the header line
    if ($line =~ /\S/)
    {
        last;
    }
}


# header line
$line =~ s/[\r\n]+//g;
$line =~ s/\"//g;
@array = split /\t/, $line;    # skip empty headers at and
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
        $field =~ s/[ ._]+Peak (\S+)$/ Peak $second_word/i;
    }
    $field =~ s/[ ._]+$//;

    $header_col_array[$col] = $field;
}
$num_header_cols = @header_col_array;


# scan for peak height and peak area
$peak_height_flag = 0;
$peak_area_flag   = 0;

for ($col = 0; $col < $num_header_cols; $col++)
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
    
    # lipidomics
    if ($field =~ /^Area,/i)
    {
         $peak_area_flag = 1;
         next;
    }
    if ($field =~ /^Height,/i)
    {
         $peak_height_flag = 1;
         next;
    }
}


# uh oh, no data columns found
# conform them so that Area or Height start at the beginning or end
if ($peak_height_flag == 0 && $peak_area_flag == 0)
{
    printf STDERR "WARNING -- non-standard Height/Area nomenclature, conforming with best guess\n";

    for ($col = 0; $col < $num_header_cols; $col++)
    {
        $field = $header_col_array[$col];

        if ($field =~ /([^A-Za-z0-9]*Height[^A-Za-z0-9]*)/i)
        {
            $height_area_str = $1;
        
            # conform the sample header
            # string may contain () or [], so escape it with \Q \E
            $field =~ s/\Q$height_area_str\E/_/;
            
            # deal with inserted _ at beginning/end
            if ($field =~ /^(_+)/)
            {
                $count_underscores_new  = length $1;
                $count_underscores_orig = 0;
                
                if ($header_col_array[$col] =~ /^(_+)/)
                {
                    $count_underscores_orig = length $1;
                }
                
                if ($count_underscores_new != $count_underscores_orig)
                {
                    $field =~ s/^_//;
                }
            }
            if ($field =~ /(_+)$/)
            {
                $count_underscores_new  = length $1;
                $count_underscores_orig = 0;
                
                if ($header_col_array[$col] =~ /(_+)$/)
                {
                    $count_underscores_orig = length $1;
                }
                
                if ($count_underscores_new != $count_underscores_orig)
                {
                    $field =~ s/_$//;
                }
            }

            # deal with dangling ( or ]
            if ($field =~ /\]$/ && !($field =~ /\[/))
            {
                $field =~ s/\]$//;
            }
            if ($field =~ /\)$/ && !($field =~ /\(/))
            {
                $field =~ s/\)$//;
            }
            
            $field .= ' Peak height';

            $header_col_array[$col] = $field;
            $peak_height_flag = 1;
        }
        if ($field =~ /([^A-Za-z0-9]*Area[^A-Za-z0-9]*)/i)
        {
            $height_area_str = $1;
        
            # conform the sample header
            # string may contain () or [], so escape it with \Q \E
            $field =~ s/\Q$height_area_str\E/_/;

            # deal with inserted _ at beginning/end
            if ($field =~ /^(_+)/)
            {
                $count_underscores_new  = length $1;
                $count_underscores_orig = 0;
                
                if ($header_col_array[$col] =~ /^(_+)/)
                {
                    $count_underscores_orig = length $1;
                }
                
                if ($count_underscores_new != $count_underscores_orig)
                {
                    $field =~ s/^_//;
                }
            }
            if ($field =~ /(_+)$/)
            {
                $count_underscores_new  = length $1;
                $count_underscores_orig = 0;
                
                if ($header_col_array[$col] =~ /(_+)$/)
                {
                    $count_underscores_orig = length $1;
                }
                
                if ($count_underscores_new != $count_underscores_orig)
                {
                    $field =~ s/_$//;
                }
            }

            # deal with dangling ( or ]
            if ($field =~ /\]$/ && !($field =~ /\[/))
            {
                $field =~ s/\]$//;
            }
            if ($field =~ /\)$/ && !($field =~ /\(/))
            {
                $field =~ s/\)$//;
            }

            $field .= ' Peak area';

            $header_col_array[$col] = $field;
            $peak_area_flag = 1;
        }
    }
}



$name_col = $header_col_hash{'LipidIon'};



$first_abundance_col = 9E99;
for ($col = 0; $col < $num_header_cols; $col++)
{
    $field = $header_col_array[$col];

    if ($field =~ / Peak height$/i ||
        $field =~ / Peak area$/i ||
        $field =~ /^Area,/i ||
        $field =~ /^Height,/i)
    {
        $sample_col_hash{$col} = 1;

        if ($col < $first_abundance_col)
        {
            $first_abundance_col = $col;
        }
    }
    else
    {
        if ($field =~ /\S/)
        {
            $metadata_col_hash{$col} = 1;
        }
    }
}


printf STDERR "First abundance col: %s\n", $first_abundance_col;


@sample_col_array = sort {$a<=>$b} keys %sample_col_hash;
#$num_sample_cols  = @sample_col_array;


if (!defined($name_col))
{
    printf STDERR "ABORT -- can't find 'LipidIon' column\n";
    exit(1);
}


#$rowid_col = $header_col_hash{'row ID'};
#if (!defined($rowid_col))
#{
#    printf STDERR "No row identifier column found, creating arbitrary identifiers\n";
#}


@sample_col_array   = sort {$a <=> $b} keys %sample_col_hash;
@metadata_col_array = sort {$a <=> $b} keys %metadata_col_hash;


# read in file, figure out which rows/columns are pos/neg later
$row = 0;
while(defined($line=<INFILE>))
{
    $line =~ s/[\r\n]+//g;
    $line =~ s/\"//g;

    @array = split /\t/, $line, -1;    # don't skip empty fields at and

    # clean up fields
    for ($col = 0; $col < @array; $col++)
    {
        $array[$col] =~ s/^\s+//;
        $array[$col] =~ s/\s+$//;
        $array[$col] =~ s/\s+/ /g;
        
        $array[$col] = reformat_sci($array[$col]);
    }
    
    # store data
    for ($col = 0; $col < $num_header_cols; $col++)
    {
        $field = $array[$col];
        
        if (!defined($field))
        {
            $field = '';
        }
        
        # set missing or non-numeric abundance data to zero
        if (defined($sample_col_hash{$col}) &&
            !is_number($field))
        {
            $field = '0';
        }
        
        $row_data_array[$row][$col] = $field;
    }

    $row++;
}
$num_rows = $row;


foreach $col (@sample_col_array)
{
    $col_counts_hash{$col}{pos} = 0;
    $col_counts_hash{$col}{neg} = 0;
}


# flag pos/neg rows
# count pos/neg for each sample
for ($row = 0; $row < $num_rows; $row++)
{
    $lipid_ion = $row_data_array[$row][$name_col];

    $pos_neg = '';
    $ion_str = '';
    
    if ($lipid_ion =~ /(\+[A-Za-z0-9-]+)$/)
    {
        $pos_neg = 'pos';
        $ion_str = $1;

        $pos_neg_row_hash{$row} = $pos_neg;
    }
    elsif ($lipid_ion =~ /(\-[A-Za-z0-9-]+)$/)
    {
        $pos_neg = 'neg';
        $ion_str = $1;

        $pos_neg_row_hash{$row} = $pos_neg;
    }
    
    # overwrite with known pos/neg when possible
    if ($ion_str ne '')
    {
        if (defined($known_ion_charge_hash{$ion_str}))
        {
            $pos_neg = $known_ion_charge_hash{$ion_str};
            $pos_neg_row_hash{$row} = $pos_neg;
        }
        else
        {
            if (!defined($guess_ion_charge_hash{$ion_str}))
            {
                $guess_ion_charge_hash{$ion_str} = $pos_neg;
                
                printf STDERR "WARNING -- guessing unknown ion charge:\t%s\t%s\n",
                    $ion_str, $pos_neg;
            }
        }
    }
    
    # count pos/neg for each sample
    # so that we can check for un-merged pos/neg columns later
    if ($pos_neg eq 'pos' || $pos_neg eq 'neg')
    {
        foreach $col (@sample_col_array)
        {
            $field = $row_data_array[$row][$col];
            
            if (is_number($field) && $field != 0)
            {
                $col_counts_hash{$col}{$pos_neg} += 1;
            }
        }
    }
}


%pos_col_hash = ();
$neg_col_hash = ();
foreach $col (@sample_col_array)
{
    $num_pos = $col_counts_hash{$col}{pos};
    $num_neg = $col_counts_hash{$col}{neg};
    
    if ($num_pos && $num_neg == 0)
    {
        $pos_col_hash{$col} = 1;
    }
    if ($num_neg && $num_pos == 0)
    {
        $neg_col_hash{$col} = 1;
    }
}

@pos_col_array = sort {$a <=> $b} keys %pos_col_hash;
@neg_col_array = sort {$a <=> $b} keys %neg_col_hash;

$num_pos_cols = @pos_col_array;
$num_neg_cols = @neg_col_array;

$merge_flag = 0;
if ($num_pos_cols && $num_neg_cols)
{
    printf STDERR "WARNING -- unmerged POS/NEG columns\n";
    
    if ($num_pos_cols != $num_neg_cols)
    {
        printf STDERR "#POS/#NEG columns differ, leave unmerged:\t%d/%d\n",
            $num_pos_cols, $num_neg_cols;
    }
    else
    {
        $merge_flag = 1;

        printf STDERR "WARNING -- assuming same input column order for pos/neg\n";
        
        if ($peak_height_flag)
        {
            for ($i = 0; $i < @pos_col_array; $i++)
            {
                $col_pos = $pos_col_array[$i];
                $col_neg = $neg_col_array[$i];
            
                $header_pos = $header_col_array[$col_pos];
                $header_neg = $header_col_array[$col_neg];
                
                $name_pos   = '';
                $metric_pos = '';
                $name_neg   = '';
                $metric_neg = '';
                
                if ($header_pos =~ /(.*?)\s*(Peak height|Peak area)$/)
                {
                    $name_pos   = $1;
                    $metric_pos = $2;
                }
                if ($header_neg =~ /(.*?)\s*(Peak height|Peak area)$/)
                {
                    $name_neg   = $1;
                    $metric_neg = $2;
                }
                
                # remove any pos/neg that may already be in the sample names
                # assume existing pos/neg are correctly labeled as such
                #
                # ^pos _pos_ _pos$
                # ^neg _neg_ _neg$
                #
                if (!($name_pos =~ s/(^|[^A-Za-z0-9]+)pos([^A-Za-z0-9]+|$)/$2/i))
                {
                    $name_pos =~ s/pos$//i;
                }
                if (!($name_neg =~ s/(^|[^A-Za-z0-9]+)neg([^A-Za-z0-9]+|$)/$2/i))
                {
                    $name_neg =~ s/neg$//i;
                }
                
                if ($metric_pos ne '' &&
                    lc($metric_pos) eq lc($metric_neg))
                {
                    # in CICPT #3501, NEG has the smaller -# postfix,
                    # so use NEG first, then POS
                    #
                    # 'neg' also comes before 'pos' alphabetically
                    #
                    $name_merged = $name_neg;
                    if ($name_neg ne $name_pos)
                    {
                        $name_merged = sprintf "%s__%s", $name_neg, $name_pos;
                    }
                    
                    # store col pair assignments for merging later
                    $paired_pos_neg_hash{$col_pos}{name}   = $name_merged;
                    $paired_pos_neg_hash{$col_pos}{metric} = $metric_neg;
                    $paired_pos_neg_hash{$col_pos}{m_col}  = $col_neg;
                    $paired_pos_neg_hash{$col_neg}{name}   = $name_merged;
                    $paired_pos_neg_hash{$col_neg}{metric} = $metric_neg;
                    $paired_pos_neg_hash{$col_neg}{m_col}  = $col_pos;
                    
                    printf STDERR "PAIRING:\t%d\t%d\t%s\t%s\n",
                        $col_neg, $col_pos,
                        $metric_neg,
                        $name_merged;
                }
            }
        }
    }
}



# print new header lines
%skip_paired_col_hash = ();
$header_pos = '';
$header_neg = '';
foreach $col (@metadata_col_array)
{
    $header = $header_col_array[$col];
    
    if ($header_pos ne '') { $header_pos .= "\t"; }
    if ($header_neg ne '') { $header_neg .= "\t"; }
    
    $header_pos .= $header;
    $header_neg .= $header;
}
foreach $col (@sample_col_array)
{
    $header = $header_col_array[$col];

    # skip paired columns we've already dealt with
    if (defined($skip_paired_col_hash{$col}))
    {
        next;
    }

    # merge paired columns into a single header
    if (defined($paired_pos_neg_hash{$col}))
    {
        $name_merged = $paired_pos_neg_hash{$col}{name};
        $metric      = $paired_pos_neg_hash{$col}{metric};

        $header = $name_merged . ' ' . $metric;

        # skip matched column later, since we're dealing with it here
        $col_match = $paired_pos_neg_hash{$col}{m_col};
        $skip_paired_col_hash{$col_match} = 1;
    }
    
    if ($header_pos ne '') { $header_pos .= "\t"; }
    if ($header_neg ne '') { $header_neg .= "\t"; }
    
    $header_pos .= 'pos_' . $header;
    $header_neg .= 'neg_' . $header;
}
print OUTFILE_POS "$header_pos\n";
print OUTFILE_NEG "$header_neg\n";



# output data, split into pos and neg files
$line_pos = '';
$line_neg = '';
for ($row = 0; $row < $num_rows; $row++)
{
    $line_pos = '';
    $line_neg = '';
    foreach $col (@metadata_col_array)
    {
        $field = $row_data_array[$row][$col];
    
        if ($line_pos ne '') { $line_pos .= "\t"; }
        if ($line_neg ne '') { $line_neg .= "\t"; }
    
        $line_pos .= $field;
        $line_neg .= $field;
    }
    foreach $col (@sample_col_array)
    {
        $field = $row_data_array[$row][$col];

        # skip paired columns we've already dealt with
        if (defined($skip_paired_col_hash{$col}))
        {
            next;
        }

        # merge paired columns into a single field
        if (defined($paired_pos_neg_hash{$col}))
        {
            $col_match = $paired_pos_neg_hash{$col}{m_col};
            
            # if current column is missing/zero, use the matched col instead
            if ($field == 0)
            {
                $field = $row_data_array[$row][$col_match];
            }

            # skip matched column later, since we're dealing with it here
            $skip_paired_col_hash{$col_match} = 1;
        }
    
        if ($line_pos ne '') { $line_pos .= "\t"; }
        if ($line_neg ne '') { $line_neg .= "\t"; }
    
        $line_pos .= $field;
        $line_neg .= $field;
    }
    
    # output positive ions to positive file
    # output negative ions to negative file
    # skip uncertain charge lines
    #
    $pos_neg = $pos_neg_row_hash{$row};
    if (!defined($pos_neg))
    {
        $pos_neg = '';
    }
    if ($pos_neg eq 'pos')
    {
        print OUTFILE_POS "$line_pos\n";
    }
    elsif ($pos_neg eq 'neg')
    {
        print OUTFILE_NEG "$line_neg\n";
    }
    else
    {
        $lipid_ion = $row_data_array[$row][$name_col];

        printf STDERR "Skipping unknown ion charge:\t%s\n", $lipid_ion;
    }
}



close OUTFILE_POS;
close OUTFILE_NEG;
