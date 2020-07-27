#!/usr/bin/perl -w

# 2020-07-27:    output metadata in (mostly) original column order
#                make sure newly added flag metadata columns are last
#
# 2020-07-24:    add Non-Spikein Identified Flag column order
#                begin adding support for El-MAVEN
#
# 2020-07-17:    deal with inconsistent case between pos/neg sample headers
#
# 2020-07-13:    deal with changes to expected column header formats
#                more aggressively strip _pos and _neg related stuff


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
    #  > 1E7, it will also automatically display it set to only 2 digits
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


sub read_in_file
{
    my $filename = $_[0];
    my $pos_neg  = $_[1];
    my $line;
    my @array = ();
    my @header_col_array  = ();
    my %header_col_hash   = ();
    my %sample_col_hash   = ();
    my %metadata_col_hash = ();
    my @sample_col_array  = ();
    my @metadata_col_array = ();
    my $sample_col;
    my $metadata_col;
    my $sample;
    my $header;
    my $field;

    open INFILE, "$filename" or die "can't open $filename\n";

    # header line
    $line = <INFILE>;
    $line =~ s/[\r\n]+//g;
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

        if ($field =~ /\S/)
        {        
            $header_col_hash{$field} = $i;
            $header_col_array[$i]    = $field;
            
            $global_concat_header_array[$global_concat_header_count] = $field;
            $global_concat_header_count++;
        }
    }

    # Excel gets very upset if first field is "ID", change it
    if ($header_col_array[0] =~ /^id$/i)
    {
        $header_col_array[0]      = 'Index';
        $header_col_hash{'Index'} = 0;
    }
    
    $row_id_col = $header_col_hash{'row ID'};
    
    # maybe it is an El-MAVEN file
    if (!defined($row_id_col))
    {
        $row_id_col = $header_col_hash{'groupId'};
    }
    
    if (!defined($row_id_col))
    {
        printf STDERR "ABORT -- can't find \'row ID\' column in file %s\n",
            $filename;
        exit(1);
    }
    
    $global_row_id_str = $header_col_array[$row_id_col];

    # categorize columns
    for ($col = 0; $col < @header_col_array; $col++)
    {
        $field = $header_col_array[$col];

        if ($field =~ /^IRON /i ||
            $field =~ /(^|[^A-Za-z0-9]+)(pos|neg)([^A-Za-z0-9]+|$)/i)
        {
            $sample_col_hash{$col} = $field;
        }
        else
        {
            $metadata_col_hash{$col} = $field;
        }
    }
    
    @sample_col_array   = sort {$a<=>$b} keys %sample_col_hash;
    @metadata_col_array = sort {$a<=>$b} keys %metadata_col_hash;

    $all_pos_start_flag = 1;
    $all_neg_start_flag = 1;
    foreach $sample_col (@sample_col_array)
    {
        $sample = $header_col_array[$sample_col];
        $sample =~ s/^IRON //;

        if ($pos_neg eq 'pos')
        {
            # ^pos _pos_ _pos$
            #
            if (!($sample =~ s/(^|[^A-Za-z0-9]+)pos([^A-Za-z0-9]+|$)/$2/i))
            {
                $all_pos_start_flag = 0;
            }
            $all_neg_start_flag = 0;
        }
        elsif ($pos_neg eq 'neg')
        {
            if (!($sample =~ s/(^|[^A-Za-z0-9]+)neg([^A-Za-z0-9]+|$)/$2/i))
            {
                $all_neg_start_flag = 0;
            }
            $all_pos_start_flag = 0;
        }
    }
    
    while(defined($line=<INFILE>))
    {
        $line =~ s/[\r\n]+//g;

        @array = split /\t/, $line;

        for ($i = 0; $i < @array; $i++)
        {
            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;
            $array[$i] =~ s/\s+/ /g;
        }
        
        $row_id = $array[$row_id_col];
        if (is_number($row_id))
        {
            $row_id = sprintf "%s_%05d", $pos_neg, $row_id;
        }
        else
        {
            $row_id = sprintf "%s_%s", $pos_neg, $row_id;
        }
        $array[$row_id_col] = $row_id;
        
        $global_row_id_hash{$row_id} = 1;
        
        # store data
        foreach $sample_col (@sample_col_array)
        {
            $sample = $header_col_array[$sample_col];
            $sample =~ s/^IRON //;
            
            # strip pos/neg from sample name
            if ($all_pos_start_flag)
            {
                $sample =~ s/(^|[^A-Za-z0-9]+)pos([^A-Za-z0-9]+|$)/$2/i;
            }
            elsif ($all_neg_start_flag)
            {
                $sample =~ s/(^|[^A-Za-z0-9]+)neg([^A-Za-z0-9]+|$)/$2/i;
            }

            $sample_lc = lc $sample;
            $sample_lc_to_orig_hash{$sample_lc}{$sample} = 1;
            
            $global_data_hash{$row_id}{$sample_lc} = $array[$sample_col];
            $global_sample_hash{$sample_lc} = 1;
            
            foreach $metadata_col (@metadata_col_array)
            {
                $header = $header_col_array[$metadata_col];
                
                $global_data_hash{$row_id}{$header} = $array[$metadata_col];
                $global_metadata_hash{$header} = 1;
            }
        }
    }
}

$filename_pos = shift;
$filename_neg = shift;

if ($filename_pos =~ /neg/i &&
    !($filename_pos =~ /pos/i))
{
    printf STDERR "ABORT -- check file names to remove pos/neg ambiguity\n";
    printf STDERR "required input order: iron_pos iron_neg\n";
    exit(1);
}
if ($filename_neg =~ /pos/i &&
    !($filename_neg =~ /neg/i))
{
    printf STDERR "ABORT -- check file names to remove pos/neg ambiguity\n";
    printf STDERR "required input order: iron_pos iron_neg\n";
    exit(1);
}
if (!($filename_pos =~ /pos/i) ||
    !($filename_neg =~ /neg/i))
{
    printf STDERR "ABORT -- check file names to remove pos/neg ambiguity\n";
    printf STDERR "required input order: iron_pos iron_neg\n";
    exit(1);
}

@global_concat_header_array = ();
$global_concat_header_count = 0;
$global_row_id_str = '';

read_in_file($filename_pos, 'pos');
read_in_file($filename_neg, 'neg');


# global header order
$index = 0;
@global_header_array = ();
@global_sample_array   = sort keys %global_metadata_hash;
#@global_metadata_array = sort keys %global_metadata_hash;
@global_sample_array   = sort keys %global_sample_hash;
@global_row_id_array   = sort keys %global_row_id_hash;
%seen_header_hash = ();

# put row identifier first, regardless of its original order
foreach $header (@global_concat_header_array)
{
    if ($header eq $global_row_id_str)
    {
        if (!defined($seen_header_hash{$header}))
        {
            $global_header_array[$index]       = $header;
            $global_header_array_print[$index] = $header;
            $index++;
        
            $seen_header_hash{$header} = 1;
        }
    }
}
# add other metadata, skipping flag columns added by the pipeline
foreach $header (@global_concat_header_array)
{
    if (!defined($seen_header_hash{$header}) &&
        defined($global_metadata_hash{$header}))
    {
        if ($header eq 'Spikein Flag' ||
            $header eq 'Identified Flag' ||
            $header eq 'Non-Spikein Identified Flag')
        {
            next;
        }
    
        $global_header_array[$index]       = $header;
        $global_header_array_print[$index] = $header;
        $index++;
        
        $seen_header_hash{$header} = 1;
    }
}
# add in final flag columns, so we can use them to know the rest are samples
foreach $header (@global_concat_header_array)
{
    if (!defined($seen_header_hash{$header}) &&
        defined($global_metadata_hash{$header}))
    {
        if ($header eq 'Spikein Flag' ||
            $header eq 'Identified Flag' ||
            $header eq 'Non-Spikein Identified Flag')
        {
            $global_header_array[$index]       = $header;
            $global_header_array_print[$index] = $header;
            $index++;
        
            $seen_header_hash{$header} = 1;
        }
    }
}


foreach $header (@global_sample_array)
{
    @temp_array = sort keys %{$sample_lc_to_orig_hash{$header}};
    
    $header_chosen = $temp_array[0];
    
    # \p{Uppercase} matching a single Unicode uppercase character
    # [A-Z] won't include unicode uppercase, so we have to get fancier
    $count_uc_best = () = $header_chosen =~ m/\p{Lu}/g;

    # pick conflicting header case by choosing one with more uppercase
    if (@temp_array > 1)
    {
        foreach $header_case (@temp_array)
        {
            $count_uc = () = $header_case =~ m/\p{Uppercase}/g;

#            printf STDERR "Conflicting case:\t%s\t%d\t%d\n",
#                $header_case, $count_uc_best, $count_uc;

            if ($count_uc > $count_uc_best)
            {
                $header_chosen = $header_case;
                $count_uc_best = $count_uc;
            }
        }
    }

    $global_header_array[$index]       = $header;
    $global_header_array_print[$index] = $header_chosen;
    $index++;

    $seen_header_hash{$header} = 1;
}


# print header
for ($i = 0; $i < @global_header_array; $i++)
{
    $header = $global_header_array_print[$i];
    
    if ($header =~ /^IRON /)
    {
        $header =~ s/^IRON\s+//;
    }
    
    if ($i)
    {
        print "\t";
    }
    
    print "$header";
}
print "\n";


# print merged lines
foreach $row_id (@global_row_id_array)
{
    for ($col = 0; $col < @global_header_array; $col++)
    {
        $header = $global_header_array[$col];
        
        $value = $global_data_hash{$row_id}{$header};
        
        if (!defined($value))
        {
            $value = '';
        }
        
#        $value = reformat_sci($value);
        
        if ($col)
        {
            print "\t";
        }
        
        print "$value";
    }
    print "\n";
}
