#!/usr/bin/perl -w

# 2022-06-16:  better handling of () and {} near pos/neg in sample names
# 2022-06-16:  better handling of [] near pos/neg in sample names
# 2022-06-16:  bugfix: stop attempting to interpret sample names as REGEX
# 2022-06-07:  add support for single input file
# 2022-04-11:  comment out metabolmoics_qc.pl fix, fixed in metabolomics_qc.pl
# 2022-04-11:  fix missing pos/neg sample names to work with metabolomics_qc.pl
# 2022-03-09:  change debug filter mean message to print pos/neg fractional rows
# 2022-03-09:  use "Non-heavy identified flag" column to ignore points
#              during pos/neg scaling training
# 2022-03-04:  globally shift all pos/neg so that pos/neg have equal means;
#              we may want this for main ion filtering of lipidomics data
# 2022-02-23:  handle merging sample-specific pos/neg metadata separately
# 2022-02-22:  fill in missing pos/neg in sample names when 100% missing
# 2021-10-19:  bugfix auto-shortened sample names
# 2021-08-19:  better support for samples missing pos/neg in their names
# 2021-08-11:  support auto-shortening when blank samples break naming convention
# 2021-08-06:  export auto-shortened sample names to sample name mapping table
# 2021-08-06:  export sample name mapping, warn on unmatched / duped pairs
# 2021-04-06:  check for and uniquify non-unique row IDs
# 2021-04-05:  fix El-MAVEN support, groupID is *NOT* a unique row identifier
# 2021-03-30:  don't convert ' ' to '_' in " Peak height" and " Peak area"
# 2021-03-30:  strip .raw from end of sample names
# 2021-01-06:  handle recent changes to strip_metabolomics_columns.pl headers
# 2020-09-25:  force row identifier header to always print as "row ID"
# 2020-09-25:  minor code clean up
# 2020-09-24:  clean up extra/dangling hyphens in sample names
# 2020-09-10:  clean up extra/dangling underscores and spaces in sample names
# 2020-09-08:  rename Spikein Flag header to Potential Spikein Flag
# 2020-08-26:  attempt to deal with more pos/neg sample name typos
# 2020-08-26:  add Number of non-zero peaks, Percent pre-gap peaks,
#              Percent non-zero peaks columns
# 2020-07-27:  output metadata in (mostly) original column order
#              make sure newly added flag metadata columns are last
#
# 2020-07-24:  add Non-Spikein Identified Flag column order
#              begin adding support for El-MAVEN
#
# 2020-07-17:  deal with inconsistent case between pos/neg sample headers
#
# 2020-07-13:  deal with changes to expected column header formats
#              more aggressively strip _pos and _neg related stuff


use Scalar::Util qw(looks_like_number);
use File::Basename;
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
    my %cols_to_merge_hash   = ();
    my %metadata_col_hash = ();
    my %actual_sample_col_hash = ();
    my @tomerge_col_array  = ();
    my @metadata_col_array = ();
    my $tomerge_col;
    my $metadata_col;
    my $tomerge;
    my $header;
    my $field;
    my $count_non_pos = 0;
    my $count_non_neg = 0;

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
    
    $row_id_col    = $header_col_hash{'row ID'};
    $nhid_flag_col = $header_col_hash{'Non-heavy identified flag'};
    
    if (!defined($row_id_col))
    {
        printf STDERR "No row identifier column found, creating arbitrary identifiers\n";
    }

    if (defined($row_id_col))
    {    
        $global_row_id_str = $header_col_array[$row_id_col];
    }
    else
    {
        $global_row_id_str = 'row ID';
    }

    # categorize columns
    for ($col = 0; $col < @header_col_array; $col++)
    {
        $field = $header_col_array[$col];

        if ($field =~ /^IRON /i ||
            $field =~ /(^|[^A-Za-z0-9]+)(pos|neg)([^A-Za-z0-9]+|$)/i ||
            $field =~ /(pos|neg)$/i)
        {
            $cols_to_merge_hash{$col} = $field;

            # keep track of which columns should be actual sample data
            if ($field =~ /^IRON /i ||
                $field =~ /([^A-Za-z0-9]*(Height|Area)[^A-Za-z0-9]*)/i)
            {
                $actual_sample_col_hash{$col} = 1;
            }
        }
        else
        {
            $metadata_col_hash{$col} = $field;
        }
    }
    
    @tomerge_col_array   = sort {$a<=>$b} keys %cols_to_merge_hash;
    @metadata_col_array = sort {$a<=>$b} keys %metadata_col_hash;

    $all_pos_start_flag = 1;
    $all_neg_start_flag = 1;
    foreach $tomerge_col (@tomerge_col_array)
    {
        $tomerge = $header_col_array[$tomerge_col];
        $tomerge =~ s/^IRON //;

        # ^pos _pos_ _pos$
        #
        if (!($tomerge =~ s/(^|[^\]\)\}A-Za-z0-9]+)pos([^\[\(\{A-Za-z0-9]+|$)/$2/i ||
              $tomerge =~ s/pos$//i))
        {
            $all_pos_start_flag = 0;
            $count_non_pos++;
        }

        if (!($tomerge =~ s/(^|[^\]\)\}A-Za-z0-9]+)neg([^\[\(\{A-Za-z0-9]+|$)/$2/i ||
              $tomerge =~ s/neg$//i))
        {
            $all_neg_start_flag = 0;
            $count_non_neg++;
        }
    }

    
    # uh oh, we are missing pos/neg in all samples
    # assume that we gave the files in the correct order,
    # prepend pos/neg to the ends
    if ($all_pos_start_flag == 0 &&
        $all_neg_start_flag == 0 &&
        $count_non_pos == @tomerge_col_array &&
        $count_non_neg == @tomerge_col_array)
    {
        printf STDERR "WARNING -- missing pos/neg from sample names, prepending %s\n",
            $pos_neg;

        if ($pos_neg =~ /pos/)
        {
            $prepended_pos_flag = 1;
        }
        if ($pos_neg =~ /neg/)
        {
            $prepended_neg_flag = 1;
        }

        foreach $tomerge_col (@tomerge_col_array)
        {
            $header_col_array[$tomerge_col] = sprintf "%s_%s",
                $pos_neg, $header_col_array[$tomerge_col];
        }
        
        if ($pos_neg eq 'pos')
        {
            $all_pos_start_flag = 1;
        }
        if ($pos_neg eq 'neg')
        {
            $all_neg_start_flag = 1;
        }
    }
    
    
    $row = 0;
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

        # row ID exists
        if (defined($row_id_col))
        {        
            $row_id = $array[$row_id_col];
        }
        # create our own
        else
        {
            $row_id = $row + 1;
        }
        
        if (is_number($row_id))
        {
            $row_id = sprintf "%s_%05d", $pos_neg, $row_id;
        }
        else
        {
            $row_id = sprintf "%s_%s", $pos_neg, $row_id;
        }

        if (defined($global_row_id_hash{$row_id}))
        {
            print STDERR "WARNING -- row ID indentifiers not unique, appending row number\n";
            
            $row_id .= sprintf "_%05d", $row + 1;
        }

        # replace original row ID, if it exists, with new one
        if (defined($row_id_col))
        {
            $array[$row_id_col] = $row_id;
        }
        
        $global_row_id_hash{$row_id} = 1;
        
        # store data
        foreach $tomerge_col (@tomerge_col_array)
        {
            $tomerge          = $header_col_array[$tomerge_col];
            $tomerge          =~ s/^((?:pos|neg)[^A-Za-z0-9]*)IRON /$1/;
            $tomerge          =~ s/^IRON //;
            $tomerge_origname = $tomerge;
            
            # strip pos/neg from sample name
            if ($all_pos_start_flag)
            {
                $tomerge =~ s/(^|[^\]\)\}A-Za-z0-9]+)pos([^\[\(\{A-Za-z0-9]+|$)/$2/i;
                $tomerge =~ s/pos$//i;

                # clean up underscores, etc.
                $tomerge =~ s/[_ ]+/_/g;
                $tomerge =~ s/\-+/\-/g;
                $tomerge =~ s/^[_ -]//;
                $tomerge =~ s/[_ -]$//;

            }
            elsif ($all_neg_start_flag)
            {
                $tomerge =~ s/(^|[^\]\)\}A-Za-z0-9]+)neg([^\[\(\{A-Za-z0-9]+|$)/$2/i;
                $tomerge =~ s/neg$//i;

                # clean up underscores, etc.
                $tomerge =~ s/[_ ]+/_/g;
                $tomerge =~ s/\-+/\-/g;
                $tomerge =~ s/^[_ -]//;
                $tomerge =~ s/[_ -]$//;
            }

            # strip .raw from end of sample name
            $tomerge =~ s/\.raw(?=(\s|[_ ]|$))//;
            
            # restore spaces in Peak height/area after conforming
            $tomerge =~ s/_Peak_(height|area)/ Peak $1/ig;


            $tomerge_lc = lc $tomerge;
            $tomerge_lc_to_origcase_hash{$tomerge_lc}{$tomerge} = 1;
            
            if ($pos_neg eq 'pos')
            {
                $tomerge_lc_to_origpos_hash{$tomerge_lc}{$tomerge_origname} = 1;
            }
            elsif ($pos_neg eq 'neg')
            {
                $tomerge_lc_to_origneg_hash{$tomerge_lc}{$tomerge_origname} = 1;
            }
            
            $global_data_hash{$row_id}{$tomerge_lc} = $array[$tomerge_col];
            $global_tomerge_hash{$tomerge_lc} = 1;
            
            # keep track of which headers are actual sample data, rather
            # than sample-related metadata
            if (defined($actual_sample_col_hash{$tomerge_col}))
            {
                $actual_sample_lc_hash{$tomerge_lc} = 1;
                
                # check to see if all data is logged or not
                # any value > 100 is unlikely to have been logged
                $value = $array[$tomerge_col];
                if (defined($value) && is_number($value) &&
                    $value > 100)
                {
                    $all_logged_flag = 0;
                }
            }
            
            foreach $metadata_col (@metadata_col_array)
            {
                $header = $header_col_array[$metadata_col];
                
                $global_data_hash{$row_id}{$header} = $array[$metadata_col];
                $global_metadata_hash{$header} = 1;
            }
            
            $row_id_pos_neg_hash{$row_id} = $pos_neg;
        }
        
        $row++;
    }
}



# begin main()

$log_of_two       = log(2.0);
$filename_pos     = '';
$filename_neg     = '';
$outname_user     = '';    # optional, user-provided sample table filename
$equal_means_flag = 1;     # shift pos/neg log2 means to be equal

$syntax_error_flag = 0;
$num_files         = 0;
for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        if ($field =~ /^--equal-means/)
        {
            $equal_means_flag = 1;
        }
        elsif ($field =~ /^--no-equal-means/)
        {
            $equal_means_flag = 0;
        }
        else
        {
            printf "ABORT -- unknown option %s\n", $field;
            $syntax_error_flag = 1;
        }
    }
    else
    {
        if ($num_files == 0)
        {
            $filename_pos = $field;
            $num_files++;
        }
        elsif ($num_files == 1)
        {
            $filename_neg = $field;
            $num_files++;
        }
        elsif ($num_files == 2)
        {
            $outname_user = $field;
            $num_files++;
        }
    }
}


if ($syntax_error_flag || $num_files < 2)
{
    $program_name = basename($0);

    printf STDERR "Usage: $program_name [options] pos_data.txt neg_data.txt [output_name.txt]\n";
    printf STDERR "\n";
    printf STDERR "  Options:\n";
    printf STDERR "    --equal-means     scale pos/neg means to be equal (default)\n";
    printf STDERR "    --no-equal-means  disable pos/neg mean scaling\n";

    exit(1);
}


$single_file_mode = '';
if ($filename_pos eq 'NULL' || $filename_neg eq 'NULL')
{
    if ($filename_pos ne 'NULL')
    {
        $single_file_mode = 'pos';
    }
    if ($filename_neg ne 'NULL')
    {
        $single_file_mode = 'neg';
    }
    
    # ABORT -- they both can't be NULL
    if ($filename_pos eq $filename_neg)
    {
        printf STDERR "ABORT -- at least one input filename must be not \"NULL\"\n";
        exit(2);
    }
}


if ($filename_pos =~ /neg/i &&
    !($filename_pos =~ /pos/i))
{
    printf STDERR "WARNING -- pos/neg file names may be swapped\n";
}
if ($filename_neg =~ /pos/i &&
    !($filename_neg =~ /neg/i))
{
    printf STDERR "WARNING -- pos/neg file names may be swapped\n";
}
if (($single_file_mode ne 'neg' && !($filename_pos =~ /pos/i)) ||
    ($single_file_mode ne 'pos' && !($filename_neg =~ /neg/i)))
{
    printf STDERR "WARNING -- pos/neg file names do not contain pos/neg\n";
}

@global_concat_header_array = ();
$global_concat_header_count = 0;
$global_row_id_str = '';

$all_logged_flag = 1;
$prepended_pos_flag = 0;
$prepended_neg_flag = 0;

if ($single_file_mode ne 'neg')
{
    read_in_file($filename_pos, 'pos');
}
if ($single_file_mode ne 'pos')
{
    read_in_file($filename_neg, 'neg');
}

@global_header_array     = ();
#@global_metadata_array  = sort keys %global_metadata_hash;
@global_tomerge_array    = sort keys %global_tomerge_hash;
@global_row_id_array     = sort keys %global_row_id_hash;
@actual_sample_lc_array  = sort keys %actual_sample_lc_hash;
$num_actual_samples      = @actual_sample_lc_array;


# global header order
$index = 0;
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
        # current column headers
        if ($header =~ /^Number of non-zero peaks$/i ||
            $header =~ /^Percent pre-gap peaks$/i    ||
            $header =~ /^Percent non-zero peaks$/i   ||
            $header =~ /^Heavy-labeled flag$/i       ||
            $header =~ /^Identified flag$/i          ||
            $header =~ /^Non-heavy identified flag$/i)

        {
            next;
        }

        # older, deprecated column headers
        if ($header =~ /^Potential Spikein Flag$/i   ||
            $header =~ /^Non-Spikein Identified Flag$/i)
        {
            next;
        }
    
        $global_header_array[$index]       = $header;
        $global_header_array_print[$index] = $header;
        $index++;
        
        $seen_header_hash{$header} = 1;
    }
}

# sample-related metadata, such as Rt[] columns
foreach $header (@global_tomerge_array)
{
    # skip actual sample column, we'll handle those later
    if (defined($actual_sample_lc_hash{$header}))
    {
        next;
    }

    @temp_array = sort keys %{$tomerge_lc_to_origcase_hash{$header}};
    
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
    $chosen_header_case_hash{$header}  = $header_chosen;
    $index++;

    $seen_header_hash{$header} = 1;
}


# add in final flag columns, so we can use them to know the rest are samples
foreach $header (@global_concat_header_array)
{
    if (!defined($seen_header_hash{$header}) &&
        defined($global_metadata_hash{$header}))
    {
        if ($header =~ /^Number of non-zero peaks$/i  ||
            $header =~ /^Percent pre-gap peaks$/i     ||
            $header =~ /^Percent non-zero peaks$/i    ||
            $header =~ /^Heavy-labeled flag$/i        ||
            $header =~ /^Potential Spikein Flag$/i    ||
            $header =~ /^Identified flag$/i           ||
            $header =~ /^Non-heavy identified flag$/i ||
            $header =~ /^Non-Spikein Identified Flag$/i)
        {
            $global_header_array[$index]       = $header;
            $global_header_array_print[$index] = $header;
            $index++;
        
            $seen_header_hash{$header} = 1;
        }
    }
}


# all the actual sample data
foreach $header (@actual_sample_lc_array)
{
    @temp_array = sort keys %{$tomerge_lc_to_origcase_hash{$header}};
    
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
    $chosen_header_case_hash{$header}  = $header_chosen;
    $index++;

    $seen_header_hash{$header} = 1;
}


# only shift pos/neg to be equal means if data is logged and both present
$pos_shift = 0;
$neg_shift = 0;
if ($equal_means_flag && $single_file_mode eq '')
{
    foreach $row_id (@global_row_id_array)
    {
        # use Non-heavy identified flag to exclude rows from training
        if (defined($nhid_flag_col))
        {
            $nhid_flag    =
                $global_data_hash{$row_id}{'Non-heavy identified flag'};

            $use_row_flag = 0;
            if (is_number($nhid_flag) && $nhid_flag >= 1)
            {
                $use_row_flag = 1;
            }

            if ($use_row_flag == 0)
            {
                next;
            }
        }
    
        $pos_neg = $row_id_pos_neg_hash{$row_id};
    
        for ($col = 0; $col < @global_header_array; $col++)
        {
            $header = $global_header_array[$col];
            $value  = $global_data_hash{$row_id}{$header};

            # sum the pos/neg log2 values over all data rows
            if (defined($actual_sample_lc_hash{$header}) &&
                defined($value) && is_number($value) && $value > 0)
            {
                if (!defined($pos_neg_stats_hash{$pos_neg}))
                {
                    $pos_neg_stats_hash{$pos_neg}{sum}   = 0;
                    $pos_neg_stats_hash{$pos_neg}{count} = 0;
                }

                # data already logged
                if ($all_logged_flag)
                {
                    $pos_neg_stats_hash{$pos_neg}{sum} += $value;
                }
                # log2 transform prior to summing
                else
                {
                    $pos_neg_stats_hash{$pos_neg}{sum} += log($value) /
                                                          $log_of_two;
                }

                $pos_neg_stats_hash{$pos_neg}{count} += 1;
            }
        }
    }

    # pos/neg stats
    $pos_avg = 0;
    $neg_avg = 0;
    if (defined($pos_neg_stats_hash{pos}))
    {
        $pos_avg = $pos_neg_stats_hash{pos}{sum} /
                   $pos_neg_stats_hash{pos}{count};
    }
    if (defined($pos_neg_stats_hash{neg}))
    {
        $neg_avg = $pos_neg_stats_hash{neg}{sum} /
                   $pos_neg_stats_hash{neg}{count};
    }

    $pos_neg_avg_avg = 0.5 * ($pos_avg + $neg_avg);
    $pos_shift       = $pos_neg_avg_avg - $pos_avg;
    $neg_shift       = $pos_neg_avg_avg - $neg_avg;
    $pos_scale       = 2 ** $pos_shift;
    $neg_scale       = 2 ** $neg_shift;

    if (defined($nhid_flag_col))
    {
        printf STDERR "Using 'Non-heavy identified flag' to filter mean calculations:  %f  %f\n",
            $pos_neg_stats_hash{pos}{count} / $num_actual_samples,
            $pos_neg_stats_hash{neg}{count} / $num_actual_samples;
    }

    if ($all_logged_flag)
    {
        printf STDERR "Shifting log pos/neg mean values to be equal:\t%f\t%f\n",
            $pos_shift, $neg_shift;
    }
    else
    {
        printf STDERR "Scaling pos/neg geometric mean values to be equal:\t%f\t%f\n",
            $pos_shift, $neg_shift;
    }
}


# match sample pairs
foreach $header (@global_tomerge_array)
{
    @temp_array = sort keys %{$tomerge_lc_to_origcase_hash{$header}};

    @pos_tomerge_matches = sort keys %{$tomerge_lc_to_origpos_hash{$header}};
    @neg_tomerge_matches = sort keys %{$tomerge_lc_to_origneg_hash{$header}};
    $count_pos_matches  = @pos_tomerge_matches;
    $count_neg_matches  = @neg_tomerge_matches;
    
    if ($count_pos_matches > 1)
    {
        $dupe_str = join "\t", @pos_tomerge_matches;

        printf STDERR "WARNING -- duplicate tomerges:\t%s\n", $dupe_str;
    }
    if ($count_neg_matches > 1)
    {
        $dupe_str = join "\t", @neg_tomerge_matches;

        printf STDERR "WARNING -- duplicate tomerges:\t%s\n", $dupe_str;
    }
    
    # warn on unmatched pairs
    if ($single_file_mode eq '')
    {
        foreach $pos_tomerge (@pos_tomerge_matches)
        {
            if ($count_neg_matches == 0)
            {
                printf STDERR "WARNING -- unmatched pos/neg tomerge:\t%s\n",
                    $pos_tomerge;
            }
        }
        foreach $neg_tomerge (@neg_tomerge_matches)
        {
            if ($count_pos_matches == 0)
            {
                printf STDERR "WARNING -- unmatched pos/neg tomerge:\t%s\n",
                    $neg_tomerge;
            }
        }
    }
}



# guess at common prefix and/or suffix
%common_prefix_hash = ();
%common_suffix_hash = ();
for ($i = 0; $i < @actual_sample_lc_array; $i++)
{
    $tomerge_1     = $actual_sample_lc_array[$i];
    $tomerge_1_rev = reverse $tomerge_1;

    for ($j = $i+1; $j < @actual_sample_lc_array; $j++)
    {
        $tomerge_2     = $actual_sample_lc_array[$j];
        $tomerge_2_rev = reverse $tomerge_2;
        
        # find common prefix, Perl "dark magic"
        $xor = "$tomerge_1" ^ "$tomerge_2";
        $xor =~ /^\0*/;
        $len = $+[0];
        $common_prefix = substr $tomerge_1, 0, $len;

        # find common suffix, Perl "dark magic"
        $xor = "$tomerge_1_rev" ^ "$tomerge_2_rev";
        $xor =~ /^\0*/;
        $len = $+[0];
        $common_suffix = substr $tomerge_1_rev, 0, $len;
        $common_suffix = reverse $common_suffix;
        
        # strip them to begin/end with a non- alphanumeric character
        $common_prefix =~ s/([^A-Za-z0-9]+)[A-Za-z0-9]+$/$1/;
        $common_suffix =~ s/^[A-Za-z0-9]+([^A-Za-z0-9]+)/$1/;
        
        # blank them if they still don't begin/end in a non- alphanumeric
        if (!($common_prefix =~ /[^A-Za-z0-9]$/))
        {
            $common_prefix = '';
        }
        if (!($common_suffix =~ /^[^A-Za-z0-9]/))
        {
            $common_suffix = '';
        }
        
        if ($common_prefix ne '')
        {
            $common_prefix_hash{$common_prefix} = 1;
        }
        if ($common_suffix ne '')
        {
            $common_suffix_hash{$common_suffix} = 1;
        }
    }
}

@common_prefix_array = sort keys %common_prefix_hash;
@common_suffix_array = sort keys %common_suffix_hash;
$count_tomerge_total = @actual_sample_lc_array;
$common_prefix_all   = '';
$common_suffix_all   = '';

$count_blank = 0;
foreach $tomerge (@actual_sample_lc_array)
{
    if ($tomerge =~ /processing_bla?n?k\d*([^A-Za-z0-9]|$)/i ||
        $tomerge =~ /(^|[^A-Za-z0-9])blank\d*([^A-Za-z0-9]|$)/i)
    {
        $count_blank++;
    }
}

foreach $common_prefix (@common_prefix_array)
{
    $count = 0;

    foreach $tomerge (@actual_sample_lc_array)
    {
        if ($tomerge =~ /^\Q$common_prefix\E/)
        {
            $count++;
        }
    }
    
    if ($count >= $count_tomerge_total - $count_blank)
    {
        $common_prefix_all = $common_prefix;
        last;
    }
}

foreach $common_suffix (@common_suffix_array)
{
    $count = 0;

    foreach $tomerge (@actual_sample_lc_array)
    {
        if ($tomerge =~ /\Q$common_suffix\E$/)
        {
            $count++;
        }
    }
    
    if ($count >= $count_tomerge_total - $count_blank)
    {
        $common_suffix_all = $common_suffix;
        last;
    }
}


# output sample mapping table
# take first of duplicate samples (should never occur)
#
if (defined($outname_user) && $outname_user =~ /\S/)
{
    $outfile_name = $outname_user;
}
else
{
    $outfile_name = 'pipeline_metabolomics_sample_table.txt';
}

open OUTFILE, ">$outfile_name" or die "ABORT -- cannot open $outfile_name for writing\n";

printf OUTFILE "%s",   'Sample';
printf OUTFILE "\t%s", 'SamplePOS';
printf OUTFILE "\t%s", 'SampleNEG';
printf OUTFILE "\t%s", 'SampleAutoShortened';
printf OUTFILE "\n";


foreach $header (@actual_sample_lc_array)
{
    @temp_array = sort keys %{$tomerge_lc_to_origcase_hash{$header}};

    @pos_tomerge_matches = sort keys %{$tomerge_lc_to_origpos_hash{$header}};
    @neg_tomerge_matches = sort keys %{$tomerge_lc_to_origneg_hash{$header}};
    $count_pos_matches  = @pos_tomerge_matches;
    $count_neg_matches  = @neg_tomerge_matches;

    $tomerge_case = $chosen_header_case_hash{$header};
    $pos_tomerge  = '';
    $neg_tomerge  = '';

    $tomerge_short = $tomerge_case;
    $tomerge_short =~ s/^\Q$common_prefix_all\E//i;
    $tomerge_short =~ s/\Q$common_suffix_all\E$//i;
    
    if ($count_pos_matches)
    {
        $pos_tomerge = $pos_tomerge_matches[0];
    }
    if ($count_neg_matches)
    {
        $neg_tomerge = $neg_tomerge_matches[0];
    }
    
    ## strip inserted pos/neg so that they match original scaling factor files
    #if ($prepended_pos_flag)
    #{
    #    $pos_tomerge =~ s/^pos_//;
    #}
    #if ($prepended_neg_flag)
    #{
    #    $neg_tomerge =~ s/^neg_//;
    #}
    
    if ($pos_tomerge eq '')
    {
        $pos_tomerge = 'NULL';
    }
    if ($neg_tomerge eq '')
    {
        $neg_tomerge = 'NULL';
    }
    
    printf OUTFILE "%s\t%s\t%s\t%s\n",
        $tomerge_case, $pos_tomerge, $neg_tomerge, $tomerge_short;
}
close OUTFILE;



# insert missing row ID column
if (!defined($row_id_col))
{
    print "row ID\t";
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
    # insert missing row ID column
    if (!defined($row_id_col))
    {
        print "$row_id\t";
    }

    for ($col = 0; $col < @global_header_array; $col++)
    {
        $header = $global_header_array[$col];
        $value  = $global_data_hash{$row_id}{$header};
        
        if (!defined($value))
        {
            $value = '';
        }
        
        #$value = reformat_sci($value);

        # shift pos/neg sample data
        $pos_neg = $row_id_pos_neg_hash{$row_id};
        if ($equal_means_flag &&
            defined($actual_sample_lc_hash{$header}) &&
            is_number($value))
        {
            # data is logged, shift the values
            if ($all_logged_flag)
            {
                if ($pos_neg eq 'pos')
                {
                    $value += $pos_shift;
                }
                elsif ($pos_neg eq 'neg')
                {
                    $value += $neg_shift;
                }
            }
            # data is unlogged, scale the values
            #
            # WARNING: negative values will likely behave strangely
            else
            {
                if ($pos_neg eq 'pos')
                {
                    $value *= $pos_scale;
                }
                elsif ($pos_neg eq 'neg')
                {
                    $value *= $neg_scale;
                }
            }
        }


        if ($col)
        {
            print "\t";
        }
        
        print "$value";
    }
    print "\n";
}
