#!/usr/bin/perl -w

# 2021-08-06:  export auto-shortened sample names to sample name mapping table
# 2021-08-06:  export sample name mapping, warn on unmatched / duped pairs
# 2021-04-06:  check for and uniquify non-unique row IDs
# 2021-04-05:  fix El-MAVEN support, groupID is *NOT* a unique row identifier
# 2021-03-30:  don't convert ' ' to '_' in " Peak height" and " Peak area"
# 2021-03-30:  strip .raw from end of sample names
# 2021-01-06:  handle recent changes to strip_metabolomics_columns.pl headers
# 2020-09-25:  force row identifier header to always print as "row ID"
# 2020-09-25:  minor code clean up
# 2020-09-24:  clean up extra/dangling hypens in sample names
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
            if (!($sample =~ s/(^|[^A-Za-z0-9]+)pos([^A-Za-z0-9]+|$)/$2/i ||
                  $sample =~ s/pos$//i))
            {
                $all_pos_start_flag = 0;
            }
            $all_neg_start_flag = 0;
        }
        elsif ($pos_neg eq 'neg')
        {
            if (!($sample =~ s/(^|[^A-Za-z0-9]+)neg([^A-Za-z0-9]+|$)/$2/i ||
                  $sample =~ s/neg$//i))
            {
                $all_neg_start_flag = 0;
            }
            $all_pos_start_flag = 0;
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
        foreach $sample_col (@sample_col_array)
        {
            $sample          = $header_col_array[$sample_col];
            $sample          =~ s/^IRON //;
            $sample_origname = $sample;
            
            # strip pos/neg from sample name
            if ($all_pos_start_flag)
            {
                $sample =~ s/(^|[^A-Za-z0-9]+)pos([^A-Za-z0-9]+|$)/$2/i;
                $sample =~ s/pos$//i;

                # clean up underscores, etc.
                $sample =~ s/[_ ]+/_/g;
                $sample =~ s/\-+/\-/g;
                $sample =~ s/^[_ -]//;
                $sample =~ s/[_ -]$//;
            }
            elsif ($all_neg_start_flag)
            {
                $sample =~ s/(^|[^A-Za-z0-9]+)neg([^A-Za-z0-9]+|$)/$2/i;
                $sample =~ s/neg$//i;

                # clean up underscores, etc.
                $sample =~ s/[_ ]+/_/g;
                $sample =~ s/\-+/\-/g;
                $sample =~ s/^[_ -]//;
                $sample =~ s/[_ -]$//;
            }

            # strip .raw from end of sample name
            $sample =~ s/\.raw(?=(\s|[_ ]|$))//;
            
            # restore spaces in Peak height/area after conforming
            $sample =~ s/_Peak_(height|area)/ Peak $1/ig;


            $sample_lc = lc $sample;
            $sample_lc_to_origcase_hash{$sample_lc}{$sample} = 1;
            
            if ($all_pos_start_flag)
            {
                $sample_lc_to_origpos_hash{$sample_lc}{$sample_origname} = 1;
            }
            elsif ($all_neg_start_flag)
            {
                $sample_lc_to_origneg_hash{$sample_lc}{$sample_origname} = 1;
            }
            
            $global_data_hash{$row_id}{$sample_lc} = $array[$sample_col];
            $global_sample_hash{$sample_lc} = 1;
            
            foreach $metadata_col (@metadata_col_array)
            {
                $header = $header_col_array[$metadata_col];
                
                $global_data_hash{$row_id}{$header} = $array[$metadata_col];
                $global_metadata_hash{$header} = 1;
            }
        }
        
        $row++;
    }
}



# begin main()

$filename_pos   = shift;
$filename_neg   = shift;
$outname_user   = shift;    # optional, user-provided sample table filename

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


# warn on unmatched pairs
foreach $header (@global_sample_array)
{
    @temp_array = sort keys %{$sample_lc_to_origcase_hash{$header}};

    @pos_sample_matches = sort keys %{$sample_lc_to_origpos_hash{$header}};
    @neg_sample_matches = sort keys %{$sample_lc_to_origneg_hash{$header}};
    $count_pos_matches  = @pos_sample_matches;
    $count_neg_matches  = @neg_sample_matches;
    
    if ($count_pos_matches > 1)
    {
        $dupe_str = join "\t", @pos_sample_matches;

        printf STDERR "WARNING -- duplicate samples:\t%s\n", $dupe_str;
    }
    if ($count_neg_matches > 1)
    {
        $dupe_str = join "\t", @neg_sample_matches;

        printf STDERR "WARNING -- duplicate samples:\t%s\n", $dupe_str;
    }
    
    foreach $pos_sample (@pos_sample_matches)
    {
        if ($count_neg_matches == 0)
        {
            printf STDERR "WARNING -- unmatched pos/neg sample:\t%s\n",
                $pos_sample;
        }
    }
    foreach $neg_sample (@neg_sample_matches)
    {
        if ($count_pos_matches == 0)
        {
            printf STDERR "WARNING -- unmatched pos/neg sample:\t%s\n",
                $neg_sample;
        }
    }
}


foreach $header (@global_sample_array)
{
    @temp_array = sort keys %{$sample_lc_to_origcase_hash{$header}};
    
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



# guess at common prefix and/or suffix
%common_prefix_hash = ();
%common_suffix_hash = ();
for ($i = 0; $i < @global_sample_array; $i++)
{
    $sample_1     = $global_sample_array[$i];
    $sample_1_rev = reverse $sample_1;

    for ($j = $i+1; $j < @global_sample_array; $j++)
    {
        $sample_2     = $global_sample_array[$j];
        $sample_2_rev = reverse $sample_2;
        
        # find common prefix, Perl "dark magic"
        $xor = "$sample_1" ^ "$sample_2";
        $xor =~ /^\0*/;
        $len = $+[0];
        $common_prefix = substr $sample_1, 0, $len;

        # find common suffix, Perl "dark magic"
        $xor = "$sample_1_rev" ^ "$sample_2_rev";
        $xor =~ /^\0*/;
        $len = $+[0];
        $common_suffix = substr $sample_1_rev, 0, $len;
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
$count_sample_total  = @global_sample_array;
$common_prefix_all   = '';
$common_suffix_all   = '';

$count = 0;
foreach $common_prefix (@common_prefix_array)
{
    foreach $sample (@global_sample_array)
    {
        if ($sample =~ /^$common_prefix/)
        {
            $count++;
        }
    }
    
    if ($count == $count_sample_total)
    {
        $common_prefix_all = $common_prefix;
        last;
    }
}

$count = 0;
foreach $common_suffix (@common_suffix_array)
{
    foreach $sample (@global_sample_array)
    {
        if ($sample =~ /$common_suffix$/)
        {
            $count++;
        }
    }
    
    if ($count == $count_sample_total)
    {
        $common_suffix_all = $common_suffix;
        last;
    }
}



# output sample mapping table
# take first of duplicate samples (should never occur)
#
if (defined($outname_user))
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

foreach $header (@global_sample_array)
{
    @temp_array = sort keys %{$sample_lc_to_origcase_hash{$header}};

    @pos_sample_matches = sort keys %{$sample_lc_to_origpos_hash{$header}};
    @neg_sample_matches = sort keys %{$sample_lc_to_origneg_hash{$header}};
    $count_pos_matches  = @pos_sample_matches;
    $count_neg_matches  = @neg_sample_matches;

    $sample_case = $chosen_header_case_hash{$header};
    $pos_sample  = '';
    $neg_sample  = '';

    $sample_short = $sample_case;
    $sample_short =~ s/^$common_prefix_all//i;
    $sample_short =~ s/$common_suffix_all$//i;
    
    if ($count_pos_matches)
    {
        $pos_sample = $pos_sample_matches[0];
    }
    if ($count_neg_matches)
    {
        $neg_sample = $neg_sample_matches[0];
    }
    
    printf OUTFILE "%s\t%s\t%s\t%s\n",
        $sample_case, $pos_sample, $neg_sample, $sample_short;
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
