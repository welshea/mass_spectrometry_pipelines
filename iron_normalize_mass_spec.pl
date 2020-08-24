#!/usr/bin/perl -w

use Scalar::Util qw(looks_like_number);

# 2020-08-24 remove --bg background subtraction option, as it is a bad idea
# 2020-08-03: output scaling factor lines in correct TMT label N/C order
# 2020-07-27: hack in support for El-MAVEN metabolomics samples
#
# 2020-02-28: added --no-log2 --log2 --unlog2 flags
#             removed old --unlog flag and reimplemented it for --unlog2
#             prioritize ModificationID over RowIdentifier column for id
#              this is to catch cases where cleanup is done in the wrong order
#
# 2020-03-12: skip "Intensity L" and "Intensity H" summary columns
#             warning, will not correctly handle case where there are real
#              samples named H and L together with non-H/L sample names
#
# 2020-03-16: fixed support for TMT-###[NC] columns with text in front of them
# 2020-03-20: don't use .mzXML for identifiying mzmine output anymore
#             re-enable pre-pending IRON to metabolomics output
#
# 2020-04-22: added support for Skyline " Total Area" sample headers

# all functions are effectively macros, since there are no local variables

# we assume that there can be no more than 2 pools per plex

sub is_number
{
    # use what Perl thinks is a number first
    if (looks_like_number($_[0]))
    {
        # Perl treats infinities as numbers, Excel does not
        if ($_[0] =~ /^[+-]*inf/)
        {
            return 0;
        }
        
        return 1;
    }

    # Perl cannot handle American comma separators within long numbers.
    # Excel does, so we have to check for it.
    # Excel doesn't handle European dot separators, at least not when it is
    #  set to the US locale (my test environment).  I am going to leave this
    #  unsupported for now.
    #
    # return ($_[0] =~ /^([+-]?)[0-9]+(,\d\d\d)*([Ee]([+-]?[0-9]+))?$/);
}


sub cmp_scale_lines
{
    my $str_a = '';
    my $str_b = '';
    my $test_a;
    my $test_b;

    # sort normally
    $str_a = $a;
    $str_b = $b;
    
    # strip GlobalScale: from beginning
    $str_a =~ s/^GlobalScale:\s+//ig;
    $str_b =~ s/^GlobalScale:\s+//ig;
    
    # make sure header line is always first
    $test_a = ($str_a =~ /SampleID/ && $str_a =~ /PresentSample/);
    $test_b = ($str_b =~ /SampleID/ && $str_b =~ /PresentSample/);
    if ($test_a && $test_b == 0) { return -1; }
    if ($test_b && $test_a == 0) { return  1; }

    # sort N/C aware if TMT sample label
    if ($str_a =~ /^([^\t]*TMT-\d+)(N|C){0,1}\t/i)
    {
        # assume C if not specified, lowercase it for sorting
        $str_a = $1 . 'c';
        
        # make sure N is uppercase and C is lowercase, so we sort N before C
        if (defined($2) && $2 =~ /N/i)
        {
            $str_a =~ s/c$/N/;
        }
    }
    if ($str_b =~ /^([^\t]*TMT-\d+)(N|C){0,1}\t/i)
    {
        # assume C if not specified, lowercase it for sorting
        $str_b = $1 . 'c';
        
        # make sure N is uppercase and C is lowercase, so we sort N before C
        if (defined($2) && $2 =~ /N/i)
        {
            $str_b =~ s/c$/N/;
        }
    }
    
    return $str_a cmp $str_b;
}


sub read_in_data_file
{
    $infile = $_[0];
    open INFILE, "$infile" or die "ABORT -- can't open $infile\n";

    # read in header line
    $line = <INFILE>;
    $line =~ s/[\r\n]+//g;
    $line =~ s/\"//;

    @array = split /\t/, $line;
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\s+/ /g;
    }
    
    $id_col = '';
    if ($id_col eq '')
    {
        for ($i = 0; $i < @array; $i++)
        {
            # select the identifier column
            if ($array[$i] =~ /^ModificationID$/i)
            {
                $id_col = $i;
                last;
            }
        }
    }
    if ($id_col eq '')
    {
        for ($i = 0; $i < @array; $i++)
        {
            # select the identifier column
            if ($array[$i] =~ /^RowIdentifier$/i)
            {
                $id_col = $i;
                last;
            }
        }
    }
    if ($id_col eq '')
    {
        for ($i = 0; $i < @array; $i++)
        {
            # select the identifier column
            if ($array[$i] =~ /^row ID$/i)
            {
                $id_col = $i;
                last;
            }
        }
    }
    if ($id_col eq '')
    {
        for ($i = 0; $i < @array; $i++)
        {
            # select the identifier column
            if ($array[$i] =~ /^id$/i)
            {
                $id_col = $i;
                last;
            }
        }
    }
    if ($id_col eq '')
    {
        for ($i = 0; $i < @array; $i++)
        {
            # select the identifier column
            if ($array[$i] =~ /^Index$/i)
            {
                $id_col = $i;
                last;
            }
        }
    }
    if ($id_col eq '')
    {
        $id_col = 0;
    }
    
    printf STDERR "Row identifier column:\t%s\n", $array[$id_col];

    @sample_array = 0;
    $num_samples = 0;
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\s+/ /g;
        
        $field = $array[$i];

        $count = 0;
        if ($field =~ /^Intensity/i)
        {
            # skip summary columns
            if ($field =~ /^Intensity$/i)
            {
                next;
            }
            if ($field =~ /^Intensity H$/i)
            {
                next;
            }
            if ($field =~ /^Intensity L$/i)
            {
                next;
            }
            
            # skip individual fractions
            if ($field =~ /___\d+$/)
            {
                next;
            }
            
            $count++;
        }
        
        if ($field =~ /^Intensity/i)
        {
            # skip summary columns
            if ($field =~ /^Intensity$/i)
            {
                next;
            }
            
            # skip H and L columns only if they aren't the only ones
            if ($count)
            {
                if ($field =~ /^Intensity H$/i)
                {
                    next;
                }
                if ($field =~ /^Intensity L$/i)
                {
                    next;
                }
            }
            
            # skip individual fractions
            if ($field =~ /___\d+$/)
            {
                next;
            }

            $sample_to_file_col_hash{$field} = $i;
            $sample_array[$num_samples++] = $field;
        }
        
        # check for mzMine abundances
        if ($field =~ / Peak height$/i ||
            $field =~ / Peak area$/i)
        {
            $sample_to_file_col_hash{$field} = $i;
            $sample_array[$num_samples++] = $field;

            $global_metabolomics_flag = 1;
        }
        
        # check for Skyline peak areas
        if ($field =~ / Total Area$/i)
        {
            $sample_to_file_col_hash{$field} = $i;
            $sample_array[$num_samples++] = $field;
        }
        

        # check to see if it is cleaned up TMT column headers
        if ($field =~ /^TMT/i)
        {
            # skip individual fractions
            if ($field =~ /___\d+$/)
            {
                next;
            }

            $sample_to_file_col_hash{$field} = $i;
            $sample_array[$num_samples++] = $field;
        }
        elsif ($field =~ /TMT\-\d+[NC]*$/i)
        {
            # skip individual fractions
            if ($field =~ /___\d+$/)
            {
                next;
            }

            $sample_to_file_col_hash{$field} = $i;
            $sample_array[$num_samples++] = $field;
        }
    }
    
    # didn't find any, maybe it is metabolomics with _pos or _neg
    if ($num_samples == 0)
    {
      for ($i = 0; $i < @array; $i++)
      {
        $field = $array[$i];

        $count = 0;

        if ($field =~ /(^|[^A-Za-z0-9]+)(pos|neg)([^A-Za-z0-9]+|$)/i)
        {
            $sample_to_file_col_hash{$field} = $i;
            $sample_array[$num_samples++] = $field;

            $global_metabolomics_flag = 1;
        }
      }
    }

    $num_header_cols = @array;
    $header_line = join "\t", @array;

    # read in the rest of the data
    $row = 0;
    while(defined($line=<INFILE>))
    {
        $line =~ s/[\r\n]+//g;
        $line =~ s/\"//;

        @array = split /\t/, $line;

        for ($i = 0; $i < @array; $i++)
        {
            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;
            $array[$i] =~ s/\s+/ /g;
        }
        
        for ($i = 0; $i < @array; $i++)
        {
            $orig_data_array[$row][$i] = $array[$i];
        }
        
        $row++;
    }
    $num_rows = $row;

    # unlog2 the sample data
    if ($unlog2_flag)
    {
        for ($row = 0; $row < $num_rows; $row++)
        {
            for ($i = 0; $i < $num_samples; $i++)
            {
                $sample = $sample_array[$i];
                $col    = $sample_to_file_col_hash{$sample};
                
                $value  = $orig_data_array[$row][$col];
                
                if (defined($value) && is_number($value))
                {
                    $orig_data_array[$row][$col] = 2.0 ** $value;
                }
            }
        }
    }
}


# store separate condensed data arrays
sub store_condensed_data
{
    for ($row = 0; $row < $num_rows; $row++)
    {
        for ($i = 0; $i < $num_samples; $i++)
        {
            $sample = $sample_array[$i];
            $col    = $sample_to_file_col_hash{$sample};

            $sample_to_condensed_col_hash{$sample} = $i;

            $value = $orig_data_array[$row][$col];

            if (defined($value) && is_number($value) && $value > 0)
            {
                $condensed_data_array[$row][$i] = $value;
            }            
        }
    }
}


$log2 = log(2.0);
sub iron_samples
{
    # open IRON input file for writing
    open INPUT_FOR_IRON, ">$iron_input_name" or die "ABORT -- can't open $iron_input_name\n";

    # print header
    printf INPUT_FOR_IRON "%s", 'ProbeID';
    for ($i = 0; $i < $num_samples; $i++)
    {
        # strip Intensity from the beginning
        $sample = $sample_array[$i];
        $sample =~ s/^Intensity\s+//i;
        
        # strip mzMine stuff
        $sample =~ s/\.mzX?ML[^.]+$//i;
        $sample =~ s/ Peak \S+$//i;

        # strip Skyline stuff
        $sample =~ s/ Total Area$//i;
        
        printf INPUT_FOR_IRON "\t%s", $sample;
    }
    printf INPUT_FOR_IRON "\n";

    for ($row = 0; $row < $num_rows; $row++)
    {
        $identifier = $orig_data_array[$row][$id_col];
        printf INPUT_FOR_IRON "%s", $identifier;

        for ($i = 0; $i < $num_samples; $i++)
        {
            $sample = $sample_array[$i];
            $col = $sample_to_condensed_col_hash{$sample};

            $value = $condensed_data_array[$row][$col];
            if (!defined($value))
            {
                $value = '';
            }

            printf INPUT_FOR_IRON "\t%s", $value;
        }
        printf INPUT_FOR_IRON "\n";
    }
    close INPUT_FOR_IRON;

    # find the median sample
    if ($force_ref_sample =~ /\S/)
    {
        $median_sample = $force_ref_sample;
    }
    else
    {
        $exclusions_string = '';
        $spikeins_string   = '';

        if ($exclusions_flag)
        {
            $exclusions_string = sprintf "-x \"%s\"", $exclusions_filename;
        }
        if ($spikeins_flag)
        {
            $spikeins_string = sprintf "-S \"%s\"", $spikeins_filename;
        }

        $cmd_string = sprintf "findmedian --pearson --spreadsheet %s %s \"%s\" 2>/dev/null > \"%s\"",
            $exclusions_string, $spikeins_string,
            $iron_input_name, $findmedian_output_name;

        `$cmd_string`;
        $cmd_string = sprintf "tail -2 \"%s\" | head -1 | cut -f 4",
            $findmedian_output_name;
        $median_sample = `$cmd_string`;
        $median_sample =~ s/[\r\n]+//;
    }
    
    # re-run to save findmedian results
    $exclusions_string = '';
    $spikeins_string   = '';
    if ($exclusions_flag)
    {
        $exclusions_string = sprintf "-x \"%s\"", $exclusions_filename;
    }
    if ($spikeins_flag)
    {
        $spikeins_string = sprintf "-S \"%s\"", $spikeins_filename;
    }

    $cmd_string = sprintf "findmedian --pearson --spreadsheet %s %s %s 2>/dev/null > %s",
        $exclusions_string, $spikeins_string,
        $iron_input_name, $findmedian_output_name;

    # run IRON
    printf STDERR "Running IRON:\t%s\n", $filename;
    printf STDERR "Median sample %s:\t%s\n", $filename, $median_sample;

    $exclusions_string = '';
    $spikeins_string   = '';
    $bg_string         = '';
    if ($exclusions_flag)
    {
        $exclusions_string = sprintf "-x \"%s\"", $exclusions_filename;
    }
    if ($spikeins_flag)
    {
        $spikeins_string = sprintf "-S \"%s\"", $spikeins_filename;
    }
    if ($bg_flag)
    {
        $bg_string = '--bg-global';
    }

    $cmd_string = sprintf "iron_generic --proteomics --norm-iron=\"%s\" %s %s %s \"%s\" -o \"%s\" |& grep -P \"^GlobalScale\"",
        $median_sample, $exclusions_string, $spikeins_string, $bg_string,
        $iron_input_name, $iron_output_name;

    @scale_lines = `$cmd_string`;
    @scale_lines = sort cmp_scale_lines @scale_lines;

    # output scaling factor file
    open OUTPUT_SCALING, ">$scales_output_name" or die "can't open scales output file $scales_output_name";
    for ($i = 0; $i < @scale_lines; $i++)
    {
        $line = $scale_lines[$i];
        $line =~ s/[\r\n]+//g;
        
        # strip GlobalScale: from beginning
        $line =~ s/^GlobalScale:\s+//ig;
        
        printf OUTPUT_SCALING "%s\n", $line;
    }
    close OUTPUT_SCALING;

    # open IRON output for reading
    open OUTPUT_FOR_IRON, "$iron_output_name" or die "ABORT -- can't open $iron_output_name\n";
    
    # skip header line
    $line = <OUTPUT_FOR_IRON>;

    # read in the rest of the data into a new log2 array
    $row = 0;
    while(defined($line=<OUTPUT_FOR_IRON>))
    {
        $line =~ s/[\r\n]+//g;
        $line =~ s/\"//;

        @array = split /\t/, $line;

        for ($i = 0; $i < @array; $i++)
        {
            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;
            $array[$i] =~ s/\s+/ /g;
        }
        
        for ($i = 0; $i < $num_samples; $i++)
        {
            $value = $array[$i+1];
            
            # paranoia for log taking
            if (is_number($value) && $value > 0)
            {
                $value = log($value) / $log2;
            }
            else
            {
                $value = '';
            }
            
            $iron_normalized_array[$row][$i] = $value;
        }
        
        $row++;
    }
    close OUTPUT_FOR_IRON;
}


sub output_final_data
{
    # print original header
    # skip original Intensity columns
    $header_line_new = '';
    @header_array_orig = split /\t/, $header_line;
    for ($i = 0; $i < @header_array_orig; $i++)
    {
        $field = $header_array_orig[$i];
        
        # skip original Intensity columns
        if (defined($sample_to_file_col_hash{$field}))
        {
            next;
        }
        
        if ($header_line_new ne '')
        {
            $header_line_new .= "\t";
        }
        
        $header_line_new .= $field;
    }

    printf "%s", $header_line_new;

    # print header for new samples
    for ($i = 0; $i < $num_samples; $i++)
    {
        $sample = $sample_array[$i];
        
        if ($strip_sample_names_flag)
        {
            # strip Intensity from beginning
            if ($sample =~ /^Intensity\s+/)
            {
                $sample =~ s/^Intensity\s+//i;
            }

            # strip mzMine stuff, add IRON
            elsif ($sample =~ / Peak height$/i ||
                   $sample =~ / Peak area$/i)
            {
                $sample =~ s/\.mzX?ML[^.]+$//i;
                $sample =~ s/ Peak \S+$//i;
                
            }

            # strip Skyline stuff
            elsif ($sample =~ / Total Area$/i)
            {
                $sample =~ s/ Total Area$//i;
            }
            
            # TMT data
#            elsif ($sample =~ /^TMT/)
#            {
#            }

            # we're going to merge POS and NEG together later
            # prepend IRON so that it is easy to detect sample columns
            if ($global_metabolomics_flag)
            {
                $sample = 'IRON ' . $sample;
            }
        }

        printf "\t%s", $sample;
    }
    printf "\n";

    # print the rest of the table
    for ($row = 0; $row < $num_rows; $row++)
    {
        @array = ();

        # count number present
        $count = 0;
        for ($i = 0; $i < $num_samples; $i++)
        {
            $value = $iron_normalized_array[$row][$i];
            
            if (defined($value) && is_number($value))
            {
                $count++;
            }
        }
        
        # skip rows with no data
        if ($count == 0)
        {
            next;
        }

        # print original columns
        $col_new = 0;
        for ($col = 0; $col < $num_header_cols; $col++)
        {
            $field = $header_array_orig[$col];

            if (defined($sample_to_file_col_hash{$field}))
            {
                next;
            }

            $value = $orig_data_array[$row][$col];
            
            $array[$col_new++] = $value;
        }

        $line_new = join "\t", @array;
        printf "%s", $line_new;
        
        # print new IRON'd columns
        for ($i = 0; $i < $num_samples; $i++)
        {
            $value = $iron_normalized_array[$row][$i];

            if ($no_log2_flag && defined($value) && is_number($value))
            {
                $value = 2.0 ** $value;
            }
            
            printf "\t%s", $value;
        }
        printf "\n";
    }
}


# read in command line arguments
$force_ref_sample        = '';
$num_files               = 0;
$syntax_error_flag       = 0;
$exclusions_flag         = 0;
$spikeins_flag           = 0;
$bg_flag                 = 0;
$strip_sample_names_flag = 1;
$unlog2_flag             = 0;    # unlog2 the input data
$no_log2_flag            = 0;   # do not log2 the output data
$global_metabolomics_flag = 0;

for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        if ($field =~ /^--iron-exclusions\=(\S+)/)
        {
            $exclusions_flag = 1;
            $exclusions_filename = $1;
            
            printf STDERR "Using exclusions file: %s\n", $exclusions_filename;
        }
        elsif ($field =~ /^--iron-spikeins\=(\S+)/)
        {
            $spikeins_flag = 1;
            $spikeins_filename = $1;
            
            printf STDERR "Using spikeins file: %s\n", $spikeins_filename;
        }
        # this appears to be a horrible idea,
        # even for metabolomics, which may have a lot of background
        # elsif ($field =~ /^--bg/)
        #{
        #    $bg_flag = 1;
        #    printf STDERR "Using RMA background subtraction\n";
        #}
        elsif ($field =~ /^--no-strip-sample-names/)
        {
            $strip_sample_names_flag = 0;
            printf STDERR "Leaving sample names untouched\n";
        }
        elsif ($field =~ /^--no-log2$/)
        {
            $no_log2_flag = 1;
        }
        elsif ($field =~ /^--log2$/)
        {
            $no_log2_flag = 0;
        }
        elsif ($field =~ /^--unlog2$/)
        {
            $unlog2_flag = 1;
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
            $filename = $field;
            $num_files++;
        }
        # treat further arguments as the reference sample
        else
        {
            $force_ref_sample = $field
        }
    }
}


$process_id = $$;

$scales_output_name =  $filename;
$scales_output_name =~ s/\.txt$//i;
$findmedian_output_name = sprintf "%s_findmedian.txt", $scales_output_name;
$scales_output_name = sprintf "%s_scaling_factors.txt", $scales_output_name;

$iron_input_name  = sprintf "temp_iron_input_%s.txt", $process_id;
$iron_output_name = sprintf "temp_iron_output_%s.txt", $process_id;

read_in_data_file($filename);
store_condensed_data();
iron_samples();
output_final_data();

`rm \"$iron_input_name\"`;
`rm \"$iron_output_name\"`;
