#!/usr/bin/perl -w

# treat any options after the filename as the pool channels
# treat them as free text, don't try to get smart,
# the user will need to specify the channel names as they rae
# in the input data file
#
# Don't forget that current file format is ex: TMT-126, not just 126
#
# 2020-09-22: respect --no-debatch flag when combined with --comp-pool
# 2020-09-22: --comp-pool uses comp pool for scaling, instead of norm pools
# 2020-08-03: output scaling factors in proper TMT N/C label order
# 2020-07-15: change TMT-126 to TMT-126C
# 2020-03-06: replace | with ~ in RowIdentifier column
# 2020-02-28: added --no-log2 --log2 --unlog2 flags
#             prioritize ModificationID over RowIdentifier column for id
#              this is to catch cases where cleanup is done in the wrong order
# 2020-02-24: added --leave-ratios flag to NOT scale ratios into abundances
#


use Scalar::Util qw(looks_like_number);


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
    
    $id_col = -1;
    
    if ($id_col == -1)
    {
        for ($i = 0; $i < @array; $i++)
        {
            if ($array[$i] eq 'ModificationID')
            {
                $id_col = $i;
                last;
            }
        }
    }
    if ($id_col == -1)
    {
        for ($i = 0; $i < @array; $i++)
        {
            if ($array[$i] eq 'RowIdentifier')
            {
                $id_col = $i;
                last;
            }
        }
    }
    
    if ($id_col == -1)
    {
        die "ABORT -- row identifier column not found\n";
    }

    $num_samples = 0;
    for ($i = 0; $i < @array; $i++)
    {
        if ($array[$i] =~ /TMT\-\d\d\d/)
        {
            $array[$i] =~ s/[^A-Za-z0-9-]+/_/g;
            $array[$i] =~ s/_+/_/g;
            $array[$i] =~ s/^_+//;
            $array[$i] =~ s/_+$//;
            
            $sample_name = $array[$i];
            
            if ($array[$i] =~ /(.*?)_(TMT-\d\d\d[NC]*)/)
            {
                $tmt_plex = $1;
                $channel  = $2;

                $sample_to_file_col_hash{$sample_name} = $i;
                
                $tmt_plex_hash{$tmt_plex}{$channel} = $sample_name;

                $channel_hash{$channel} = 1;
            }
        }
    }
    $num_header_cols = @array;
    $header_line = join "\t", @array;

    @sample_array = sort keys %sample_to_file_col_hash;
    @channel_array = sort keys %channel_hash;
    @tmt_plex_array = sort keys %tmt_plex_hash;

    $num_samples  = @sample_array;
    $num_channels = @channel_array;
    $num_plexes   = @tmt_plex_array;
    
    # generate channel string to @channel_array index lookup table
    for ($i = 0; $i < $num_channels; $i++)
    {
        $channel_string = $channel_array[$i];
        $channel_string_to_index_hash{$channel_string} = $i;
    }

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
            if ($i == $id_col)
            {
                # HACK -- Evince search is broken if identifier contains |
                # replace | with ~
                $array[$i] =~ s/\|/\~/g;
            }
        
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
                $transposed_log_array[$i][$row] = log($value);
            }            
        }
    }
}


# crudely scale each log sample to be roughly normalized
# assume they all have roughly the same distributions
# scale each 6-plex independently
sub normalize_crude
{
    @temp_averages = ();

    for ($p = 0; $p < $num_plexes; $p++)
    {
        $tmt_plex = $tmt_plex_array[$p];

        for ($ch = 0; $ch < $num_channels; $ch++)
        {
            $sample = $tmt_plex_hash{$tmt_plex}{$channel_array[$ch]};
            $col = $sample_to_condensed_col_hash{$sample};

            $avg = 0;
            $n   = 0;
            
            for ($row = 0; $row < $num_rows; $row++)
            {
                $value = $transposed_log_array[$col][$row];
                if (defined($value))
                {
                    $avg += $value;
                    $n++;
                }
            }
            
            if ($n)
            {
                $avg /= $n;
            }
            
            $temp_averages[$ch] = $avg;
        }

        $avg = 0;
        for ($ch = 0; $ch < $num_channels; $ch++)
        {
            $avg += $temp_averages[$ch];
        }
        $avg /= $num_channels;
        
        for ($ch = 0; $ch < $num_channels; $ch++)
        {
            $sample = $tmt_plex_hash{$tmt_plex}{$channel_array[$ch]};
            $col = $sample_to_condensed_col_hash{$sample};

            $offset = $avg - $temp_averages[$ch];

            $scale = 1.0;
            if ($temp_averages[$ch])
            {
                $scale = $avg / $temp_averages[$ch];
            }
            
#   #        printf STDERR "SCALE\t%s\t%s\n",
#   #            $channel_array[$ch], $scale;

            for ($row = 0; $row < $num_rows; $row++)
            {
                $value = $transposed_log_array[$col][$row];

                if (defined($value))
                {
                    # shift the log values
                    $transposed_log_array[$col][$row] += $offset;
                    
                    # scale the unlogged values
                    $condensed_data_array[$row][$col] *= $scale;
                }
            }
        }
    }
}


# scan through the 6-plexes to automatically identify the pool channels
# assume the two most similar channels are the pools
sub identify_pools
{
    @fixed_pool_array = sort keys %fixed_pool_hash;
    $num_fixed_pools = @fixed_pool_array;
    
    # support "variable", since default is now TMT-126C
    if (defined($fixed_pool_hash{'variable'}))
    {
        @fixed_pool_array = ();
        $num_fixed_pools = 0;
    }

    # fixed pool channel(s)
    if ($num_fixed_pools > 0)
    {
        for ($fp = 0; $fp < $num_fixed_pools; $fp++)
        {
            $channel_string = $fixed_pool_array[$fp];
            $ch = $channel_string_to_index_hash{$channel_string};

            for ($p = 0; $p < $num_plexes; $p++)
            {
                $plex_pool_channels[$p][$fp] = $ch;
            }
        }
        
        $num_pools = $num_fixed_pools;
    }


    # calculate RMSD of the fixed-position pools
    if ($num_fixed_pools > 1)
    {
      for ($p = 0; $p < $num_plexes; $p++)
      {
        $tmt_plex = $tmt_plex_array[$p];

        $sample1 = $tmt_plex_hash{$tmt_plex}{$fixed_pool_array[0]};
        $sample2 = $tmt_plex_hash{$tmt_plex}{$fixed_pool_array[1]};

        $s1 = $sample_to_condensed_col_hash{$sample1};
        $s2 = $sample_to_condensed_col_hash{$sample2};

        $rmsd = 0;
        $n    = 0;

        for ($i = 0; $i < $num_rows; $i++)
        {
            $value1 = $transposed_log_array[$s1][$i];
            $value2 = $transposed_log_array[$s2][$i];

            # both points exist
            if (defined($value1) && defined($value2))
            {
                $diff = $value2 - $value1;
                $rmsd += $diff * $diff;
                $n++;
            }
        }
        
        if ($n)
        {
            $rmsd = sqrt($rmsd / $n);
        }
        
        $plex_pool_channels[$p][2] = $rmsd;

        if (1)
        {
            $ch1 = $plex_pool_channels[$p][0];
            $ch2 = $plex_pool_channels[$p][1];
            
            $pool1 = $channel_array[$ch1];
            $pool2 = $channel_array[$ch2];

            $pool1 =~ s/^TMT-//;
            $pool2 =~ s/^TMT-//;

            printf STDERR "POOL\t%s\t%s\t%s\t%f\n",
                $tmt_plex,
                $pool1,
                $pool2,
                $rmsd;
        }
      }
    }
    if ($num_fixed_pools == 1)
    {
        printf STDERR "POOL\t%s\t%s\n",
            'fixed_pool_ch', $fixed_pool_array[0];
    }
    
    
    # otherwise, assume 2 pools per plex, automatically guess them
    if ($num_fixed_pools == 0)
    {
      $num_pools = 2;
    
      for ($p = 0; $p < $num_plexes; $p++)
      {
        $tmt_plex = $tmt_plex_array[$p];
        
        # [0] [1] are the pools channels
        # [2] is the best RMSD
        $plex_pool_channels[$p][0] = '';
        $plex_pool_channels[$p][1] = '';
        $plex_pool_channels[$p][2] = '';
        
        $best_rmsd = 9E99;
        for ($ch1 = 0; $ch1 < $num_channels; $ch1++)
        {
            $sample1 = $tmt_plex_hash{$tmt_plex}{$channel_array[$ch1]};
            $s1 = $sample_to_condensed_col_hash{$sample1};

            for ($ch2 = $ch1+1; $ch2 < $num_channels; $ch2++)
            {
                $sample2 = $tmt_plex_hash{$tmt_plex}{$channel_array[$ch2]};
                $s2 = $sample_to_condensed_col_hash{$sample2};

                $rmsd = 0;
                $n    = 0;

                for ($i = 0; $i < $num_rows; $i++)
                {
                    $value1 = $transposed_log_array[$s1][$i];
                    $value2 = $transposed_log_array[$s2][$i];

                    # both points exist
                    if (defined($value1) && defined($value2))
                    {
                        $diff = $value2 - $value1;
                        $rmsd += $diff * $diff;
                        $n++;
                    }
                }
                
                if ($n)
                {
                    $rmsd = sqrt($rmsd / $n);
                }
                
                if ($rmsd < $best_rmsd)
                {
                    $best_rmsd = $rmsd;

                    $plex_pool_channels[$p][0] = $ch1;
                    $plex_pool_channels[$p][1] = $ch2;
                    $plex_pool_channels[$p][2] = $best_rmsd;
                }
            }
        }

        if (1)
        {
            $ch1 = $plex_pool_channels[$p][0];
            $ch2 = $plex_pool_channels[$p][1];
            
            $pool1 = $channel_array[$ch1];
            $pool2 = $channel_array[$ch2];

            $pool1 =~ s/^TMT-//;
            $pool2 =~ s/^TMT-//;

            printf STDERR "POOL\t%s\t%s\t%s\t%f\n",
                $tmt_plex,
                $pool1,
                $pool2,
                $best_rmsd;
        }
      }
    }
}


# print warnings to STDERR if the pools look too different
sub check_outlier_pools
{
    if ($num_fixed_pools == 1)
    {
        return;
    }

    # scan for plexes with outlier pools
    $avg = 0;
    $sd  = 0;

    for ($p = 0; $p < $num_plexes; $p++)
    {
        $avg += $plex_pool_channels[$p][2];
    }
    $avg /= $num_plexes;

    for ($p = 0; $p < $num_plexes; $p++)
    {
        $diff = $plex_pool_channels[$p][2] - $avg;
        $sd += $diff * $diff;
    }
    $sd = sqrt($sd / $num_plexes);

    for ($p = 0; $p < $num_plexes; $p++)
    {
        $diff = $plex_pool_channels[$p][2] - $avg;

        if ($diff > 3 * $sd)
        {
            $tmt_plex = $tmt_plex_array[$p];

            $ch1 = $plex_pool_channels[$p][0];
            $ch2 = $plex_pool_channels[$p][1];
        
            $pool1 = $channel_array[$ch1];
            $pool2 = $channel_array[$ch2];

            $pool1 =~ s/^TMT-//;
            $pool2 =~ s/^TMT-//;

            printf STDERR "WARNPOOL\t%s\t%s\t%s\t%.4f\t%.4f\n",
                $tmt_plex,
                $pool1,
                $pool2,
                $plex_pool_channels[$p][2],
                $avg;
        }
    }
}


# normalize pan-pool cohort for estimating pseudo-abundance scales
sub iron_pools
{
    # open IRON input file for writing
    open INPUT_FOR_IRON, ">$iron_input_name" or die "ABORT -- can't open $iron_input_name\n";

    # print header stuff
    printf INPUT_FOR_IRON "%s", 'ProbeID';
    for ($p = 0; $p < $num_plexes; $p++)
    {
        $tmt_plex = $tmt_plex_array[$p];
        
        # use computational pool instead of real pool
        if ($comp_pool_flag)
        {
            printf INPUT_FOR_IRON "\tcomp_pool_%s", $tmt_plex;
        }
        else
        {
            for ($i = 0; $i < $num_pools; $i++)
            {
                $pool_ch   = $plex_pool_channels[$p][$i];
                $pool_samp = $tmt_plex_hash{$tmt_plex}{$channel_array[$pool_ch]};
        
                printf INPUT_FOR_IRON "\t%s", $pool_samp;
            }
        }
    }
    printf INPUT_FOR_IRON "\n";

    # output the rows
    for ($row = 0; $row < $num_rows; $row++)
    {
        $identifier = $orig_data_array[$row][$id_col];
        printf INPUT_FOR_IRON "%s", $identifier;

        for ($p = 0; $p < $num_plexes; $p++)
        {
            $tmt_plex = $tmt_plex_array[$p];

            # use computational pool instead of real pool
            if ($comp_pool_flag)
            {
                $log_avg = 0;
                $n = 0;

                for ($ch = 0; $ch < $num_channels; $ch++)
                {
                    $sample = $tmt_plex_hash{$tmt_plex}{$channel_array[$ch]};
                    $col = $sample_to_condensed_col_hash{$sample};
                
                    $value = $condensed_data_array[$row][$col];
                    if (defined($value))
                    {
                        # data is unlogged, so best to take geometric mean
                        $log_avg += log($value);
                        $n++;
                    }
                }
                
                $value = '';
                
                if ($n)
                {
                    $value = exp($log_avg / $n);
                }

                printf INPUT_FOR_IRON "\t%s", $value;
            }
            else
            {
                for ($i = 0; $i < $num_pools; $i++)
                {
                    $pool_ch = $plex_pool_channels[$p][$i];
                    $pool_samp =
                        $tmt_plex_hash{$tmt_plex}{$channel_array[$pool_ch]};
                    $pool_s = $sample_to_condensed_col_hash{$pool_samp};

                    $value = $condensed_data_array[$row][$pool_s];

                    if (!defined($value))
                    {
                        $value = '';
                    }

                    printf INPUT_FOR_IRON "\t%s", $value;
                }
            }
        }
        printf INPUT_FOR_IRON "\n";
    }
    close INPUT_FOR_IRON;


    # IRON the pools
    if ($no_iron_flag == 0)
    {
        # find the median pool
        $cmd_string = sprintf "findmedian --pearson --spreadsheet \"%s\" 2>/dev/null | tail -2 | head -1 | cut -f 4",
            $iron_input_name;
        $median_sample = `$cmd_string`;
        $median_sample =~ s/[\r\n]+//;


        # IRON the pools
        printf STDERR "Running IRON:\tpools\n";

        if ($exclusions_flag)
        {
          $cmd_string = sprintf "iron_generic --proteomics --iron-exclusions=\"%s\" --norm-iron=\"%s\" \"%s\" -o \"%s\" |& grep -P \"^GlobalScale\"",
              $exclusions_filename, $median_sample, $iron_input_name, $iron_output_name;
        }
        else
        {
          $cmd_string = sprintf "iron_generic --proteomics --norm-iron=\"%s\" \"%s\" -o \"%s\" |& grep -P \"^GlobalScale\"",
              $median_sample, $iron_input_name, $iron_output_name;
        }

        `$cmd_string`;
        # `cp -a $iron_output_name debug_iron_pools.txt`;

        @scale_lines = `$cmd_string`;
        @scale_lines = sort cmp_scale_lines @scale_lines;

        # output scaling factor file
        $log2 = log(2.0);
        open OUTPUT_SCALING, ">>$scales_output_name" or die "ABORT -- can't open scales output file $scales_output_name\n";
        for ($i = 0; $i < @scale_lines; $i++)
        {
            $line = $scale_lines[$i];
            $line =~ s/[\r\n]+//g;

            # strip GlobalScale: from beginning
            $line =~ s/^GlobalScale:\s+//ig;
            
            # only print the header line once
            if ($i == 0)
            {
                printf OUTPUT_SCALING "%s",   'SampleName';
                printf OUTPUT_SCALING "\t%s", 'Plex';
                printf OUTPUT_SCALING "\t%s", 'Scale';
                printf OUTPUT_SCALING "\t%s", 'Log2Scale';
                printf OUTPUT_SCALING "\t%s", 'AbsLog2Scale';
                printf OUTPUT_SCALING "\t%s", 'TrainingSet';
                printf OUTPUT_SCALING "\t%s", 'PresentBoth';
                printf OUTPUT_SCALING "\t%s", 'PresentSample';
                printf OUTPUT_SCALING "\t%s", 'PresentDataset';
                printf OUTPUT_SCALING "\t%s", 'FractionTrain';
                printf OUTPUT_SCALING "\n";
            
                next;
            }
            
            @array = split /\t/, $line;
            
            $sample          = $array[0];
            $scale           = $array[1];
            $training        = $array[3];
            $present_both    = $array[4];
            $present_sample  = $array[5];
            $present_dataset = $array[6];
            $frac_train      = $array[7];

            $plex            = 'PoolNorm';
            $log2_scale      = log($scale) / $log2;

            printf OUTPUT_SCALING "%s",   $sample;
            printf OUTPUT_SCALING "\t%s", $plex;
            printf OUTPUT_SCALING "\t%s", $scale;
            printf OUTPUT_SCALING "\t%s", $log2_scale;
            printf OUTPUT_SCALING "\t%s", abs($log2_scale);
            printf OUTPUT_SCALING "\t%s", $training;
            printf OUTPUT_SCALING "\t%s", $present_both;
            printf OUTPUT_SCALING "\t%s", $present_sample;
            printf OUTPUT_SCALING "\t%s", $present_dataset;
            printf OUTPUT_SCALING "\t%s", $frac_train;
            printf OUTPUT_SCALING "\n";
        }
        close OUTPUT_SCALING;
    }
    # copy input to output if we didn't run IRON
    else
    {
        `cp $iron_input_name $iron_output_name`;
    }


    # open IRON output for reading
    open OUTPUT_FOR_IRON, "$iron_output_name" or die "ABORT -- can't open $iron_output_name\n";
    $line = <OUTPUT_FOR_IRON>;
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
        
        $max_value = -9E99;
        $row_avg = 0;
        $n = 0;
        for ($i = 1; $i < @array; $i++)
        {
            $value = $array[$i];

            if (defined($value) && is_number($value) && $value > 0)
            {
                if ($value > $max_value)
                {
                    $max_value = $value;
                }
                
                $row_avg += log($value);
                $n++;
            }
        }
        
#        $row_pool_max_values[$row] = $max_value;
        
        if ($n)
        {
            $row_avg /= $n;
        }
        
        $row_pool_avg_values[$row] = exp($row_avg);
        
        $row++;
    }
    close OUTPUT_FOR_IRON;
}


sub iron_samples
{
    # output temp plex for iron
    for ($p = 0; $p < $num_plexes; $p++)
    {
        $tmt_plex = $tmt_plex_array[$p];
        
        $pool_ch1 = $plex_pool_channels[$p][0];
        $pool_samp1 = $tmt_plex_hash{$tmt_plex}{$channel_array[$pool_ch1]};
        $pool_s1 = $sample_to_condensed_col_hash{$pool_samp1};

        # HACK -- average the single pool with itself, so I don't have to
        # write more special case code for it
        if ($num_pools == 1)
        {
            $pool_ch2   = $pool_ch1;
            $pool_samp2 = $pool_samp1;
            $pool_s2    = $pool_s1;
        }
        else
        {
            $pool_ch2 = $plex_pool_channels[$p][1];
            $pool_samp2 = $tmt_plex_hash{$tmt_plex}{$channel_array[$pool_ch2]};
            $pool_s2 = $sample_to_condensed_col_hash{$pool_samp2};
        }
        
        # open IRON input file for writing
        open INPUT_FOR_IRON, ">$iron_input_name" or die "ABORT -- can't open $iron_input_name\n";

        # print header
        printf INPUT_FOR_IRON "%s", 'ProbeID';
        for ($ch = 0; $ch < $num_channels; $ch++)
        {
            $pool_samp2 = $tmt_plex_hash{$tmt_plex}{$channel_array[$pool_ch2]};
            printf INPUT_FOR_IRON "\t%s", $tmt_plex_hash{$tmt_plex}{$channel_array[$ch]}
        }
        printf INPUT_FOR_IRON "\t%s", '__AvgPool__';
        printf INPUT_FOR_IRON "\n";

        for ($row = 0; $row < $num_rows; $row++)
        {
            $identifier = $orig_data_array[$row][$id_col];
        
            $value1 = $condensed_data_array[$row][$pool_s1];
            $value2 = $condensed_data_array[$row][$pool_s2];

            $avg = '';
            if (defined($value1) && defined($value2))
            {
#                $avg = 0.5 * ($value1 + $value2);
                $avg = sqrt($value1 * $value2);
            }
            elsif (defined($value1))
            {
                $avg = $value1;
            }
            elsif (defined($value2))
            {
                $avg = $value2;
            }
            
            printf INPUT_FOR_IRON "%s", $identifier;
            for ($ch = 0; $ch < $num_channels; $ch++)
            {
                $sample = $tmt_plex_hash{$tmt_plex}{$channel_array[$ch]};
                $col = $sample_to_condensed_col_hash{$sample};

                $value = $condensed_data_array[$row][$col];
                if (!defined($value))
                {
                    $value = '';
                }

                printf INPUT_FOR_IRON "\t%s", $value;
            }
            printf INPUT_FOR_IRON "\t%s", $avg;
            printf INPUT_FOR_IRON "\n";
        }
        close INPUT_FOR_IRON;

        # IRON the samples
        if ($no_iron_flag == 0)
        {
            if (0)
            {
                # find the median sample
                $cmd_string = sprintf "findmedian --pearson --spreadsheet \"%s\" 2>/dev/null | tail -2 | head -1 | cut -f 4",
                    $iron_input_name;
                $median_sample = `$cmd_string`;
                $median_sample =~ s/[\r\n]+//;
            }
            else
            {
                $median_sample = '__AvgPool__';
            }

            # run IRON
            printf STDERR "Running IRON:\t%s\n", $tmt_plex;
#   #       printf STDERR "Median sample %s:\t%s\n", $filename, $median_sample;

            if ($exclusions_flag)
            {
              $cmd_string = sprintf "iron_generic --proteomics --iron-exclusions=\"%s\" --norm-iron=\"%s\" \"%s\" -o \"%s\" |& grep -P \"^GlobalScale\"",
                  $exclusions_filename, $median_sample, $iron_input_name, $iron_output_name;
            }
            else
            {
              $cmd_string = sprintf "iron_generic --proteomics --norm-iron=\"%s\" \"%s\" -o \"%s\" |& grep -P \"^GlobalScale\"",
                  $median_sample, $iron_input_name, $iron_output_name;
            }

            @scale_lines = `$cmd_string`;
            @scale_lines = sort cmp_scale_lines @scale_lines;

            # output scaling factor file
            $log2 = log(2.0);
            open OUTPUT_SCALING, ">>$scales_output_name" or die "ABORT -- can't open scales output file $scales_output_name\n";
            for ($i = 0; $i < @scale_lines; $i++)
            {
                $line = $scale_lines[$i];
                $line =~ s/[\r\n]+//g;

                # strip GlobalScale: from beginning
                $line =~ s/^GlobalScale:\s+//ig;
                
                if ($i == 0)
                {
                    if (0)
                    {
                        printf OUTPUT_SCALING "%s",   'SampleName';
                        printf OUTPUT_SCALING "\t%s", 'Plex';
                        printf OUTPUT_SCALING "\t%s", 'Scale';
                        printf OUTPUT_SCALING "\t%s", 'Log2Scale';
                        printf OUTPUT_SCALING "\t%s", 'AbsLog2Scale';
                        printf OUTPUT_SCALING "\t%s", 'TrainingSet';
                        printf OUTPUT_SCALING "\t%s", 'PresentBoth';
                        printf OUTPUT_SCALING "\t%s", 'PresentSample';
                        printf OUTPUT_SCALING "\t%s", 'PresentDataset';
                        printf OUTPUT_SCALING "\t%s", 'FractionTrain';
                        printf OUTPUT_SCALING "\n";
                    }
                
                    next;
                }
                
                # skip average pool
                if ($line =~ /^__AvgPool__/)
                {
                    next;
                }
                
                @array = split /\t/, $line;
                
                $sample          = $array[0];
                $scale           = $array[1];
                $training        = $array[3];
                $present_both    = $array[4];
                $present_sample  = $array[5];
                $present_dataset = $array[6];
                $frac_train      = $array[7];

                $plex            = $tmt_plex;
                $log2_scale      = log($scale) / $log2;
#   #           $sample          = sprintf "%s_%s", $plex, $sample;

                printf OUTPUT_SCALING "%s",   $sample;
                printf OUTPUT_SCALING "\t%s", $plex;
                printf OUTPUT_SCALING "\t%s", $scale;
                printf OUTPUT_SCALING "\t%s", $log2_scale;
                printf OUTPUT_SCALING "\t%s", abs($log2_scale);
                printf OUTPUT_SCALING "\t%s", $training;
                printf OUTPUT_SCALING "\t%s", $present_both;
                printf OUTPUT_SCALING "\t%s", $present_sample;
                printf OUTPUT_SCALING "\t%s", $present_dataset;
                printf OUTPUT_SCALING "\t%s", $frac_train;
                printf OUTPUT_SCALING "\n";
            }
            close OUTPUT_SCALING;
        }
        # copy input to output if we didn't run IRON
        else
        {
            `cp $iron_input_name $iron_output_name`;
        }


        # open IRON output for reading
        open OUTPUT_FOR_IRON, "$iron_output_name" or die "ABORT -- can't open $iron_output_name\n";
        
        # skip header line
        $line = <OUTPUT_FOR_IRON>;

        # read in the rest of the data into the condensed data array
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
            
            for ($i = 1, $ch = 0; $i <= $num_channels; $i++, $ch++)
            {
                $sample = $tmt_plex_hash{$tmt_plex}{$channel_array[$ch]};
                $col = $sample_to_condensed_col_hash{$sample};
                
                $value = $array[$i];

                if (defined($value) && is_number($value) && $value > 0)
                {
                    $condensed_data_array[$row][$col] = $value;
                }
            }
            
            $row++;
        }
        close OUTPUT_FOR_IRON;
    }
}


sub correct_abundances
{
    if ($comp_pool_flag)
    {
        print STDERR "De-batching with computational pools\n";
    }

    # overwrite condensed data with log2 ratio derived pseudo-abundances
    for ($p = 0; $p < $num_plexes; $p++)
    {
        $tmt_plex = $tmt_plex_array[$p];

        $pool_ch1 = $plex_pool_channels[$p][0];
        $pool_samp1 = $tmt_plex_hash{$tmt_plex}{$channel_array[$pool_ch1]};
        $pool_s1 = $sample_to_condensed_col_hash{$pool_samp1};


        # HACK -- average the single pool with itself, so I don't have to
        # write more special case code for it
        if ($num_pools == 1)
        {
            $pool_ch2   = $pool_ch1;
            $pool_samp2 = $pool_samp1;
            $pool_s2    = $pool_s1;
        }
        else
        {
            $pool_ch2 = $plex_pool_channels[$p][1];
            $pool_samp2 = $tmt_plex_hash{$tmt_plex}{$channel_array[$pool_ch2]};
            $pool_s2 = $sample_to_condensed_col_hash{$pool_samp2};
        }

        
        for ($row = 0; $row < $num_rows; $row++)
        {
            $avg = '';

            # use computational pool instead of real pool
            if ($comp_pool_flag)
            {
                $log_avg = 0;
                $n = 0;

                for ($ch = 0; $ch < $num_channels; $ch++)
                {
                    $sample = $tmt_plex_hash{$tmt_plex}{$channel_array[$ch]};
                    $col = $sample_to_condensed_col_hash{$sample};
                
                    $value = $condensed_data_array[$row][$col];
                    if (defined($value))
                    {
                        # data is unlogged, so best to take geometric mean
                        $log_avg += log($value);
                        $n++;
                    }
                }
                
                if ($n)
                {
                    $avg = exp($log_avg / $n);
                }
            }
            else
            {
                $value1 = $condensed_data_array[$row][$pool_s1];
                $value2 = $condensed_data_array[$row][$pool_s2];

                if (defined($value1) && defined($value2))
                {
                    # unlogged space, so take geometric mean
                    $avg = sqrt($value1 * $value2);
                }
                elsif (defined($value1))
                {
                    $avg = $value1;
                }
                elsif (defined($value2))
                {
                    $avg = $value2;
                }
            }
            
            for ($ch = 0; $ch < $num_channels; $ch++)
            {
                $sample = $tmt_plex_hash{$tmt_plex}{$channel_array[$ch]};
                $col = $sample_to_condensed_col_hash{$sample};
                
                # no average to divide by, undefine it
                if ($avg eq '')
                {
                    undef($condensed_data_array[$row][$col]);
                }

                $value = $condensed_data_array[$row][$col];
                if (defined($value))
                {
                    # HACK -- leave as ratios, rather than abundances
                    if ($leave_ratios_flag)
                    {
                        $condensed_data_array[$row][$col] = ($value / $avg);
                    }
                    # scale ratios back into abundances
                    else
                    {
                        $condensed_data_array[$row][$col] =
                            $row_pool_avg_values[$row] * ($value / $avg);
                    }
                }
            }
        }
    }
}


sub output_final_data
{
    %bad_row_hash = ();

    # overwrite original data with post-processed condensed data
    for ($row = 0; $row < $num_rows; $row++)
    {
        $blank_count = 0;

        for ($i = 0; $i < $num_samples; $i++)
        {
            $sample = $sample_array[$i];
            $col    = $sample_to_file_col_hash{$sample};
        
            $value = $condensed_data_array[$row][$i];

            if (defined($value) && is_number($value))
            {
                # print log2 value
                if ($no_log2_flag == 0)
                {
                    $orig_data_array[$row][$col] = log($value) / log(2.0);
                }
            }
            else
            {
                $orig_data_array[$row][$col] = '';
                $blank_count++;
            }
        }

        # flag rows that are all-blank data
        if ($blank_count == $num_samples)
        {
            $bad_row_hash{$row} = 1;
        }
    }
    


    # rearrange columns to put normalized abundances last

    @temp_sample_array = keys %sample_to_file_col_hash;
    %sample_col_hash = ();
    foreach $sample (@temp_sample_array)
    {
        $col = $sample_to_file_col_hash{$sample};
        $sample_col_hash{$col} = 1;
    }
    @sample_col_array = sort {$a<=>$b} keys %sample_col_hash;

    
    @new_col_order_array = ();
    # fill with non-samples first
    for ($i = 0, $col=0; $i < $num_header_cols; $i++)
    {
        if (!defined($sample_col_hash{$i}))
        {
            $new_col_order_array[$i] = $col++;
        }
    }
    # then put the samples at the end
    foreach $i (@sample_col_array)
    {
        $new_col_order_array[$i] = $col++;
    }
    
    # print header
    @header_array = split /\t/, $header_line;
    @header_array_new = ();
    for ($i = 0; $i < $num_header_cols; $i++)
    {
        $col = $new_col_order_array[$i];

        $header_array_new[$col] = $header_array[$i];
    }
    $header_line_new = join "\t", @header_array_new;
    printf "%s\n", $header_line_new;
    

    # print the rest of the table
    for ($row = 0; $row < $num_rows; $row++)
    {
        if (defined($bad_row_hash{$row}))
        {
            printf STDERR "EMPTYROW:\t%s\n",
                $orig_data_array[$row][$id_col];

            next;
        }

        @array = ();

        for ($col = 0; $col < $num_header_cols; $col++)
        {
            $value = $orig_data_array[$row][$col];
            
            $array[$col] = $value;
        }

        # reorder the columns
        @reordered_array = ();
        for ($i = 0; $i < $num_header_cols; $i++)
        {
            $col = $new_col_order_array[$i];

            $reordered_array[$col] = $array[$i];
        }

        $line_new = join "\t", @reordered_array;
        printf "%s\n", $line_new;
    }

    # output file containing the row-scales
    $log2 = log(2.0);
    open OUTFILE_ROW_SCALES, ">$row_scales_output_name" or die "ABORT -- can't open output row scales file $row_scales_output_name\n";
    printf OUTFILE_ROW_SCALES "%s\t%s\t%s\n", 'Identifier', 'RowScale', 'Log2Scale';
    for ($row = 0; $row < $num_rows; $row++)
    {
        if (defined($bad_row_hash{$row}))
        {
            next;
        }

        $sample = $orig_data_array[$row][$id_col];
        $scale  = $row_pool_avg_values[$row];
        
        printf OUTFILE_ROW_SCALES "%s\t%s\t%f\n",
            $sample, $scale, log($scale) / $log2;
    }
}


# read in command line arguments
$num_files         = 0;
$num_fixed_pools   = 0;
$exclusions_flag   = 0;
$no_debatch_flag   = 0;
$no_iron_flag      = 0;
$leave_ratios_flag = 0;
$unlog2_flag       = 0;    # unlog2 the input data
$no_log2_flag      = 0;    # do not log2 the output data
$comp_pool_flag    = 0;    # use computational pool instead of real pool

for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        if ($field =~ /^--iron-exclusions\=(\S+)/)
        {
            $exclusions_flag = 1;
            $exclusions_filename = $1;
            
            printf STDERR "Using exclusion file: %s\n", $exclusions_filename;
        }
        elsif ($field =~ /^--no-debatch$/)
        {
            $no_debatch_flag = 1;
        }
        elsif ($field =~ /^--debatch$/)
        {
            $no_debatch_flag = 0;
        }
        elsif ($field =~ /^--no-iron$/)
        {
            $no_iron_flag = 1;
        }
        elsif ($field =~ /^--iron$/)
        {
            $no_iron_flag = 0;
        }
        elsif ($field =~ /^--leave-ratios$/)
        {
            $leave_ratios_flag = 1;
        }
        elsif ($field =~ /^--comp-pool$/)
        {
            $comp_pool_flag = 1;
        }
        elsif ($field =~ /^--no-comp-pool$/)
        {
            $comp_pool_flag = 0;
        }
        else
        {
            printf "ABORT -- unknown option %s\n", $field;
            exit(1);
        }
    }
    else
    {
        if ($num_files == 0)
        {
            $filename = $field;
            $num_files++;
        }
        # treat any options after the filename as the pool channels
        # treat them as free text, don't try to get smart,
        # the user will need to specify the channel names as they rae
        # in the input data file
        #
        # Don't forget that current file format is ex: TMT-126C, not just 126
        else
        {
            $fixed_pool_hash{$field} = 1;
        }
    }
}


$process_id = $$;

$scales_output_name =  $filename;
$scales_output_name =~ s/\.txt$//i;
$scales_output_name = sprintf "%s_scaling_factors.txt", $scales_output_name;
`rm \"$scales_output_name\" >& /dev/null`;

$row_scales_output_name =  $filename;
$row_scales_output_name =~ s/\.txt$//i;
$row_scales_output_name = sprintf "%s_row_scales.txt", $row_scales_output_name;
`rm \"$row_scales_output_name\" >& /dev/null`;

$iron_input_name  = sprintf "temp_iron_input_%s.txt", $process_id;
$iron_output_name = sprintf "temp_iron_output_%s.txt", $process_id;

# default to using TMT-126C as the normalization channel
if (!defined(%fixed_pool_hash))
{
    $fixed_pool_hash{'TMT-126C'} = 1;
}

read_in_data_file($filename);
store_condensed_data();
#normalize_crude();
identify_pools();
check_outlier_pools();
iron_pools();
iron_samples();

if ($no_debatch_flag == 0)
{
    correct_abundances();
}

output_final_data();

`rm \"$iron_input_name\"`;
`rm \"$iron_output_name\"`;
