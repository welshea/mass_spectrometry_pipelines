#!/usr/bin/perl -w

# treat any options after the filename as the pool channels
# treat them as free text, don't try to get smart,
# the user will need to specify the channel names as they rae
# in the input data file
#
# Don't forget that current file format is ex: TMT-126, not just 126
#
# 2023-08-07: --iron-untilt overrides --no-iron flag if specified afterwards
# 2023-06-27: update is_number() to not treat NaNs as numbers
# 2023-05-18: strip UTF-8 BOM from MaxQuant 2.4 output, which broke many things
# 2022-12-09: change |& to 2>&1| to fix Ubuntu/Debian sh --> dash
# 2022-07-08: better handling of 126C default when no pool specified
# 2022-04-26: add --iron-untilt to usage statement (no longer experimental)
# 2022-04-19: strip accidental .txt from --comp-pool flags documentation
# 2022-01-24: add --comp-pool-exclusions-boost; exclude last and last -2
# 2022-01-19: ability to auto-flag and exclude dark samples from comp pool
# 2022-01-19: remove deprecated "variable" reference channel (now auto2)
# 2021-01-06: refactor scaling factor output, finish untilt support
# 2021-01-06: remove AbsLog2Scale from scaling factor output
# 2021-01-05: add experimental per-plex auto IRON reference sample picking
# 2021-12-22: add experimental --iron-untilt flag, not fully implemented
# 2021-12-17: add --comp-pool-exclusions flag to exclude samples from pool
# 2021-10-28: print Usage statement on command line error
# 2020-10-15: make sure split doesn't remove empty trailing fields
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

# any sample with a scaling factor more than +log2(6x) from that of the
# median log2 scaling factor is probably a bad sample
# this should be fairly conservative, we can maybe tighten it a bit
$log2_sf_above_median_cutoff = log(6.0) / log(2.0);


# all functions are effectively macros, since there are no local variables

# we assume that there can be no more than 2 pools per plex



sub is_number
{
    # use what Perl thinks is a number first
    # this is purely for speed, since the more complicated REGEX below should
    #  correctly handle all numeric cases
    if (looks_like_number($_[0]))
    {
        # Perl treats infinities as numbers, Excel does not.
        #
        # Perl treats NaN or NaNs, and various mixed caps, as numbers.
        # Weird that not-a-number is a number... but it is so that
        # it can do things like nan + 1 = nan, so I guess it makes sense
        #
        if ($_[0] =~ /^[-+]*(Inf|NaN)/i)
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
    $str_a =~ s/^Global(Scale|FitLine):\s+//ig;
    $str_b =~ s/^Global(Scale|FitLine):\s+//ig;
    
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


sub cmp_orig_sample_order
{
    my $col_a;
    my $col_b;

    $col_a = $sample_to_file_col_hash{$a};
    $col_b = $sample_to_file_col_hash{$b};

    if ($col_a < $col_b) { return -1; }
    if ($col_a > $col_b) { return  1; }

    return ($a cmp $b);
}


sub read_in_data_file
{
    $infile = $_[0];
    open INFILE, "$infile" or die "ABORT -- can't open $infile\n";

    # read in header line
    $line = <INFILE>;
    $line =~ s/[\r\n]+//g;
    $line =~ s/\"//;

    # remove UTF-8 byte order mark, since it causes many problems
    # remove some other BOM I've observed in the wild
    $line =~ s/(^\xEF\xBB\xBF|^\xEF\x3E\x3E\xBF)//;

    # do NOT strip off trailing empty fields !!
    @array = split /\t/, $line, -1;

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

        # do NOT strip off trailing empty fields !!
        @array = split /\t/, $line, -1;

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
    close INFILE;

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


# file should have no headers, just a list of samples
sub read_in_comp_pool_exclusions_file
{
    $infile = $_[0];

    open INFILE, "$infile" or die "ABORT -- can't open $infile\n";
    
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
        
        # Allow multiple columns of exclusions
        # Since we may have different sample names between pY, pSTY, etc.
        # Assuming different samples don't share the same name across
        #  assays, we can put them all in a single file
        foreach $sample (@array)
        {
            if ($sample =~ /\S/)
            {
                #printf STDERR "EXCLUDE0\t%s\n",
                #    $sample;

                $comp_pool_exclude_hash{$sample} = 1;
            }
        }
    }
    close INFILE;
}


# scan through the stored log2 scaling factors for each plex
# flag any sample that is much worse than the median for that plex as bad
sub auto_flag_failed_samples
{
    for ($p = 0; $p < $num_plexes; $p++)
    {
        @temp_array     = ();
        $median_log2_sf = '';

        for ($ch = 0; $ch < $num_channels; $ch++)
        {
            $log2_sf         = $log2_sf_array[$p][$ch];
            $temp_array[$ch] = $log2_sf;
        }

        @temp_array = sort {$a <=> $b} @temp_array;
    
        $num_ch_half = $num_channels >> 1;
        # odd
        if ($num_channels % 2)
        {
            $median_log2_sf = $temp_array[$num_ch_half];
        }
        # even
        else
        {
            $median_log2_sf = 0.5 *
                              ($temp_array[$num_ch_half - 1] +
                               $temp_array[$num_ch_half]);
        }

        $tmt_plex = $tmt_plex_array[$p];

        for ($ch = 0; $ch < $num_channels; $ch++)
        {
            $sample  = $tmt_plex_hash{$tmt_plex}{$channel_array[$ch]};
            $log2_sf = $log2_sf_array[$p][$ch];
            $delta   = $log2_sf - $median_log2_sf;
            
            # fudge it a little, for round off error
            if ($delta > $log2_sf_above_median_cutoff - 1E-5)
            {
                #printf STDERR "Auto-exclude dark channel from comp pool:\t%s\n",
                #    $sample;

                $comp_pool_exclude_hash{$sample} = 1;
            }
        }
    }
}


# exclude the last and 2nd to last channels
sub flag_boosting_channels
{
    @temp_array = sort cmp_orig_sample_order @sample_array;

    for ($p = 0; $p < $num_plexes; $p++)
    {
        $tmt_plex = $tmt_plex_array[$p];

        $sample = $temp_array[($p * $num_channels) - 1];
        #printf STDERR "Auto-exclude boost channel from comp pool:\t%s\n",
        #    $sample;
        $comp_pool_exclude_hash{$sample} = 1;

        $sample = $temp_array[($p * $num_channels) - 3];
        #printf STDERR "Auto-exclude boost channel from comp pool:\t%s\n",
        #    $sample;
        $comp_pool_exclude_hash{$sample} = 1;
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
    
    # support auto identification of 2 pool channels
    # "auto2" is now what used to be "variable" mode
    if (defined($fixed_pool_hash{'auto2'}))
    {
        @fixed_pool_array               = ();
        $num_fixed_pools                = 0;
        $auto_single_variable_pool_flag = 0;
    }
    # automatically pick the best sample per plex, can be different for each
    if (defined($fixed_pool_hash{'auto'}) ||
        defined($fixed_pool_hash{'auto1'}))
    {
        @fixed_pool_array               = ();
        $num_fixed_pools                = 0;
        $auto_single_variable_pool_flag = 1;
    }

    # fixed pool channel(s)
    $num_pools = 1;             # also for single auto-identified ref sample
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
    if ($num_fixed_pools == 0 && $auto_single_variable_pool_flag == 0)
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
    
    # else, automatically pick the best sample per plex, can be different
    # put it off for later, within the sample normalization routine
}


# print warnings to STDERR if the pools look too different
sub check_outlier_pools
{
    if ($num_fixed_pools == 1 || $auto_single_variable_pool_flag)
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
                    
                    if (defined($comp_pool_exclude_hash{$sample}))
                    {
                        #printf STDERR "EXCLUDE1\t%s\n",
                        #    $sample;

                        next;
                    }
                
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
          $cmd_string = sprintf "iron_generic --proteomics --iron-exclusions=\"%s\" --norm-iron=\"%s\" \"%s\" -o \"%s\" 2>&1| grep -P \"^Global(Scale|Fit)\"",
              $exclusions_filename, $median_sample, $iron_input_name, $iron_output_name;
        }
        else
        {
          $cmd_string = sprintf "iron_generic --proteomics --norm-iron=\"%s\" \"%s\" -o \"%s\" 2>&1| grep -P \"^Global(Scale|Fit)\"",
              $median_sample, $iron_input_name, $iron_output_name;
        }
        
        if ($iron_untilt_flag)
        {
            $cmd_string =~ s/--proteomics/--rnaseq/g;
        }

        `$cmd_string`;
        # `cp -a $iron_output_name debug_iron_pools.txt`;

        @scale_lines_unsorted = `$cmd_string`;
        @scale_lines = sort cmp_scale_lines @scale_lines_unsorted;

        # output scaling factors file
        $log2 = log(2.0);

        if ($scales_output_init == 0)
        {
            open OUTPUT_SCALING, ">$scales_output_name" or die "ABORT -- can't open scales output file $scales_output_name\n";
        }
        else
        {
            open OUTPUT_SCALING, ">>$scales_output_name" or die "ABORT -- can't open scales output file $scales_output_name\n";
        }

        for ($i = 0; $i < @scale_lines; $i++)
        {
            $line = $scale_lines[$i];
            $line =~ s/[\r\n]+//g;

            # strip GlobalScale: from beginning
            $line =~ s/^Global(Scale|FitLine):\s+//ig;

            # strip Fit from headers
            # for compatability with older scripts
            if ($i == 0)
            {
                $line =~ s/Fit//g;
            }

            # do NOT strip off trailing empty fields !!
            @array = split /\t/, $line, -1;
            
            # only print the header line once
            if ($i == 0)
            {
                if ($scales_output_init == 0)
                {
                    $scales_output_init = 1;
                    
                    printf OUTPUT_SCALING "%s",   $array[0];
                    printf OUTPUT_SCALING "\t%s", 'Plex';
                    for ($j = 1; $j < @array; $j++)
                    {
                        printf OUTPUT_SCALING "\t%s", $array[$j];
                    }
                    printf OUTPUT_SCALING "\n";
                    

                    # store which column has the log2 scaling factor
                    for ($j = 0; $j < @array; $j++)
                    {
                        if ($array[$j] =~ /Log2/i &&
                            $array[$j] =~ /Scale/i)
                        {
                            $global_log2_sf_col = $j;
                        }
                    }
                }
            
                next;
            }
            
            printf OUTPUT_SCALING "%s",   $array[0];
            printf OUTPUT_SCALING "\t%s", 'PoolNorm';
            for ($j = 1; $j < @array; $j++)
            {
                printf OUTPUT_SCALING "\t%s", $array[$j];
            }
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

        # do NOT strip off trailing empty fields !!
        @array = split /\t/, $line, -1;

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
        $tmt_plex     = $tmt_plex_array[$p];
        %temp_ch_hash = ();

        if ($auto_single_variable_pool_flag == 0)
        {
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
        }
        
        # open IRON input file for writing
        open INPUT_FOR_IRON, ">$iron_input_name" or die "ABORT -- can't open $iron_input_name\n";

        # print header
        printf INPUT_FOR_IRON "%s", 'ProbeID';
        for ($ch = 0; $ch < $num_channels; $ch++)
        {
            $sample = $tmt_plex_hash{$tmt_plex}{$channel_array[$ch]};
            printf INPUT_FOR_IRON "\t%s", $sample;
            
            $temp_ch_hash{$sample} = $ch;
        }
        if ($auto_single_variable_pool_flag == 0)
        {
            printf INPUT_FOR_IRON "\t%s", '__AvgPool__';

            # $temp_ch_hash{'__AvgPool__'} = $num_channels;
        }
        printf INPUT_FOR_IRON "\n";

        for ($row = 0; $row < $num_rows; $row++)
        {
            $identifier = $orig_data_array[$row][$id_col];
            
            if ($auto_single_variable_pool_flag == 0)
            {
                $value1 = $condensed_data_array[$row][$pool_s1];
                $value2 = $condensed_data_array[$row][$pool_s2];

                $avg = '';
                if (defined($value1) && defined($value2))
                {
#                    $avg = 0.5 * ($value1 + $value2);
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
            if ($auto_single_variable_pool_flag == 0)
            {
                printf INPUT_FOR_IRON "\t%s", $avg;
            }
            
            printf INPUT_FOR_IRON "\n";
        }
        close INPUT_FOR_IRON;

        # IRON the samples
        if ($no_iron_flag == 0)
        {
            # automatically pick the best sample per plex, can be different
            if ($auto_single_variable_pool_flag)
            {
                # find the median sample
                $cmd_string = sprintf "findmedian --pearson --spreadsheet \"%s\" 2>/dev/null | tail -2 | head -1 | cut -f 4",
                    $iron_input_name;
                $median_sample = `$cmd_string`;
                $median_sample =~ s/[\r\n]+//;
                
                $pool_ch1   = $temp_ch_hash{$median_sample};
                $pool_samp1 = $tmt_plex_hash{$tmt_plex}{$channel_array[$pool_ch1]};
                $pool_s1    = $sample_to_condensed_col_hash{$pool_samp1};
                
                $plex_pool_channels[$p][0] = $pool_ch1;
                $plex_pool_channels[$p][1] = '';
                $plex_pool_channels[$p][2] = '';

                printf STDERR "POOL\t%s\t%s\n", $tmt_plex, $pool_samp1;
            }
            # use the fixed or averaged pools
            else
            {
                $median_sample = '__AvgPool__';
            }

            # run IRON
            printf STDERR "Running IRON:\t%s\n", $tmt_plex;
#   #       printf STDERR "Median sample %s:\t%s\n", $filename, $median_sample;

            if ($exclusions_flag)
            {
              $cmd_string = sprintf "iron_generic --proteomics --iron-exclusions=\"%s\" --norm-iron=\"%s\" \"%s\" -o \"%s\" 2>&1| grep -P \"^Global(Scale|FitLine)\"",
                  $exclusions_filename, $median_sample, $iron_input_name, $iron_output_name;
            }
            else
            {
              $cmd_string = sprintf "iron_generic --proteomics --norm-iron=\"%s\" \"%s\" -o \"%s\" 2>&1| grep -P \"^Global(Scale|FitLine)\"",
                  $median_sample, $iron_input_name, $iron_output_name;
            }

            if ($iron_untilt_flag)
            {
                $cmd_string =~ s/--proteomics/--rnaseq/g;
            }

            @scale_lines_unsorted = `$cmd_string`;
            @scale_lines = sort cmp_scale_lines @scale_lines_unsorted;

            # output scaling factors file
            $log2 = log(2.0);
            if ($scales_output_init == 0)
            {
                open OUTPUT_SCALING, ">$scales_output_name" or die "ABORT -- can't open scales output file $scales_output_name\n";
            }
            else
            {
                open OUTPUT_SCALING, ">>$scales_output_name" or die "ABORT -- can't open scales output file $scales_output_name\n";
            }

            # output scaling factors file
            for ($i = 0; $i < @scale_lines; $i++)
            {
                $line = $scale_lines[$i];
                $line =~ s/[\r\n]+//g;

                # strip GlobalScale: from beginning
                $line =~ s/^Global(Scale|FitLine):\s+//ig;

                # strip Fit from headers
                # for compatability with older scripts
                if ($i == 0)
                {
                    $line =~ s/Fit//g;
                }

                # do NOT strip off trailing empty fields !!
                @array = split /\t/, $line, -1;

                # only print the header line once
                if ($i == 0)
                {
                    if ($scales_output_init == 0)
                    {
                        $scales_output_init = 1;
                        
                        printf OUTPUT_SCALING "%s",   $array[0];
                        printf OUTPUT_SCALING "\t%s", 'Plex';
                        for ($j = 1; $j < @array; $j++)
                        {
                            printf OUTPUT_SCALING "\t%s", $array[$j];
                        }
                        printf OUTPUT_SCALING "\n";


                        # store which column has the log2 scaling factor
                        for ($j = 0; $j < @array; $j++)
                        {
                            if ($array[$j] =~ /Log2/i &&
                                $array[$j] =~ /Scale/i)
                            {
                                $global_log2_sf_col = $j;
                            }
                        }
                    }
                
                    next;
                }
                
                # skip average pool
                if ($line =~ /^__AvgPool__/)
                {
                    next;
                }
                
                printf OUTPUT_SCALING "%s",   $array[0];
                printf OUTPUT_SCALING "\t%s", $tmt_plex;
                for ($j = 1; $j < @array; $j++)
                {
                    printf OUTPUT_SCALING "\t%s", $array[$j];
                }
                printf OUTPUT_SCALING "\n";
            }
            close OUTPUT_SCALING;


            # store scaling factors
            # must use original unsorted lines, since that lines up with
            #  all the other bookkeeping
            $k = 0;
            for ($i = 1; $i < @scale_lines; $i++)
            {
                $line = $scale_lines_unsorted[$i];
                $line =~ s/[\r\n]+//g;

                # strip GlobalScale: from beginning
                $line =~ s/^Global(Scale|FitLine):\s+//ig;

                # skip average pool
                if ($line =~ /^__AvgPool__/)
                {
                    next;
                }

                # do NOT strip off trailing empty fields !!
                @array = split /\t/, $line, -1;
                
                # store the log2 scaling factors for later
                $log2_sf               = $array[$global_log2_sf_col];
                $log2_sf_array[$p][$k] = $log2_sf;
                $k++;
            }
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

            # do NOT strip off trailing empty fields !!
            @array = split /\t/, $line, -1;

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

                    if (defined($comp_pool_exclude_hash{$sample}))
                    {
                        #printf STDERR "EXCLUDE2\t%s\n",
                        #    $sample;

                        next;
                    }
                
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
    # do NOT strip off trailing empty fields !!
    @header_array = split /\t/, $header_line, -1;
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
$auto_single_variable_pool_flag = 0;
$unlog2_flag       = 0;    # unlog2 the input data
$no_log2_flag      = 0;    # do not log2 the output data
$comp_pool_flag               = 0;    # use comp pool instead of real pool
$comp_pool_exclude_flag       = 0;
$comp_pool_exclude_filename   = '';
$comp_pool_exclude_dark_flag  = 0;
$comp_pool_exclude_boost_flag = 0;
$iron_untilt_flag  = 0;    # --rnaseq flag

$error_flag = 0;
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
        elsif ($field =~ /^--comp-pool-exclusions\=(\S+)/)
        {
            # $comp_pool_flag = 1;
            $comp_pool_exclude_flag = 1;
            $comp_pool_exclude_filename = $1;

            printf STDERR "Using comp pool sample exclusion file: %s\n",
                          $comp_pool_exclude_filename;
        }
        elsif ($field =~ /^--comp-pool-exclusions-dark$/)
        {
            # $comp_pool_flag = 1;
            $comp_pool_exclude_flag = 1;
            $comp_pool_exclude_filename = '';
            $comp_pool_exclude_dark_flag = 1;

            printf STDERR "Auto-detecting comp pool dark channel exclusions\n";
        }
        elsif ($field =~ /^--comp-pool-exclusions-boost$/)
        {
            # $comp_pool_flag = 1;
            $comp_pool_exclude_flag = 1;
            $comp_pool_exclude_boost_flag = 1;
            $comp_pool_exclude_filename = '';

            printf STDERR "Excluding boost and boost-blank channels from comp pool\n",
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
        elsif ($field =~ /^--iron-untilt$/)
        {
            $iron_untilt_flag = 1;
            $no_iron_flag     = 0;
        }
        else
        {
            printf "ABORT -- unknown option %s\n", $field;

            $error_flag = 1;
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


if (!defined($filename))
{
    printf STDERR "ABORT -- no input file specified\n";

    $error_flag = 1;
}

if ($error_flag)
{
    print STDERR "Usage: automate_tmt.pl [options] maxquant_output_file.txt [IRON ref channels] > iron_output.txt\n";
    print STDERR "\n";
    print STDERR "  Options:\n";
    print STDERR "    --iron             normalize within-plex prior to other calcuations [default]\n";
    print STDERR "    --iron-untilt      account for relative dynamic range\n";
    print STDERR "    --no-iron          do not normalize within each plex first\n";
    print STDERR "    --debatch          use reference channel for cross-plex normalization [default]\n";
    print STDERR "    --no-debatch       do not perform cross-plex normalization\n";
    print STDERR "    --leave-ratios     leave cross-plex normalized data as log2 ratios\n";
    print STDERR "    --no-leave-ratios  scale cross-plex normalized log2 ratios back into abundances [default]\n";
    print STDERR "    --comp-pool        use all-channel geometric mean for cross-plex debatching\n";
    print STDERR "    --no-comp-pool     do not create a computational reference pool for cross-plex debatching [default]\n";
    print STDERR "    --comp-pool-exclusions-dark\n";
    print STDERR "                       auto-excludes dark samples from computational pool\n";
    print STDERR "    --comp-pool-exclusions-boost\n";
    print STDERR "                       excludes boosting channels (N, N-2) from comp pool\n";
    print STDERR "    --comp-pool-exclusions=filename.txt\n";
    print STDERR "                       load comp pool sample exclusions from tab-delimited file\n";
    print STDERR "    --iron-exclusions=filename.txt\n";
    print STDERR "                       exclude row identifiers from IRON training\n";
    print STDERR "\n";
    print STDERR "    --no-iron --no-debatch will leave the output abundances unchanged\n";
    print STDERR "\n";
    print STDERR "  Output:\n";
    print STDERR "    Normalized data is output to STDOUT.\n";
    print STDERR "\n";
    print STDERR "    Several other files are generated, extending the original input filename:\n";
    print STDERR "\n";
    print STDERR "      *_scaling_factors.txt  various statistics from the IRON normalizations\n";
    print STDERR "      *_findmedian.txt       results from the various findmedian runs\n";

    exit(1);
}


$process_id = $$;

$scales_output_init =  0;
$scales_output_name =  $filename;
$scales_output_name =~ s/\.txt$//i;
$scales_output_name =  sprintf "%s_scaling_factors.txt", $scales_output_name;
`rm \"$scales_output_name\" >& /dev/null`;

$row_scales_output_name =  $filename;
$row_scales_output_name =~ s/\.txt$//i;
$row_scales_output_name = sprintf "%s_row_scales.txt", $row_scales_output_name;
`rm \"$row_scales_output_name\" >& /dev/null`;

$iron_input_name  = sprintf "temp_iron_input_%s.txt", $process_id;
$iron_output_name = sprintf "temp_iron_output_%s.txt", $process_id;

# default to using TMT-126C as the normalization channel
@fixed_pool_array = sort keys %fixed_pool_hash;
if (@fixed_pool_array == 0)
{
    $fixed_pool_hash{'TMT-126C'} = 1;
}

%comp_pool_exclude_hash = ();

read_in_data_file($filename);
store_condensed_data();

if ($comp_pool_flag && $comp_pool_exclude_flag &&
    $comp_pool_exclude_filename ne '')
{
    read_in_comp_pool_exclusions_file($comp_pool_exclude_filename);
}

#normalize_crude();
identify_pools();
check_outlier_pools();

# must come before iron_pools() now, to support auto ref sample picking
iron_samples();

# scan scaling factors for samples to auto-exclude from comp pool
# must come after iron_samples() and before iron_pools()
if ($comp_pool_flag && $comp_pool_exclude_flag &&
    $comp_pool_exclude_dark_flag)
{
    auto_flag_failed_samples();
}
if ($comp_pool_flag && $comp_pool_exclude_flag &&
    $comp_pool_exclude_boost_flag)
{
    flag_boosting_channels();
}
foreach $sample (sort keys %comp_pool_exclude_hash)
{
    # only print excluded samples that exist in this file
    if (defined($sample_to_file_col_hash{$sample}))
    {
        printf STDERR "Excluding sample from comp pool:\t%s\n",
            $sample;
    }
}



# must come after iron_samples() now, to support auto ref sample picking
iron_pools();


if ($no_debatch_flag == 0)
{
    correct_abundances();
}

output_final_data();

`rm \"$iron_input_name\"`;
`rm \"$iron_output_name\"`;
