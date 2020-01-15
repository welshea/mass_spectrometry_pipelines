#!/usr/bin/perl -w

# treat any options after the filename as the pool channels
# treat them as free text, don't try to get smart,
# the user will need to specify the channel names as they rae
# in the input data file
#
# Don't forget that current file format is ex: TMT-126, not just 126

# Syntax: rollup_tmt.pl [options] log2_normalized_data.txt regrouped_peptides.txt [pool channel]
#  data must have already run through the IRON normalization pipeline
#  abundances are in log2 units already

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


sub read_in_data_file
{
    $infile = $_[0];
    open INFILE, "$infile" or die "ABORT -- can't open $infile\n";

    # read in header line
    $line = <INFILE>;
    $line =~ s/[\r\n]+//g;
    $line =~ s/\"//;

    @array = split /\t/, $line;

    $num_samples = 0;
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\s+/ /g;
        
        if ($array[$i] eq 'RowIdentifier')
        {
            $id_col = $i;
        }
        elsif ($array[$i] eq 'Sequence')
        {
            $data_seq_col = $i;
        }
        else
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
    }
    $num_header_cols = @array;
    $header_line = join "\t", @array;
    
    if (!defined($data_seq_col))
    {
        printf "ABORT -- cannot find Sequence column in data file %s\n",
            $infile;
        exit(1);
    }

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
            $orig_data_array[$row][$i] = $array[$i];
        }
        
        $seq                      = $array[$data_seq_col];
        $seq_row_hash{$seq}{$row} = 1;
#        $row_seq_array[$row]      = $seq;
        
        $row++;
    }
    $num_rows = $row;
    
    close INFILE;
}


sub read_in_group_file
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
        
        $group_header_col_hash{$array[$i]} = $i;
    }
    
    $group_seq_col = $group_header_col_hash{'Sequence'};
    
    if ($use_maxquant_unique_flag)
    {
        $group_group_col = $group_header_col_hash{'Group_Maxquant'};
    }
    elsif ($use_maxquant_razor_flag)
    {
        $group_group_col = $group_header_col_hash{'Group_Maxquant_Razor'};
    }
    else
    {
        $group_group_col = $group_header_col_hash{'Group_Regrouped'};
    }
    
    # we might be reading a filtered summarized file
    if (!defined($group_group_col))
    {
        $group_group_col = $group_header_col_hash{'ProteinGroup'};
    }

    
    if (!defined($group_seq_col))
    {
        printf "ABORT -- cannot find Sequence column in data file %s\n",
            $infile;
        exit(1);
    }
    if (!defined($group_group_col))
    {
        printf "ABORT -- cannot find Group_* column in data file %s\n",
            $infile;
        exit(1);
    }
    
    %group_seq_hash = ();
    %group_row_hash = ();

    # read in the groups
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
        
        $seq       = $array[$group_seq_col];
        $group_str = $array[$group_group_col];
        
        # WARNING -- maxquant can assign a peptide to multiple groups
        $group_seq_hash{$group_str}{$seq} = 1;
    }
    close INFILE;
    
    # map rows to groups
    foreach $group (keys %group_seq_hash)
    {
        foreach $seq (keys %{$group_seq_hash{$group}})
        {
            foreach $row (keys %{$seq_row_hash{$seq}})
            {
                $group_row_hash{$group}{$row} = 1;
            }
        }
    }
}


sub read_in_group_summary_file
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
        
        $summary_header_col_hash{$array[$i]} = $i;
    }
    
    $summary_group_col = $summary_header_col_hash{'ProteinGroup'};
    
    if (!defined($summary_group_col))
    {
        printf "ABORT -- cannot find ProteinGroup column in group summary file %s\n",
            $infile;
        exit(1);
    }

    # assemble header line
    @array_new = ();
    $j = 0;
    for ($i = 0; $i < @array; $i++)
    {
        if ($i != $summary_group_col)
        {
            $array_new[$j++] = $array[$i];
        }
    }
    $summary_header_line = join "\t", @array_new;

    # read in the groups
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
        
        $group = $array[$summary_group_col];
        
        if ($group =~ /\S/)
        {
            # assemble new line
            @array_new = ();
            $j = 0;
            for ($i = 0; $i < @array; $i++)
            {
                if ($i != $summary_group_col)
                {
                    $array_new[$j++] = $array[$i];
                }
            }
            $line_new = join "\t", @array_new;

            $group_annotation_hash{$group} = $line_new;
        }
    }
    close INFILE;
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



# scan through the 6-plexes to automatically identify the pool channels
# assume the two most similar channels are the pools
sub identify_pools
{
    @fixed_pool_array = sort keys %fixed_pool_hash;
    $num_fixed_pools = @fixed_pool_array;
    
    # support "auto", since default is now TMT-126
    if (defined($fixed_pool_hash{'auto'}))
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



# normalize pan-pool cohort for estimating pseudo-abundance scales
sub average_pools
{
    @row_pool_avg_values = ();

    for ($row = 0; $row < $num_rows; $row++)
    {
        @value_array = ();
        $avg         = 0;
        $sd          = 0;
        $count       = 0;
    
        for ($p = 0; $p < $num_plexes; $p++)
        {
            $tmt_plex = $tmt_plex_array[$p];

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
                
                if (is_number($value))
                {
                    $value_array[$count++] = $value;
                }
            }
        }
        
        for ($i = 0; $i < $count; $i++)
        {
            $value = $value_array[$i];
            $avg += $value;
        }

        if ($count)
        {
            $avg /= $count;
        }

        for ($i = 0; $i < $count; $i++)
        {
            $value = $value_array[$i];
            $diff  = $avg - $value;
            $sd   += $diff * $diff;
        }
        if ($count)
        {
            $sd = sqrt($sd / $count);
        }
        
        if ($sd > 1e-5)
        {
            printf STDERR "ABORT -- pools not identical for row %s %s %s\n",
                $row, $avg, $sd;
            
            exit(2);
        }
        
        $row_pool_avg_values[$row] = $avg;
    }
}



sub rollup_pool
{
    @group_array = sort keys %group_row_hash;
    
    foreach $group (@group_array)
    {
        $avg   = 0;
        $count = 0;

        foreach $row (keys %{$group_row_hash{$group}})
        {
            $value = $row_pool_avg_values[$row];
            
            $avg += $value;
            $count++;
        }
        
        if ($count)
        {
            $avg /= $count;
        }
        
        if ($count)
        {
            $group_pool_rollup_hash{$group} = $avg;
        }
    }
}



sub rollup_samples
{
    for ($p = 0; $p < $num_plexes; $p++)
    {
        $tmt_plex = $tmt_plex_array[$p];

        for ($ch = 0; $ch < $num_channels; $ch++)
        {
            $sample = $tmt_plex_hash{$tmt_plex}{$channel_array[$ch]};
            $col    = $sample_to_condensed_col_hash{$sample};
            
            foreach $group (@group_array)
            {
                $avg_log2_ratio = 0;
                $count          = 0;
            
                foreach $row (keys %{$group_row_hash{$group}})
                {
                    # all pools should be identical, assume that they are
                    # if they aren't, average_pools() will have already
                    #  issues warnings to this effect
                    $pool_avg = $row_pool_avg_values[$row];

                    $value = $condensed_data_array[$row][$col];

                    if (defined($value) && defined($pool_avg))
                    {
                        $avg_log2_ratio += $value - $pool_avg;
                        $count++;
                        
#                        printf "FOOBAR\t%s\t%s\t%s\t%s\t%s\n",
#                            $group, $p, $col, $row, $value - $pool_avg;
                    }
                }
                
                if ($count)
                {
                    $group_pool_rollup = $group_pool_rollup_hash{$group};
                
                    $avg_log2_ratio /= $count;
                
                    $group_col_rollup_hash{$group}[$col] =
                        $group_pool_rollup + $avg_log2_ratio;
                    
#                    printf "FOOBAR\t%s\t%s\t%s\n",
#                        $group, $group_pool_rollup, $avg_log2_ratio;
                }
            }
        }
    }
}



sub output_rollup
{
    @array = split /\t/, $summary_header_line;
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] = '';
    }
    $group_annotation_blank = join "\t", @array;
    $num_summary_fields = @array;


    # print headers
    printf "%s",   'ProteinGroup';
    printf "\t%s", 'NumPeptidesObserved';
    printf "\t%s", $summary_header_line;
    
    for ($p = 0; $p < $num_plexes; $p++)
    {
        $tmt_plex = $tmt_plex_array[$p];

        for ($ch = 0; $ch < $num_channels; $ch++)
        {
            $sample = $tmt_plex_hash{$tmt_plex}{$channel_array[$ch]};
            
            printf "\t%s", $sample;
        }
    }
    printf "\n";


    foreach $group (@group_array)
    {
        @row_array = keys %{$group_row_hash{$group}};
        
        $num_peptides = @row_array;
        
        printf "%s",   $group;
        printf "\t%d", $num_peptides;
        
        $group_annotation = $group_annotation_hash{$group};
        
        # missing group annotation
        if (!defined($group_annotation))
        {
            printf STDERR "WARNING -- missing %s group annotation\n", $group;
            
            printf "\t%s", $group_annotation_blank;
        }
        else
        {
            @array = split /\t/, $group_annotation;
            
            if (@array ne $num_summary_fields)
            {
                @array_new = '';
                $j = 0;
                # make sure we truncate any lines that are too long
                for ($i = 0; $i < $num_summary_fields; $i++)
                {
                    $array_new[$j++] = $array_new[$i];
                }
                # make sure we pad any lines that are too short
                for (; $i < $num_summary_fields; $i++)
                {
                    $array_new[$j++] = '';
                }
                
                $group_annotation = join "\t", @array_new;
            }
            
            printf "\t%s", $group_annotation;
        }

        for ($p = 0; $p < $num_plexes; $p++)
        {
            $tmt_plex = $tmt_plex_array[$p];

            for ($ch = 0; $ch < $num_channels; $ch++)
            {
                $sample = $tmt_plex_hash{$tmt_plex}{$channel_array[$ch]};
                $col    = $sample_to_condensed_col_hash{$sample};
                
                $value = '';
                if (defined($group_col_rollup_hash{$group}))
                {
                    $value = $group_col_rollup_hash{$group}[$col];
                    
                    if (!defined($value))
                    {
                        $value = '';
                    }
                }
                
                printf "\t%s", $value;
            }
        }
        
        printf "\n";
    }
}


sub output_original_data
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
                $orig_data_array[$row][$col] = $value;
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
}


# read in command line arguments
$num_files = 0;
$num_fixed_pools = 0;
$syntax_error_flag = 0;
$use_maxquant_unique_flag = 0;
$use_maxquant_razor_flag  = 0;

for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        if ($field eq '--regroup' ||
            $field eq '--regrouped')
        {
            $use_maxquant_razor_flag  = 0;
            $use_maxquant_unique_flag = 0;
        }
        elsif ($field eq '--maxquant-razor')
        {
            $use_maxquant_razor_flag = 1;
        }
        elsif ($field eq '--maxquant-unique')
        {
            $use_maxquant_unique_flag = 1;
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
        elsif ($num_files == 1)
        {
            $grouped_peptides_filename = $field;
            $num_files++;
        }
        elsif ($num_files == 2)
        {
            $group_summary_filename = $field;
            $num_files++;
        }
        # treat any options after the filename as the pool channels
        # treat them as free text, don't try to get smart,
        # the user will need to specify the channel names as they rae
        # in the input data file
        #
        # Don't forget that current file format is ex: TMT-126, not just 126
        else
        {
            $fixed_pool_hash{$field} = 1;
        }
    }
}


# default to using TMT-126 as the normalization channel
if (!defined(%fixed_pool_hash))
{
    $fixed_pool_hash{'TMT-126'} = 1;
}

read_in_data_file($filename);
read_in_group_file($grouped_peptides_filename);
read_in_group_summary_file($group_summary_filename);
store_condensed_data();
identify_pools();
average_pools();
#iron_pools();
#iron_samples();
rollup_pool();
rollup_samples();

#output_original_data();
output_rollup();
