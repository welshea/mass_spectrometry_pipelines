#!/usr/bin/perl -w


# 2022-06-07:  add support for single input file
# 2022-04-11:  attempt to auto-prepend pos/neg auto-prepended during merging
# 2022-03-10:  changelog edits
# 2022-03-09:  re-score after initial scoring to update median sample
# 2021-08-19:  bugfix: print some errors to STDERR instead of STDOUT
# 2021-08-13:  print samples sorted by worst score first
# 2021-08-13:  fix typo in median sample indexing for odd N


use Scalar::Util qw(looks_like_number);
use POSIX;
use File::Basename;



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


# use on array of merged sample names
sub cmp_sample
{
    my $sample_a = $a;
    my $sample_b = $b;
    my $value_a;
    my $value_b;
    
    # HACK -- sort on score (worst first), if they exit at this point
    #
    # We want to bring the worst samples to the top before printing
    #
    $value_a = $sample_merged_hash{$sample_a}{merged}{score};
    $value_b = $sample_merged_hash{$sample_b}{merged}{score};
    if (defined($value_a) && defined($value_b))
    {
        if ($value_a < $value_b) { return -1; }
        if ($value_a > $value_b) { return  1; }
    }

    $value_a = $sample_merged_hash{$sample_a}{merged}{log2scale};
    $value_b = $sample_merged_hash{$sample_b}{merged}{log2scale};
    if ($value_a < $value_b) { return -1; }
    if ($value_a > $value_b) { return  1; }
    
    $value_a = $sample_merged_hash{$sample_a}{merged}{mean_dist};
    $value_b = $sample_merged_hash{$sample_b}{merged}{mean_dist};
    if ($value_a < $value_b) { return -1; }
    if ($value_a > $value_b) { return  1; }

    $value_a = $sample_merged_hash{$sample_a}{merged}{present};
    $value_b = $sample_merged_hash{$sample_b}{merged}{present};
    if ($value_a > $value_b) { return -1; }
    if ($value_a < $value_b) { return  1; }

    $value_a = $sample_merged_hash{$sample_a}{merged}{mean_log2};
    $value_b = $sample_merged_hash{$sample_b}{merged}{mean_log2};
    if ($value_a > $value_b) { return -1; }
    if ($value_a < $value_b) { return  1; }
    
    return $sample_a cmp $sample_b;
}


# call after reading in the sample name table file
sub read_in_scaling_factors_file
{
    my $filename = $_[0];
    my $pos_neg  = $_[1];
    my $line;
    my @array = ();
    my @header_col_array  = ();
    my %header_col_hash   = ();
    my $sample;
    my $header;
    my $field;

    my $sampleid_col;
    my $log2scale_col;
    my $present_sample_col;
    my $present_dataset_col;
    my $sampleid;
    my $sample_merged;
    my $log2scale;
    my $present_sample;
    my $present_dataset;

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
        }
    }
 
    $sampleid_col        = $header_col_hash{'SampleID'};
    $log2scale_col       = $header_col_hash{'Log2Scale'};
    $present_sample_col  = $header_col_hash{'PresentSample'};
    $present_dataset_col = $header_col_hash{'PresentDataset'};
    
    if (!defined($sampleid_col))
    {
        printf STDERR "ABORT -- scaling factors SampleID not found in file %s\n",
            $filename;
        exit(2);
    }
    if (!defined($log2scale_col))
    {
        printf STDERR "ABORT -- scaling factors Log2Scale not found in file %s\n",
            $filename;
        exit(2);
    }
    if (!defined($present_sample_col))
    {
        printf STDERR "ABORT -- scaling factors PresentSample not found in file %s\n",
            $filename;
        exit(2);
    }
    if (!defined($present_dataset_col))
    {
        printf STDERR "ABORT -- scaling factors PresentDataset not found in file %s\n",
            $filename;
        exit(2);
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

        $sampleid        = $array[$sampleid_col];
        $log2scale       = $array[$log2scale_col];
        $present_sample  = $array[$present_sample_col];
        $present_dataset = $array[$present_dataset_col];
        
        $global_present_dataset = $present_dataset;

        
        if ($pos_neg =~ /pos/i)
        {
            $sample_merged = $global_pos_map_hash{$sampleid};
            
            # maybe pos/neg were auto-prepended elsewhere
            if (!defined($sample_merged))
            {
                $sample_merged = $sampleid;
                $sample_merged = 'pos_' . $sample_merged;
                $sample_merged = $global_pos_map_hash{$sample_merged};
            }

            if (!defined($sample_merged))
            {
                printf STDERR "ABORT -- pos scaling factors %s not mapped to merged sample name\n",
                    $sampleid;

                exit(3);
            }
        }
        if ($pos_neg =~ /neg/i)
        {
            $sample_merged = $global_neg_map_hash{$sampleid};

            # maybe pos/neg were auto-prepended elsewhere
            if (!defined($sample_merged))
            {
                $sample_merged = $sampleid;
                $sample_merged = 'neg_' . $sample_merged;
                $sample_merged = $global_neg_map_hash{$sample_merged};
            }

            if (!defined($sample_merged))
            {
                printf STDERR "ABORT -- neg scaling factors %s not mapped to merged sample name\n",
                    $sampleid;

                exit(3);
            }
        }
        
        if (defined($sample_merged))
        {
            $sample_merged_hash{$sample_merged}{$pos_neg}{log2scale} =
                $log2scale;
            $sample_merged_hash{$sample_merged}{$pos_neg}{present} =
                $present_sample;
        }
    }
    close INFILE;
}


# call after reading in the sample name table file
sub read_in_findmedian_file
{
    my $filename = $_[0];
    my $pos_neg  = $_[1];
    my $line;
    my @array = ();
    my @header_col_array  = ();
    my %header_col_hash   = ();
    my $sample;
    my $header;
    my $field;

    my $sample_col;
    my $mean_dist_col;
    my $mean_log2_col;

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
        }
    }

    $sample_col    = $header_col_hash{'SampleName'};
    $mean_dist_col = $header_col_hash{'MeanDistance'};
    $mean_log2_col = $header_col_hash{'MeanLog2Abundance'};
    
    if (!defined($sample_col))
    {
        printf STDERR "ABORT -- findmedian SampleName not found in file %s\n",
            $filename;
        exit(2);
    }
    if (!defined($mean_dist_col))
    {
        printf STDERR "ABORT -- findmedian MeanDistance not found in file %s\n",
            $filename;
        exit(2);
    }
    if (!defined($mean_log2_col))
    {
        printf STDERR "ABORT -- findmedian MeanLog2Abundance not found in file %s\n",
            $filename;
        exit(2);
    }

    while(defined($line=<INFILE>))
    {
        # skip non-sample rows
        if (!($line =~ /^Score/))
        {
            next;
        }

        $line =~ s/[\r\n]+//g;

        @array = split /\t/, $line;

        for ($i = 0; $i < @array; $i++)
        {
            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;
            $array[$i] =~ s/\s+/ /g;
        }

        $sample    = $array[$sample_col];
        $mean_dist = $array[$mean_dist_col];
        $mean_log2 = $array[$mean_log2_col];

        if ($pos_neg =~ /pos/i)
        {
            $sample_merged = $global_pos_map_hash{$sample};

            # maybe pos/neg were auto-prepended elsewhere
            if (!defined($sample_merged))
            {
                $sample_merged = $sample;
                $sample_merged = 'pos_' . $sample_merged;
                $sample_merged = $global_pos_map_hash{$sample_merged};
            }

            if (!defined($sample_merged))
            {
                printf STDERR "ABORT -- pos findmedian %s not mapped to merged sample name\n",
                    $sample;

                exit(3);
            }
        }
        if ($pos_neg =~ /neg/i)
        {
            $sample_merged = $global_neg_map_hash{$sample};

            # maybe pos/neg were auto-prepended elsewhere
            if (!defined($sample_merged))
            {
                $sample_merged = $sample;
                $sample_merged = 'neg_' . $sample_merged;
                $sample_merged = $global_neg_map_hash{$sample_merged};
            }

            if (!defined($sample_merged))
            {
                printf STDERR "ABORT -- neg findmedian %s not mapped to merged sample name\n",
                    $sample;

                exit(3);
            }
        }
        
        if (defined($sample_merged))
        {
            $sample_merged_hash{$sample_merged}{$pos_neg}{mean_dist} =
                $mean_dist;
            $sample_merged_hash{$sample_merged}{$pos_neg}{mean_log2} =
                $mean_log2;
        }
    }
    close INFILE;
}


# call before reading in pos/neg scaling factor files
sub read_in_sample_table_file
{
    my $filename = $_[0];
    my $line;
    my @array = ();
    my @header_col_array  = ();
    my %header_col_hash   = ();
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
        }
    }
 
    $sample_col  = $header_col_hash{'Sample'};
    $pos_col     = $header_col_hash{'SamplePOS'};
    $neg_col     = $header_col_hash{'SampleNEG'};
    $short_col   = $header_col_hash{'SampleAutoShortened'};
    
    if (!defined($sample_col))
    {
        printf STDERR "ABORT -- sample table Sample not found in file %s\n",
            $filename;
        exit(2);
    }
    if (!defined($pos_col))
    {
        printf STDERR "ABORT -- sample table SamplePOS not found in file %s\n",
            $filename;
        exit(2);
    }
    if (!defined($neg_col))
    {
        printf STDERR "ABORT -- sample table SampleNEG not found in file %s\n",
            $filename;
        exit(2);
    }
    if (!defined($short_col))
    {
        printf STDERR "ABORT -- sample table SampleAutoShortened not found in file %s\n",
            $filename;
        exit(2);
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

        $sample = $array[$sample_col];
        $pos    = $array[$pos_col];
        $neg    = $array[$neg_col];
        $short  = $array[$short_col];

        # map pos/neg sample names to merged sample names
        if ($single_file_mode ne 'neg')
        {
            $global_pos_map_hash{$pos} = $sample;
            $global_sample_map_hash{$sample}{pos}   = $pos;
        }
        if ($single_file_mode ne 'pos')
        {
            $global_neg_map_hash{$neg} = $sample;
            $global_sample_map_hash{$sample}{neg}   = $neg;
        }
        
        # all sample name mappings, mapped from merged sample name
        $global_sample_map_hash{$sample}{short} = $short;
    }
    close INFILE;
}


sub score_qc
{
    my $sample;
    my $sample_pos;
    my $sample_neg;
    my $sample_short;
    my $scale_pos;
    my $scale_neg;
    my $log2scale_pos;
    my $log2scale_neg;
    my $present_pos;
    my $present_neg;
    my $mean_dist_pos;
    my $mean_dist_neg;
    my $mean_log2_pos;
    my $mean_log2_neg;
    my $scale_merged;
    my $log2scale_merged;
    my $present_merged;
    my $mean_dist_merged;
    my $mean_log2_merged;
    my $value_merged;
    my $score;

    my $log2scale_adj;
    my $present_adj;
    my $mean_dist_adj;
    my $mean_log2_adj;
    
    my @sample_array;
    my $best_value;
    my $sample_median;
    my $i;
    my $numerator;
    my $denominator;
    
    @sample_array = sort keys %sample_merged_hash;
    
    # combine pos/neg metrics together
    foreach $sample (@sample_array)
    {
        $log2scale_pos = $sample_merged_hash{$sample}{pos}{log2scale};
        $present_pos   = $sample_merged_hash{$sample}{pos}{present};
        $mean_dist_pos = $sample_merged_hash{$sample}{pos}{mean_dist};
        $mean_log2_pos = $sample_merged_hash{$sample}{pos}{mean_log2};

        $log2scale_neg = $sample_merged_hash{$sample}{neg}{log2scale};
        $present_neg   = $sample_merged_hash{$sample}{neg}{present};
        $mean_dist_neg = $sample_merged_hash{$sample}{neg}{mean_dist};
        $mean_log2_neg = $sample_merged_hash{$sample}{neg}{mean_log2};

        # missing pos sample
        if (!defined($log2scale_pos))
        {
            $log2scale_merged = $log2scale_neg;
            $mean_dist_merged = $mean_dist_neg;
            $mean_log2_merged = $mean_log2_neg;

            # this isn't going to be scaled correctly...
            # leaving out the missing counts would be too small
            # assume roughly equal counts between them, so double it
            $present_merged   = 2 * $present_neg;
        }
        # missing neg sample
        elsif (!defined($log2scale_neg))
        {
            $log2scale_merged = $log2scale_pos;
            $mean_dist_merged = $mean_dist_pos;
            $mean_log2_merged = $mean_log2_pos;

            # this isn't going to be scaled correctly...
            # leaving out the missing counts would be too small
            # assume roughly equal counts between them, so double it
            $present_merged   = 2 * $present_pos;
        }
        # assume we have both pos and neg
        else
        {
            $log2scale_merged = 0.5 * ($log2scale_pos + $log2scale_neg);
            $mean_dist_merged = 0.5 * ($mean_dist_pos + $mean_dist_neg);
            $present_merged   = $present_pos + $present_neg;
            $mean_log2_merged = 0.5 * ($mean_log2_pos + $mean_log2_neg);
        }

        
        $sample_merged_hash{$sample}{merged}{log2scale} = $log2scale_merged;
        $sample_merged_hash{$sample}{merged}{present}   = $present_merged;
        $sample_merged_hash{$sample}{merged}{mean_dist} = $mean_dist_merged;
        $sample_merged_hash{$sample}{merged}{mean_log2} = $mean_log2_merged;
    }

    # sort samples mainly on merged log2scale
    # pick median sample as the actual median sample in the sorted array
    @sample_array = sort cmp_sample @sample_array;
    if (@sample_array % 2)
    {
        $sample_median = $sample_array[@sample_array >> 1];
    }
    # earlier of the two middle samples
    else
    {
        $sample_median = $sample_array[(@sample_array >> 1) - 1];
    }
    
    # adjust merged metrics relative to median sample
    for ($i = 0; $i < @sample_array; $i++)
    {
        $sample = $sample_array[$i];
    
        $sample_merged_hash{$sample}{merged}{log2scale_adj} =
            $sample_merged_hash{$sample}{merged}{log2scale} -
            $sample_merged_hash{$sample_median}{merged}{log2scale};

        $sample_merged_hash{$sample}{merged}{mean_log2_adj} =
            $sample_merged_hash{$sample}{merged}{mean_log2} -
            $sample_merged_hash{$sample_median}{merged}{mean_log2};

        $sample_merged_hash{$sample}{merged}{present_adj}   =
            $sample_merged_hash{$sample}{merged}{present}   /
            $sample_merged_hash{$sample_median}{merged}{present};

        $sample_merged_hash{$sample}{merged}{mean_dist_adj} =
            $sample_merged_hash{$sample}{merged}{mean_dist} /
            $sample_merged_hash{$sample_median}{merged}{mean_dist};
        

        $log2scale_adj = $sample_merged_hash{$sample}{merged}{log2scale_adj};
        $present_adj   = $sample_merged_hash{$sample}{merged}{present_adj};
        $mean_dist_adj = $sample_merged_hash{$sample}{merged}{mean_dist_adj};
        $mean_log2_adj = $sample_merged_hash{$sample}{merged}{mean_log2_adj};
        
        $score = 1.0 / ($mean_dist_adj * pow(2, $log2scale_adj));

        # $score = pow($present_adj / $mean_dist_adj, 1.0 / 2.0) /
        #          pow(2, $log2scale_adj);


        $sample_merged_hash{$sample}{merged}{score} = $score;
    }
}


sub print_qc
{
    my $sample;
    my $sample_pos;
    my $sample_neg;
    my $sample_short;

    my $score;
    my $log2scale_adj;
    my $present_adj;
    my $mean_dist_adj;
    my $mean_log2_adj;

    my $scale_pos;
    my $scale_neg;
    my $log2scale_pos;
    my $log2scale_neg;
    my $present_pos;
    my $present_neg;
    my $mean_dist_pos;
    my $mean_dist_neg;
    my $mean_log2_pos;
    my $mean_log2_neg;

    my @sample_array;
    

    @sample_array = sort cmp_sample keys %sample_merged_hash;

    printf "%s",   'Sample';
    printf "\t%s", 'SamplePOS';
    printf "\t%s", 'SampleNEG';
    printf "\t%s", 'SampleAutoShortened';

    printf "\t%s", 'ScoreRelative';
    printf "\t%s", 'ScaleRelative';
    printf "\t%s", 'Log2ScaleRelative';
    printf "\t%s", 'DistanceRelative';
    printf "\t%s", 'PresentRelative';
    printf "\t%s", 'AbundanceRelative';
    printf "\t%s", 'Log2AbundanceRelative';

    printf "\t%s", 'ScalePOS';
    printf "\t%s", 'ScaleNEG';
    printf "\t%s", 'Log2ScalePOS';
    printf "\t%s", 'Log2ScaleNEG';
    printf "\t%s", 'MeanDistancePOS';
    printf "\t%s", 'MeanDistanceNEG';
    printf "\t%s", 'PresentPOS';
    printf "\t%s", 'PresentNEG';
    printf "\t%s", 'PresentDataset';
    printf "\t%s", 'MeanLog2AbundancePOS';
    printf "\t%s", 'MeanLog2AbundanceNEG';
    printf "\n";

    foreach $sample (@sample_array)
    {
        $sample_pos = 'NULL';
        $sample_neg = 'NULL';
        
        if ($single_file_mode ne 'neg')
        {
            $sample_pos    = $global_sample_map_hash{$sample}{pos};
        }
        if ($single_file_mode ne 'pos')
        {
            $sample_neg    = $global_sample_map_hash{$sample}{neg};
        }

        $sample_short  = $global_sample_map_hash{$sample}{short};

        $score         = $sample_merged_hash{$sample}{merged}{score};
        $log2scale_adj = $sample_merged_hash{$sample}{merged}{log2scale_adj};
        $present_adj   = $sample_merged_hash{$sample}{merged}{present_adj};
        $mean_dist_adj = $sample_merged_hash{$sample}{merged}{mean_dist_adj};
        $mean_log2_adj = $sample_merged_hash{$sample}{merged}{mean_log2_adj};

        $log2scale_pos = $sample_merged_hash{$sample}{pos}{log2scale};
        $mean_dist_pos = $sample_merged_hash{$sample}{pos}{mean_dist};
        $present_pos   = $sample_merged_hash{$sample}{pos}{present};
        $mean_log2_pos = $sample_merged_hash{$sample}{pos}{mean_log2};

        $log2scale_neg = $sample_merged_hash{$sample}{neg}{log2scale};
        $mean_dist_neg = $sample_merged_hash{$sample}{neg}{mean_dist};
        $present_neg   = $sample_merged_hash{$sample}{neg}{present};
        $mean_log2_neg = $sample_merged_hash{$sample}{neg}{mean_log2};

        # if ($sample_pos eq '')        { $sample_pos    = ''; }
        # if ($sample_neg eq '')        { $sample_neg    = ''; }
        if (!defined($log2scale_pos)) { $log2scale_pos = ''; }
        if (!defined($log2scale_neg)) { $log2scale_neg = ''; }
        if (!defined($mean_dist_pos)) { $mean_dist_pos = ''; }
        if (!defined($mean_dist_neg)) { $mean_dist_neg = ''; }
        if (!defined($present_pos))   { $present_pos   = ''; }
        if (!defined($present_neg))   { $present_neg   = ''; }
        if (!defined($mean_log2_pos)) { $mean_log2_pos = ''; }
        if (!defined($mean_log2_neg)) { $mean_log2_neg = ''; }

        $scale_pos = '';
        $scale_neg = '';
        if (is_number($log2scale_pos))
        {
            $scale_pos = pow(2, $log2scale_pos);
        }
        if (is_number($log2scale_neg))
        {
            $scale_neg = pow(2, $log2scale_neg);
        }

        printf "%s",   $sample;
        printf "\t%s", $sample_pos;
        printf "\t%s", $sample_neg;
        printf "\t%s", $sample_short;

        printf "\t%f", $score;
        printf "\t%f", pow(2, $log2scale_adj);
        printf "\t%f", $log2scale_adj;
        printf "\t%f", $mean_dist_adj;
        printf "\t%f", $present_adj;
        printf "\t%f", pow(2, $mean_log2_adj);
        printf "\t%f", $mean_log2_adj;

        printf "\t%s", $scale_pos;
        printf "\t%s", $scale_neg;
        printf "\t%s", $log2scale_pos;
        printf "\t%s", $log2scale_neg;
        printf "\t%s", $mean_dist_pos;
        printf "\t%s", $mean_dist_neg;
        printf "\t%s", $present_pos;
        printf "\t%s", $present_neg;
        printf "\t%s", $global_present_dataset;
        printf "\t%s", $mean_log2_pos;
        printf "\t%s", $mean_log2_neg;

        printf "\n";
    }
}


# begin main()

$filename_table  = shift;   # _sample_table.txt
$filename_sf_pos = shift;   # _pos_cleaned_scaling_factors.txt
$filename_sf_neg = shift;   # _neg_cleaned_scaling_factors.txt
$filename_fm_pos = shift;   # _pos_cleaned_findmedian.txt
$filename_fm_neg = shift;   # _neg_cleaned_findmedian.txt


if (!defined($filename_table) ||
    !defined($filename_sf_pos) ||
    !defined($filename_sf_neg) ||
    !defined($filename_fm_pos) ||
    !defined($filename_fm_neg))
{
    $program_name = basename($0);
    
    print "Usage: $program_name sample_table.txt pos_sf.txt neg_sf.txt pos_fm.txt neg_fm.txt\n";
    exit(1);
}


$single_file_mode  = '';
if (($filename_sf_pos eq 'NULL' && $filename_fm_pos eq 'NULL') ||
    ($filename_sf_neg eq 'NULL' && $filename_fm_neg eq 'NULL'))
{
    if ($filename_sf_pos ne 'NULL')
    {
        $single_file_mode = 'pos';
    }
    if ($filename_sf_neg ne 'NULL')
    {
        $single_file_mode = 'neg';
    }
    
    # ABORT -- they both can't be NULL
    if ($filename_sf_pos eq $filename_sf_neg)
    {
        printf STDERR "ABORT -- at least one input filename must be not \"NULL\"\n";
        exit(2);
    }
}


read_in_sample_table_file($filename_table);

if ($single_file_mode ne 'neg')
{
    read_in_scaling_factors_file($filename_sf_pos, 'pos');
    read_in_findmedian_file($filename_fm_pos, 'pos');
}

if ($single_file_mode ne 'pos')
{
    read_in_scaling_factors_file($filename_sf_neg, 'neg');
    read_in_findmedian_file($filename_fm_neg, 'neg');
}

score_qc();
score_qc();    # score again with updated median sample
print_qc();
