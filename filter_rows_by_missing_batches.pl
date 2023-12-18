#!/usr/bin/perl -w


# 2023-06-27:  update is_number() to not treat NaNs as numbers


use Scalar::Util qw(looks_like_number);
use File::Basename;


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


#$batch_file_name = '';
#$data_file_name  = '';
#$min_values      = shift;    # default of 2

$keep_all_cols_flag = 0;
$num_not_options    = 0;
$error_flag         = 0;
for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];
    
    if ($field =~ /^-/ && $field ne '-')
    {
        if ($field =~ /^--keep-all-cols/)
        {
            $keep_all_cols_flag = 1;
        }
        else
        {
            printf STDERR "ABORT -- unknown option %s\n", $field;
        
            $error_flag = 1;
        }
    }
    else
    {
        if ($num_not_options == 0)
        {
            $batch_file_name = $field;
        }
        elsif ($num_not_options == 1)
        {
            $data_file_name  = $field;
        }
        elsif ($num_not_options == 2)
        {
            $min_value = $field;
            
            if (!is_number($min_value))
            {
                printf STDERR "ABORT -- non-numeric minimum number of values per each\n";
            
                $error_flag = 1;
            }
        }

        $num_not_options++;
    }
}


if ($error_flag ||
    !defined($batch_file_name) ||
    !defined($data_file_name))
{
    $program_name = basename($0);

    printf STDERR "Usage: $program_name [options] combat_batch_input_file data_file [min # value per each]\n";
    printf STDERR "\n";
    printf STDERR "  Options:\n";
    printf STDERR "    --keep-all-cols    keep all columns, even if not in batch input file\n";
    printf STDERR "\n";
    printf STDERR "  First column is sample name\n";
    printf STDERR "  Batch column header is given as \"Batch\"\n";
    printf STDERR "  Any columns after Batch are treated as co-variates\n";
    printf STDERR "\n";
    printf STDERR "Output:\n";
    printf STDERR "  Rows with < 2 values per batch + co-variate combination are pruned\n";
    printf STDERR "  Only samples present in both files add to batch counts\n";
    printf STDERR "\n";
    printf STDERR "  By default, columns not present in the batch input file are removed\n";
    printf STDERR "  Use the --keep-all-cols flag to keep all columns\n";
    
    exit(1);
}

if (!defined($min_values))
{
    $min_values = 2;
}


open BATCH_FILE, "$batch_file_name" or die "can't open $batch_file_name\n";
open DATA_FILE,  "$data_file_name"  or die "can't open $data_file_name\n";


# read in the data header first, so we know which samples are present
$header_line_data = <DATA_FILE>;
$header_line_data =~ s/[\r\n]+//g;
@array = split /\t/, $header_line_data;
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
        
    $data_col_headers[$i] = $array[$i];
    
    # assume first column is row identifier, all others are sample names
    if ($i > 0)
    {
        $seen_sample_hash{$array[$i]} = 1;
    }
}
$header_line_data = join "\t", @array;


# batch metadata header line
$header_line_batch = <BATCH_FILE>;
$header_line_batch =~ s/[\r\n]+//g;
@array = split /\t/, $header_line_batch;
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    
    $header_col_hash{$array[$i]} = $i;
    $header_col_array[$i] = $array[$i];
}
$header_line_batch = join "\t", @array;


$batch_col = $header_col_hash{'Batch'};
if (!defined($batch_col))
{
    $batch_col = $header_col_hash{'batch'};
}
if (!defined($batch_col))
{
    $batch_col = $header_col_hash{'BATCH'};
}
if (!defined($batch_col))
{
    print STDERR "ABORT -- Batch column header not found in $batch_file_name\n";

    exit(1);
}


# read in batch assignments
# assume first column is sample name
while(defined($line=<BATCH_FILE>))
{
    $line =~ s/[\r\n]+//g;
    
    @array = split /\t/, $line;

    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
    }
    
    $sample = $array[0];
    
    # tack on additional covariates
    $batch = '';
    $ok_flag = 1;
    for ($i = $batch_col; $i < @header_col_array; $i++)
    {
        $header = $header_col_array[$i];
        $value  = $array[$i];
        
        if ($header =~ /\S/)
        {
            if ($value =~ /\S/)
            {
                if ($batch ne '')
                {
                    $batch .= '|||';
                }

                $batch .= sprintf "%s: %s", $header, $value;
            }
            else
            {
                $ok_flag = 0;
            }
        }
    }
    
    if ($ok_flag == 0)
    {
        print STDERR "ABORT -- $sample is not assigned a batch in batch file\n";
        exit(2);
    }
    
    #printf STDERR "%s\t%s\n", $sample, $batch;
    
    # only store batch info if we've seen the sample
    if (defined($seen_sample_hash{$sample}))
    {
        $sample_batch_hash{$sample} = $batch;
        $seen_batch_hash{$batch} = 1;
    }
}



@batch_array = sort keys %seen_batch_hash;
$num_batches = @batch_array;
$batch_string = join "\n", @batch_array;
printf STDERR "Found %s batches:\n%s\n", $num_batches, $batch_string;


%skip_col_hash = ();
if (@data_col_headers)
{
    printf "%s", $data_col_headers[0];
}
# by default, skip columns not present in batch file
for ($i = 1; $i < @data_col_headers; $i++)
{
    $sample = $data_col_headers[$i];
    $batch  = $sample_batch_hash{$sample};

    if (!defined($batch))
    {
        printf STDERR "WARNING -- skipping column not present in batch file: %s\n",
            $sample;

        if ($keep_all_cols_flag == 0)
        {
            $skip_col_hash{$i} = 1;

            next;
        }
    }
    
    printf "\t%s", $sample;
}
printf "\n";

# assume first column is row identifier, all other columns are samples
$num_total = 0;
$num_kept  = 0;
while(defined($line=<DATA_FILE>))
{
    $num_total++;

    $line =~ s/[\r\n]+//g;
    
    @array = split /\t/, $line;

    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
    }
    
    # clear batch counts
    foreach $batch (@batch_array)
    {
        $temp_batch_counts_hash{$batch} = 0;
    }
    
    # count observations per batch
    $happy_batch_count = 0;
    for ($i = 1; $i < @array; $i++)
    {
        $sample = $data_col_headers[$i];
        $batch  = $sample_batch_hash{$sample};
        
        if (!defined($batch))
        {
            next;
        }
        
        $value = $array[$i];
        
        if (is_number($value) && $value > 0)
        {
            $temp_batch_counts_hash{$batch} += 1;
        }
    }
    
    # count number of happy batches
    foreach $batch (@batch_array)
    {
        if ($temp_batch_counts_hash{$batch} >= $min_values)
        {
            $happy_batch_count++;
        }
    }

    $num_batches = @batch_array;

    # printf STDERR "FOOBAR\t%s\t%s\n", $happy_batch_count, $num_batches;
    
    # print row if all batches have enough samples
    if ($happy_batch_count == @batch_array)
    {
        if (@data_col_headers)
        {
            printf "%s", $array[0];
        }
        # by default, skip columns not present in batch file
        for ($i = 1; $i < @array; $i++)
        {
            if (defined($skip_col_hash{$i}))
            {
                next;
            }

            printf "\t%s", $array[$i];
        }
        printf "\n";
        
        $num_kept++;
    }
}

printf STDERR "Kept\t%d\tout of\t%d\trows\n", $num_kept, $num_total;
