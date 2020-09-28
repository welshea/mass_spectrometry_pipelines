#!/usr/bin/perl -w

# Run on merged file after the cleanup part of the pipeline has already been
# run.  Assume abundances are unlogged abundances.  Also assume that Peak
# Height, etc. has already been stripped from the sample names.  Assume we'll
# need to use the last of the columns added by the pipeline to detect where
# the samples start (the pipeline already has logic to place these columns
# right before the start of the samples).

# Assume row identifiers start with pos_ or neg_ so that we can tell which
# rows come from which ion mode.  The pipeline should already take care of
# this.

# Data will probably be El-MAVEN output, but it could still be MZmine or
# something else.  If it is El-MAVEN, the spike-in column will probably
# not be useful (all zeroes), but we'll check for it anyways, just in case
# we can use it to NOT normalize the spike-ins.

# Assume first column is row identifier


use Scalar::Util qw(looks_like_number);


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


$merged_file_name     = shift;
$scales_pos_file_name = shift;
$scales_neg_file_name = shift;

open MERGED,     $merged_file_name     or die "ABORT -- can't open file $merged_file_name\n";
open SCALES_POS, $scales_pos_file_name or die "ABORT -- can't open file $scales_pos_file_name\n";
open SCALES_NEG, $scales_neg_file_name or die "ABORT -- can't open file $scales_neg_file_name\n";


# assume column order is SampleID, Scale, Log2Scale
$line = <SCALES_POS>;
while (defined($line=<SCALES_POS>))
{
    $line =~ s/[\r\n]+//g;
    @array = split /\t/, $line;

    $sampleid  = $array[0];
    $scale     = $array[1];
    # $log2scale = $array[2];

    if (defined($sampleid)  && $sampleid  =~ /[0-9]/ &&
        defined($scale)     && $scale     =~ /[0-9]/)
    {
        # strip sampleid of pos/neg
        # apply same clean up used in merge script
        if ($sampleid =~ s/(^|[^A-Za-z0-9]+)pos([^A-Za-z0-9]+|$)/$2/i ||
            $sampleid =~ s/pos$//i)
        {
            # clean up underscores, etc.
            $sampleid =~ s/[_ ]+/_/g;
            $sampleid =~ s/\-+/\-/g;
            $sampleid =~ s/^[_ -]//;
            $sampleid =~ s/[_ -]$//;
        }
    
        $scale_hash{$sampleid}{pos} = $scale;
    }
}
close SCALES_POS;

$line = <SCALES_NEG>;
while (defined($line=<SCALES_NEG>))
{
    $line =~ s/[\r\n]+//g;
    @array = split /\t/, $line;

    $sampleid  = $array[0];
    $scale     = $array[1];
    # $log2scale = $array[2];

    if (defined($sampleid)  && $sampleid  =~ /[0-9]/ &&
        defined($scale)     && $scale     =~ /[0-9]/)
    {
        # strip sampleid of pos/neg
        # apply same clean up used in merge script
        if ($sampleid =~ s/(^|[^A-Za-z0-9]+)neg([^A-Za-z0-9]+|$)/$2/i ||
            $sampleid =~ s/neg$//i)
        {
            # clean up underscores, etc.
            $sampleid =~ s/[_ ]+/_/g;
            $sampleid =~ s/\-+/\-/g;
            $sampleid =~ s/^[_ -]//;
            $sampleid =~ s/[_ -]$//;
        }
    
        $scale_hash{$sampleid}{neg} = $scale;
    }
}
close SCALES_NEG;




$line = <MERGED>;
$line =~ s/[\r\n]+//g;
@header_col_array_merged = split /\t/, $line;
for ($col = 0; $col < @header_col_array_merged; $col++)
{
    $field = $header_col_array_merged[$col];
    
    $header_col_hash_merged{$field} = $col;
}

# use last of the pipeline generated columns to determine start of sample data
$first_abundance_col = $header_col_hash_merged{'Non-Spikein Identified Flag'};
if (defined($first_abundance_col))
{
    $first_abundance_col += 1;
}
else
{
    printf STDERR "ABORT -- cannot find Non-Spikein Identified Flag column in file %s\n",
        $merged_file_name;
    
    exit(1);
}

#printf STDERR "First abundance col: %s\n", $first_abundance_col;

print "$line\n";

# normalize the data
$log2_const = log(2.0);
while (defined($line=<MERGED>))
{
    $line =~ s/[\r\n]+//g;
    @array = split /\t/, $line;
    
    $rowid = $array[0];

    if (!defined($rowid) || !($rowid =~ /^(pos|neg)_/))
    {
        print STDERR "WARNING -- bad rowid:\t$rowid\n";
    
        print "$line\n";
    
        next;
    }
    
    $pos_neg = 'pos';
    if ($rowid =~ /^neg_/)
    {
        $pos_neg = 'neg';
    }
    
    for ($col = $first_abundance_col; $col < @array; $col++)
    {
        $sampleid = $header_col_array_merged[$col];

        undef($scale);
        if (defined($scale_hash{$sampleid}))
        {
            $scale = $scale_hash{$sampleid}{$pos_neg};
        }
        
        if (!defined($scale))
        {
            print STDERR "WARNING -- undefined scaling factor:\t$rowid\t$sampleid\n";
        }
        
        $value = $array[$col];
        
        if (!defined($value))
        {
            $value = '';
        }
        
        if (is_number($value))
        {
            if ($value > 0)
            {
                # overwrite with normalized data
                $array[$col] = log($scale * $value) / $log2_const;
            }
            else
            {
                $array[$col] = '';
            }
        }
    }
    
    $line_new = join "\t", @array;
    
    print "$line_new\n";
}
close MERGED;
