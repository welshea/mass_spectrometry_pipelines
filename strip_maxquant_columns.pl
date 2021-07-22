#!/usr/bin/perl -w

# 2021-07-22:  add re-capitalization code from reformat_modification_sites.pl
# 2020-06-15:  fix broken TMT handling introduced in heavy/light changes
# 2020-04-07:  remove even more SILAC-related summary columns and ratios
# 2020-03-12:  skip "Intensity L" and "Intensity H" summary columns
#              warning, will not correctly handle case where there are real
#               samples named H and L together with non-H/L sample names

$filename = shift;

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


    # damn it, Maxquant keeps changing case between versions...
    # this is *ROYALLY*screwing everything up
    #
    # I am just going to conform all the capitalization here
    # If a word is >= 3 long, and the first letter of a word isn't
    #  already capitalized, capitalize it
    #
    # 2021-07-22 Maxquant broke capitalization yet again, this time with
    #  Desthiobiotin vs. DesthioBiotin
    # HACK: conform DesthioBiotin to Desthiobiotin,
    #  since DesthioBiotin is poorly capitalized anyways
    #  they've spelled in wrong in the past too, so fix that typo as well
    #
    $array[$i] =~ s/DesthioBiotin/Desthiobiotin/ig;
    $array[$i] =~ s/Dethiobiotin/Desthiobiotin/ig;
    @word_array = split /\s+/, $array[$i];
    for ($j = 0; $j < @word_array; $j++)
    {
        $c = substr $word_array[$j], 0, 1;
        if (length $word_array[$j] >= 3 && $c =~ /[a-z]/)
        {
            $c = uc($c);
            $word_array[$j] =~ s/^([a-z])/$c/;
        }
    }
    $array[$i] = join ' ', @word_array;

    
    $field = $array[$i];
    $header_col_hash{$field} = $i;
    $header_col_array[$i] = $field;
}


# count Intensity and iBaq columns
# single sample only has one of either/each
# multiple samples should exclude them, since they are the sum of all samples
$count_intensity        = 0;
$count_intensity_non_hl = 0;
$count_ibaq             = 0;
$count_ibaq_non_hl      = 0;
$count_reporter         = 0;
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];
    
    if ($field =~ /^Intensity\b/i)
    {
        $count_intensity++;
        
        if (!($field =~ /^Intensity$/i) &&
            !($field =~ /^Intensity H$/i) &&
            !($field =~ /^Intensity L$/i))
        {
            $count_intensity_non_hl++;
        }
    }
    if ($field =~ /^iBAQ\b/i)
    {
        $count_ibaq++;

        if (!($field =~ /^iBAQ$/i) &&
            !($field =~ /^iBAQ H$/i) &&
            !($field =~ /^iBAQ L$/i))
        {
            $count_ibaq_non_hl++;
        }
    }
    if ($field =~ /^Reporter intensity \d/i)
    {
        $count_reporter++;
    }
}


# flag columns to remove, as they clutter up the file
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];
    
    if ($field =~ /^[A-Z] Count$/i)
    {
        $col_to_remove_hash{$col} = 1;
    }

    if ($field =~ /^Identification type /i)
    {
        $col_to_remove_hash{$col} = 1;
    }

    if ($field =~ /^Fraction /i)
    {
        $col_to_remove_hash{$col} = 1;
    }

    if ($field =~ /^Slice /i)
    {
        $col_to_remove_hash{$col} = 1;
    }

    if ($field =~ /^Experiment /i)
    {
        $col_to_remove_hash{$col} = 1;
    }

    if ($count_intensity_non_hl > 0)
    {
        if ($field =~ /^Intensity$/i ||
            $field =~ /^Intensity H$/i ||
            $field =~ /^Intensity L$/i)
        {
            $col_to_remove_hash{$col} = 1;
        }
    }
    if ($count_ibaq_non_hl > 0)
    {
        if ($field =~ /^iBAQ$/i ||
            $field =~ /^iBAQ H$/i ||
            $field =~ /^iBAQ L$/i)
        {
            $col_to_remove_hash{$col} = 1;
        }
    }
    # TMT experiment, all Intensity columns are summary columns
    if ($count_reporter > 0)
    {
        if ($field =~ /^Intensity/i)
        {
            $col_to_remove_hash{$col} = 1;
        }
    }

    # protein-level stuff

    if ($field =~ /^Peptides \w+/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^Razor \+ unique peptides \w+/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^Unique peptides \w+/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^Sequence coverage \w+/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^LFQ intensity \w+/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^MS\/MS Count \w+/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^Peptide IDs$/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^Peptide is razor$/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^Mod. peptide IDs$/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^Evidence IDs$/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^MS\/MS IDs$/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^Best MS\/MS$/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^Oxidation \(M\) site IDs$/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^Oxidation \(M\) site positions$/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^Reporter intensity count \d+/)
    {
        $col_to_remove_hash{$col} = 1;
    }

    if ($field =~ /^Occupancy/i)
    {
        $col_to_remove_hash{$col} = 1;
    }

    # they're all NaN's anyways...
    if ($field =~ /^Ratio mod\/base/i)
    {
        $col_to_remove_hash{$col} = 1;
    }

    # H/L ratio stuff, ratios will need to be recalculated after renorm
    if ($field =~ /^Ratio H\/L/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    
    # more sample-specific stuff we don't want
    if ($field =~ /^Localization Prob \w+/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^Score Diff \w+/i)
    {
        $col_to_remove_hash{$col} = 1;
    }
    if ($field =~ /^PEP \w+/i)
    {
        $col_to_remove_hash{$col} = 1;
    }

    if ($field =~ /^Score \w+/i)
    {
        # Score Diff sample
        if ($field =~ /^Score Diff/i)
        {
            if ($field =~ /^Score Diff \w+/i)
            {
                $col_to_remove_hash{$col} = 1;
            }
        }
        # Score sample, where sample != "Diff"
        else
        {
            $col_to_remove_hash{$col} = 1;
        }
    }

    
    # blank column headers
    if ($field eq '')
    {
        $col_to_remove_hash{$col} = 1;
    }
}


# deal with even more H/L summary column crap
$silac_flag             = 0;
$count_remaining_hl     = 0;
$count_remaining_non_hl = 0;
for ($col = 0; $col < @header_col_array; $col++)
{
    # skip columns we are already removing
    if (defined($col_to_remove_hash{$col}))
    {
        next;
    }

    $field = $header_col_array[$col];

    if ($field =~ /^Intensity/)
    {
        if ($field =~ /^Intensity H /i ||
            $field =~ /^Intensity L /i)
        {
            $count_remaining_hl++;
        }
        else
        {
            $count_remaining_non_hl++;
        }
    }
}

if ($count_remaining_hl &&
    $count_remaining_hl >= 2 * $count_remaining_non_hl)
{
    $silac_flag = 1;
}


# remove even more SILAC summary columns
if ($silac_flag)
{
  for ($col = 0; $col < @header_col_array; $col++)
  {
    # skip columns we are already removing
    if (defined($col_to_remove_hash{$col}))
    {
        next;
    }

    $field = $header_col_array[$col];

    if ($field =~ /^Intensity/)
    {
        if (!($field =~ /^Intensity H /i) &&
            !($field =~ /^Intensity L /i))
        {
            $col_to_remove_hash{$col} = 1;
        }
    }
  }
}


# deal with ___# issues
# newer versions of Maxquant fail to report the summarized columns
%cols_to_sum_hash = ();
%summarized_field_hash = ();
@summarized_field_array = ();
$num_summarized_fields = 0;
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];

    # skip columns we are removing
    if (defined($col_to_remove_hash{$col}))
    {
        next;
    }
    
    # skip types of columns we don't want
    if ($field =~ / intensity count /i)
    {
        next;
    }


    if ($field =~ /___[123]$/)
    {
        $root = $field;
        $root =~ s/\s*___[123]$//;

        if (!defined($header_col_hash{$root}))
        {
            $cols_to_sum_hash{$col} = 1;
            
            # store new headers in original input order
            if (!defined($summarized_field_hash{$root}))
            {
                $summarized_field_array[$num_summarized_fields++] = $root;
            }

            $summarized_field_hash{$root} = 1;
        }

        $col_to_remove_hash{$col} = 1;
    }
}

# print header line
$print_flag = 0;
for ($col = 0; $col < @header_col_array; $col++)
{
    # deal with ___# issues
    # newer versions of Maxquant fail to report the summarized columns
    if (defined($col_to_remove_hash{$col}))
    {
        next;
    }
    
    if ($print_flag == 1)
    {
        print "\t";
    }
    else
    {
        $print_flag = 1;
    }

    print $header_col_array[$col];
}

# deal with ___# issues
# newer versions of Maxquant fail to report the summarized columns
@cols_to_sum_array = sort {$a<=>$b} keys %cols_to_sum_hash;
foreach $field (@summarized_field_array)
{
    if ($print_flag == 1)
    {
        print "\t";
    }
    else
    {
        $print_flag = 1;
    }

    print $field;
}

print "\n";


while(defined($line=<INFILE>))
{
    $line =~ s/[\r\n]+//g;

    @array = split /\t/, $line;

    $print_flag = 0;

    for ($col = 0; $col < @array; $col++)
    {
        if (defined($col_to_remove_hash{$col}))
        {
            next;
        }

        $field = $array[$col];
        $field =~ s/^\s+//;
        $field =~ s/\s+$//;
        $field =~ s/\s+/ /g;
        
        # clean accession junk
        $header = $header_col_array[$col];
        if ($header =~ /^Proteins$/i ||
            $header =~ /^Leading .*?protein$/i ||
            $header =~ /^Protein IDs$/i ||
            $header =~ /^Majority protein IDs$/i)
        {
#print STDERR "FOOBAR   $field\n";
#            $field =~ s/sp\|[^\|\;\,]+\|([^\|\;\,]+)/$1/g;
#            $field =~ s/ENSEMBL:([A-Za-z0-9-\.]+)/$1/ig;
#            $field =~ s/REFSEQ:([A-Za-z0-9\.]+)/$1/ig;
        }
        
    
        if ($print_flag == 1)
        {
            print "\t";
        }
        else
        {
            $print_flag = 1;
        }
        
        print $field;
    }


    # deal with ___# issues
    # newer versions of Maxquant fail to report the summarized columns
    %sums_hash = ();
    foreach $col (@cols_to_sum_array)
    {
        $value = $array[$col];
        
        if (defined($value) && $value =~ /[0-9]/)
        {
            $header = $header_col_array[$col];
            $root   = $header;
            $root   =~ s/\s*___[123]$//;

            if (!defined($sums_hash{$root}))
            {
                $sums_hash{$root} = 0;
            }
            
            $sums_hash{$root} += $value;
        }
    }
    
    foreach $header (@summarized_field_array)
    {
        if ($print_flag == 1)
        {
            print "\t";
        }
        else
        {
            $print_flag = 1;
        }
        
        $value = $sums_hash{$header};
        
        if (defined($value) && $value =~ /[0-9]/)
        {
            print $value;
        }
    }

    print "\n";
}
