#!/usr/bin/perl -w

#take the iBAQ data and move it into the Intensity columns

$filename = shift;

open INFILE, "$filename" or die "can't open $filename\n";

# header line
$line = <INFILE>;
@array = split /\t/, $line;
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    $array[$i] =~ s/\s+/ /g;
    
    $field = $array[$i];
#   $header_col_hash{$field} = $i;
    $header_col_array[$i] = $field;
}


# flag columns to remove, as they clutter up the file
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $array[$col];

    if ($field =~ /^Intensity\s(.*)$/i)
    {
        $sample = $1;
        
        if (defined($sample_col_hash{$sample}{intensity}))
        {
            printf "ABORT -- multiple sample columns for %s\n", $field;
            exit(1);
        }
        
        $sample_col_hash{$sample}{intensity} = $col;
        $sample_col_ok_status_hash{$col} = 0;	# default to not ok
    }

    if ($field =~ /^iBAQ\s(.*)$/i)
    {
        $sample = $1;

        if (defined($sample_col_hash{$sample}{ibaq}))
        {
            printf "ABORT -- multiple sample columns for %s\n", $field;
            exit(1);
        }

        $sample_col_hash{$sample}{ibaq} = $col;
        $sample_col_ok_status_hash{$col} = 0;	# default to not ok
    }
}


# identify samples that have both iBAQ and Intensity data
foreach $sample (sort keys %sample_col_hash)
{
    $intensity_col = $sample_col_hash{$sample}{intensity};
    $ibaq_col      = $sample_col_hash{$sample}{ibaq};

    if (defined($intensity_col) && defined($ibaq_col))
    {
        $sample_col_ok_status_hash{$intensity_col} = 1;
        
        $intensity_col_to_ibaq_col_hash{$intensity_col} = $ibaq_col;
    }
}



# remove all iBAQ columns, replace Intensity data with iBAQ data

# print header line
$print_flag = 0;
for ($col = 0; $col < @header_col_array; $col++)
{
    # skip sample columns we don't want to print
    if (defined($sample_col_ok_status_hash{$col}) &&
        $sample_col_ok_status_hash{$col} == 0)
    {
        next;
    }
    
    if ($print_flag == 0)
    {
        print "$header_col_array[$col]";
    }
    else
    {
        print "\t$header_col_array[$col]";
    }
    
    $print_flag = 1;
}
print "\n";


while(defined($line=<INFILE>))
{
    @array = split /\t/, $line;
    for ($col = 0; $col < @array; $col++)
    {
        $array[$col] =~ s/^\s+//;
        $array[$col] =~ s/\s+$//;
    }

    $print_flag = 0;
    for ($col = 0; $col < @array; $col++)
    {
        # skip sample columns we don't want to print
        if (defined($sample_col_ok_status_hash{$col}) &&
            $sample_col_ok_status_hash{$col} == 0)
        {
            next;
        }
        
        $field = $array[$col];
        

        # replace Intensity data with iBAQ data in place
        $ibaq_col = $intensity_col_to_ibaq_col_hash{$col};
        if (defined($ibaq_col))
        {
            $field = $array[$ibaq_col];
        }


        if ($print_flag == 0)
        {
            print "$field";
        }
        else
        {
            print "\t$field";
        }
    
        $print_flag = 1;
    }
    print "\n";
}
