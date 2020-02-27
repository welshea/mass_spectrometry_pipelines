#!/usr/bin/perl -w

sub is_number
{
    return $_[0] =~ /^([+-]?)(?=[0-9]|\.[0-9])[0-9]*(\.[0-9]*)?([Ee]([+-]?[0-9]+))?$/;
}

# choose longest field, then highest in regular sort order
sub compare_fields
{
    my $len1 = length $a;
    my $len2 = length $b;
    my $ok_flag1 = 0;
    my $ok_flag2 = 0;
    
    # check for alphanumeric
    if ($a =~ /[A-Za-z0-9]/)
    {
        $ok_flag1 = 1;
    }
    if ($b =~ /[A-Za-z0-9]/)
    {
        $ok_flag2 = 1;
    }
    
    # also check for + symbol
    if ($a =~ /\+/)
    {
        $ok_flag1 = 1;
    }
    if ($b =~ /\+/)
    {
        $ok_flag2 = 1;
    }
    
    
    # put blank/--- fields last
    if ($ok_flag1 && !$ok_flag2) { return -1; }
    if ($ok_flag2 && !$ok_flag1) { return  1; }
    
    # longest first, regardless of if they are numbers
    if ($len1 > $len2) { return -1; }
    if ($len1 < $len2) { return  1; }
    
    # highest numerically
    if (is_number($a) && is_number($b))
    {
        # keep smallest
        if ($global_sort_header =~ /^Localization Prob/i ||
            $global_sort_header =~ /^Score/i ||
            $global_sort_header =~ /Probability/i)
        {
            if ($a < $b) { return -1; }
            if ($a > $b) { return  1; }
        }
        # keep biggest
        else
        {
            # use absolute values, in case of Mass Error and whatnot
            if (abs($a) > abs($b)) { return -1; }
            if (abs($a) < abs($b)) { return  1; }

            # take positive in case of equal positive and negative
            if ($a > $b) { return -1; }
            if ($a < $b) { return  1; }
        }
        
        return 0;
    }
    
    # highest alphabetically
    if ($a gt $b) { return -1; }
    if ($a lt $b) { return  1; }
    
    return 0;
}



$num_files = 0;

# read in command line arguments
for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        {
            printf "ABORT -- unknown option %s\n", $field;
            $syntax_error_flag = 1;
        }
    }
    else
    {
        $filename = $field;
        $filename_array[$num_files] = $filename;
#        $filename_to_index_hash{$filename} = $num_files;

        $num_files++;
    }
}


# print syntax error message
if ($num_files == 0 || $syntax_error_flag)
{
    printf "Syntax: merge_proteomics.pl [options] tab_delimited_text_files*.txt\n";
    exit(1);
}


# find common filename prefix
$common_prefix = $filename_array[0];
for ($i = 0; $i < $num_files; $i++)
{
    $filename = $filename_array[$i];
    
    # XOR the two strings, using "" to force string XOR
    $xor = "$common_prefix" ^ "$filename";
    
    # common prefix will be null'd out
    $xor =~ /^\0*/;
    $prefix_length = $+[0];    # length of match
    
    $common_prefix = substr $common_prefix, 0, $prefix_length;
    
    # only go up to the first _ or -
#    $common_prefix =~ s/(_).*/$1/;
#    $common_prefix =~ s/(-).*/$1/;
}


# remove common prefix, store in mapping table
for ($i = 0; $i < $num_files; $i++)
{
    $filename = $filename_array[$i];
    
    $filename_stripped = $filename;
    $filename_stripped =~ s/^$common_prefix//;
#    $filename_stripped =~ s/\.txt$//;
#    $filename_stripped =~ s/_blessed//;
#    $filename_stripped =~ s/_reannotated//;
#    $filename_stripped =~ s/_annotated//;
    
    $filename_stripped_hash{$filename} = $filename_stripped;
}


# find common filename suffix after removing common prefix
$filename = $filename_array[0];
$common_suffix = $filename_stripped_hash{$filename};
for ($i = 0; $i < $num_files; $i++)
{
    $filename = $filename_array[$i];
    $filename = $filename_stripped_hash{$filename};
    
    # XOR the two strings, using "" to force string XOR
    $xor = "$common_suffix" ^ "$filename";

    # common suffix will be null'd out
    $xor =~ /(\0*)$/;
    $xor = $1;
    $suffix_length = length $xor;

    $offset        = length($common_suffix) - $suffix_length;
    $common_suffix = substr $common_suffix, $offset, $suffix_length;
}


# remove common suffix, store in mapping table
for ($i = 0; $i < $num_files; $i++)
{
    $filename = $filename_array[$i];
    
    $filename_stripped =  $filename_stripped_hash{$filename};
    $filename_stripped =~ s/$common_suffix$//;

    $filename_stripped_hash{$filename} = $filename_stripped;
}


# read in files
for ($f = 0; $f < $num_files; $f++)
{
    $filename = $filename_array[$f];
    $num_lines = 0;

    open INFILE, "$filename" or die "can't open $filename\n";
    while(defined($line=<INFILE>))
    {
        $line =~ s/[\r\n]+//g;
        $line =~ s/\"//g;
        
        # skip comment lines and blank lines
        if ($line =~ /^#/ || !($line =~ /\S/))
        {
            next;
        }
    
        @array = split /\t/, $line;
    
        for ($i = 0; $i < @array; $i++)
        {
            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;

            # fix Dethiobiotin typos
            $array[$i] =~ s/Dethiobiotin/Desthiobiotin/g;
            
            $value = $array[$i];
            
            
            # header line
            if ($num_lines == 0)
            {
                # headers to skip
                if ($value =~ /^Best Identification /i ||
                    $value =~ /^Best Localization /i ||
                    $value =~ /^Diagnostic peak/i ||
                    $value =~ /^id$/i ||
                    $value =~ /^Mod\. Peptide IDs/i ||
                    $value =~ /^Peptide IDs/i ||
                    $value =~ /^Protein Group IDs/i)
                {
                    next;
                }

                if ($value =~ /^Ratio /i || $value =~ /^Intensity$/i ||
                    $value =~ /^Intensity__/i || $value =~ /^Occupancy/i)
                {
                    next;
                }
            
                # keep first occurrence
                if (!defined($header_to_col_hash{$filename}{$value}))
                {
                    $header_to_col_hash{$filename}{$value} = $i;
                }

                $col_to_header_hash{$filename}{$i} = $value;
                $header_to_multiple_col_hash{$filename}{$value}{$i} = 1;
                
                # flag this file as being Heavy/Light
                # must use lookahead/lookbehind due to \b including _ 
                if ($value =~ /(?<![A-Za-z0-9])H\/L(?![A-Za-z0-9])/)
                {
                    $is_heavy_light_hash{$filename} = 1;
                }
            }
            else
            {
                $data_hash{$filename}{data}[$num_lines][$i] = $value;
            }
        }

        $num_lines++;
    }
    $data_hash{$filename}{num_lines} = $num_lines;
    close INFILE;
}


# find ModificationID columns and map them to rows
for ($f = 0; $f < @filename_array; $f++)
{
    $filename = $filename_array[$f];
    
    @col_array = sort {$a<=>$b} keys
                 %{$header_to_multiple_col_hash{$filename}{'ModificationID'}};
    if (!defined(@col_array) || @col_array == 0)
    {
        @col_array = sort {$a<=>$b} keys
                 %{$header_to_multiple_col_hash{$filename}{'RowIdentifier'}};
    }

    if (!defined(@col_array) || @col_array == 0)
    {
        die "ABORT -- ModificationID column not found: $filename\n";
    }
    
    $num_lines = $data_hash{$filename}{num_lines};

    for ($r = 1; $r < $num_lines; $r++)
    {
        %temp_hash = ();
        foreach $col (@col_array)
        {
            $value = $data_hash{$filename}{data}[$r][$col];
            
            if ($value =~ /[A-Za-z0-9]/)
            {
                $temp_hash{$value} = 1;
            }
        }
        
        foreach $value (keys %temp_hash)
        {
            $modificationid_row_hash{$filename}{$value} = $r;
            $modificationid_seen_hash{$value} = 1;
        }
    }
}


# find common headers and assign them a common ordering
for ($f = 0; $f < @filename_array; $f++)
{
    $filename = $filename_array[$f];

    foreach $header (sort keys %{$header_to_col_hash{$filename}})
    {
        # skip H/L, H, L headers in heavy/light experiments
        # since they do not contain sample IDs, and are thus falsely in
        # common
        #
        # must use lookahead/lookbehind due to \b including _ 
        #
        # we only want to keep Ratio H/L Normalized separate within Ratio*
        #
        if ($header =~ /(?<![A-Za-z0-9])H\/L[ _]+?Normalized$/)
        {
            printf STDERR "%s\n", $header;
            next;
        }
        if ($header =~ /(?<![A-Za-z0-9])(H|L)$/)
        {
            printf STDERR "%s\n", $header;
            next;
        }

        # make darn sure we don't summarize measurements
        if ($header =~ /^Ratio /i || $header =~ /^Intensity /i ||
            $header =~ /^Intensity$/i || $header =~ /^Intensity__/i ||
            $header =~ /^IRON / ||
            $header =~ /^Occupancy/i ||
            $header =~ /TMT-/)
        {
            next;
        }

        if (!defined($header_info_hash{$header}))
        {
            $header_info_hash{$header}{count} = 0;
            $header_info_hash{$header}{rank}  = 0;
        }

        $header_info_hash{$header}{count} += 1;
        $header_info_hash{$header}{rank}  +=
            $header_to_col_hash{$filename}{$header} + 1;
    }
}

$max_count = 0;
foreach $header (keys %header_info_hash)
{
    $count = $header_info_hash{$header}{count};

    if ($count > $max_count)
    {
        $max_count = $count;
    }
}

foreach $header (keys %header_info_hash)
{
    $count = $header_info_hash{$header}{count};
    
    if ($count == $max_count)
    {
        $rank = $header_info_hash{$header}{rank} / $max_count;
        $common_header_hash{$header} = $rank;
    }
}

@common_header_array = sort {$common_header_hash{$a} <=>
                             $common_header_hash{$b}}
                            keys %common_header_hash;

$final_num_cols = 0;


# order merged common headers first
foreach $header (@common_header_array)
{
    for ($f = 0; $f < @filename_array; $f++)
    {
        $filename = $filename_array[$f];
        $filename_stripped = $filename_stripped_hash{$filename};
        $header_new = sprintf "%s_%s", $header, $filename_stripped;
        
        @col_array = sort {$a<=>$b} keys
                      %{$header_to_multiple_col_hash{$filename}{$header}};
        
        foreach $col (@col_array)
        {
            $final_order_hash{$final_num_cols}{filename}    = $filename;
            $final_order_hash{$final_num_cols}{header}      = $header_new;
            $final_order_hash{$final_num_cols}{header_orig} = $header;
            $final_order_hash{$final_num_cols}{col}         = $col;

            $final_num_cols++;
        }
    }
}


# not in common, not anything else that I particularly care about
for ($f = 0; $f < @filename_array; $f++)
{
    $filename = $filename_array[$f];
    
    @col_array_orig = sort {$a<=>$b} keys %{$col_to_header_hash{$filename}};
    
    foreach $i (@col_array_orig)
    {
        $header = $col_to_header_hash{$filename}{$i};
        $filename_stripped = $filename_stripped_hash{$filename};
        $header_new = sprintf "%s_%s", $header, $filename_stripped;
        
        if (defined($common_header_hash{$header}))
        {
            next;
        }
        
        if ($header =~ /^Localization Prob/i || $header =~ /^Score/i ||
            $header =~ /^PEP/i)
        {
            next;
        }
        
        # we need to skip ALL the individual intensities, so they aren't
        #  duplicated or summarized, etc.
        if ($header =~ /^Ratio /i || $header =~ /^Intensity /i ||
            $header =~ /^Intensity$/i || $header =~ /^Intensity__/i ||
            $header =~ /^IRON / ||
            $header =~ /^Occupancy/i ||
            $header =~ /TMT-/)
        {
            next;
        }

        
        $final_order_hash{$final_num_cols}{filename}    = $filename;
        $final_order_hash{$final_num_cols}{header}      = $header_new;
        $final_order_hash{$final_num_cols}{header_orig} = $header;
        $final_order_hash{$final_num_cols}{col}         = $i;

        $final_num_cols++;
    }
}


# not in common, Localization Probe
for ($f = 0; $f < @filename_array; $f++)
{
    $filename = $filename_array[$f];
    
    @col_array_orig = sort {$a<=>$b} keys %{$col_to_header_hash{$filename}};
    
    foreach $i (@col_array_orig)
    {
        $header = $col_to_header_hash{$filename}{$i};

        $header_new = $header;
        if (1 || defined($is_heavy_light_hash{$filename}))
        {
            $filename_stripped = $filename_stripped_hash{$filename};
            $header_new = sprintf "%s_%s", $header, $filename_stripped;
        }

        if (defined($common_header_hash{$header}))
        {
            next;
        }
        
        if (!($header =~ /^Localization Prob/))
        {
            next;
        }
        
        $final_order_hash{$final_num_cols}{filename}    = $filename;
        $final_order_hash{$final_num_cols}{header}      = $header_new;
        $final_order_hash{$final_num_cols}{header_orig} = $header;
        $final_order_hash{$final_num_cols}{col}         = $i;

        $final_num_cols++;
    }
}


# not in common, Score Diff
for ($f = 0; $f < @filename_array; $f++)
{
    $filename = $filename_array[$f];
    
    @col_array_orig = sort {$a<=>$b} keys %{$col_to_header_hash{$filename}};
    
    foreach $i (@col_array_orig)
    {
        $header = $col_to_header_hash{$filename}{$i};
        $header_new = $header;
        if (1 || defined($is_heavy_light_hash{$filename}))
        {
            $filename_stripped = $filename_stripped_hash{$filename};
            $header_new = sprintf "%s_%s", $header, $filename_stripped;
        }

        if (defined($common_header_hash{$header}))
        {
            next;
        }
        
        if (!($header =~ /^Score Diff/))
        {
            next;
        }
        
        $final_order_hash{$final_num_cols}{filename}    = $filename;
        $final_order_hash{$final_num_cols}{header}      = $header_new;
        $final_order_hash{$final_num_cols}{header_orig} = $header;
        $final_order_hash{$final_num_cols}{col}         = $i;
        
        $final_num_cols++;
    }
}


# not in common, PEP
for ($f = 0; $f < @filename_array; $f++)
{
    $filename = $filename_array[$f];
    
    @col_array_orig = sort {$a<=>$b} keys %{$col_to_header_hash{$filename}};
    
    foreach $i (@col_array_orig)
    {
        $header = $col_to_header_hash{$filename}{$i};
        $header_new = $header;
        if (1 || defined($is_heavy_light_hash{$filename}))
        {
            $filename_stripped = $filename_stripped_hash{$filename};
            $header_new = sprintf "%s_%s", $header, $filename_stripped;
        }

        if (defined($common_header_hash{$header}))
        {
            next;
        }
        
        if (!($header =~ /^PEP /))
        {
            next;
        }
        
        $final_order_hash{$final_num_cols}{filename}    = $filename;
        $final_order_hash{$final_num_cols}{header}      = $header_new;
        $final_order_hash{$final_num_cols}{header_orig} = $header;
        $final_order_hash{$final_num_cols}{col}         = $i;
        
        $final_num_cols++;
    }
}


# not in common, Score
for ($f = 0; $f < @filename_array; $f++)
{
    $filename = $filename_array[$f];
    
    @col_array_orig = sort {$a<=>$b} keys %{$col_to_header_hash{$filename}};
    
    foreach $i (@col_array_orig)
    {
        $header = $col_to_header_hash{$filename}{$i};
        $header_new = $header;
        if (1 || defined($is_heavy_light_hash{$filename}))
        {
            $filename_stripped = $filename_stripped_hash{$filename};
            $header_new = sprintf "%s_%s", $header, $filename_stripped;
        }

        if (defined($common_header_hash{$header}))
        {
            next;
        }
        
        if (!($header =~ /^Score /))
        {
            next;
        }

        if ($header =~ /^Score Diff/)
        {
            next;
        }
        
        $final_order_hash{$final_num_cols}{filename}    = $filename;
        $final_order_hash{$final_num_cols}{header}      = $header_new;
        $final_order_hash{$final_num_cols}{header_orig} = $header;
        $final_order_hash{$final_num_cols}{col}         = $i;
        
        $final_num_cols++;
    }
}


# not in common, Occupancy
for ($f = 0; $f < @filename_array; $f++)
{
    $filename = $filename_array[$f];
    
    @col_array_orig = sort {$a<=>$b} keys %{$col_to_header_hash{$filename}};
    
    foreach $i (@col_array_orig)
    {
        $header = $col_to_header_hash{$filename}{$i};
        $header_new = $header;
        if (1 || defined($is_heavy_light_hash{$filename}))
        {
            $filename_stripped = $filename_stripped_hash{$filename};
            $header_new = sprintf "%s_%s", $header, $filename_stripped;
        }

        if (defined($common_header_hash{$header}))
        {
            next;
        }
        
        if (!($header =~ /^Occupancy /))
        {
            next;
        }
        
        $final_order_hash{$final_num_cols}{filename}    = $filename;
        $final_order_hash{$final_num_cols}{header}      = $header_new;
        $final_order_hash{$final_num_cols}{header_orig} = $header;
        $final_order_hash{$final_num_cols}{col}         = $i;
        
        $final_num_cols++;
    }
}


# not in common, Ratio (not Normalized)
for ($f = 0; $f < @filename_array; $f++)
{
    $filename = $filename_array[$f];
    
    @col_array_orig = sort {$a<=>$b} keys %{$col_to_header_hash{$filename}};
    
    foreach $i (@col_array_orig)
    {
        $header = $col_to_header_hash{$filename}{$i};
        $header_new = $header;
        if (1 || defined($is_heavy_light_hash{$filename}))
        {
            $filename_stripped = $filename_stripped_hash{$filename};
            $header_new = sprintf "%s_%s", $header, $filename_stripped;
        }

        if (defined($common_header_hash{$header}))
        {
            next;
        }
        
        if (!($header =~ /^Ratio /))
        {
            next;
        }
        if (!($header =~ /Normalized$/))
        {
            next;
        }
        
        $final_order_hash{$final_num_cols}{filename}    = $filename;
        $final_order_hash{$final_num_cols}{header}      = $header_new;
        $final_order_hash{$final_num_cols}{header_orig} = $header;
        $final_order_hash{$final_num_cols}{col}         = $i;
        
        $final_num_cols++;
    }
}


# not in common, Ratio (Normalized)
for ($f = 0; $f < @filename_array; $f++)
{
    $filename = $filename_array[$f];
    
    @col_array_orig = sort {$a<=>$b} keys %{$col_to_header_hash{$filename}};
    
    foreach $i (@col_array_orig)
    {
        $header = $col_to_header_hash{$filename}{$i};
        $header_new = $header;
        if (1 || defined($is_heavy_light_hash{$filename}))
        {
            $filename_stripped = $filename_stripped_hash{$filename};
            $header_new = sprintf "%s_%s", $header, $filename_stripped;
        }

        if (defined($common_header_hash{$header}))
        {
            next;
        }
        
        if (!($header =~ /^Ratio / && $header =~ /Normalized$/))
        {
            next;
        }
        
        $final_order_hash{$final_num_cols}{filename}    = $filename;
        $final_order_hash{$final_num_cols}{header}      = $header_new;
        $final_order_hash{$final_num_cols}{header_orig} = $header;
        $final_order_hash{$final_num_cols}{col}         = $i;
        
        $final_num_cols++;
    }
}


# summarize common headers, right before we print out intensities
# order merged common headers first
foreach $header (@common_header_array)
{
    # keep RowIdentifier and ModificationID as-is, without _summary
    $header_new = $header;
    if (!($header =~ /ModificationID/) &&
        !($header =~ /RowIdentifier/))
    {
        $header_new = sprintf "%s_summary", $header;
    }

    $final_order_hash{$final_num_cols}{filename}    = '';
    $final_order_hash{$final_num_cols}{header}      = $header_new;
    $final_order_hash{$final_num_cols}{header_orig} = $header;
    $final_order_hash{$final_num_cols}{col}         = '';
            
    $summary_final_col_hash{$final_num_cols}        = $header;

    $final_num_cols++;
}


# not in common, Intensity
for ($f = 0; $f < @filename_array; $f++)
{
    $filename = $filename_array[$f];
    
    @col_array_orig = sort {$a<=>$b} keys %{$col_to_header_hash{$filename}};
    
    foreach $i (@col_array_orig)
    {
        $header = $col_to_header_hash{$filename}{$i};
        $header_new = $header;
        if (1 || defined($is_heavy_light_hash{$filename}))
        {
            $filename_stripped = $filename_stripped_hash{$filename};
            $header_new = sprintf "%s_%s", $header, $filename_stripped;
        }

        if (defined($common_header_hash{$header}))
        {
            next;
        }
        
        if (!($header =~ /^Intensity /i))
        {
            next;
        }
        
        $final_order_hash{$final_num_cols}{filename}    = $filename;
        $final_order_hash{$final_num_cols}{header}      = $header_new;
        $final_order_hash{$final_num_cols}{header_orig} = $header;
        $final_order_hash{$final_num_cols}{col}         = $i;
        
        $final_num_cols++;
    }
}


# not in common, TMT-
for ($f = 0; $f < @filename_array; $f++)
{
    $filename = $filename_array[$f];
    
    @col_array_orig = sort {$a<=>$b} keys %{$col_to_header_hash{$filename}};
    
    foreach $i (@col_array_orig)
    {
        $header = $col_to_header_hash{$filename}{$i};
        $header_new = $header;

        if (!($header =~ /TMT-/))
        {
            next;
        }

        if (defined($common_header_hash{$header}))
        {
            next;
        }

        if (1 || defined($is_heavy_light_hash{$filename}))
        {
            $filename_stripped = $filename_stripped_hash{$filename};

            # HACK -- TMT pipelines expect TMT channel at very end
            if ($header =~ /(.*?)_(TMT-[0-9]{3}[NC]?)$/)
            {
                $sample  = $1;
                $channel = $2;

                $header_new = sprintf "%s_%s_%s",
                    $sample, $filename_stripped, $channel;
            }
            else
            {
                $header_new = sprintf "%s_%s", $header, $filename_stripped;
            }
        }

        $final_order_hash{$final_num_cols}{filename}    = $filename;
        $final_order_hash{$final_num_cols}{header}      = $header_new;
        $final_order_hash{$final_num_cols}{header_orig} = $header;
        $final_order_hash{$final_num_cols}{col}         = $i;
        
        $final_num_cols++;
    }
}


# not in common, IRON
for ($f = 0; $f < @filename_array; $f++)
{
    $filename = $filename_array[$f];
    
    @col_array_orig = sort {$a<=>$b} keys %{$col_to_header_hash{$filename}};
    
    foreach $i (@col_array_orig)
    {
        $header = $col_to_header_hash{$filename}{$i};
        $header_new = $header;
        if (1 || defined($is_heavy_light_hash{$filename}))
        {
            $filename_stripped = $filename_stripped_hash{$filename};
            $header_new = sprintf "%s_%s", $header, $filename_stripped;
        }

        if (defined($common_header_hash{$header}))
        {
            next;
        }
        
        if (!($header =~ /^IRON /))
        {
            next;
        }
        
        $final_order_hash{$final_num_cols}{filename}    = $filename;
        $final_order_hash{$final_num_cols}{header}      = $header_new;
        $final_order_hash{$final_num_cols}{header_orig} = $header;
        $final_order_hash{$final_num_cols}{col}         = $i;
        
        $final_num_cols++;
    }
}


# print final headers
$final_header_string = '';
for ($i = 0; $i < $final_num_cols; $i++)
{
    $header = $final_order_hash{$i}{header};
    
    if ($i)
    {
        $final_header_string .= "\t";
    }
    
    $final_header_string .= $header;
}
printf "%s\n", $final_header_string;


@modificationid_array = sort keys %modificationid_seen_hash;


# print final values
foreach $modificationid (@modificationid_array)
{
    $final_line = '';

    for ($c = 0; $c < $final_num_cols; $c++)
    {
      # summarize several fields together
      if (defined($summary_final_col_hash{$c}))
      {
        $num_fields = 0;
        @temp_array = ();
        
        $header = $summary_final_col_hash{$c};

        for ($f = 0; $f < @filename_array; $f++)
        {
            $filename = $filename_array[$f];
            $r = $modificationid_row_hash{$filename}{$modificationid};
            
            if (!defined($r))
            {
                next;
            }

            @data_array = @{$data_hash{$filename}{data}[$r]};

            @col_array = sort {$a<=>$b} keys
                          %{$header_to_multiple_col_hash{$filename}{$header}};
            
            foreach $col (@col_array)
            {
                $value = $data_array[$col];
                
                if (defined($value) && $value =~ /[A-Za-z0-9+]/)
                {
                    $temp_array[$num_fields++] = $value;
                }
            }
        }
        
        $value = '---';
        
        $global_sort_header = $header;

        if ($num_fields)
        {
            @temp_array = sort compare_fields @temp_array;
            $value = $temp_array[0];
        }
      }
      # regular data field
      else
      {
        $filename = $final_order_hash{$c}{filename};
        $r = $modificationid_row_hash{$filename}{$modificationid};
        
        if (!defined($r))
        {
            $value = '---';
        }
        else
        {
            @data_array = @{$data_hash{$filename}{data}[$r]};

            $col = $final_order_hash{$c}{col};
            $value = $data_array[$col];
        }
      }
      
      # convert zero intensity into missing values (---)
      if ($final_order_hash{$c}{header_orig} =~ /^Intensity / ||
          $final_order_hash{$c}{header_orig} =~ /TMT- /)
      {
          if (!defined($value) ||
              (is_number($value) && $value == 0))
          {
              $value = '---';
          }
      }

      if ($c)
      {
        $final_line .= "\t";
      }

      if (defined($value))
      {
        $final_line .= $value;
      }
    }

    printf "%s\n", $final_line;
}
