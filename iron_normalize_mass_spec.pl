#!/usr/bin/perl -w

use Scalar::Util qw(looks_like_number);
use File::Basename;

# 2026-02-20:  correctly handle enclosing double quotes
# 2026-01-08:  skip normalization of common RNA-Seq annotation headers
# 2025-04-01:  detect ORIEN Avatar SL###### format sample identifiers
# 2025-01-17:  preserve transcript* fields
# 2024-09-12:  preserve Symbol field
# 2024-09-12:  remove leftover STDERR debugging
# 2024-09-09:  rescale output to global raw mean, instead of per-row raw mean
# 2024-02-28:  add --keep-zero-rows flag
# 2024-02-27:  if no sample columns detected, default to all columns
# 2024-02-27:  add --iron-ignore-low to untilt flags, --no-ignore-low option
# 2023-11-13:  support iTRAQ-4 and iTRAQ-8
# 2023-10-05:  add --iron-ignore-low to untilt flags, for IRON >= v2.3.0
# 2023-08-07:  --iron-untilt overrides --no-iron flag if specified afterwards
# 2023-08-07:  change --norm-none to --no-iron for consistency with other .pl
# 2023-06-27:  update is_number() to not treat NaNs as numbers
# 2023-05-25:  refactor how stripping sample names is handled
# 2023-05-25:  improved pos/neg detection from merge_metabolomics_pos_neg.pl
# 2023-05-25:  support alternative lipidomics sample headers
# 2023-05-18:  strip UTF-8 BOM from MaxQuant 2.4 output, which broke many things
# 2023-04-18:  support --first-ch flag
# 2023-04-10:  support --boost and --last-ch flags
# 2023-03-24:  add preliminary FragPipe support
# 2023-01-11:  rescale normalized rows so log2 normalized mean == log2 raw mean
# 2022-12-09:  change |& to 2>&1| to fix Ubuntu/Debian sh --> dash
# 2022-04-26:  add --iron-untilt to usage statement (no longer experimental)
# 2022-04-21:  begin adding support for Proteome Discoverer
# 2022-04-21:  remove all double quotes, not just the first one per line (oops)
# 2022-04-21:  make sure split doesn't remove empty trailing fields
# 2021-12-22:  add experimental --iron-untilt flag
# 2021-08-08:  print usage statement on command line error
# 2021-04-14:  begin adding lipidomics support
# 2020-09-24:  add --norm-none flag (has some use for metabolomics pipeline)
# 2020-08-26:  attempt to deal with more pos/neg sample labeling typos
# 2020-08-24:  remove --bg background subtraction option, as it is a bad idea
# 2020-08-03:  output scaling factor lines in correct TMT label N/C order
# 2020-07-27:  hack in support for El-MAVEN metabolomics samples
#
# 2020-02-28:  added --no-log2 --log2 --unlog2 flags
#              removed old --unlog flag and reimplemented it for --unlog2
#              prioritize ModificationID over RowIdentifier column for id
#               this is to catch cases where cleanup is done in the wrong order
#
# 2020-03-12:  skip "Intensity L" and "Intensity H" summary columns
#              warning, will not correctly handle case where there are real
#               samples named H and L together with non-H/L sample names
#
# 2020-03-16:  fixed support for TMT-###[NC] columns with text in front of them
# 2020-03-20:  don't use .mzXML for identifiying mzmine output anymore
#              re-enable pre-pending IRON to metabolomics output
#
# 2020-04-22:  added support for Skyline " Total Area" sample headers

# all functions are effectively macros, since there are no local variables

# we assume that there can be no more than 2 pools per plex

$headers_lc_to_skip_hash{'ensembl_gene_id'} = 1;
$headers_lc_to_skip_hash{'gencodeid'} = 1;
$headers_lc_to_skip_hash{'gene'} = 1;
$headers_lc_to_skip_hash{'order_in_genome'} = 1;
$headers_lc_to_skip_hash{'kept_in_filtered_gtf'} = 1;
$headers_lc_to_skip_hash{'exist_in_ensembl'} = 1;
$headers_lc_to_skip_hash{'hgnc_gene_id'} = 1;
$headers_lc_to_skip_hash{'entrez_gene_id'} = 1;
$headers_lc_to_skip_hash{'symbol'} = 1;
$headers_lc_to_skip_hash{'alias'} = 1;
$headers_lc_to_skip_hash{'type'} = 1;
$headers_lc_to_skip_hash{'type_gene'} = 1;
$headers_lc_to_skip_hash{'location'} = 1;
$headers_lc_to_skip_hash{'description'} = 1;
$headers_lc_to_skip_hash{'genetype'} = 1;
$headers_lc_to_skip_hash{'genename'} = 1;
$headers_lc_to_skip_hash{'transcript'} = 1;


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


sub read_in_data_file
{
    $infile = $_[0];
    open INFILE, "$infile" or die "ABORT -- can't open $infile\n";

    # read in header line
    $line = <INFILE>;
    $line =~ s/[\r\n]+//g;

    # remove UTF-8 byte order mark, since it causes many problems
    # remove some other BOM I've observed in the wild
    $line =~ s/(^\xEF\xBB\xBF|^\xEF\x3E\x3E\xBF)//;

    @array = split /\t/, $line, -1;
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\s+/ /g;

        # handle enclosing quotes
        if ($array[$i] =~ /^\".*\"$/)
        {
            $array[$i] =~ s/^\"(.*)\"$/$1/;
            $array[$i] =~ s/\"\"/\"/g;
        }

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
        # check for lipidomics abundances
        if ($field =~ / Peak height$/i ||
            $field =~ / Peak area$/i ||
            $field =~ /^(Area|Height),/i ||
            $field =~ /^(Area|Height)\[.*\]$/i)
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

        # Proteome Discoverer
        if ($field =~ /^Abundances \(Grouped\):\s*/)
        {
            $sample_to_file_col_hash{$field} = $i;
            $sample_array[$num_samples++] = $field;
        }
        

        # check to see if it is cleaned up TMT column headers
        if ($field =~ /^(TMT|iTRAQ)/i)
        {
            # skip individual fractions
            if ($field =~ /___\d+$/)
            {
                next;
            }

            $sample_to_file_col_hash{$field} = $i;
            $sample_array[$num_samples++] = $field;
        }
        elsif ($field =~ /(TMT|iTRAQ)\-\d+[NC]*$/i)
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
        
        # skip LipidSearch Rt columns, since they can contain pos/neg
        if ($field =~ /^Rt\[/)
        {
            next;
        }

        if ($field =~ /(^|[^A-Za-z0-9]+)(pos|neg)([^A-Za-z0-9]+|$)/i ||
            $field =~ /[^A-Za-z0-9](pos|neg)[0-9]{0,1}(?:\]*)$/i ||
            $field =~ /^(pos|neg)[^A-Za-z0-9]/)
        {
            $sample_to_file_col_hash{$field} = $i;
            $sample_array[$num_samples++] = $field;

            $global_metabolomics_flag = 1;
        }
      }
    }

    # still didn't find any, maybe it is FragPipe " Intensity$" headers
    if ($num_samples == 0)
    {
        for ($i = 0; $i < @array; $i++)
        {
            $field = $array[$i];

            if ($field =~ /\s+Intensity$/i &&
                !($field =~ /\s+MaxLFQ\s+Intensity$/i))
            {
                $sample_to_file_col_hash{$field} = $i;
                $sample_array[$num_samples++] = $field;
                
                $intensity_at_end_flag = 1;
            }
        }
    }


    # still didn't find any, maybe it is ORIEN Avatar SLIDs
    if ($num_samples == 0)
    {
        for ($i = 0; $i < @array; $i++)
        {
            $field = $array[$i];

            if ($field =~ /^SL[0-9]{6}$/)
            {
                $sample_to_file_col_hash{$field} = $i;
                $sample_array[$num_samples++] = $field;
                
                $intensity_at_end_flag = 1;
            }
        }
    }
    
    # still didn't find any, print warning and assume all columns are data
    if ($num_samples == 0)
    {
        printf STDERR "WARNING -- no sample columns detected, using all cols after identifier\n";
    
        for ($i = $id_col + 1; $i < @array; $i++)
        {
            $field    = $array[$i];
            $field_lc = lc $field;
            
            if (defined($headers_lc_to_skip_hash{$field_lc}) ||
                $field =~ /^Accession/i)
            {
                printf STDERR "Treating %s column as non-data\n", $field;
                next;
            }

            $sample_to_file_col_hash{$field} = $i;
            $sample_array[$num_samples++] = $field;
                
            $intensity_at_end_flag = 1;
        }
    }
    
    # strip header stuff from sample names
    @sample_strip_array = @sample_array;
    if ($strip_sample_names_flag && $num_samples)
    {
        for ($i = 0; $i < $num_samples; $i++)
        {
            $sample = $sample_array[$i];
            
            # strip Intensity from beginning
            if ($sample =~ /^Intensity\s+/)
            {
                $sample =~ s/^Intensity\s+//i;
            }
            
            # strip Intensity from end
            elsif ($intensity_at_end_flag)
            {
                $sample =~ s/\s+Intensity$//i;
            }

            # strip mzMine stuff, add IRON
            elsif ($sample =~ / Peak height$/i ||
                   $sample =~ / Peak area$/i)
            {
                $sample =~ s/\.mzX?ML[^.]+$//i;
                $sample =~ s/ Peak \S+$//i;
                
            }
            
            # strip lipidomics stuff
            elsif ($sample =~ /^Height,/i ||
                   $sample =~ /^Area,/i)
            {
                $sample =~ s/^(Area|Height),\s*//i;
            }
            elsif ($sample =~ /^(?:Area|Height)\[(.*)\]$/i)
            {
                $sample =~ s/^(?:Area|Height)\[(.*)\]$/$1/i;
            }

            # strip Skyline stuff
            elsif ($sample =~ / Total Area$/i)
            {
                $sample =~ s/ Total Area$//i;
            }

            # strip Proteome Discoverer stuff
            elsif ($sample =~ /^Abundances \(Grouped\):\s*/i)
            {
                $sample =~ s/^Abundances \(Grouped\):\s*//i;
            }
            
            $sample_strip_array[$i] = $sample;
        }
    }


    $num_header_cols = @array;
    $header_line = join "\t", @array;

    # read in the rest of the data
    $row = 0;
    while(defined($line=<INFILE>))
    {
        $line =~ s/[\r\n]+//g;

        # be sure to not skip trailing blank fields
        @array = split /\t/, $line, -1;

        for ($i = 0; $i < @array; $i++)
        {
            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;
            $array[$i] =~ s/\s+/ /g;

            # handle enclosing quotes
            if ($array[$i] =~ /^\".*\"$/)
            {
                $array[$i] =~ s/^\"(.*)\"$/$1/;
                $array[$i] =~ s/\"\"/\"/g;
            }

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
        $sample = $sample_strip_array[$i];
        $sample =~ s/^Intensity\s+//i;

        # strip Intensity from the emd
        if ($intensity_at_end_flag)
        {
            $sample =~ s/\s+Intensity$//i;
        }
        
        # strip mzMine stuff
        $sample =~ s/\.mzX?ML[^.]+$//i;
        $sample =~ s/ Peak \S+$//i;

        # strip Skyline stuff
        $sample =~ s/ Total Area$//i;
        
        # strip Proteome Discoverer stuff
        $sample =~ s/^Abundances \(Grouped\):\s*//i;
        
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
    elsif ($boost_flag)
    {
        $median_sample = $sample_strip_array[$num_samples - 1];
    }
    elsif ($first_flag)
    {
        $median_sample = $sample_strip_array[0];
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
    
    # override --proteomics or --rnaseq by putting it afterwards
    if ($ignore_low_flag)
    {
        # use for all mass spec data, where values < 1 are junk
        $ignore_low_str = '--iron-ignore-low';
    }
    else
    {
        # use for TPM/FPKM RNA-Seq data, values 0.00001 < x < 1 are OK
        $ignore_low_str = '--iron-no-ignore-low';
    }

    # normalize data as usual
    $cmd_string = sprintf "iron_generic --proteomics --norm-iron=\"%s\" %s %s %s %s \"%s\" -o \"%s\" 2>&1| grep -P \"^Global(Scale|FitLine)\"",
        $median_sample, $exclusions_string, $spikeins_string, $bg_string,
        $ignore_low_str, $iron_input_name, $iron_output_name;

    if ($iron_untilt_flag)
    {
        $cmd_string =~ s/--proteomics/--rnaseq/g;
    }


    @scale_lines = `$cmd_string`;
    @scale_lines = sort cmp_scale_lines @scale_lines;

    # output scaling factor file
    open OUTPUT_SCALING, ">$scales_output_name" or die "can't open scales output file $scales_output_name";
    for ($i = 0; $i < @scale_lines; $i++)
    {
        $line = $scale_lines[$i];
        $line =~ s/[\r\n]+//g;
        
        # strip GlobalScale: from beginning
        $line =~ s/^Global(Scale|FitLine):\s+//ig;
        
        printf OUTPUT_SCALING "%s\n", $line;
    }
    close OUTPUT_SCALING;

    # we didn't actually want to normalize the data, we only iron'd it to
    # get the scaling factors
    if ($no_iron_flag)
    {
        printf STDERR "--no-iron; overwriting iron output with unnormalized data\n";
    
        $cmd_string = sprintf "iron_generic --proteomics --norm-none \"%s\" -o \"%s\"",
            $iron_input_name, $iron_output_name;
        
        `$cmd_string`;
    }


    # for rescaling
    $mean_raw  = 0;
    $n_raw     = 0;
    $mean_norm = 0;
    $n_norm    = 0;
    $offset    = 0;


    # open IRON output for reading
    open OUTPUT_FOR_IRON, "$iron_output_name" or die "ABORT -- can't open $iron_output_name\n";
    
    # skip header line
    $line = <OUTPUT_FOR_IRON>;


    # read in the rest of the data into a new log2 array
    $row = 0;
    while(defined($line=<OUTPUT_FOR_IRON>))
    {
        $line =~ s/[\r\n]+//g;

        @array = split /\t/, $line, -1;

        for ($i = 0; $i < @array; $i++)
        {
            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;
            $array[$i] =~ s/\s+/ /g;

            # handle enclosing quotes
            if ($array[$i] =~ /^\".*\"$/)
            {
                $array[$i] =~ s/^\"(.*)\"$/$1/;
                $array[$i] =~ s/\"\"/\"/g;
            }

            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;
            $array[$i] =~ s/\s+/ /g;
        }
        
        # mean of original log2 values
        for ($i = 0; $i < $num_samples; $i++)
        {
            $sample = $sample_array[$i];
            $col    = $sample_to_file_col_hash{$sample};
            
            $value  = $orig_data_array[$row][$col];
            
            if (defined($value) && is_number($value) && $value > 0)
            {
                $mean_raw += log($value);
                $n_raw++;
            }
        }

        for ($i = 0; $i < $num_samples; $i++)
        {
            $value = $array[$i+1];
            
            if (defined($value) && is_number($value) && $value > 0)
            {
                $mean_norm += log($value);
                $n_norm++;
            }
        }

        for ($i = 0; $i < $num_samples; $i++)
        {
            $value = $array[$i+1];
            
            # paranoia for log taking
            if (is_number($value) && $value > 0)
            {
                $iron_normalized_array[$row][$i] = $value;
            }
            else
            {
                $iron_normalized_array[$row][$i] = '';
            }
        }
        
        $row++;
    }
    $num_rows = $row;
    close OUTPUT_FOR_IRON;


    # amount to shift normalized data so it has same log2 mean as raw
    $offset = 0;
    if ($n_raw)
    {
        $mean_raw  /= $n_raw  * $log2;
    }
    if ($n_norm)
    {
        $mean_norm /= $n_norm * $log2;
    }
    if ($mean_raw > 0 && $mean_norm > 0)
    {
        $offset = $mean_raw - $mean_norm;
    }
    
    
    # do not scale the data after all
    if ($rescale_flag == 0)
    {
        $offset = 0;
    }
    

    # rescale the data
    for ($row = 0; $row < $num_rows; $row++)
    {
        for ($i = 0; $i < $num_samples; $i++)
        {
            $value = $iron_normalized_array[$row][$i];

            # paranoia for log taking
            if (is_number($value) && $value > 0)
            {
                $value = log($value) / $log2 + $offset;
            }
            else
            {
                $value = '';
            }
            
            $iron_normalized_array[$row][$i] = $value;
        }
    }
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
        $sample = $sample_strip_array[$i];

        # we're going to merge POS and NEG together later
        # prepend IRON so that it is easy to detect sample columns
        if ($global_metabolomics_flag)
        {
            $sample = 'IRON ' . $sample;
        }

        printf "\t%s", $sample;
    }
    printf "\n";

    # print the rest of the table
    for ($row = 0; $row < $num_rows; $row++)
    {
        @array = ();

        if ($keep_zero_rows_flag == 0)
        {
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
            
            # in case the line removed trailing blank fields
            if (!defined($value))
            {
                $value = '';
            }
            
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
$unlog2_flag             = 0;   # unlog2 the input data
$no_log2_flag            = 0;   # do not log2 the output data
$no_iron_flag            = 0;   # do not normalize the data at all
$iron_untilt_flag        = 0;   # --rnaseq flag
$boost_flag              = 0;   # assume last sample is boosting channel
$first_flag              = 0;   # assume first sample is pool
$ignore_low_flag         = 1;   # ignore values < 1 during training
$keep_zero_rows_flag     = 0;   # keep all-zero rows
$global_metabolomics_flag = 0;
$rescale_flag            = 1;   # rescale normalized data to mean raw data


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
        elsif ($field =~ /^--no-iron$/)
        {
            $no_iron_flag = 1;
        }
        elsif ($field =~ /^--iron-untilt$/)
        {
            $iron_untilt_flag = 1;
            $no_iron_flag = 0;
        }
        elsif ($field =~ /^--boost$/ ||
               $field =~ /^--last-ch$/)
        {
            $boost_flag = 1;
            $first_flag = 0;
        }
        elsif ($field =~ /^--first-ch$/)
        {
            $first_flag = 1;
            $boost_flag = 0;
        }
        # include values 0.00001 < x < 1 during normalization training
        # used for TPM/FPKM RNA-Seq data
        elsif ($field =~ /^--no-ignore-low/)
        {
            $ignore_low_flag = 0;
        }
        elsif ($field =~ /^--keep-zero-rows/)
        {
            $keep_zero_rows_flag = 1;
        }
        elsif ($field =~ /^--rescale/)
        {
            $rescale_flag = 1;
        }
        elsif ($field =~ /^--no-rescale/)
        {
            $rescale_flag = 0;
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

if ($syntax_error_flag || $num_files == 0)
{
    $program_name = basename($0);

    printf STDERR "Usage: $program_name [options] input_file.txt [reference_sample]\n";
    printf STDERR "\n";
    printf STDERR "  Options:\n";
    printf STDERR "    --iron-untilt                 account for relative dynamic range\n";
    printf STDERR "    --iron-exclusions=filename    identifiers to exclude from training\n";
    printf STDERR "    --iron-spikeins=filename      list of spikein identifiers\n";
    printf STDERR "\n";
    printf STDERR "    --rescale                     rescale output to raw log2 mean (default)\n";
    printf STDERR "    --no-rescale                  do not rescale normalized output\n";
    printf STDERR "\n";
    printf STDERR "    --keep-zero-rows              keep rows with all-zero data\n";
    printf STDERR "    --log2                        output log2 abundances [default]\n";
    printf STDERR "    --no-strip-sample-names       keep original full sample headers\n";
    printf STDERR "    --no-log2                     output unlogged abundances\n";
    printf STDERR "    --no-ignore-low               *must use* for RNA-Seq TPM/FPKM data\n";
    printf STDERR "    --no-iron                     disable normalization\n";
    printf STDERR "    --unlog2                      exponentiate input data, pow(2, value)\n";
    printf STDERR "\n";
    printf STDERR "    --boost                       use last  sample for normalization\n";
    printf STDERR "    --last-ch                     use last  sample for normalization\n";
    printf STDERR "    --first-ch                    use first sample for normalization\n";

    exit(1);
}


$process_id = $$;

$scales_output_name =  $filename;
$scales_output_name =~ s/\.txt$//i;
$findmedian_output_name = sprintf "%s_findmedian.txt", $scales_output_name;
$scales_output_name = sprintf "%s_scaling_factors.txt", $scales_output_name;

$iron_input_name  = sprintf "temp_iron_input_%s.txt", $process_id;
$iron_output_name = sprintf "temp_iron_output_%s.txt", $process_id;

# set if samples have Intensity at the end of the column headers
$intensity_at_end_flag = 0;

read_in_data_file($filename);
store_condensed_data();
iron_samples();
output_final_data();

`rm \"$iron_input_name\"`;
`rm \"$iron_output_name\"`;
