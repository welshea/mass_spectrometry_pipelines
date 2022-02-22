#!/usr/bin/perl -w

# 2021-02-22:  detect whether each input file is csv or tab-delimited
# 2021-12-08:  fix broken spaces in filename support after lipidomics changes
# 2021-11-30:  add support for lipidomics
# 2021-08-19:  change default back to leaving heavy unscaled
# 2021-08-19:  --scale-heavy --no-scale-heavy to --heavy-tracer --heavy-spikein
# 2021-08-19:  add --norm-none flag to disable normalization
# 2021-08-18:  add --scale-heavy --no-scale-heavy flags
# 2021-08-18:  add --ppm flag to set m/z ppm tolerance
# 2021-08-18:  implement new --scale-heavy and --no-scale-heavy flags
# 2021-08-17:  changed default behavior to scale heavy labeled metabolites;
#              it was too easy to mess up and forget to scale them in tracer
#              experiments
# 2021-08-11:  add QC script, delete sample, scaling factor, findmedian files
# 2021-08-06:  support naming of output merged sample table
# 2021-07-23:  add _log2 to iron output filenames to indicate it is log2
# 2021-06-14:  add annotate_metabolomics.pl script
# 2021-05-27:  add default output_root_name
# 2021-01-06:  add support for strip_metabolomics_columns.pl command options

use Scalar::Util qw(looks_like_number);
use File::Spec;
use File::Basename;


$keep_single_pregap_flag   = 0;
$discard_unidentified_flag = 0;
$discard_heavy_flag        = 0;
$scale_heavy_flag          = 0;    # control normalization of heavy rows
$mz_tol_ppm                = '';   # '' means no --ppm argument specified
$norm_none_flag            = 0;    # disable normalization

$syntax_error_flag         = 0;
$num_files                 = 0;


# this correctly doesn't follow symlinked scripts, so you get the right path
# various other methods return the destination linked path instead,
#  which we don't want here
$script_path   = dirname(File::Spec->rel2abs(__FILE__));


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




# begin main()

for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        if ($field =~ /^--keep-single-pregap$/)
        {
            $keep_single_pregap_flag = 1;
        }
        elsif ($field =~ /^--discard-single-pregap$/)
        {
            $keep_single_pregap_flag = 0;
        }
        elsif ($field =~ /^--discard-unidentified$/)
        {
            $discard_unidentified_flag = 1;
        }
        elsif ($field =~ /^--discard-heavy$/)
        {
            $discard_heavy_flag = 1;
        }
        elsif ($field =~ /^--heavy-tracer$/)
        {
            $scale_heavy_flag = 1;
        }
        elsif ($field =~ /^--heavy-spikein$/)
        {
            $scale_heavy_flag = 0;
        }
        elsif ($field =~ /^--norm-none$/)
        {
            $norm_none_flag = 1;
        }
        # override default PPM tolerance
        elsif ($field eq '--ppm' ||
               $field =~ /^--ppm\=/)
        {
            $ppm_arg = '';
            
            if ($field eq '--ppm')
            {
                if ($i + 1 < @ARGV)
                {
                    $i++;
                    $ppm_arg = $ARGV[$i];
                }
            }
            elsif ($field =~ /^--ppm=(.*)/)
            {
                $ppm_arg = $1;
            }
            
            if (is_number($ppm_arg))
            {
                $mz_tol_ppm = $ppm_arg;
            }
            else
            {
                printf STDERR "ABORT -- you must specify a numeric value after --ppm\n";
                $syntax_error_flag = 1;
            }
        }
        else
        {
            printf STDERR "ABORT -- unknown option %s\n", $field;
            $syntax_error_flag = 1;
        }
    }
    else
    {
        if ($num_files == 0)
        {
            $pos_csv_filename = $field;
            $num_files++;
        }
        elsif ($num_files == 1)
        {
            $neg_csv_filename = $field;
            $num_files++;
        }
        elsif ($num_files == 2)
        {
            $output_root_name = $field;
            $num_files++;
        }
    }
}


if ($syntax_error_flag ||
    !defined($pos_csv_filename) || !defined($neg_csv_filename))
{
    printf STDERR "Usage: generate_metabolomics_glue_script.pl [options] mzmine_pos.csv mzmine_neg.csv output_prefix\n";
    printf STDERR "\n";
    printf STDERR "Options:\n";
    printf STDERR "    --heavy-spikein            heavy rows are spikeins, leave unscaled (default)\n";
    printf STDERR "    --heavy-tracer             heavy rows are biological, normalize them\n";
    printf STDERR "    --norm-none                disable normalization; use on targeted panels\n";
    printf STDERR "    --ppm N                    override default m/z PPM tolerance\n";
    printf STDERR "\n";
    printf STDERR "    --discard-heavy            discard heavy labeled rows\n";
    printf STDERR "    --discard-unidentified     discard unidentified rows\n";
    printf STDERR "\n";
    printf STDERR "  options which use the MZmine \"row number of detected peaks\" column:\n";
    printf STDERR "    --discard-single-pregap    discard pre gap-filled single-hit rows (default)\n";
    printf STDERR "    --keep-single-pregap       keep pre gap-filled single-hit rows\n";


    exit(1);
}


if (!defined($output_root_name))
{
    $output_root_name = 'metabolomics_pipeline';
}


#$process_id = $$;


$pos_unidentified_output_name = sprintf "%s_pos_unidentified.txt",
                                   $output_root_name;
$pos_spikeins_output_name     = sprintf "%s_pos_spikeins.txt",
                                   $output_root_name;
$pos_cleaned_filename         = sprintf "%s_pos_cleaned.txt",
                                   $output_root_name;
$pos_iron_filename            = sprintf "%s_pos_iron_log2.txt",
                                   $output_root_name;
$neg_unidentified_output_name = sprintf "%s_neg_unidentified.txt",
                                   $output_root_name;
$neg_spikeins_output_name     = sprintf "%s_neg_spikeins.txt",
                                   $output_root_name;
$neg_cleaned_filename         = sprintf "%s_neg_cleaned.txt",
                                   $output_root_name;
$neg_iron_filename            = sprintf "%s_neg_iron_log2.txt",
                                   $output_root_name;
$merged_filename              = sprintf "%s_iron_log2_merged.txt",
                                   $output_root_name;
$sample_table_filename        = sprintf "%s_sample_table.txt",
                                   $output_root_name;
$sample_qc_filename           = sprintf "%s_sample_qc_table.txt",
                                   $output_root_name;


$pos_sf_filename              = sprintf "%s_pos_cleaned_scaling_factors.txt",
                                   $output_root_name;
$neg_sf_filename              = sprintf "%s_neg_cleaned_scaling_factors.txt",
                                   $output_root_name;
$pos_fm_filename              = sprintf "%s_pos_cleaned_findmedian.txt",
                                   $output_root_name;
$neg_fm_filename              = sprintf "%s_neg_cleaned_findmedian.txt",
                                   $output_root_name;


# options to pass to strip_metabolomics_columns.pl
$strip_options_str = '';
if ($keep_single_pregap_flag)
{
    $strip_options_str .= ' --keep-single-pregap';
}
else
{
    $strip_options_str .= ' --discard-single-pregap';
}
if ($discard_unidentified_flag)
{
    $strip_options_str .= ' --discard-unidentified';
}
if ($discard_heavy_flag)
{
    $strip_options_str .= ' --discard-heavy';
}
if ($scale_heavy_flag)
{
    $strip_options_str .= ' --heavy-tracer';
}
else
{
    $strip_options_str .= ' --heavy-spikein';
}


# detect lipidomics data
$lipidomics_flag = 0;
$cmd_str = "head -1 \"$pos_csv_filename\"";
$result_str = `$cmd_str`;
if ($result_str =~ /LipidIon/ ||
    $result_str =~ /FattyAcid/)
{
    $lipidomics_flag = 1;
}


# read in each file, determine whether it is comma or tab-delimited
$count_comma_pos = 0;
$count_comma_neg = 0;
$count_tab_pos   = 0;
$count_tab_neg   = 0;

open INFILE, "$pos_csv_filename" or die "can't open file $pos_csv_filename\n";
while (defined($line=<INFILE>))
{
    $line_keep_comma = $line;
    $line_keep_tab   = $line;
    
    $line_keep_comma =~ s/[^,]//g;
    $line_keep_tab   =~ s/[^\t]//g;
    
    $count_comma_pos += length $line_keep_comma;
    $count_tab_pos   += length $line_keep_tab;
}
close INFILE;

open INFILE, "$neg_csv_filename" or die "can't open file $neg_csv_filename\n";
while (defined($line=<INFILE>))
{
    $line_keep_comma = $line;
    $line_keep_tab   = $line;
    
    $line_keep_comma =~ s/[^,]//g;
    $line_keep_tab   =~ s/[^\t]//g;
    
    $count_comma_neg += length $line_keep_comma;
    $count_tab_neg   += length $line_keep_tab;
}
close INFILE;


$csv_command_pos = 'cat';
$csv_command_neg = 'cat';
if ($count_comma_pos > $count_tab_pos)
{
    $csv_command_pos = 'csv2tab_not_excel.pl';
}
if ($count_comma_pos > $count_tab_pos)
{
    $csv_command_neg = 'csv2tab_not_excel.pl';
}


$cmd_str_strip_pos = sprintf "%s \"%s\" | %s%s - \"%s\" \"%s\" > \"%s\"",
                         $csv_command_pos,
                         $pos_csv_filename,
                         'strip_metabolomics_columns.pl',
                         $strip_options_str,
                         $pos_unidentified_output_name,
                         $pos_spikeins_output_name,
                         $pos_cleaned_filename;

$cmd_str_strip_neg = sprintf "%s \"%s\" | %s%s - \"%s\" \"%s\" > \"%s\"",
                         $csv_command_neg,
                         $neg_csv_filename,
                         'strip_metabolomics_columns.pl',
                         $strip_options_str,
                         $neg_unidentified_output_name,
                         $neg_spikeins_output_name,
                         $neg_cleaned_filename;

$norm_extra_options_str = '';
if ($norm_none_flag)
{
    # disable normalization
    $norm_extra_options_str = '--norm-none';
}

# do not treat heavy labeled rows as spike-ins, normalize them
if ($scale_heavy_flag)
{
    $cmd_str_iron_pos  = sprintf "%s %s --iron-exclusions=\"%s\" \"%s\" > \"%s\"",
                             'iron_normalize_mass_spec.pl',
                             $norm_extra_options_str,
                             $pos_unidentified_output_name,
                             $pos_cleaned_filename,
                             $pos_iron_filename;

    $cmd_str_iron_neg  = sprintf "%s %s --iron-exclusions=\"%s\" \"%s\" > \"%s\"",
                             'iron_normalize_mass_spec.pl',
                             $norm_extra_options_str,
                             $neg_unidentified_output_name,
                             $neg_cleaned_filename,
                             $neg_iron_filename;
}
# treat heavy labeled rows as spike-ins, do not scale them
else
{
    $cmd_str_iron_pos  = sprintf "%s %s --iron-exclusions=\"%s\" --iron-spikeins=\"%s\" \"%s\" > \"%s\"",
                             'iron_normalize_mass_spec.pl',
                             $norm_extra_options_str,
                             $pos_unidentified_output_name,
                             $pos_spikeins_output_name,
                             $pos_cleaned_filename,
                             $pos_iron_filename;

    $cmd_str_iron_neg  = sprintf "%s %s --iron-exclusions=\"%s\" --iron-spikeins=\"%s\" \"%s\" > \"%s\"",
                             'iron_normalize_mass_spec.pl',
                             $norm_extra_options_str,
                             $neg_unidentified_output_name,
                             $neg_spikeins_output_name,
                             $neg_cleaned_filename,
                             $neg_iron_filename;
}


$annotate_options_str = '';
if ($mz_tol_ppm ne '')
{
    $annotate_options_str = '--ppm ' . $mz_tol_ppm;
}

$cmd_str_merge     = sprintf "%s \"%s\" \"%s\" %s | %s %s \"%s/%s\" - > \"%s\"",
                         'merge_metabolomics_pos_neg.pl',
                         $pos_iron_filename,
                         $neg_iron_filename,
                         $sample_table_filename,
                         'annotate_metabolomics.pl',
                         $annotate_options_str,
                         $script_path,
                         'metabolite_database_latest.txt',
                         $merged_filename;

$cmd_str_qc        = sprintf "%s \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" > \"%s\"",
                         'metabolomics_qc.pl',
                         $sample_table_filename,
                         $pos_sf_filename,
                         $neg_sf_filename,
                         $pos_fm_filename,
                         $neg_fm_filename,
                         $sample_qc_filename;

$cmd_str_rm        = sprintf "rm \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" \"%s\"",
                         $pos_unidentified_output_name,
                         $pos_spikeins_output_name,
                         $pos_cleaned_filename,
                         $pos_iron_filename,
                         $neg_unidentified_output_name,
                         $neg_spikeins_output_name,
                         $neg_cleaned_filename,
                         $neg_iron_filename,
                         $sample_table_filename,
                         $pos_sf_filename,
                         $neg_sf_filename,
                         $pos_fm_filename,
                         $neg_fm_filename,


print "$cmd_str_strip_pos\n";
print "$cmd_str_strip_neg\n";
print "$cmd_str_iron_pos\n";
print "$cmd_str_iron_neg\n";
print "$cmd_str_merge\n";
print "$cmd_str_qc\n";
print "$cmd_str_rm\n";
