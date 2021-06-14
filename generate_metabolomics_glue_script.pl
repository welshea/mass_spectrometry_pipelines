#!/usr/bin/perl -w

# 2021-06-14:  add annotate_metabolomics.pl script
# 2021-05-27:  add default output_root_name
# 2021-01-06:  add support for strip_metabolomics_columns.pl command options

use File::Spec;
use File::Basename;


$keep_single_pregap_flag   = 0;
$discard_unidentified_flag = 0;
$discard_heavy_flag        = 0;

$syntax_error_flag         = 0;
$num_files                 = 0;


# this correctly doesn't follow symlinked scripts, so you get the right path
# various other methods return the destination linked path instead,
#  which we don't want here
$script_path   = dirname(File::Spec->rel2abs(__FILE__));


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
$pos_iron_filename            = sprintf "%s_pos_iron.txt",
                                   $output_root_name;
$neg_unidentified_output_name = sprintf "%s_neg_unidentified.txt",
                                   $output_root_name;
$neg_spikeins_output_name     = sprintf "%s_neg_spikeins.txt",
                                   $output_root_name;
$neg_cleaned_filename         = sprintf "%s_neg_cleaned.txt",
                                   $output_root_name;
$neg_iron_filename            = sprintf "%s_neg_iron.txt",
                                   $output_root_name;
$merged_filename              = sprintf "%s_iron_merged.txt",
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


$cmd_str_strip_pos = sprintf "%s \"%s\" | %s%s - \"%s\" \"%s\" > \"%s\"",
                         'csv2tab_not_excel.pl',
                         $pos_csv_filename,
                         'strip_metabolomics_columns.pl',
                         $strip_options_str,
                         $pos_unidentified_output_name,
                         $pos_spikeins_output_name,
                         $pos_cleaned_filename;

$cmd_str_strip_neg = sprintf "%s \"%s\" | %s%s - \"%s\" \"%s\" > \"%s\"",
                         'csv2tab_not_excel.pl',
                         $neg_csv_filename,
                         'strip_metabolomics_columns.pl',
                         $strip_options_str,
                         $neg_unidentified_output_name,
                         $neg_spikeins_output_name,
                         $neg_cleaned_filename;

$cmd_str_iron_pos  = sprintf "%s --iron-exclusions=\"%s\" --iron-spikeins=\"%s\" \"%s\" > \"%s\"",
                         'iron_normalize_mass_spec.pl',
                         $pos_unidentified_output_name,
                         $pos_spikeins_output_name,
                         $pos_cleaned_filename,
                         $pos_iron_filename;

$cmd_str_iron_neg  = sprintf "%s --iron-exclusions=\"%s\" --iron-spikeins=\"%s\" \"%s\" > \"%s\"",
                         'iron_normalize_mass_spec.pl',
                         $neg_unidentified_output_name,
                         $neg_spikeins_output_name,
                         $neg_cleaned_filename,
                         $neg_iron_filename;

$cmd_str_merge     = sprintf "%s \"%s\" \"%s\" | %s \"%s/%s\" - > \"%s\"",
                         'merge_metabolomics_pos_neg.pl',
                         $pos_iron_filename,
                         $neg_iron_filename,
                         'annotate_metabolomics.pl',
                         $script_path,
                         'metabolite_database_latest.txt',
                         $merged_filename;

$cmd_str_rm        = sprintf "rm \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" \"%s\"",
                         $pos_unidentified_output_name,
                         $pos_spikeins_output_name,
                         $pos_cleaned_filename,
                         $pos_iron_filename,
                         $neg_unidentified_output_name,
                         $neg_spikeins_output_name,
                         $neg_cleaned_filename,
                         $neg_iron_filename;


print "$cmd_str_strip_pos\n";
print "$cmd_str_strip_neg\n";
print "$cmd_str_iron_pos\n";
print "$cmd_str_iron_neg\n";
print "$cmd_str_merge\n";
print "$cmd_str_rm\n";
