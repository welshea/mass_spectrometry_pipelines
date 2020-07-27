#!/usr/bin/perl -w


$pos_csv_filename = shift;
$neg_csv_filename = shift;
$output_root_name = shift;

if (!defined($pos_csv_filename) || !defined($neg_csv_filename) ||
    !defined($output_root_name))
{
    printf STDERR "Usage: automate_metabolomics.pl mzmine_pos.csv mzmine_neg.csv output_prefix\n";

    exit(1);
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


$cmd_str_strip_pos = sprintf "%s \"%s\" | %s - \"%s\" \"%s\" > \"%s\"",
                         'csv2tab_not_excel.pl',
                         $pos_csv_filename,
                         'strip_metabolomics_columns.pl',
                         $pos_unidentified_output_name,
                         $pos_spikeins_output_name,
                         $pos_cleaned_filename;

$cmd_str_strip_neg = sprintf "%s \"%s\" | %s - \"%s\" \"%s\" > \"%s\"",
                         'csv2tab_not_excel.pl',
                         $neg_csv_filename,
                         'strip_metabolomics_columns.pl',
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

$cmd_str_merge     = sprintf "%s \"%s\" \"%s\" > \"%s\"",
                         'merge_metabolomics_pos_neg.pl',
                         $pos_iron_filename,
                         $neg_iron_filename,
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
