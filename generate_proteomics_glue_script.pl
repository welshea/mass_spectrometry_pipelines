#!/usr/bin/perl -w

# 2020-02-28: clean phosphosites before cleaning TMT
#              this allows clean TMT to detect phospho data and use
#              the ModificationID column as the RowIdentifier

use File::Spec;
use File::Basename;

$input_filename = shift;
$species = shift;    # must be 'human' or 'mouse' for the time being
$autodetect_filename = shift;

if (!defined($species))
{
    $species = '';
}

# this correctly doesn't follow symlinked scripts, so you get the right path
# various other methods return the destination linked path instead,
#  which we don't want here
$script_path   = dirname(File::Spec->rel2abs(__FILE__));
$file_path     = dirname(File::Spec->rel2abs($input_filename));
$filename_base = basename($input_filename);

#$pid = $$;
#$tmp_path = sprintf "/tmp/proteomics_pipeline_%s/", $pid;

#`mkdir $tmp_path`;

printf STDERR "script_path: %s\n", $script_path;
printf STDERR "file_path:   %s\n", $file_path;
printf STDERR "basename:    %s\n", $filename_base;

if (!defined($autodetect_filename))
{
    $cmd_str = sprintf "\"%s/%s\" \"%s/%s\"",
        $script_path, 'autodetect_ms_data.pl',
        $file_path, $filename_base;
}
else
{
    $cmd_str = sprintf "cat \"%s\"", $autodetect_filename;
}

$autodetect_output = `$cmd_str`;

@autodetect_line_array = split /\n/, $autodetect_output;
foreach $line (@autodetect_line_array)
{
    $line =~ s/[\r\n]//g;
    
    @array = split /\t/, $line;
    
    if (@array == 2)
    {
        $key   = $array[0];
        $value = $array[1];
        
        $autodetect_hash{$key} = $value;
    }
}

@autodetect_key_array = sort keys %autodetect_hash;

foreach $key (@autodetect_key_array)
{
    $value = $autodetect_hash{$key};
    
    printf STDERR "%s\t%s\n", $key, $value;
}


# initial command which will pipe the data into additional future commands
$pipeline_preprocess_str = sprintf "\"%s/%s\" \"%s/%s\"",
    $script_path, 'strip_maxquant_columns.pl',
    $file_path, $filename_base;


if (defined($autodetect_hash{Modification}) &&
            $autodetect_hash{Modification} eq 'yes')
{
    $pipeline_preprocess_str .= sprintf " | \"%s/%s\" -",
        $script_path, 'reformat_modification_sites.pl';
}


# clean up TMT output
if (defined($autodetect_hash{TMT}) &&
            $autodetect_hash{TMT} ne 'no')
{
    $pipeline_preprocess_str .= sprintf " | \"%s/%s\" -",
        $script_path, 'clean_tmt.pl';
}


$human_flag = 0;
$mouse_flag = 0;
if ($species =~ /Mouse/i) { $mouse_flag = 1; }
if ($species =~ /Human/i) { $human_flag = 1; }


# default to uniprot human annotation
$annotation_species = 'human';
$annotation_source  = 'uniprot';
$annotation_file    = 'ipi_uniprot_annotation_human.txt';

# use refseq instead
if (defined($autodetect_hash{RefSeq}) &&
            $autodetect_hash{RefSeq} eq 'yes')
{
    $annotation_source = 'refseq';
    
    # default to human
    $annotation_file   = 'human_refseq_table.txt';
}

# use mouse annotation
if ($mouse_flag && $human_flag == 0)
{
    $annotation_species = 'mouse';

    # default to uniprot
    $annotation_file = 'ipi_uniprot_annotation_mouse.txt';
    
    if ($annotation_source eq 'refseq')
    {
        $annotation_file = 'mouse_refseq_table.txt';
    }
}

# use human + mouse annotation
if ($mouse_flag && $human_flag)
{
    $annotation_species = 'human_plus_mouse';

    # default to uniprot
    $annotation_file = 'ipi_uniprot_annotation_human_mouse.txt';
    
    if ($annotation_source eq 'refseq')
    {
        $annotation_file = 'human_mouse_refseq_table.txt';
    }
}

printf STDERR "species:     %s\n", $annotation_species;



# run the annotation script
$pipeline_preprocess_str .= sprintf " | \"%s/%s\" \"%s/%s\" -",
        $script_path, 'reannotate_proteomics.pl',
        $script_path, $annotation_file;


# output to cleaned_reannotated.txt
$pipeline_preprocess_str .= sprintf " > \"%s/%s\"",
    $file_path, 'pipeline_cleaned_reannotated.txt';



# deal with intensity/ibaq stuff
$pipeline_ibaq_str = '';
if (defined($autodetect_hash{Rollup}))
{
    # single ibaq file, rename iBAQ columns to Intensity columns
    if ($autodetect_hash{Rollup} eq 'ibaq')
    {
        $pipeline_ibaq_str = sprintf "\"%s/%s\" \"%s/%s\" | sed 's/^iBAQ/Intensity/i' | \"%s/%s\" > \"%s/%s\"",
            $script_path, 'transpose',
            $file_path, 'pipeline_cleaned_reannotated.txt',
            $script_path, 'transpose',
            $file_path, 'pipeline_orig_ibaq.txt';
    }
    # separate intensity and ibaq files
    elsif ($autodetect_hash{Rollup} eq 'intensity_and_ibaq')
    {
        $pipeline_ibaq_str = sprintf "\"%s/%s\" \"%s/%s\" | grep -v -i '^iBAQ' | \"%s/%s\" > \"%s/%s\"",
            $script_path, 'transpose',
            $file_path, 'pipeline_cleaned_reannotated.txt',
            $script_path, 'transpose',
            $file_path, 'pipeline_orig_intensity.txt';

        $pipeline_ibaq_str .= sprintf "; \"%s/%s\" \"%s/%s\" > \"%s/%s\"",
            $script_path, 'replace_maxquant_intensity_with_ibaq.pl',
            $file_path, 'pipeline_cleaned_reannotated.txt',
            $file_path, 'pipeline_orig_ibaq.txt';
    }
    # copy to intensity file unchanged
    else
    {
        $pipeline_ibaq_str .= sprintf "cp -a \"%s/%s\" \"%s/%s\"",
            $file_path, 'pipeline_cleaned_reannotated.txt',
            $file_path, 'pipeline_orig_intensity.txt';
    }
}
# copy to intensity file unchanged
else
{
    $pipeline_ibaq_str .= sprintf "cp -a \"%s/%s\" \"%s/%s\"",
        $file_path, 'pipeline_cleaned_reannotated.txt',
        $file_path, 'pipeline_orig_intensity.txt';
}
# remove the original cleaned file
$pipeline_ibaq_str .= sprintf "; rm \"%s/%s\"",
    $file_path, 'pipeline_cleaned_reannotated.txt';



# normalize the data

# default pipeline
$pipeline_norm_str = sprintf "\"%s/%s\" \"%s/%s\" > \"%s/%s\"",
    $script_path, 'iron_normalize_mass_spec.pl',
    $file_path, 'pipeline_orig_intensity.txt',
    $file_path, 'pipeline_iron_intensity.txt';

if (defined($autodetect_hash{TMT}) &&
            $autodetect_hash{TMT} ne 'no')
{
    if ($autodetect_hash{TMT} eq 'multi' &&
        defined($autodetect_hash{TMT_Channel}) &&
                $autodetect_hash{TMT_Channel} =~ /TMT/)
    {
        $pipeline_norm_str = sprintf "\"%s/%s\" \"%s/%s\" %s > \"%s/%s\"",
            $script_path, 'automate_tmt.pl',
            $file_path, 'pipeline_orig_intensity.txt',
            $autodetect_hash{TMT_Channel},
            $file_path, 'pipeline_iron_intensity.txt';
    }
    # else use default normalization
}


# handle ibaq stuff
if (defined($autodetect_hash{Rollup}))
{
    # single ibaq file, rename iBAQ columns to Intensity columns
    if ($autodetect_hash{Rollup} eq 'ibaq')
    {
        $pipeline_norm_str =~ s/_intensity\.txt/_ibaq\.txt/g;
    }
    # separate intensity and ibaq files
    elsif ($autodetect_hash{Rollup} eq 'intensity_and_ibaq')
    {
        $pipeline_norm_ibaq_str = $pipeline_norm_str;
        $pipeline_norm_ibaq_str =~ s/_intensity\.txt/_ibaq\.txt/g;
        
        $pipeline_norm_str .= '; ' . $pipeline_norm_ibaq_str;
    }
}



printf "%s\n", $pipeline_preprocess_str;
printf "\n";
printf "%s\n", $pipeline_ibaq_str;
printf "\n";
printf "%s\n", $pipeline_norm_str;