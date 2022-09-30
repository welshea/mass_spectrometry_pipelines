#!/usr/bin/perl -w

# 2022-06-11:  add --iron-auto-ch to auto-pick each ref channel per-plex
# 2022-04-06:  add --iron-untilt to address rare dynamic range compression
# 2022-02-03:  bugfix: honor auto-detected species
# 2022-01-24:  support passing new comp pool exclusions flags
# 2022-01-19:  auto-exclude dark samples from injection replicates comp pool
# 2022-01-19:  add some support for auto reference channels
# 2021-10-28:  support --(no-)iron --(no-)debatch --(no-)comp-pool
# 2021-09-08:  print usage statement to STDERR instead of STDOUT
# 2021-07-23:  add _log2 to iron output filenames to indicate it is log2
# 2021-05-27:  change default output_root_name
# 2021-05-21:  add new autodetect_ms_data.pl arguments
# 2020-12-21:  add missing # in front of --comp-pool changelog line below...
# 2020-12-18:  change behavior/syntax to that of metabolomics glue script
#              use --comp-pool for TMT 100% injection replicates
# 2020-09-10:  use different annotation files
# 2020-02-28:  clean phosphosites before cleaning TMT
#               this allows clean TMT to detect phospho data and use
#               the ModificationID column as the RowIdentifier

use File::Spec;
use File::Basename;


$input_filename    = '';
$num_non_options   = 0;
$boost_flag        = 0;
$no_debatch_flag   = 0;
$no_iron_flag      = 0;
$comp_pool_flag                 = 0;    # use comp pool instead of real pool
$comp_pool_exclude_flag         = 0;
$comp_pool_exclude_filename     = '';
$comp_pool_exclude_dark_flag    = 0;
$comp_pool_exclude_boost_flag   = 0;
$iron_untilt_flag               = 0;
$auto_single_variable_pool_flag = 0;    # different auto reference per plex


$syntax_error_flag = 0;
for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        if ($field =~ /^--boost$/ ||
            $field =~ /^--last-ch$/)
        {
            $boost_flag = 1;
        }

        elsif ($field =~ /^--no-debatch$/)
        {
            $no_debatch_flag = 1;
        }
        elsif ($field =~ /^--debatch$/)
        {
            $no_debatch_flag = 0;
        }

        elsif ($field =~ /^--no-iron$/)
        {
            $no_iron_flag = 1;
        }
        elsif ($field =~ /^--iron$/)
        {
            $no_iron_flag = 0;
        }

        elsif ($field =~ /^--comp-pool$/)
        {
            $comp_pool_flag = 1;
        }
        elsif ($field =~ /^--no-comp-pool$/)
        {
            $comp_pool_flag = 0;
        }
        elsif ($field =~ /^--comp-pool-exclusions\=(\S+)/)
        {
            # $comp_pool_flag = 1;
            $comp_pool_exclude_flag = 1;
            $comp_pool_exclude_filename = $1;
        }
        elsif ($field =~ /^--comp-pool-exclusions-dark$/)
        {
            # $comp_pool_flag = 1;
            $comp_pool_exclude_flag = 1;
            $comp_pool_exclude_filename = '';
            $comp_pool_exclude_dark_flag = 1;
        }
        elsif ($field =~ /^--comp-pool-exclusions-boost$/)
        {
            # $comp_pool_flag = 1;
            $comp_pool_exclude_flag = 1;
            $comp_pool_exclude_boost_flag = 1;
            $comp_pool_exclude_filename = '';
        }
        elsif ($field =~ /^--iron-untilt$/)
        {
            $iron_untilt_flag = 1;
            $no_iron_flag     = 0;
        }
        elsif ($field =~ /^--iron-auto-ch$/)
        {
            $auto_single_variable_pool_flag = 1;
            $no_debatch_flag = 1;
        }
        else
        {
            printf "ABORT -- unknown option %s\n", $field;
            $syntax_error_flag = 1;
        }
    }
    else
    {
        if ($num_non_options == 0)
        {
            $input_filename = $field;
            $num_non_options++;
        }
        elsif ($num_non_options == 1)
        {
            $output_root_name = $field;
            $num_non_options++;
        }
        # must be 'human', 'mouse', or 'human_and_mouse'
        elsif ($num_non_options == 2)
        {
            $species = $field;
            $num_non_options++;
        }
        elsif ($num_non_options == 3)
        {
            $autodetect_filename = $field;
            $num_non_options++;
        }
    }
}


# force --no-debatch if different ref channels are used per-plex
if ($auto_single_variable_pool_flag)
{
    $no_debatch_flag = 1;

    printf STDERR "cross-plex de-batching disabled due to --iron-auto-ch\n";
}



if ($input_filename eq '')
{
    $syntax_error_flag = 1;
}

if ($syntax_error_flag)
{
    print STDERR "Usage: generate_proteomics_glue_script.pl [options] maxquant_output.txt output_root_name [[species] autodetect.txt] > run_proteomics.sh\n";
    print STDERR "\n";
    print STDERR "  Options:\n";
    print STDERR "    --boost         use highest channel for normalization\n";
    print STDERR "    --last-ch       use highest channel for normalization\n";
    print STDERR "\n";
    print STDERR "    --iron          apply IRON normalization (default)\n";
    print STDERR "    --iron-auto-ch  auto-pick different reference channel per-plex\n";
    print STDERR "                    *** forces --no-debatch ***\n";
    print STDERR "    --iron-untilt   account for relative dynamic range\n";
    print STDERR "    --no-iron       disable IRON normalization\n";
    print STDERR "\n";
    print STDERR "    --debatch       cross-plex de-batch TMT data (default)\n";
    print STDERR "    --no-debatch    disable cross-plex de-batching\n";
    print STDERR "\n";
    print STDERR "    --comp-pool     average all plex channels for cross-plex de-batching\n";
    print STDERR "    --no-comp-pool  use reference channel for de-batching (default)\n";
    print STDERR "    --comp-pool-exclusions-dark\n";
    print STDERR "                       auto-excludes dark samples from computational pool\n";
    print STDERR "    --comp-pool-exclusions-boost\n";
    print STDERR "                       excludes boosting channels (N, N-2) from comp pool\n";
    print STDERR "    --comp-pool-exclusions=filename.txt\n";
    print STDERR "                       load comp pool sample exclusions from tab-delimited file\n";
    
    exit(1);
}


# optional
if (!defined($output_root_name))
{
    $output_root_name = 'proteomics_pipeline';
}
if (!defined($species))
{
    $species = '';
}

# this correctly doesn't follow symlinked scripts, so you get the right path
# various other methods return the destination linked path instead,
#  which we don't want here
$script_path   = dirname(File::Spec->rel2abs(__FILE__));
#$file_path     = dirname(File::Spec->rel2abs($input_filename));
#$filename_base = basename($input_filename);

#$pid = $$;
#$tmp_path = sprintf "/tmp/proteomics_pipeline_%s/", $pid;

#`mkdir $tmp_path`;

#printf STDERR "script_path: %s\n", $script_path;
#printf STDERR "file_path:   %s\n", $file_path;
#printf STDERR "basename:    %s\n", $filename_base;

if (!defined($autodetect_filename))
{
    if ($boost_flag)
    {
        $cmd_str = sprintf "%s --boost \"%s\" \"%s\"",
            'autodetect_ms_data.pl',
            $input_filename,
            $species;
    }
    else
    {
        $cmd_str = sprintf "%s \"%s\" \"%s\"",
            'autodetect_ms_data.pl',
            $input_filename,
            $species;
    }
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


# set TMT channel to "auto" if --iron-auto-ch was enabled
if ($auto_single_variable_pool_flag)
{
    $autodetect_hash{TMT_Channel} = "auto";
}

@autodetect_key_array = sort keys %autodetect_hash;

foreach $key (@autodetect_key_array)
{
    $value = $autodetect_hash{$key};
    
    printf STDERR "%s\t%s\n", $key, $value;
}


# initial command which will pipe the data into additional future commands
$pipeline_preprocess_str = sprintf "%s \"%s\"",
    'strip_maxquant_columns.pl',
    $input_filename;


if (defined($autodetect_hash{Modification}) &&
            $autodetect_hash{Modification} eq 'yes')
{
    $pipeline_preprocess_str .= sprintf " \\\n  | %s -",
        'reformat_modification_sites.pl';
}


# clean up TMT output
if (defined($autodetect_hash{TMT}) &&
            $autodetect_hash{TMT} ne 'no')
{
    $pipeline_preprocess_str .= sprintf " \\\n  | %s -",
        'clean_tmt.pl';
}


$human_flag = 0;
$mouse_flag = 0;
if (defined($autodetect_hash{Species}))
{
    $species = $autodetect_hash{Species};
}
if ($species =~ /Mouse/i) { $mouse_flag = 1; }
if ($species =~ /Human/i) { $human_flag = 1; }


# default to uniprot human annotation
$annotation_species = 'human';
$annotation_source  = 'uniprot';

$annotation_file  = 'merged_protein_annotations_human.txt';

# use refseq instead
if (defined($autodetect_hash{RefSeq}) &&
            $autodetect_hash{RefSeq} eq 'yes')
{
    $annotation_source = 'refseq';
    
    # default to human
    $annotation_file   = 'merged_protein_annotations_human.txt';
    
    # override with combined IPI + Uniprot + Refseq + some Ensembl
    $annotation_file   = 'merged_protein_annotations_human.txt';
}

# use mouse annotation
if ($mouse_flag && $human_flag == 0)
{
    $annotation_species = 'mouse';

    # default to uniprot
    $annotation_file = 'merged_protein_annotations_mouse.txt';
    
    if ($annotation_source eq 'refseq')
    {
        $annotation_file = 'merged_protein_annotations_mouse.txt';
    }

    # override with combined IPI + Uniprot + Refseq + some Ensembl
    $annotation_file = 'merged_protein_annotations_mouse.txt';
}

# use human + mouse annotation
if ($mouse_flag && $human_flag)
{
    $annotation_species = 'human_plus_mouse';

    # default to uniprot
    $annotation_file = 'merged_protein_annotations_human_mouse.txt';
    
    if ($annotation_source eq 'refseq')
    {
        $annotation_file = 'merged_protein_annotations_human_mouse.txt';
    }

    # override with combined IPI + Uniprot + Refseq + some Ensembl
    $annotation_file = 'merged_protein_annotations_human_mouse.txt';
}

#printf STDERR "species:     %s\n", $annotation_species;



# run the annotation script
$pipeline_preprocess_str .= sprintf " \\\n  | %s \"%s/%s\" -",
        'reannotate_proteomics.pl',
        $script_path, $annotation_file;


# output to cleaned_reannotated.txt
$pipeline_preprocess_str .= sprintf " \\\n  > %s%s",
    $output_root_name, '_cleaned_reannotated.txt';



# deal with intensity/ibaq stuff
$pipeline_ibaq_str = '';
if (defined($autodetect_hash{Rollup}))
{
    # single ibaq file, rename iBAQ columns to Intensity columns
    if ($autodetect_hash{Rollup} eq 'ibaq')
    {
        $pipeline_ibaq_str = sprintf "%s \"%s%s\" \\\n  | sed 's/^iBAQ/Intensity/i' \\\n  | %s \\\n  > \"%s%s\"",
            'transpose',
            $output_root_name, '_cleaned_reannotated.txt',
            'transpose',
            $output_root_name, '_orig_ibaq.txt';
    }
    # separate intensity and ibaq files
    elsif ($autodetect_hash{Rollup} eq 'intensity_and_ibaq')
    {
        $pipeline_ibaq_str = sprintf "%s \"%s%s\" \\\n  | grep -v -i '^iBAQ' \\\n  | %s \\\n  > \"%s%s\"",
            'transpose',
            $output_root_name, '_cleaned_reannotated.txt',
            'transpose',
            $output_root_name, '_orig_intensity.txt';

        $pipeline_ibaq_str .= sprintf "\n%s \"%s%s\" \\\n  > \"%s%s\"",
            'replace_maxquant_intensity_with_ibaq.pl',
            $output_root_name, '_cleaned_reannotated.txt',
            $output_root_name, '_orig_ibaq.txt';
    }
    # copy to intensity file unchanged
    else
    {
        $pipeline_ibaq_str .= sprintf "cp -a \"%s%s\" \"%s%s\"",
            $output_root_name, '_cleaned_reannotated.txt',
            $output_root_name, '_orig_intensity.txt';
    }
}
# copy to intensity file unchanged
else
{
    $pipeline_ibaq_str .= sprintf "cp -a \"%s%s\" \"%s%s\"",
        $output_root_name, '_cleaned_reannotated.txt',
        $output_root_name, '_orig_intensity.txt';
}
# remove the original cleaned file
$pipeline_ibaq_str .= sprintf "\nrm \"%s%s\"",
    $output_root_name, '_cleaned_reannotated.txt';



# normalize the data

# default pipeline
$options_str = '';
if ($no_iron_flag)
{
    $options_str = '--norm-none';
}
if ($iron_untilt_flag)
{
    $options_str = '--iron-untilt';
}
$pipeline_norm_str = sprintf "%s %s \"%s%s\" \\\n  > \"%s%s\"",
    'iron_normalize_mass_spec.pl',
    $options_str,
    $output_root_name, '_orig_intensity.txt',
    $output_root_name, '_iron_log2_intensity.txt';

if (defined($autodetect_hash{TMT}) &&
            $autodetect_hash{TMT} ne 'no')
{
    $options_str = '';
    if ($no_iron_flag)
    {
        $options_str .= ' --no-iron';
    }
    elsif ($iron_untilt_flag)
    {
        $options_str .= ' --iron-untilt';
    }
    
    if ($no_debatch_flag)
    {
        $options_str .= ' --no-debatch';
    }
    if ($comp_pool_flag)
    {
        $options_str .= ' --comp-pool';
    }
    if ($comp_pool_flag && $comp_pool_exclude_flag &&
        $comp_pool_exclude_filename ne '')
    {
        $options_str .= sprintf " --comp-pool-exclusions=\"%s\"",
            $comp_pool_exclude_filename;
    }
    if ($comp_pool_flag && $comp_pool_exclude_flag &&
        $comp_pool_exclude_dark_flag)
    {
        $options_str .= ' --comp-pool-exclusions-dark';
    }
    if ($comp_pool_flag && $comp_pool_exclude_flag &&
        $comp_pool_exclude_boost_flag)
    {
        $options_str .= ' --comp-pool-exclusions-boost';
    }

    $options_str =~ s/^\s+//;

    # multiple plexes
    if ($autodetect_hash{TMT} eq 'multi' &&
        defined($autodetect_hash{TMT_Channel}) &&
                ($autodetect_hash{TMT_Channel} =~ /TMT/ ||
                 $autodetect_hash{TMT_Channel} eq 'auto'))
    {
        $pipeline_norm_str = sprintf "%s %s \"%s%s\" %s \\\n  > \"%s%s\"",
            'automate_tmt.pl',
            $options_str,
            $output_root_name, '_orig_intensity.txt',
            $autodetect_hash{TMT_Channel},
            $output_root_name, '_iron_log2_intensity.txt';
    }
    # multiple plexes, 100% injection replicates
    elsif ($autodetect_hash{TMT} eq 'injection' &&
           defined($autodetect_hash{TMT_Channel}) &&
           ($autodetect_hash{TMT_Channel} =~ /TMT/ ||
            $autodetect_hash{TMT_Channel} eq 'auto'))
    {
        $pipeline_norm_str = sprintf "%s %s --comp-pool --comp-pool-exclusions-dark \"%s%s\" %s \\\n  > \"%s%s\"",
            'automate_tmt.pl',
            $options_str,
            $output_root_name, '_orig_intensity.txt',
            $autodetect_hash{TMT_Channel},
            $output_root_name, '_iron_log2_intensity.txt';

        # remove 2nd, redundant, --comp-pool flag
        if ($comp_pool_flag)
        {
            $pipeline_norm_str =~ s/ --comp-pool\b//;
        }
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
        
        $pipeline_norm_str .= "\n" . $pipeline_norm_ibaq_str;
    }
}



printf "%s\n", $pipeline_preprocess_str;
printf "\n";
printf "%s\n", $pipeline_ibaq_str;
printf "\n";
printf "%s\n", $pipeline_norm_str;
