Brief overview of running the proteomics and metabolomics pipelines:


  Usage:
    generate_proteomics_glue_script.pl
      input_maxquant_file.txt output_filenames_prefix [[species] autodetect.txt]
       > run_stuff.sh

      # -- or --

    generate_metabolomics_glue_script.pl
      POS_mzmine_output.csv NEG_mzmine_output.csv output_filenames_prefix
       > run_stuff.sh

    # -- followed by: --
    sh run_stuff.sh


    If the input data is lipidomics, split_lipidomics_into_pos_neg.pl must
    be run first, followed by generate_metabolomics_glue_script.pl on the
    resulting separate pos and neg files.

    If the experiment was only assayed in positive or negative ion mode,
    not both, then replace all occurrences of the missing pos/neg file name(s)
    with "NULL", to indicate that no data is available for this injection.


  Species and edited auto-detect parameters file arguments are optional and
  usually not needed.

  See additional documentation for each glue script command further below,
  where the scripts are documented in alphabetical order.


Species will default to auto-detection based on the presence of the species
within the sequence annotation provided by MaxQuant.  If species auto-
detection fails (no or little species annotation present), or the provided
species argument doesn't contain human or mouse, the species will default to
human.  If the wrong species is used, the wrong annotation file will be used,
and the output Target_Species_Flag column will be incorrect as well.  Use
"human_and_mouse", or any string containing both human and mouse within it
somewhere (ie: "mouse-human", "human-bear-pig-mouse") for a mouse human-
xenograph experiment.  This should be auto-detected correctly most of the
time as well, so you generally shouldn't need to manually override the
species auto-detection.

Sometimes, the auto-detection (guessing) of all the various parameters fails
to choose the correct TMT channel for normalization (TMT-126, when the pool
was loaded in the last channel instead).  You'll need to do one of the
following:

 A) re-run generate_proteomics_glue_script.pl with --last-ch or --boost flags
 B) edit (or sed -i) the resulting script to replace the incorrect channel
 C) Copy/paste the auto-detect results into a file and re-generate the script.
    Since I was lazy in my command line parsing, you *MUST* specify the
    species in order to then read in the autodetect.txt (or whatever you named
    the file containing the edited auto-detect output) parameters.

    "generate_proteomics_glue_script.pl input.txt output_prefix species autodetext.txt > run_stuff.sh"

 Generally, --last-ch or --boost is going to be the easiest/best solution.




Documentation for individual programs within the pipeline:



findmedian, iron_generic
  see http://gene.moffitt.org/libaffy/



annotate_metabolomics.pl

  Annotate MZMine (and eventually other software's) metabolomics data with
  external identifiers and updated metabolite names.  Automatically removes
  untrustworthy MZMine peak identifications with m/z error > 10 PPM. Warns of
  unmatched rows, different m/z, and retention times > +/-1 minute from
  reference.

  Usage: annotate_metabolomics.pl identifier_mapping.txt cleaned_mzmine.txt

  Options:
      --ppm N          override default m/z PPM tolerance

  Mapped fields are | (vertical bar) delimited, rather than ; delimited.
  Multiple | delimited subfields have 1:1 correspondence across mapped
  fields.  Original name ("identity") usually exhibits 1:1 correspondence
  with newly mapped name ("Identity Mapped"), but not always, such as when
  multiple original names map to the same metabolite.  There are currently
  no 1:many original:new mappings that have been observed, but many:1 can
  ocasionally occur in older datasets (this is intended behavior).

  Mapping Type field legend:
    1A/B:  conformed name,     original + pos/neg m/z
    2A/B:  conformed synonyms, original + pos/neg m/z
    3A/B:  conformed name,     all rows
    4A/B:  conformed synonyms, all rows
    5A/B:  fuzzy,              original + pos/neg m/z
    9X:    identified row not matched to search database
           (either a bad identification, or metabolite not in database)

    A/B indicate whether there was only a single "reasonable" candidate
    (A) within the m/z range, or multiple "reasonable" candidates (B).



autodetect_ms_data.pl

  Scan MaxQuant output file to determine (or guess) various parameters to feed
  into generate_proteomics_glue_script.pl:

  Usage: autodetect_ms_data.pl [options] maxquant_output.txt [species]

    Options:
      --boost     use highest channel for normalization
      --last-ch   use highest channel for normalization

    Output:
      Modification  [yes/no]         does the data contain modification sites
      RefSeq        [yes/no]         is the data mainly matched against RefSeqs
      TMT           [single/multi/injection]
                                     1 plex, >=2 plexes, >=2 all injection reps
      TMT_Channel   [auto/TMT-????]  reference channel for IRON normalization
      Rollup        [ibaq/intensity/intensity_and_ibaq]
                                     which type of rollup columns to keep
      Species       [human/mouse/human_and_mouse]



automate_tmt.pl

  Normalization of multiple TMT plexes within the same file.

  Usage:
    automate_tmt.pl [options] maxquant_output_file.txt [IRON ref channels] > iron_output.txt

    If IRON ref channel is left unspecified, it defaults to TMT-126C.

    If "variable" is given as the IRON ref channel, things will work a bit
    differently than usual....  The algorithm will assume two replicate
    reference channels and will automatically guess which ones they are by
    finding the two most similar channels within each TMT plex, separately.
    This is a *deprecated* feature, implemented for a single pilot experiment
    a long time ago.

  Options:
    --iron             normalize within-plex prior to other calcuations [default]
    --iron-untilt      account for relative dynamic range
    --no-iron          do not normalize within each plex first
    --debatch          use reference channel for cross-plex normalization [default]
    --no-debatch       do not perform cross-plex normalization
    --leave-ratios     leave cross-plex normalized data as log2 ratios
    --no-leave-ratios  scale cross-plex normalized log2 ratios back into abundances [default]

    --comp-pool        use all-channel geometric mean for cross-plex debatching
    --no-comp-pool     do not create a computational reference pool for cross-plex debatching [default]
    --comp-pool-exclusions-dark
                       auto-excludes dark samples from computational pool
    --comp-pool-exclusions-boost
                       excludes boosting channels (N, N-2) from comp pool
    --comp-pool-exclusions=filename.txt
                       load comp pool sample exclusions from tab-delimited file
    --iron-exclusions=filename.txt
                       exclude row identifiers from IRON training

    --no-iron --no-debatch will leave the output abundances unchanged

  Output:
    Normalized data is output to STDOUT.

    Several other files are generated, extending the original input filename:

      *_scaling_factors.txt  various statistics from the IRON normalizations
      *_findmedian.txt       results from the various findmedian runs
    


clean_tmt.pl

  Reformat TMT data, removing more unwanted columns, rename channel headers

  Usage:
    clean_tmt.pl maxquant_output_file.txt > output.txt

  Output:
    Original file, but with sample columns renamed and various columns
    removed.



csv2tab_not_excel.pl

  Convert CSV text files to tab-delimited text.

  Usage;
    csv2tab_not_excel.pl input.csv > output_tab_delimited.txt

  Output:
    Simple tab-delimited output.  Fields containing commas and/or double-
    quotes are not re-escaped, nor are embedded tabs, newlines, or carriage
    returns supported.  It is slower than python's build-in functions, but is
    more robust to various errors in CSV file formats, and produces a properly
    cleaned (no escaped fields) simple tab-delimited output.

    See comments within script for additonal discussion and rationale.



generate_metabolomics_glue_script.pl

  Generate glue script to launch other scripts involved in the metabolomics
  pipeline.

  Usage: generate_metabolomics_glue_script.pl [options] mzmine_pos.csv mzmine_neg.csv
           [output_prefix [lipidsearch_summary.txt]]

  Options:
      --heavy-spikein            heavy rows are spikeins, leave unscaled (default)
      --heavy-tracer             heavy rows are biological, normalize them
      --no-log2                  disable log2 transform of output data
      --norm-none                disable normalization; use on targeted panels
      --ppm N                    override default m/z PPM tolerance

      --discard-heavy            discard heavy labeled rows
      --discard-unidentified     discard unidentified rows

    options which use the MZmine "row number of detected peaks" column:
      --discard-single-pregap    discard pre gap-filled single-hit rows (default)
      --keep-single-pregap       keep pre gap-filled single-hit rows

    Specifying a LipidSearch summary file will use the summary file to add an
    additional pipeline output column to flag rows identified as main ions by
    LipidSearch.  This is in addition to the BBSR method of lipidomics main ion
    assignment, which uses summed normalized log2 abundances instead of unlogged
    unnormalized abundances.

    MZmine output files for positive and negative ion mode acquisitions must
    contain only pos or neg mode data.  The file names themselves must also
    contain pos or neg, as appropriate, within their names, otherwise the
    program will abort due to failed sanity checking.  This overly paranoid
    sanity checking is there due to various errors in files I have received in
    the past, errors which this will catch.

    If sample names contain non-alphanumeric- separated pos/neg in their names,
    the pos/neg will be removed from the sample names in order to facilitate
    merging pos/neg data together.  If pos/neg sample names conflict with the
    file mode specified, or conflicting pos/neg sample names occur within a
    file, pos/neg will not be stripped from any sample names within that file
    prior to merging.  No warning messages are currently issued when most such
    conflicts are detected (I should probably add some).

    If the experiment was only assayed in positive or negative ion mode,
    not both, then replace all occurrences of the missing pos/neg file name(s)
    with "NULL", to indicate that no data is available for this injection.


  Output:
    .sh file to glue the various scripts together.

    Running the script will generate IRON normalized, merged pos/neg
    metabolites, with new columns added to denote spike-ins, identified, and
    non-spikein identified rows.



generate_proteomics_glue_script.pl

  Generate glue script to launch other scripts involved in the proteomics
  pipeline.  Makes a few assumptions if no autodetect.txt file is specified.
  Assumes human if no species is given.

  Usage: generate_proteomics_glue_script.pl [options] maxquant_output.txt output_root_name [[species] autodetect.txt] > run_proteomics.sh

  Options:
    --boost         use highest channel for normalization
    --last-ch       use highest channel for normalization

    --iron          apply IRON normalization (default)
    --iron-auto-ch  auto-pick different reference channel per-plex
                    *** forces --no-debatch ***
    --iron-untilt   account for relative dynamic range
    --no-iron       disable IRON normalization

    --debatch       cross-plex de-batch TMT data (default)
    --no-debatch    disable cross-plex de-batching

    --comp-pool     average all plex channels for cross-plex de-batching
    --no-comp-pool  use reference channel for de-batching (default)
    --comp-pool-exclusions-dark
                       auto-excludes dark samples from computational pool
    --comp-pool-exclusions-boost
                       excludes boosting channels (N, N-2) from comp pool
    --comp-pool-exclusions=filename.txt
                       load comp pool sample exclusions from tab-delimited file

  Output:
    .sh file to glue the various scripts together.

    Running the script will generate IRON normalized, proteins or peptides or
    modification sites, depending on the input file, with various columns
    removed, reformatted, or reannotated depending on the input data type.

  

iron_normalize_mass_spec.pl

  General purpose mass spec normalization, all samples vs. single reference.
  Can be used on single-plex TMT data as well.  It will not give any warnings
  if run on TMT data containing multiple plexes, but the output will likely
  not be normalized as desired (no between-plex "batch correction" will have
  been performed).  Will also work on metabolomics data.

  Usage:
    iron_normalize_mass_spec.pl [options] input_data.pl > normalized_output.txt

  Options:
    --iron-untilt                 account for relative dynamic range
    --iron-exclusions=filename    identifiers to exclude from training
    --iron-spikeins=filename      list of spikein identifiers

    --log2                        output log2 abundances [default]
    --no-strip-sample-names       keep original full sample headers
    --no-log2                     output unlogged abundances
    --norm-none                   disable normalization
    --unlog2                      exponentiate input data, pow(2, value)

  Output:
    Assuming the input file format is recognized and it is therefore able
    to correctly identify which columns are data, then the input data columns
    will be replaced with normalized values, stripping various header text
    from the sample names unless --no-strip-sample-names is specified.



lipidomics_assign_main_ion.pl

  Assign lipid isomer groups based on retention times, flag the most abundant
  ion within each lipid isomer gorup as the main ion.  If an optional
  LipidSearch summary file is provided, an additional column of LipidSearch
  main ion assignments will be added.

  Usage:
    lipidomics_assign_main_ion.pl iron_log2_merged.txt [lipidsearch_summary.txt]



merge_cleaned_modification_site_files.pl

  Merge multiple different cleaned modification site proteomics files together

  Usage:
    merge_cleaned_modification_site_files.pl clean1.txt clean2.txt [...] > output.txt

    Attempt to merge multiple different cleaned modification site files
    together and generate various consensus columns (such as annotation, PEP
    score, etc.) where possible.

    This tool needs a good bit more work to support other types of mass spec
    data.  This is probably not so much for general use at the moment.

  Output:
    Various columns are duplicated/conformed/renamed to merge the various
    files together as best as possible.



merge_metabolomics_pos_neg.pl

  Merge positive and negative ion mode tab-delimited metabolomics data files
  together

  Usage:
    merge_metabolomics_pos_neg.pl mzmine_pos_file.txt mzmine_neg_file.txt [output_sample_table.txt]

    !!! NOTE -- files must be tab-delimited, not the original csv files !!!

  Output:
    Merged positive and negative data files, auto-conforming the sample
    identifiers and stripping their pos/neg strings in order to merge
    them together.  The auto-conforming process is reasonbly robust to
    differences in spacing and/or punctiation surrounding the pos/neg
    substrings, however other differences in spacing or punctuation elsewhere
    in the sample names will not be auto-conformed, and will result in
    separate pos/neg columns in the merged output.  It is up to the user to
    make sure that sample naming conventions are the same between the pos and
    neg sample names.

    A sample name mapping table is exported to a default file name of
    "pipeline_metabolomics_sample_table.txt".  An optional 3rd agument
    overrides the default output sample table file name with the string
    given in the argument.

    If the experiment was only assayed in positive or negative ion mode,
    not both, then replace all occurrences of the missing pos/neg file name(s)
    with "NULL", to indicate that no data is available for this injection.



reannotate_proteomics.pl

  Add annotation columns based on various proteomics accession columns

  Usage:
    reannotate_proteomics.pl accession_annotation.txt data_file.txt > reannotated.txt

  Output:
    New columns will be added containing the additional annotation, as well
    as flags for all-reverse, contaminant, and cow contaminants (if any of
    the accession fields contain the text "bos taurus").



reformat_modification_sites.pl

  Reposition incorrectly reported modification site positions in the
  Maxquant modified sequence string, as well as conform how they are
  represented.

  Usage:
    reformat_modification_sites.pl maxquant_output.txt > cleaned.txt

    Supports all of the common modifications that we assay; will need to be
    edited to support new ones that have not been requested yet:

      pS/T/Y
      Desthiobiotin / ATP
      Acetyl (K)
      GlyGly (K)
      Oxidation (M)

  Output:
    Modification sequences are re-generated with the modification in the
    corrected position and conformed to an older style of site format.
    Various new columns, including a ModificationID column, are also generated,
    which can be used to better match modification sites between files.



replace_maxquant_intensity_with_ibaq.pl

  Move iBAQ rollup columns into Intensity columns, so other scripts will work

  Usage:
    replace_maxquant_intensity_with_ibaq.pl input.txt > output.txt

  Output:
    replaces the Intensity columns with iBAQ columns, deleting the iBAQ columns



scan_oct.pl

  Scan positive ion MZmine2 metabolomics output for OCT contamination

  Usage:

    scan_oct.pl mzmine_pos.csv [output_oct_rows.txt [output_oct_metrics.txt]]

    Both CSV and tab-delimited input is accepted.
    (+) ion mode files yield best results, since OCT doesn't show up in (-).
    Any (-) ion mode rows detected as OCT are likely false positives.

    If output file names are not specified, default file names will be used:
      oct_detected_rows.txt oct_sample_metrics.txt

    Output files contain rows identified as potential OCT, and OCT abundance-
    related metrics, respectively.



split_lipidomics_into_pos_neg.pl

    Split single tab-delimited lipidomics file into separate pos/neg files.

    Usage: split_lipidomics_into_pos_neg.pl [options] tab_delimited.txt
             [outfile_pos outfile_neg]

      Options:
        --all-adducts     keep all adducts, including unobservable

    Output:
      outfile_pos and outfile_neg default to lipidomics_split_pos.txt and
      lipidomics_split_neg.txt if no output filenames are specified.



strip_maxquant_columns.pl

  Remove columns from MaxQuant data that we don't use

  Usage:
    strip_maxquant_columns.pl maxquant_file.txt > output.txt

  Output:
    various uninformative or problematic (too wide) columns are removed



strip_metabolomics_columns.pl

  Remove MZMine columns that we don't use

  Usage: strip_mzmine_columns.pl [options] mzmine_tab_delimited.txt [unidentified.txt spikeins.txt]

  Options:
      --heavy-spikein            treat heavy rows as spike-ins
      --heavy-tracer             treat heavy rows as biological
      --ppm N                    override default m/z PPM tolerance

      --discard-heavy            discard heavy labeled rows
      --discard-unidentified     discard unidentified rows

    options which use the MZmine "row number of detected peaks" column:
      --discard-single-pregap    discard pre gap-filled single-hit rows (default)
      --keep-single-pregap       keep pre gap-filled single-hit rows

  Output:
    Various uninformative columns are removed.
    "one-hit wonder" rows will also be removed, unless --keep-single-pregap
     is also specified.



transpose

  Transpose row/col tab-delimited text file to col/row.

  Usage:
    (to compile executable from source code: gcc -O3 transpose.c -o transpose)

    transpose input.txt > output.txt
      or
    transpose input.txt output.txt

  Output:
    Rows/columns are transposed.  Missing data at the ends of rows and columns
    is preserved, if possible.

    If a second argument is given as the output file, the program will allow
    for the output file to overwrite the input file, since the input file is
    fully buffered into memory and closed before the output file is written.
    The program does not check for this or prevent it.  USE THE FOLLOWING
    SYNTAX WITH CAUTION:

      transpose input_to_be_overwritten.txt input_to_be_overwritten.txt



Annotation files:

The following text files are soft links to files on the Moffitt unix cluster,
since Git cannot handle versioning large files efficiently.  I usually update
these files once per year, around the first of the new year.

human_mouse_refseq_table.txt              Human+Mouse NCBI accessions
human_refseq_table.txt                    Human NCBI accessions
ipi_uniprot_annotation_human.txt          Human IPI & Uniprot accessions
ipi_uniprot_annotation_human_mouse.txt    Human+Mouse IPI & Uniprot accessions
ipi_uniprot_annotation_mouse.txt          Mouse IPI & Uniprot accessions
mouse_refseq_table.txt                    Mouse NCBI accessions



Scripts for re-grouping maxquant peptides:

  The following are a work in progress, and not yet supported for general use.

  regroup_maxquant_peptides.pl
    create new gene-level peptide groupings

  filter_regrouped_maxquant_peptides.pl  
    filter new peptide groupings to remove less useful groupings

  summarize_regrouped_maxquant_peptides.pl
    summarize/annotate new peptide groupings

  rollup_normalized_tmt.pl
    use new peptide groupings to re-rollup peptides into gene-level values
