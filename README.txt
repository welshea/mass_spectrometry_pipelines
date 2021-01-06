Documentation for individual programs within the pipeline:



findmedian, iron_generic
  see http://gene.moffitt.org/libaffy/




autodetect_ms_data.pl

  Scan MaxQuant output file to determine (or guess) various parameters to feed
  into generate_proteomics_glue_script.pl:

  Usage:
    autodetect_ms_data.pl maxquant_output_file.txt > autodetect.txt

  Output:
    Modification  [yes/no]         does the data contain modification sites
    RefSeq        [yes/no]         is the data mainly matched against RefSeqs
    TMT           [single/multi/injection]
                                   1 plex, >=2 plexes, >=2 all injection reps
    TMT_Channel   [auto/TMT-????]  reference channel for IRON normalization
    Rollup        [ibaq/intensity/intensity_and_ibaq]
                                   which type of rollup columns to keep



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
    --no-iron          do not normalize within each plex first
    --debatch          use reference channel for cross-plex normalization [default]
    --no-debatch       do not perform cross-plex normalization
    --leave-ratios     leave cross-plex normalized data as log2 ratios
    --no-leave-ratios  scale cross-plex normalized log2 ratios back into abundances [default]
    --comp-pool        use all-channel geometric mean for cross-plex debatching
    --no-comp-pool     do not create a computational reference pool for cross-plex debatching [default]
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

  Syntax:
    clean_tmt.pl maxquant_output_file.txt > output.txt

  Output:
    Original file, but with sample columns renamed and various columns
    removed.




csv2tab_not_excel.pl

  Convert CSV text files to tab-delimited text.

  Syntax;
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

  Usage: generate_metabolomics_glue_script.pl [options] mzmine_pos.csv mzmine_neg.csv output_prefix

  Options:
      --discard-heavy            discard heavy labeled rows
      --discard-unidentified     discard unidentified rows

    options which use the MZmine "row number of detected peaks" column:
      --discard-single-pregap    discard pre gap-filled single-hit rows (default)
      --keep-single-pregap       keep pre gap-filled single-hit rows

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

  Output:
    .sh file to glue the various scripts together.

    Running the script will generate IRON normalized, merged pos/neg
    metabolites, with new columns added to denote spike-ins, identified, and
    non-spikein identified rows.




generate_proteomics_glue_script.pl

  Generate glue script to launch other scripts involved in the proteomics
  pipeline.  Makes a few assumptions if no autodetect.txt file is specified.
  Assumes human if no species is given.

  Usage:
    generate_proteomics_glue_script.pl maxquant_output.txt output_root_name [[human/mouse/human_plus_mouse] autodetect.txt] > run_proteomics.sh

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
    --iron-exclusions=exclusions.txt
      list of row identifiers to exclude from findmedian and IRON training
    --iron-spieins=spikeins.txt
      list of row identifiers to exclude from calculations and normalization;
      similar to --iron-exlusions but also leaves them un-normalized
    --no-strip-sample-names
      do not strip various strings from sample names, such as:
        Intensity, Peak height, Peak area, .mzXML, Total Area
    --log2   (default)
      output log2 abundances
    --no-log2
      output unlogged abundances
    --unlog2
      un-log2 the input data prior to normalization

  Output:
    Assuming the input file format is recognized and it is therefore able
    to correctly identify which columns are data, then the input data columns
    will be replaced with normalized values, stripping various header text
    from the sample names unless --no-strip-sample-names is specified.




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
    merge_metabolomics_pos_neg.pl mzmine_pos_file.txt mzmine_neg_file.txt

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
