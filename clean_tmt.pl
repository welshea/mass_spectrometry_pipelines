#!/usr/bin/perl -w

# 2022-05-06: fixed sample channel detection bug introduced on 2022-04-21
# 2022-04-21: make sure split doesn't remove empty trailing fields
# 2022-04-21: begin adding support for Proteome Discoverer
# 2022-04-21: remove all double quotes, not just the first one per line (oops)
# 2021-11-02: preserve and reformat run# in single-plex injection replicates
# 2020-08-03: finish support for TMT16, improve unsupported #channels error
# 2020-07-09: add support for TMT16, add letters to the ends of ALL channels
# 2020-02-28: replace | with ~ in RowIdentifier

sub bless_delimiter_space
{
    my $text = $_[0];

    $text =~ s/\//\|/g;
    $text =~ s/\;/\|/g;
    $text =~ s/,/\|/g;
    $text =~ s/\s+/\|/g;
    $text =~ s/\|+/\|/g;
    $text =~ s/^\|//;
    $text =~ s/\|$//;
    $text =~ s/\|/ /g;
    
    return $text;
}

# to be used with sort function
sub compare_accession
{
    my $value_a = $a;
    my $value_b = $b;

    # not REV__
    if (!($value_a =~ /REV__/) &&   $value_b =~ /REV__/)  {return -1;}
    if (  $value_a =~ /REV__/  && !($value_b =~ /REV__/)) {return 1;}

    # not CON__
    if (!($value_a =~ /CON__/) &&   $value_b =~ /CON__/)  {return -1;}
    if (  $value_a =~ /CON__/  && !($value_b =~ /CON__/)) {return 1;}

    # not H-INV
    if (!($value_a =~ /\b_*HIT/) &&   $value_b =~ /\b_*HIT/)  {return -1;}
    if (  $value_a =~ /\b_*HIT/  && !($value_b =~ /\b_*HIT/)) {return 1;}

    # not ENSEMBL
    if (!($value_a =~ /\b_*ENS/) &&   $value_b =~ /\b_*ENS/)  {return -1;}
    if (  $value_a =~ /\b_*ENS/  && !($value_b =~ /\b_*ENS/)) {return 1;}
    
    # best refseqs
    if (  $value_a =~ /\bNM_/  && !($value_b =~ /\bNM_/)) {return -1;}
    if (!($value_a =~ /\bNM_/) &&   $value_b =~ /\bNM_/)  {return 1;}
    if (  $value_a =~ /\bNR_/  && !($value_b =~ /\bNR_/)) {return -1;}
    if (!($value_a =~ /\bNR_/) &&   $value_b =~ /\bNR_/)  {return 1;}
    if (  $value_a =~ /\bNP_/  && !($value_b =~ /\bNP_/)) {return -1;}
    if (!($value_a =~ /\bNP_/) &&   $value_b =~ /\bNP_/)  {return 1;}
    
    # worse refseqs
    if (  $value_a =~ /\bXM_/  && !($value_b =~ /\bXM_/)) {return -1;}
    if (!($value_a =~ /\bXM_/) &&   $value_b =~ /\bXM_/)  {return 1;}
    if (  $value_a =~ /\bXR_/  && !($value_b =~ /\bXR_/)) {return -1;}
    if (!($value_a =~ /\bXR_/) &&   $value_b =~ /\bXR_/)  {return 1;}
    if (  $value_a =~ /\bXP_/  && !($value_b =~ /\bXP_/)) {return -1;}
    if (!($value_a =~ /\bXP_/) &&   $value_b =~ /\bXP_/)  {return 1;}

    # worse refseqs
    if (  $value_a =~ /\bYM_/  && !($value_b =~ /\bYM_/)) {return -1;}
    if (!($value_a =~ /\bYM_/) &&   $value_b =~ /\bYM_/)  {return 1;}
    if (  $value_a =~ /\bYR_/  && !($value_b =~ /\bYR_/)) {return -1;}
    if (!($value_a =~ /\bYR_/) &&   $value_b =~ /\bYR_/)  {return 1;}
    if (  $value_a =~ /\bYP_/  && !($value_b =~ /\bYP_/)) {return -1;}
    if (!($value_a =~ /\bYP_/) &&   $value_b =~ /\bYP_/)  {return 1;}

    # accession, alphabetical
    return ($value_a cmp $value_b);
}



# rename MaxQuant sample column headers to conform to the TMT pipeline
$channel_map_table[6][0]   = '126C';	# 126C
$channel_map_table[6][1]   = '127C';
$channel_map_table[6][2]   = '128C';
$channel_map_table[6][3]   = '129C';
$channel_map_table[6][4]   = '130C';
$channel_map_table[6][5]   = '131C';

$channel_map_table[10][0]  = '126C';	# 126C
$channel_map_table[10][1]  = '127N';
$channel_map_table[10][2]  = '127C';
$channel_map_table[10][3]  = '128N';
$channel_map_table[10][4]  = '128C';
$channel_map_table[10][5]  = '129N';
$channel_map_table[10][6]  = '129C';
$channel_map_table[10][7]  = '130N';
$channel_map_table[10][8]  = '130C';
$channel_map_table[10][9]  = '131N';	# 131N

$channel_map_table[11][0]  = '126C';	# 126C
$channel_map_table[11][1]  = '127N';
$channel_map_table[11][2]  = '127C';
$channel_map_table[11][3]  = '128N';
$channel_map_table[11][4]  = '128C';
$channel_map_table[11][5]  = '129N';
$channel_map_table[11][6]  = '129C';
$channel_map_table[11][7]  = '130N';
$channel_map_table[11][8]  = '130C';
$channel_map_table[11][9]  = '131N';
$channel_map_table[11][10] = '131C';

$channel_map_table[16][0]  = '126C';	# 126C
$channel_map_table[16][1]  = '127N';
$channel_map_table[16][2]  = '127C';
$channel_map_table[16][3]  = '128N';
$channel_map_table[16][4]  = '128C';
$channel_map_table[16][5]  = '129N';
$channel_map_table[16][6]  = '129C';
$channel_map_table[16][7]  = '130N';
$channel_map_table[16][8]  = '130C';
$channel_map_table[16][9]  = '131N';
$channel_map_table[16][10] = '131C';
$channel_map_table[16][11] = '132N';
$channel_map_table[16][12] = '132C';
$channel_map_table[16][13] = '133N';
$channel_map_table[16][14] = '133C';
$channel_map_table[16][15] = '134N';	# 134N


$pep_cutoff  = 0.05;
#$pep_cutoff_site = 0.20;	# 0.20 resulted in bad PCA
#score_cutoff = 75;

#$pep_cutoff   =  9999;
#$pep_cutoff   =  0.2;		# discard if > $pep_cutoff
#$score_cutoff = -9999;
$mass_err_cutoff = 5;		# discard if > $mass_err_cutoff


$infile = shift;
$skip_trypsin_opt = shift;

$skip_trypsin_flag = 1;
if (defined($skip_trypsin_opt))
{
    if ($skip_trypsin_opt =~ /keep/i)
    {
        $skip_trypsin_flag = 0;
    }
    elsif ($skip_trypsin_opt =~ /omit/i)
    {
        $skip_trypsin_flag = 1;
    }
    elsif ($skip_trypsin_opt =~ /skip/i)
    {
        $skip_trypsin_flag = 1;
    }
    elsif ($skip_trypsin_opt =~ /discard/i)
    {
        $skip_trypsin_flag = 1;
    }
}


open INFILE, "$infile" or die "can't open $infile\n";

# read in header line
$line = <INFILE>;
$line =~ s/[\r\n]+//g;
$line =~ s/\"//g;

@array = split /\t/, $line;

for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    $array[$i] =~ s/\s+/ /g;
    
    $header_to_col_hash{$array[$i]} = $i;

    # idpicker uses accession column as identifier column
    if ($array[$i] eq 'Accession' || $array[$i] eq 'Protein IDs')
    {
        $id_col = $i;
    }
    elsif ($array[$i] eq 'Description')
    {
        $desc_col = $i;
        $array[$i] = 'Description IDPicker';
    }
    elsif ($array[$i] eq 'Gene Id')
    {
        # this is usually actually the gene symbol, *not* the Entrez GeneID
        $symbol_col = $i;
        
        $array[$i] = 'Symbol IDPicker';
    }
    # clean up messed up sample names, IRON treats '/' as path delimiters
    else
    {
        if ($array[$i] =~ /TMT\-\d\d\d/)
        {
            $array[$i] =~ s/[^A-Za-z0-9-]+/_/g;
            $array[$i] =~ s/_+/_/g;
            $array[$i] =~ s/^_+//;
            $array[$i] =~ s/_+$//;
            
            # flag channel summary columns for removal later
            if ($array[$i] =~ /^TMT\-\d\d\d\w*/ ||
                $array[$i] =~ /^[^_]+_TMT\-\d\d\d\w*/)
            {
                $strip_col_flags{$i} = 1;
            }
        }
    }
    
    $header_array[$i] = $array[$i];
}


# check to see if it has modification sites, use a different ID column
$site_flag = 0;
for ($i = 0; $i < @array; $i++)
{
    # idpicker uses accession column as identifier column
    if ($array[$i] eq 'ModificationID')
    {
        $id_col = $i;
        $site_flag = 1;
    }
    if ($array[$i] =~ /Localization prob/i)
    {
        $site_flag = 1;
    }
}


$peptide_flag = '';
$seq_col   = $header_to_col_hash{'Sequence'};
$nterm_col = $header_to_col_hash{'N-term cleavage window'};
$cterm_col = $header_to_col_hash{'C-term cleavage window'};
if (defined($seq_col) && defined($nterm_col) && defined($cterm_col))
{
    $peptide_flag = 'maxquant';
}

# hack to detect IDPicker peptides
if (!defined($id_col) && $peptide_flag eq '' &&
    defined($seq_col))
{
    $peptide_flag = 'idpicker';
}


# it might be an un-blessed MaxQuant file
if (!defined($id_col) && $peptide_flag eq 'maxquant')
{
    for ($i = 0; $i < @array; $i++)
    {
        if (!defined($id_col) && $array[$i] =~ /^Sequence$/i)
        {
            $id_col = $i;
        }
    }
}
if (!defined($id_col))
{
    for ($i = 0; $i < @array; $i++)
    {
        if (!defined($id_col) && $array[$i] =~ /^Leading proteins$/i)
        {
            $id_col = $i;
        }
    }
}
if (!defined($id_col))
{
    for ($i = 0; $i < @array; $i++)
    {
        if (!defined($id_col) && $array[$i] =~ /^Proteins$/i)
        {
            $id_col = $i;
        }
    }
}
if (!defined($id_col))
{
    for ($i = 0; $i < @array; $i++)
    {
        if (!defined($id_col) && $array[$i] =~ /^Protein$/i)
        {
            $id_col = $i;
        }
    }
}


# MaxQuant stuff
$pep_col      = $header_to_col_hash{'PEP'};
$mass_err_col = $header_to_col_hash{'Mass Error [ppm]'};
#$score_col   = $header_to_col_hash{'Score'};

# columns to be removed, since they can be huge and cause Excel to puke
if (defined($header_to_col_hash{'Evidence IDs'}))
{
    $strip_col_flags{$header_to_col_hash{'Evidence IDs'}} = 1;
}
if (defined($header_to_col_hash{'MS/MS IDs'}))
{
    $strip_col_flags{$header_to_col_hash{'MS/MS IDs'}}    = 1;
}

# Proteome Discoverer
for ($i = 0; $i < @array; $i++)
{
    $header = $array[$i];

    if ($header =~ /^Found in Sample/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}}    = 1;
    }
    if ($header =~ /CV \[%\]/)
    {
        $strip_col_flags{$header_to_col_hash{$header}}    = 1;
    }
}


# if there is just one plex, maxquant DOES NOT print the individual
# plex columns, and instead only prints the summary columns
# so, we have to detect the format change and dynamically adjust
$multi_plex_flag = 0;
$max_channel = -9E99;		# figure out how many channels we have
$min_channel =  9E99;		# figure out how many channels we have
foreach $header (keys %header_to_col_hash)
{
    if ($header =~ /^Reporter intensity (\d+) .*?$/i)
    {
        $channel = $1;
        
        if ($channel > $max_channel)
        {
            $max_channel = $channel;
        }
        if ($channel < $min_channel)
        {
            $min_channel = $channel;
        }

        $multi_plex_flag = 1;
    }
    
    if ($header =~ /^Reporter intensity (\d+)$/i)
    {
        $channel = $1;
        
        if ($channel > $max_channel)
        {
            $max_channel = $channel;
        }
        if ($channel < $min_channel)
        {
            $min_channel = $channel;
        }
    }
}


# recent Maxquant nonsense where they changed the channel numbering
#  convention from 0 to N-1 to 1 to N.
#
# ** DAMN IT, THEY CHANGED IT BACK AGAIN !! **
#  current code should still handle it, though
$max_channel -= $min_channel;


# Uh oh, this isn't a MaxQuant formatted TMT file
%channel_order_hash = ();
$num_seen_channels = 0;
if ($min_channel == 9E99)
{
    $max_channel = -9E99;

    # Proteome Discoverer
    for ($i = 0; $i < @array; $i++)
    {
        $header = $array[$i];
    
        if ($header =~ /^Abundances \(Grouped\):\s*([0-9NC]+)/)
        {
            $channel = $1;

            if (!defined($channel_order_hash{$channel}))
            {
                $channel_order_hash{$channel} = $num_seen_channels;
                $num_seen_channels++;
            }
        }
    }
    $max_channel = $num_seen_channels - 1;

    if ($max_channel >= 0)
    {
        $min_channel = 0;
    }
}


# fill channel string map for appropriate plex size
if (defined($channel_map_table[$max_channel+1]))
{
    for ($i = 0; $i <= $max_channel; $i++)
    {
        $channel_map[$i] = $channel_map_table[$max_channel+1][$i];
    }
}
# uh oh, we haven't seen this one, just use the original numbers
else
{
    printf STDERR "WARNING -- unsupported TMT-# plex; using original channel numbers";

    for ($i = 0; $i <= $max_channel; $i++)
    {
        $channel_map[$i] = $i;
    }
}


# other not-so-useful MaxQuant columns to be removed
# Damn it Maxquant! Stop changing the case of your column headers!!
foreach $header (keys %header_to_col_hash)
{
#    if ($header =~ /^Peptides\b/i)
#    {
#        $strip_col_flags{$header_to_col_hash{$header}} = 1;
#    }
#    if ($header =~ /^Razor \+ unique peptides\b/i)
#    {
#        $strip_col_flags{$header_to_col_hash{$header}} = 1;
#    }
#    if ($header =~ /^Unique peptides\b/i)
#    {
#        $strip_col_flags{$header_to_col_hash{$header}} = 1;
#    }
    if ($header =~ /^Fraction\b/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
    if ($header =~ /^Intensity /i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
    if ($header =~ /^Intensity$/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
    if ($header =~ /^Intensity__/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
    if ($header =~ /^iBAQ /i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
    if ($header =~ /^iBAQ$/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
    if ($header =~ /^iBAQ__/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
    
    # in the case of a single plex, Maxquant doesn't print individual
    # sample data, so we have to pull it out of the summary columns....
    if ($multi_plex_flag)
    {
      # summary columns, not actual sample data
      if ($header =~ /^Reporter intensity \d+$/i)
      {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
      }
      # summary columns, not actual sample data
      if ($header =~ /^Reporter intensity not corrected \d+$/i)
      {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
      }
      # summary columns, not actual sample data
      if ($header =~ /^Reporter intensity corrected \d+$/i)
      {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
      }
      # summary columns, not actual sample data
      if ($header =~ /^Reporter intensity count \d+$/i)
      {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
      }
    }

    # We don't want to use the un-corrected signals
    # Even though we don't really know what sort of correction MaxQuant might
    #  be doing, we can hope that it is reasonable and non-destructive...
    if ($header =~ /^Reporter intensity not corrected/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
    # we don't want to use the designated corrected columns, if they exist,
    #  either, since Bin thinks they may be bleed-over correction, which in
    #  our experience does more harm than good
    if ($header =~ /^Reporter intensity corrected/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }

    # skip counts
    if ($header =~ /^Reporter intensity count/i ||
        $header =~ / Count$/i ||
        $header =~ /^Experiment /i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
    
    # more per-run summaries
    if ($header =~ /^Sequence coverage [\w+].*?\[\%\]$/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }

    # misc stuff at the end
    if ($header =~ /^Peptide IDs$/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
    if ($header =~ /^Peptide is razor$/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
    if ($header =~ /^Mod. peptide IDs$/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
    if ($header =~ /^Best MS\/MS$/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
    if ($header =~ /^Oxidation \(M\) site/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
    if ($header =~ /^Occupancy/i)
    {
        $strip_col_flags{$header_to_col_hash{$header}} = 1;
    }
}


@intensities_col_hash = ();
foreach $header (keys %header_to_col_hash)
{
    # Damn it Maxquant! Stop changing the case of your column headers!!
    
    # previous blessing script may have already addressed 1-plex issues
    # so don't require this check to only be applied to multi-plex
    if ($header =~ /^Reporter intensity (\d+) (.*?)$/i)
    {
        $channel_number = $1;
        $plex           = $2;
        
        # recent Maxquant nonsense where they changed the channel numbering
        #  convention from 0 to N-1 to 1 to N.
        $channel_number -= $min_channel;
        
        $channel        = $channel_map[$channel_number];

        # reformat channel to the convention used in the lung squamous
        #  protegenomics project
        # assume that there are purely technical replicates of a single
        #  biological sample, so call them different TMTs to prevent them
        #  from being averaged together later in the pipeline
        if ($plex =~ /^TMTrun(\d+)$/i)
        {
            $replicate = $1;
#            $plex = sprintf "TMT_%d", $replicate;
#            $plex = sprintf "TMT%d", $replicate;
            $plex = sprintf "Plex1_run%d", $replicate;
        }

        elsif ($plex =~ /^Run\s*(\d+)$/i)
        {
            $replicate = $1;
#            $plex = sprintf "TMT_%d", $replicate;
#            $plex = sprintf "TMT%d", $replicate;
            $plex = sprintf "Plex1_run%d", $replicate;
        }
        
        $sample_name = sprintf "%s_TMT-%s",
            $plex, $channel;
        
        $array[$header_to_col_hash{$header}] = $sample_name;
        
        $intensities_col_hash{$header_to_col_hash{$header}} = $sample_name;
    }
    # use the summary columns instead
    if ($multi_plex_flag == 0 && $header =~ /^Reporter intensity (\d+)$/i)
    {
        $channel_number = $1;
        $plex           = 'Plex1';

        # recent Maxquant nonsense where they changed the channel numbering
        #  convention from 0 to N-1 to 1 to N.
        $channel_number -= $min_channel;

        $channel        = $channel_map[$channel_number];
        
        # reformat channel to the convention used in the lung squamous
        #  protegenomics project
        # assume that there are purely technical replicates of a single
        #  biological sample, so call them different TMTs to prevent them
        #  from being averaged together later in the pipeline
        if ($plex =~ /^TMTrun(\d+)$/i)
        {
            $replicate = $1;
#            $plex = sprintf "TMT_%d", $replicate;
#            $plex = sprintf "TMT%d", $replicate;
            $plex = sprintf "Plex1_run%d", $replicate;
        }

        $sample_name = sprintf "%s_TMT-%s", $plex, $channel;
        
        $array[$header_to_col_hash{$header}] = $sample_name;
        
        $intensities_col_hash{$header_to_col_hash{$header}} = $sample_name;
    }
    
    # Proteome Discoverer
    #
    # I only have single-plex data so far
    # I don't know what multiple plex data looks like yet
    #
    if ($header =~ /^Abundances \(Grouped\):\s*([0-9NC]+)/)
    {
        $channel_orig   = $1;
        $channel_number = $channel_order_hash{$channel_orig};

        $plex           = 'Plex1';
        $channel        = $channel_map[$channel_number];

        $sample_name = sprintf "%s_TMT-%s", $plex, $channel;
        
        $array[$header_to_col_hash{$header}] = $sample_name;
        
        $intensities_col_hash{$header_to_col_hash{$header}} = $sample_name;
    }
}
@intensities_col_array = sort keys %intensities_col_hash;


# insert new row identifier column
printf "%s\t", 'RowIdentifier';
printf "%s\t", 'OrigIdentifier';

# print header, minus channel summary columns
$i = 0;
while (defined($strip_col_flags{$i}))
{
    $i++;
}
printf "%s", $array[$i];
$i++;

# insert RefType column
printf "\t%s", 'RefType';

for (; $i < @array; $i++)
{
    if (!defined($strip_col_flags{$i}))
    {
        printf "\t%s", $array[$i];
    }
}
printf "\n";


# clean the rest of the data
$row = 0;
while(defined($line=<INFILE>))
{
    $line =~ s/[\r\n]+//g;

    @array = split /\t/, $line, -1;

    for ($i = 0; $i < @array; $i++)
    {
        # remove mis-escaped quotes
        if ($array[$i] =~ /^\"/ && $array[$i] =~ /\"$/)
        {
            $array[$i] =~ s/^\"//;
            $array[$i] =~ s/\"$//;
        }

        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\s+/ /g;
    }
    
    $orig_row_id = $array[$id_col];

    # MaxQuant often neglects to carry over REV_ and CON_ from one column
    # to another, so we have to scan them all, then add them where needed
    %maxquant_con_hash = ();
    %maxquant_rev_hash = ();
    for ($i = 0; $i < @array; $i++)
    {
        # clean up GenBank junk
        if ((defined($id_col) && $i == $id_col &&
             !($header_array[$i] =~ /ModificationID/i)) ||
            $header_array[$i] =~ /^Protein/i ||
            $header_array[$i] =~ /protein$/i ||
            $header_array[$i] =~ /\bprotein ids\b/i ||
            $header_array[$i] =~ /^Protein Accession/i)
        {
            $blessed_string = $array[$i];
            $blessed_string =~ s/\,/;/g;
            $blessed_string =~ s/\;+/;/;
            @split_array = split /\;/, $blessed_string;
            
            for ($j = 0; $j < @split_array; $j++)
            {
                $split_array[$j] =~ s/\|+$//;
                $split_array[$j] =~ s/^\|+//;

                # XXX_ at the beginning indicates a reverse match
                $reverse_flag = 0;
                if ($split_array[$j] =~ /^XXX_/)
                {
                    $reverse_flag = 1;
                }

                $split_array[$j] =~ s/^.*?gi\|\d+\|ref\|//;
                
                if ($reverse_flag)
                {
                    $split_array[$j] = sprintf "REV__%s", $split_array[$j];
                }
            }
            
            @split_array = sort compare_accession @split_array;

            $array[$i] = join ';', @split_array;
        }

        # clean up MaxQuant junk
        # holy crap is this bad
        if ((defined($id_col) && $i == $id_col &&
             !($header_array[$i] =~ /ModificationID/i)) ||
            $header_array[$i] =~ /^Protein/i ||
            $header_array[$i] =~ /protein$/i ||
            $header_array[$i] =~ /\bprotein ids\b/i ||
            $header_array[$i] =~ /^Protein Accession/i ||
            $header_array[$i] =~ /^Leading Proteins/i)
        {

            @split_array = split /;/, $array[$i];
            
#            $field = $array[$i];
            
            for ($j = 0; $j < @split_array; $j++)
            {
                # clean up sp| junk
#                $split_array[$j] =~ s/sp\|([^\|]+)\|[^\|]+$/$1/;
#                $split_array[$j] =~ s/sp\|([^\|]+)$/$1/;
                
                # clean up REFSEQ junk
#                $split_array[$j] =~ s/REFSEQ://g;

                # clean up H-INV junk
#                $split_array[$j] =~ s/H-INV://g;

                # clean up ENSEMBL junk
#                $split_array[$j] =~ s/ENSEMBL://g;
                
                if ($split_array[$j] =~ /CON__/)
                {
                    $split_array[$j] =~ s/CON__//g;

                    $maxquant_con_hash{$split_array[$j]} = 1;
                }
                if ($split_array[$j] =~ /REV__/)
                {
                    $split_array[$j] =~ s/REV__//g;

                    $maxquant_rev_hash{$split_array[$j]} = 1;
                }
            }
            
            $field_new = join ';', @split_array;
            
            $array[$i] = $field_new;
        }
        
        # strip Unmappable text
        if (defined($symbol_col) && $i == $symbol_col)
        {
            if ($array[$i] =~ /^Unmapped_/)
            {
                $array[$i] = '';
            }
        }
        
        # strip concatenated junk
        if (defined($desc_col) && $i == $desc_col)
        {
            $array[$i] =~ s/[^|]+\|[^|]+\|//;
            $array[$i] =~ s/^PREDICTED:\s+//;
#            $array[$i] =~ s/\s+\[Home sapiens]$//i;
        }
    }


    # add back in missing REV_ and CON_
    for ($i = 0; $i < @array; $i++)
    {
        if ((defined($id_col) && $i == $id_col &&
             !($header_array[$i] =~ /ModificationID/i)) ||
            $header_array[$i] =~ /^Protein/i ||
            $header_array[$i] =~ /protein$/i ||
            $header_array[$i] =~ /\bprotein ids\b/i ||
            $header_array[$i] =~ /^Protein Accession/i ||
            $header_array[$i] =~ /^Leading Proteins/i)
        {
            @split_array = split /;/, $array[$i];

            for ($j = 0; $j < @split_array; $j++)
            {
                $accession = $split_array[$j];

                if (defined($maxquant_con_hash{$accession}))
                {
                    $accession = 'CON__' . $accession;
                }
                if (defined($maxquant_rev_hash{$accession}))
                {
                    $accession = 'REV__' . $accession;
                }
                
                $split_array[$j] = $accession;
            }

            $field_new = join ';', @split_array;
            
            $array[$i] = $field_new;
        }
    }

    if ($skip_trypsin_flag)
    {
        # skip Bovine trypsin
        if (defined($id_col) && $array[$id_col] =~ /TRY1_BOVIN/)
        {
            next;
        }
        $skip_flag = 0;
        for ($i = 0; $i < @array; $i++)
        {
            if ($header_array[$i] =~ /^Protein/i ||
                $header_array[$i] =~ /protein$/i ||
                $header_array[$i] =~ /^Protein Accession/i)
            {
                if ($array[$i] =~ /TRY1_BOVIN/)
                {
                    $skip_flag = 1;
                    last;
                }
            }
        }
        if ($skip_flag)
        {
            next;
        }
    }

    
    # skip poor quality MaxQuant assignments
    # PCA behaves very poorly with p < 0.20; works pretty well for p < 0.05
#    if ($site_flag == 1 &&
#        defined($pep_col) && $array[$pep_col] > $pep_cutoff_site)
#    {
#        next;
#    }
#    if ($site_flag == 0 &&
#        defined($pep_col) && $array[$pep_col] > $pep_cutoff)
#    {
#        next;
#    }
    if (defined($pep_col) && $array[$pep_col] > $pep_cutoff)
    {
        next;
    }
    if (defined($mass_err_col) && $array[$mass_err_col] > $mass_err_cutoff)
    {
        next;
    }
    

    # skip lines with no intensities
    if (@intensities_col_array)
    {
        $has_non_zero_intensities_flag = 0;
        
        foreach $col (@intensities_col_array)
        {
            if ($array[$col] =~ /[1-9]/)
            {
                $has_non_zero_intensities_flag = 1;
                last;
            }
        }
        
        if ($has_non_zero_intensities_flag == 0)
        {
            next;
        }
    }
    

    # insert new peptide sequence
    if ($peptide_flag eq 'maxquant')
    {
        $seq    = $array[$seq_col];
        $nterm  = $array[$nterm_col];
        $cterm  = $array[$cterm_col];

        $seq_window = sprintf "%s.%s.%s", $nterm, $seq, $cterm;
        
        $seq_window =~ s/^_+/_/;
        $seq_window =~ s/_+$/_/;
        
        printf "%s\t", $seq_window;
    }
    elsif ($peptide_flag eq 'idpicker')
    {
        $seq    = $array[$seq_col];
        
        printf "%s\t", $seq;
    }
    # clone the id_col
    elsif (defined($id_col))
    {
        $value = $array[$id_col];
        
        # HACK -- Evince search is broken if the identifier has | in it
        # replace | with ~
        $value =~ s/\|/\~/g;
        
        printf "%s\t", $value;
    }
    else
    {
        printf STDERR "ERROR - something bad happened on row %d\n", $row;
        printf "%s\t", $row;

        $row++;
    }
    
    printf "%s\t", $orig_row_id;

    
    $i = 0;
    while (defined($strip_col_flags{$i}))
    {
        $i++;
    }
    printf "%s", $array[$i];
    $i++;
    
    # insert RefType column
    if ($peptide_flag ne '')
    {
      printf "\t%s", 'Peptide';
    }
    elsif (defined($id_col))
    {
      $accession = $array[$id_col];
      if ($accession =~ /REV__/)
      {
        printf "\t%s", '__REV';
      }
      elsif ($accession =~ /NP_/)
      {
        printf "\t%s", 'NP';
      }
      elsif ($accession =~ /XP_/)
      {
        printf "\t%s", 'XP';
      }
      elsif ($accession =~ /YP_/)
      {
        printf "\t%s", 'YP';
      }
      else
      {
        printf "\t%s", '__Other';
      }
    }
    # assume it is a peptide
    else
    {
      printf "\t%s", '__Other';
    }
    
    for (; $i < @array; $i++)
    {
        if (!defined($strip_col_flags{$i}))
        {
            printf "\t%s", $array[$i];
        }
    }
    printf "\n";
}
