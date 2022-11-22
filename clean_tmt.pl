#!/usr/bin/perl -w

# 2022-11-22: zero-pad Plex numbers in large experiments
# 2022-11-22: sort samples into plex + label order
# 2022-11-21: add support for multiple Proteome Discoverer plexes
# 2022-08-17: support TMT-18
# 2022-05-06: fixed sample channel detection bug introduced on 2022-04-21
# 2022-04-21: make sure split doesn't remove empty trailing fields
# 2022-04-21: begin adding support for Proteome Discoverer
# 2022-04-21: remove all double quotes, not just the first one per line (oops)
# 2021-11-02: preserve and reformat run# in single-plex injection replicates
# 2020-08-03: finish support for TMT16, improve unsupported #channels error
# 2020-07-09: add support for TMT16, add letters to the ends of ALL channels
# 2020-02-28: replace | with ~ in RowIdentifier, otherwise Evince is upset


use POSIX;

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


sub cmp_renamed_header_cols
{
    my $header_a    = $header_array_orig_order[$a];
    my $header_b    = $header_array_orig_order[$b];
    my $plex_a      = '';
    my $plex_b      = '';
    my $ch_a        = '';
    my $ch_b        = '';
    my $ch_a_digits = '';
    my $ch_b_digits = '';

    # may contain a _run1 or _run2 after the plex
    if ($header_a =~ /^Plex([0-9]+).*_TMT-([0-9CN]+)$/)
    {
        $plex_a = $1;
        $ch_a   = $2;

        $ch_a_digits = $ch_a;
        $ch_a_digits =~ s/[^0-9]//g;
    }
    if ($header_b =~ /^Plex([0-9]+).*_TMT-([0-9CN]+)$/)
    {
        $plex_b = $1;
        $ch_b   = $2;

        $ch_b_digits = $ch_b;
        $ch_b_digits =~ s/[^0-9]//g;
    }

    # sort non-samples by original order
    if ($plex_a eq '' && $plex_b eq '')
    {
        return ($a <=> $b);
    }
    
    # put samples last
    if ($plex_a ne '' && $plex_b eq '') { return  1; }
    if ($plex_a eq '' && $plex_b ne '') { return -1; }
    
    # sort samples by plex
    if ($plex_a < $plex_b) { return -1; }
    if ($plex_a > $plex_b) { return  1; }

    # sort samples by channel
    if ($ch_a_digits < $ch_b_digits) { return -1; }
    if ($ch_a_digits > $ch_b_digits) { return  1; }
    
    # N comes before C or empty
    if (  $ch_a =~ /N$/i  && !($ch_b =~ /N$/i)) { return -1; }
    if (!($ch_a =~ /N$/i) &&   $ch_b =~ /N$/i)  { return  1; }

    # shouldn't ever trigger, but if it does, leave them in place
    return ($a <=> $b);
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

$channel_map_table[18][0]  = '126C';	# 126C
$channel_map_table[18][1]  = '127N';
$channel_map_table[18][2]  = '127C';
$channel_map_table[18][3]  = '128N';
$channel_map_table[18][4]  = '128C';
$channel_map_table[18][5]  = '129N';
$channel_map_table[18][6]  = '129C';
$channel_map_table[18][7]  = '130N';
$channel_map_table[18][8]  = '130C';
$channel_map_table[18][9]  = '131N';
$channel_map_table[18][10] = '131C';
$channel_map_table[18][11] = '132N';
$channel_map_table[18][12] = '132C';
$channel_map_table[18][13] = '133N';
$channel_map_table[18][14] = '133C';
$channel_map_table[18][15] = '134N';
$channel_map_table[18][16] = '134C';
$channel_map_table[18][17] = '135N';	# 135N

$plex_size_last_ch_hash{6}  = '131C';   # 6 and 11 both have 131C last
$plex_size_last_ch_hash{10} = '131N';
$plex_size_last_ch_hash{11} = '131C';   # 6 and 11 both have 131C last
$plex_size_last_ch_hash{16} = '134N';
$plex_size_last_ch_hash{18} = '135N';


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

@array = split /\t/, $line;

for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s*\"(.*)\"$/$1/;
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
%channel_order_hash  = ();
@channel_order_array = ();
$temp_plex_str_hash  = ();
$multi_plex_flag = 0;
$max_channel = -9E99;		# figure out how many channels we have
$min_channel =  9E99;		# figure out how many channels we have
$pd_grouped_flag = 0;
$num_seen_channels = 0;         # currently only used for Proteome Discoverer
foreach $header (@header_array)
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
    
    # Proteome Discoverer
    if ($header =~ /^Abundances \(Grouped\):\s*([0-9A-Za-z]+)/)
    {
        $ch_str   = $1;

        $plex_num = 1;
        if ($header =~ /^Abundances \(Grouped\):\s*([0-9A-Za-z]+),\s*([0-9]+)/)
        {
            $plex_num = $2;
        }

        # add missing C labels
        # 126 always leaves it off, is probably missing for all 6-plex ch as well
        if (!($ch_str =~ /[A-Za-z]$/i))
        {
            $ch_str .= 'C';
        }
        

        ## DEBUG -- skip a channel for testing purposes
        ## don't forget to comment this out after debugging
        #if ($ch_str eq '128N') { next; }
        #if ($ch_str eq '134N') { next; }


        if (!defined($channel_order_hash{$ch_str}))
        {
            $channel_order_hash{$ch_str} = $num_seen_channels;
            $channel_order_array[$num_seen_channels] = $ch_str;
            $num_seen_channels++;
        }

        $temp_plex_str_hash{$plex_num} = 1;
        
        $pd_grouped_flag = 1;
    }
}


# Proteome Discoverer, determine number of channels, whether multiple plexes
if ($pd_grouped_flag)
{
    @temp_array = keys %temp_plex_str_hash;
    if (@temp_array > 1)
    {
        $multi_plex_flag = 1;
    }
    
    
    # sanity check for missing channels at the end
    $last_label = $channel_order_array[$num_seen_channels - 1];
    
    # scan labels for N to differentiate TMT-6 from TMT-11
    $n_flag = 0;
    foreach $ch_str (@channel_order_array)
    {
        if ($ch_str =~ /N$/i)
        {
            $n_flag = 1;
        }
    }
    
    # scan plex sizes for channel in reverse order
    $smallest_size   = 0;
    @plex_size_array = sort {$b<=>$a} keys %plex_size_last_ch_hash;
    foreach $plex_size (@plex_size_array)
    {
        # early exit, all remaining plexs are too small
        if ($plex_size < $num_seen_channels)
        {
            last;
        }
    
        # skip 6-plex, since there are no N-labels
        if ($n_flag && $plex_size == 6)
        {
            next;
        }

        # scan each plex size to see if contains the last channel
        # scan all channels, in case the true last channel is missing
        for ($ch = $plex_size - 1; $ch >= 0; --$ch)
        {
            $ch_lookup = $channel_map_table[$plex_size][$ch];
            
            #printf STDERR "DEBUG\t%s\t%s\t%s\n",
            #    $plex_size, $last_label, $ch_lookup;
            
            if ($last_label eq $ch_lookup)
            {
                $smallest_size = $plex_size;
                
                last;
            }
        }
    }
    
    if ($smallest_size)
    {
        $max_channel = $smallest_size - 1;
        $min_channel = 0;
        
        # re-fill channel order, with missing channels inserted
        %channel_order_hash  = ();
        @channel_order_array = ();
        for ($i = 0; $i < $smallest_size; $i++)
        {
            $ch_str = $channel_map_table[$smallest_size][$i];
            $channel_order_hash{$ch_str} = $i;
            $channel_order_array[$i]     = $ch_str;
        }
    }
    
    #printf STDERR "DEBUG\t%s\t%s\t%g\n",
    #    $max_channel+1, $num_seen_channels, $smallest_size;
}


# recent Maxquant nonsense where they changed the channel numbering
#  convention from 0 to N-1 to 1 to N.
#
# ** DAMN IT, THEY CHANGED IT BACK AGAIN !! **
#  current code should still handle it, though
$max_channel -= $min_channel;



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
foreach $header (sort keys %header_to_col_hash)
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
    if ($header =~ /^Abundances \(Grouped\):\s*([0-9NC]+)/)
    {
        $channel_orig = $1;

        # add missing C labels
        # 126 always leaves it off, is probably missing for all 6-plex ch as well
        if (!($channel_orig =~ /[A-Za-z]$/i))
        {
            $channel_orig .= 'C';
        }

        $channel_number = $channel_order_hash{$channel_orig};


        $plex_num = 1;
        if ($header =~ /^Abundances \(Grouped\):\s*([0-9A-Za-z]+),\s*([0-9]+)/)
        {
            $plex_num = $2;
        }

        $plex           = 'Plex' . $plex_num;
        $channel        = $channel_map[$channel_number];

        $sample_name = sprintf "%s_TMT-%s", $plex, $channel;
        
        $array[$header_to_col_hash{$header}] = $sample_name;
        
        $intensities_col_hash{$header_to_col_hash{$header}} = $sample_name;
    }
}
@intensities_col_array = sort {$a<=>$b} keys %intensities_col_hash;


# Scan for missing channels
#
# It turns out there weren't any missing channels after all,
# so all this new missing channel code is going to waste.
# Oh well, it's good to be extra paranoid when the samples aren't arranged
# in the correct order within the output file...
#
# I haven't written any code to insert blank columns at the end yet,
# I'll only implement that if I ever see any missing channel warnings.
#
%seen_sample_hash     = ();
%seen_plex_hash       = ();
%missing_sample_hash  = ();
@missing_sample_array = ();
$num_missing_samples  = 0;
if ($pd_grouped_flag)
{
    foreach $col (@intensities_col_array)
    {
        $header = $array[$col];

        if ($header =~ /Plex([0-9]+)_TMT-([0-9CN]+)/)
        {
            $plex   = $1;
            $ch_str = $2;

            $sample_name = sprintf "Plex%s_TMT-%s", $plex, $ch_str;

            $seen_sample_hash{$sample_name} = 1;
            $seen_plex_hash{$plex} = 1;
        }
    }

    @seen_plex_array = sort {$a<=>$b} keys %seen_plex_hash;
    foreach $plex (@seen_plex_array)
    {
        for ($i = 0; $i <= $max_channel; $i++)
        {
            $ch_str      = $channel_map_table[$max_channel+1][$i];
            $sample_name = sprintf "Plex%s_TMT-%s", $plex, $ch_str;

            if (!defined($seen_sample_hash{$sample_name}))
            {
                $missing_sample_hash{$sample_name} = 1;
                
                $missing_sample_array[$num_missing_samples++] =
                    $sample_name;
            }
        }
    }
    
    foreach $sample_name (@missing_sample_array)
    {
        printf "Missing sample:\t%s\n", $sample_name;
    }
}


# determine maximum width to zero-pad Plex numbers in large experiments
$max_plex = 1;
for ($i = 0; $i < @array; $i++)
{
    $header = $array[$i];

    if ($header =~ /^Plex([0-9]+)(.*_TMT-[0-9CN]+)$/)
    {
        $plex    = $1;
        $end_str = $2;
        
        if ($plex > $max_plex)
        {
            $max_plex = $plex;
        }
    }
}
$max_plex_digits = floor(log10($max_plex)) + 1;

# zero-pad Plex numbers in large experiments
if ($max_plex_digits > 1)
{
    for ($i = 0; $i < @array; $i++)
    {
        $header = $array[$i];

        if ($header =~ /^Plex([0-9]+)(.*_TMT-[0-9CN]+)$/)
        {
            $plex_num = $1;
            $end_str  = $2;
            
            $plex_new = sprintf "Plex%0*d", $max_plex_digits, $plex_num;
            
            $array[$i] =~ s/^Plex$plex_num/$plex_new/;
        }
    }
}


# get ready to reorder the headers
@header_array_orig_order     = @array;
@header_col_new_order_array  = ();
for ($i = 0; $i < @header_array_orig_order; $i++)
{
    @header_col_new_order_array[$i] = $i;
}

@header_col_new_order_array =
    sort cmp_renamed_header_cols @header_col_new_order_array;


# insert new row identifier column
printf "%s\t", 'RowIdentifier';
printf "%s\t", 'OrigIdentifier';

# print header, minus channel summary columns
$i   = 0;
$col = $header_col_new_order_array[$i];
while (defined($strip_col_flags{$col}))
{
    $i++;
    $col = $header_col_new_order_array[$i];
}
printf "%s", $header_array_orig_order[$col];
$i++;
$col = $header_col_new_order_array[$i];

# insert RefType column
printf "\t%s", 'RefType';

for (; $i < @array; $i++)
{
    $col = $header_col_new_order_array[$i];

    if (!defined($strip_col_flags{$col}))
    {
        printf "\t%s", $header_array_orig_order[$col];
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
        $array[$i] =~ s/^\s*\"(.*)\"$/$1/;
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
             !($header_array_orig_order[$i] =~ /ModificationID/i)) ||
            $header_array_orig_order[$i] =~ /^Protein/i ||
            $header_array_orig_order[$i] =~ /protein$/i ||
            $header_array_orig_order[$i] =~ /\bprotein ids\b/i ||
            $header_array_orig_order[$i] =~ /^Protein Accession/i)
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
             !($header_array_orig_order[$i] =~ /ModificationID/i)) ||
            $header_array_orig_order[$i] =~ /^Protein/i ||
            $header_array_orig_order[$i] =~ /protein$/i ||
            $header_array_orig_order[$i] =~ /\bprotein ids\b/i ||
            $header_array_orig_order[$i] =~ /^Protein Accession/i ||
            $header_array_orig_order[$i] =~ /^Leading Proteins/i)
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
             !($header_array_orig_order[$i] =~ /ModificationID/i)) ||
            $header_array_orig_order[$i] =~ /^Protein/i ||
            $header_array_orig_order[$i] =~ /protein$/i ||
            $header_array_orig_order[$i] =~ /\bprotein ids\b/i ||
            $header_array_orig_order[$i] =~ /^Protein Accession/i ||
            $header_array_orig_order[$i] =~ /^Leading Proteins/i)
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
            if ($header_array_orig_order[$i] =~ /^Protein/i ||
                $header_array_orig_order[$i] =~ /protein$/i ||
                $header_array_orig_order[$i] =~ /^Protein Accession/i)
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
    $col = $header_col_new_order_array[$i];
    while (defined($strip_col_flags{$col}))
    {
        $i++;
        $col = $header_col_new_order_array[$i];
    }
    printf "%s", $array[$col];
    $i++;
    $col = $header_col_new_order_array[$i];
    
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
        $col = $header_col_new_order_array[$i];

        if (!defined($strip_col_flags{$col}))
        {
            printf "\t%s", $array[$col];
        }
    }
    printf "\n";
}
