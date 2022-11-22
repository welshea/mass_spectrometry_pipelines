#!/usr/bin/perl -w

# Changelog:
#
# 2022-11-21: add support for multiple Proteome Discoverer plexes
# 2022-11-21: remove first level of enclosing double quotes from fields
# 2022-08-17: support TMT-18
# 2022-04-21: begin adding Proteome Discoverer support
# 2022-01-19: default 100% injection replicates to "auto" reference channel
# 2021-08-12: change deprecated "if (defined(%plex_hash))" to "if (%plex_hash)"
# 2021-07-28: respect species argument override, autodetect if non-mouse/human
# 2021-05-21: added --boost and --last-ch flags to force last channel norm
#             print Usage statement
#             use last argument as species, default to human
#             added autodetection of species based on SwissProt annotation
# 2020-12-18: detect multiple injection replicates of a single plex
# 2020-09-10: 2020-08-03 #channels test broke non-TMT data
# 2020-08-03: fix pY last channel for TMT16, abort on unsupported #channels
# 2020-07-09: add support for TMT16, add letters to the ends of ALL channels
#
# 2020-02-14: better internal handling of ___#, no actual change in results
#
# 2020-01-22: dealt with new edge case where both Reporter intensity
#             and Reporter intensity Samplename were present for single plex
#             (this never occured before...)

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

$intensity_flag  = 0;
$ibaq_flag       = 0;
$refseq_flag     = 0;
$mod_flag        = 0;
$psty_flag       = 0;
$py_flag         = 0;
$boost_flag      = 0;
$tmt_flag        = 0;
$multi_plex_flag = 0;
$species         = '';

$filename = '';
$num_files = 0;
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
        # treat non-option following filename as the reference species
        elsif ($species eq '')
        {
            $species = $field
        }
    }
}

if ($filename eq '')
{
    $syntax_error_flag = 1;
}

if ($syntax_error_flag)
{
    printf "Usage: autodetect_ms_data.pl [options] maxquant_output.txt [species]\n";
    printf "\n";
    printf "  Options:\n";
    printf "    --boost     use highest channel for normalization\n";
    printf "    --last-ch   use highest channel for normalization\n";
    printf "\n";
    printf "  Output:\n";
    printf "    Modification  [yes/no]         does the data contain modification sites\n";
    printf "    RefSeq        [yes/no]         is the data mainly matched against RefSeqs\n";
    printf "    TMT           [single/multi/injection]\n";
    printf "                                   1 plex, >=2 plexes, >=2 all injection reps\n";
    printf "    TMT_Channel   [auto/TMT-????]  reference channel for IRON normalization\n";
    printf "    Rollup        [ibaq/intensity/intensity_and_ibaq]\n";
    printf "                                   which type of rollup columns to keep\n";
    printf "    Species       [human/mouse/human_and_mouse]\n";

    exit(1);
}


open INFILE, "$filename" or die "can't open input file $filename\n";

$line = <INFILE>;
$line =~ s/[\r\n]+//g;
#$line =~ s/\"//g;

@header_array = split /\t/, $line;

for ($i = 0; $i < @header_array; $i++)
{
    $header_array[$i] =~ s/^\s*\"(.*)\"$/$1/;
    $header_array[$i] =~ s/^\s+//;
    $header_array[$i] =~ s/\s+$//;
    $header_array[$i] =~ s/\s+/ /g;
    
#    $header_to_col_hash{$header_array[$i]} = $i;
}

$max_channel = -9E99;		# figure out how many channels we have
$min_channel =  9E99;		# figure out how many channels we have
for($i = 0; $i < @header_array; $i++)
{
    $header = $header_array[$i];
    
    if ($header =~ /^Reporter intensity (\d+) (.*?)(___[123])*$/i)
    {
        $channel = $1;
        $plex    = $2;
        
        if ($channel > $max_channel)
        {
            $max_channel = $channel;
        }
        if ($channel < $min_channel)
        {
            $min_channel = $channel;
        }

        $tmt_flag        = 1;
        
        $plex_hash{$plex} = 1;
    }
    
    # deal with ___# bug
    if ($header =~ /^Reporter intensity (\d+)(___[123])*$/i)
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
        
        $tmt_flag = 1;
    }

    # Proteome Discoverer
    if ($header =~ /^Abundances \(Grouped\):\s*/)
    {
        $plex_num = 1;
        if ($header =~ /^Abundances \(Grouped\):\s*([0-9A-Za-z]+),\s*([0-9]+)/)
        {
            $plex_num = $2;
        }

        $plex             = 'Plex' . $plex_num;
        $plex_hash{$plex} = 1;

        $tmt_flag         = 1;
        $intensity_flag   = 1;
    }

    if ($header =~ /^Intensity \S/i)
    {
        $intensity_flag = 1;
    }
    
    if ($header =~ /^iBAQ\s/i)
    {
        $ibaq_flag = 1;
    }
    
    if ($header =~ /^Amino Acid$/i)
    {
        $mod_flag = 1;
        $aa_col = $i;
    }
    
    if ($header =~ /^Phospho \(STY\)/i)
    {
        $psty_flag = 1;
    }
}


# Uh oh, this isn't a MaxQuant formatted TMT file
if ($tmt_flag && $max_channel < 0)
{
    %channel_order_hash = ();
    $num_seen_channels = 0;
    $max_channel = -9E99;
    if ($min_channel == 9E99)
    {
        # Proteome Discoverer
        for ($i = 0; $i < @header_array; $i++)
        {
            $header = $header_array[$i];
        
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
}

# check to see if we have multiple plexes
if (%plex_hash)
{
    @plex_array = sort keys %plex_hash;
    
    if (@plex_array > 1)
    {
        $multi_plex_flag = 1;
    }
}

# if we have multiple plexes, check to see if they are 100% injection reps
# assume injection replicates end in case-insensitive run#
$injection_plex_flag = 0;
if ($multi_plex_flag)
{
    %temp_hash  = ();
    
    foreach $plex (@plex_array)
    {
        $plex =~ s/[-_. ]*run[-_. ]*[0-9]+$//i;
        $temp_hash{$plex} = 1;
    }
    
    @temp_array = sort keys %temp_hash;
    if (@temp_array == 1)
    {
        $injection_plex_flag = 1;
    }
}


# scan for details
%accession_hash = ();
$count_human = 0;
$count_mouse = 0;
$count_rows = 0;
while(defined($line=<INFILE>))
{
    $line =~ s/[\r\n]+//g;

    @array = split /\t/, $line;
    
    # check for likely SwissProt accession annotation, guess species
    if ($line =~ /(_human|Homo sapiens)/i)
    {
        $count_human++;
    }
    if ($line =~ /(_mouse|Mus musculus)/i)
    {
        $count_mouse++;
    }
    if ($line =~ /\S/)
    {
        $count_rows++;
    }

    for ($i = 0; $i < @array; $i++)
    {
        # remove mis-escaped quotes
        if ($array[$i] =~ /^\"/ && $array[$i] =~ /\"$/)
        {
            $array[$i] =~ s/^\"//;
            $array[$i] =~ s/\"$//;
        }

        $array[$i] =~ s/^\s*\"(.*)\"$/$1/;
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\s+/ /g;
    }
    
    for ($i = 0; $i < @array; $i++)
    {
        # clean up GenBank junk
        if ($header_array[$i] =~ /^Protein/i ||
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
            
            $array[$i] = join ';', @split_array;
        }

        # clean up MaxQuant junk
        # holy crap is this bad
        if ($header_array[$i] =~ /^Protein/i ||
            $header_array[$i] =~ /protein$/i ||
            $header_array[$i] =~ /\bprotein ids\b/i ||
            $header_array[$i] =~ /^Protein Accession/i ||
            $header_array[$i] =~ /^Leading Proteins/i)
        {

            @split_array = split /;/, $array[$i];
            
            for ($j = 0; $j < @split_array; $j++)
            {
                # clean up sp| junk
                $split_array[$j] =~ s/sp\|([^\|]+)\|[^\|]+$/$1/;
                $split_array[$j] =~ s/sp\|([^\|]+)$/$1/;
                
                # clean up REFSEQ junk
                $split_array[$j] =~ s/REFSEQ://g;

                # clean up new, different, ref| junk
                $split_array[$j] =~ s/ref\|//g;

                # clean up H-INV junk
                $split_array[$j] =~ s/H-INV://g;

                # clean up ENSEMBL junk
                $split_array[$j] =~ s/ENSEMBL://g;
                
                if ($split_array[$j] =~ /CON__/)
                {
                    $split_array[$j] =~ s/CON__//g;
                }
                if ($split_array[$j] =~ /REV__/)
                {
                    $split_array[$j] =~ s/REV__//g;
                }
            }
            
            $field_new = join ';', @split_array;

            foreach $accession (@split_array)
            {
                $accession_hash{$accession} = 1;
            }
            
            $array[$i] = $field_new;
        }
        
        if (defined($aa_col) && $i == $aa_col)
        {
            $field = $array[$i];
        
            if (!defined($psty_counts_hash{$field}))
            {
                $psty_counts_hash{$field} = 0;
            }
            
            $psty_counts_hash{$field} += 1;
        }
    }
}

@accession_array = sort keys %accession_hash;

$non_refseq_count = 0;
$refseq_count = 0;
foreach $accession (@accession_array)
{
    if ($accession =~ /^NP_/ || $accession =~ /^XP_/ ||
        $accession =~ /^YP_/)
    {
        $refseq_count++;
    }
    else
    {
        $non_refseq_count++;
    }
    
#    print "$accession\n";
}

if ($refseq_count > $non_refseq_count)
{
    $refseq_flag = 1;
}


$pst_count = 0;
$py_count  = 0;

if (defined($psty_counts_hash{'S'}))
{
    $pst_count += $psty_counts_hash{'S'};
}
if (defined($psty_counts_hash{'T'}))
{
    $pst_count += $psty_counts_hash{'T'};
}
if (defined($psty_counts_hash{'Y'}))
{
    $py_count += $psty_counts_hash{'Y'};
}

if ($py_count > $pst_count)
{
    $py_flag = 1;
}



# recent Maxquant nonsense where they changed the channel numbering
#  convention from 0 to N-1 to 1 to N.
#
# ** DAMN IT, THEY CHANGED IT BACK AGAIN !! **
#  current code should still handle it, though

if ($tmt_flag)
{
    $max_channel -= $min_channel;

    if ($tmt_flag && !defined($channel_map_table[$max_channel+1]))
    {
        printf STDERR "ABORT -- unsupported number of channels: %d\n",
            $max_channel+1;
        exit(1);
    }
}

$tmt_channel = 'none';
if ($tmt_flag)
{
    $tmt_channel = 'auto';
}
if ($tmt_flag && $multi_plex_flag)
{
    $tmt_channel = 'TMT-' .
                   $channel_map_table[$max_channel+1][0]
}
if ($tmt_flag && ($boost_flag || $py_flag))
{
    $tmt_channel = 'TMT-' . $channel_map_table[$max_channel+1][$max_channel]
}
if ($tmt_flag && $injection_plex_flag)
{
    $tmt_channel = 'auto';
}

if ($ibaq_flag)
{
    $rollup = 'ibaq';
    
    if ($intensity_flag)
    {
        $rollup = 'intensity_and_ibaq';
    }
}
else
{
    $rollup = 'intensity';
}

# override intensity/ibaq flag
# any Intensity and iBAQ columns in TMT data are summary columns
if ($tmt_flag)
{
    $rollup = 'intensity';
}


$tmt_type = 'no';
if ($tmt_flag)
{
    $tmt_type = 'single';
}
if ($multi_plex_flag)
{
    $tmt_type = 'multi';
}
if ($injection_plex_flag)
{
    $tmt_type = 'injection';
}

$mod_string = 'no';
if ($mod_flag) { $mod_string = 'yes'; }
$refseq_string = 'no';
if ($refseq_flag) { $refseq_string = 'yes'; }


# guess species
if ($count_rows &&
    !($species =~ /Human/i || $species =~ /Mouse/i))
{
    $fraction_human = $count_human / $count_rows;
    $fraction_mouse = $count_mouse / $count_rows;
    
    ## printf STDERR "%s\t%s\n", $fraction_human, $fraction_mouse;
    
    # default to human
    $species = 'human';
    
    if ($fraction_human > 0.1)
    {
        $species = 'human';
    }
    if ($fraction_mouse > 0.1)
    {
        $species = 'mouse';
    }
    if ($fraction_human > 0.1 && $fraction_mouse > 0.1)
    {
        $species = 'human_and_mouse';
    }
}
# use the given species override
elsif ($species =~ /Human/i || $species =~ /Mouse/i)
{
    $species = $species;
}
# default to human
else
{
    $species = 'human';
}


printf "Modification\t%s\n",  $mod_string;
printf "RefSeq\t%s\n",        $refseq_string;
printf "TMT\t%s\n",           $tmt_type;
printf "TMT_Channel\t%s\n",   $tmt_channel;
printf "Rollup\t%s\n",        $rollup;
printf "Species\t%s\n",       $species;
