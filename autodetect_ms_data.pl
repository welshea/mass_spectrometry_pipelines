#!/usr/bin/perl -w

# Changelog:
#
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

$intensity_flag  = 0;
$ibaq_flag       = 0;
$refseq_flag     = 0;
$mod_flag        = 0;
$psty_flag       = 0;
$py_flag         = 0;
$tmt_flag        = 0;
$multi_plex_flag = 0;

$filename = shift;

open INFILE, "$filename" or die "can't open input file $filename\n";

$line = <INFILE>;
$line =~ s/[\r\n]+//g;
$line =~ s/\"//;

@header_array = split /\t/, $line;

for ($i = 0; $i < @header_array; $i++)
{
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


# check to see if we have multiple plexes
if (defined(%plex_hash))
{
    @plex_array = sort keys %plex_hash;
    
    if (@plex_array > 1)
    {
        $multi_plex_flag = 1;
    }
}


# scan for details
%accession_hash = ();
while(defined($line=<INFILE>))
{
    $line =~ s/[\r\n]+//g;

    @array = split /\t/, $line;

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
if ($tmt_flag && $py_flag)
{
    $tmt_channel = 'TMT-' . $channel_map_table[$max_channel+1][$max_channel]
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

$mod_string = 'no';
if ($mod_flag) { $mod_string = 'yes'; }
$refseq_string = 'no';
if ($refseq_flag) { $refseq_string = 'yes'; }

printf "Modification\t%s\n",  $mod_string;
printf "RefSeq\t%s\n",        $refseq_string;
printf "TMT\t%s\n",           $tmt_type;
printf "TMT_Channel\t%s\n",   $tmt_channel;
printf "Rollup\t%s\n",        $rollup;
