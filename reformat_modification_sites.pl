#!/usr/bin/perl -w

# 2022-05-06:  support missing Positions Within Proteins column (old MaxQuant)
# 2021-10-28:  add support for Biotin-HPDP
# 2021-08-30:  deal with missing accessions/positions in various fields

use POSIX;

sub cmp_phospho_probs
{
    my $value1;
    my $value2;

    # put observed sites first
    $value1 = defined($global_populated_sites_hash{$a});
    $value2 = defined($global_populated_sites_hash{$b});
    
    if ($value1 > $value2) { return -1; };
    if ($value2 > $value1) { return 1; };


    # then by probabilities
    $value1 = $global_sites_prob_hash{$a};
    $value2 = $global_sites_prob_hash{$b};
    
    if ($value1 > $value2) { return -1; };
    if ($value2 > $value1) { return 1; };
    

    # some probabilities, due to reported significant digits, can be
    # broken with score diffs

    $value1 = $global_sites_diff_hash{$a};
    $value2 = $global_sites_diff_hash{$b};

    if ($value1 > $value2) { return -1; };
    if ($value2 > $value1) { return 1; };
    
    return ($a <=> $b);
}


sub fix_phosphopeptide
{
    my $mod_type_char        = $_[0];
    my $site                 = $_[1];
    my $probSTY              = $_[2];
    my $diffSTY              = $_[3];
    my $positions_prot       = $_[4];
    my %mod_hash             = ();
    my @sites_list           = ();
    my @temp_array           = ();
    my @temp_array2          = ();
    my @much_greater_array   = ();
    my @greater_array        = ();
    my @equal_array          = ();
    my @multi_protein_positions_array = ();
    my $num_potential_sites  = 0;
    my $seq_len;
    my $aa_index;
    my $aa_index2;
    my $fix_flag;
    my $site_to_remove;
    my $seq_new;
    my $prob;
    my $min_prob;
    my $temp;
    my $i;
    my $j;
    my $c;
    my $old_c;
    my $value1;
    my $value2;
    my $seq_unmodified;

    my $position_prot;
    my $protein_peptide_pos_diff;
    my $positions_prot_order;

    %global_populated_sites_hash = ();
    %global_sites_prob_hash = ();
    %global_sites_diff_hash = ();


    # extract site probabilities
    $seq_unmodified = '';
    $temp = $probSTY;
    $temp =~ s/[\(\)]/\|/g;
    @temp_array = split /\|/, $temp;
    $aa_index = 0;
    for ($i = 0; $i < @temp_array; $i++)
    {
        # it is a probability
        $prob = $temp_array[$i];
        if ($prob =~ /^[0-9\.]+$/)
        {
            # floor it, just in case
            # smallest I've observed in Maxquant output
            if ($prob < 0.001)
            {
                $prob = 0.001;
            }
            
            $global_sites_prob_hash{$aa_index} = $prob;
        }
        # it is sequence
        else
        {
            for ($j = 0; $j < length $prob; $j++)
            {
                $c = substr $prob, $j, 1;
                if ($c =~ /[A-Z]/)
                {
                    $aa_index++;
                    
                    $seq_unmodified .= $c;
                }
            }
        }
    }

    
    # fill in missing potential phosphosites using Score Diffs
    # assign probability of 0.001
    $temp = $diffSTY;
    $temp =~ s/[\(\)]/\|/g;
    @temp_array = split /\|/, $temp;
    $aa_index = 0;
    for ($i = 0; $i < @temp_array; $i++)
    {
        # it is a score
        $prob = $temp_array[$i];
        if ($prob =~ /^[-0-9\.]+$/)
        {
            if (!defined($global_sites_prob_hash{$aa_index}))
            {
                $global_sites_prob_hash{$aa_index} = 0.001;
            }
            
            $global_sites_diff_hash{$aa_index} = $prob;
        }
        # it is sequence
        else
        {
            for ($j = 0; $j < length $prob; $j++)
            {
                $c = substr $prob, $j, 1;
                if ($c =~ /[A-Z]/)
                {
                    $aa_index++;
                }
            }
        }
    }
    
    # number of potential STY in the sequence
    @temp_array = keys %global_sites_prob_hash;
    $num_potential_sites = @temp_array;

    # store the probability of the queried phosphosite
    $global_return_phosphosite_probability =
        $global_sites_prob_hash{$site};



    # MaxQuant got rid of the modified sequence column, and hosed numSTY, etc.
    # re-generate them using the probabilities
    $global_return_num_sites = 0;
    @temp_array = keys %global_sites_prob_hash;
    $num_potential_sites = @temp_array;
    for ($i = 0; $i < $num_potential_sites; $i++)
    {
        $prob = $global_sites_prob_hash{$temp_array[$i]};
        
        $global_return_num_sites += $prob;
    }
    $global_return_num_sites = floor($global_return_num_sites + 0.5);


    
    # populate reported phosphosite
    $global_populated_sites_hash{$site} = 1;
    
    # fill in remaining phosphosites
    @temp_array = sort cmp_phospho_probs keys %global_sites_prob_hash;
    for ($i = 0; $i < $global_return_num_sites; $i++)
    {
        $global_populated_sites_hash{$temp_array[$i]} = 1;
    }
    
    # generate new modified sequence
    # remember, all of the probability coordinates, etc. are base-1
    for ($i = 1; $i <= length $seq_unmodified; $i++)
    {
        $c = substr $seq_unmodified, $i-1, 1;
        
        if (defined($global_populated_sites_hash{$i}))
        {
            $seq_new .= $mod_type_char;
        }
        
        $seq_new .= $c;
    }


    # calculate cumulative probability
    $prob = 1.0;
    $aa_index = 0;
    for ($i = 0; $i < length $seq_new; $i++)
    {
        $c = substr $seq_new, $i, 1;
        if ($c =~ /[A-Z]/)
        {
            $aa_index++;
        }
        elsif ($c eq $mod_type_char)
        {
            $prob *= $global_sites_prob_hash{$aa_index + 1};
        }
    }
    
    # Regardless of individual probabilities, if global_return_num_sites == num_potential_sites,
    #  the probability is 100% that it is the correct assignment.
    # This should never be otherwise in the data, but I'll check anyways.
    if ($global_return_num_sites == $num_potential_sites)
    {
        $prob = 1.0;
    }
    
    # Floor the probability, since multiple low-confidence phosphosites
    # could result in a vanishingly small probability, although I have
    # not actually seen this occur.  Best to handle it, just in case.
    if ($prob < 0.001)
    {
        $prob = 0.001;
    }
    $global_return_cumulative_prob = sprintf "%.3f", $prob;


    # create probability signature

    # Emperically, the noise threshold is around 0.03 from sampling
    # So, floor data at 0.05, and later treat differences of <= 0.05 as equal
    #
    # Further empirical results, using 0.05, shows that it works very well,
    # and that *at least* 0.046 is required 
    #
    foreach $aa_index (keys %global_sites_prob_hash)
    {
        if ($global_sites_prob_hash{$aa_index} < 0.05)
        {
            $global_sites_prob_hash{$aa_index} = 0.05;
        }
    }
    
    $global_return_prob_rank = join ",",
               sort cmp_phospho_probs keys %global_sites_prob_hash;
    @temp_array = split ",", $global_return_prob_rank;
    $global_return_prob_rank = '';
    for ($i = 0; $i < @temp_array; $i++)
    {
        $aa_index  = $temp_array[$i];

        if ($i)
        {
            $aa_index2 = $temp_array[$i-1];
            
            # compare only if there are multiple alternative combinations
            if ($global_return_num_sites !=  $num_potential_sites &&
                ($global_sites_prob_hash{$aa_index2} -
                 $global_sites_prob_hash{$aa_index} > 0.05))
            {
                $global_return_prob_rank .= '>';
                
                # indicate when probabilities are quite different
                if ($global_sites_prob_hash{$aa_index2} -
                    $global_sites_prob_hash{$aa_index} >= 0.333)
                {
                    $global_return_prob_rank .= '>';
                }
            }
            else
            {
                $global_return_prob_rank .= ',';
            }
        }
        
        $global_return_prob_rank .= $aa_index;
    }
    
    # resort equalities to be sure that they are increasing site order
    @much_greater_array = split />>/, $global_return_prob_rank;
    for ($i = 0; $i < @much_greater_array; $i++)
    {
        @greater_array = split />/, $much_greater_array[$i];

        for ($j = 0; $j < @greater_array; $j++)
        {
           @equal_array = split /,/, $greater_array[$j];
           @equal_array = sort cmp_phospho_probs @equal_array;
           
           $greater_array[$j] = join ",", @equal_array;
        }
        
        $much_greater_array[$i] = join ">", @greater_array;
    }
    $global_return_prob_rank = join ">>", @much_greater_array;
    
    # Keep only the first global_return_num_sites+1 probs
    # report the +1 site as 'n'
    $temp = '';
    $aa_index = 0;
    $old_c = '';
    for ($i = 0; $i < length $global_return_prob_rank; $i++)
    {
        $c = substr $global_return_prob_rank, $i, 1;

        if ($c =~ /[^>,]/)
        {
            # don't increment if we're still in the same number
            if (!($i && (substr $global_return_prob_rank, $i-1, 1) =~ /[^>,]/))
            {
                $aa_index++;
            }
            
            if ($aa_index >= $global_return_num_sites + 1)
            {
                $c = 'n';
            }
            
            if ($old_c ne 'n')
            {
                $temp .= $c;
            }
        }
        else
        {
            if ($aa_index >= $global_return_num_sites + 1)
            {
                last;
            }

            $temp .= $c;
        }
        
        $old_c = $c;
    }
    # tack on a 'x' if there is no 'n' at the end
    # this will prevent Excel from potentially destroying the data on import
    if (!($temp =~ /n/))
    {
        $temp .= ',x';
    }

    $global_return_prob_rank = $temp;
    
    

    # fix protein positions
    $positions_prot =~ s/,/;/g;
    $positions_prot =~ s/\s+/;/g;
    $positions_prot =~ s/;+/;/g;
    $positions_prot =~ s/^;//;
    $positions_prot =~ s/;$//;
    @multi_protein_positions_array = split /;/, $positions_prot;
    for ($i = 0; $i < @multi_protein_positions_array; $i++)
    {
        $position_prot = $multi_protein_positions_array[$i];
        
        $protein_peptide_pos_diff = $position_prot - $site;
        
        @temp_array  = ($global_return_prob_rank =~ m/([^>,]+)/g);
        @temp_array2 = ($global_return_prob_rank =~ m/([>,]+)/g);
        
        $positions_prot_order = $temp_array[0] + $protein_peptide_pos_diff;
        
        for ($j = 0; $j < @temp_array2; $j++)
        {
            $positions_prot_order .= $temp_array2[$j];
            
            if ($temp_array[$j+1] =~ /[0-9]/)
            {
                $positions_prot_order .= $temp_array[$j+1] +
                                         $protein_peptide_pos_diff;
            }
            else
            {
                $positions_prot_order .= $temp_array[$j+1];
            }
        }
        
        $multi_protein_positions_array[$i] = $positions_prot_order;
        
#        printf "FOOBAR\t%s\t%s\t%s\t%s\n",
#                $site, $position_prot, $global_return_prob_rank,
#                $positions_prot_order;
    }
    $global_return_phosphosites_order = join ';', @multi_protein_positions_array;
#    printf "FOOBLAH\t%s\n", $global_return_positions_prot;
    
    if ($fix_flag)
    {
#        printf STDERR "Prob:    %.3f\n", $prob;
#       printf STDERR "Probs:   %s\n", $global_return_prob_rank;
    }
    
    return $seq_new;
}



#$pep_cutoff  = 0.05;
#score_cutoff = 75;

#$pep_cutoff   =  9999;
#$pep_cutoff   =  0.2;		# discard if > $pep_cutoff
$pep_cutoff  = 0.05;		# empirically much better than PEP < 0.20
$score_cutoff = -9999;
$mass_err_cutoff = 5;		# discard if > $mass_err_cutoff

$filename = shift;


if (!defined($filename))
{
    printf "Syntax: program.pl \"Phospho (STY) Sites.txt\"\n";
    exit (1);
}



# Read in Header
open INFILE, "$filename" or die "can't open $filename\n";

# skip down to first non-blank, non-comment line
while(defined($line=<INFILE>))
{
    $line =~ s/[\r\n]+//g;
    $flag = 1;

    if (!($line =~ /\S/)) { $flag = 0; }
    if ($line =~ /^#/) { $flag = 0; }

    if ($flag) { last; }
}

# header line
$line =~ s/[\r\n]+//;
$line =~ s/\"//g;
@array = split /\t/, $line;
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    $array[$i] =~ s/\s+/ /g;

    # damn it, Maxquant keeps changing case between versions...
    # this is *ROYALLY*screwing everything up
    #
    # I am just going to conform all the capitalization here
    # If a word is >= 3 long, and the first letter of a word isn't
    #  already capitalized, capitalize it
    #
    # 2021-07-22 Maxquant broke capitalization yet again, this time with
    #  Desthiobiotin vs. DesthioBiotin
    # HACK: conform DesthioBiotin to Desthiobiotin,
    #  since DesthioBiotin is poorly capitalized anyways
    #  they've spelled in wrong in the past too, so fix that typo as well
    #
    $array[$i] =~ s/DesthioBiotin/Desthiobiotin/ig;
    $array[$i] =~ s/Dethiobiotin/Desthiobiotin/ig;
    @word_array = split /\s+/, $array[$i];
    for ($j = 0; $j < @word_array; $j++)
    {
        $c = substr $word_array[$j], 0, 1;
        if (length $word_array[$j] >= 3 && $c =~ /[a-z]/)
        {
            $c = uc($c);
            $word_array[$j] =~ s/^([a-z])/$c/;
        }
    }
    $array[$i] = join ' ', @word_array;
    
    $header1_hash{$array[$i]} = $i;
    $header_array[$i] = $array[$i];
}
$header_line = join "\t", @array;
$num_cols_1 = @array;

# latest MaxQuant does *NOT* have a modified sequence window
# I'll need to parse the new formats and create one myself
# The "Sequence window" is bogus now, too, since it is too wide
# We'll have to derive our own from the Phospho (STY) Proabilities column
#
# Sequence window, Modification window, Peptide Window coverage
#  all have no bearing on the Phospho (STY) Probabilities, nor the number
#  or position of the modification.  They can even be missing!!
#
#$modified_sequence_col = $header1_hash{''};

#$charge_col            = $header1_hash{'Charge'};
$pep_col               = $header1_hash{'PEP'};
$score_col             = $header1_hash{'Score'};
$mass_err_col          = $header1_hash{'Mass Error [ppm]'};


# these go together
# more complete than Leading/Positions, butempty for REV_ accessions
$proteins_col          = $header1_hash{'Proteins'};
$positions_within_col  = $header1_hash{'Positions Within Proteins'};

# these go together
# has data for REV_ accessions, but can be missing accessions
$leading_prot_col      = $header1_hash{'Leading Proteins'};
$positions_prot_col    = $header1_hash{'Positions'};


$position_pep_col      = $header1_hash{'Position in Peptide'};
$window_col            = $header1_hash{'Sequence Window'};


# latest MaxQuant often is *blank* in this column, when it should have a number
# when it *does* have a number, the number is often wrong
# sometimes the number is different for the different accessions!
# well crap, now we can't even use that for figuring out how to parse the
# Phospho (STY) Proabilities correctly...
#
#$numSTY_col            = $header1_hash{'Number of Phospho (STY)'};

$probSTY_col           = $header1_hash{'Phospho (STY) Probabilities'};
$diffSTY_col           = $header1_hash{'Phospho (STY) Score Diffs'};

$reverse_col           = $header1_hash{'Reverse'};
$contaminant_col       = $header1_hash{'Contaminant'};

if (!defined($contaminant_col))
{
    $contaminant_col   = $header1_hash{'Potential Contaminant'};
}

# columns to be removed, since they can be huge and cause Excel to puke
if (defined($header1_hash{'Evidence IDs'}))
{
    $cols_to_remove_hash{$header1_hash{'Evidence IDs'}} = 'Evidence IDs';
}
if (defined($header1_hash{'MS/MS IDs'}))
{
    $cols_to_remove_hash{$header1_hash{'MS/MS IDs'}}    = 'MS/MS IDs';
}

# columns to be removed, since they clutter up the file with no information
foreach $header (sort keys %header1_hash)
{
    if ($header =~ /^Ratio mod\/base/i)
    {
        $cols_to_remove_hash{$header1_hash{$header}} = $header;
    }
}


$modification_type_char = 'p';	# phosphorylation

# might be some other modification
if (!defined($probSTY_col))
{
#    $numSTY_col            = $header1_hash{'Number of Oxidation (M)'};
    $probSTY_col           = $header1_hash{'Oxidation (M) Probabilities'};
    $diffSTY_col           = $header1_hash{'Oxidation (M) Score Diffs'};
    $modification_type_char     = 'o';	# oxidation
}
# might be some other modification
if (!defined($probSTY_col))
{
#    $numSTY_col            = $header1_hash{'Number of Desthiobiotin-ATP'};
    $probSTY_col           = $header1_hash{'Desthiobiotin-ATP Probabilities'};
    $diffSTY_col           = $header1_hash{'Desthiobiotin-ATP Score Diffs'};
    $modification_type_char     = 'd';	# desthiobiotin
}
# looks like they renamed the column headers for Desthiobiotin yet again...
if (!defined($probSTY_col))
{
    $probSTY_col           = $header1_hash{'ATP Probabilities'};
    $diffSTY_col           = $header1_hash{'ATP Score diffs'};
    $modification_type_char     = 'd';	# desthiobiotin
}
# or it could be a new typo that Maxquant introduced...
if (!defined($probSTY_col))
{
#    $numSTY_col            = $header1_hash{'Number of Dethiobiotin-ATP'};
    $probSTY_col           = $header1_hash{'Dethiobiotin-ATP Probabilities'};
    $diffSTY_col           = $header1_hash{'Dethiobiotin-ATP Score Diffs'};
    $modification_type_char     = 'd';	# desthiobiotin
}
# might be some other modification
if (!defined($probSTY_col))
{
#    $numSTY_col            = $header1_hash{'Number of Acetyl (K)'};
    $probSTY_col           = $header1_hash{'Acetyl (K) Probabilities'};
    $diffSTY_col           = $header1_hash{'Acetyl (K) Score Diffs'};
    $modification_type_char     = 'a';	# acetylation
}
# might be some other modification
if (!defined($probSTY_col))
{
#    $numSTY_col            = $header1_hash{'Number of GlyGly (K)'};
    $probSTY_col           = $header1_hash{'GlyGly (K) Probabilities'};
    $diffSTY_col           = $header1_hash{'GlyGly (K) Score Diffs'};
    $modification_type_char     = 'g';	# glygly
}
# might be some other modification
if (!defined($probSTY_col))
{
#    $numSTY_col            = $header1_hash{'Number of Biotin-HPDP'};
    $probSTY_col           = $header1_hash{'Biotin-HPDP Probabilities'};
    $diffSTY_col           = $header1_hash{'Biotin-HPDP Score Diffs'};
    $modification_type_char     = 'b';	# biotin
}

#if (!defined($modified_sequence_col))
#{
#    die "Sequence window column not found in file $filename\n";
#}
#if (!defined($charge_col))
#{
#    die "Charge column not found in file $filename\n";
#}
if (!defined($pep_col))
{
    die "PEP column not found in file $filename\n";
}
if (!defined($score_col))
{
    die "Score column not found in file $filename\n";
}
if (!defined($mass_err_col))
{
    die "Mass Error [ppm] column not found in file $filename\n";
}

if (!defined($positions_prot_col))
{
    die "Positions column not found in file $filename\n";
}
if (!defined($leading_prot_col))
{
    die "Leading Proteins column not found in file $filename\n";
}

#if (!defined($numSTY_col))
#{
#    die "Number of Phospho (STY) column not found in file $filename\n";
#}
if (!defined($position_pep_col))
{
    die "Position in peptide column not found in file $filename\n";
}
if (!defined($window_col))
{
    die "Sequence Window column not found in file $filename\n";
}
if (!defined($probSTY_col))
{
    die "Phospho (STY) Probabilities column not found in file $filename\n";
}
if (!defined($diffSTY_col))
{
    die "Phospho (STY) Score Diffs column not found in file $filename\n";
}
if (!defined($proteins_col))
{
    die "Proteins column not found in file $filename\n";
}
#if (!defined($positions_within_col))
#{
#    die "Positions Within Proteins column not found in file $filename\n";
#}

if (!defined($reverse_col))
{
    die "Reverse column not found in file $filename\n";
}
if (!defined($contaminant_col))
{
    die "Contaminant column not found in file $filename\n";
}


# scan for Intensity * columns
$has_intensities_flag = 0;
foreach $header (sort keys %header1_hash)
{
    if ($header =~ /^Intensity\s+[^\s]+/)
    {
        $has_intensities_flag = 1;

        $intensities_col_hash{$header1_hash{$header}} = 1;
    }
}
@intensities_col_array = sort {$a<=>$b} keys %intensities_col_hash;

# uh oh, there aren't any
# check to see if it is TMT data with >1 plex
$maxquant_tmt_multi_flag = 0;
if ($has_intensities_flag == 0)
{
    foreach $header (sort keys %header1_hash)
    {
        if ($header =~ /^Reporter intensity\s+[0-9]+\s+.*?$/i)
        {
            $has_intensities_flag = 1;

            $intensities_col_hash{$header1_hash{$header}} = 1;

            $maxquant_tmt_multi_flag = 1;
        }
    }
    @intensities_col_array = sort {$a<=>$b} keys %intensities_col_hash;
}


# uh oh, there aren't any
# must be a poorly formatted single TMT output
$maxquant_tmt_single_flag = 0;
if ($has_intensities_flag == 0)
{
    foreach $header (sort keys %header1_hash)
    {
        if ($header =~ /^Reporter intensity\s+[0-9]+/i)
        {
            $has_intensities_flag = 1;

            $intensities_col_hash{$header1_hash{$header}} = 1;

            $maxquant_tmt_single_flag = 1;
        }
    }
    @intensities_col_array = sort {$a<=>$b} keys %intensities_col_hash;
}


%seen_line_hash = ();


# remove any columns from the header that we don't want
@temp_array_i = split /\t/, $header_line;
@temp_array_j = ();
for ($i = 0, $j = 0; $i < @temp_array_i; $i++)
{
    if (!defined($cols_to_remove_hash{$i}))
    {
        $temp_array_j[$j] = $temp_array_i[$i];
        
        # fix potential maxquant single TMT header crap
        # conform them to be the same as >1 plex output
        if ($maxquant_tmt_single_flag &&
            defined($intensities_col_hash{$i}))
        {
            $temp_array_j[$j] .= ' TMT01';
        }
        
        $j++;
    }
}
$header_line_new = join "\t", @temp_array_j;

$header_line_new = sprintf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                            $header_line_new,
                            'Modification Type',
                            'Inferred Number of Sites',
                            'Inferred Modified Sequence',
                            'Modification Probability Order',
                            'Modification Site Probability Order',
                            'Modification Site Probability',
                            'Cumulative Modification Probability',
                            'Unmodified Sequence',
                            'CIC_Identifier',
                            'ModificationID';

printf "%s\n", $header_line_new;


# Read in file1 data
while(defined($line=<INFILE>))
{
    # skip blank lines
    if (!($line =~ /[A-Za-z0-9]/))
    {
        next;
    }

    $line =~ s/[\r\n]+//g;
    $line =~ s/\"//g;

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
    	
    	if (!($array[$i] =~ /\S/))
    	{
    	    $array[$i] = '---';
    	}
    }
    
    # check for missing fields at end and insert blanks
    for ($i = @array; $i < $num_cols_1; $i++)
    {
    	$array[$i] = '---';
    }

    # Clean up MaxQuant accession junk
    # Holy crap is this bad
    #
    # MaxQuant often neglects to carry over REV_ and CON_ from one column
    # to another, so we have to scan them all, then add them where needed
    %maxquant_con_hash = ();
    %maxquant_rev_hash = ();
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
            
#            $field = $array[$i];
            
            for ($j = 0; $j < @split_array; $j++)
            {
                # clean up sp| junk
                $split_array[$j] =~ s/sp\|([^\|]+)\|[^\|]+$/$1/;
                $split_array[$j] =~ s/sp\|([^\|]+)$/$1/;
                
                # clean up REFSEQ junk
                $split_array[$j] =~ s/REFSEQ://g;

                # clean up H-INV junk
                $split_array[$j] =~ s/H-INV://g;

                # clean up ENSEMBL junk
                $split_array[$j] =~ s/ENSEMBL://g;
                
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
    
#    $modified_sequence = $array[$modified_sequence_col];
#    $charge            = $array[$charge_col];

    $pep               = $array[$pep_col];
    $score             = $array[$score_col];
    $mass_err          = $array[$mass_err_col];

    $proteins          = $array[$proteins_col];
    if (defined($positions_within_col))
    {
        $positions_within = $array[$positions_within_col];
    }
    else
    {
        $positions_within = '';
    }
    
    $leading_prot      = $array[$leading_prot_col];
    $positions_prot    = $array[$positions_prot_col];

#    $numSTY            = $array[$numSTY_col];
    $position_pep      = $array[$position_pep_col];
    $window            = $array[$window_col];
    $probSTY           = $array[$probSTY_col];
    $diffSTY           = $array[$diffSTY_col];

    $reverse           = $array[$reverse_col];
    $contaminant       = $array[$contaminant_col];

    # skip reverse sequences and contaminants
    if ($reverse =~ /\+/)
    {
#        next;
    }
    if ($contaminant =~ /\+/)
    {
#        next;
    }
    
    # skip low confidence peptide assignments
    if (!($pep =~ /[0-9]/) || $pep > $pep_cutoff) { next; }
    if (!($score =~ /[0-9]/) || $score < $score_cutoff) { next; }
    if (!($mass_err =~ /[0-9]/) || abs($mass_err) > $mass_err_cutoff) { next; }

    # skip charges that aren't 2
#    if (!defined($charge) || $charge ne '2')
#    {
#        next;
#    }

    # skip lines with no ratios
#    if (!defined($ratio) || !($ratio =~ /[0-9]/))
#    {
#        next;
#    }

    # skip peptides with no intensities
    if ($has_intensities_flag)
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
    

    # count the number of accessions within each accession/position pair
    @temp_array       = split /;/, $proteins;
    $proteins_count   = 0;
    foreach $accession (@temp_array)
    {
        if ($accession =~ /[A-Za-z0-9]/)
        {
            $proteins_count++;
        }
    }
    @temp_array       = split /;/, $leading_prot;
    $leading_count    = 0;
    foreach $accession (@temp_array)
    {
        if ($accession =~ /[A-Za-z0-9]/)
        {
            $leading_count++;
        }
    }

    # pick one of the two paired sources of accessions/positions
    # choose Proteins / Positions Within Proteins first
    # use Leading Proteins / Positions if it has more
    #
    # Positions Within Proteins is missing from older MaxQuant output
    # Use Positions column when that is the only choice available
    #
    $proteins_chosen  = $proteins;
    $positions_chosen = $positions_within;
    if ($leading_count > $proteins_count ||
        $positions_within eq '')
    {
        $proteins_chosen  = $leading_prot;
        $positions_chosen = $positions_prot;
    }


    $modified_sequence = fix_phosphopeptide($modification_type_char,
                                            $position_pep,
                                            $probSTY, $diffSTY,
                                            $positions_chosen);

    $unmodified_sequence = $modified_sequence;
    $unmodified_sequence =~ s/[^A-Z]//g;

    $query_field       = sprintf "%s;%s;%s;%s;%s;%s",
                                 $modified_sequence,
                                 $proteins_chosen,
                                 $window,
                                 $position_pep,
                                 $global_return_num_sites,
                                 $global_return_prob_rank;

#    $query_field       = sprintf "%s;%s;%s;%s;%s;%s;%s",
#                                 $modified_sequence,
#                                 $proteins_chosen,
#                                 $window,
#                                 $position_pep,
#                                 $charge,
#                                 $global_return_num_sites,
#                                 $global_return_prob_rank;

    $accession_str = $proteins_chosen;
    $accession_str =~ s/\|/\;/g;
#    $accession_str =~ s/\:/\;/g;
    $accession_str =~ s/,/\;/g;
    $accession_str =~ s/\s+/\;/g;
    $accession_str =~ s/\;+/\;/g;
    $accession_str =~ s/^\;//;
    $accession_str =~ s/\;$//;
    
    @accession_array = split /\;/, $accession_str;
    %contam_accession_hash = ();
    %reversed_accession_hash = ();
    %local_accession_hash = ();

    
    # store reversed accessions for later removal
    # store positions for each accession
    @positions_array = split /\;/, $positions_chosen;
    for ($i = 0; $i < @accession_array; $i++)
    {
        $accession = $accession_array[$i];

        # remove FASTA junk
        # strangely, \b> does NOT match at the beginnig of the string !!
        #  so, I have to use negative lookbehind instead... perl bug ??
#        $accession =~ s/\b>//g;
        $accession =~ s/(?<![A-Za-z0-9])>//g;

        # strip -# from end of accession
        $accession =~ s/-\d+$//g;

        # uppercase accession for later lookup
        $accession =~ tr/a-z/A-Z/;
        
        # remove CON__
        $contam_flag = 0;
        if ($accession =~ /(?<![A-Za-z0-9])CON__/)
        {
            $accession =~ s/(?<![A-Za-z0-9])CON__//g;
            $contam_flag = 1;
        }

        # remove any leading/trailing whitespace
        $accession =~ s/^\s+//;
        $accession =~ s/\s+$//;

        # remove REV__
        $reverse_flag = 0;
        if ($accession =~ /(?<![A-Za-z0-9])REV__/)
        {
            $accession =~ s/(?<![A-Za-z0-9])REV__//g;
            $reverse_flag = 1;
        }
        
        # deal with :
        if (1 && $accession =~ /\:/)
        {
            @accession_array2 = split /\:/, $accession;
            
            foreach $accession2 (@accession_array2)
            {
                if ($contam_flag)
                {
                    $contam_accession_hash{$accession2} = 1;
                }
                if ($reverse_flag)
                {
                    $reversed_accession_hash{$accession2} = 1;
                }
                $local_accession_hash{$accession2} = $positions_array[$i];
            }
        }
        else
        {
            if ($contam_flag)
            {
                $contam_accession_hash{$accession} = 1;
            }
            if ($reverse_flag)
            {
                $reversed_accession_hash{$accession} = 1;
            }

            $local_accession_hash{$accession} = $positions_array[$i];
        }
    }

    # store new accessions
    @new_accession_array = ();
    @new_positions_array = ();
    $num_new_accessions = 0;
    foreach $accession (sort keys %local_accession_hash)
    {
        # skip reversed accessions
        if (defined($reversed_accession_hash{$accession}))
        {
#            printf STDERR "REVERSE\t%s\n", $accession;
            next;
        }

        $new_accession_array[$num_new_accessions] = $accession;
        $new_positions_array[$num_new_accessions] =
                             $local_accession_hash{$accession};

        $num_new_accessions++;
    }

    $proteins_new  = join ';', @new_accession_array;
    $positions_new = join ';', @new_positions_array;
    
    # no non-reversed accessions, include the reversed accessions back in...
    if ($num_new_accessions == 0)
    {
        foreach $accession (sort keys %local_accession_hash)
        {
            # must look up position here, since we add REV__ below
            $new_positions_array[$num_new_accessions] =
                                 $local_accession_hash{$accession};

            # rename reversed accessions
            if (defined($reversed_accession_hash{$accession}))
            {
#                printf STDERR "REVERSE\t%s\n", $accession;
                 $accession = 'REV__' . $accession;
            }

            $new_accession_array[$num_new_accessions] = $accession;
            
            $num_new_accessions++;
        }
    }

    $proteins_new  = join ';', @new_accession_array;
    $positions_new = join ';', @new_positions_array;

    # abbreviations for ModificationID
    $mod_type = $modification_type_char;
    $abbrev = 'xx';
    if ($mod_type eq 'p') { $abbrev = 'ph'; }
    if ($mod_type eq 'o') { $abbrev = 'ox'; }
    if ($mod_type eq 'd') { $abbrev = 'de'; }
    if ($mod_type eq 'g') { $abbrev = 'gl'; }
    if ($mod_type eq 'a') { $abbrev = 'ac'; }
    if ($mod_type eq 'b') { $abbrev = 'bh'; }   # Biotin-HPDP

    if    ($mod_type eq 'p') { $mod_type = 'Phosphorylation (STY)'; }
    elsif ($mod_type eq 'o') { $mod_type = 'Oxidation (M)'; }
    elsif ($mod_type eq 'd') { $mod_type = 'Desthiobiotin-ATP'; }
    elsif ($mod_type eq 'g') { $mod_type = 'GlyGly (K)'; }
    elsif ($mod_type eq 'a') { $mod_type = 'Acetyl (K)'; }
    elsif ($mod_type eq 'b') { $mod_type = 'Biotin-HPDP'; }
    else                     { $mod_type = 'Unknown' };

    $phosphosite_id = sprintf "%s_%s:%s",
        $abbrev, $proteins_new, $positions_new;
    
    # HACK -- heavy labeled peptides
#    if ($modified_sequence =~ /\(si\)/)
#    {
#        $phosphosite_id = sprintf "%s:Heavy", $phosphosite_id;
#        $modified_sequence =~ s/\(si\)//g,
#    }

    # remove any columns from the header that we don't want
    @temp_array_j = ();
    for ($i = 0, $j = 0; $i < @array; $i++)
    {
        if (!defined($cols_to_remove_hash{$i}))
        {
            $temp_array_j[$j++] = $array[$i];
        }
    }
    $line_new = join "\t", @temp_array_j;

    $line_print = sprintf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                                  $line_new,
                                  $mod_type,
                                  $global_return_num_sites,
                                  $modified_sequence,
                                  $global_return_prob_rank,
                                  $global_return_phosphosites_order,
                                  $global_return_phosphosite_probability,
                                  $global_return_cumulative_prob,
                                  $unmodified_sequence,
                                  $query_field,
                                  $phosphosite_id;
    

    if (!defined($seen_line_hash{$line_print}))
    {
        printf "%s\n", $line_print;
    }
    
    $seen_line_hash{$line_print} = 1;
}
close INFILE;
