#!/usr/bin/perl -w

# 2020-10-15  keep track of non- species of interest contaminants as well


# TODO --
#  add trypsin exclusion support
#  add unique + razor maxquant group inclusion/exclusion criteria

# can't pipe input, since we need to read through the file twice

sub bless_delimiter_bar_no_space
{
    my $text = $_[0];

    $text =~ s/\;/\|/g;
    $text =~ s/\/\//\|/g;
    $text =~ s/,/\|/g;
#    $text =~ s/\s+/\|/g;
    $text =~ s/\|+/\|/g;
    $text =~ s/^\|//;
    $text =~ s/\|$//;
    
    return $text;
}

sub bless_delimiter_bar
{
    my $text = $_[0];

    $text =~ s/\;/\|/g;
    $text =~ s/\/\//\|/g;
    $text =~ s/,/\|/g;
    $text =~ s/\s+/\|/g;
    $text =~ s/\|+/\|/g;
    $text =~ s/^\|//;
    $text =~ s/\|$//;
    
    return $text;
}

sub bless_delimiter_space
{
    my $text = $_[0];

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

    # not H-INV
    if (!($value_a =~ /\b_*HIT/) &&   $value_b =~ /\b_*HIT/)  {return -1;}
    if (  $value_a =~ /\b_*HIT/  && !($value_b =~ /\b_*HIT/)) {return 1;}
    if (!($value_a =~ /H\-INV/) &&   $value_b =~ /H\-INV/)  {return -1;}
    if (  $value_a =~ /H\-INV/  && !($value_b =~ /H\-INV/)) {return 1;}

    # uniprot first
    if ($value_a =~ /_/ && !($value_b =~ /_/)) {return -1;}
    if ($value_b =~ /_/ && !($value_a =~ /_/)) {return  1;}

    # accession, alphabetical
    return ($value_a cmp $value_b);
}


$i = 0;
$annotation_headers_to_keep[$i++] = 'Species';
$annotation_headers_to_keep[$i++] = 'Species_flag';
$annotation_headers_to_keep[$i++] = 'Species_contaminant_flag';
$annotation_headers_to_keep[$i++] = 'Other_species_flag';
$annotation_headers_to_keep[$i++] = 'Unknown_species_flag';
$annotation_headers_to_keep[$i++] = 'Non-Species_contaminant_flag';
$annotation_headers_to_keep[$i++] = 'Reverse_flag';
$annotation_headers_to_keep[$i++] = 'Trypsin_flag';
$annotation_headers_to_keep[$i++] = 'Accession_Condensed';
#$annotation_headers_to_keep[$i++] = 'GeneID';
#$annotation_headers_to_keep[$i++] = 'NumGenes';
#$annotation_headers_to_keep[$i++] = 'WarningLevel';
#$annotation_headers_to_keep[$i++] = 'Symbol';
#$annotation_headers_to_keep[$i++] = 'Alias';
#$annotation_headers_to_keep[$i++] = 'Type';
#$annotation_headers_to_keep[$i++] = 'Location';
#$annotation_headers_to_keep[$i++] = 'Description';


# read in command line arguments
$num_files = 0;
$syntax_error_flag = 0;
$use_maxquant_unique_flag = 0;
$use_maxquant_razor_flag  = 0;
$filter_flag = 0;

for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        if ($field eq '--filter')
        {
            $filter_flag = 1;
        }
        elsif ($field eq '--regrouped' ||
               $field eq '--regroup')
        {
            $use_maxquant_razor_flag  = 0;
            $use_maxquant_unique_flag = 0;
        }
        elsif ($field eq '--maxquant-razor')
        {
            $use_maxquant_razor_flag = 1;
        }
        elsif ($field eq '--maxquant-unique')
        {
            $use_maxquant_unique_flag = 1;
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
    }
}

if ($syntax_error_flag)
{
    exit(1);
}

open INFILE, "$filename" or die "can't open file $filename\n";

$line = <INFILE>;
$line =~ s/[\r\n]+//g;
@array = split /\t/, $line;
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    $array[$i] =~ s/\"//g;
    
    $header_col_hash{$array[$i]} = $i;
}

# new headers from new maxquant pipeline
$sequence_col        = $header_col_hash{'Sequence'};
$accession_short_col = $header_col_hash{'Accession_Condensed'};
$accession_full_col  = $header_col_hash{'Accession_Full'};
$trypsin_col         = $header_col_hash{'Trypsin_flag'};
$regroup_col         = $header_col_hash{'Group_Regrouped'};

if ($use_maxquant_unique_flag)
{
    $group_col       = $header_col_hash{'Group_Maxquant'};
}
elsif ($use_maxquant_razor_flag)
{
    $group_col       = $header_col_hash{'Group_Maxquant_Razor'};
}
else
{
    $group_col       = $header_col_hash{'Group_Regrouped'};
}

if (!defined($sequence_col))
{
    print STDERR "ABORT -- missing sequence column\n";
    exit(1);
}
if (!defined($accession_short_col))
{
    print STDERR "ABORT -- missing condensed accession column\n";
    exit(1);
}
if (!defined($accession_full_col))
{
    print STDERR "ABORT -- missing full accession column\n";
    exit(1);
}
if (!defined($regroup_col))
{
    print STDERR "ABORT -- missing regrouped group column\n";
    exit(1);
}
if (!defined($group_col))
{
    print STDERR "ABORT -- missing group column\n";
    exit(1);
}
if (!defined($trypsin_col))
{
    print STDERR "ABORT -- missing trypsin column\n";
    exit(1);
}


# older headers from IDP pipeline
$geneid_col    = $header_col_hash{'GeneID'};
$numgenes_col  = $header_col_hash{'NumGenes'};
$accession_col = $header_col_hash{'Accession_Full'};

if (!defined($accession_col))
{
    $accession_col = $header_col_hash{'Accession_Protein'};
}
if (!defined($accession_col))
{
    $accession_col = $header_col_hash{'Accession'};
}
if (!defined($accession_col))
{
    $accession_col = $header_col_hash{'Accession_RNA'};
}
if (!defined($accession_col))
{
    $accession_col = $header_col_hash{'Accession_Genomic'};
}


if (!defined($geneid_col))
{
    print STDERR "ABORT -- missing GeneID column\n";
    exit(1);
}
if (!defined($numgenes_col))
{
    print STDERR "ABORT -- missing NumGenes column\n";
    exit(1);
}
if (!defined($accession_col))
{
    print STDERR "ABORT -- missing accession column\n";
    exit(1);
}



$global_worst_n  = -9E9;
$global_worst_n2 = -9E9;

while(defined($line=<INFILE>))
{
    $line =~ s/[\r\n]+//g;

    @array = split /\t/, $line;
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\"//g;
    }

    $sequence        = $array[$sequence_col];
#    $accession_short = $array[$accession_short_col];
#    $accession_full  = $array[$accession_full_col ];
    $group           = $array[$group_col];
    $regroup         = $array[$regroup_col];

    # strip g: from beginning, since code assumes it won't be there
    if ($use_maxquant_razor_flag  == 0 &&
        $use_maxquant_unique_flag == 0)
    {
        $group =~ s/^g://;
    }
    $regroup =~ s/^g://;

    # for now, use the peptide sequence as the probeid
    # will replace with spectra identifier later once we add spectra support
    $probeid         = $sequence;

#    $geneid_str      = $array[$geneid_col];
    $numgenes        = $array[$numgenes_col];
    $accession_str   = $array[$accession_col];

    # create a geneid string out of the group name
    $geneid_str = $regroup;
    $geneid_str = bless_delimiter_space($geneid_str);

    $accession_str = bless_delimiter_space($accession_str);
    @temp_accession_array = split / /, $accession_str;

    %temp_geneid_hash = ();
    @temp_array = split / /, $geneid_str;
    foreach $geneid (@temp_array)
    {
        if ($geneid =~ /[A-Za-z0-9]/)
        {
            $temp_geneid_hash{$geneid} = 1;
        }
    }

    
    $n  = @temp_array;
    $n2 = $numgenes;

    if ($n > $global_worst_n)
    {
        $global_worst_n = $n;
    }
    if ($n2 > $global_worst_n2)
    {
        $global_worst_n2 = $n;
    }

    foreach $geneid (@temp_array)
    {
        $geneid_probeid_hash{$geneid}{$probeid} = 1;
    
        if (!defined($geneid_hash{$geneid}))
        {
            $geneid_hash{$geneid}{best_n}  = 9E9;
            $geneid_hash{$geneid}{best_n2} = 9E9;
        }
        
        if ($n2 <= $geneid_hash{$geneid}{best_n2})
        {
            $geneid_hash{$geneid}{best_n2} = $n2;
        
            if ($n < $geneid_hash{$geneid}{best_n})
            {
                $geneid_hash{$geneid}{best_n} = $n;
            }
        }
        
        # store spectra identifiers for each gene
        $geneid_hash{$geneid}{groups}{$group} = 1;

        $group_hash{$group}{geneids}{$geneid} = 1;
    }

    $group_hash{$group}{n}  = $n;
    $group_hash{$group}{n2} = $n2;
    $group_hash{$group}{probeids}{$probeid} = 1;
    
    # store group accessions
    foreach $accession (@temp_accession_array)
    {
        $group_hash{$group}{accessions}{$accession} = 1;
    }
    
    # store annotation line for later
    $probeid_annotation_line_hash{$probeid} = $line;
}

@group_array = sort keys %group_hash;
@geneid_array = sort keys %geneid_hash;


# count number of spectra for each geneid
foreach $geneid (@geneid_array)
{
    @temp_array = keys %{$geneid_probeid_hash{$geneid}};
    
    $geneid_spectra_counts_all{$geneid} = @temp_array;
}


# store maximum number of members in a group for each geneid
foreach $group (@group_array)
{
    @temp_array = keys %{$group_hash{$group}{geneids}};
    $num_members = @temp_array;
    
    @temp_array = keys %{$group_hash{$group}{geneids}};
    foreach $geneid (@temp_array)
    {
        if (!defined($geneid_max_group_membership{$geneid}))
        {
            $geneid_max_group_membership{$geneid} = $num_members;
        }
        elsif ($num_members > $geneid_max_group_membership{$geneid})
        {
            $geneid_max_group_membership{$geneid} = $num_members;
        }
    }
}


# score groups
foreach $group (@group_array)
{
    @temp_array = keys %{$group_hash{$group}{probeids}};
    $num_members = @temp_array;
    
    $group_hash{$group}{num_members} = $num_members;
    
    $n = $group_hash{$group}{n};
    $n2 = $group_hash{$group}{n2};

    @temp_array = keys %{$group_hash{$group}{accessions}};
    
    # hack for unannotated accessions
    if ($n < 1 || $n2 < 1)
    {
        $num_accessions = @temp_array;

        if ($n < 1)
        {
#            $n = $global_worst_n + $num_accessions;
            $n = $num_accessions;
        }
        if ($n2 < 1)
        {
#            $n2 = $global_worst_n2 + $num_accessions;
            $n2 = $num_accessions;
        }
    }
    
    $group_hash{$group}{score} = ($num_members * $num_members) / $n;
}


%geneid_assigned_hash = ();
%group_kept_hash      = ();


# keep all groups with exactly one gene
foreach $group (@group_array)
{
    if ($group_hash{$group}{n} == 1)
    {
        @temp_array = keys %{$group_hash{$group}{geneids}};

        # skip truly one-hit wonders, not represented in any other group
        $geneid = $temp_array[0];
        $num_spectra = $geneid_spectra_counts_all{$geneid};
        if ($num_spectra < 2)
        {
            $one_hit_wonder_hash{$group} = 1;
            
#            printf STDERR "One-hit wonder:\t%s\n", $geneid;
#            next;
        }
        
        # skip geneids with no group that has >=2 members
##        if ($geneid_max_group_membership{$geneid} < 2)
##        {
##            printf STDERR "Poor membership:\t%s\n", $geneid;
##            next;
##        }

        $group_kept_hash{$group} = 1;

        foreach $geneid (@temp_array)
        {
            $geneid_assigned_hash{$geneid} = 1;
        }
    }
}


# add in 2-gene groups, for those not already assigned
# this is to be sure that we include gene1 + gene1-gene2 readthroughs
# however, 1-peptide multi-gene groupings are highly suspect
#  in that there usually isn't any obvious functional or symbol pattern
#  so, we're only going to add in those with >= 2 peptides here
foreach $geneid (@geneid_array)
{
    if (!defined($geneid_assigned_hash{$geneid}))
    {
        @temp_array = keys %{$geneid_hash{$geneid}{groups}};

        $best_score = -9E99;
        foreach $group (@temp_array)
        {
            if ($group_hash{$group}{n} == 2 &&
                $group_hash{$group}{num_members} >= 2)
            {
                $score = $group_hash{$group}{score};

                if ($score > $best_score)
                {
                    $best_score = $score;
                }
            }
        }
        
        if ($best_score > 0)
        {
            foreach $group (@temp_array)
            {
                if ($group_hash{$group}{n} == 2 &&
                    $group_hash{$group}{num_members} >= 2 &&
                    $group_hash{$group}{score} == $best_score)
                {
                    $group_kept_hash{$group} = 1;
                }
            }
        }
    }
}


# re-flag kept geneids
foreach $group (keys %group_kept_hash)
{
    @temp_array = keys %{$group_hash{$group}{geneids}};

    foreach $geneid (@temp_array)
    {
        $geneid_assigned_hash{$geneid} = 1;
    }
}


# add in groups with the best score for each gene
foreach $geneid (@geneid_array)
{
    # skip truly one-hit wonders, not represented in any other group
    if ($geneid_spectra_counts_all{$geneid} < 2)
    {
#        next;
    }

    # skip geneids with no group that has >=2 members
#    if ($geneid_max_group_membership{$geneid} < 2)
#    {
#        printf STDERR "Poor membership:\t%s\n", $geneid;
#        next;
#    }

    if (1 || !defined($geneid_assigned_hash{$geneid}))
    {
        @temp_array = keys %{$geneid_hash{$geneid}{groups}};
        
        $best_score = -9E99;
        foreach $group (@temp_array)
        {
            $score = $group_hash{$group}{score};

            if ($score > $best_score)
            {
                $best_score = $score;
            }
        }

        if ($best_score > 0)
        {
            foreach $group (@temp_array)
            {
                if ($group_hash{$group}{score} == $best_score)
                {
                    $group_kept_hash{$group} = 1;
                }
            }
        }
    }
}

# re-flag kept geneids
foreach $group (keys %group_kept_hash)
{
    @temp_array = keys %{$group_hash{$group}{geneids}};

    foreach $geneid (@temp_array)
    {
        $geneid_assigned_hash{$geneid} = 1;
    }
}

if (0)
{
  foreach $geneid (@geneid_array)
  {
    if (!defined($geneid_assigned_hash{$geneid}))
    {
        printf STDERR "UNASSIGNED\t$geneid\n";
    }
  }
}


@group_kept_array = sort keys %group_kept_hash;
%seen_probeid_hash = ();


# assign group numbers to geneid combinations
#for ($i = 0; $i < @group_kept_array; $i++)
#{
#    $group = $group_kept_array[$i];
#    $group_renumbered_hash{$group} = $i + 1;
#}


# assemble group annotation stuff
%group_annotation_hash = ();
%seen_annotation_hash = ();
foreach $group (@group_array)
{
    @temp_array = keys %{$group_hash{$group}{probeids}};
    foreach $probeid (@temp_array)
    {
        $line = $probeid_annotation_line_hash{$probeid};
        
        @array = split /\t/, $line;
        
        for ($i = 0; $i < @annotation_headers_to_keep; $i++)
        {
            $field = '';
            $header = $annotation_headers_to_keep[$i];
            $col = $header_col_hash{$header};
            
            if (defined($col))
            {
                $field = $array[$col];
                
                if (!defined($field))
                {
                    $field = '';
                }
            }
            
            if ($field ne '')
            {
                $field = bless_delimiter_bar_no_space($field);
                @temp_array2 = split /\|/, $field;
                
                foreach $value (@temp_array2)
                {
                    $group_annotation_hash{$group}{$header}{$value} = 1;
                    $seen_annotation_hash{$header} = 1;
                }
            }
        }
    }
}


# output header line
printf "%s",   'ProteinGroup';
printf "\t%s", 'NumPeptidesPreNorm';
printf "\t%s", 'OneHitWonder';
for ($i = 0; $i < @annotation_headers_to_keep; $i++)
{
    $header = $annotation_headers_to_keep[$i];
    if (defined($seen_annotation_hash{$header}))
    {
        printf "\t$header";
    }
}
printf "\t%s", 'Accession';
print "\n";


# print *all* groups if we told it not to filter them
@groups_to_print_array = @group_kept_array;
if ($filter_flag == 0)
{
    @groups_to_print_array = @group_array;
}


foreach $group (@groups_to_print_array)
{
    $one_hit_wonder_flag = 0;
    if (defined($one_hit_wonder_hash{$group}))
    {
        $one_hit_wonder_flag = 1;
    }

    @temp_array = keys %{$group_hash{$group}{probeids}};
    foreach $probeid (@temp_array)
    {
        $seen_probeid_hash{$probeid} = 1;
#        $probeid_group_hash{$probeid}{$group} = 1;
    }
    
    $num_peptides = @temp_array;

    @temp_array = sort compare_accession keys
        %{$group_hash{$group}{accessions}};
    $accession_str = join ' ', @temp_array;
    
    if ($use_maxquant_razor_flag  == 0 &&
        $use_maxquant_unique_flag == 0)
    {
        print 'g:';
    }
    # skip multi-group maxquant groups, since maxquant does not keep them
    else
    {
        if ($group =~ /;/)
        {
            next;
        }
    }
    
    print $group;

    printf "\t%s", $num_peptides;
    printf "\t%s", $one_hit_wonder_flag;

    # print more annotation
    for ($i = 0; $i < @annotation_headers_to_keep; $i++)
    {
        $header = $annotation_headers_to_keep[$i];
        if (defined($seen_annotation_hash{$header}))
        {
            $value = '';
            
            if (defined($group_annotation_hash{$group}) &&
                defined($group_annotation_hash{$group}{$header}))
            {
                @temp_array =
                    sort keys %{$group_annotation_hash{$group}{$header}};
                
                $value = join '; ', @temp_array;
            }
            
            # replace multiple flags with single flag
            if ($header =~ /_flag/)
            {
                if ($value eq '0; 1')
                {
                    $value = 1;
                }
            }
            
            printf "\t$value";
        }
    }

    printf "\t%s", $accession_str;
    
    print "\n";
}
