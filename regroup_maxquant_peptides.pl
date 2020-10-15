#!/usr/bin/perl -w


# TODO -- write a separate script to keep only the best groups per geneid

# 2020-10-15  fixed bug where unknown species were annotated as known
#             keep track of non- species of interest contaminants as well

sub cmp_razor
{
    my $razor_a    = $maxquant_group_hash{$a}{'razor'};
    my $score_a    = $maxquant_group_hash{$a}{'score'};
    my $unique_a   = $maxquant_group_hash{$a}{'unique'};
    my $peptides_a = $maxquant_group_hash{$a}{'peptides'};

    my $razor_b    = $maxquant_group_hash{$b}{'razor'};
    my $score_b    = $maxquant_group_hash{$b}{'score'};
    my $unique_b   = $maxquant_group_hash{$b}{'unique'};
    my $peptides_b = $maxquant_group_hash{$b}{'peptides'};

    # razor + unique
    if ($razor_a > $razor_b) { return -1; }
    if ($razor_a < $razor_b) { return  1; }

    # score
    if ($score_a > $score_b) { return -1; }
    if ($score_a < $score_b) { return  1; }

    # this number can be larger than razor + unique
    # I do not know what this means, if the missing peptides are crap or not
    # so far, this has not made a difference in best scoring razor group
    if ($peptides_a > $peptides_b) { return -1; }
    if ($peptides_a < $peptides_b) { return  1; }

#    if ($unique_a > $unique_b) { return -1; }
#    if ($unique_a < $unique_b) { return  1; }
    
    # only need to print when debugging, or if I ever encounter
    #  final razor group selections that are tied
    # printf STDERR "TIED_POTENTIAL_RAZOR_GROUPS:\t%s\t%s\n", $a, $b;
    
    return ($a<=>$b);
}


sub cmp_accessions
{
    my $accession_a = $a;
    my $accession_b = $b;

    # best refseqs
    if (  $accession_a =~ /^NM_/  && !($accession_b =~ /^NM_/)) {return -1;}
    if (!($accession_a =~ /^NM_/) &&   $accession_b =~ /^NM_/)  {return 1;}
    if (  $accession_a =~ /^NR_/  && !($accession_b =~ /^NR_/)) {return -1;}
    if (!($accession_a =~ /^NR_/) &&   $accession_b =~ /^NR_/)  {return 1;}
    if (  $accession_a =~ /^NP_/  && !($accession_b =~ /^NP_/)) {return -1;}
    if (!($accession_a =~ /^NP_/) &&   $accession_b =~ /^NP_/)  {return 1;}
    
    # uniprot first
    if ($accession_a =~ /_/ && !($accession_b =~ /_/)) {return -1;}
    if ($accession_b =~ /_/ && !($accession_a =~ /_/)) {return  1;}

    # accession, alphabetical
    return ($accession_a cmp $accession_b);
}


# we need to exclude any peptides that hit both human and digest trypsin
$trypsin_hash{TRYP_PIG}   = 'P00761';
$trypsin_hash{P00761}     = 'TRYP_PIG';
$trypsin_hash{TRY1_BOVIN} = 'P00760';
$trypsin_hash{P00760}     = 'TRY1_BOVIN';

# ENSBTAP00000038253 is E1B991_BOVIN E1B991_BOVIN KRT2 Keratin 2 (Bos taurus)
# not sure what I want to do with this, it is junk from extra search file

$annotation_filename  = shift;   # ipi_uniprot_annotation_human.txt
$mq_proteins_filename = shift;   # protein level, has off-species annotation
$mq_peptides_filename = shift;   # peptides.txt


open ANNOTATION, "$annotation_filename" or die "can't open $annotation_filename\n";
open PROTEINS, "$mq_proteins_filename" or die "can't open $mq_proteins_filename\n";
open PEPTIDES, "$mq_peptides_filename" or die "can't open $mq_peptides_filename\n";




# header line
$line = <ANNOTATION>;
$line =~ s/[\r\n]+//g;
@array = split /\t/, $line;
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    $array[$i] =~ s/\s+/ /g;
    
    $field = $array[$i];
    $annotation_header_col_hash{$field} = $i;
#    $annotation_header_col_array[$i] = $field;
}

$annotation_accession_col = $annotation_header_col_hash{'Accession'};
$annotation_geneid_col    = $annotation_header_col_hash{'GeneID'};
$annotation_prot_col      = $annotation_header_col_hash{'Accession_Protein'};

if (!defined($annotation_accession_col))
{
    die "ABORT -- cannot find Accession column in annotation file\n";
}
if (!defined($annotation_geneid_col))
{
    die "ABORT -- cannot find GeneID column in annotation file\n";
}
if (!defined($annotation_prot_col))
{
    die "ABORT -- cannot find Accession_Protein column in annotation file\n";
}


# map proteins to geneids, map uniprot to swissprot when possible
while(defined($line=<ANNOTATION>))
{
    $line =~ s/[\r\n]+//g;

    @array = split /\t/, $line;

    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\s+/ /g;
    }
    
    $accession = $array[$annotation_accession_col];
    $geneid    = $array[$annotation_geneid_col];
    $prot      = $array[$annotation_prot_col];

    $stripped  = $accession;
    $stripped  =~ s/(\-|\.)\d+$//;
    
    if ($geneid =~ /[0-9]/)
    {
        $accession_geneid_hash{$accession} = $geneid;
        $accession_geneid_hash{$stripped}  = $geneid;
    }

    # map uniprot <-> swissprot, where possible
    if (!($accession =~ /[^_]+_[^_]+/) &&
        $prot =~ /\b[^_]{3,}_[^_]{3,}\b/)
    {
        @prot_array = split /\s+/, $prot;
        
        foreach $temp_accession (@prot_array)
        {
            if ($temp_accession =~ /^([^_]{3,})_[^_]{3,}$/)
            {
                $uniprot   = $accession;
                $swissprot = $temp_accession;
            
                $accession_root = $1;
                
                if (1 || $accession_root ne $accession)
                {
                    $uniprot_to_sp_hash{$uniprot} = $swissprot;
                    $sp_to_uniprot_hash{$swissprot} = $uniprot;
                }
                
                # assume that the annotation file contains only
                #  the species of interest
                $species_accession_hash{$uniprot}   = 1;
                $species_accession_hash{$swissprot} = 1;
            }
        }
    }
    
    # assume that the annotation file contains only
    #  the species of interest
    $species_accession_hash{$accession}   = 1;
}
close ANNOTATION;
printf STDERR "Finished reading in annotation file\n";


# header line
$line = <PROTEINS>;
if (defined($line))
{
    $line =~ s/[\r\n]+//g;
    @array = split /\t/, $line;
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\s+/ /g;
    
        $field = $array[$i];
        $prot_header_col_hash{$field} = $i;
    }
    #$header_line = join "\t", @array;

}


$prot_fa_col       = $prot_header_col_hash{'Fasta headers'};
$prot_peptides_col = $prot_header_col_hash{'Peptides'};
$prot_razor_col    = $prot_header_col_hash{'Razor + unique peptides'};
$prot_unique_col   = $prot_header_col_hash{'Unique peptides'};
$prot_group_col    = $prot_header_col_hash{'id'};
$prot_score_col    = $prot_header_col_hash{'Score'};

if (!defined($prot_fa_col))
{
   die "ABORT -- can't find Fasta headers column in protein file\n";
}
if (!defined($prot_peptides_col))
{
   die "ABORT -- can't find Peptides column in protein file\n";
}
if (!defined($prot_razor_col))
{
   die "ABORT -- can't find Razor + unique peptides column in protein file\n";
}
if (!defined($prot_unique_col))
{
   die "ABORT -- can't find Unique peptides column in protein file\n";
}
if (!defined($prot_group_col))
{
   die "ABORT -- can't find id column in protein file\n";
}
if (!defined($prot_score_col))
{
   die "ABORT -- can't find Score column in protein file\n";
}


while(defined($line=<PROTEINS>))
{
    # not a useful proteins file, skip it with a warning
    if (!defined($prot_fa_col))
    {
        print STDERR "WARNING -- can't find Fasta headers line in protein file\n";
        last;
    }

    $line =~ s/[\r\n]+//g;

    @array = split /\t/, $line;

    for ($i = 0; $i < @array; $i++)
    {
        $array_len[$i] = length $array[$i];
    
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\s+/ /g;
    }

    $fasta    = $array[$prot_fa_col];
    $peptides = $array[$prot_peptides_col];
    $razor    = $array[$prot_razor_col];
    $unique   = $array[$prot_unique_col];
    $group    = $array[$prot_group_col];
    $score    = $array[$prot_score_col];
    
    $maxquant_group_hash{$group}{'peptides'} = $peptides;
    $maxquant_group_hash{$group}{'razor'}    = $razor;
    $maxquant_group_hash{$group}{'unique'}   = $unique;
    $maxquant_group_hash{$group}{'score'}    = $score;
    
    # parse out additional annotation
    $fasta =~ s/^>//;
    @array2 = split />/, $fasta;
    for ($i = 0; $i < @array2; $i++)
    {
        $field = $array2[$i];
        $field =~ s/[ ;]+$//g;
        
        $accession = '';
        $uniprot   = '';
        $swissprot = '';

        $desc_str  = '';
        $symbol    = '';
        $species   = '';

        # sp|
        if ($field =~ /sp\|([^| ]+)\|([^| ]+)\s+(.*)/)
        {
            $uniprot   = $1;
            $swissprot = $2;
            $desc_str  = $3;

            $stripped  = $uniprot;
            $stripped  =~ s/(\-|\.)\d+$//;

            # extract gene symbol
            if ($desc_str =~ /GN=(.*)/)
            {
                $symbol = $1;

                # strip remaining fields
                $symbol =~ s/\s+[A-Za-z]+\=.*//;
            }

            # extract species
            if ($desc_str =~ /OS=(.*)/)
            {
                $species = $1;

                # strip remaining fields
                $species =~ s/\s+[A-Za-z]+\=.*//;
            }
            
            # check for likely truncated fields, set them to blank
            # poor use of database with fixed width fields...
            if ($i == @array2 - 1 &&
                $array_len[$prot_fa_col] == 256)
            {
                if (!($desc_str =~ /[A-Z]\=\S/))
                {
                    $desc_str = '';
                }
                if (!($desc_str =~ /GN\=\S/))
                {
                    $species = '';
                }
                if (!($desc_str =~ /PE\=\S/))
                {
                    $symbol = '';
                }
            }

            # strip remaining fields
            $desc_str =~ s/\s*[A-Za-z]+\=.*//;
            
            # infer species from swissprot
            if ($species eq '' && $swissprot =~ /_([A-Za-z]+)/)
            {
                $species = $1;
                $species =~ s/HUMAN/Homo sapiens/;
                $species =~ s/MOUSE/Mus musculus/;
                $species =~ s/BOVIN/Bos taurus/;
            }
            
            # swissprot accession isn't actually swissprot
            if ($uniprot  eq $swissprot ||
                $stripped eq $swissprot)
            {
                $swissprot = '';
            }
            # already obsolete uniprot entry, replaced by another
            elsif (!($swissprot =~ /_[A-Z]$/))
            {
               # TODO
            }
            
#            printf "%s\t%s\t%s\t%s\t%s\n",
#                $uniprot,
#                $swissprot,
#                $species,
#                $symbol,
#                $desc_str;
        }
        # ENSEMBL: or REFSEQ:
        # there can be an accession that preceeds SWISS-PROT: or TREMBL:
        # some :accession;accession can happen
        elsif ($field =~ /\b[A-Za-z]+:(\S+)/)
        {
            $accession = $1;
            
            if ($field =~ /\b[A-Za-z]+:\S+\s+(.*)/)
            {
                $desc_str = $1;
                
                # (species)
                if ($desc_str =~ /\s*\(([A-Za-z ]+)\)\s*/)
                {
                    $species = $1;
                    
                    # check for "- Bos taurus (Bovine)", etc.
                    if ($desc_str =~ /\s*-*\s*Homo sapiens \(Human\)/i)
                    {
                        $species = 'Homo sapiens';
                        $desc_str =~ s/\s*-*\s*Homo sapiens \(Human\)\s*//i;
                    }
                    elsif ($desc_str =~ /\s*-*\s*Bos taurus \(Bovine\)/i)
                    {
                        $species = 'Bos taurus';
                        $desc_str =~ s/\s*-*\s*Bos taurus \(Bovine\)\s*//i;
                    }
                    elsif ($desc_str =~ /\s*-*\s*Mus musculus \(Mouse\)/i)
                    {
                        $species = 'Mus musculus';
                        $desc_str =~ s/\s*-*\s*Mus musculus \(Mouse\)\s*//i;
                    }
                    else
                    {
                        # strip the species
                        $desc_str =~ s/^\s*\(([A-Za-z ]+)\)\s*//;
                    }
                }
                # other format
                elsif ($desc_str =~ /[A-Za-z_]\=\S+\s/)
                {
                    # taxid
                    if ($desc_str =~ /Tax_Id\=(\S+)\s/i)
                    {
                        $species = $1;
                        
                        if ($species == 9606)
                        {
                            $species = 'Homo sapiens';
                        }
                        elsif ($species == 10090)
                        {
                            $species = 'Mus musculus';
                        }
                        elsif ($species == 9913)
                        {
                            $species = 'Bos taurus';
                        }
                    }
                    # symbol, assume it doesn't have a space in it...
                    # this is NOT an entirely safe assumption, but due to
                    #  the crappy formatting, we are forced to make it here,
                    #  which may lead to corrupt gene symbols and descriptions
                    if ($desc_str =~ /Gene_Symbol\=(\S+)\s/i)
                    {
                        $symbol = $1;
                    }
                }
            }

            # check for likely truncated fields, set them to blank
            # poor use of database with fixed width fields...
            if ($i == @array2 - 1 &&
                $array_len[$prot_fa_col] == 256)
            {
                if (!($desc_str =~ /Gene_Symbol\=/))
                {
                    $species = '';
                }

                $desc_str = '';
            }

            # strip key=value pairs
            $desc_str =~ s/\s*[A-Za-z_]+\=\S+\s*//g;
            
#            printf "%s\t%s\t%s\t%s\n",
#                $accession, $species, $symbol, $desc_str;
        }
        # unknown format
        elsif ($field =~ /(\S+)\s+(.*)/)
        {
            $accession = $1;
            $desc_str  = $2;
        }

        # store uniprot:sp mappings
        if ($uniprot   =~ /[A-Za-z0-9]/ &&
            $swissprot =~ /[A-Za-z0-9]+_[A-Za-z0-9]+/)
        {
            if (!defined($uniprot_to_sp_hash{$uniprot}) &&
                !defined($sp_to_uniprot_hash{$swissprot}))
            {
                $uniprot_to_sp_hash{$uniprot}   = $swissprot;
                $sp_to_uniprot_hash{$swissprot} = $uniprot;
            }
        }
        
        # store annotation, do not overwrite existing annotation
        if ($accession =~ /[A-Za-z0-9]/)
        {
            # deal with multiple accessions per accession
            @temp_array = split /;/, $accession;
            
            foreach $temp_accession (@temp_array)
            {
                if ($symbol ne '' &&
                    !defined($annotation_hash{$temp_accession}{symbol}))
                {
                    $annotation_hash{$temp_accession}{symbol} = $symbol;
                }
                if ($species ne '' &&
                    !defined($annotation_hash{$temp_accession}{species}))
                {
                    $annotation_hash{$temp_accession}{species} = $species;
                }
                if ($desc_str ne '' &&
                    !defined($annotation_hash{$temp_accession}{desc}))
                {
                    $annotation_hash{$temp_accession}{desc} = $desc;
                }
            }
        }
        if ($uniprot =~ /[A-Za-z0-9]/)
        {
            if ($symbol ne '' &&
                !defined($annotation_hash{$uniprot}{symbol}))
            {
                $annotation_hash{$uniprot}{symbol} = $symbol;
            }
            if ($species ne '' &&
                !defined($annotation_hash{$uniprot}{species}))
            {
                $annotation_hash{$uniprot}{species} = $species;
            }
            if ($desc_str ne '' &&
                !defined($annotation_hash{$uniprot}{desc}))
            {
                $annotation_hash{$uniprot}{desc} = $desc;
            }

            if ($symbol ne '' &&
                !defined($annotation_hash{$stripped}{symbol}))
            {
                $annotation_hash{$stripped}{symbol} = $symbol;
            }
            if ($species ne '' &&
                !defined($annotation_hash{$stripped}{species}))
            {
                $annotation_hash{$stripped}{species} = $species;
            }
            if ($desc_str ne '' &&
                !defined($annotation_hash{$stripped}{desc}))
            {
                $annotation_hash{$stripped}{desc} = $desc;
            }
        }
        if ($swissprot =~ /[A-Za-z0-9]/)
        {
            if ($symbol ne '' &&
                !defined($annotation_hash{$swissprot}{symbol}))
            {
                $annotation_hash{$swissprot}{symbol} = $symbol;
            }
            if ($species ne '' &&
                !defined($annotation_hash{$swissprot}{species}))
            {
                $annotation_hash{$swissprot}{species} = $species;
            }
            if ($desc_str ne '' &&
                !defined($annotation_hash{$swissprot}{desc}))
            {
                $annotation_hash{$swissprot}{desc} = $desc;
            }
        }
    }
}
close PROTEINS;


@accession_array = sort keys %annotation_hash;

# guess species of interest
%species_count_hash = ();
foreach $accession (@accession_array)
{
    $species = $annotation_hash{$accession}{species};
    
    if (defined($species))
    {
        if (!defined($species_count_hash{$species}))
        {
            $species_count_hash{$species} = 0;
        }
        $species_count_hash{$species} += 1;
    }
}

@species_array = keys %species_count_hash;
if (@species_array)
{
    @species_array = sort
                     {-$species_count_hash{$a} <=> -$species_count_hash{$b}}
                     @species_array;
    
    foreach $species (@species_array)
    {
        printf STDERR "%s\t%s\n", $species, $species_count_hash{$species};
    }
}

$species_majority = '';
if (@species_array)
{
    $species_majority = $species_array[0];
}



# header line
$line = <PEPTIDES>;
$line =~ s/[\r\n]+//g;
@array = split /\t/, $line;
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    $array[$i] =~ s/\s+/ /g;
    
    $field = $array[$i];
    $header_col_hash{$field} = $i;
#    $header_col_array[$i] = $field;
}
#$header_line = join "\t", @array;

$seq_col            = $header_col_hash{'Sequence'};
$proteins_col       = $header_col_hash{'Proteins'};
$leading_col        = $header_col_hash{'Leading razor protein'};
$id_col             = $header_col_hash{'id'};
$maxquant_group_col = $header_col_hash{'Protein group IDs'};
$peptides_pep_col   = $header_col_hash{'PEP'};

if (!defined($seq_col))
{
    die "ABORT -- cannot find Sequence column in peptides file\n";
}
if (!defined($proteins_col))
{
    die "ABORT -- cannot find Proteins column in peptides file\n";
}
if (!defined($leading_col))
{
    die "ABORT -- cannot find Leading razor protein column in peptides file\n";
}
if (!defined($id_col))
{
    die "ABORT -- cannot find id column in peptides file\n";
}
if (!defined($maxquant_group_col))
{
    die "ABORT -- cannot find Protein group IDs column in peptides file\n";
}
if (!defined($peptides_pep_col))
{
    die "ABORT -- cannot find PEP column in peptides file\n";
}



#$header_line_new  = $header_line;
$header_line_new  = 'Sequence';
$header_line_new .= "\t" . 'Accession_Condensed';
$header_line_new .= "\t" . 'Accession_Full';
$header_line_new .= "\t" . 'Group_Regrouped';
$header_line_new .= "\t" . 'Group_Maxquant';
$header_line_new .= "\t" . 'Group_Maxquant_Razor';
$header_line_new .= "\t" . 'MQ Unique';
$header_line_new .= "\t" . 'MQ Razor+Unique';
$header_line_new .= "\t" . 'Score';
$header_line_new .= "\t" . 'Species';
$header_line_new .= "\t" . 'Species_flag';
$header_line_new .= "\t" . 'Species_contaminant_flag';
$header_line_new .= "\t" . 'Other_species_flag';
$header_line_new .= "\t" . 'Unknown_species_flag';
$header_line_new .= "\t" . 'Non-Species_contaminant_flag';
$header_line_new .= "\t" . 'Reverse_flag';
$header_line_new .= "\t" . 'Trypsin_flag';


%missed_hash = ();

printf "%s\n", $header_line_new;

while(defined($line=<PEPTIDES>))
{
    $line =~ s/[\r\n]+//g;

    @array = split /\t/, $line;

    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\s+/ /g;
    }

    $seq            = $array[$seq_col];
    $proteins       = $array[$proteins_col];
    $leading        = $array[$leading_col];
    $maxquant_group = $array[$maxquant_group_col];
#    $pep_score      = $array[$peptides_pep_col];
    
    $proteins =~ s/\s+//;
    $leading  =~ s/\s+//;
    
    # try to deduce human, contaminant, other species, reverse, etc.
    
    # real human proteins are always in the Proteins col
    # they are almost all sp| and have HUMAN in them
    #  sometimes there is a bug and the sp| is left out of CON__ entries...
    # ignore the Potential contaminant column, it doesn't seem to work well
    #
    # Proteins field never contains any reverse hits
    
    $species_flag         = 0;
    $species_con_flag     = 0;
    $other_species_flag   = 0;
    $unknown_species_flag = 0;
    $non_species_con_flag = 0;
    $reverse_flag         = 0;
    $trypsin_flag         = 0;
    
    %condensed_accession_hash = ();
    %full_accession_hash = ();
    %tmp_species_hash = ();
    
    
    # first, build SwissProt <-> Uniprot mapping
    # ** WARNING ** the maxquant data may be woefully out of date, and
    #  have incorrect mappings; don't overwrite existing mappings
    #  only use it to fill in (old, out of date) missing mappings
    @subfields = split /;/, $proteins;
    for ($i = 0; $i < @subfields; $i++)
    {
        $uniprot   = '';
        $swissprot = '';

        $subfield = $subfields[$i];

        # handle full uniprot/sp accession format
        if ($subfield =~ /(sp|tr)\|([A-Za-z0-9\.-]+)\|([A-Za-z0-9\._-]+)/)
        {
            $uniprot   = $2;
            $swissprot = $3;
            $stripped  = $uniprot;
            $stripped  =~ s/(\-|\.)\d+$//;

            if (!defined($uniprot_to_sp_hash{$uniprot}) &&
                !defined($sp_to_uniprot_hash{$swissprot}))
            {
                $uniprot_to_sp_hash{$uniprot}   = $swissprot;
                $sp_to_uniprot_hash{$swissprot} = $uniprot;
            }
            # it is out of date, try to correct the gene symbol
            elsif (defined($uniprot_to_sp_hash{$uniprot}))
            {
                $swissprot_old = $swissprot;
                $swissprot     = $uniprot_to_sp_hash{$uniprot};
                
                if ($swissprot_old ne $swissprot)
                {
                    $outdated_sp_hash{$swissprot_old} = $swissprot;
                }
            }

            $full_accession_hash{$uniprot} = 1;
            $full_accession_hash{$swissprot} = 1;


            $species = '';
            if (defined($annotation_hash{$uniprot}{species}))
            {
                $species = $annotation_hash{$uniprot}{species};
            }
            elsif (defined($annotation_hash{$stripped}{species}))
            {
                $species = $annotation_hash{$stripped}{species};
            }
            elsif (defined($annotation_hash{$swissprot}{species}))
            {
                $species = $annotation_hash{$swissprot}{species};
            }
            elsif (defined($annotation_hash{$swissprot_old}{species}))
            {
                $species = $annotation_hash{$swissprot_old}{species};
            }
            elsif ($swissprot =~ /_([A-Za-z]+)/)
            {
                $species = $1;
                $species =~ s/HUMAN/Homo sapiens/;
                $species =~ s/MOUSE/Mus musculus/;
                $species =~ s/BOVIN/Bos taurus/;
            }
            
            if ($species ne '')
            {
                $tmp_species_hash{$species} = 1;
            }

            if ($species eq '')
            {
                $unknown_species_flag = 1;

                next;
            }
            if ($species ne $species_majority)
            {
                $other_species_flag = 1;

                next;
            }
            
            $species_accession_hash{$uniprot}   = 1;
            $species_accession_hash{$swissprot} = 1;
            $species_accession_hash{$stripped}  = 1;
            
            # store the accession in our growing list
            # store SwissProt, since it is more descriptive
            #
            # if the swissprot doesn't map to a geneid, but the uniprot
            #  does, then keep the uniprot instead of the swissprot
            if (!defined($accession_geneid_hash{$swissprot}) &&
                 (defined($accession_geneid_hash{$uniprot}) ||
                  defined($accession_geneid_hash{$stripped})))
            {
                $condensed_accession_hash{$uniprot} = 1;
            }
            else
            {
                $condensed_accession_hash{$swissprot} = 1;
            }

            $species_flag = 1;
        }


        # add in a few known species that slipped through

        # clean database stuff
        $subfield =~ s/(sp|tr)\|[^\|\;\,]+\|([^\|\;\,]+)/$2/g;
        $subfield =~ s/ENSEMBL:([A-Za-z0-9-\.]+)/$1/ig;
        $subfield =~ s/REFSEQ:([A-Za-z0-9\.]+)/$1/ig;
        $subfield =~ s/^CON__//;

        $stripped = $subfield;
        $stripped =~ s/(\-|\.)\d+$//;

        $species = '';
        if (defined($annotation_hash{$subfield}{species}))
        {
            $species = $annotation_hash{$subfield}{species};
        }
        elsif (defined($annotation_hash{$stripped}{species}))
        {
            $species = $annotation_hash{$stripped}{species};
        }
        elsif ($swissprot =~ /_([A-Za-z]+)/)
        {
            $species = $1;
            $species =~ s/HUMAN/Homo sapiens/;
            $species =~ s/MOUSE/Mus musculus/;
            $species =~ s/BOVIN/Bos taurus/;
        }

        # do not flag as non-species, as we skipped obsolete mappings before
        if ($species ne '' && $species eq $species_majority)
        {
            $species_accession_hash{$uniprot}  = 1;
            $species_accession_hash{$stripped} = 1;
            $species_flag                      = 1;
        }
    }

    
    # strip database and CON__, only add if they aren't already included
    # anything without HUMAN in it will be assumed to be non-human later
    for ($i = 0; $i < @subfields; $i++)
    {
        $uniprot   = '';
        $swissprot = '';

        $subfield = $subfields[$i];

        # clean database stuff
        $subfield =~ s/(sp|tr)\|[^\|\;\,]+\|([^\|\;\,]+)/$2/g;
        $subfield =~ s/ENSEMBL:([A-Za-z0-9-\.]+)/$1/ig;
        $subfield =~ s/REFSEQ:([A-Za-z0-9\.]+)/$1/ig;
        
        $was_con_flag = 0;
        if ($subfield =~ /^CON__/)
        {
            $was_con_flag = 1;
        }
        
        $subfield =~ s/^CON__//;
        $stripped = $subfield;
        $stripped =~ s/(\-|\.)\d+$//;

        # assume these are species of interest
        # we need to eat the outdated sp mappings here
        if (defined($uniprot_to_sp_hash{$subfield}) ||
            defined($sp_to_uniprot_hash{$subfield}) ||
            defined($outdated_sp_hash{$subfield}) ||
            defined($outdated_sp_hash{$stripped}) ||
            defined($species_accession_hash{$subfield}) ||
            defined($species_accession_hash{$stripped}))
        {
            if (defined($uniprot_to_sp_hash{$subfield}))
            {
                $uniprot   = $subfield;
                $swissprot = $uniprot_to_sp_hash{$subfield};
            }
            elsif (defined($sp_to_uniprot_hash{$subfield}))
            {
                $swissprot = $subfield;
                $uniprot   = $sp_to_uniprot_hash{$subfield};
            }
            # we don't have a swissprot <-> uniprot mapping for this one
            else
            {
                $uniprot   = $subfield;
                $swissprot = $subfield;
            }

            $full_accession_hash{$uniprot}   = 1;
            $full_accession_hash{$swissprot} = 1;

            $species = '';
            if (defined($annotation_hash{$uniprot}{species}))
            {
                $species = $annotation_hash{$uniprot}{species};
            }
            elsif (defined($annotation_hash{$stripped}{species}))
            {
                $species = $annotation_hash{$stripped}{species};
            }
            elsif (defined($annotation_hash{$swissprot}{species}))
            {
                $species = $annotation_hash{$swissprot}{species};
            }
            elsif ($swissprot =~ /_([A-Za-z]+)/)
            {
                $species = $1;
                $species =~ s/HUMAN/Homo sapiens/;
                $species =~ s/MOUSE/Mus musculus/;
                $species =~ s/BOVIN/Bos taurus/;
            }

            # if the swissprot doesn't map to a geneid, but the uniprot
            #  does, then keep the uniprot instead of the swissprot
            if (!defined($accession_geneid_hash{$swissprot}) &&
                 (defined($accession_geneid_hash{$uniprot}) ||
                  defined($accession_geneid_hash{$stripped})))
            {
                $condensed_accession_hash{$uniprot} = 1;
            }
            else
            {
                $condensed_accession_hash{$swissprot} = 1;
            }

            if ($species ne '')
            {
                $tmp_species_hash{$species} = 1;
            }

            if ($species eq '')
            {
                next;
            }
            elsif ($species ne $species_majority)
            {
                $other_species_flag = 1;
            
                next;
            }

            if ($was_con_flag)
            {
                $species_con_flag = 1;
                $species_con_hash{$uniprot}   = 1;
                $species_con_hash{$stripped}  = 1;
                $species_con_hash{$swissprot} = 1;
            }

            # Maxquant can screw up and not report the sp| stuff sometimes
            # we can catch some of these if we saw them earlier
            $species_flag = 1;
            $species_accession_hash{$uniprot}   = 1;
            $species_accession_hash{$swissprot} = 1;
        }

        # not already handled, generally not species of interest
        elsif (!defined($full_accession_hash{$subfield}))
        {
            $species = '';
            if (defined($annotation_hash{$subfield}{species}))
            {
                $species = $annotation_hash{$subfield}{species};
            }
            elsif (defined($annotation_hash{$stripped}{species}))
            {
                $species = $annotation_hash{$stripped}{species};
            }
            elsif ($subfield =~ /_([A-Za-z]+)/)
            {
                $species = $1;
                $species =~ s/HUMAN/Homo sapiens/;
                $species =~ s/MOUSE/Mus musculus/;
                $species =~ s/BOVIN/Bos taurus/;
            }
            # This shouldn't happen, unless the accessions are obsolete.
            #
            # The only missing hits so far are Q9UE12 and ENSP00000377550
            #  both of which are obsolete accessions in the mysterious extra
            #  contaminant database.
            # Both also contain the current _HUMAN identifier in the same field.
            if ($species ne '' && $species eq $species_majority)
            {
                if (!defined($missed_hash{$subfield}))
                {
                    printf STDERR "MISSED\t%s\t%s\n", $species, $subfield;
                }

                $missed_hash{$subfield} = 1;
                
                $species_flag = 1;
            }
            else
            {
                if ($species eq '')
                {
                    $unknown_species_flag = 1;
                }
                else
                {
                    $other_species_flag = 1;
                }
            }

            if (defined($trypsin_hash{$subfield}))
            {
                $trypsin_flag = 1;
                
                # replace uniprot trypsin accession with swissprot accession
                $swissprot = $trypsin_hash{$subfield};
                if ($swissprot =~ /^([^_]{3,})_[^_]{3,}$/)
                {
                    $subfield = $swissprot;
                }
            }
            
            $condensed_accession_hash{$subfield} = 1;
            $full_accession_hash{$subfield}      = 1;
            
            if ($was_con_flag)
            {
                $non_species_con_flag = 1;
                # $unmapped_con_hash{$subfield} = 1;
            }
        }

        if ($species ne '')
        {
            $tmp_species_hash{$species} = 1;
        }
    }
    
    # leading proteins can contain REV__ hits
    # they should not contain any other hits not already in Proteins
    @subfields = split /;/, $leading;
    for ($i = 0; $i < @subfields; $i++)
    {
        $subfield = $subfields[$i];

        # clean database stuff
        $subfield =~ s/(sp|tr)\|[^\|\;\,]+\|([^\|\;\,]+)/$2/g;
        $subfield =~ s/ENSEMBL:([A-Za-z0-9-\.]+)/$1/ig;
        $subfield =~ s/REFSEQ:([A-Za-z0-9\.]+)/$1/ig;
        
        if ($subfield =~ /^REV__/)
        {
            $reverse_flag = 1;
            
            # only include a reverse hit if not already a non-reverse protein
            if ($species_flag == 0 && $other_species_flag == 0)
            {
                $condensed_accession_hash{$subfield} = 1;
                $full_accession_hash{$subfield}      = 1;
            }
        }
    }
    
    # remove outdated swissprot from condensed accession list
    foreach $accession (keys %condensed_accession_hash)
    {
        $updated = $outdated_sp_hash{$accession};
        
        if (defined($updated) &&
            $accession ne $updated &&
            defined($condensed_accession_hash{$updated}))
        {
            delete $condensed_accession_hash{$accession};
        }
    }
    
    $condensed_accession_str = join ';', sort cmp_accessions keys %condensed_accession_hash;
    $full_accession_str      = join ';', sort cmp_accessions keys %full_accession_hash;
    $group                   = $condensed_accession_str;

    $species_str = '';
    @tmp_species_array =
        sort {-$species_count_hash{$a} <=> -$species_count_hash{$b}}
             keys %tmp_species_hash;
    if (@tmp_species_array)
    {
        $species_str = join ';', @tmp_species_array;
    }
    
    %temp_geneid_hash = ();
    foreach $accession (keys %full_accession_hash)
    {
        $geneid_str = $accession_geneid_hash{$accession};
        
        if (defined($geneid_str) && $geneid_str =~ /[0-9]/)
        {
            $geneid_str =~ s/\s+/\;/g;
            @temp_geneid_array = split ';', $geneid_str;
        
            foreach $geneid (@temp_geneid_array)
            {
                $temp_geneid_hash{$geneid} = 1;
            }
        }
        # try stripping off the isoform and/or version from the end
        else
        {
            $stripped = $accession;
            $stripped =~ s/(\-|\.)\d+$//;

            $geneid_str = $accession_geneid_hash{$stripped};

            if (defined($geneid_str) && $geneid_str =~ /[0-9]/)
            {
                $geneid_str =~ s/\s+/\;/g;
                @temp_geneid_array = split ';', $geneid_str;
        
                foreach $geneid (@temp_geneid_array)
                {
                    $temp_geneid_hash{$geneid} = 1;
                }
            }
        }
    }
    @temp_geneid_array = sort keys %temp_geneid_hash;
    if (@temp_geneid_array)
    {
        $group = join ';', @temp_geneid_array;
    }
    
    # Excel Remove Duplicates *does not remove all duplicates* (!!) if the
    #  list is large, contains mixed numbers and text, and some of the fields
    #  start with a digit.  This is a long standing Excel bug.
    #
    # http://superuser.com/questions/572226/remove-duplicates-feature-does-not-remove-all-duplicates
    #
    # Tacking on a non-number to the beginning avoids the issue entirely.
    $group = 'g:' . $group;

if (0)
{
    printf "%s\t%s\t%s", $condensed_accession_str, $proteins, $leading;
    printf "\t%s\t%s\t%s\t%s\t%s",
        $species_flag, $species_con_flag,
        $other_species_flag, $reverse_flag, $trypsin_flag;
    printf "\n";
}

    @mq_group_array = split /;/, $maxquant_group;

    # check to see if any of the groups are missing from the proteins file
    foreach $group (@mq_group_array)
    {
        # initialize them to horrible values
        if (!defined($maxquant_group_hash{$group}))
        {
            printf STDERR
              "WARNING -- protein group %s missing from proteinGroups file\n",
              $group;
        
            $maxquant_group_hash{$group}{'razor'}    = 1;    # at least 1
            $maxquant_group_hash{$group}{'score'}    = -9E9; # garbage
            $maxquant_group_hash{$group}{'unique'}   = 0;    # assume no unique
            $maxquant_group_hash{$group}{'peptides'} = 1;    # at least 1
        }
    }

    # pick a razor group
    @mq_group_array_sorted = sort cmp_razor (@mq_group_array);
#    $maxquant_group     = join ';', @mq_group_array;
    $razor_group        = $mq_group_array_sorted[0];
    
    # print tie warning
    if (@mq_group_array_sorted > 1)
    {
        $a = $mq_group_array_sorted[0];
        $b = $mq_group_array_sorted[1];
        if (cmp_razor($a, $b) == 0)
        {
            printf STDERR "WARNING -- TIED_RAZOR_GROUPS:\t%s\t%s\n",
                $a, $b;
        }
    }

    # create group strings
    $razor_str = '';
    foreach $group (@mq_group_array)
    {
        $razor = $maxquant_group_hash{$group}{'razor'};
        if (!defined($razor))
        {
            $razor = 0;
        }
        
        if ($razor_str ne '')
        {
            $razor_str .= ';'
        }
        
        $razor_str .= $razor;
    }
    $unique_str = '';
    foreach $group (@mq_group_array)
    {
        $unique = $maxquant_group_hash{$group}{'unique'};
        if (!defined($unique))
        {
            $unique = 0;
        }
        
        if ($unique_str ne '')
        {
            $unique_str .= ';'
        }
        
        $unique_str .= $unique;
    }
    $score_str = '';
    foreach $group (@mq_group_array)
    {
        $score = $maxquant_group_hash{$group}{'score'};
        if (!defined($score))
        {
            $score = 0;
        }
        
        if ($score_str ne '')
        {
            $score_str .= ';'
        }
        
        $score_str .= $score;
    }
    

#    printf "%s", $line;
    printf "%s", $seq;
    printf "\t%s", $condensed_accession_str;
    printf "\t%s", $full_accession_str;
    printf "\t%s", $group;
    printf "\t%s", $maxquant_group;
    printf "\t%s", $razor_group;
    printf "\t%s", $unique_str;
    printf "\t%s", $razor_str;
    printf "\t%s", $score_str;
    printf "\t%s", $species_str;
    printf "\t%s", $species_flag;
    printf "\t%s", $species_con_flag;
    printf "\t%s", $other_species_flag;
    printf "\t%s", $unknown_species_flag;
    printf "\t%s", $non_species_con_flag;
    printf "\t%s", $reverse_flag;
    printf "\t%s", $trypsin_flag;
    printf "\n";
}
close PEPTIDES;
