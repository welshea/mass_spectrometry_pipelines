#!/usr/bin/perl -w


# 2026-03-31:  expand Usage statement
# 2025-07-17:  add Selenocysteine stuff
# 2024-09-03:  extract activity annotation


use Scalar::Util qw(looks_like_number);
use File::Basename;


# sort numeric parts as numbers, not strings
sub cmp_args_alphanumeric
{
    my @array_a = split /([0-9]+)/, $_[0];
    my @array_b = split /([0-9]+)/, $_[1];
    my $count_a = @array_a;
    my $count_b = @array_b;
    my $min_count;
    my $i;
    my $j;
    
    $min_count = $count_a;
    if ($count_b < $min_count)
    {
        $min_count = $count_b;
    }
    
    for ($i = 0; $i < $min_count; $i += 2)
    {
        # even fields sort alphabetically
        if ($array_a[$i] lt $array_b[$i]) { return -1; }
        if ($array_a[$i] gt $array_b[$i]) { return  1; }
        
        # odd fields sort numerically
        $j = $i + 1;
        if ($j < $min_count)
        {
            if ($array_a[$j] < $array_b[$j]) { return -1; }
            if ($array_a[$j] > $array_b[$j]) { return  1; }
        }
    }

    # sort shorter remaining portion first
    if ($count_a < $count_b) { return -1; }
    if ($count_a > $count_b) { return  1; }

    # this shouldn't ever trigger
    return $_[0] cmp $_[1];
}


# assume that annotation hash has been filled, preferably pruned already
sub cmp_args_go_source
{
    my $acc_sp   = $_[0];
    my $source_a = $_[1];
    my $source_b = $_[2];
    my @array_a;
    my @array_b;
    my $string_a;
    my $string_b;
    my $len_a;
    my $len_b;
    my $compare;
    
    @array_a  = sort keys %{$annotation_hash{$acc_sp}{'Activity'}{$source_a}};
    @array_b  = sort keys %{$annotation_hash{$acc_sp}{'Activity'}{$source_b}};

    $string_a = join "\|", @array_a;
    $string_b = join "\|", @array_b;

    $len_a    = length $string_a;
    $len_b    = length $string_b;


    # sort shorter merged terms first
    if ($len_a < $len_b)     { return -1; }
    if ($len_a > $len_b)     { return  1; }

    # sort smaller number of terms first
    if (@array_a < @array_b) { return -1; }
    if (@array_a > @array_b) { return  1; }
    
    # sort terms alphanumerically
    $compare = cmp_args_alphanumeric($string_a, $string_b);
    if ($compare)            { return $compare; }

    # sort source alphanumerically
    $compare = cmp_args_alphanumeric($source_a, $source_b);
    return $compare;
}


sub is_number
{
    # use what Perl thinks is a number first
    # this is purely for speed, since the more complicated REGEX below should
    #  correctly handle all numeric cases
    if (looks_like_number($_[0]))
    {
        # Perl treats infinities as numbers, Excel does not.
        #
        # Perl treats NaN or NaNs, and various mixed caps, as numbers.
        # Weird that not-a-number is a number... but it is so that
        # it can do things like nan + 1 = nan, so I guess it makes sense
        #
        if ($_[0] =~ /^[-+]*(Inf|NaN)/i)
        {
            return 0;
        }
        
        return 1;
    }

    # optional + or - sign at beginning
    # then require either:
    #  a number followed by optional comma stuff, then optional decimal stuff
    #  mandatory decimal, followed by optional digits
    # then optional exponent stuff
    #
    # Perl cannot handle American comma separators within long numbers.
    # Excel does, so we have to check for it.
    # Excel doesn't handle European dot separators, at least not when it is
    #  set to the US locale (my test environment).  I am going to leave this
    #  unsupported for now.
    #
    if ($_[0] =~ /^([-+]?)([0-9]+(,[0-9]{3,})*\.?[0-9]*|\.[0-9]*)([Ee]([-+]?[0-9]+))?$/)
    {
        # current REGEX can treat '.' as a number, check for that
        if ($_[0] eq '.')
        {
            return 0;
        }
        
        return 1;
    }
    
    return 0;
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


# sort numeric parts as numbers, not strings
# case-insensitive
sub cmp_args_alphanumeric_i
{
    my @array_a = split /([0-9]+)/, lc $_[0];
    my @array_b = split /([0-9]+)/, lc $_[1];
    my $count_a = @array_a;
    my $count_b = @array_b;
    my $min_count;
    my $i;
    my $j;
    
    $min_count = $count_a;
    if ($count_b < $min_count)
    {
        $min_count = $count_b;
    }
    
    for ($i = 0; $i < $min_count; $i += 2)
    {
        # even fields sort alphabetically
        if ($array_a[$i] lt $array_b[$i]) { return -1; }
        if ($array_a[$i] gt $array_b[$i]) { return  1; }
        
        # odd fields sort numerically
        $j = $i + 1;
        if ($j < $min_count)
        {
            if ($array_a[$j] < $array_b[$j]) { return -1; }
            if ($array_a[$j] > $array_b[$j]) { return  1; }
        }
    }

    # sort shorter remaining portion first
    if ($count_a < $count_b) { return -1; }
    if ($count_a > $count_b) { return  1; }

    # this shouldn't ever trigger
    return $_[0] cmp $_[1];
}


sub cmp_header
{
    my $order_a = $header_order_hash{$a};
    my $order_b = $header_order_hash{$b};
    
    if (defined($order_a) && defined($order_b))
    {
        return $order_a <=> $order_b;
    }
    
    if (defined($order_a)) { return -1; }
    if (defined($order_b)) { return  1; }

    return $a cmp $b;
}


sub cmp_accession
{
    my $acc_sp_a = $a;
    my $acc_sp_b = $b;
    my $acc_to_use_a;
    my $acc_to_use_b;
    my $symbol_a;
    my $symbol_b;
    my $score_a;
    my $score_b;
    my $desc_a;
    my $desc_b;
    my $type_a;    # SwissProt, Canonical, Non-nanonical, Deprecated
    my $type_b;
    my $parent_a;
    my $parent_b;
    my $comparison;
    my @temp_array;
    
    $type_a = 0;
    if (defined($canonical_to_sp_hash{$acc_sp_a}))
    {
        $type_a = 1;
    }
    elsif (defined($noncanonical_to_sp_hash{$acc_sp_a}))
    {
        $type_a = 2;
    }
    elsif (defined($deprecated_to_sp_hash{$acc_sp_a}))
    {
        $type_a = 3;
    }

    $type_b = 0;
    if (defined($canonical_to_sp_hash{$acc_sp_b}))
    {
        $type_b = 1;
    }
    elsif (defined($noncanonical_to_sp_hash{$acc_sp_b}))
    {
        $type_b = 2;
    }
    elsif (defined($deprecated_to_sp_hash{$acc_sp_b}))
    {
        $type_b = 3;
    }
    

    $acc_to_use_a = $acc_sp_a;
    $acc_to_use_b = $acc_sp_b;

    if ($type_a == 1)
    {
        @temp_array   = sort { cmp_args_alphanumeric_i($a, $b) }
                        keys %{$canonical_to_sp_hash{$acc_sp_a}};
        $acc_to_use_a = $temp_array[0];
    }
    if ($type_a == 2)
    {
        @temp_array   = sort { cmp_args_alphanumeric_i($a, $b) }
                        keys %{$noncanonical_to_sp_hash{$acc_sp_a}};
        $acc_to_use_a = $temp_array[0];
    }
    if ($type_a == 3)
    {
        @temp_array   = sort { cmp_args_alphanumeric_i($a, $b) }
                        keys %{$deprecated_to_sp_hash{$acc_sp_a}};
        $acc_to_use_a = $temp_array[0];
    }

    if ($type_b == 1)
    {
        @temp_array   = sort { cmp_args_alphanumeric_i($a, $b) }
                        keys %{$canonical_to_sp_hash{$acc_sp_b}};
        $acc_to_use_b = $temp_array[0];
    }
    if ($type_b == 2)
    {
        @temp_array   = sort { cmp_args_alphanumeric_i($a, $b) }
                        keys %{$noncanonical_to_sp_hash{$acc_sp_b}};
        $acc_to_use_b = $temp_array[0];
    }
    if ($type_b == 3)
    {
        @temp_array   = sort { cmp_args_alphanumeric_i($a, $b) }
                        keys %{$deprecated_to_sp_hash{$acc_sp_b}};
        $acc_to_use_b = $temp_array[0];
    }


    $symbol_a = '';
    $symbol_b = '';
    if (defined($annotation_hash{$acc_to_use_a}) &&
        defined($annotation_hash{$acc_to_use_a}{'Symbol'}))
    {
        @temp_array = sort { cmp_args_alphanumeric_i($a, $b) }
                      keys %{$annotation_hash{$acc_to_use_a}{'Symbol'}};
        $symbol_a = join '|', @temp_array;
    }
    if (defined($annotation_hash{$acc_to_use_b}) &&
        defined($annotation_hash{$acc_to_use_b}{'Symbol'}))
    {
        @temp_array = sort { cmp_args_alphanumeric_i($a, $b) }
                      keys %{$annotation_hash{$acc_to_use_b}{'Symbol'}};
        $symbol_b = join '|', @temp_array;
    }
    

    # sort missing symbols last
    if ($symbol_a ne '' && $symbol_b eq '') { return -1; }
    if ($symbol_a eq '' && $symbol_b ne '') { return  1; }
    

    # sort by symbol
    if ($symbol_a ne '' && $symbol_b ne '')
    {
        $comparison = cmp_args_alphanumeric_i($symbol_a, $symbol_b);
        if ($comparison)
        {
            return $comparison;
        }
    }


    # then by description if we don't have symbols for either
    if ($symbol_a eq '' && $symbol_b eq '')
    {
        $desc_a = '';
        $desc_b = '';
        if (defined($annotation_hash{$acc_to_use_a}) &&
            defined($annotation_hash{$acc_to_use_a}{'Description'}))
        {
            @temp_array = sort { cmp_args_alphanumeric_i($a, $b) }
                          keys %{$annotation_hash{$acc_to_use_a}{'Description'}};
            $desc_a = join '|', @temp_array;
        }
        if (defined($annotation_hash{$acc_to_use_b}) &&
            defined($annotation_hash{$acc_to_use_b}{'Description'}))
        {
            @temp_array = sort { cmp_args_alphanumeric_i($a, $b) }
                          keys %{$annotation_hash{$acc_to_use_b}{'Description'}};
            $desc_b = join '|', @temp_array;
        }


        # sort missing descriptions last
        if ($desc_a ne '' && $desc_b eq '') { return -1; }
        if ($desc_a eq '' && $desc_b ne '') { return  1; }


        # sort by description
        if ($desc_a ne '' && $desc_b ne '')
        {
            $comparison = cmp_args_alphanumeric_i($desc_a, $desc_b);
            if ($comparison)
            {
                return $comparison;
            }
        }
    }

    
    # then by accession type
    if ($type_a < $type_b) { return -1; }
    if ($type_a > $type_b) { return  1; }
    
    
    # then by completeness of isoform annotation (without mapping back to sp),
    # but only if both have gene symbols, since we're only sorting on
    # completness of annotation in order to pick the the "best" accession
    # per-symbol
    if ($symbol_a ne '' && $symbol_b ne '')
    {
        $score_a = 0;
        $score_b = 0;
        if (defined($annotation_hash{$acc_sp_a}))
        {
            @temp_array = keys %{$annotation_hash{$acc_sp_a}};
            $score_a = @temp_array;
        }
        if (defined($annotation_hash{$acc_sp_b}))
        {
            @temp_array = keys %{$annotation_hash{$acc_sp_b}};
            $score_b = @temp_array;
        }
        if ($score_a > $score_b) { return -1; }
        if ($score_a < $score_b) { return  1; }
    }
    
    
    # then by SwissProt parent
    $parent_a = '';
    $parent_b = '';
    if (defined($annotation_hash{$acc_to_use_a}) &&
        defined($annotation_hash{$acc_to_use_a}{'ParentSwissProt'}))
    {
        @temp_array = sort { cmp_args_alphanumeric_i($a, $b) }
            keys %{$annotation_hash{$acc_to_use_a}{'ParentSwissProt'}};
        $parent_a = join '|', @temp_array;
    }
    if (defined($annotation_hash{$acc_to_use_b}) &&
        defined($annotation_hash{$acc_to_use_b}{'ParentSwissProt'}))
    {
        @temp_array = sort { cmp_args_alphanumeric_i($a, $b) }
            keys %{$annotation_hash{$acc_to_use_b}{'ParentSwissProt'}};
        $parent_b = join '|', @temp_array;
    }

    # sort missing parents last
    if ($parent_a ne '' && $parent_b eq '') { return -1; }
    if ($parent_a eq '' && $parent_b ne '') { return  1; }

    # sort by parents
    if ($parent_a ne '' && $parent_b ne '')
    {
        $comparison = cmp_args_alphanumeric_i($parent_a, $parent_b);
        if ($comparison)
        {
            return $comparison;
        }
    }
    
    
    # then alphanumeric by accession
    return cmp_args_alphanumeric_i($acc_sp_a, $acc_sp_b);
}


sub print_usage_statement
{
    $program_name = basename($0);

    printf STDERR "Usage: $program_name [options] UP#_#.dat\n";
    printf STDERR "\n";
    printf STDERR "  Map as many UniProt accessions/isoforms to their parent SwissProt accession\n";
    printf STDERR "  as possible, classify them as SwissProt/Canonical/Non-canonical/Deprecated,\n";
    printf STDERR "  and annotate them with metadata we are interested in.\n";
    printf STDERR "\n";
    printf STDERR "  Isoform extracellular/helix metadata is cloned from the canonical isoform\n";
    printf STDERR "  and is not isoform-specific.  Sub-cellular locations are merged from all\n";
    printf STDERR "  isoforms combined.\n";
    printf STDERR "\n";
    printf STDERR "  Options:\n";
    printf STDERR "    --best-per-symbol     keep only most-annotated SwissProt per gene symbol\n";
    printf STDERR "\n";
    printf STDERR "  DAT file download location:\n";
    printf STDERR "    http://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/\n";
    printf STDERR "      reference_proteomes/Eukaryota/\n";
    printf STDERR "\n";
    printf STDERR "        Human: UP000005640/UP000005640_9606.dat.gz\n";
    printf STDERR "        Mouse: UP000000589/UP000000589_10090.dat.gz\n";
    printf STDERR "\n";
    printf STDERR "    see download site ../reference_proteomes/README for UP# species table\n";
}


# Example:
# ID   OR2T8_HUMAN             Reviewed;         312 AA.
# AC   A6NH00;
# DE   RecName: Full=Olfactory receptor 2T8;
# GN   Name=OR2T8; Synonyms=OR2T8P;
# DR   RefSeq; NP_001005522.1; NM_001005522.1.
# DR   RefSeq; XP_011542480.1; XM_011544178.2.
# DR   RefSeq; XP_016856661.1; XM_017001172.1.
# DR   Ensembl; ENST00000641945; ENSP00000493286; ENSG00000177462.
# DR   HGNC; HGNC:15020; OR2T8.
# DR   GeneID; 343172; -.
# FT   TOPO_DOM        1..26
# FT                   /note="Extracellular"
#
# ID   FM25C_HUMAN             Reviewed;          89 AA.
# AC   B3EWG5; B2RV02; Q5VTM1;
#
# DR   RefSeq; NP_001160060.1; NM_001166588.2. [Q9NR50-2]
# DR   RefSeq; NP_001248347.1; NM_001261418.1. [Q9NR50-3]
# DR   RefSeq; NP_065098.1; NM_020365.4. [Q9NR50-1]
#
# (for determining which isoform is canonical)
# CC       Name=1;
# CC         IsoId=Q5H8A4-1; Sequence=Displayed;
#
# (WARNING -- can line-wrap mid-phrase, contain Note=, other stuff)
# CC   -!- SUBCELLULAR LOCATION: Preautophagosomal structure membrane
# CC       {ECO:0000269|PubMed:20810658, ECO:0000305|PubMed:21893048}; Single-pass
# CC       type III membrane protein {ECO:0000305|PubMed:9794794}. Endoplasmic
# CC       reticulum membrane {ECO:0000269|PubMed:24837458}; Single-pass type III
# CC       membrane protein {ECO:0000305|PubMed:9794794}. Cell membrane,
# CC       sarcolemma, T-tubule {ECO:0000269|PubMed:9794794}. Note=Also detected
# CC       near the junctional sarcoplasmic reticulum (PubMed:9794794).
# CC       Concentrates at perinuclear structures (PubMed:21893048).
# CC       {ECO:0000269|PubMed:21893048, ECO:0000269|PubMed:9794794}.
# CC   -!- (other stuff...)
# CC   -------------------------------------------------------------------------
# CC   Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms
# CC   Distributed under the Creative Commons Attribution (CC BY 4.0) License
# CC   -------------------------------------------------------------------------
#
# (Isoforms can cause serious parsing issues, since they use the isoform Name,
# instead of the accession, so we'd need to build a map from the CC Name=
# and CC IsoId fields that are *after* the SUBCELLULAR LOCATION lines,
# which just is out of order, and not something I want to deal with right now.
# For now, skip parsing any subcellular location that is isoform-specific.
#
# CC   -!- SUBCELLULAR LOCATION: [Isoform PDE9A1]: Cell projection, ruffle
# CC       membrane {ECO:0000269|PubMed:17090334}. Cytoplasm, perinuclear region
# CC       {ECO:0000269|PubMed:17090334}. Golgi apparatus
# CC       {ECO:0000269|PubMed:17090334}. Endoplasmic reticulum
# CC       {ECO:0000269|PubMed:17090334}. Cell membrane, sarcolemma
# CC       {ECO:0000269|PubMed:25799991}.
# CC   -!- SUBCELLULAR LOCATION: [Isoform PDE9A2]: Cell projection, ruffle
# CC       membrane {ECO:0000269|PubMed:17090334}. Cytoplasm, perinuclear region
# CC       {ECO:0000269|PubMed:17090334}.
# CC   -!- SUBCELLULAR LOCATION: [Isoform PDE9A3]: Cytoplasm
# CC       {ECO:0000269|PubMed:17090334}. Endoplasmic reticulum
# CC       {ECO:0000269|PubMed:17090334}.
# CC   -!- SUBCELLULAR LOCATION: [Isoform PDE9A17]: Cytoplasm
# CC       {ECO:0000269|PubMed:17090334}. Endoplasmic reticulum
# CC       {ECO:0000269|PubMed:17090334}.



$header_order_hash{'Accession'}         =  0;  # manually printed first
$header_order_hash{'ParentSwissProt'}   =  1;
$header_order_hash{'AccTypeUniProt'}    =  2;  # manually after ParentSwissProt
$header_order_hash{'Length'}            =  3;
$header_order_hash{'Accession_Protein'} =  4;
$header_order_hash{'Accession_RNA'}     =  5;
$header_order_hash{'Ensembl'}           =  6;
$header_order_hash{'HGNC'}              =  7;
$header_order_hash{'GeneID'}            =  8;
$header_order_hash{'Symbol'}            =  9;
$header_order_hash{'Alias'}             = 10;
$header_order_hash{'ExCellLoop'}        = 11;
$header_order_hash{'InCellLoop'}        = 12;
$header_order_hash{'TMHelix'}           = 13;
$header_order_hash{'SubCellLoc'}        = 14;
$header_order_hash{'Activity'}          = 15;
$header_order_hash{'SelenoAA'}          = 16;    # contained Selenocyteine
$header_order_hash{'SelenoRelated'}     = 17;    # binds to / pathway
$header_order_hash{'Description'}       = 18;


$opt_best_per_symbol = 0;

# read in command line arguments
$syntax_error_flag = 0;
$num_files         = 0;
for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        if ($field eq '--best-per-symbol')
        {
            $opt_best_per_symbol = 1;
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
        elsif ($num_files == 1)
        {
            $filename_brenda = $field;
            $num_files++;
        }
    }
}

# default to stdin if no filename given
if ($num_files == 0)
{
    $filename = '-';
    $num_files = 1;
}


# print syntax error message
if ($num_files == 0 || $syntax_error_flag)
{
    print_usage_statement();
    exit(1);
}


# read in BRENDA data, if specified
if (defined($filename_brenda))
{
    open BRENDA, "$filename_brenda" or
                 die "can't open input file $filename_brenda\n";

    $in_ok_record_flag = 0;
    while(defined($line=<BRENDA>))
    {
        if ($line =~ /^ID\s/)
        {
            $line =~ s/[\r\n]+//g;

            $id   = $line;
            $id   =~ s/^ID\s+//;
            
            if ($id =~ /^([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)/)
            {
                # overwrite, since it could have (deleted, transferred, etc.)
                # after the identifier
                $id = $1;

                $in_ok_record_flag = 1;
            }
            # only the initial "spontaneous" entry
            else
            {
                printf STDERR "Skipping:\t%s\n", $id;
            }
        }

        if ($in_ok_record_flag && $line =~ /^RN\s/)
        {
            $line =~ s/[\r\n]+//g;
            
            $name = $line;
            $name =~ s/^RN\s+//;
            
            if ($name =~ /[A-Za-z0-9]/)
            {
                # assume only a single name per ID
                $brenda_name_hash{$id} = $name;
            }
        }
    }
    close BRENDA;
}

open INFILE, "$filename" or die "ABORT -- cannot open file $filename\n";

%annotation_hash = ();

$within_entry_flag = 0;
$prev_ft_type      = '';
$subcell_str       = '';
$prev_cc_bang      = '';

while(defined($line=<INFILE>))
{
    $line =~ s/[\r\n]+//g;
    
    $line_strip = $line;
    $line_strip =~ s/\.$//;
    
    # start new entry
    if ($line_strip =~ /^ID\s+(\S+)/)
    {
        $acc_sp = $1;
        $within_entry_flag = 1;
        $prev_ft_type      = '';
        $prev_cc_bang      = '';
        $subcell_str       = '';

        $length = '';
        if ($line_strip =~ /([0-9]+)\s+AA/)
        {
            $length = $1;
        }

        $annotation_hash{$acc_sp}{'Accession_Protein'}{$acc_sp} = 1;
        $annotation_hash{$acc_sp}{'ParentSwissProt'}{$acc_sp}   = 1;

        $annotation_major_keys_hash{'Accession_Protein'} = 1;
        $annotation_major_keys_hash{'ParentSwissProt'}   = 1;

        if ($length ne '')
        {
            $annotation_hash{$acc_sp}{'Length'}{$length} = 1;
            $annotation_major_keys_hash{'Length'}        = 1;
        }
    }
    
    # end an entry
    elsif ($line_strip =~ /^\/\//)
    {
        $within_entry_flag = 0;
        $prev_ft_type      = '';
        $prev_cc_bang      = '';
        $subcell_str       = '';
    }
    
    # can have multiple comma-delimited accessions
    # also many:many mappings (!!)
    #
    # GAL3A_HUMAN
    #   CC         IsoId=P0DPI2-1, P30042-1; Sequence=Displayed;
    # GAL3B_HUMAN
    #   CC         IsoId=A0A0B4J2D5-1, P30042-1; Sequence=Displayed;
    #
    # Isoforms after the first one may not actually still exist.
    # These appear to be deprecated, just like any accessions that come
    # after the first accession in the AC lines.  Replace the old deprecated
    # accessions that don't have -# isoform numbers with ones that do.
    #
    # CC       Name=Non-muscle; Synonyms=MLC3nm, LC17A, LC17-nm;
    # CC         IsoId=P60660-1, P16475-1; Sequence=Displayed;
    # CC       Name=Smooth muscle; Synonyms=MLC3sm, LC17B, LC17-sm;
    # CC         IsoId=P60660-2, P24572-2; Sequence=VSP_009735;
    #
    elsif ($within_entry_flag &&
           $line_strip =~ /^CC\s+IsoId\s*=\s*([^;]+).*Sequence\s*\=\s*(.*)/i)
    {
        $isoform_str = $1;
        $seq_str     = $2;
        
        @isoform_array = split /\s*,\s*/, $isoform_str;
        
        for ($i = 0; $i < @isoform_array; $i++)
        {
            $isoform =  $isoform_array[$i];
            $isoform =~ s/^\s+//;
            $isoform =~ s/\s+$//;

            # primary accessions
            if ($i == 0)
            {
                # canonical isoform
                if ($seq_str =~ /^Displayed/)
                {
                    $sp_to_canonical_hash{$acc_sp}{$isoform} = 1;
                    $canonical_to_sp_hash{$isoform}{$acc_sp} = 1;

                    $annotation_hash{$acc_sp}{'Accession_Protein'}{$isoform}
                        = 1;
                }

                # non-canonical isoform
                #   example: Q9NQ94-3 isn't otherwise annotated elsewhere
                else
                {
                    $annotation_hash{$isoform}{'Accession_Protein'}{$isoform}
                        = 1;
                    $noncanonical_to_sp_hash{$isoform}{$acc_sp} = 1;
                }
            }
            # deprecated
            else
            {
                $isoform_trunc =  $isoform;
                $isoform_trunc =~ s/\-[0-9]+$//;
                
                # delete the deprecated entry that doesn't have -# number
                delete $deprecated_to_sp_hash{$isoform_trunc};
                delete $sp_to_deprecated_hash{$acc_sp}{$isoform_trunc};
                delete $annotation_hash{$isoform_trunc};
            
                # replace it with the entry for which we have -# number
                $deprecated_to_sp_hash{$isoform}{$acc_sp} = 1;
                $sp_to_deprecated_hash{$acc_sp}{$isoform} = 1;
                $annotation_hash{$isoform}{'Accession_Protein'}{$isoform}
                    = 1;

                #printf STDERR "FOOBAR\t%s\n", $isoform;
            }
        }
    }
    
    # Uniprot accession(s)
    elsif ($within_entry_flag &&
           $line_strip =~ /^AC\s+(\S.*)/)
    {
        $acc_str = $1;
        $acc_str =~ s/\.$//;
        $acc_str = bless_delimiter_bar($acc_str);
        
        @array = split /\|/, $acc_str;
        
        foreach $acc_uniprot (@array)
        {
            if ($acc_uniprot =~ /[A-Za-z0-9]/)
            {
                # store first accession for filling in missing isoforms later
                if (!defined($first_acc_hash{$acc_sp}))
                {
                    $first_acc_hash{$acc_sp} = $acc_uniprot;
                    $annotation_hash{$acc_sp}{'Accession_Protein'}{$acc_uniprot}
                        = 1;
                }
                
                # store remaining deprecated accessions
                else
                {
                    $sp_to_deprecated_hash{$acc_sp}{$acc_uniprot} = 1;
                    $deprecated_to_sp_hash{$acc_uniprot}{$acc_sp} = 1;

                    $annotation_hash{$acc_uniprot}{'Accession_Protein'}{$acc_uniprot}
                        = 1;
                }
            
                $annotation_major_keys_hash{'Accession_Protein'} = 1;

                #printf STDERR "FOOBAR\t%s\n", $acc_uniprot;
            }
        }
    }
    
    # gene symbols and aliases
    # symbols and aliases can be on different lines, or on same line
    #
    # Example line wrap:
    #
    # GN   Name=FITM1 {ECO:0000255|HAMAP-Rule:MF_03229,
    # GN   ECO:0000312|HGNC:HGNC:33714};
    # GN   Synonyms=FIT1 {ECO:0000255|HAMAP-Rule:MF_03229};
    #
    # so far, I don't see any need to scan ahead/behind to un-wrap lines,
    # since it doesn't appear to matter for anything I want to store
    #
    # In the past, SwissProt/UniProtKB merged identical sequences from
    # multiple genes into a single entry representing multiple genes.
    # This was a bad idea.  They appear to have finally realized this and
    # started to demerge them in 2011, but have only de-merged *some* genes
    # over time, but still not all of them.  There are still 74 merged
    # Human SwissProt accessions (covering 205 genes) as of 2021_04.
    # 10 years is a long time to have not finished cleaning up Human....
    # As of 2024-04-22, there are still 71 merged Human SwissProt accessions.
    # ~1 demerge per year over 3 years is pretty sad :-(
    #
    # https://www.uniprot.org/help/redundancy
    #   "When different genes in the same species give rise to the same protein
    #    sequence, they were merged in a single UniProtKB/Swiss-Prot record and
    #    the gene names listed in the gene name subsection. See for example the
    #    entry for human histone H3.1 (P68431). However, we tend to demerge
    #    many of these entries, and for newly annotated proteins, generate
    #    separate sequence entries in case of multiple genes coding for
    #    identical protein sequences, e.g. P08409."
    #
    # If newly added genes aren't being merged, then they really should go back
    # and demerge *ALL* of the merged genes in order to fix it properly.
    # Alas, they don't seem to be too inclined to do that yet :-(
    #
    if ($within_entry_flag &&
        $line_strip =~ /^GN\s+(\S.*)/)
    {
        $acc_str = $1;
        $acc_str =~ s/\.$//;
        $acc_str = bless_delimiter_bar($acc_str);
        
        @array = split /\|/, $acc_str;
        
        foreach $field (@array)
        {
            if ($field =~ /^Name\s*\=\s*(\S.*)/)
            {
                $symbol = $1;
                
                $annotation_hash{$acc_sp}{'Symbol'}{$symbol} = 1;
                $annotation_major_keys_hash{'Symbol'} = 1;

                #printf STDERR "FOOBAR\t%s\n", $symbol;
            }
            elsif ($field =~ /^Synonyms*\s*\=\s*(\S.*)/)
            {
                $alias_str = $1;
                $alias_str = bless_delimiter_bar($alias_str);
                
                @array2 = split /\|/, $alias_str;
                
                foreach $alias (@array2)
                {
                    if ($alias =~ /[A-Za-z0-9]/)
                    {
                        $annotation_hash{$acc_sp}{'Alias'}{$alias} = 1;
                        $annotation_major_keys_hash{'Alias'} = 1;
                    }

                    #printf STDERR "FOOBAR\t%s\n", $alias;
                }
            }
        }
    }
    
    # description
    # *MUST* have exactly 3 spaces after ^DE, otherwise,
    # it would also includes other additional further-indented descriptions
    elsif ($within_entry_flag &&
           $line_strip =~ /^DE   RecName:\s*Full\s*\=\s*(\S.*)/)
    {
        $description =  $1;
        $description =~ s/\;$//;
        $description =~ s/\{[^{}]+\}//g;
        $description =~ s/^\s+//;
        $description =~ s/\s+$//;
        
        $annotation_hash{$acc_sp}{'Description'}{$description} = 1;
        $annotation_major_keys_hash{'Description'} = 1;

        #printf STDERR "FOOBAR\t%s\n", $description;
    }

    # transcript/protein accessions
    elsif ($within_entry_flag &&
           $line_strip =~ /^DR\s+EMBL;*\s(\S.*)/)
    {
        $accession_str = $1;
        
        # only keep non-genomic regions
        if (!($accession_str =~ /genomic/i))
        {
            @temp_array = split /\s*\;\s*/, $accession_str;
            
            $accession_rna = '';
            if (@temp_array)
            {
                $accession_rna = $temp_array[0];
            }

            $accession_prot = '';
            if (@temp_array > 1)
            {
                $accession_prot = $temp_array[1];
            }
            
            if ($accession_rna =~ /[A-Za-z0-9]/)
            {
                $annotation_hash{$acc_to_use}{'Accession_RNA'}{$accession_rna} = 1;
            }

            if ($accession_prot =~ /[A-Za-z0-9]/)
            {
                $annotation_hash{$acc_to_use}{'Accession_Protein'}{$accession_prot}
                    = 1;
            }
        }
    }
    
    # RefSeq
    elsif ($within_entry_flag &&
           $line_strip =~ /^DR\s+RefSeq;*\s(\S.*)/)
    {
        $accession_str = $1;

        # deal with isoforms, to some extent
        # assume there is only one isoform given
        # we'll finish dealing with isoforms more later
        #
        # must extract isoform before converting delimiters,
        # since it will turn spaces into |
        $isoform = '';
        if ($accession_str =~ s/\s*\.\s*\[([^\[\] ]+)\]\s*//)
        {
            $isoform = $1;
        }

        # store annotation under the non-canonical isoform instead
        $acc_to_use = $acc_sp;
        if ($isoform ne '' && !defined($canonical_to_sp_hash{$isoform}))
        {
            $acc_to_use = $isoform;

            # save non-canonical isoforms to deal with more later
            $noncanonical_to_sp_hash{$isoform}{$acc_sp} = 1;
        }

        $accession_str = bless_delimiter_bar($accession_str);
        
        @array = split /\|/, $accession_str;
        
        foreach $accession (@array)
        {
            $first3 = substr $accession, 0, 3;
            
            # protein
            if ($first3 =~ /P_$/)
            {
                $annotation_hash{$acc_to_use}{'Accession_Protein'}{$accession} = 1;
                $annotation_major_keys_hash{'Accession_Protein'} = 1;
            }
            # RNA
            elsif ($first3 =~ /M_$/)
            {
                $annotation_hash{$acc_to_use}{'Accession_RNA'}{$accession} = 1;
                $annotation_major_keys_hash{'Accession_RNA'} = 1;
            }

            #printf STDERR "FOOBAR\t%s\t%s\n", $acc_to_use, $accession;
        }
    }

    # Ensembl
    #
    # examples: ENSMUSP00000129828 (11-digits), ENSMUS = mouse
    #           MGP_C3HHeJ_P...
    #
    # species prefix is usually ENS[A-Z]{3}, or just ENS for human
    # sometimes prefix is MGP_[A-Za-z0-9]{2,}_ instead of ENS...
    #   MGP_AJ_ (shortest for now)   MGP_129S1SvImJ_ (longest for now)
    #
    # there is a human gene symbol ENST869, so require all 11 digits
    #
    # AGAP_ENSANGG[0-9]{11} are symbols for Anopheles gambiae str. PEST
    #
    elsif ($within_entry_flag &&
           $line_strip =~ /^DR\s+Ensembl;*\s(\S.*)/)
    {
        $accession_str = $1;

        # deal with isoforms, to some extent
        # assume there is only one isoform given
        # we'll finish dealing with isoforms more later
        #
        # must extract isoform before converting delimiters,
        # since it will turn spaces into |
        $isoform = '';
        if ($accession_str =~ s/\s*\.\s*\[([^\[\] ]+)\]\s*//)
        {
            $isoform = $1;
        }

        # store annotation under the non-canonical isoform instead
        $acc_to_use = $acc_sp;
        if ($isoform ne '' && !defined($canonical_to_sp_hash{$isoform}))
        {
            $acc_to_use = $isoform;

            # save non-canonical isoforms to deal with more later
            $noncanonical_to_sp_hash{$isoform}{$acc_sp} = 1;
        }

        @array = $accession_str =~
            m/(?:^|[^-A-Za-z0-9])((?:ENS(?:[A-Z]{3})*|MGP_[A-Za-z0-9]{2,}_)(?:[EGPRT]|FM|GT)[0-9]{11}(?:\.[0-9]+(?:_par_[A-Z])*)*)/ig;
        
        foreach $accession (@array)
        {
            # protein
            if ($accession =~ /P[0-9]{11}/)
            {
                $annotation_hash{$acc_to_use}{'Accession_Protein'}{$accession} = 1;
                $annotation_major_keys_hash{'Accession_Protein'} = 1;
            }
            # RNA
            elsif ($accession =~ /T[0-9]{11}/)
            {
                $annotation_hash{$acc_to_use}{'Accession_RNA'}{$accession} = 1;
                $annotation_major_keys_hash{'Accession_RNA'} = 1;
            }
            elsif ($accession =~ /G[0-9]{11}/)
            {
                $annotation_hash{$acc_to_use}{'Ensembl'}{$accession} = 1;
                $annotation_major_keys_hash{'Ensembl'} = 1;
            }
        
            #printf STDERR "FOOBAR\t%s\t%s\n", $acc_to_use, $accession;
        }
    }

    # HGNC id and symbol
    elsif ($within_entry_flag &&
           $line_strip =~ /^DR\s+HGNC;*\s(\S.*)/)
    {
        $hgnc_str = $1;

        $hgnc_str = $1;
        $hgnc_str =~ s/\.$//;
        $hgnc_str = bless_delimiter_bar($hgnc_str);
        
        @array = split /\|/, $hgnc_str;
        
        $hgnc_id = $array[0];
        $symbol  = $array[1];

        # must check before we store first HGNC
        # there can be multiple HGNC rows per accession
        #
        #
        # Disable this.
        #
        # While NCBI has good mappings to Ensembl, and both HGNC and
        # NCBI have better gene symbols than Ensembl, Ensembl does *NOT*
        # have good mappins to *them*.  It appears to be sequence
        # homology based, without any sanity checking or manual curation??
        # At any rate, the mappings are often to both a gene and its
        # readthrough, and sometimes to multiple highly related family
        # members.
        #
        # So, since the Ensembl --> HGNC/NCBI mappings themselves aren't
        # trustworthy, don't use the HGNC gene symbols, because they
        # will sometimes be from the wrong gene(s).
        #
        if (0 && $symbol =~ /[A-Za-z0-9]/)
        {
            # Overwrite any existing symbols with HGNC symbol
            if (defined($annotation_hash{$acc_sp}) &&
                !defined($annotation_hash{$acc_sp}{'HGNC'}) &&
                defined($annotation_hash{$acc_sp}{'Symbol'}) &&
                !defined($annotation_hash{$acc_sp}{'Symbol'}{$symbol}))
            {
                # store the old symbols as aliases
                @array = sort keys %{$annotation_hash{$acc_sp}{'Symbol'}};
                foreach $alias (@array)
                {
                    $annotation_hash{$acc_sp}{'Alias'}{$alias} = 1;
                    $annotation_major_keys_hash{'Alias'} = 1;
                }
            
                # purge the old symbols
                delete $annotation_hash{$acc_sp}{'Symbol'};
            }
        
            $annotation_hash{$acc_sp}{'Symbol'}{$symbol} = 1;
            $annotation_major_keys_hash{'Symbol'} = 1;
        }
        
        if ($hgnc_id =~ /[A-Za-z0-9]/)
        {
            $annotation_hash{$acc_sp}{'HGNC'}{$hgnc_id} = 1;
            $annotation_major_keys_hash{'HGNC'} = 1;
        }
        
        #printf STDERR "FOOBAR\t%s\t%s\t%s\n", $acc_sp, $hgnc_id, $symbol;
    }

    # Entrez GeneID
    elsif ($within_entry_flag &&
           $line_strip =~ /^DR\s+GeneID;*\s(\S.*)/)
    {
        $entrez_str = $1;

        $entrez_str = $1;
        $entrez_str =~ s/\.$//;
        $entrez_str = bless_delimiter_bar($entrez_str);
        
        @array = split /\|/, $entrez_str;
        
        $entrez_id  = $array[0];

        if ($entrez_id =~ /[0-9]/)
        {
            $annotation_hash{$acc_sp}{'GeneID'}{$entrez_id} = 1;
            $annotation_major_keys_hash{'GeneID'} = 1;

            #printf STDERR "FOOBAR\t%s\t%s\n", $acc_sp, $entrez_id;
        }
    }
    
    # domain definition line
    elsif ($within_entry_flag &&
           $line_strip =~ /^FT\s+(\S+)\s+([0-9.]+)/)
    {
        $prev_ft_type = $1;
        $range        = $2;
        
        # start..end
        if ($range =~ /^\s*([0-9]+)\.+([0-9]+)/)
        {
            $ft_start  = $1;
            $ft_end    = $2;
            $ft_length = abs($ft_end - $ft_start) + 1;
        }
        # single amino acid
        else
        {
            $ft_start  = $range;
            $ft_end    = $range;
            $ft_length = 1;
        }

        #printf STDERR "FOOBAR\t%s\t%s\n",
        #    $acc_sp, $line_strip;
    }
    
    # Extracellular loop
    elsif ($within_entry_flag &&
           $prev_ft_type eq 'TOPO_DOM' &&
           $line_strip =~ /^FT\s+\/note\s*\=\s*\"Extracellular/)
    {
        if (!defined($annotation_hash{$acc_sp}{'ExCellLoop'}))
        {
            $annotation_hash{$acc_sp}{'ExCellLoop'}{'Count'}  = 0;
            $annotation_hash{$acc_sp}{'ExCellLoop'}{'Ranges'} = '';
            $annotation_major_keys_hash{'ExCellLoop'}         = 1;
        }
        
        $count = $annotation_hash{$acc_sp}{'ExCellLoop'}{'Count'};
        $annotation_hash{$acc_sp}{'ExCellLoop'}{'Array'}[$count] = $ft_length;
        $annotation_hash{$acc_sp}{'ExCellLoop'}{'Count'} += 1;
        
        # Assume they are given in sequence order, since they appear to be.
        $range_str = $ft_start . '--' . $ft_end;
        if ($annotation_hash{$acc_sp}{'ExCellLoop'}{'Ranges'} ne '')
        {
            $annotation_hash{$acc_sp}{'ExCellLoop'}{'Ranges'} .=
                ':' . $range_str;
        }
        else
        {
            $annotation_hash{$acc_sp}{'ExCellLoop'}{'Ranges'} = $range_str;
        }
    
        #printf STDERR "FOOBAR\t%s\t%s\t%s\t%s\t%s\n",
        #    $acc_sp, $ft_start, $ft_end, $ft_length,
        #    $annotation_hash{$acc_sp}{'ExCellLoop'}{'Ranges'};

        #printf STDERR "FOOBAR\t%s\t%s\n",
        #    $acc_sp, $line_strip;
    }


    # Intracellular loop (anything not Extracellular)
    elsif ($within_entry_flag &&
           $prev_ft_type eq 'TOPO_DOM' &&
           $line_strip =~ /^FT\s+\/note\s*\=\s*\"([^";]+)/)
    {
        $note = $1;
        $note =~ s/\s+$//;
    
        if (!defined($annotation_hash{$acc_sp}{'InCellLoop'}))
        {
            $annotation_hash{$acc_sp}{'InCellLoop'}{'Count'}  = 0;
            $annotation_major_keys_hash{'InCellLoop'}         = 1;
        }
        
        $annotation_hash{$acc_sp}{'InCellLoop'}{'Count'}       += 1;
        $annotation_hash{$acc_sp}{'InCellLoop'}{'Type'}{$note}  = 1;
    
        #printf STDERR "FOOBAR\t%s\t%s\n",
        #    $acc_sp, $note;
    }

    # Trans-membrane regions
    # AJAP1_HUMAN is line-wrapped in its note= line, so no end "
    elsif ($within_entry_flag &&
           $prev_ft_type eq 'TRANSMEM' &&
           $line_strip =~ /^FT\s+\/note\s*\=\s*\"(.*)/)
    {
        $note = $1;
        
        # "Helical"
        # "Discontinuously helical"
        #
        if ($note =~ /helical/i)
        {
            if (!defined($annotation_hash{$acc_sp}{'TMHelix'}))
            {
                $annotation_hash{$acc_sp}{'TMHelix'}{'TotalCount'}  = 0;
                $annotation_hash{$acc_sp}{'TMHelix'}{'AnchorCount'} = 0;
                $annotation_major_keys_hash{'TMHelix'}              = 1;
            }

            $annotation_hash{$acc_sp}{'TMHelix'}{'TotalCount'} += 1;

            # "Signal-anchor for type II membrane protein"
            # "Anchor for type IV membrane protein"
            # "Signal-anchor"
            #
            if ($line_strip =~ /anchor/i)
            {
                $annotation_hash{$acc_sp}{'TMHelix'}{'AnchorCount'} += 1;
            }

            #printf STDERR "FOOBAR\t%s\t%s\t%s\n",
            #    $acc_sp, $prev_ft_type, $note;
        }

        # "Beta stranded"
        # Very few proteins annotated, not going to bother with it:
        #   CO9_HUMAN,   GSDMA_HUMAN
        #   VDAC1_HUMAN, VDAC2_HUMAN, VDAC3_HUMAN
    }
    
    # various biological stuff, first line of line-wrap
    # can have isoforms, which we'll include in the field name
    #
    # Ex: CC   -!- SUBCELLULAR LOCATION: [Isoform PDE9A1]:
    #
    # Ex: CC   -!- SUBCELLULAR LOCATION: [Gasdermin-A, N-terminal]:
    # (this one doesn't even have any IsoIds)
    #
    # (Sigh... the isoform can even line wrap...)
    # CC   -!- SUBCELLULAR LOCATION: [Processed sterol regulatory element-binding
    # CC       protein 1]: Nucleus {ECO:0000269|PubMed:11477106,
    # CC       ECO:0000269|PubMed:32322062}.
    #
    # Give up on the idea of parsing out the isoform here.
    # Merge all entries together and deal with them later.
    #
    elsif ($within_entry_flag &&
           $line =~ /^CC\s*-\!\-\s*([^:]+)\s*:\s*(.*)/)
    {
        $prev_cc_bang = $1;
        $rest_of_line = $2;
        
        ## clean up spaces and terminal :
        #$prev_cc_bang =~ s/^\s+//;
        #$prev_cc_bang =~ s/\s*:\s*$//;
        
        ## skip stuff involving isoforms, etc.
        #if ($rest_of_line =~ /^\[/)
        #{
        #    $prev_cc_bang = '';
        #
        #    next;
        #}
        
        # begin new subcellular location string
        if ($prev_cc_bang eq 'SUBCELLULAR LOCATION')
        {
            # start new string
            if ($subcell_str eq '')
            {
                $subcell_str = $rest_of_line;
            }
            # append to new ||| delimited entry of current string
            else
            {
                $subcell_str .= '|||' . $rest_of_line;
            }
        }
        
        #printf STDERR "FOOBAR\t|%s|\n", $prev_cc_bang;
        #printf STDERR "FOOBAR\t%s\t%s\t|%s|\n",
        #    $acc_sp, $prev_cc_bang, $rest_of_line;
    }

    # SUBCELLULAR LOCATION
    # line wrap, indent must be >= 4 spaces
    #
    elsif ($within_entry_flag &&
           $prev_cc_bang eq 'SUBCELLULAR LOCATION' &&
           $line =~ /^CC\s{4,}(.*)/)
    {
        $rest_of_line = $1;
        
        # NOTE -- this can insert too many spaces, clean them up later
        $subcell_str .= ' ' . $rest_of_line;

        #printf STDERR "FOOBAR\t|%s|\n", $rest_of_line;
    }
    
    # store CC -!- fields we've accumulated
    elsif ($within_entry_flag &&
           $prev_cc_bang ne '' &&
           !($line =~ /^CC/))
    {
        if ($subcell_str ne '')
        {
            @subentry_array = split /\|\|\|/, $subcell_str;

            foreach $subentry (@subentry_array)
            {
                # clean up spaces
                $subentry =~ s/^\s+//;
                $subentry =~ s/\s+$//;
                $subentry =~ s/\s+/ /g;
                
                # remove spaces after hyphens, probably from line wraps
                $subentry =~ s/\-\s+/\-/g;
                
                # isoform, etc.
                if ($subentry =~ s/^\[([^]]+)\]\s*:\s*//)
                {
                    $isoform_str  = $1;
                    $isoform_str .= '';    # shut up Perl warning
                }

                # strip out stuff we don't want
                $subentry =~ s/\s*\{[^{}]+\}\s*//g;

                # remove spaces after delimiters
                $subentry =~ s/\s*([;.])+\s*/$1/g;

                # remove Note= from end
                $subentry =~ s/Note\s*\=.*//;

                # . delimits different subcellular locations
                # ; delimits location followed by protein class
                @major_array = split /\./, $subentry;

                foreach $field_major (@major_array)
                {
                    # strip out the ; delimited protein class
                    $field_major =~ s/\;.*//;

                    # only keep the highest level subcellular location
                    $field_major =~ s/\,.*//;

                    #printf STDERR "FOOBAR\t%s\t%s\n",
                    #    $acc_sp, $field_major;

                    $annotation_hash{$acc_sp}{'SubCellLoc'}{$field_major} = 1;
                    $annotation_major_keys_hash{'SubCellLoc'}             = 1;
                }

                #printf STDERR "FOOBAR\t%s\t%s\n",
                #$acc_sp, $subentry;
            }
        }
        
        $prev_cc_bang = '';
    }

    # GO F: activity
    elsif ($within_entry_flag &&
           $line_strip =~ /^DR\s+GO;*\s.*F:(.*)/)
    {
        $go_activity =  $1;
        
        ## only include UniProt and GO_Central sources
        ## UniProtKB also includes UniProtKB-KW, UniProt-UniRule, etc.
        #if ($go_activity =~ /:(GO_Central)/i ||
        #    $go_activity =~ /:(UniProtKB(?:\-[A-Za-z0-9]+)*)/i)
        if ($go_activity =~ /:([^\. ]+)/i)
        {
            $source = $1;

            # skip sources we don't want
            if (!($source eq 'ARUK-UCL'))
            {
                $go_activity =~ s/;.*//;

                if ($go_activity =~ /\sactivity$/)
                {
                    $go_activity =~ s/\s+activity$//;

                    $annotation_hash{$acc_sp}{'Activity'}{$source}{$go_activity} = 1;
                    $annotation_major_keys_hash{'Activity'} = 1;
                }
            }
        }
    }

    # BRENDA
    elsif ($within_entry_flag &&
           $line_strip =~ /^DR\s+BRENDA;*\s+([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)/)
    {
        $brenda_id   = $1;
        $brenda_name = $brenda_name_hash{$brenda_id};
        
        # store if
        if (defined($brenda_name) && $brenda_name =~ /[A-Za-z0-9]/)
        {
            $annotation_hash{$acc_sp}{'Activity'}{'BRENDA'}{$brenda_name} = 1;
            $annotation_major_keys_hash{'Activity'} = 1;

            #printf STDERR "FOOBAR\t%s\t%s\n", $acc_sp, $brenda_name;
        }
    }

    # DR selenium/seleno
    if ($within_entry_flag &&
        $line_strip =~ /^DR\s+.*?selen(ium|oprot|ocys)/i)
    {
        if ($line_strip =~ /^DR\s+([^;]+);\s*([^;]+)/)
        {
            $seleno_db  = $1;
            $seleno_acc = $2;

            # skip "family" proteins, such as some TCDB 9.B.87.*.* proteins
            if (!($line_strip =~ /family/i))
            {
                $seleno_str = $seleno_db . '~' . $seleno_acc;

                $annotation_hash{$acc_sp}{'SelenoRelated'}{$seleno_str} = 1;
                $annotation_major_keys_hash{'SelenoRelated'} = 1;
            }
        }
    }

    # Selenocysteine itself
    if ($within_entry_flag &&
        $prev_ft_type eq 'NON_STD' &&
        $line_strip =~ /^FT\s+.*?Selenocysteine/i)
    {
        $annotation_hash{$acc_sp}{'SelenoAA'}   = 1;
        $annotation_major_keys_hash{'SelenoAA'} = 1;
    }
}


# Fill in missing isoforms.
# IsoId fields are only specified when there are >= 2 isoforms,
# otherwise, the first ACC listed is the canonical accession.
# No isoform -# is specified in the ACC field, so we infer it to be "1".
#
foreach $acc_sp (sort keys %first_acc_hash)
{
    # append "-1" to end of first ACC, since it is isoform 1 of 1
    $isoform = $first_acc_hash{$acc_sp} . '-1';

    if (!defined($sp_to_canonical_hash{$acc_sp}))
    {
        $sp_to_canonical_hash{$acc_sp}{$isoform} = 1;
        $canonical_to_sp_hash{$isoform}{$acc_sp} = 1;
    }
}


# clean up many:many canonical mappings
%isoform_to_purge_hash = ();
foreach $acc_sp (sort keys %sp_to_canonical_hash)
{
    foreach $isoform (sort keys %{$sp_to_canonical_hash{$acc_sp}})
    {
        # check to see if it maps to multiple SwissProt entries
        @temp_array = sort keys %{$canonical_to_sp_hash{$isoform}};
        
        # P30042-1 --> GAL3A_HUMAN, GAL3B_HUMAN
        #   only instance of many:1 in release 2021_04
        #   was de-merged in release 2018_07
        #
        if (@temp_array >= 2)
        {
            # Remove the problematic accession from canonical mappings,
            # since we don't know which it should really be assigned to.
            #
            # For P30042-1, it was entirely deleted in 2018_07,
            # leaving the original A0A0B4J2D5 and P0DPI2 mappings,
            # so this appears to be a reasonable thing to do.
            #
            # NOTE -- accession will remain in Accession_Protein hashes,
            # in case we encounter it at some point, so it can be back-mapped
            # to the many genes it mapped to.
            #
            foreach $acc_sp2 (@temp_array)
            {
                $isoform_to_purge_hash{$isoform}{$acc_sp}  = 1;
                $isoform_to_purge_hash{$isoform}{$acc_sp2} = 1;
            }
        }
    }
}


# purge offending canonical linkages,
# but leave any offending isoform annotation intact
foreach $isoform (sort keys %isoform_to_purge_hash)
{
    foreach $acc_sp (sort keys %{$isoform_to_purge_hash{$isoform}})
    {
        delete $annotation_hash{$acc_sp};
        delete $sp_to_canonical_hash{$acc_sp}{$isoform};
        delete $canonical_to_sp_hash{$isoform}{$acc_sp};
        
        @temp_array = keys %{$sp_to_canonical_hash{$acc_sp}};
        if (@temp_array == 0)
        {
            delete $sp_to_canonical_hash{$acc_sp};
        }

        @temp_array = keys %{$canonical_to_sp_hash{$isoform}};
        if (@temp_array == 0)
        {
            delete $canonical_to_sp_hash{$isoform};
        }

        printf STDERR "PURGE_MULTI_CANONICAL\t%s\t%s\n",
            $isoform, $acc_sp;
    }
}


# purge 1:many non-canonical:sp mappings
# this happens once we add in IsoId isoforms that aren't otherwise annotated
foreach $isoform (sort keys %noncanonical_to_sp_hash)
{
    @temp_array = sort keys %{$noncanonical_to_sp_hash{$isoform}};
    
    if (@temp_array >= 2)
    {
        $temp_str = join '|', @temp_array;
        printf STDERR "PURGE_NONCANONICAL\t%s\t%s\n",
            $isoform, $temp_str;
        
        delete $annotation_hash{$isoform};
        delete $noncanonical_to_sp_hash{$isoform};

        @temp_array = keys %{$noncanonoical_isoform_hash{$isoform}};
        if (@temp_array == 0)
        {
            delete $noncanonoical_isoform_hash{$isoform};
        }
    }
}


# check for presence in both canonical and noncanonical
# there should be only 1:1 sp:canonical mappings left at this point
foreach $isoform (sort keys %noncanonical_to_sp_hash)
{
    if (defined($canonical_to_sp_hash{$isoform}))
    {

        @temp_array  = sort keys %{$canonical_to_sp_hash{$isoform}};
        $temp_str    = join '|', @temp_array;

        @temp_array2 = sort keys %{$noncanonical_to_sp_hash{$isoform}};
        $temp_str2   = join '|', @temp_array2;

        printf STDERR "PURGE_CANON_NONCANON\t%s\t%s|%s\n",
            $isoform, $temp_str, $temp_str2;


        ## leave the canonical version, only remove the non-canonical version
        ##
        #foreach $acc_sp (@temp_array)
        #{
        #    delete $sp_to_canonical_hash{$acc_sp};
        #}
        #delete $canonical_to_sp_hash{$isoform};


        delete $noncanonical_to_sp_hash{$isoform};
        @temp_array2 = keys %{$noncanonoical_isoform_hash{$isoform}};
        if (@temp_array2 == 0)
        {
            delete $noncanonoical_isoform_hash{$isoform};
        }

        delete $annotation_hash{$isoform};
    }
}


# clean up deprecated accessions
@acc_dep_array = sort keys %deprecated_to_sp_hash;
foreach $acc_dep (@acc_dep_array)
{
    @temp_array = sort keys %{$deprecated_to_sp_hash{$acc_dep}};
    
    # delete any deprecated accessions that map to multiple current ones,
    # since we can't uniquely map those
    if (@temp_array >= 2)
    {
        foreach $acc_sp (@temp_array)
        {
            delete $annotation_hash{$acc_dep};
            delete $sp_to_deprecated_hash{$acc_sp}{$acc_dep};
            
            # delete key entirely if it is empty underneath it
            @test_array = keys %{$sp_to_deprecated_hash{$acc_sp}};
            if (@test_array == 0)
            {
                delete $sp_to_deprecated_hash{$acc_sp};
            }
        }
    
        delete $deprecated_to_sp_hash{$acc_dep};
        
        @temp_array = keys %{$deprecated_to_sp_hash{$acc_dep}};
        if (@temp_array == 0)
        {
            delete $deprecated_to_sp_hash{$acc_dep};
        }
    }
}


# remove any deprecated accessions that are found to still be in use
# this shouldn't happen anymore, now that much of it is dealt with elsewhere
foreach $isoform (sort keys %canonical_to_sp_hash)
{
    $isoform_trunc =  $isoform;
    $isoform_trunc =~ s/\-[0-9]+$//;

    if (defined($deprecated_to_sp_hash{$isoform_trunc}))
    {
        delete $deprecated_to_sp_hash{$isoform_trunc};
    
        printf STDERR "DEDEPRECATED\t%s\n", $isoform;
    }
}
foreach $isoform (sort keys %noncanonical_to_sp_hash)
{
    $isoform_trunc =  $isoform;
    $isoform_trunc =~ s/\-[0-9]+$//;
    
    if (defined($deprecated_to_sp_hash{$isoform_trunc}))
    {
        delete $deprecated_to_sp_hash{$isoform_trunc};

        printf STDERR "DEDEPRECATED\t%s\n", $isoform;
    }
}


@acc_sp_array = sort { cmp_args_alphanumeric_i($a, $b) }
                keys %annotation_hash;


# deal with GO terms that complete substrings of other GO terms
foreach $acc_sp (@acc_sp_array)
{
    @source_array = keys %{$annotation_hash{$acc_sp}{'Activity'}};

    foreach $source (@source_array)
    {
        %bad_hash = ();
        @go_array = keys %{$annotation_hash{$acc_sp}{'Activity'}{$source}};

        if (@go_array > 1)
        {
            @go_array = sort { (length $a <=> length $b) } @go_array;
            
            for ($i = 0; $i < @go_array; $i++)
            {
                $go_string_shorter = $go_array[$i];
                
                for ($j = $i + 1; $j < @go_array; $j++)
                {
                    $go_string_longer = $go_array[$j];
                    
                    # longer string completely contains shorter string
                    if ($go_string_longer =~ /\Q$go_string_shorter\E/i)
                    {
                        # store first-encountered in bad hash
                        if (!defined($bad_hash{$j}))
                        {
                            $bad_hash{$j} = $i;
                        }
                    }
                }
            }
            
            # remove flagged entries from annotation hash
            @bad_array = sort { $a <=> $b } keys %bad_hash;
            if (@bad_array)
            {
                foreach $j (@bad_array)
                {
                    delete $annotation_hash{$acc_sp}{'Activity'}{$source}{$go_array[$j]};

                    #printf STDERR "FOOBAR  %s  %s  %s ||| %s\n",
                    #    $acc_sp, $source, $go_array[$j],
                    #    $go_array[$bad_hash{$j}];
                }
            }
        }
    }
}


# infer missing LOC1######## symbols for large Entrez GeneID
#
# A0A0G2JNH3_HUMAN   LOC102723532
# DUTL_HUMAN         LOC100506422  (HGNC symbol: LINC03106)
# YP021_HUMAN        LOC107987235
#
# For DUTL_HUMAN, LINC03106 was likely assigned after the 2021_04
# Uniprot release, and NBCI Gene search works correctly for
# LOC100506422 (even though it isn't listed as a synonym...).
# The other two still don't have "proper" gene symbols as of
# 2024-04-17.  So, filling missing symbols with LOC# in these cases
# appears to be very reasonable.  I still wonder why NCBI doesn't
# list LOC100506422 as an alias for LINC03106, though....
#
foreach $acc_sp (@acc_sp_array)
{
    if (defined($annotation_hash{$acc_sp}) &&
        !defined($annotation_hash{$acc_sp}{'Symbol'}) &&
        defined($annotation_hash{$acc_sp}{'GeneID'}))
    {
        @subfield_array = sort { cmp_args_alphanumeric_i($a, $b) }
            keys %{$annotation_hash{$acc_sp}{'GeneID'}};

        %temp_hash = ();
        foreach $geneid (@subfield_array)
        {
            if (is_number($geneid) && $geneid >= 100000000)
            {
                $symbol = sprintf "LOC%s", $geneid;
                $temp_hash{$symbol} = 1;
            }
        }

        @subfield_array = sort { cmp_args_alphanumeric_i($a, $b) }
            keys %temp_hash;

        $symbol = join '|', @subfield_array;

        if ($symbol =~ /[A-Za-z0-9]/)
        {
            $annotation_hash{$acc_sp}{'Symbol'}{$symbol} = 1;

            printf STDERR "INFER_SYMBOL\t%s\t%s\n", $acc_sp, $symbol;
        }
    }
}


# populate ParentSwissProt after cleaning up
# separate lookups from other --> SwissProt,
# rather than SwissProt --> other, just in case I missed some book keeping
foreach $isoform (sort keys %canonical_to_sp_hash)
{
    foreach $acc_sp (sort keys %{$canonical_to_sp_hash{$isoform}})
    {
        if (defined($annotation_hash{$isoform}))
        {
            $annotation_hash{$isoform}{'ParentSwissProt'}{$acc_sp} = 1;
        }
    }
}
foreach $isoform (sort keys %noncanonical_to_sp_hash)
{
    foreach $acc_sp (sort keys %{$noncanonical_to_sp_hash{$isoform}})
    {
        if (defined($annotation_hash{$isoform}))
        {
            $annotation_hash{$isoform}{'ParentSwissProt'}{$acc_sp} = 1;
        }
    }
}
foreach $isoform (sort keys %deprecated_to_sp_hash)
{
    foreach $acc_sp (sort keys %{$deprecated_to_sp_hash{$isoform}})
    {
        if (defined($annotation_hash{$isoform}))
        {
            $annotation_hash{$isoform}{'ParentSwissProt'}{$acc_sp} = 1;
        }
    }
}


# re-create and re-sort after populating/cleaning
@acc_sp_array = sort cmp_accession keys %annotation_hash;


# fill in Selenocysteine-containing proteins that aren't annotated as
# being selenocysteine-related
foreach $acc_sp (@acc_sp_array)
{
    $header     = 'SelenoAA';
    $acc_to_use = $acc_sp;

    # use base accession for most isoform fields
    
    # non-canonical isoforms
    if (defined($noncanonical_to_sp_hash{$acc_sp}) &&
        $header ne 'Accession_Protein' &&
        $header ne 'Accession_RNA' &&
        !($header eq 'Ensembl' &&
          defined($annotation_hash{$acc_sp}{'Ensembl'})))
    {
        @temp_array = sort keys
            %{$noncanonical_to_sp_hash{$acc_sp}};

        $acc_to_use    = $temp_array[0];
    }
    # deprecated accessions
    elsif (defined($deprecated_to_sp_hash{$acc_sp}) &&
           $header ne 'Accession_Protein' &&
           $header ne 'Accession_RNA' &&
           !($header eq 'Ensembl' &&
             defined($annotation_hash{$acc_sp}{'Ensembl'})))
    {
        @temp_array = sort keys
            %{$deprecated_to_sp_hash{$acc_sp}};

        $acc_to_use    = $temp_array[0];
    }
    
    # use canonical SwissProt accession if Ensembl ENSG is missing
    # already handled immediately above, but leave code here for now
    if (0 && $header eq 'Ensembl' &&
        !defined($annotation_hash{$acc_to_use}{'Ensembl'}))
    {
        if (defined($noncanonical_to_sp_hash{$acc_sp}))
        {
            @temp_array = sort keys
                %{$noncanonical_to_sp_hash{$acc_sp}};

            $acc_temp   = $temp_array[0];
            
            if (defined($annotation_hash{$acc_temp}{'Ensembl'}))
            {
                $acc_to_use = $acc_temp;
            }
        }
        elsif (defined($deprecated_to_sp_hash{$acc_sp}))
        {
            @temp_array = sort keys
                %{$deprecated_to_sp_hash{$acc_sp}};

            $acc_temp   = $temp_array[0];
            
            if (defined($annotation_hash{$acc_temp}{'Ensembl'}))
            {
                $acc_to_use = $acc_temp;
            }
        }
    }

    # add amino acid sequence to list of evidence sources
    if (defined($annotation_hash{$acc_to_use}) &&
        defined($annotation_hash{$acc_to_use}{'SelenoAA'}))
    {
        $annotation_hash{$acc_sp}{'SelenoAA'} = 1;

        $seleno_str = 'AA' . '~' . $acc_to_use;
        
        $annotation_hash{$acc_sp}{'SelenoRelated'}{$seleno_str} = 1;
        $annotation_hash{$acc_to_use}{'SelenoRelated'}{$seleno_str} = 1;
    }
}


# remove any SelenoRelated proteins only supported by 60S ribosome,
# since a lot of ribosomal proteins are effectively false-positive in Reactome
# acc_to_use stuff already handled in the loop just above here
#
# also remove any that are only supported by TCDB, which has many "family"
# members that aren't actually related to selenocysteine processes
foreach $acc_sp (@acc_sp_array)
{
    if (defined($annotation_hash{$acc_sp}) &&
        defined($annotation_hash{$acc_sp}{'SelenoRelated'}))
    {
        @temp_array = sort keys %{$annotation_hash{$acc_sp}{'SelenoRelated'}};
        
        if (@temp_array == 1 &&
            ($temp_array[0] =~ /R-HSA-2408557/ ||
             $temp_array[0] =~ /^TCDB/))
        {
            delete $annotation_hash{$acc_sp}{'SelenoRelated'};
        }
    }
}


# short table of only "best" SwissProt per gene symbol
%acc_to_skip_hash = ();
if ($opt_best_per_symbol)
{
    %symbol_acc_hash = ();

    # we're fine with non-canonical isoforms missing symbols and geneid here,
    # since they will use the SwissProt entry for those later
    foreach $acc_sp (@acc_sp_array)
    {
        # skip non-canonical isoforms
        if (defined($noncanonical_to_sp_hash{$acc_sp}))
        {
            $acc_to_skip_hash{$acc_sp} = 1;
        
            next;
        }
        
        # skip any accessions without a symbol
        if (!defined($annotation_hash{$acc_sp}{'Symbol'}))
        {
            $acc_to_skip_hash{$acc_sp} = 1;

            next;
        }

        @temp_array = keys %{$annotation_hash{$acc_sp}};
        $score      = @temp_array;

        @temp_array = sort { cmp_args_alphanumeric_i($a, $b) }
                      keys %{$annotation_hash{$acc_sp}{'Symbol'}};
        $symbol     = join '|', @temp_array;
        
        $symbol_acc_hash{$symbol}{$acc_sp} = $score;
    }
    
    @symbol_array = sort { cmp_args_alphanumeric_i($a, $b) }
                    keys %symbol_acc_hash;

    # flag accessions to skip
    foreach $symbol (@symbol_array)
    {
        $best_score = -999;
        @temp_array = sort { cmp_args_alphanumeric_i($a, $b) }
                      keys %{$symbol_acc_hash{$symbol}};

        # determine best score per symbol
        foreach $acc_sp (@temp_array)
        {
            $score = $symbol_acc_hash{$symbol}{$acc_sp};
            
            if ($score > $best_score)
            {
                $best_score = $score;
            }
        }
        
        # skip any accessions worse than best score
        $count_best = 0;
        foreach $acc_sp (@temp_array)
        {
            $score = $symbol_acc_hash{$symbol}{$acc_sp};
            
            if ($score < $best_score)
            {
                $acc_to_skip_hash{$acc_sp} = 1;
            }
            else
            {
                $count_best++;
            }
        }
        
        # warn of any remaining ties
        if ($count_best >= 2)
        {
            printf STDERR "TIED_SWISSPROT_SYMBOL\t%s\n", $symbol;
        }
    }
}


# output table

# print header line
@header_array = sort cmp_header keys %annotation_major_keys_hash;

$header_line  = 'Accession';
foreach $header (@header_array)
{
    # foreach points to original within the array,
    # modifying the foreach variable modifies the value within the array (!!)
    $header_print = $header;

    # special case ExCellLoop columns
    if ($header eq 'ExCellLoop')
    {
        $header_line .=
            sprintf "\t%s:%s\t%s:%s\t%s:%s\t%s:%s\t%s:%s",
                'ExCellLoop', 'Count',
                'ExCellLoop', 'SumLength',
                'ExCellLoop', 'MinLength',
                'ExCellLoop', 'MaxLength',
                'ExCellLoop', 'Ranges';
    }
    # special case InCellLoop columns
    elsif ($header eq 'InCellLoop')
    {
        $header_line .=
            sprintf "\t%s:%s\t%s:%s",
                'InCellLoop', 'Count',
                'InCellLoop', 'Type';
    }
    # special case TMHelix columns
    elsif ($header eq 'TMHelix')
    {
        $header_line .=
            sprintf "\t%s:%s\t%s:%s",
                'TMHelix', 'TotalCount',
                'TMHelix', 'AnchorCount',
    }
    # insert PotentialCellSurface column after SubCellLoc
    elsif ($header eq 'SubCellLoc')
    {
        $header_line .=
            sprintf "\t%s\t%s", $header, 'PotentialCellSurface';
    }
    else
    {
        # HACK -- rename Length to LengthCanonical,
        #         since we don't have lengths for the isoforms
        if ($header eq 'Length')
        {
            $header_print = 'LengthCanonical';
        }
    
        $header_line .= "\t" . $header_print;


        # insert AccTypeUniProt column
        if ($header eq 'ParentSwissProt')
        {
            $header_line .= "\t" . 'AccTypeUniProt';
        }
    }
}
print "$header_line\n";


foreach $acc_sp (@acc_sp_array)
{
    # skip any accessions we don't want to print
    if (defined($acc_to_skip_hash{$acc_sp}))
    {
        next;
    }

    $line = $acc_sp;
    
    $acc_type = 'SwissProt';
    if (defined($noncanonical_to_sp_hash{$acc_sp}))
    {
        $acc_type = 'Non-canonical';
    }
    elsif (defined($deprecated_to_sp_hash{$acc_sp}))
    {
        $acc_type = 'Deprecated';
    }

    foreach $header (@header_array)
    {
        $acc_to_use = $acc_sp;

        # use base accession for most isoform fields
        
        # non-canonical isoforms
        if (defined($noncanonical_to_sp_hash{$acc_sp}) &&
            $header ne 'Accession_Protein' &&
            $header ne 'Accession_RNA' &&
            !($header eq 'Ensembl' &&
              defined($annotation_hash{$acc_sp}{'Ensembl'})))
        {
            @temp_array = sort keys
                %{$noncanonical_to_sp_hash{$acc_sp}};

            $acc_to_use    = $temp_array[0];
        }
        # deprecated accessions
        elsif (defined($deprecated_to_sp_hash{$acc_sp}) &&
               $header ne 'Accession_Protein' &&
               $header ne 'Accession_RNA' &&
               !($header eq 'Ensembl' &&
                 defined($annotation_hash{$acc_sp}{'Ensembl'})))
        {
            @temp_array = sort keys
                %{$deprecated_to_sp_hash{$acc_sp}};

            $acc_to_use    = $temp_array[0];
        }
        
        # use canonical SwissProt accession if Ensembl ENSG is missing
        # already handled immediately above, but leave code here for now
        if (0 && $header eq 'Ensembl' &&
            !defined($annotation_hash{$acc_to_use}{'Ensembl'}))
        {
            if (defined($noncanonical_to_sp_hash{$acc_sp}))
            {
                @temp_array = sort keys
                    %{$noncanonical_to_sp_hash{$acc_sp}};

                $acc_temp   = $temp_array[0];
                
                if (defined($annotation_hash{$acc_temp}{'Ensembl'}))
                {
                    $acc_to_use = $acc_temp;
                }
            }
            elsif (defined($deprecated_to_sp_hash{$acc_sp}))
            {
                @temp_array = sort keys
                    %{$deprecated_to_sp_hash{$acc_sp}};

                $acc_temp   = $temp_array[0];
                
                if (defined($annotation_hash{$acc_temp}{'Ensembl'}))
                {
                    $acc_to_use = $acc_temp;
                }
            }
        }


        # special case ExCellLoop columns
        if ($header eq 'ExCellLoop')
        {
            $count = $sum_length = $min_length = $max_length =
                     $range_str  = '---';
            
            if (defined($annotation_hash{$acc_to_use}{'ExCellLoop'}))
            {
                 $count     =
                     $annotation_hash{$acc_to_use}{'ExCellLoop'}{'Count'};
                 $range_str =
                     $annotation_hash{$acc_to_use}{'ExCellLoop'}{'Ranges'};

                 @subfield_array = sort { $a <=> $b }
                   @{$annotation_hash{$acc_to_use}{'ExCellLoop'}{'Array'}};
                 $min_length = $subfield_array[0];
                 $max_length = $subfield_array[$count-1];

                 $sum_length = 0;
                 foreach $length (@subfield_array)
                 {
                     $sum_length += $length;
                 }
            }

            $loop_str = sprintf "\t%s\t%s\t%s\t%s\t%s",
                            $count,
                            $sum_length, $min_length, $max_length,
                            $range_str;
            
            $line .= $loop_str;
        }

        # special case InCellLoop columns
        #
        # "Exoplasmic loop" appears to be mostly specific to flippases,
        # some of which can be located on the PM, while others are in other
        # organelles, such as golgi.  Flippase DNF1/DNF2 in yeast
        # (AT10B_HUMAN ortholog) is known to be in golgi in yeast, and is
        # annotated as "Exoplasmic loop".  So, these may or may not be
        # truly extracellular.  Either way, we probably don't want to be
        # designing antigens vs. flippases, since I'm thinking that would
        # cause nasty side effects.  So, safe to not count as
        # extracellular?
        #
        elsif ($header eq 'InCellLoop')
        {
            $total_count = $type_str = '---';
            
            if (defined($annotation_hash{$acc_to_use}{'InCellLoop'}))
            {
                 $total_count =
                     $annotation_hash{$acc_to_use}{'InCellLoop'}{'Count'};
                 
                 @temp_array = sort { cmp_args_alphanumeric_i($a, $b) }
                   keys %{$annotation_hash{$acc_to_use}{'InCellLoop'}{'Type'}};
                 $type_str = join '|', @temp_array;
            }

            $loop_str = sprintf "\t%s\t%s",
                            $total_count, $type_str;
            
            $line .= $loop_str;
        }

        # special case TMHelix columns
        elsif ($header eq 'TMHelix')
        {
            $total_count = $anchor_count = '---';
            
            if (defined($annotation_hash{$acc_to_use}{'TMHelix'}))
            {
                 $total_count  =
                   $annotation_hash{$acc_to_use}{'TMHelix'}{'TotalCount'};
                 $anchor_count =
                   $annotation_hash{$acc_to_use}{'TMHelix'}{'AnchorCount'};
            }

            $loop_str = sprintf "\t%s\t%s",
                            $total_count, $anchor_count;
            
            $line .= $loop_str;
        }

        # insert PotentialCellSurface column after SubCellLoc
        elsif ($header eq 'SubCellLoc')
        {
            $value_str = '---';

            if (defined($annotation_hash{$acc_to_use}{$header}))
            {
                @subfield_array = sort { cmp_args_alphanumeric_i($a, $b) }
                    keys %{$annotation_hash{$acc_to_use}{$header}};
                $value_str = join '|', @subfield_array;
            }
            
            # insert PotentialCellSurface column after SubCellLoc
            $potential_flag = 0;
            if (defined($annotation_hash{$acc_to_use}{'ExCellLoop'}) ||
                $value_str =~ /cell /i)
            {
                $potential_flag = 1;
            }

            $line .= "\t" . $value_str . "\t" . "$potential_flag";
        }

        # special case Activity
        elsif ($header eq 'Activity')
        {
            @source_array = sort { cmp_args_go_source($acc_to_use, $a, $b) }
                            keys %{$annotation_hash{$acc_to_use}{'Activity'}};

            $value_str      = '---';
            
            # maybe require at least 2 sources ??
            if (@source_array >= 1)
            {
                $source         = $source_array[0];
                @subfield_array = sort keys
                    %{$annotation_hash{$acc_to_use}{'Activity'}{$source}};
                $value_str      = join '|', @subfield_array;
            }

            $line .= "\t" . $value_str;
        }

        # special case SelenoRelated
        elsif ($header eq 'SelenoRelated')
        {
            @source_array = sort { cmp_args_alphanumeric($a, $b) }
                            keys %{$annotation_hash{$acc_to_use}{'SelenoRelated'}};

            $value_str      = '---';
            
            if (@source_array >= 1)
            {
                $value_str = join '|', @source_array;
            }

            $line .= "\t" . $value_str;
        }

        # special case SelenoAA
        elsif ($header eq 'SelenoAA')
        {
            $value_str = '---';

            if (defined($annotation_hash{$acc_to_use}{$header}))
            {
                $value_str = 'U';
            }

            $line .= "\t" . $value_str;
        }

        # all other fields
        else
        {
            $value_str = '---';

            if (defined($annotation_hash{$acc_to_use}{$header}))
            {
                @subfield_array = sort { cmp_args_alphanumeric_i($a, $b) }
                    keys %{$annotation_hash{$acc_to_use}{$header}};
                $value_str = join '|', @subfield_array;
            }

            $line .= "\t" . $value_str;


            # insert AccTypeUniProt column
            if ($header eq 'ParentSwissProt')
            {
                $line .= "\t" . $acc_type;
            }
        }
    }
    
    print "$line\n";
    
    # clone line into canonical isoform
    if ($opt_best_per_symbol == 0 &&
        defined($sp_to_canonical_hash{$acc_sp}))
    {
        @isoform_array = sort keys %{$sp_to_canonical_hash{$acc_sp}};
        
        foreach $isoform (@isoform_array)
        {
            # Be paranoid about canonical/non-canonical book keeping.
            # Don't clone if it was somehow also stored as non-canonical...
            if (defined($isoform) && $isoform ne '' &&
                !defined($noncanonical_to_sp_hash{$isoform}))
            {
                $line_isoform = $line;
                $line_isoform =~ s/^\S+/$isoform/;
                $line_isoform =~ s/\tSwissProt/\tCanonical/;

                print "$line_isoform\n";
            }
        }
    }
}
