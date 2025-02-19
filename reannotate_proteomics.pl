#!/usr/bin/perl -w

# 2024-08-19  delimit extra annotation fields by | instead of spaces
# 2024-08-19  strip --- from CRAPome fields prior to sorting
# 2024-02-28  annotate proteogenomics with RNA identifiers
# 2024-02-27  map more proteogenomics identifiers to extra fields
# 2024-01-16  detect CRAPomeScore and output new CRAPomeScoreMax column
# 2024-01-03  check for _SHEEP and _LYSEN contaminants (found in CRAPome)
# 2024-01-03  support additional accession column headers
# 2023-12-15  check for _BOVIN, _PIG, _GRIFR contaminants
# 2023-12-07  bugfix accession sorting in output sub-fields
# 2023-11-07  bugfix NumGenes to better handle fusions/readthroughs
# 2023-08-23  detect FragPipe reverse sequences
# 2023-08-23  flag proteogenomics mutant-only peptides
# 2023-08-23  bugfix proteogeonomics 2-character chromosomes
# 2023-08-07  more support for ENST-based FragPipe proteogenomics
# 2023-08-03  add support for ENST-based Gencode proteogenomics annotation
# 2023-06-27  update is_number() to not treat NaNs as numbers
# 2023-05-18  strip UTF-8 BOM from MaxQuant 2.4 output, which broke many things
# 2023-04-07  add equivalent accessions (RefSeq, ENSP, etc.) to output
# 2023-04-07  fix ENS* sort order
# 2023-03-17  better sorting of isoforms, fragments, Ensembl accessions
# 2022-12-23  add in additional annotation columns at beginning, if present
# 2022-08-03  Has_Contaminant_Flag, All_Reverse_Flag now [0,1] instead of [0-3]
# 2021-08-30  move amino acid column, insert new position column
# 2021-08-30  generate new ModificationID from new sorted accession order
# 2021-08-30  remove columns-to-remove code; all in another script now
# 2021-03-19  add derived Target_Species_Flag column
# 2021-03-02  support various CPTAC expression files
# 2020-06-09  strip ref|, etc. from accessions, so that reverse/contaminate tracking works again!!
#
# 2020-06-03  strip sp|, etc. from accessions, so that reverse/contaminate tracking works again!!
#             always print --- for empty annotation, rather than sometimes blanks
#             2nd update of the day: when stripping, keep uniprot part instead of swissprot, so reverse/con tracking works properly elsewhere

use Scalar::Util qw(looks_like_number);
use POSIX;

$header_remap_hash{'ProbeID'}       = 'ProbeID';
$header_remap_hash{'ProbeSet'}      = 'ProbeID';
$header_remap_hash{'PROBESET'}      = 'ProbeID';
$header_remap_hash{'GeneID'}        = 'GeneID';
$header_remap_hash{'GENEID'}        = 'GeneID';
$header_remap_hash{'Accession'}     = 'Accession';
$header_remap_hash{'ACC'}           = 'Accession';
$header_remap_hash{'Symbol'}        = 'Symbol';
$header_remap_hash{'SYMBOL'}        = 'Symbol';
$header_remap_hash{'Alias'}         = 'Alias';
$header_remap_hash{'ALIAS'}         = 'Alias';
$header_remap_hash{'Synonym'}       = 'Alias';
$header_remap_hash{'Description'}   = 'Description';
$header_remap_hash{'DESCP'}         = 'Description';
$header_remap_hash{'Annotation'}    = 'Description';

$header_remap_hash{'Accession_RNA'}      = 'Accession_RNA';
$header_remap_hash{'Accession_Protein'}  = 'Accession_Protein';
$header_remap_hash{'RefSeq'}             = 'Accession_Protein';
$header_remap_hash{'RefSeq_NT'}          = 'Accession_RNA';
$header_remap_hash{'REFSEQ_INFERRED'}    = 'Accession_Protein';
$header_remap_hash{'REFSEQ_MODEL'}       = 'Accession_Protein';
$header_remap_hash{'REFSEQ_PREDICTED'}   = 'Accession_Protein';
$header_remap_hash{'REFSEQ_PROVISIONAL'} = 'Accession_Protein';
$header_remap_hash{'REFSEQ_REVIEWED'}    = 'Accession_Protein';
$header_remap_hash{'REFSEQ_VALIDATE'}    = 'Accession_Protein';

# support CPTAC expression files
$header_remap_hash{'Acetylsite'}  = 'Accession';
$header_remap_hash{'Phosphosite'} = 'Accession';
$header_remap_hash{'NCBIGeneID'}  = 'GeneID';
$header_remap_hash{'Gene'}        = 'Symbol';


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


sub bless_delimiter_bar_only
{
    my $text = $_[0];

    #$text =~ s/\;/\|/g;
    #$text =~ s/\/\//\|/g;
    #$text =~ s/,/\|/g;
    #$text =~ s/\s+/\|/g;
    $text =~ s/\s+\|/\|/g;
    $text =~ s/\|\s+/\|/g;
    $text =~ s/\|+/\|/g;
    $text =~ s/^\|//;
    $text =~ s/\|$//;
    
    return $text;
}


sub bless_delimiter_space
{
    my $text = $_[0];

    $text =~ s/\;/\|/g;
    $text =~ s/\/\//\|/g;
    $text =~ s/,/\|/g;
    $text =~ s/\s+/\|/g;
    $text =~ s/\|+/\|/g;
    $text =~ s/^\|//;
    $text =~ s/\|$//;
    $text =~ s/\|/ /g;
    
    return $text;
}


sub score_symbol_sketchiness
{
    my @symbol_array;
    my $symbol_str = $_[0];
    my $symbol;
    my $score = -999;
    my $temp;
    
    @symbol_array = split /\|/, $symbol_str;
    
    foreach $symbol (@symbol_array)
    {
        if ($symbol =~ /^FAM\d+[A-Za-z]+/)
        {
            $temp = -1;
            if ($temp > $score) { $score = $temp; }
            next;
        }

        if ($symbol =~ /^C[XY0-9]+orf\d+/i)
        {
            $temp = -2;
            if ($temp > $score) { $score = $temp; }
            next;
        }
        
        if ($symbol =~ /^LOC\d+$/)
        {
            $temp = -3;
            if ($temp > $score) { $score = $temp; }
            next;
        }
        
        if ($symbol =~ /^RIK\d+$/)
        {
            $temp = -4;
            if ($temp > $score) { $score = $temp; }
            next;
        }

        if ($symbol =~ /^FLJ\d+$/)
        {
            $temp = -5;
            if ($temp > $score) { $score = $temp; }
            next;
        }

        if ($symbol =~ /^KIAA\d+$/ || $symbol =~ /^mKIAA\d+$/)
        {
            $temp = -6;
            if ($temp > $score) { $score = $temp; }
            next;
        }
        
        if ($symbol =~ /^MGC\d+$/)
        {
            $temp = -7;
            if ($temp > $score) { $score = $temp; }
            next;
        }
        
        if ($symbol =~ /\d+\.\d+$/)
        {
            $temp = -8;
            if ($temp > $score) { $score = $temp; }
            next;
        }

        if ($symbol =~ /[0-9]{4}/)
        {
            $temp = -9;
            if ($temp > $score) { $score = $temp; }
            next;
        }

        if ($symbol =~ /^DKFZp\d+[A-Za-z]+\d+/i)
        {
            $temp = -10;
            if ($temp > $score) { $score = $temp; }
            next;
        }
        
        if ($symbol =~ /[A-Za-z0-9]/)
        {
            $score = 0;
            last;
        }
    }
    
    return $score;
}


sub filter_description_species
{
    my $description = $_[0];

    $description =~ s/^H\.sapiens\b//;
    $description =~ s/^Homo sapiens\b//;
    $description =~ s/^PREDICTED:\s+H\.sapiens\b/PREDICTED:/;
    $description =~ s/^PREDICTED:\s+Homo sapiens\b/PREDICTED:/;
    $description =~ s/^M\.musculus\b//;
    $description =~ s/^Mus musculus\b//;
    $description =~ s/^PREDICTED:\s+M\.musculus\b/PREDICTED:/;
    $description =~ s/^PREDICTED:\s+Mus musculus\b/PREDICTED:/;
    $description =~ s/^\s+//;
    $description =~ s/\s+$//;
    
    return $description;
}


sub filter_description_general
{
    my $description = $_[0];

    $description =~ s/\s+/ /g;
    $description =~ s/_predicted//ig;
    $description =~ s/\bgene\b//ig;

    $description =~ s/\bcDNA clone.*3\'//ig;
    $description =~ s/\bcDNA clone.*5\'//ig;
    $description =~ s/\bcDNA clone\b//ig;
    $description =~ s/\bIMAGE:\d+\b//g;
    $description =~ s/\bMGC:\d+\b//g;
    $description =~ s/\(pseudogene\)//ig;
    $description =~ s/[^A-Za-z0-9]pseudogene//ig;
    $description =~ s/pseudogene//ig;
    $description =~ s/pseuodgene//ig;
    $description =~ s/hypothetical protein//ig;
    $description =~ s/[^A-Za-z0-9\-]protein\b//ig;
    $description =~ s/putative//ig;
    $description =~ s/predicted//ig;

    $description =~ s/family member//ig;
    $description =~ s/family, member//ig;
    $description =~ s/\bfamily\b//ig;
    $description =~ s/\bsubfamily\b//ig;
    $description =~ s/\bmember\b//ig;
    $description =~ s/-containing\b//ig;
    $description =~ s/\bcontaining\b//ig;
    $description =~ s/\bcontains\b//ig;
    $description =~ s/\bdomains\b//ig;
    $description =~ s/\bdomain\b//ig;
    $description =~ s/-like\b//ig;
    $description =~ s/\bweakly similar to\b//ig;
    $description =~ s/\bmoderately similar to\b//ig;
    $description =~ s/\bhighly similar to\b//ig;
    $description =~ s/\bsimilar to\b//ig;
    $description =~ s/\bsimilar\b//ig;
    $description =~ s/\bwith sequence similarity to\b//ig;
    $description =~ s/\bwith sequence similarity\b//ig;
    $description =~ s/\bwith similarity to\b//ig;
    $description =~ s/\bwith similarity\b//ig;
    $description =~ s/\bsimilarity to\b//ig;
    $description =~ s/\bsimilarity\b//ig;
    $description =~ s/\bhomolog\b//ig;
    $description =~ s/(?:^|\|\s+)Isoform of [A-Za-z0-9_.-]+[,;| ]*//ig;
    $description =~ s/[,;| ]*Isoform [A-Za-z0-9_.-]+ of\s+//ig;
    $description =~ s/\bisoform\b//ig;
    $description =~ s/\(fragment\)//ig;
    $description =~ s/\(\)//ig;

    $description =~ s/\bpartial cds\b//ig;
    $description =~ s/\bmRNA sequence\b//ig;
    $description =~ s/\bnon-coding RNA\b//ig;
    $description =~ s/\bmiscRNA\b//ig;
    $description =~ s/\bmRNA\b//ig;

    $description =~ s/[^A-Za-z0-9]//g;

    return $description;
}


# pass it gene_id
sub filter_description_accession_symbol
{
    my $description = $_[0];
    my $accession_str  = $_[1];
    my $symbol_str  = $_[2];
    my $accession;
    my $symbol;
    my $test_symbol1;
    my $test_symbol2;
    my @symbol_array;
    my @accession_array;

    @accession_array = split /\|/, $accession_str;
    # remove accessions from description
    foreach $accession (@accession_array)
    {
        $description =~ s/\b\Q$accession\E\b//ig;
    }
    
    @symbol_array = split /\|/, $symbol_str;
    # remove gene symbols from description
    foreach $symbol (@symbol_array)
    {
        $description =~ s/\b\Q$symbol\E\b//ig;
    }

    $description = filter_description_general($description);

    # check to see if filtered description is a reduced gene symbol
    foreach $symbol (@symbol_array)
    {
        # FAM proteins
        $test_symbol1 = $symbol;
        $test_symbol2 = 'FAM' . $description;
        if (lc($test_symbol1) eq lc($test_symbol2))
        {
            $description = '';
            last;
        }
        # -like
        $test_symbol2 = 'FAM' . $description . 'L';
        if (lc($test_symbol1) eq lc($test_symbol2))
        {
            $description = '';
            last;
        }
        # pseudogene
        $test_symbol2 = 'FAM' . $description;
        $test_symbol1 =~ s/p//ig;
        $test_symbol2 =~ s/p//ig;
        if (lc($test_symbol1) eq lc($test_symbol2))
        {
            $description = '';
            last;
        }

        # FAM "linked"
        $test_symbol1 = $symbol;
        $test_symbol2 = 'FAM' . $description;
        $test_symbol2 =~ s/linked//ig;
        if (lc($test_symbol1) eq lc($test_symbol2))
        {
            $description = '';
            last;
        }
        # FAM "pseudogene"
        $test_symbol2 = 'FAM' . $description;
        $test_symbol1 =~ s/p//ig;
        $test_symbol2 =~ s/p//ig;
        if (lc($test_symbol1) eq lc($test_symbol2))
        {
            $description = '';
            last;
        }

        # FAM "CC"
        $test_symbol1 = $symbol;
        $test_symbol2 = 'FAM' . $description;
        $test_symbol2 =~ s/chemokineCCmotif//ig;
        if (lc($test_symbol1) eq lc($test_symbol2))
        {
            $description = '';
            last;
        }

        # FAM "OS"
        $test_symbol1 = $symbol;
        $test_symbol2 = 'FAM' . $description;
        $test_symbol2 =~ s/oppositestrandnonproteincoding/OS/ig;
        $test_symbol2 =~ s/oppositestrand/OS/ig;
        if (lc($test_symbol1) eq lc($test_symbol2))
        {
            $description = '';
            last;
        }

        # FAM "OS", but already has FAM
        $test_symbol1 = $symbol;
        $test_symbol2 = $description;
        $test_symbol2 =~ s/oppositestrandnonproteincoding/OS/ig;
        $test_symbol2 =~ s/oppositestrand/OS/ig;
        if (lc($test_symbol1) eq lc($test_symbol2))
        {
            $description = '';
            last;
        }

        # "domain"
        $test_symbol1 = $symbol;
        $test_symbol1 =~ s/d//ig;
        $test_symbol2 = $description;
        $test_symbol2 =~ s/d//ig;
        if (lc($test_symbol1) eq lc($test_symbol2))
        {
            $description = '';
            last;
        }

        # "domain containing"
        $test_symbol1 = $symbol;
        $test_symbol1 =~ s/dc//ig;
        $test_symbol2 = $description;
        $test_symbol2 =~ s/dc//ig;
        if (lc($test_symbol1) eq lc($test_symbol2))
        {
            $description = '';
            last;
        }

        # "-like"
        $test_symbol1 = $symbol;
        $test_symbol1 =~ s/l//ig;
        $test_symbol2 = $description;
        $test_symbol2 =~ s/l//ig;
        if (lc($test_symbol1) eq lc($test_symbol2))
        {
            $description = '';
            last;
        }

        # "pseudogene"
        $test_symbol1 = $symbol;
        $test_symbol1 =~ s/p//ig;
        $test_symbol2 = $description;
        $test_symbol2 =~ s/p//ig;
        if (lc($test_symbol1) eq lc($test_symbol2))
        {
            $description = '';
            last;
        }

        # ribosomal
#        $test_symbol1 = $symbol;
#        $test_symbol2 = $description;
#        $test_symbol2 =~ s/ribosomal/RP/ig;
#        if (lc($test_symbol1) eq lc($test_symbol2))
#        {
#            $description = '';
#            last;
#        }
    }
    
    return $description;
}


# to be used with sort function
sub compare_accession
{
    my $accession_a   = $a;
    my $accession_b   = $b;
    my %hash_a;
    my %hash_b;
    my $geneid_a = '';
    my $geneid_b = '';
    my $symbol_a = '';
    my $symbol_b = '';
    my $alias_a  = '';
    my $alias_b  = '';
    my $type_a   = '';
    my $type_b   = '';
    my $description_a = '';
    my $description_b = '';
    my $value_a;
    my $value_b;
    my $temp_a;
    my $temp_b;
    my $sub_description;
    my $temp_length;
    my $debug = 0;

    if (defined($accession_annotation_hash{$accession_a}))
    {
        %hash_a = %{$accession_annotation_hash{$accession_a}};
    }
    if (defined($accession_annotation_hash{$accession_b}))
    {
        %hash_b = %{$accession_annotation_hash{$accession_b}};
    }

    if (%hash_a)
    {
        $geneid_a      = $hash_a{gene_id};
        $symbol_a      = $hash_a{symbol};
        $alias_a       = $hash_a{alias};
        $type_a        = $hash_a{type};
        $description_a = $hash_a{description};
    }
    if (%hash_b)
    {
        $geneid_b      = $hash_b{gene_id};
        $symbol_b      = $hash_b{symbol};
        $alias_b       = $hash_b{alias};
        $type_b        = $hash_b{type};
        $description_b = $hash_b{description};
    }

    if (!defined($geneid_a))       { $geneid_a = ''; }
    if (!defined($geneid_b))       { $geneid_b = ''; }
    if (!defined($symbol_a))       { $symbol_a = '~'; }
    if (!defined($symbol_b))       { $symbol_b = '~'; }
    if (!defined($type_a))         { $type_a   = ''; }
    if (!defined($type_b))         { $type_b   = ''; }
    #if (!defined($description_a)) { $description_a = '~'; }
    #if (!defined($description_b)) { $description_b = '~'; }


    # has geneid assigned
    if ($geneid_a =~ /[0-9]/ && !($geneid_b =~ /[0-9]/)) { return -1; }
    if ($geneid_b =~ /[0-9]/ && !($geneid_a =~ /[0-9]/)) { return  1; }


    # has symbol
    if ($symbol_a =~ /[A-Za-z0-9]/ && !($symbol_b =~ /[A-Za-z0-9]/)) { return -1; }
    if ($symbol_b =~ /[A-Za-z0-9]/ && !($symbol_a =~ /[A-Za-z0-9]/)) { return  1; }

    
    # best refseqs
    if (  $accession_a =~ /^NM_/  && !($accession_b =~ /^NM_/)) {return -1;}
    if (!($accession_a =~ /^NM_/) &&   $accession_b =~ /^NM_/)  {return  1;}
    if (  $accession_a =~ /^NR_/  && !($accession_b =~ /^NR_/)) {return -1;}
    if (!($accession_a =~ /^NR_/) &&   $accession_b =~ /^NR_/)  {return  1;}
    if (  $accession_a =~ /^NP_/  && !($accession_b =~ /^NP_/)) {return -1;}
    if (!($accession_a =~ /^NP_/) &&   $accession_b =~ /^NP_/)  {return  1;}

    
    # symbol sketchiness
    $value_a = score_symbol_sketchiness($symbol_a);
    $value_b = score_symbol_sketchiness($symbol_b);
    if ($value_a > $value_b) {return -1;}
    if ($value_a < $value_b) {return  1;}


    # worse refseqs
    if (  $accession_a =~ /^XM_/  && !($accession_b =~ /^XM_/)) {return -1;}
    if (!($accession_a =~ /^XM_/) &&   $accession_b =~ /^XM_/)  {return  1;}
    if (  $accession_a =~ /^XR_/  && !($accession_b =~ /^XR_/)) {return -1;}
    if (!($accession_a =~ /^XR_/) &&   $accession_b =~ /^XR_/)  {return  1;}
    if (  $accession_a =~ /^XP_/  && !($accession_b =~ /^XP_/)) {return -1;}
    if (!($accession_a =~ /^XP_/) &&   $accession_b =~ /^XP_/)  {return  1;}
    if (  $accession_a =~ /^YP_/  && !($accession_b =~ /^YP_/)) {return -1;}
    if (!($accession_a =~ /^YP_/) &&   $accession_b =~ /^YP_/)  {return  1;}
    if (  $accession_a =~ /^NG_/  && !($accession_b =~ /^NG_/)) {return -1;}
    if (!($accession_a =~ /^NG_/) &&   $accession_b =~ /^NG_/)  {return  1;}


    # Ensembl accessions
    if (  $accession_a =~ /^ENS[A-Z]{0,3}T/  && !($accession_b =~ /^ENS[A-Z]{0,3}T/)) {return -1;}
    if (!($accession_a =~ /^ENS[A-Z]{0,3}T/) &&   $accession_b =~ /^ENS[A-Z]{0,3}T/)  {return  1;}
    if (  $accession_a =~ /^ENS[A-Z]{0,3}P/  && !($accession_b =~ /^ENS[A-Z]{0,3}P/)) {return -1;}
    if (!($accession_a =~ /^ENS[A-Z]{0,3}P/) &&   $accession_b =~ /^ENS[A-Z]{0,3}P/)  {return  1;}
    if (  $accession_a =~ /^ENS[A-Z]{0,3}G/  && !($accession_b =~ /^ENS[A-Z]{0,3}G/)) {return -1;}
    if (!($accession_a =~ /^ENS[A-Z]{0,3}G/) &&   $accession_b =~ /^ENS[A-Z]{0,3}G/)  {return  1;}


    # non-NG genomic annotations are often flat out *WRONG* -- *SO WRONG*
    #  like, WTF? kind of wrong...
    $value_a = $value_b = 0;
    $value_a = ($type_a =~ /genomic/);
    $value_b = ($type_b =~ /genomic/);
    if ($value_a > $value_b) { return -1; }
    if ($value_b > $value_a) { return  1; }


    # presence of annotation
    $value_a = 0;
    if ($description_a =~ /[A-Za-z0-9]/)
    {
        $value_a = 1;
    }
    $value_b = 0;
    if ($description_b =~ /[A-Za-z0-9]/)
    {
        $value_b = 1;
    }
    if ($value_a > $value_b) {return -1;}
    if ($value_a < $value_b) {return  1;}


    # absence of pseudogene
    $value_a = 0;
    if ($type_a        =~ /pseudo/i ||
        $description_a =~ /pseudogene/i)
    {
        $value_a = 1;
    }
    $value_b = 0;
    if ($type_b        =~ /pseudo/i ||
        $description_b =~ /pseudogene/i)
    {
        $value_b = 1;
    }
    if ($value_a < $value_b) {return -1;}
    if ($value_a > $value_b) {return  1;}


    # absence of hyphen (fusions, anti-sense transcripts, etc.)
    $value_a = 0;
    if ($symbol_a =~ /\-/)
    {
        $value_a = 1;
    }
    $value_b = 0;
    if ($symbol_b =~ /\-/)
    {
        $value_b = 1;
    }
    if ($value_a < $value_b) {return -1;}
    if ($value_a > $value_b) {return  1;}


    # protein-coding in type
    $value_a = $value_b = 0;
    $value_a = ($type_a =~ /protein-coding/);
    $value_b = ($type_b =~ /protein-coding/);
    if ($value_a > $value_b) { return -1; }
    if ($value_b > $value_a) { return  1; }


    # rRNA in type
    $value_a = $value_b = 0;
    $value_a = ($type_a =~ /rRNA/);
    $value_b = ($type_b =~ /rRNA/);
    if ($value_a > $value_b) { return -1; }
    if ($value_b > $value_a) { return  1; }


    # tRNA in type
    $value_a = $value_b = 0;
    $value_a = ($type_a =~ /tRNA/);
    $value_b = ($type_b =~ /tRNA/);
    if ($value_a > $value_b) { return -1; }
    if ($value_b > $value_a) { return  1; }


    # snRNA in type
    $value_a = $value_b = 0;
    $value_a = ($type_a =~ /snRNA/);
    $value_b = ($type_b =~ /snRNA/);
    if ($value_a > $value_b) { return -1; }
    if ($value_b > $value_a) { return  1; }


    # snoRNA in type
    $value_a = $value_b = 0;
    $value_a = ($type_a =~ /snoRNA/);
    $value_b = ($type_b =~ /snoRNA/);
    if ($value_a > $value_b) { return -1; }
    if ($value_b > $value_a) { return  1; }


    # ncRNA in type
    $value_a = $value_b = 0;
    $value_a = ($type_a =~ /ncRNA/);
    $value_b = ($type_b =~ /ncRNA/);
    if ($value_a > $value_b) { return -1; }
    if ($value_b > $value_a) { return  1; }


    # other in type
    $value_a = $value_b = 0;
    $value_a = ($type_a =~ /other/);
    $value_b = ($type_b =~ /other/);
    if ($value_a > $value_b) { return -1; }
    if ($value_b > $value_a) { return  1; }


    # -AS#/-OS#/-OT/-IT# symbol ending, or antisense / opposite strand
    #  in description
    $value_a   = $value_b  = 0;
    $value_a   = ($symbol_a =~ /-AS[0-9]?\b/);
    $value_b   = ($symbol_b =~ /-AS[0-9]?\b/);
    $value_a   = ($symbol_a =~ /-OS[0-9]?\b/);
    $value_b   = ($symbol_b =~ /-OS[0-9]?\b/);
    $value_a   = ($symbol_a =~ /-OT[0-9]?\b/);
    $value_b   = ($symbol_b =~ /-OT[0-9]?\b/);
    $value_a   = ($symbol_a =~ /-IT[0-9]?\b/);
    $value_b   = ($symbol_b =~ /-IT[0-9]?\b/);
    $value_a   = ($symbol_a =~ /\bAS[0-9]?-/);
    $value_b   = ($symbol_b =~ /\bAS[0-9]?-/);
    $value_a   = ($symbol_a =~ /\bOS[0-9]?-/);
    $value_b   = ($symbol_b =~ /\bOS[0-9]?-/);
    $value_a   = ($symbol_a =~ /\bOT[0-9]?-/);
    $value_b   = ($symbol_b =~ /\bOT[0-9]?-/);
    $value_a   = ($symbol_a =~ /\bIT[0-9]?-/);
    $value_b   = ($symbol_b =~ /\bIT[0-9]?-/);
    ## can't trust Alias, so don't check them; see example NCOA7
    $value_a  += ($description_a =~ /\bantisense\b/i);
    $value_b  += ($description_b =~ /\bantisense\b/i);
    $value_a  += ($description_a =~ /\bopposite strand\b/i);
    $value_b  += ($description_b =~ /\bopposite strand\b/i);
    $value_a  += ($description_a =~ /\boverlapping transcript\b/i);
    $value_b  += ($description_b =~ /\boverlapping transcript\b/i);
    $value_a  += ($description_a =~ /\bintronic transcript\b/i);
    $value_b  += ($description_b =~ /\bintronic transcript\b/i);
    if ($value_a < $value_b) { return -1; }
    if ($value_b < $value_a) { return  1; }

    
    # alphabetical by symbol
    if ($symbol_a lt $symbol_b) {return -1;}
    if ($symbol_a gt $symbol_b) {return  1;}


    # absence of fragment
    $value_a = 0;
    if ($description_a =~ /\(fragment\)/i)
    {
        $value_a = 1;
    }
    $value_b = 0;
    if ($description_b =~ /\(fragment\)/i)
    {
        $value_b = 1;
    }
    if ($value_a < $value_b) {return -1;}
    if ($value_a > $value_b) {return  1;}


    # absence of isoform
    $value_a = 0;
    if ($description_a =~ /\bIsoform\b/i)
    {
        $value_a = 1;
    }
    $value_b = 0;
    if ($description_b =~ /\bIsoform\b/i)
    {
        $value_b = 1;
    }
    if ($value_a < $value_b) {return -1;}
    if ($value_a > $value_b) {return  1;}
    

    # annotation length, after filtering out unuseful stuff
    $temp_a = filter_description_accession_symbol($description_a,
                                                  $accession_a, $symbol_a);
    $temp_b = filter_description_accession_symbol($description_b,
                                                  $accession_b, $symbol_b);
    # measure longest sub-description
    $value_a = 0;
    foreach $description (split "\|", $temp_a)
    {
        $description =~ s/^\s+//;
        $description =~ s/\s+$//;
        
        $temp_length = length $description;
        if ($temp_length > $value_a)
        {
            $value_a = $temp_length;
        }
    }
    $value_b = 0;
    foreach $description (split "\|", $temp_b)
    {
        $description =~ s/^\s+//;
        $description =~ s/\s+$//;
        
        $temp_length = length $description;
        if ($temp_length > $value_b)
        {
            $value_b = $temp_length;
        }
    }
    if ($value_a > $value_b) {return -1;}
    if ($value_a < $value_b) {return  1;}
    

    # presence of similarity terms
    $value_a = 0;
    if ($description_a =~ /-like\b/ ||
        $description_a =~ /\bsimilar\b/ ||
        $description_a =~ /\bsimilarity\b/ ||
        $description_a =~ /\bhomolog\b/)
    {
        $value_a = 1;
    }
    $value_b = 0;
    if ($description_b =~ /-like\b/ ||
        $description_b =~ /\bsimilar\b/ ||
        $description_b =~ /\bsimilarity\b/ ||
        $description_b =~ /\bhomolog\b/)
    {
        $value_b = 1;
    }
    if ($value_a > $value_b) {return -1;}
    if ($value_a < $value_b) {return  1;}


    # presence of domain terms
    $value_a = 0;
    if ($description_a =~ /-containing\b/ ||
        $description_a =~ /\bcontaining\b/ ||
        $description_a =~ /\bcontains\b/ ||
        $description_a =~ /\bdomain\b/ ||
        $description_a =~ /\bdomains\b/)
    {
        $value_a = 1;
    }
    $value_b = 0;
    if ($description_b =~ /-containing\b/ ||
        $description_b =~ /\bcontaining\b/ ||
        $description_b =~ /\bcontains\b/ ||
        $description_b =~ /\bdomain\b/ ||
        $description_b =~ /\bdomains\b/)
    {
        $value_b = 1;
    }
    if ($value_a > $value_b) {return -1;}
    if ($value_a < $value_b) {return  1;}


    # presence of family terms
    $value_a = 0;
    if ($description_a =~ /family member\b/ ||
        $description_a =~ /family, member\b/ ||
        $description_a =~ /\bfamily\b/ ||
        $description_a =~ /\bsubfamily\b/ ||
        $description_a =~ /\bmember\b/)
    {
        $value_a = 1;
    }
    $value_b = 0;
    if ($description_b =~ /family member\b/ ||
        $description_b =~ /family, member\b/ ||
        $description_b =~ /\bfamily\b/ ||
        $description_b =~ /\bsubfamily\b/ ||
        $description_b =~ /\bmember\b/)
    {
        $value_b = 1;
    }
    if ($value_a > $value_b) {return -1;}
    if ($value_a < $value_b) {return  1;}
    

    # description, alphabetical
    if (lc($description_a) lt lc($description_b)) {return -1;}
    if (lc($description_a) gt lc($description_b)) {return  1;}
    

    # shorter symbol is better
    $value_a = length $symbol_a;
    $value_b = length $symbol_b;
    if ($value_a < $value_b) {return -1;}
    if ($value_a > $value_b) {return  1;}
    

    # presence of alias
    $value_a = 0;
    if ($alias_a =~ /[A-Za-z0-9]/)
    {
        $value_a = 1;
    }
    $value_b = 0;
    if ($alias_b =~ /[A-Za-z0-9]/)
    {
        $value_b = 1;
    }
    if ($value_a > $value_b) {return -1;}
    if ($value_a < $value_b) {return  1;}


    # sketchiness of alias
    $value_a = score_symbol_sketchiness($alias_a);
    $value_b = score_symbol_sketchiness($alias_b);
    if ($value_a > $value_b) {return -1;}
    if ($value_a < $value_b) {return  1;}
    

    # alias, alphabetical
    if ($alias_a ne '' && $alias_b ne '')
    {
        if (lc($alias_a) lt lc($alias_b)) {return -1;}
        if (lc($alias_a) gt lc($alias_b)) {return  1;}
    }
    

    # sp|tr _SPECIES accession
    $value_a = 0;
    if ($accession_a =~ /\b[^_]{3,}_[^_]{3,}\b/)
    {
        $value_a = 1;
    }
    $value_b = 0;
    if ($accession_b =~ /\b[^_]{3,}_[^_]{3,}\b/)
    {
        $value_b = 1;
    }
    if ($value_a > $value_b) {return -1;}
    if ($value_a < $value_b) {return  1;}


    # accession, alphabetical
    return ($accession_a cmp $accession_b);
}



sub print_probeid_annotation
{
    my $line_old      = $_[0];
    my $accession_str = $_[1];

    $accession_str = bless_delimiter_bar($accession_str);
    
    %seen_accession_hash = ();
    @accession_array = ();
    $num_accessions = 0;

    foreach $accession (split /\|/, $accession_str)
    {
        # uppercase accession for later looks
        $accession =~ tr/a-z/A-Z/;

        # some non-standard accessions use periods as delimiters
        @temp_array = $accession =~ /\./g;
        if (@temp_array == 1)
        {
            $accession =~ s/\.\S+//g;
        }

        if ($accession =~ /[A-Za-z0-9]/ &&
            !defined($seen_accession_hash{$accession}))
        {
            $accession_array[$num_accessions++] = $accession;
            $seen_accession_hash{$accession} = 1;
        }
    }
    
    # strip off -# if they aren't annotated, but could be after stripping
    %temp_accession_hash = ();
    foreach $accession (@accession_array)
    {
        if (!defined($accession_annotation_hash{$accession}))
        {
            # strip -# from end of accession
            $accession_stripped = $accession;
            $accession_stripped =~ s/(\d+)-\d+$/$1/g;

#if ($accession_stripped ne $accession)
#{
#printf STDERR "%s\t%s\n", $accession, $accession_stripped;
#}
            
            # we can annotate it with *something* if we strip the isoform !!
            # strip it -- better than not annotating it at all
            if (defined($accession_annotation_hash{$accession_stripped}))
            {
                $temp_accession_hash{$accession_stripped} = 1;
            }
            # otherwise, keep the isoform as-is
            else
            {
                $temp_accession_hash{$accession} = 1;
            }
        }
        else
        {
            $temp_accession_hash{$accession} = 1;
        }
    }
    @accession_array = sort keys %temp_accession_hash;

    # check for unknown accessions
    foreach $accession (@accession_array)
    {
        if (!defined($accession_annotation_hash{$accession}))
        {
            $accession_annotation_hash{$accession}{accession_rna} = '';
            $accession_annotation_hash{$accession}{accession_protein} = '';
            $accession_annotation_hash{$accession}{gene_id} = '';
            $accession_annotation_hash{$accession}{symbol} = '';
            $accession_annotation_hash{$accession}{alias} = '';
            $accession_annotation_hash{$accession}{type} = '';
            $accession_annotation_hash{$accession}{location} = '';
            $accession_annotation_hash{$accession}{description} = '';
        }
    }
    
    if ($num_accessions)
    {
        @accession_array = sort compare_accession @accession_array;
        $accession_str = join " ", @accession_array;
    }
    
    if (!($accession_str =~ /[A-Za-z0-9]/))
    {
        $accession_str = '---';
    }



    # create other fields from accession array

    $accession_rna_str = '';
    $accession_protein_str = '';
    $gene_id_str = '';
    $symbol_str = '';
    $alias_str = '';
    $type_str = '';
    $location_str = '';
    $description_str = '';
    $crapome_str = '';
    
    # initialize extra annotation strings
    foreach $col (@annotation_extra_col_array)
    {
        # skip CRAPome, since we handle it separtaely
        if (defined($accession_crapome_col) &&
            $col == $accession_crapome_col)
        {
            next;
        }
    
        $row_extra_str_by_col_hash{$col} = '';
    }

    %seen_accession_rna_hash = ();
    %seen_accession_protein_hash = ();
    %seen_gene_id_hash = ();
    %seen_symbol_hash = ();
    %seen_alias_hash = ();
    %seen_type_hash = ();
    %seen_location_hash = ();
    %seen_description_hash = ();
    %seen_extra_hash = ();
    %seen_crapome_hash = ();
    %symbol_with_geneid_hash = ();

    # flag accessions with pseudogenes
    %pseudo_hash = ();
    $count_pseudo = 0;
    $count_non_pseudo = 0;
    foreach $accession (@accession_array)
    {
        $type              = $accession_annotation_hash{$accession}{type};

        # remove NULL fields
        $type              =~ s/^null$//i;

        $type              = bless_delimiter_bar($type);
        
        $non_pseudo_flag = 0;
        $pseudo_flag     = 0;
        foreach $field (split /\|/, $type)
        {
            if ($field =~ /[A-Za-z0-9]/)
            {
                if ($field =~ /_pseudo/)
                {
                    $pseudo_flag = 1;
                }
                else
                {
                    $non_pseudo_flag = 1;
                }
            }
        }
        
        if ($pseudo_flag && $non_pseudo_flag == 0)
        {
            $pseudo_hash{$accession} = 1;
            $count_pseudo++;
        }
        else
        {
            $not_pseudo_hash{$accession} = 1;
            $count_non_pseudo++;
        }
    }


    $crapome_max = 0;
    foreach $accession (@accession_array)
    {
        # skip pseudogenes, unless there are no non-pseudogene hits
        if ($count_non_pseudo &&
            !defined($not_pseudo_hash{$accession}) &&
            defined($pseudo_hash{$accession}))
        {
            # printf STDERR "SKIP: %s\n", $accession;
            next;
        }

        $accession_rna     = $accession_annotation_hash{$accession}{accession_rna};
        $accession_protein = $accession_annotation_hash{$accession}{accession_protein};
        $gene_id           = $accession_annotation_hash{$accession}{gene_id};
        $symbol            = $accession_annotation_hash{$accession}{symbol};
        $alias             = $accession_annotation_hash{$accession}{alias};
        $type              = $accession_annotation_hash{$accession}{type};
        $location          = $accession_annotation_hash{$accession}{location};
        $description       = $accession_annotation_hash{$accession}{description};
        $crapome           = $accession_annotation_hash{$accession}{crapome};
        
#        if (!defined($accession_rna))     { $accession_rna     = ''; }
#        if (!defined($accession_protein)) { $accession_protein = ''; }
#        if (!defined($gene_id))           { $gene_id           = ''; }
#        if (!defined($symbol))            { $symbol            = ''; }
#        if (!defined($alias))             { $alias             = ''; }
#        if (!defined($type))              { $type              = ''; }
#        if (!defined($location))          { $location          = ''; }
#        if (!defined($description))       { $description       = ''; }
         if (!defined($crapome))           { $crapome           = ''; }

        # remove NULL fields
        $accession_rna     =~ s/^null$//i;
        $accession_protein =~ s/^null$//i;
        $gene_id           =~ s/^null$//i;
        $symbol            =~ s/^null$//i;
        $alias             =~ s/^null$//i;
        $type              =~ s/^null$//i;
        $location          =~ s/^null$//i;
        $description       =~ s/^null$//i;
        $crapome           =~ s/^null$//i;
        $crapome           =~ s/^\-\-\-$//i;

        $accession_rna     = bless_delimiter_bar($accession_rna);
        $accession_protein = bless_delimiter_bar($accession_protein);
        $gene_id           = bless_delimiter_bar($gene_id);
        $symbol            = bless_delimiter_bar($symbol);
        $alias             = bless_delimiter_bar($alias);
        $type              = bless_delimiter_bar($type);
        $crapome           = bless_delimiter_bar($crapome);
        
        foreach $field (split /\|/, $accession_rna)
        {
            if ($field =~ /[A-Za-z0-9]/ &&
                !defined($seen_accession_rna_hash{$field}))
            {
                if ($accession_rna_str =~ /\S/)
                {
                    $accession_rna_str .= ' ';
                }
                $accession_rna_str .= $field;

                $seen_accession_rna_hash{$field} = 1;
            }
        }

        foreach $field (split /\|/, $accession_protein)
        {
            if ($field =~ /[A-Za-z0-9]/ &&
                !defined($seen_accession_protein_hash{$field}))
            {
                if ($accession_protein_str =~ /\S/)
                {
                    $accession_protein_str .= ' ';
                }
                $accession_protein_str .= $field;

                $seen_accession_protein_hash{$field} = 1;
            }
        }

        foreach $field (split /\|/, $gene_id)
        {
            if ($field =~ /[A-Za-z0-9]/ &&
                !defined($seen_gene_id_hash{$field}))
            {
                if ($gene_id_str =~ /\S/)
                {
                    $gene_id_str .= ' ';
                }
                $gene_id_str .= $field;

                $seen_gene_id_hash{$field} = 1;
            }
        }

        foreach $field (split /\|/, $symbol)
        {
            if ($field =~ /[A-Za-z0-9]/ &&
                !defined($seen_symbol_hash{$field}))
            {
                if ($symbol_str =~ /\S/)
                {
                    $symbol_str .= ' ';
                }
                $symbol_str .= $field;

                $seen_symbol_hash{$field} = 1;
                
                if ($gene_id =~ /[0-9]/)
                {
                    $symbol_with_geneid_hash{$field} = 1;
                }
            }
        }

        foreach $field (split /\|/, $alias)
        {
            if ($field =~ /[A-Za-z0-9]/ &&
                !defined($seen_alias_hash{$field}))
            {
                if ($alias_str =~ /\S/)
                {
                    $alias_str .= ' ';
                }
                $alias_str .= $field;

                $seen_alias_hash{$field} = 1;
            }
        }

        foreach $field (split /\|/, $type)
        {
            if ($field =~ /[A-Za-z0-9]/ &&
                !defined($seen_type_hash{$field}))
            {
                if ($type_str =~ /\S/)
                {
                    $type_str .= '|';
                }
                $type_str .= $field;

                $seen_type_hash{$field} = 1;
            }
        }

        foreach $field (split /\!/, $location)
        {
            if ($field =~ /[A-Za-z0-9]/ &&
                !defined($seen_location_hash{$field}))
            {
                if ($location_str =~ /\S/)
                {
                    $location_str .= '!';
                }
                $location_str .= $field;

                $seen_location_hash{$field} = 1;
            }
        }

        foreach $field (split /\|/, $description)
        {
            if ($field =~ /[A-Za-z0-9]/ &&
                !defined($seen_description_hash{$field}))
            {
                if ($description_str =~ /\S/)
                {
                    $description_str .= '|';
                }
                $description_str .= $field;

                $seen_description_hash{$field} = 1;
            }
        }

        # sort CRAPome sub-fields in reverse order
        @temp_array = sort { $b <=> $a } split /\|/, $crapome;
        foreach $field (@temp_array)
        {
            if ($field =~ /[A-Za-z0-9]/ &&
                is_number($field) &&
                $field > 0 &&
                !defined($seen_crapome_hash{$field}))
            {
                if ($crapome_str =~ /\S/)
                {
                    $crapome_str .= '|';
                }
                $crapome_str .= $field;
                
                if ($field > $crapome_max)
                {
                    $crapome_max = $field;
                }

                $seen_crapome_hash{$field} = 1;
            }
        }
        
        # handle extra annotation columns
        foreach $col (@annotation_extra_col_array)
        {
            # skip CRAPome, since we handle it separately
            if (defined($accession_crapome_col) &&
                $col == $accession_crapome_col)
            {
                next;
            }
        
            $value_str  = $row_extra_str_by_col_hash{$col};
            $value_orig = $accession_annotation_extra_hash{$accession}{$col};
            
            if (!defined($value_orig))
            {
                $value_orig = '';
            }

            $value_orig =~ s/^null$//i;
            $value_orig =  bless_delimiter_bar_only($value_orig);

            foreach $field (split /\|/, $value_orig)
            {
                if ($field =~ /\S/ &&
                    !(defined($seen_extra_hash{$col}) &&
                      defined($seen_extra_hash{$col}{$field})))
                {
                    if ($value_str =~ /\S/)
                    {
                        $value_str .= '|';
                    }
                    $value_str .= $field;

                    $seen_extra_hash{$col}{$field} = 1;
                }
            }
            
            # store extra annotation strings as we grow them
            $row_extra_str_by_col_hash{$col} = $value_str;
        }
    }

    
    if (!($accession_rna_str =~ /[A-Za-z0-9]/))
    {
        $accession_rna_str = '---';
    }
    if (!($accession_protein_str =~ /[A-Za-z0-9]/))
    {
        $accession_protein_str = '---';
    }
    if (!($gene_id_str =~ /[A-Za-z0-9]/))
    {
        $gene_id_str = '---';
    }
    if (!($symbol_str =~ /[A-Za-z0-9]/))
    {
        $symbol_str = '---';
    }
    if (!($alias_str =~ /[A-Za-z0-9]/))
    {
        $alias_str = '---';
    }
    if (!($description_str =~ /[A-Za-z0-9]/))
    {
        $description_str = '---';
    }
    if (!($location_str =~ /[A-Za-z0-9]/))
    {
        $location_str = '---';
    }
    if (!($type_str =~ /[A-Za-z0-9]/))
    {
        $type_str = '---';
    }
    if (!($crapome_str =~ /[A-Za-z0-9]/))
    {
        $crapome_str = '---';
    }

    $accession_rna_str     = bless_delimiter_space($accession_rna_str);
    $accession_protein_str = bless_delimiter_space($accession_protein_str);
    $gene_id_str           = bless_delimiter_space($gene_id_str);
    $symbol_str            = bless_delimiter_space($symbol_str);
    $alias_str             = bless_delimiter_space($alias_str);
    $type_str              = bless_delimiter_bar($type_str);
    $type_str              =~ s/\;/\|/g;
    $crapome_str           = bless_delimiter_space($crapome_str);

    # clean up extra annotation
    foreach $col (@annotation_extra_col_array)
    {
        # skip CRAPome, since we handle it separately
        if (defined($accession_crapome_col) &&
            $col == $accession_crapome_col)
        {
            next;
        }
    
        $value_str = $row_extra_str_by_col_hash{$col};
    
        if (!($value_str =~ /\S/))
        {
            $value_str = '---';
        }

        $value_str = bless_delimiter_bar_only($value_str);

        $row_extra_str_by_col_hash{$col} = $value_str;
    }

    # remove symbols without geneids if we already have some
    # no need to sort the temp array, we're not using its values here
    @temp_array = keys %symbol_with_geneid_hash;
    if (@temp_array)
    {
        $symbol_str_new = '';

        @symbol_array = split / /, $symbol_str;
        foreach $symbol (@symbol_array)
        {
            if (defined($symbol_with_geneid_hash{$symbol}))
            {
                if ($symbol_str_new ne '')
                {
                    $symbol_str_new .= ' ';
                }
                $symbol_str_new .= $symbol;
            }
            else
            {
                delete $seen_symbol_hash{$symbol};
            }
        }
        
        $symbol_str = $symbol_str_new;
    }


    # pad | delimited fields
    $description_str =~ s/\s+\|/\|/g;
    $description_str =~ s/\|\s+/\|/g;
    $description_str =~ s/\|+/\|/g;
    $description_str =~ s/\|/ \|\|\| /g;
    
    # remove +2 from all_reverse, if there winds up being one that isn't REV__
    if ($all_reverse >= 2 &&
        ($accession_rna_str     =~ /[A-Za-z0-9]/ ||
         $accession_protein_str =~ /[A-Za-z0-9]/))
    {
        $all_reverse -= 2;
    }

    # count symbols, ignoring any fusions that are already otherwise present
    $count_symbols = 0;
    if ($symbol_str =~ /[A-Za-z0-9]/)
    {
        %fusion_hash     = ();
        %non_fusion_hash = ();
        @symbol_array = split / /, $symbol_str;
        
        foreach $field (@symbol_array)
        {
            # potential fusion gene
            # some can have multiple fusions and/or hyphens, check both ends
            #
            # hyphen in individual gene: ERVK3-1 ZNF8-ERVK3-1
            # was missed before: CD302 LY75-CD302
            # not fusions: AS-, -AS, -AS#, -#, -OT#, -DT
            #
            # problems: ABCF2-H2BE1, C1QTNF3-AMACR, C8orf44-SGK3,
            #           CCZ1P-OR7E38P
            # weird:    CH17-340M24.3, CH507-42P11.6, CTB-178M22.2,
            #           CTB-1I21.1, CTB-49A3.2, CTC-338M12.4, CTD-3080P12.3
            #
            # XXX, XX#(X|#), X#(X|#)
            #
            # must have hyphen in middle, tolerate up to 2 alphanum in middle
            #
            # RB1 is a real symbol, so we may miss some potential real fusions
            #
            if ($field =~
                /\b[A-Za-z]([A-Za-z]{2}|[A-Za-z]?[0-9][A-Za-z0-9])[^-]*\-([A-Za-z0-9]{1,2}[0-9]?\-)*[A-Za-z]([A-Za-z]{2}|[A-Za-z]?[0-9][A-Za-z0-9])[^-]*\b/)
            {
                $fusion_hash{$field} = 1;
            }
            # non_fusion gene, count it
            else
            {
                $non_fusion_hash{$field} = 1;
                $count_symbols++;
            }
        }
        
        # decide whether to add potential fusion genes to counts
        foreach $fusion (keys %fusion_hash)
        {
            $fusion_flag = 0;
            
            # fusion that includes a non-fusion gene
            # we already counted the non-fusion, so don't count the fusion
            foreach $field (keys %non_fusion_hash)
            {
                # match within the fusion, make sure it passes the same
                # starting character rules as were used in the fusion matching
                if ($fusion =~ /\b\Q$field\E\b/ &&
                    $fusion =~ /^[A-Za-z]([A-Za-z]{2}|[A-Za-z]?[0-9][A-Za-z0-9])/)
                {
                    $fusion_flag = 1;
                    last;
                }
            }

            # didn't match any non-fusions, count it
            if ($fusion_flag == 0)
            {
                $count_symbols++;
            }
        }
    }
    
    $warning_level = 0;			# all is well
    if ($count_symbols < 1)		# no hits
    {
        $warning_level = 1;
    }
    if ($count_symbols > 1)		# multiple genes
    {
        $warning_level = 5;
    }
    
    # assign a warning level of 2 if the genes are close to each other
    # 3 if on same arm
    # 4 if on same chromosome
    # 5 if different chromosomes
    if ($warning_level == 5)
    {
        $location_str =~ s/\|/\!/g;
        $location_str =~ s/\!+/\!/g;
        @location_array = split /\!/, $location_str;
        
        # not just a single locus
        # the following will mostly only work for Human
        #    Mouse would require alphabetical ranges... not going there yet...
        # some Mouse may still work well enough anyways?
        if (@location_array > 1)
        {
            $num_locations = @location_array;
            $orig_num_locations = $num_locations;

            # strip sub-bands
            for ($i = 0; $i < $orig_num_locations; $i++)
            {
                $location_array[$i] =~ s/\.[0-9]+//g;
            }

            # handle ranges
            for ($i = 0; $i < $orig_num_locations; $i++)
            {
                $field = $location_array[$i];
                
                if ($field =~ /^([^-]+)\-([^-]+)$/)
                {
                    $location1 = $1;
                    $location2 = $2;
                    
                    # hack to deal with centromeres
                    $location1 =~ s/cen/cen0/g;
                    $location2 =~ s/cen/cen0/g;
                    
                    # hack to deal with terminus
                    $location1 =~ s/([a-z]+)ter/\Q$1\E999/g;
                    $location2 =~ s/([a-z]+)ter/\Q$1\E999/g;
                    
#                    printf STDERR "DEBUG RANGE0\t%s\t%s\t%s\n",
#                                   $field, $location1, $location2;
                    
                    # mostly Human-specific nomenclature here
                    # TO-DO: handle arms: cen, pter, qter
                    if ($location1 =~ /^([A-Z0-9]+)([^A-Z0-9]+)([0-9]+)$/)
                    {
                        $chromosome  = $1;
                        $arm1        = $2;
                        $band1       = $3;
                        
#                        printf STDERR "DEBUG RANGE1\t%s\t%s\t%s\n",
#                                   $chromosome, $arm1, $band1;

                        if ($location2 =~ /^([^A-Z0-9]+)([0-9]+)$/)
                        {
                            $arm2        = $1;
                            $band2       = $2;

                            # replace centromere with 0 on same arm
                            if ($arm1 =~ /^cen/)
                            {
                                $arm1  = $arm2;
                                $band1 = 0;
                            }
                            # replace centromere with 0 on same arm
                            if ($arm2 =~ /^cen/)
                            {
                                $arm2  = $arm1;
                                $band2 = 0;
                            }
                            
#                            printf STDERR "DEBUG RANGE2\t%s\t%s\t%s\n",
#                                   $chromosome, $arm2, $band2;
                            
                            # on the same arm
                            if ($arm1 eq $arm2 &&
                                is_number($band1) && is_number($band2))
                            {
                                # same band, no band range
                                if ($band1 eq $band2)
                                {
                                    $new_location = sprintf "%s%s%d",
                                        $chromosome, $arm1, $band1;

                                    $location_array[$i] = $new_location;

                                    next;
                                }
                            
                                $start = $band1;
                                $end   = $band2;
                                
                                if ($band1 > $band2)
                                {
                                    $start = $band2;
                                    $end   = $band1;
                                }
                                
                                for ($j = $start; $j <= $end; $j++)
                                {
                                    $new_location = sprintf "%s%s%d",
                                        $chromosome, $arm1, $j;

#                                    printf STDERR "DEBUG RANGE3\t%s\n",
#                                        $new_location;

                                    # overwrite original location
                                    if ($j == $start)
                                    {
                                        $location_array[$i] = $new_location;
                                    }
                                    # add new location
                                    else
                                    {
                                        $location_array[$num_locations++] =
                                            $new_location;
                                    }
                                }
                            }
                        }
                    }
                }
            }


            # generate hash of final locations
            %location_hash = ();
            for ($i = 0; $i < $num_locations; $i++)
            {
                $location_hash{$location_array[$i]} = 1;

#                printf STDERR "DEBUG FIELD\t%s\n",
#                               $location_array[$i];
            }
            
            @new_location_array = sort keys %location_hash;
            $num_locations = @new_location_array;
            $new_location_str   = join '!', @new_location_array;

            # all genes within same band
            if (@new_location_array == 1)
            {
                $warning_level = 2;
            }
            # multiple locations, check to see if all are contiguous
            else
            {
#               printf STDERR "DEBUG OLD_LOC\t%s\n", $location_str;
#               printf STDERR "DEBUG TMP_LOC\t%s\n", $new_location_str;

                $start = 9E99;
                $end   = -9E99;
                
                %chromosome_hash = ();
                %arm_hash        = ();
                %band_hash       = ();
                for ($i = 0; $i < $num_locations; $i++)
                {
                    $field = $new_location_array[$i];
                
                    # mostly Human-specific nomenclature here
                    if ($field =~ /^([A-Z0-9]+)([^A-Z0-9]+)$/)
                    {
                        $chromosome = $1;
                        $arm        = $2;

                        $chromosome_hash{$chromosome} = 1;
                        $arm_hash{$arm}               = 1;
                    }
                    if ($field =~ /^([A-Z0-9]+)([^A-Z0-9]+)([0-9]+)$/)
                    {
                        $chromosome = $1;
                        $arm        = $2;
                        $band       = $3;
                        
                        if ($band < $start) { $start = $band; }
                        if ($band > $end)   { $end   = $band; }
                        
                        $chromosome_hash{$chromosome} = 1;
                        $arm_hash{$arm}               = 1;
                        $band_hash{$band}             = 1;

#                        printf STDERR "DEBUG BAND %s %s %s %s\n", $field, $band, $start, $end;
                    }
                }
                
                @chromosome_array = sort keys %chromosome_hash;
                @arm_array        = sort keys %arm_hash;
                @band_array       = sort {$a <=> $b} keys %band_hash;
                $num_bands = @band_array;

                if (@chromosome_array > 1)
                {
                    $warning_level = 5;
                }
                elsif (@arm_array > 1)
                {
                    $warning_level = 4;
                }
                # we don't count an arm without a band as being part of range
                elsif ($num_bands != $num_locations)
                {
                    $warning_level = 3;
                }
                else
                {
                    $warning_level = 3;
                    
                    $max_dist = 0;
                    for ($i = 1; $i < $num_bands; $i++)
                    {
                        $dist = $band_array[$i] - $band_array[$i-1];
                        if ($dist > $max_dist) { $max_dist = $dist; }
                    }
                
                    # allow for slop of +/- 2, due to loc annotation issues
                    if ($max_dist <= 2)
                    {
                        $warning_level = 2;
                        
                        $new_location_str = sprintf "FOOBAR %s%s%d-%d",
                            $chromosome_array[0], $arm_array[0],
                            $start, $end;

#                        printf STDERR "DEBUG OLD_LOC\t%s\n", $location_str;
#                        printf STDERR "DEBUG NEW_LOC\t%s\t%d\t%d\n", $new_location_str, $width, $num_bands;
                    }
                }
            }
            
        }
        # single locus, genes are clearly near each other
        else
        {
            $warning_level = 2;
        }
    }
    
    if ($warning_level == 0) { $stats_hash{0} += 1; }
    if ($warning_level == 1) { $stats_hash{1} += 1; }
    if ($warning_level == 2) { $stats_hash{2} += 1; }
    if ($warning_level == 3) { $stats_hash{3} += 1; }
    if ($warning_level == 4) { $stats_hash{4} += 1; }
    if ($warning_level == 5) { $stats_hash{5} += 1; }
    
    
    # prune expanded accessions back down to original accessions
    %accession_pruned_hash = ();
    $accession_pruned_str = '---';
    foreach $accession (keys %seen_accession_protein_hash)
    {
        if (defined($seen_accession_hash{$accession}))
        {
            $accession_pruned_hash{$accession} = 1;
        }
    }
    @accession_pruned_array =
        sort compare_accession keys %accession_pruned_hash;

    # tack on additional accessions annotated as equivalent
    # must keep these separate, for the modification site reordering below
    @accession_less_pruned_array = @accession_pruned_array;
    $count_pruned                = @accession_pruned_array;
    foreach $accession (sort compare_accession
                        keys %seen_accession_protein_hash)
    {
        if (defined($seen_accession_hash{$accession}) ||
            defined($accession_annotation_hash{$accession}))
        {
            if (!defined($accession_pruned_hash{$accession}))
            {
                $accession_less_pruned_array[$count_pruned++] = $accession;

                $accession_pruned_hash{$accession} = 1;
            }
        }
    }


    # this is just for printing, so use the less pruned array
    if (@accession_less_pruned_array)
    {
        $accession_pruned_str = join ' ', @accession_less_pruned_array;
    }


    if (!($accession_rna_str =~ /[A-Za-z0-9]/))
    {
        $accession_rna_str = '---';
    }
    if (!($accession_protein_str =~ /[A-Za-z0-9]/))
    {
        $accession_protein_str = '---';
    }
    if (!($gene_id_str =~ /[A-Za-z0-9]/))
    {
        $gene_id_str = '---';
    }
    if (!($symbol_str =~ /[A-Za-z0-9]/))
    {
        $symbol_str = '---';
    }
    if (!($alias_str =~ /[A-Za-z0-9]/))
    {
        $alias_str = '---';
    }
    if (!($description_str =~ /[A-Za-z0-9]/))
    {
        $description_str = '---';
    }
    if (!($location_str =~ /[A-Za-z0-9]/))
    {
        $location_str = '---';
    }
    if (!($type_str =~ /[A-Za-z0-9]/))
    {
        $type_str = '---';
    }


    # flag whether it is likely not a real target species protein or not
    $target_species_flag = 1;
    if ($all_reverse ||
        ($has_contam && !($gene_id_str =~ /[0-9]/)))
    {
        $target_species_flag = 0;
    }


    # re-order ModificationID fields to match new sorted accession order
    $mod_pos_str = '';
    if (defined($modificationid_col))
    {
        @array = split /\t/, $line_old;
    
        $modificationid = $array[$modificationid_col];
        %temp_seen_hash = ();
        @index_array    = ();
        $count_new      = 0;

        # mod_acc1;acc2:pos1;pos2
        if ($modificationid =~ /([^_]+)_([^:]+):(.*)/)
        {
            $mod_type    = $1;
            $mod_acc_str = $2;
            $mod_pos_str = $3;
            
            @accession_mod_array = split ';', $mod_acc_str;
            @position_mod_array  = split ';', $mod_pos_str;

            # scan through new accession order
            foreach $accession_new (@accession_pruned_array)
            {
                # scan ModificationID accessions, add it to new list
                for ($i = 0; $i < @accession_mod_array; $i++)
                {
                    $accession_mod = $accession_mod_array[$i];

                    # exact match, assume we've already cleaned the versions, etc.
                    if ($accession_mod eq $accession_new)
                    {
                        $index_array[$count_new] = $i;
                        $count_new++;
                        
                        $temp_seen_hash{$i} = 1;

                        last;
                    }
                }
            }

            # add in remaining ModificationID accessions
            # this should be mostly CON_ and REV_
            for ($i = 0; $i < @accession_mod_array; $i++)
            {
                if (!defined($temp_seen_hash{$i}))
                {
                    $index_array[$count_new] = $i;
                    $count_new++;

                    $temp_seen_hash{$i} = 1;
                }
            }
            
            # generate new accession portion
            $modificationid_new = $mod_type . '_';
            for ($i = 0; $i < $count_new; $i++)
            {
                $index         = $index_array[$i];
                $accession_mod = $accession_mod_array[$index];
                
                if ($i)
                {
                    $modificationid_new .= ';';
                }
                
                $modificationid_new .= $accession_mod;
            }

            # generate new positions portion
            $mod_pos_str = '';
            for ($i = 0; $i < $count_new; $i++)
            {
                $index        = $index_array[$i];
                $position_mod = $position_mod_array[$index];
                
                if ($i)
                {
                    $mod_pos_str .= ';';
                }
                
                $mod_pos_str .= $position_mod;
            }
            $modificationid_new .= ':' . $mod_pos_str;
            
            
            # overwrite old ModificationID with new
            if ($modificationid ne $modificationid_new)
            {
                $array[$modificationid_col] = $modificationid_new;
                $line_old = join "\t", @array;
            
                # printf STDERR "FOOBAR\t%s\t%s\n",
                #    $modificationid, $modificationid_new;
            }
        }
    }
    
    # remove amino acid column, since we're going to move it
    if (defined($amino_acid_col))
    {
        @array     = split /\t/, $line_old;
        @array_new = ();
        
        $j = 0;
        for ($i = 0; $i < @array; $i++)
        {
            if ($i != $amino_acid_col)
            {
                $array_new[$j++] = $array[$i];
            }
        }
        $line_old = join "\t", @array_new;
    }

    # floor flags to 1
    # the extra level of detail from 1-3 just adds confusion
    if ($has_contam  > 1) { $has_contam  = 1; }
    if ($all_reverse > 1) { $all_reverse = 1; }

    $line_new  = "$has_contam";
    $line_new .= "\t$all_reverse";
    $line_new .= "\t$cow_contam";
    $line_new .= "\t$target_species_flag";

    # insert proteogenomics column(s)
    if ($proteogenomics_flag)
    {
        # unique to mutant
        if ($mutant_flag && $ensp_flag == 0)
        {
            $line_new .= "\t" . '1';
        }
        # not unique to mutant, leave blank
        else
        {
            $line_new .= "\t" . '0';
        }
    }

    $line_new .= "\t$accession_pruned_str";

    # insert CRAPome fields
    if (defined($accession_crapome_col))
    {
        $line_new .= "\t$crapome_str";
        $line_new .= "\t$crapome_max";
    }

    # insert proteogenomics column(s)
    if ($proteogenomics_flag)
    {
        $line_new .= "\t$accession_rna_str";
    }

    # insert extra annotation fields
    foreach $col (@annotation_extra_col_array)
    {
        # skip CRAPome, since we handle it separately
        if (defined($accession_crapome_col) &&
            $col == $accession_crapome_col)
        {
            next;
        }
    
        $value_str = $row_extra_str_by_col_hash{$col};
        $line_new .= "\t$value_str";
    }

    $line_new .= "\t$gene_id_str";
    $line_new .= "\t$count_symbols";
    $line_new .= "\t$warning_level";
    $line_new .= "\t$symbol_str";
    if (defined($amino_acid_col))
    {
        $line_new .= "\t" . $array[$amino_acid_col];
    }
    if (defined($modificationid_col))
    {
        $line_new .= "\t" . $mod_pos_str;
    }
    $line_new .= "\t$alias_str";
    $line_new .= "\t$type_str";
    $line_new .= "\t$location_str";
    $line_new .= "\t$description_str";

    if (1 || $line_new =~ /[A-Za-z0-9]/)
    {
        $line_new = sprintf "%s\t%s", $line_old, $line_new;
        printf "%s\n", $line_new;
    }
}




# Begin main()

$strip_flag          = 0;
$proteogenomics_flag = 0;

$num_files = 0;
for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        if ($field eq '--strip' || $field eq '--prune')
        {
            $strip_flag = 1;
        }
        # only the data content can tell us if it is proteogenomics or not,
        # so we need a flag to tell us ahead of time for inserting columns
        elsif ($field eq '--proteogenomics')
        {
            $proteogenomics_flag = 1;
        }
        else
        {
            printf "ABORT -- unknown option %s\n", $field;
            printf "Syntax: program accession_annotation.txt data.txt\n";
            exit(1);
        }
    }
    else
    {
        if ($num_files == 0)
        {
            $accession_annotation_filename = $field;
            $num_files++;
        }
        elsif ($num_files == 1)
        {
            $data_filename = $field;
            $num_files++;
        }
    }
}


if (!defined($data_filename) || !defined($accession_annotation_filename))
{
    printf "Syntax: program accession_annotation.txt data.txt\n";
    exit(1);
}

open ANNOTATION, "$accession_annotation_filename" or die "can't open $accession_annotation_filename";
open DATA, "$data_filename" or die "can't open $data_filename";


# accession annotation header line
$line = <ANNOTATION>;
$line =~ s/[\r\n]+//;

# remove UTF-8 byte order mark, since it causes many problems
# remove some other BOM I've observed in the wild
$line =~ s/(^\xEF\xBB\xBF|^\xEF\x3E\x3E\xBF)//;

@array = split /\s+/, $line;
for ($i = 0; $i < @array; $i++)
{
    $field = $array[$i];
    $field =~ s/^\s+//;
    $field =~ s/\s+$//;
    
    if ($field =~ /^\(/)
    {
        last;
    }

    if (defined($header_remap_hash{$field}))
    {
        $field = $header_remap_hash{$field};
    }
    
    $accession_header_hash{$field} = $i;
    
    if ($field =~ /[A-Za-z0-9]/)
    {
        $annotation_col_to_header_hash{$i} = $field;
    }
}

$accession_accession_col         = $accession_header_hash{'Accession'};
$accession_accession_rna_col     = $accession_header_hash{'Accession_RNA'};
$accession_accession_protein_col = $accession_header_hash{'Accession_Protein'};
$accession_gene_id_col           = $accession_header_hash{'GeneID'};
$accession_symbol_col            = $accession_header_hash{'Symbol'};
$accession_alias_col             = $accession_header_hash{'Alias'};
$accession_type_col              = $accession_header_hash{'Type'};
$accession_location_col          = $accession_header_hash{'Location'};
$accession_description_col       = $accession_header_hash{'Description'};
$accession_crapome_col           = $accession_header_hash{'CRAPomeScore'};


# support merged Gencode ENST annotation
$map_with_additional_columns_flag = 0;
if ($line =~ /^ensembl_transcript_id/ &&
    $line =~ /type_transcript/)
{
    $accession_accession_col         = $accession_header_hash{'ensembl_transcript_id'};
    $accession_accession_rna_col     = $accession_header_hash{'ensembl_transcript_id'};
    $accession_accession_protein_col = $accession_header_hash{'protein_id'};
    $accession_gene_id_col           = $accession_header_hash{'entrez_gene_id'};
    $accession_symbol_col            = $accession_header_hash{'symbol'};
    $accession_alias_col             = $accession_header_hash{'alias'};
    $accession_type_col              = $accession_header_hash{'type_gene'};
    $accession_location_col          = $accession_header_hash{'location'};
    $accession_description_col       = $accession_header_hash{'description'};

    $map_with_additional_columns_flag = 1;
}


if (!defined($accession_accession_col))
{
    printf STDERR "ABORT -- missing accession Accession column\n";
    exit(1);
}
if (!defined($accession_accession_rna_col))
{
    printf STDERR "ABORT -- missing accession Accession_RNA column\n";
    exit(1);
}
if (!defined($accession_accession_protein_col))
{
    printf STDERR "ABORT -- missing accession Accession_Protein column\n";
    exit(1);
}
if (!defined($accession_gene_id_col))
{
    printf STDERR "ABORT -- missing accession GeneID column\n";
    exit(1);
}
if (!defined($accession_symbol_col))
{
    printf STDERR "ABORT -- missing accession Symbol column\n";
    exit(1);
}
if (!defined($accession_alias_col))
{
    printf STDERR "ABORT -- missing accession Alias column\n";
    exit(1);
}
if (!defined($accession_type_col))
{
    printf STDERR "ABORT -- missing accession Type column\n";
    exit(1);
}
if (!defined($accession_location_col))
{
    printf STDERR "ABORT -- missing accession Location column\n";
    exit(1);
}
if (!defined($accession_description_col))
{
    printf STDERR "ABORT -- missing accession Description column\n";
    exit(1);
}


# flag regular annotation columns
$annotation_regular_col_hash{$accession_accession_col} = 1;
$annotation_regular_col_hash{$accession_accession_rna_col} = 1;
$annotation_regular_col_hash{$accession_accession_protein_col} = 1;
$annotation_regular_col_hash{$accession_gene_id_col} = 1;
$annotation_regular_col_hash{$accession_symbol_col} = 1;
$annotation_regular_col_hash{$accession_alias_col} = 1;
$annotation_regular_col_hash{$accession_type_col} = 1;
$annotation_regular_col_hash{$accession_location_col} = 1;
$annotation_regular_col_hash{$accession_description_col} = 1;

# detect extra annotation columns that we didn't originally support
%annotation_extra_col_hash = ();
foreach $col (sort {$a<=>$b} keys %annotation_col_to_header_hash)
{
    if (defined($annotation_regular_col_hash{$col}))
    {
        next;
    }

    $header = $annotation_col_to_header_hash{$col};
    $annotation_extra_col_hash{$col} = $header;
}
@annotation_extra_col_array = sort {$a<=>$b} keys %annotation_extra_col_hash;
$i = 0;
foreach $col (@annotation_extra_col_array)
{
    $annotation_extra_header_array[$i++] = $annotation_extra_col_hash{$col};
}


while(defined($line=<ANNOTATION>))
{
    $line =~ s/[\r\n]+//;
    @array = split /\t/, $line, -1;   # don't auto-strip empty descriptions...
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
    }
    
    $accession           = $array[$accession_accession_col];
    $accession_rna       = $array[$accession_accession_rna_col];
    $accession_protein   = $array[$accession_accession_protein_col];
    $gene_id             = $array[$accession_gene_id_col];
    $symbol              = $array[$accession_symbol_col];
    $alias               = $array[$accession_alias_col];
    $type                = $array[$accession_type_col];
    $location            = $array[$accession_location_col];
    $description         = $array[$accession_description_col];

    $crapome             = '';
    if (defined($accession_crapome_col))
    {
        $crapome         = $array[$accession_crapome_col];
    }

    # some non-standard accessions use periods as delimiters
    @temp_array = $accession =~ /\./g;
    if (@temp_array == 1)
    {
        $accession =~ s/\.\S+//g;
    }
    
    # uppercase accession for later lookup
    $accession =~ tr/a-z/A-Z/;

    $accession_rna     = bless_delimiter_bar($accession_rna);
    $accession_protein = bless_delimiter_bar($accession_protein);
    $gene_id           = bless_delimiter_bar($gene_id);
    $symbol            = bless_delimiter_bar($symbol);
    $alias             = bless_delimiter_bar($alias);
    $type              = bless_delimiter_bar($type);
    $crapome           = bless_delimiter_bar($crapome);

    # don't load the annotation if it has no geneid
    # since some isoforms aren't annotated right,
    # and having present -- but blank -- annotation will screw up rules later
    #
    # *** I'm not sure if this is true anymore, and I need to keep everything
    # *** due to HERV annotations now.
#    if (!($gene_id =~ /[0-9]/))
#    {
#        next;
#    }

    # strip rna/prot distinction from type
    $type =~ s/rna_//g;
    $type =~ s/prot_//g;

    $accession_annotation_hash{$accession}{accession_rna}     = $accession_rna;
    $accession_annotation_hash{$accession}{accession_protein} = $accession_protein;
    $accession_annotation_hash{$accession}{gene_id}           = $gene_id;
    $accession_annotation_hash{$accession}{symbol}            = $symbol;
    $accession_annotation_hash{$accession}{alias}             = $alias;
    $accession_annotation_hash{$accession}{type}              = $type;
    $accession_annotation_hash{$accession}{location}          = $location;
    $accession_annotation_hash{$accession}{description}       = $description;
    
    if (defined($accession_crapome_col))
    {
        $accession_annotation_hash{$accession}{crapome}           = $crapome;
    }

    # store RNA/Protein lookups as well, in case we need them
    if ($map_with_additional_columns_flag)
    {
      @accession_array = split /\|/, $accession_rna;
      foreach $accession_temp (@accession_array)
      {
        if ($accession_temp =~ /[A-Za-z0-9]/ &&
            !($accession_temp =~ /^N\/*A$/i))
        {
          $accession_annotation_hash{$accession_temp}{accession_rna}     = $accession_rna;
          $accession_annotation_hash{$accession_temp}{accession_protein} = $accession_protein;
          $accession_annotation_hash{$accession_temp}{gene_id}           = $gene_id;
          $accession_annotation_hash{$accession_temp}{symbol}            = $symbol;
          $accession_annotation_hash{$accession_temp}{alias}             = $alias;
          $accession_annotation_hash{$accession_temp}{type}              = $type;
          $accession_annotation_hash{$accession_temp}{location}          = $location;
          $accession_annotation_hash{$accession_temp}{description}       = $description;

          if (defined($accession_crapome_col))
          {
            $accession_annotation_hash{$accession}{crapome}              = $crapome;
          }
        }

      }
      @accession_array = split /\|/, $accession_protein;
      foreach $accession_temp (@accession_array)
      {
        if ($accession_temp =~ /[A-Za-z0-9]/ &&
            !($accession_temp =~ /^N\/*A$/i))
        {
          $accession_annotation_hash{$accession_temp}{accession_rna}     = $accession_rna;
          $accession_annotation_hash{$accession_temp}{accession_protein} = $accession_protein;
          $accession_annotation_hash{$accession_temp}{gene_id}           = $gene_id;
          $accession_annotation_hash{$accession_temp}{symbol}            = $symbol;
          $accession_annotation_hash{$accession_temp}{alias}             = $alias;
          $accession_annotation_hash{$accession_temp}{type}              = $type;
          $accession_annotation_hash{$accession_temp}{location}          = $location;
          $accession_annotation_hash{$accession_temp}{description}       = $description;

          if (defined($accession_crapome_col))
          {
            $accession_annotation_hash{$accession}{crapome}              = $crapome;
          }
        }
      }
    }
    
    # store extra annotation, by col
    # if we store by col, we don't need to worry about header collisions
    foreach $col (@annotation_extra_col_array)
    {
        # skip CRAPome, since we handle it separately
        if (defined($accession_crapome_col) &&
            $col == $accession_crapome_col)
        {
            next;
        }
    
        $value = $array[$col];
        
        if (!defined($value) || !($value =~ /\S/) ||
            $value eq '---')
        {
            $value = '';
        }

        $accession_annotation_extra_hash{$accession}{$col} = $value;

        if ($map_with_additional_columns_flag)
        {
          @accession_array = split /\|/, $accession_rna;
          foreach $accession_temp (@accession_array)
          {
            if ($accession_temp =~ /[A-Za-z0-9]/ &&
                !($accession_temp =~ /^N\/*A$/i))
            {
                $accession_annotation_extra_hash{$accession_temp}{$col} =
                    $value;
            }
          }
          @accession_array = split /\|/, $accession_protein;
          foreach $accession_temp (@accession_array)
          {
            if ($accession_temp =~ /[A-Za-z0-9]/ &&
                !($accession_temp =~ /^N\/*A$/i))
            {
                $accession_annotation_extra_hash{$accession_temp}{$col} =
                    $value;
            }
          }
        }
    }
}
close ANNOTATION;



# skip down to first non-blank line
while(defined($line=<DATA>))
{
    $line =~ s/[\r\n]+//g;
    $line =~ s/\"//g;
    $flag = 1;

    if (!($line =~ /\S/)) { $flag = 0; }
#    if ($line =~ /^#/) { $flag = 0; }

    if ($flag) { last; }
    
    # print non-skipped lines unaltered
#    printf "%s\n", $line;
}


# data header line
$line =~ s/[\r\n]+//;
$line =~ s/\"//g;
@array = split /\t/, $line;
$num_accession_cols = 0;
@accession_col_array = ();
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;

    $field = $array[$i];

    if ($field =~ /\S/)
    {
        if (defined($header_remap_hash{$field}))
        {
            $field = $header_remap_hash{$field};
        }
        
        if ($field =~ /^Accession/i ||
            $field =~ /^Proteins[_ ]/i ||
            $field =~ /^Uniprot[_ ]/i ||
            $field =~ /^SwissProt[_ ]/i ||
            $field =~ /^RefSeq[_ ]/i ||
            $field =~ /^Ensemble[_ ]/i ||
            $field eq 'Uniprot' ||
            $field eq 'prey_ac' ||
            ($field =~ /^Protein/i && !($field =~ /^Protein\s*Group\s*ID/i)) ||
            $field =~ /^Leading\s*Protein/i ||
            $field =~ /Leading razor protein/i ||
            $field =~ /\bprotein ids*\b/i ||
            $field =~ /^Majority Protein IDs*$/i ||
            $field =~ /^Parent Protein$/i ||
            $field =~ /^Mapped Proteins*$/i)
        {
            $accession_col_array[$num_accession_cols++] = $i;
        }

        if ($field eq 'Protein' ||
            $field eq 'Proteins' ||
            $field =~ /^Leading\s*Protein/i ||
            $field =~ /^Majority\s*Protein/i)
        {
            $conrev_col_hash{$i} = 1;
        }
        
        if ($field =~ /^ModificationID$/i)
        {
            $modificationid_col = $i;
        }
        
        if ($field =~ /^Amino Acid$/i)
        {
            $amino_acid_col = $i;
        }

        $header_col_hash_data{$field} = $i;
    }
}
#$header_line_data = join "\t", @array;
$num_header_cols = @array;


$reverse_col = $header_col_hash_data{'Reverse'};
$contam_col  = $header_col_hash_data{'Contaminant'};

if (!defined($reverse_col))
{
    $reverse_col = $header_col_hash_data{'Reverse_flag'};
}
if (!defined($contam_col))
{
    $contam_col = $header_col_hash_data{'Species_contaminant_flag'};
}

# remove amino acid column from header, since we're going to move it
if (defined($amino_acid_col))
{
    @array_new = ();
    
    $j = 0;
    for ($i = 0; $i < @array; $i++)
    {
        if ($i != $amino_acid_col)
        {
            $array_new[$j++] = $array[$i];
        }
    }

    $header_line_data_new = join "\t", @array_new;
}
else
{
    $header_line_data_new = join "\t", @array;
}


printf "%s",   $header_line_data_new;
printf "\t%s", 'Has_Contaminant_Flag';
printf "\t%s", 'All_Reverse_Flag';
printf "\t%s", 'Cow_Contaminant_Flag';
printf "\t%s", 'Target_Species_Flag';

if ($proteogenomics_flag)
{
    printf "\t%s", 'Unique_To_Mutant';
}

printf "\t%s", 'Accession_Protein';


# insert CRAPome fields
if (defined($accession_crapome_col))
{
    printf "\tCRAPomeScore\tCRAPomeScoreMax";
}

if ($proteogenomics_flag)
{
    printf "\t%s", 'Accession_RNA';
}

# insert extra annotation fields
foreach $header (@annotation_extra_header_array)
{
    # skip CRAPomeScore, since we handle it separately
    if ($header eq 'CRAPomeScore')
    {
        next;
    }

    printf "\t%s", $header;
}


printf "\t%s", 'GeneID';
printf "\t%s", 'NumGenes';
printf "\t%s", 'WarningLevel';
printf "\t%s", 'Symbol';
if (defined($amino_acid_col))
{
    printf "\t%s", $array[$amino_acid_col];
}
if (defined($modificationid_col))
{
    printf "\t%s", 'Positions Reannotated';
}
printf "\t%s", 'Alias';
printf "\t%s", 'Type';
printf "\t%s", 'Location';
printf "\t%s", 'Description';
printf "\n";


if ($num_accession_cols == 0)
{
#    printf STDERR "ABORT -- No accession columns present\n";
#    exit(1);
    printf STDERR "WARNING -- No accession columns present\n";
}


# read in data
while(defined($line=<DATA>))
{
    $line =~ s/[\r\n]+//g;
    $line =~ s/\"//g;

    @array = split /\t/, $line;

    # clean up line
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        
        if (!($array[$i] =~ /\S/))
        {
#            $array[$i] = '---';
            $array[$i] = '';
        }
    }
    # pad out missing fields at end of line
    for ($i = @array; $i < $num_header_cols; $i++)
    {
#        $array[$i] = '---';
        $array[$i] = '';
    }


    %contam_accession_hash   = ();
    %reversed_accession_hash = ();
    %local_accession_hash    = ();
    %potential_conrev_accession_hash = ();
    
    $reverse_orig = 0;
    $contam_orig  = 0;
    
    $mutant_flag  = 0;
    $ensp_flag    = 0;
    
    if (defined($reverse_col))
    {
        if ($array[$reverse_col] eq '+' ||
            $array[$reverse_col] eq '1')
        {
            $reverse_orig = 1;
        }
    }
    if (defined($contam_col))
    {
        if ($array[$contam_col] eq '+' ||
            $array[$contam_col] eq '1')
        {
            $contam_orig = 1;
        }
    }

    $cow_pig_shroom_flag = 0;
    foreach $accession_col (@accession_col_array)
    {
        $accession_str = $array[$accession_col];

        # check SwissProt for contaminants before we strip SwissProt out
        if ($accession_str =~ /_(BOVIN|PIG|GRIFR|SHEEP|LYSEN)\b/i)
        {
            $cow_pig_shroom_flag = 1;
        }

        # handle database identifier junk
        # otherwise, REV_ CON_ mapping will break on sp|stuff, etc.
#        $accession_str =~ s/sp\|[^\|\;\,]+\|([^\|\;\,]+)/$1/g;
        $accession_str =~ s/sp\|([^\|\;\,]+)\|([^\|\;\,]+)/$1/g;
        $accession_str =~ s/ENSEMBL:([A-Za-z0-9-\.]+)/$1/ig;
        $accession_str =~ s/REFSEQ:([A-Za-z0-9\.]+)/$1/ig;
        $accession_str =~ s/ref\|([A-Za-z0-9\.]+)/$1/ig;

        # strip SwissProt alias from uniprot fasta identifiers
        #  only the uniprot identifier part is flagged as con/rev elsewhere
        #  including the SwissProt identifier results in con/rev inclusion!
        $accession_str =~ s/>[A-Za-z0-9]+\|([A-Za-z0-9]+)\|[A-Za-z0-9]+/$1/g;
    
        $accession_str =~ s/\|/\;/g;
#        $accession_str =~ s/\:/\;/g;
        $accession_str =~ s/,/\;/g;
        $accession_str =~ s/\s+/\;/g;
        $accession_str =~ s/\;+/\;/g;
        $accession_str =~ s/^\;//;
        $accession_str =~ s/\;$//;
    
        @accession_array = split /\;/, $accession_str;

        # store reversed accessions for later removal
        for ($i = 0; $i < @accession_array; $i++)
        {
            $accession = $accession_array[$i];
            
            # some non-standard accessions use periods as delimiters
            @temp_array = $accession =~ /\./g;
            if (@temp_array == 1)
            {
                $accession =~ s/\.\S+//g;
            }
            
            # remove FASTA junk
            # strangely, \b> does NOT match at the beginnig of the string !!
            #  so, I have to use negative lookbehind instead... perl bug ??
#            $accession =~ s/\b>//g;
            $accession =~ s/(?<![A-Za-z0-9])>//g;

            # strip -# from end of accession
#            $accession =~ s/(\d+)-\d+$/$1/g;

            # uppercase accession for later lookup
#            $accession =~ tr/a-z/A-Z/;
        
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

            # remove REV__, rev_
            $reverse_flag = 0;
            if ($accession =~ /(?<![A-Za-z0-9])REV_+/i)
            {
                $accession =~ s/(?<![A-Za-z0-9])REV_+//ig;
                $reverse_flag = 1;
            }
            
            # skip empty accessions
            if (!($accession =~ /[A-Za-z0-9]/))
            {
                next;
            }


            # detect and deal with proteogenomics ENSTs
            # use \bchr, not just chr, to skip rev_chr
            if ($reverse_flag == 0 &&
                $accession =~ /\bchr[A-Za-z0-9]+_[0-9]+_[^ ]*(ENS[A-Z]{0,3}[PTG][0-9]+)/)
            {
                $accession   = $1;
                $mutant_flag = 1;
            }
            # check for \bENS, not just ENS,
            # since FragPipe can have rev_ENSP reverse sequences
            #
            # \b counts _ as a word, not a word boundary,
            # so it does what we want here
            #
            elsif ($reverse_flag == 0 &&
                   $accession =~ /\bENS[A-Z]{0,3}P[0-9]+/)
            {
                $ensp_flag   = 1;
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
                    $local_accession_hash{$accession2} = 1;
                    
                    if (defined($conrev_col_hash{$accession_col}))
                    {
                        $potential_conrev_accession_hash{$accession2} = 1;
                    }
                }
            }
            else
            {
                if ($contam_flag)
                {
                    $contam_accession_hash{$accession} = 1;
                    $potential_conrev_accession_hash{$accession} = 1;
                }
                if ($reverse_flag)
                {
                    $reversed_accession_hash{$accession} = 1;
                    $potential_conrev_accession_hash{$accession} = 1;
                }

                $local_accession_hash{$accession} = 1;

                if (defined($conrev_col_hash{$accession_col}))
                {
                    $potential_conrev_accession_hash{$accession} = 1;
                }
            }
        }
    }

    # store new accessions
    @new_accession_array = ();
    $num_new_accessions = 0;
    foreach $accession (sort compare_accession keys %local_accession_hash)
    {
        # skip reversed accessions
        if (defined($reversed_accession_hash{$accession}))
        {
#            printf STDERR "REVERSE\t%s\n", $accession;
             next;
        }

        $new_accession_array[$num_new_accessions++] = $accession;
    }
    

    # sanity check, reverse sequences aren't actually contaminants
    @temp_array = keys %contam_accession_hash;
    foreach $accession (@temp_array)
    {
        if (defined($reversed_accession_hash{$accession}))
        {
            delete $contam_accession_hash{$accession};
        }
    }
    

    # count CON__ and REV__
    $has_contam  = $contam_orig;
    $all_reverse = $reverse_orig;
    $num_total   = 0;
    $num_reverse = 0;
    $num_contam  = 0;
    foreach $accession (keys %local_accession_hash)
    {
        $num_total++;
    
        # adding cow_pig_shroom flag will overcount, but it doesn't matter
        if (defined($contam_accession_hash{$accession}) ||
            $cow_pig_shroom_flag)
        {
            $num_contam++;
        }
        if (defined($reversed_accession_hash{$accession}))
        {
            $num_reverse++;
        }
    }
    if ($num_contam)
    {
        $has_contam += 2;
    }
    # we may need to remove the +2 later, during printing
    if ($num_reverse)
    {
        $all_reverse += 2;
    }
    
    $new_accession_str = join "\|", @new_accession_array;
    
    # scan for cow contaminates
    $cow_contam = 0;
    if ($num_contam)
    {
        for ($i = 0; $i < @array; $i++)
        {
            if ($array[$i] =~ /bos taurus/i ||
                $array[$i] =~ /_BOVIN\b/i)
            {
                $cow_contam = 1;
            }
        }
    }

    
#    printf STDERR "FOOBAR\t%s\n", $new_accession_str;

    $line_new = join "\t", @array;

    print_probeid_annotation($line_new, $new_accession_str);
}
