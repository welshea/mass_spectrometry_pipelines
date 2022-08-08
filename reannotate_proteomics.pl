#!/usr/bin/perl -w

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
        # Perl treats infinities as numbers, Excel does not
        if ($_[0] =~ /^[-+]*inf/)
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
    my $accession_a = $a;
    my $accession_b = $b;
    my %hash_a = %{$accession_annotation_hash{$accession_a}};
    my %hash_b = %{$accession_annotation_hash{$accession_b}};
    my $symbol_a = $hash_a{symbol};
    my $symbol_b = $hash_b{symbol};
    my $description_a = $hash_a{description};
    my $description_b = $hash_b{description};
    my $type_a   = $hash_a{type};
    my $type_b   = $hash_b{type};
    my $value_a;
    my $value_b;
    my $temp_a;
    my $temp_b;
    my $sub_description;
    my $temp_length;
    my $debug = 0;
    
    # best refseqs
    if (  $accession_a =~ /^NM_/  && !($accession_b =~ /^NM_/)) {return -1;}
    if (!($accession_a =~ /^NM_/) &&   $accession_b =~ /^NM_/)  {return 1;}
    if (  $accession_a =~ /^NR_/  && !($accession_b =~ /^NR_/)) {return -1;}
    if (!($accession_a =~ /^NR_/) &&   $accession_b =~ /^NR_/)  {return 1;}
    if (  $accession_a =~ /^NP_/  && !($accession_b =~ /^NP_/)) {return -1;}
    if (!($accession_a =~ /^NP_/) &&   $accession_b =~ /^NP_/)  {return 1;}

    # symbol sketchiness
    $value_a = score_symbol_sketchiness($symbol_a);
    $value_b = score_symbol_sketchiness($symbol_b);
    if ($value_a > $value_b) {return -1;}
    if ($value_a < $value_b) {return 1;}

    # worse refseqs
    if (  $accession_a =~ /^XM_/  && !($accession_b =~ /^XM_/)) {return -1;}
    if (!($accession_a =~ /^XM_/) &&   $accession_b =~ /^XM_/)  {return 1;}
    if (  $accession_a =~ /^XR_/  && !($accession_b =~ /^XR_/)) {return -1;}
    if (!($accession_a =~ /^XR_/) &&   $accession_b =~ /^XR_/)  {return 1;}
    if (  $accession_a =~ /^XP_/  && !($accession_b =~ /^XP_/)) {return -1;}
    if (!($accession_a =~ /^XP_/) &&   $accession_b =~ /^XP_/)  {return 1;}
    if (  $accession_a =~ /^NG_/  && !($accession_b =~ /^NG_/)) {return -1;}
    if (!($accession_a =~ /^NG_/) &&   $accession_b =~ /^NG_/)  {return 1;}
    
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
    if ($value_a < $value_b) {return 1;}

    # absence of pseudogene
    $value_a = 0;
    if ($hash_a{type} =~ /pseudo/i ||
        $description_a =~ /pseudogene/i)
    {
        $value_a = 1;
    }
    $value_b = 0;
    if ($hash_b{type} =~ /pseudo/i ||
        $description_b =~ /pseudogene/i)
    {
        $value_b = 1;
    }
    if ($value_a < $value_b) {return -1;}
    if ($value_a > $value_b) {return 1;}

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
    if ($value_a > $value_b) {return 1;}

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
    if ($value_a > $value_b) {return 1;}


    # alphabetical
    if ($symbol_a lt $symbol_b) {return -1;}
    if ($symbol_a gt $symbol_b) {return 1;}
    
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
    if ($value_a < $value_b) {return 1;}
    
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
    if ($value_a < $value_b) {return 1;}

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
    if ($value_a < $value_b) {return 1;}

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
    if ($value_a < $value_b) {return 1;}
    
    # description, alphabetical
    if (lc($description_a) lt lc($description_b)) {return -1;}
    if (lc($description_a) gt lc($description_b)) {return 1;}
    
    # shorter symbol is better
    $value_a = length $symbol_a;
    $value_b = length $symbol_b;
    if ($value_a < $value_b) {return -1;}
    if ($value_a > $value_b) {return 1;}
    
    # presence of alias
    $value_a = 0;
    if ($hash_a{alias} =~ /[A-Za-z0-9]/)
    {
        $value_a = 1;
    }
    $value_b = 0;
    if ($hash_b{alias} =~ /[A-Za-z0-9]/)
    {
        $value_b = 1;
    }
    if ($value_a > $value_b) {return -1;}
    if ($value_a < $value_b) {return 1;}

    # sketchiness of alias
    $value_a = score_symbol_sketchiness($value_a);
    $value_b = score_symbol_sketchiness($value_b);
    if ($value_a > $value_b) {return -1;}
    if ($value_a < $value_b) {return 1;}
    
    # alias, alphabetical
    $value_a = $hash_a{alias};
    $value_b = $hash_b{alias};
    if (lc($value_a) lt lc($value_a)) {return -1;}
    if (lc($value_a) gt lc($value_b)) {return 1;}
    
    # uniprot first
    if ($accession_a =~ /_/ && !($accession_b =~ /_/)) {return -1;}
    if ($accession_b =~ /_/ && !($accession_a =~ /_/)) {return  1;}

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

    %seen_accession_rna_hash = ();
    %seen_accession_protein_hash = ();
    %seen_gene_id_hash = ();
    %seen_symbol_hash = ();
    %seen_alias_hash = ();
    %seen_type_hash = ();
    %seen_location_hash = ();
    %seen_description_hash = ();
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
        
#        if (!defined($accession_rna))     { $accession_rna     = ''; }
#        if (!defined($accession_protein)) { $accession_protein = ''; }
#        if (!defined($gene_id))           { $gene_id           = ''; }
#        if (!defined($symbol))            { $symbol            = ''; }
#        if (!defined($alias))             { $alias             = ''; }
#        if (!defined($type))              { $type              = ''; }
#        if (!defined($location))          { $location          = ''; }
#        if (!defined($description))       { $description       = ''; }

        # remove NULL fields
        $accession_rna     =~ s/^null$//i;
        $accession_protein =~ s/^null$//i;
        $gene_id           =~ s/^null$//i;
        $symbol            =~ s/^null$//i;
        $alias             =~ s/^null$//i;
        $type              =~ s/^null$//i;
        $location          =~ s/^null$//i;
        $description       =~ s/^null$//i;

        $accession_rna     = bless_delimiter_bar($accession_rna);
        $accession_protein = bless_delimiter_bar($accession_protein);
        $gene_id           = bless_delimiter_bar($gene_id);
        $symbol            = bless_delimiter_bar($symbol);
        $alias             = bless_delimiter_bar($alias);
        $type              = bless_delimiter_bar($type);
        
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

    $accession_rna_str     = bless_delimiter_space($accession_rna_str);
    $accession_protein_str = bless_delimiter_space($accession_protein_str);
    $gene_id_str           = bless_delimiter_space($gene_id_str);
    $symbol_str            = bless_delimiter_space($symbol_str);
    $alias_str             = bless_delimiter_space($alias_str);
    $type_str              = bless_delimiter_bar($type_str);
    $type_str              =~ s/\;/\|/g;


    # remove symbols without geneids if we already have some
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
            if ($field =~ /^[^-]+\-[^-]+$/)
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
        
        # add fusion genes to counts
        foreach $field (sort keys %fusion_hash)
        {
            $count_flag = 1;

            if ($field =~ /^([^-]+)\-([^-]+)$/)
            {
                $symbol1 = $1;
                $symbol2 = $2;
                
                # check to see if symbol2 looks like a real symbol or not

                # example: FOOBAR-AS
                $length2 = length $symbol2;
                if ($length2 <= 2)
                {
                    $count_flag = 1;
                }
                # example: FOOBAR-AS1
                elsif ($symbol2 =~ /^[A-Za-z][A-Za-z][0-9]+$/)
                {
                    $count_flag = 1;
                }
                # example: FOOBAR-1
                elsif ($symbol2 =~ /^[0-9]+$/)
                {
                    $count_flag = 1;
                }
                
                # do not count if already covered by a non_fusion symbol
                elsif (defined($non_fusion_hash{$symbol1}) ||
                    defined($non_fusion_hash{$symbol2}))
                {
                    $count_flag = 0;
                }
            }
            
            if ($count_flag)
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
    if (@accession_pruned_array)
    {
        $accession_pruned_str = join ' ', @accession_pruned_array;
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
    $line_new .= "\t$accession_pruned_str";
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

$strip_flag = 0;

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



while(defined($line=<ANNOTATION>))
{
    $line =~ s/[\r\n]+//;
    @array = split /\t/, $line;
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
            $field =~ /^Proteins_/i ||
            $field eq 'Proteins' ||
            $field =~ /^Uniprot_/i ||
            $field eq 'Uniprot' ||
            $field eq 'prey_ac' ||
            ($field =~ /^Protein/i && !($field =~ /^Protein\s*Group\s*ID/i)) ||
            $field =~ /^Leading\s*Protein/i ||
            $field =~ /Leading razor protein/i ||
            $field =~ /\bprotein ids\b/i ||
            $field =~ /^Protein IDs$/i ||
            $field =~ /^Majority Protein IDs$/i ||
            $field =~ /^Parent Protein$/i)
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
printf "\t%s", 'Accession_Protein';
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

    foreach $accession_col (@accession_col_array)
    {
        $accession_str = $array[$accession_col];

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

            # remove REV__
            $reverse_flag = 0;
            if ($accession =~ /(?<![A-Za-z0-9])REV__/)
            {
                $accession =~ s/(?<![A-Za-z0-9])REV__//g;
                $reverse_flag = 1;
            }
            
            # skip empty accessions
            if (!($accession =~ /[A-Za-z0-9]/))
            {
                next;
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
    foreach $accession (sort keys %local_accession_hash)
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
    
        if (defined($contam_accession_hash{$accession}))
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
            if ($array[$i] =~ /bos taurus/i)
            {
                $cow_contam = 1;
            }
        }
    }

    
#    printf STDERR "FOOBAR\t%s\n", $new_accession_str;

    $line_new = join "\t", @array;

    print_probeid_annotation($line_new, $new_accession_str);
}
