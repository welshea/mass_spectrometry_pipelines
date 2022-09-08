#!/usr/bin/perl -w


# 2022-09-08:  don't split m/z on forward slash
# 2022-07-15:  conform by removing likely abbreviations from ends
# 2022-07-15:  update is_heavy_labeled() function
# 2022-07-08:  added some formula-based sanity checks
# 2022-07-08:  rework and renumber the match groups
# 2022-06-30:  negligible speedup of certain match types
# 2022-06-30:  sort rows on name, to get more predictable annotation order
# 2022-06-30:  improved name conforming; DL/+-, sulph->sulf
# 2022-06-29:  reorder and alter various match types
# 2022-06-27:  swap 5A/B and 6A/B categories (fuzzy pos/neg is last now)
# 2022-06-27:  do not warn on mismatched bogus m/z values
# 2022-06-25:  preliminary synonym matching support
# 2022-06-25:  limit search to appropriate pos/neg ionization
# 2022-06-25:  conform various "acids", which should be "acid"
# 2022-06-24:  handle non-zero formal charges
# 2022-06-22:  don't include placeholder bad rowid in all-row scan
# 2022-06-22:  big speedups
# 2022-06-22:  begin adding HMDB adduct support
# 2022-06-22:  conform unicode, sulfid/sulfide/sulphid/sulphide
# 2022-06-08:  enable 5A/5B fuzzy matching on opposite pos/neg
# 2022-06-08:  adjust alignment method and scoring cutoff
# 2022-04-04:  conform [1-7] to [abgdezh], "Acid,*ic" to "*ic acid"
# 2022-02-15:  include PPM as 3rd field of m/z WARNING messages
# 2021-11-30:  add minimal lipidomics support
# 2021-10-27:  rename gmiddle to elocal
# 2021-10-21:  rename glocal to gmiddle
# 2021-08-17:  begin adding auto heavy label matching support
# 2021-08-17:  add --ppm flag to set m/z PPM tolerance


# some problematic ambiguous names, not much way around these:
#
#    Sorbate
#      HMDB0029581  (2E,4E)-2,4-Hexadienoic acid
#      HMDB0253115  Hexadienic acid
#      HMDB0256823  Propenylacrylic Acid
#
#    Methyl glutarate
#      HMDB0000752  Methylglutaric acid
#      HMDB0000858  Monomethyl glutaric acid

# problems with Lactose
#
# HMDB0041627  beta-Lactose    C12H22O11  milk lactose
# HMDB0035312  Hebevinoside I  C44H72O13  toxin from Hebeloma mushroom
# HMDB0000186  Alpha-Lactose   C12H22O11  milk lactose
#
# Unfortunately, HMDB0035312 contains most or all of the usual synoyms
# for Lactose, and links up the the correct KEGG and other external accessions.
# HMDB0041627, which is actually Lactose, isn't linked up to the correct
# external identifiers (linked to a KEGG identifier that isn't in any
# pathways).  HMDB0000186 links out correctly, but is missing the D-lactose
# synonym, so we don't map it to D-Lactose.  Maybe I can get HMDB to fix
# their entry for Hebevinoside I to de-lactoseify it.
#
# So, until HMDB fixes HMDB0035312, we need some way to map D-lactose
# to Alpha-Lactose (HMDB0000186), since that is the most correct.
# Implementing a formula sanity check filter should force us to map to
# both alpha- and beta- lactose instead.  Unfortunately, B-D-lactose
# aligns better to D-lactose than A-lactose does, so beta-lactose is still
# the chosen HMDB accession instead of alpha-lactose :-(
#
# I'm going to add in KEGG pathway identifiers, to select those over
# KEGG entries that aren't in pathways, but that still won't help Lactose,
# due to the synonym alignment scoring issue.


# set lib search path to directory the script is run from
use File::Basename;
use lib dirname (__FILE__);

use Scalar::Util qw(looks_like_number);
use POSIX;
use align_text;    # text string alignment module

$mz_tol_ppm   = 10;    # 10 ppm
$rt_tol       = 1.0;   # minutes
$align_method = 'overlap';

$bad_row_id      = 9E99;
$lipidomics_flag = 0;

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


# handle heavy elements as well
sub cmp_elements
{
    my $ele_a;
    my $ele_b;
    my $heavy_a;
    my $heavy_b;

    $heavy_a = '';
    $heavy_b = '';

    if ($a =~ /^\[([0-9]+)\]/)
    {
        $heavy_a = $1;
    }
    if ($b =~ /^\[([0-9]+)\]/)
    {
        $heavy_b = $1;
    }
    
    $a     =~ /([A-Za-z]+)/;
    $ele_a = $1;
    
    $b     =~ /([A-Za-z]+)/;
    $ele_b = $1;

    # first by element
    if ($ele_a ne $ele_b)
    {
        return $ele_a cmp $ele_b;
    }

    # then put heavy labeled atoms first
    if ($heavy_a ne '' && $heavy_b eq '') { return -1; }
    if ($heavy_b ne '' && $heavy_a eq '') { return  1; }

    # then sort by number of heavy
    if ($heavy_a != $heavy_b)
    {
        return $heavy_a <=> $heavy_b;
    }
   
    return $a cmp $b;
}


sub cmp_rows
{
    my $value_a;
    my $value_b;
    
    $value_a = $annotation_hash{$a}{name};
    $value_b = $annotation_hash{$b}{name};
    if ($value_a lt $value_b) { return -1; }
    if ($value_a gt $value_b) { return  1; }

    $value_a = $annotation_hash{$a}{mz};
    $value_b = $annotation_hash{$b}{mz};
    if ($value_a < $value_b) { return -1; }
    if ($value_a > $value_b) { return  1; }
    
    return ($a <=> $b);
}


sub cmp_conformed_rows
{
    my $value_a;
    my $value_b;
    
    # has KEGG
    $value_a = $annotation_hash{$a}{kegg};
    $value_b = $annotation_hash{$b}{kegg};
    if (defined($value_a) && $value_a ne '') { $value_a = 1; }
    else                                     { $value_a = 0; }
    if (defined($value_b) && $value_b ne '') { $value_b = 1; }
    else                                     { $value_b = 0; }
    if ($value_a && $annotation_hash{$a}{kegg_map} ne '') { $value_a = 2; }
    if ($value_b && $annotation_hash{$b}{kegg_map} ne '') { $value_b = 2; }
    if ($value_a > $value_b) { return -1; }
    if ($value_a < $value_b) { return  1; }
    
    # by number of interesting columns
    $value_a = $annotation_hash{$a}{col_count};
    $value_b = $annotation_hash{$b}{col_count};
    if ($value_a > $value_b) { return -1; }
    if ($value_a < $value_b) { return  1; }
    
    # by score
    $value_a = $temp_row_score_hash{$a};
    $value_b = $temp_row_score_hash{$b};
    if ($value_a > $value_b) { return -1; }
    if ($value_a < $value_b) { return  1; }
    

    # then by the usual alphabetical and m/z sort order

    $value_a = $annotation_hash{$a}{name};
    $value_b = $annotation_hash{$b}{name};
    if ($value_a lt $value_b) { return -1; }
    if ($value_a gt $value_b) { return  1; }

    $value_a = $annotation_hash{$a}{mz};
    $value_b = $annotation_hash{$b}{mz};
    if ($value_a < $value_b) { return -1; }
    if ($value_a > $value_b) { return  1; }
    
    return ($a <=> $b);
}


sub cmp_fuzzy_rows
{
    my $value_a;
    my $value_b;

    # by score
    $value_a = $temp_row_score_hash{$a};
    $value_b = $temp_row_score_hash{$b};
    if ($value_a > $value_b) { return -1; }
    if ($value_a < $value_b) { return  1; }
    
    # has KEGG
    $value_a = $annotation_hash{$a}{kegg};
    $value_b = $annotation_hash{$b}{kegg};
    if (defined($value_a) && $value_a ne '') { $value_a = 1; }
    else                                     { $value_a = 0; }
    if (defined($value_b) && $value_b ne '') { $value_b = 1; }
    else                                     { $value_b = 0; }
    if ($value_a && $annotation_hash{$a}{kegg_map} ne '') { $value_a = 2; }
    if ($value_b && $annotation_hash{$b}{kegg_map} ne '') { $value_b = 2; }
    if ($value_a > $value_b) { return -1; }
    if ($value_a < $value_b) { return  1; }
    
    # by number of interesting columns
    $value_a = $annotation_hash{$a}{col_count};
    $value_b = $annotation_hash{$b}{col_count};
    if ($value_a > $value_b) { return -1; }
    if ($value_a < $value_b) { return  1; }


    # then by the usual alphabetical and m/z sort order

    $value_a = $annotation_hash{$a}{name};
    $value_b = $annotation_hash{$b}{name};
    if ($value_a lt $value_b) { return -1; }
    if ($value_a gt $value_b) { return  1; }

    $value_a = $annotation_hash{$a}{mz};
    $value_b = $annotation_hash{$b}{mz};
    if ($value_a < $value_b) { return -1; }
    if ($value_a > $value_b) { return  1; }
    
    return ($a <=> $b);
}


# The Hill system specifies C#H#D#, not C#D#H#
#   list all elements in alphabetical order,
#   unless it contains a carbon, then list carbon then hydrogen first
#
# I check for 3-letter elements (all the Uuu's have real symbols by now),
# and elements listed multiple times, and exit early with the original
# formula if such errors are detected.
#
# I don't current check to see if the given 1- or 2- letter elements are
# valid known elements or not.  I could, but that would require a good bit
# more work than I have time for at the moment.  I'm not *quite* that paranoid
# about the formulas just yet, although part of me still worries about it...
#
sub conform_formula
{
    my $formula_orig = $_[0];
    my $formula;
    my $formula_new = '';
    my @match_array;
    my @element_array;
    my $match;
    my $heavy;
    my $element;
    my $heavy_plus_element;
    my $count;
    my %count_hash = ();
    my $has_carbon_flag = 0;
    
    $formula         = $formula_orig;
    $has_carbon_flag = 0;
    
    # element with number
    @match_array = $formula =~ m/(?:\[[0-9]+\])*[A-Z][a-z]*(?:[0-9]+)*/g;
    foreach $match (@match_array)
    {
        $match =~ /((?:\[[0-9]+\])*)([A-Za-z]+)([0-9]+)*/;

        $heavy   = $1;
        $element = $2;
        $count   = $3;

        if (!defined($heavy))
        {
            $heavy = '';
        }
        if (!defined($count))
        {
            $count = 1;
        }

        if ($element eq 'C')
        {
            $has_carbon_flag = 1;
        }

        $heavy_plus_element = $heavy . $element;

        if (length $element > 2 || defined($count_hash{$heavy_plus_element}))
        {
            printf STDERR "WARNING -- error in formula %s\n", $formula_orig;
            return $formula_orig;
        }
        
        $count_hash{$heavy_plus_element} = $count;
    }
    
    @element_array = sort cmp_elements keys %count_hash;
    
    # order all carbons first, followed by hydrogens
    if ($has_carbon_flag)
    {
        # print all carbons first
        foreach $heavy_plus_element (@element_array)
        {
            $heavy_plus_element =~ /([A-Za-z]+)/;
            $element = $1;
            
            if ($element eq 'C')
            {
                $count = $count_hash{$heavy_plus_element};

                # replace 1 count with blank
                if ($count == 1)
                {
                    $count = '';
                }

                $formula_new .= $heavy_plus_element . $count;
            }
        }
        
        # then all hydrogens that aren't D's
        foreach $heavy_plus_element (@element_array)
        {
            $heavy_plus_element =~ /([A-Za-z]+)/;
            $element = $1;
            
            if ($element eq 'H')
            {
                $count = $count_hash{$heavy_plus_element};

                # replace 1 count with blank
                if ($count == 1)
                {
                    $count = '';
                }

                $formula_new .= $heavy_plus_element . $count;
            }
        }
    }

    # order the remaining elements alphabetically, including deuterium
    foreach $heavy_plus_element (@element_array)
    {
        # skip C's and H's we've already placed first
        if ($has_carbon_flag)
        {
            $heavy_plus_element =~ /([A-Za-z]+)/;
            $element = $1;
            
            if ($element eq 'C' || $element eq 'H')
            {
                next;
            }
        }
    
        $count = $count_hash{$heavy_plus_element};
        
        # replace 1 count with blank
        if ($count == 1)
        {
            $count = '';
        }
        
        $formula_new .= $heavy_plus_element . $count;
    }
    
    return $formula_new;
}


# oh no, semicolon can appear within heavy labels
# ex: L-LYSINE (13C6, 99%; 15N2, 99%)
sub bless_delimiter_bar_metabolomics
{
    my $text = $_[0];
    my @temp_array;
    my $n;
    my $i;

    # convert semicolons only if they are not within ()
    # incrementing $i within array access for minor speed increase
    #   requires initializing things to -2
    #
    # regular expression from MichaelRushton:
    #   https://stackoverflow.com/questions/133601/can-regular-expressions-be-used-to-match-nested-patterns
    @temp_array = split /(\((?>[^()]+|(?1))*\))/, $text;
    $n = @temp_array - 2;
    for ($i = -2; $i < $n;)
    {
        $i += 2;
        $temp_array[$i] =~ tr/\;/\|/;

        # / can also be a delimiter in some much older metabolomics data
        # protect m/z
        $temp_array[$i] =~ s/\bm\/z\b/M_OvEr_Z/g;
        $temp_array[$i] =~ tr/\//\|/;
        $temp_array[$i] =~ s/M_OvEr_Z/m\/z/g;
    }
    $text = join '', @temp_array;

    # clean up delimiters
    $text =~ s/\|+/\|/g;
    $text =~ s/^\|//;
    $text =~ s/\|$//;
    
    
    # clean up spaces
    $text =~ s/^\s+//;
    $text =~ s/\s+$//;
    $text =~ s/\s+\|/\|/g;
    $text =~ s/\|\s+/\|/g;
    
    return $text;
}

sub reformat_sci
{
    my $field = $_[0];
    my $temp;

    # n.nnEn, nn.nEn, etc.
    #
    # If it looks like scientific notation, Excel will automatically
    #  format it as scientific notation.  However, if the magnitude is
    #  >= 1E7, it will also automatically display it set to only 2 digits
    #  to the right of the decimal (it is still fine internally).  If the
    #  file is re-exported to text, truncation to 3 significant digits
    #  will occur!!!
    #
    # Reformat the number to (mostly) avoid this behavior.
    #
    # Unfortunately, the >= 11 significant digits behavior still
    #  triggers, so it still truncates to 10 digits when re-exporting
    #  General format.  10 digits is still better than 3....
    #
    # The re-export truncation behavior can only be more fully avoided by
    #  manually setting the format to Numeric and specifying a large
    #  number of digits after the decimal place for numbers with
    #  fractions, or 0 digits after the decimal for whole numbers.
    #
    # Ugh.
    #
    # There is no fully fixing this brain-damagedness automatically,
    #  I can only decrease the automatic truncation of significant
    #  digits from 3 to 10 digits :(  Any precision beyond 10 digits
    #  *WILL* be lost on re-export if the format is set to General.
    #
    # NOTE -- We truncate to 16 significant digits by going through
    #         a standard IEEE double precision intermediate.
    #         However, Excel imports numbers as double precision
    #         anyways, so we aren't losing any precision that Excel
    #         wouldn't already be discarding.
    #
    if (is_number($field))
    {
          # strip commas
          $temp = $field;
          $temp =~ s/\,//g;
          
          if (abs($temp) >= 1 &&
              $temp =~ /^([-+]?[0-9]*\.*[0-9]*)[Ee]([-+]?[0-9]+)$/)
          {
#              $number   = $1 + 0;
#              $exponent = $2 + 0;
              
              $temp /= 1;

              # replace original with new scientific notation format
              $field = $temp;
          }
    }
    
    return $field;
}


sub is_heavy_labeled
{
    my $string = $_[0];
    
    if ($string =~ /\([^()]*\b13C[0-9]*\b[^()]*\)/) { return 1; }
    if ($string =~ /\([^()]*\b14N[0-9]*\b[^()]*\)/) { return 1; }
    if ($string =~ /\([^()]*\bD[0-9]*\b[^()]*\)/)   { return 1; }
    if ($string =~ /\([^()]*\bBOC\b[^()]*\)/)       { return 1; }

    if ($string =~ /\b13C[0-9]+\b/) { return 1; }
    if ($string =~ /\b14N[0-9]+\b/) { return 1; }
    if ($string =~ /\bD[0-9]+\b/)   { return 1; }
    if ($string =~ /\bBOC\b/)       { return 1; }
    
    return 0;
}


# rules notes:
#
# see https://hmdb.ca/metabolites/HMDB0001906 for Greek examples
#
# replace Greek at word boundaries, s/\bGreek\b/English/i
#   alpha           : a
#   beta            : b
#   gamma           : g
#   delta           : d
#   epsilon         : e
#   zeta            : z
#   eta             : h
#   lambda          : l

# L-aminoacid     : aminoacid        ex. L-Alanine     --> Alanine
# "anoic"         : "yric"           ex. Butanoic acid --> Butyric acid
# "anoate\b"      : "yrate"          ex. Butanaoate    --> Butyrate
# "ic acid\b"     : "ate"            ex. glutamic acid --> glutamate
# "monosomething" : "something"      ex. monophosphate --> phosphate

# strip ()
# strip spaces
# maybe strip -

# lookup table for converting numbers into ordered greek letters
$number_letter_hash{1} = 'a';
$number_letter_hash{2} = 'b';
$number_letter_hash{3} = 'g';
$number_letter_hash{4} = 'd';
$number_letter_hash{5} = 'e';
$number_letter_hash{6} = 'z';
$number_letter_hash{7} = 'h';

# common Greek and punctuation unicode seen in metabolite names
# convert them to their nearest ASCII equivalent
#
# Greek
$unicode_to_ascii_hash{"\x{0391}"} = 'A';
$unicode_to_ascii_hash{"\x{0392}"} = 'B';
$unicode_to_ascii_hash{"\x{0393}"} = 'G';
$unicode_to_ascii_hash{"\x{0394}"} = 'D';
$unicode_to_ascii_hash{"\x{0395}"} = 'E';
$unicode_to_ascii_hash{"\x{0396}"} = 'Z';
$unicode_to_ascii_hash{"\x{0397}"} = 'H';
$unicode_to_ascii_hash{"\x{039B}"} = 'L';
$unicode_to_ascii_hash{"\x{03A6}"} = 'Phi';
$unicode_to_ascii_hash{"\x{03A8}"} = 'Psi';
$unicode_to_ascii_hash{"\x{03A9}"} = 'O';
$unicode_to_ascii_hash{"\x{03B1}"} = 'a';
$unicode_to_ascii_hash{"\x{03B2}"} = 'b';
$unicode_to_ascii_hash{"\x{03B3}"} = 'g';
$unicode_to_ascii_hash{"\x{03B4}"} = 'd';
$unicode_to_ascii_hash{"\x{03B5}"} = 'e';
$unicode_to_ascii_hash{"\x{03B6}"} = 'z';
$unicode_to_ascii_hash{"\x{03B7}"} = 'h';
$unicode_to_ascii_hash{"\x{03BB}"} = 'l';
$unicode_to_ascii_hash{"\x{03C6}"} = 'phi';
$unicode_to_ascii_hash{"\x{03C8}"} = 'psi';
$unicode_to_ascii_hash{"\x{03C9}"} = 'o';


# Latin-1 supplement letters block
if (0)
{
$unicode_to_ascii_hash{"\x{00C0}"} = 'A';
$unicode_to_ascii_hash{"\x{00C1}"} = 'A';
$unicode_to_ascii_hash{"\x{00C2}"} = 'A';
$unicode_to_ascii_hash{"\x{00C3}"} = 'A';
$unicode_to_ascii_hash{"\x{00C4}"} = 'A';
$unicode_to_ascii_hash{"\x{00C5}"} = 'A';
$unicode_to_ascii_hash{"\x{00C6}"} = 'AE';
$unicode_to_ascii_hash{"\x{00C7}"} = 'C';
$unicode_to_ascii_hash{"\x{00C8}"} = 'E';
$unicode_to_ascii_hash{"\x{00C9}"} = 'E';
$unicode_to_ascii_hash{"\x{00CA}"} = 'E';
$unicode_to_ascii_hash{"\x{00CB}"} = 'E';
$unicode_to_ascii_hash{"\x{00CC}"} = 'I';
$unicode_to_ascii_hash{"\x{00CD}"} = 'I';
$unicode_to_ascii_hash{"\x{00CE}"} = 'I';
$unicode_to_ascii_hash{"\x{00CF}"} = 'I';
$unicode_to_ascii_hash{"\x{00D0}"} = 'D';
$unicode_to_ascii_hash{"\x{00D1}"} = 'N';
$unicode_to_ascii_hash{"\x{00D2}"} = 'O';
$unicode_to_ascii_hash{"\x{00D3}"} = 'O';
$unicode_to_ascii_hash{"\x{00D4}"} = 'O';
$unicode_to_ascii_hash{"\x{00D5}"} = 'O';
$unicode_to_ascii_hash{"\x{00D6}"} = 'O';
$unicode_to_ascii_hash{"\x{00D7}"} = 'x';    # multiplication
$unicode_to_ascii_hash{"\x{00D8}"} = 'O';
$unicode_to_ascii_hash{"\x{00D9}"} = 'U';
$unicode_to_ascii_hash{"\x{00DA}"} = 'U';
$unicode_to_ascii_hash{"\x{00DB}"} = 'U';
$unicode_to_ascii_hash{"\x{00DC}"} = 'U';
$unicode_to_ascii_hash{"\x{00DD}"} = 'Y';
$unicode_to_ascii_hash{"\x{00DE}"} = 'TH';
$unicode_to_ascii_hash{"\x{00DF}"} = 'ss';
$unicode_to_ascii_hash{"\x{00E0}"} = 'a';
$unicode_to_ascii_hash{"\x{00E1}"} = 'a';
$unicode_to_ascii_hash{"\x{00E2}"} = 'a';
$unicode_to_ascii_hash{"\x{00E3}"} = 'a';
$unicode_to_ascii_hash{"\x{00E4}"} = 'a';
$unicode_to_ascii_hash{"\x{00E5}"} = 'a';
$unicode_to_ascii_hash{"\x{00E6}"} = 'ae';
$unicode_to_ascii_hash{"\x{00E7}"} = 'c';
$unicode_to_ascii_hash{"\x{00E8}"} = 'e';
$unicode_to_ascii_hash{"\x{00E9}"} = 'e';
$unicode_to_ascii_hash{"\x{00EA}"} = 'e';
$unicode_to_ascii_hash{"\x{00EB}"} = 'e';
$unicode_to_ascii_hash{"\x{00EC}"} = 'i';
$unicode_to_ascii_hash{"\x{00ED}"} = 'i';
$unicode_to_ascii_hash{"\x{00EE}"} = 'i';
$unicode_to_ascii_hash{"\x{00EF}"} = 'i';
$unicode_to_ascii_hash{"\x{00F0}"} = 'd';
$unicode_to_ascii_hash{"\x{00F1}"} = 'n';
$unicode_to_ascii_hash{"\x{00F2}"} = 'o';
$unicode_to_ascii_hash{"\x{00F3}"} = 'o';
$unicode_to_ascii_hash{"\x{00F4}"} = 'o';
$unicode_to_ascii_hash{"\x{00F5}"} = 'o';
$unicode_to_ascii_hash{"\x{00F6}"} = 'o';
$unicode_to_ascii_hash{"\x{00F7}"} = '/';    # division
$unicode_to_ascii_hash{"\x{00F8}"} = 'o';
$unicode_to_ascii_hash{"\x{00F9}"} = 'u';
$unicode_to_ascii_hash{"\x{00FA}"} = 'u';
$unicode_to_ascii_hash{"\x{00FB}"} = 'u';
$unicode_to_ascii_hash{"\x{00FC}"} = 'u';
$unicode_to_ascii_hash{"\x{00FD}"} = 'y';
$unicode_to_ascii_hash{"\x{00FE}"} = 'th';
$unicode_to_ascii_hash{"\x{00FF}"} = 'y';
# OE from Latin Extended-A
$unicode_to_ascii_hash{"\x{0152}"} = 'OE';
$unicode_to_ascii_hash{"\x{0153}"} = 'oe';
}
# only the most common ones, from synonyms field
# name field can contain junk, which I may want to try to salvage later
else
{
$unicode_to_ascii_hash{"\x{00E4}"} = 'a';
$unicode_to_ascii_hash{"\x{00E8}"} = 'e';
$unicode_to_ascii_hash{"\x{00E9}"} = 'e';
$unicode_to_ascii_hash{"\x{00EF}"} = 'i';
$unicode_to_ascii_hash{"\x{00F4}"} = 'o';
$unicode_to_ascii_hash{"\x{00F6}"} = 'o';
$unicode_to_ascii_hash{"\x{00FC}"} = 'u';
}


# +/-, dashes, and quotes
$unicode_to_ascii_hash{"\x{00B1}"} = '+-';   # synonyms use (+-), not (+/-)
$unicode_to_ascii_hash{"\x{2010}"} = '-';    # name only, real (50 of them)
$unicode_to_ascii_hash{"\x{2011}"} = '-';
$unicode_to_ascii_hash{"\x{2012}"} = '-';
$unicode_to_ascii_hash{"\x{2013}"} = '-';
$unicode_to_ascii_hash{"\x{2014}"} = '-';
$unicode_to_ascii_hash{"\x{2015}"} = '-';
$unicode_to_ascii_hash{"\x{2192}"} = '-';
$unicode_to_ascii_hash{"\x{2212}"} = '-';
#$unicode_to_ascii_hash{"\x{00A8}"} = '"';    # name only, corrupt text
$unicode_to_ascii_hash{"\x{02B9}"} = "'";
$unicode_to_ascii_hash{"\x{02BA}"} = '"';
$unicode_to_ascii_hash{"\x{2018}"} = "'";
$unicode_to_ascii_hash{"\x{2019}"} = "'";
#$unicode_to_ascii_hash{"\x{201A}"} = "'";    # name only, corrupt text
$unicode_to_ascii_hash{"\x{201B}"} = "'";
$unicode_to_ascii_hash{"\x{201C}"} = '"';
$unicode_to_ascii_hash{"\x{201D}"} = '"';
$unicode_to_ascii_hash{"\x{201E}"} = '"';
$unicode_to_ascii_hash{"\x{201F}"} = '"';
$unicode_to_ascii_hash{"\x{2032}"} = "'";
$unicode_to_ascii_hash{"\x{2033}"} = '"';
$unicode_to_ascii_hash{"\x{2035}"} = "'";
$unicode_to_ascii_hash{"\x{2034}"} = "'''";
$unicode_to_ascii_hash{"\x{2036}"} = '"';
$unicode_to_ascii_hash{"\x{2037}"} = "'''";   # triple prime
$unicode_to_ascii_hash{"\x{2057}"} = '""';    # quadruple prime
$unicode_to_ascii_hash{"\x{301D}"} = '"';
$unicode_to_ascii_hash{"\x{301E}"} = '"';
$unicode_to_ascii_hash{"\x{301F}"} = '"';

# other punctuation
# should be removed, shouldn't be there in the first place
$unicode_to_ascii_hash{"\x{00AB}"} = '';   # <<;  alpha,<<gamma>>-butadiene
$unicode_to_ascii_hash{"\x{00BB}"} = '';   # >>;  alpha,<<gamma>>-butadiene
$unicode_to_ascii_hash{"\x{2020}"} = '';   # dagger;        end of HMDB0240697
$unicode_to_ascii_hash{"\x{2021}"} = '';   # double dagger; end of HMDB0240697


# superscript numbers
$unicode_to_ascii_hash{"\x{00B2}"} = '2';
$unicode_to_ascii_hash{"\x{00B3}"} = '3';
$unicode_to_ascii_hash{"\x{00B9}"} = '1';
$unicode_to_ascii_hash{"\x{2070}"} = '0';
$unicode_to_ascii_hash{"\x{2074}"} = '4';
$unicode_to_ascii_hash{"\x{2075}"} = '5';
$unicode_to_ascii_hash{"\x{2076}"} = '6';
$unicode_to_ascii_hash{"\x{2077}"} = '7';
$unicode_to_ascii_hash{"\x{2078}"} = '8';
$unicode_to_ascii_hash{"\x{2079}"} = '9';
# subscript numbers
$unicode_to_ascii_hash{"\x{2080}"} = '0';
$unicode_to_ascii_hash{"\x{2081}"} = '1';
$unicode_to_ascii_hash{"\x{2082}"} = '2';
$unicode_to_ascii_hash{"\x{2083}"} = '3';
$unicode_to_ascii_hash{"\x{2084}"} = '4';
$unicode_to_ascii_hash{"\x{2085}"} = '5';
$unicode_to_ascii_hash{"\x{2086}"} = '6';
$unicode_to_ascii_hash{"\x{2087}"} = '7';
$unicode_to_ascii_hash{"\x{2088}"} = '8';
$unicode_to_ascii_hash{"\x{2089}"} = '9';
# zero-width spaces, should be removed
$unicode_to_ascii_hash{"\x{200B}"} = '';
$unicode_to_ascii_hash{"\x{FEFF}"} = '';    # name only, Bosutinib

# U+00A0 (non-breaking space)
#
# HMDB0062476
#    GalNAc(3S)-GlcA-Gal-Gal-Xyl??
# Non-breaking space at the end
# This occurs several times, likely a copy/paste error
$unicode_to_ascii_hash{"\x{00A0}"} = '';

# U+00AC (NOT symbol)
#
# appears to be inserted junk in front of +/-
# example: HMDB0303381 (+/-)-Isobornyl acetate
#          https://foodb.ca/compounds/FDB012445
# Or part of corrupted multibyte unicode:
#     HMDB0251069
#     HMDB0250632

# bogus character, probably supposed to be alpha?
# it only occurs once in all of the HMDB compound name fields
# HMDB0242122 ?-D-galactopyranoside, ethyl
$unicode_to_ascii_hash{"\x{FFFD}"} = 'a';



sub unicode_to_ascii
{
    my $value = $_[0];
    my $len;
    my $string_new;
    
    if ($value =~ /[\x80-\xFF]/)
    {
        # first, decode the unicode string
        # into single characters, so substr works correctly
        utf8::decode($value);

        #if ($value =~ /[\x{0370}-\x{03ff}]/)
        #{
        #    $temp = $value;
        #    utf8::encode($temp);
        #    printf STDERR "$accession\t$temp\n";
        #}

        ## HACK -- {NOT}+/-
        $value =~ s/\(\x{00AC}\x{00B1}\)/\(\x{00B1}\)/g;
        
        ## HACK -- HMDB0304547
        ## corrupted omega
        ## thank you python ftfy package for confirming!
        $value =~ s/\x{0153}\x{00E2}/\x{03C9}/g;
        
        ## HACK -- HMDB0304570 HMDB0304569
        ## appears to be corrupted '
        ## thank you python ftfy package!
        $value =~ s/\x{201A}\x{00C4}\x{2264}/\'/g;
        
        ## HMDB0300900    (2E_4Z)???\decadienoyl-CoA
        ##                http://qpmf.rx.umaryland.edu/PAMDB?MetID=PAMDB001410
        ##                should probably be (2E_4Z)-decadienoyl-CoA
        ## remove offending unicode entirely, plus the following backslash
        $value =~ s/\x{00D4}\x{00F8}\x{03A9}\\//g;

        ## HMDB0251069    2,2???-(Hydroxynitrosohydrazino)bis-ethanamine
        ##                bloodexposome.org: 2,2'-(Hydroxynitrosohydrazino)...
        ## remove offending unicode entirely
        $value =~ s/\x{201A}\x{00C4}\x{00F6}\x{221A}\x{00D1}\x{221A}\x{2202}\?//g;
        $value =~ s/\x{201A}\x{00C4}\x{00F6}\x{221A}\x{00A2}\x{00AC}\x{00DF}//g;
        
        ## HMDB0250632    ... cyclic (3?5)-disulfide
        ##                bloodexposome.org: ... cyclic (35)-disulfide
        ## remove offending unicode entirely
        $value =~ s/\x{00AC}\x{00A8}\x{00AC}\x{00AE}//g;
        $value =~ s/\x{201A}\x{00E0}\x{00F6}\x{221A}\x{00FA}//g;


        $string_new = '';
        $len        = length $value;
        for ($j = 0; $j < $len; $j++)
        {
            $c_single = substr $value, $j, 1;
            $c_new = $unicode_to_ascii_hash{$c_single};

            if (defined($c_new))
            {
                $string_new .= $c_new;
            }
            else
            {
                $string_new .= $c_single;
            }
        }
        $value = $string_new;
        
        # remove leading/trailing whitespace
        $value =~ s/^\s+//;
        $value =~ s/\s+$//;

        # encode it back again, since input is multi-byte chars
        utf8::encode($value);
    }
    
    return $value;
}

sub conform_name
{
    my $name       = $_[0];
    my $name_len;
    my $half_name_len;
    my $name_half1 = '';
    my $name_half2 = '';

    my @temp_array;
    my $i;
    
    # lowercase everything
    $name = lc $name;
    
    # convert unicode to ASCII
    $name = unicode_to_ascii($name);
    
    # convert _ to space, so \w and \b won't match on them
    $name =~ s/_/ /g;
    
    # remove likely abbreviations or synonyms from end
    $name =~ s/ +\([^0-9()]+\)$//;

    # replace Greek letters at letter boundaries
    # example: HMDB0000708 Glycoursodeoxycholic acid
    $name =~ s/(?<![A-Za-z])alpha(?![A-Za-z])/a/g;
    $name =~ s/(?<![A-Za-z])beta(?![A-Za-z])/b/g;
    $name =~ s/(?<![A-Za-z])gamma(?![A-Za-z])/g/g;
    $name =~ s/(?<![A-Za-z])delta(?![A-Za-z])/d/g;
    $name =~ s/(?<![A-Za-z])epsilon(?![A-Za-z])/e/g;
    $name =~ s/(?<![A-Za-z])zeta(?![A-Za-z])/z/g;
    $name =~ s/(?<![A-Za-z])eta(?![A-Za-z])/h/g;
    
    # replace single numbers with romanized greek letters
    #    2-aminoethylphosphonate --> b-aminoethylphosphonate
    @temp_array = split /\b([1-7])\b/, $name;
    for ($i = 1; $i < @temp_array; $i += 2)
    {
        $temp_array[$i] = $number_letter_hash{$temp_array[$i]};
    }
    $name = join '', @temp_array;

    # conform acids
    # example of "acids": HMDB0001202
    #
    $name =~ s/\bacids\b/acid/g;
    $name =~ s/acids*,(.*)ic$/$1ic acid/; # reorder weird MeSH, HMDB notation
    $name =~ s/anoic/yric/g;              # Butanoic acid --> Butyric acid
    $name =~ s/anoate\b/yrate/g;          # Butanoate     --> Butyrate
    $name =~ s/ic acids*\b/ate/g;         # Glutamic acid --> Glutamate
    $name =~ s/ates\b/ate/g;              # Benzoates     --> Benzoate
    
    # strip the L- from L-aminoacids
    # only on word boundary, so we don't strip DL-aminoacid
    $name =~ s/\bl-(alanine)/$1/g;
    $name =~ s/\bl-(arginine)/$1/g;
    $name =~ s/\bl-(asparagine)/$1/g;
    $name =~ s/\bl-(aspartic acid)/$1/g;
    $name =~ s/\bl-(cysteine)/$1/g;
    $name =~ s/\bl-(glutamine)/$1/g;
    $name =~ s/\bl-(glutamic acid)/$1/g;
    $name =~ s/\bl-(glycine)/$1/g;    # L-Glycine is used, but not by itself
    $name =~ s/\bl-(histidine)/$1/g;
    $name =~ s/\bl-(isoleucine)/$1/g;
    $name =~ s/\bl-(leucine)/$1/g;
    $name =~ s/\bl-(lysine)/$1/g;
    $name =~ s/\bl-(methionine)/$1/g;
    $name =~ s/\bl-(phenylalanine)/$1/g;
    $name =~ s/\bl-(proline)/$1/g;
    $name =~ s/\bl-(serine)/$1/g;
    $name =~ s/\bl-(threonine)/$1/g;
    $name =~ s/\bl-(tryptophan)/$1/g;          # also Tryptophanamide
    $name =~ s/\bl-(tyrosine)/$1/g;
    $name =~ s/\bl-(valine)/$1/g;

    # sulfid/sulfide/sulphid/sulphide
    # HMDB0042033 Thiodiglycol is the only entry with sulfid/sulphid
    # so, replace sulfid/sulphid with sulfide, since sulfid/sulphid is odd
    # sulfide is kept over sulphide due to fewer characters
    #
    $name =~ s/sulph/sulf/g;
    $name =~ s/sulfid\b/sulfide/g;
    #$name =~ s/sulphid\b/sulfide/g;
    #$name =~ s/sulphide\b/sulfide/g;

    # "ic acid" / "ate"
    $name =~ s/\bl-(aspartate)\b/$1/g;
    $name =~ s/\bl-(glutamate)\b/$1/g;

    # artificial, non-Human amino acids or dipeptides
    $name =~ s/\bl-(cysteinylglycine)/$1/g;    # Cysteinylglycine
    $name =~ s/\bl-(homocysteine)/$1/g;        # Homocysteine
    $name =~ s/\bl-(norleucine)/$1/g;          # Norleucine
    $name =~ s/\bl-(selenomethionine)/$1/g;    # Selenomethionine
    $name =~ s/\bl-(anserine)/$1/g;            # Anserine
    $name =~ s/\bl-(homoserine)/$1/g;          # Homoserine
    $name =~ s/\bl-(allothreonine)/$1/g;       # Allothreonine
    $name =~ s/\bl-(norvaline)/$1/g;           # Norvaline
    
    # ethyl, methyl, etc.
    $name =~ s/thane/thyl/g;          # (2-Aminoethane)phosphonic acid
                                      # 2-Aminoethylphosphonate

    # mono is redundant and often left out in synonyms
    #
    # WARNING:  conform to the same string, but methyl in different places
    #   HMDB0000752  Methylglutaric acid
    #   HMDB0000858  Monomethyl glutaric acid
    #
    $name =~ s/mono(\S)/$1/g;
    
    # replace (+/-) with (+-)
    $name =~ s/\(\+\/\-\)/\(\+\-\)/g;

    # keep only D,L,DL when together with (+),(-),(+-)
    #
    # examples: D-(-)-Arabinose
    #           D-(+)-Galactosamine
    #           D-(+)-Galacturonic acid
    #           D-(+)-Glucosamine
    #
    #swap them around if +/- comes before D/L
    #$name =~ s/(?<!\w)(d|l|dl)(?!\w)(.*?)(?<!\w)\((\+|\-|\+\-)\)(?!\w)/$1$2/g;
    #$name =~ s/(?<!\w)\((\+|\-|\+\-)\)(?!\w)(.*?)(?<!\w)(d|l|dl)(?!\w)/$3$2/g;
    $name =~ s/(?<!\w)(d|l|dl)-*\((\+|\-|\+\-)\)(?!\w)/$1/g;
    $name =~ s/(?<!\w)\((\+|\-|\+\-)\)-*(d|l|dl)(?!\w)/$2/g;

    ## protect -) as in (+/-) or (-) using capital letters    
    $name =~ s/[\(\[\{]-|-[\)\]\}]/MINUS/g;       # protect minus signs
    
    
    ## condense everything that isn't a letter, number, comma, or plus/minus
    ## except when between two numbers
    ##
    $name =~ s/[^,+\w]/-/g;           # convert to hyphens
    $name =~ s/-+/-/g;                # condense multiple hyphens in a row
    $name =~ s/(^-|-$)//g;            # strip leading/trailing hyphens
    $name =~ s/([,+])-/$1/g;          # strip hyphens next to comma or plus
    $name =~ s/-([,+])/$1/g;          # strip hyphens next to comma or plus
    $name =~ s/([a-z])-/$1/g;         # strip hyphens next to letters
    $name =~ s/-([a-z])/$1/g;         # strip hyphens near to letters
    $name =~ s/,+/,/g;                # condense multiple commas in a row
    $name =~ s/\++/\+/g;              # condense multiple pluses in a row
    
    ## de-protect and condense minus signs
    ##
    ## hypothetical example: 1-(-)-galacturonate   -->   1-galacturonate
    ##                       1-(+)-galacturonate   -->   1+galacturonate
    ##                       1-galacturonate       -->   1galacturonate
    $name =~ s/MINUS/-/g;
    $name =~ s/-+/-/g;

    # check for tandem duplicate names after conforming
    $name_len = length $name;
    if ($name_len % 2 == 0)
    {
        $half_name_len = 0.5 * $name_len;
        $name_half1 = substr $name, 0, $half_name_len;
        $name_half2 = substr $name, $half_name_len, $half_name_len;
        
        if ($name_half1 eq $name_half2)
        {
            $name = $name_half1;
        }
    }

    return $name;
}



# begin main()

$syntax_error_flag         = 0;
$num_files                 = 0;

for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        # override default PPM tolerance
        if ($field eq '--ppm' ||
            $field =~ /^--ppm\=/)
        {
            $ppm_arg = '';
            
            if ($field eq '--ppm')
            {
                if ($i + 1 < @ARGV)
                {
                    $i++;
                    $ppm_arg = $ARGV[$i];
                }
            }
            elsif ($field =~ /^--ppm=(.*)/)
            {
                $ppm_arg = $1;
            }
            
            if (is_number($ppm_arg))
            {
                $mz_tol_ppm = $ppm_arg;
                
                printf STDERR "Overiding default m/z PPM tolerance of 10:\t%s\n",
                    $mz_tol_ppm;
            }
            else
            {
                printf STDERR "ABORT -- you must specify a numeric value after --ppm\n";
                $syntax_error_flag = 1;
            }
        }
        else
        {
            printf STDERR "ABORT -- unknown option %s\n", $field;
            $syntax_error_flag = 1;
        }
    }
    else
    {
        if ($num_files == 0)
        {
            $annotation_filename = $field;
            $num_files++;
        }
        elsif ($num_files == 1)
        {
            $data_filename = $field;
            $num_files++;
        }
    }
}



if (!defined($data_filename) || !defined($annotation_filename) ||
    $syntax_error_flag)
{
    print STDERR "Usage: annotate_metabolomics.pl identifier_mapping.txt cleaned_mzmine.txt\n";
    print STDERR "\n";
    print STDERR "Options:\n";
    print STDERR "    --ppm N          override default m/z PPM tolerance\n";
    exit(1);
}



open ANNOTATION, "$annotation_filename" or die "can't open $annotation_filename\n";
open DATA,       "$data_filename"       or die "can't open $data_filename\n";


# read in annotation file

# skip down to first line that has anything on it
while($line=<ANNOTATION>)
{
    # skip comment lines
    if ($line =~ /^#/)
    {
        next;
    }

    # this line isn't purely whitespace, assume it is the header line
    if ($line =~ /\S/)
    {
        last;
    }
}


# annotation header line
$line =~ s/[\r\n]+//g;
$line =~ s/\"//g;
@array = split /\t/, $line;    # skip empty headers at and
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    $array[$i] =~ s/\s+/ /g;

    $field = $array[$i];
    $annotation_header_col_hash{$field} = $i;
    # $annotation_header_col_array[$i] = $field;
}

# important annotation column headers
for ($i = 0; $i < @array; $i++)
{
    $header = $array[$i];

    # skip Conflict columns
    if ($header =~ /conflict/i)
    {
        next;
    }
    
    if (!defined($annotation_mz_col) &&
        $header =~ /m\/z/i)
    {
        $annotation_mz_col = $i;
    }
    elsif (!defined($annotation_formula_col) &&
           $header =~ /formula/i)
    {
        $annotation_formula_col = $i;
    }
    elsif (!defined($annotation_name_col) &&
           $header =~ /identity/i)
    {
        $annotation_name_col = $i;
    }
    elsif (!defined($annotation_synonym_col) &&
           $header =~ /synonym$/i)
    {
        $annotation_synonym_col = $i;
    }
    elsif (!defined($annotation_traditional_col) &&
           $header =~ /traditional_iupac$/i)
    {
        $annotation_traditional_col = $i;
    }
    elsif (!defined($annotation_kegg_map_col) &&
           $header =~ /\bkegg_map/i)
    {
        $annotation_kegg_map_col = $i;
    }
    elsif (!defined($annotation_kegg_col) &&
           $header =~ /\bkegg/i)
    {
        $annotation_kegg_col = $i;
    }
    elsif (!defined($annotation_hmdb_col) &&
           $header =~ /\bhmdb/i)
    {
        $annotation_hmdb_col = $i;
    }
    elsif (!defined($annotation_pubchem_col) &&
           $header =~ /\bpubchem/i)
    {
        $annotation_pubchem_col = $i;
    }
    elsif (!defined($annotation_pos_neg_col) &&
           $header =~ /pos.*neg/i)
    {
        $annotation_pos_neg_col = $i;
    }
    elsif (!defined($annotation_rt_col) &&
           $header =~ /retention time/i)
    {
        $annotation_rt_col = $i;
    }
}

# use name column, if present, instead of any other column for name
if (defined($annotation_header_col_hash{'name'}))
{
    $annotation_name_col = $annotation_header_col_hash{'name'};
}

$annotation_adduct_flag   = 0;
$annotation_mz_pos_col    = $annotation_header_col_hash{'m/z_pos'};
$annotation_mz_neg_col    = $annotation_header_col_hash{'m/z_neg'};
$annotation_mono_mass_col = $annotation_header_col_hash{'mono_mass'};
$annotation_acc_col       = $annotation_header_col_hash{'accession'};

if ((defined($annotation_mz_pos_col) ||
     defined($annotation_mz_neg_col)) &&
    defined($annotation_mono_mass_col) &&
    defined($annotation_acc_col))
{
    $annotation_adduct_flag = 1;
}

# replace HMDB column with accession column,
# since these only come from HMDB at the moment
if ($annotation_adduct_flag)
{
    $annotation_hmdb_col = $annotation_acc_col;
}


if (!defined($annotation_mz_col))
{
    printf STDERR "ABORT -- m/z column not found in annotation file %s\n",
                 $annotation_filename;
    exit(1);
}
if (!defined($annotation_formula_col))
{
    printf STDERR "ABORT -- formula column not found in annotation file %s\n",
                 $annotation_filename;
    exit(1);
}
if (!defined($annotation_name_col))
{
    printf STDERR "ABORT -- name/identity column not found in annotation file %s\n",
                 $annotation_filename;
    exit(1);
}
#if (!defined($annotation_kegg_col))
#{
#    printf STDERR "ABORT -- KEGG ID column not found in annotation file %s\n",
#                 $annotation_filename;
#    exit(1);
#}
#if (!defined($annotation_hmdb_col))
#{
#    printf STDERR "ABORT -- HMDB column not found in annotation file %s\n",
#                 $annotation_filename;
#    exit(1);
#}
#if (!defined($annotation_pubchem_col))
#{
#    printf STDERR "ABORT -- PubChem column not found in annotation file %s\n",
#                 $annotation_filename;
#    exit(1);
#}


# columns to exclude from count of interesting columns
%annotation_no_count_col_hash = ();
if (defined($annotation_mz_col))
{
    $annotation_no_count_col_hash{$annotation_mz_col} = 1;
}
if (defined($annotation_formula_col))
{
    $annotation_no_count_col_hash{$annotation_formula_col} = 1;
}
if (defined($annotation_name_col))
{
    $annotation_no_count_col_hash{$annotation_name_col} = 1;
}
if (defined($annotation_synonym_col))
{
    $annotation_no_count_col_hash{$annotation_synonym_col} = 1;
}
if (defined($annotation_traditional_col))
{
    $annotation_no_count_col_hash{$annotation_traditional_col} = 1;
}
if (defined($annotation_pos_neg_col))
{
    $annotation_no_count_col_hash{$annotation_pos_neg_col} = 1;
}
if (defined($annotation_rt_col))
{
    $annotation_no_count_col_hash{$annotation_rt_col} = 1;
}
if (defined($annotation_name_col))
{
    $annotation_no_count_col_hash{$annotation_name_col} = 1;
}
if (defined($annotation_mz_pos_col))
{
    $annotation_no_count_col_hash{$annotation_mz_pos_col} = 1;
}
if (defined($annotation_mz_neg_col))
{
    $annotation_no_count_col_hash{$annotation_mz_neg_col} = 1;
}
if (defined($annotation_mono_mass_col))
{
    $annotation_no_count_col_hash{$annotation_mono_mass_col} = 1;
}
if (defined($annotation_acc_col))
{
    $annotation_no_count_col_hash{$annotation_acc_col} = 1;
}
if (defined($annotation_kegg_map_col))
{
    $annotation_no_count_col_hash{$annotation_kegg_map_col} = 1;
}


# read in the annotation
@conformed_name_hash = ();
$row = 0;
while(defined($line=<ANNOTATION>))
{
    $line =~ s/[\r\n]+//g;
    $line =~ s/\"//g;

    @array = split /\t/, $line, -1;    # don't skip empty fields at and

    # clean up fields
    for ($col = 0; $col < @array; $col++)
    {
        $array[$col] =~ s/^\s+//;
        $array[$col] =~ s/\s+$//;
        $array[$col] =~ s/\s+/ /g;
        
#        $array[$col] = reformat_sci($array[$col]);
    }

    $mz          = $array[$annotation_mz_col];
    $formula     = $array[$annotation_formula_col];
    $name        = $array[$annotation_name_col];

    $kegg        = '';
    $kegg_map    = '';
    $hmdb        = '';
    $pubchem     = '';
    $synonym_str = '';
    $traditional_str = '';
    
    if (defined($annotation_synonym_col))
    {
        $synonym_str = $array[$annotation_synonym_col];
    }
    if (defined($annotation_traditional_col))
    {
        $traditional_str = $array[$annotation_traditional_col];
    }
    if (defined($annotation_kegg_col))
    {
        $kegg     = $array[$annotation_kegg_col];
    }
    if (defined($annotation_kegg_map_col))
    {
        $kegg_map = $array[$annotation_kegg_map_col];
    }
    if (defined($annotation_hmdb_col))
    {
        $hmdb     = $array[$annotation_hmdb_col];
    }
    if (defined($annotation_pubchem_col))
    {
        $pubchem  = $array[$annotation_pubchem_col];
    }
    
    # retention time sanity checks, if available
    if (defined($annotation_rt_col))
    {
        $rt  = $array[$annotation_rt_col];
    }
    else
    {
        $rt  = '';
    }
    
    if (!defined($mz))       { $mz       = ''; }
    if (!defined($formula))  { $formula  = ''; }
    if (!defined($name))     { $name     = ''; }
    if (!defined($kegg))     { $kegg     = ''; }
    if (!defined($kegg_map)) { $kegg_map = ''; }
    if (!defined($hmdb))     { $hmdb     = ''; }
    if (!defined($pubchem))  { $pubchem  = ''; }
    if (!defined($pubchem))  { $rt       = ''; }

    # skip bad rows
    if (!($name =~ /\S/))    { next; }


    # count interesting columns
    $count = 0;
    for ($col = 0; $col < @array; $col++)
    {
        if (!defined($annotation_no_count_col_hash{$col}) &&
            $array[$col] =~ /\S/)
        {
            $count++;
        }
    }

    # sanity check formula, stripped of H's
    $fsanity = '';
    if ($formula =~ /[A-Za-z0-9]/)
    {
        $fsanity = conform_formula($formula);
        
        # remove H's, since protonation state is questionable
        $fsanity =~ s/(D|H)[0-9]*//g;
    }
    
    $conformed_name_hash{$name} = conform_name($name);

    # concatenate synonym and traditional iupac fields
    if ($traditional_str ne '')
    {
        if ($synonym_str ne '')
        {
            $synonym_str .= '|' . $traditional_str;
        }
        else
        {
            $synonym_str = $traditional_str;
        }
    }
    
    @synonym_array = split /\|/, $synonym_str;
    %synonym_hash  = ();
    for ($i = 0; $i < @synonym_array; $i++)
    {
        $synonym = conform_name($synonym_array[$i]);
        $synonym_hash{$synonym} = 1;
    }
    @synonym_array = sort keys %synonym_hash;

    $mz = bless_delimiter_bar_metabolomics($mz);
    @mz_array = split /\|/, $mz;
    foreach $mz (@mz_array)
    {
        # skip bad m/z
        if (!($mz =~ /[0-9]/)) { next; }

          $annotation_hash{$row}{formula}   = $formula;
          $annotation_hash{$row}{name}      = $name;
        @{$annotation_hash{$row}{syn_arr}}  = @synonym_array;
          $annotation_hash{$row}{kegg}      = $kegg;
          $annotation_hash{$row}{kegg_map}  = $kegg_map;
          $annotation_hash{$row}{hmdb}      = $hmdb;
          $annotation_hash{$row}{pubchem}   = $pubchem;
          $annotation_hash{$row}{rt}        = $rt;
          $annotation_hash{$row}{mz}        = $mz;
          $annotation_hash{$row}{col_count} = $count;
          $annotation_hash{$row}{fsanity}   = $fsanity;
        
        # bin m/z for rapid scanning later
        # largest m/z we see is 900
        # +/-1 is ~1000 ppm, which is *way* over our 10 ppm tolerance
        # bins are too large, but are very simple to code up
        #
        # also factor in +/- 2.014552 (pos - neg) difference
        #
        $ppm_offset = $mz_tol_ppm * $mz / 1000000;
        $mz_floor   = floor($mz);
        $mz_minus1  = floor($mz - $ppm_offset);
        $mz_plus1   = floor($mz + $ppm_offset);
        
        $mz_row_bins_hash{$mz_minus1}{$row} = 1;
        $mz_row_bins_hash{$mz_floor}{$row}  = 1;
        $mz_row_bins_hash{$mz_plus1}{$row}  = 1;

        # determine pos/neg ionization mode
        $pos_flag = 0;
        $neg_flag = 0;
        if (defined($annotation_pos_neg_col))
        {
            $pos_neg = $array[$annotation_pos_neg_col];

            if ($pos_neg =~ /pos/i) { $pos_flag = 1; }
            if ($pos_neg =~ /neg/i) { $neg_flag = 1; }
        }

        if ($annotation_adduct_flag)
        {
            $pos_flag = 0;
            $neg_flag = 0;

            if ($annotation_mz_col == $annotation_mz_pos_col)
            {
                $pos_flag = 1;
            }
            if ($annotation_mz_col == $annotation_mz_neg_col)
            {
                $neg_flag = 1;
            }
        }
        
        if ($pos_flag == 1 && $neg_flag == 0)
        {
            $annotation_hash{$row}{pos_neg} = 'pos';
        }
        if ($pos_flag == 0 && $neg_flag == 1)
        {
            $annotation_hash{$row}{pos_neg} = 'neg';
        }

        if (1 || $annotation_adduct_flag == 0)
        {
            $mz_minus3 = floor($mz - 2.014552 - $ppm_offset);
            $mz_minus2 = floor($mz - 2.014552 + $ppm_offset);
            $mz_plus2  = floor($mz + 2.014552 - $ppm_offset);
            $mz_plus3  = floor($mz + 2.014552 + $ppm_offset);

            # include -2.014552 if pos
            if ($pos_flag)
            {
                $mz_row_bins_hash{$mz_minus3}{$row} = 1;
                $mz_row_bins_hash{$mz_minus2}{$row} = 1;
            }
            # include +2.014552 if neg
            if ($neg_flag)
            {
                $mz_row_bins_hash{$mz_plus2}{$row}  = 1;
                $mz_row_bins_hash{$mz_plus3}{$row}  = 1;
            }

            # if neither, include both
            if ($pos_flag == 0 && $neg_flag == 0)
            {
                $mz_row_bins_hash{$mz_minus3}{$row} = 1;
                $mz_row_bins_hash{$mz_minus2}{$row} = 1;
                $mz_row_bins_hash{$mz_plus2}{$row}  = 1;
                $mz_row_bins_hash{$mz_plus3}{$row}  = 1;
            }
        }

        $row++;
    }

    # assume one of the pos/neg mz cols was already stored
    # store the other one
    if ($annotation_adduct_flag)
    {
        $temp_mz_col = $annotation_mz_pos_col;
        $pos_flag    = 1;
        $neg_flag    = 0;
        if ($annotation_mz_col == $annotation_mz_pos_col)
        {
            $temp_mz_col = $annotation_mz_neg_col;
            $pos_flag    = 0;
            $neg_flag    = 1;
        }
        $mz = $array[$annotation_mz_neg_col];
    
        $mz = bless_delimiter_bar_metabolomics($mz);
        @mz_array = split /\|/, $mz;
        foreach $mz (@mz_array)
        {
            # skip bad m/z
            if (!($mz =~ /[0-9]/)) { next; }

              $annotation_hash{$row}{formula}   = $formula;
              $annotation_hash{$row}{name}      = $name;
            @{$annotation_hash{$row}{syn_arr}}  = @synonym_array;
              $annotation_hash{$row}{kegg}      = $kegg;
              $annotation_hash{$row}{kegg_map}  = $kegg_map;
              $annotation_hash{$row}{hmdb}      = $hmdb;
              $annotation_hash{$row}{pubchem}   = $pubchem;
              $annotation_hash{$row}{rt}        = $rt;
              $annotation_hash{$row}{mz}        = $mz;
              $annotation_hash{$row}{col_count} = $count;
              $annotation_hash{$row}{fsanity}   = $fsanity;

            if ($pos_flag == 1 && $neg_flag == 0)
            {
                $annotation_hash{$row}{pos_neg} = 'pos';
            }
            if ($pos_flag == 0 && $neg_flag == 1)
            {
                $annotation_hash{$row}{pos_neg} = 'neg';
            }
            
            # bin m/z for rapid scanning later
            # largest m/z we see is 900
            # +/-1 is ~1000 ppm, which is *way* over our 10 ppm tolerance
            # bins are too large, but are very simple to code up
            #
            # also factor in +/- 2.014552 (pos - neg) difference
            #
            $ppm_offset = $mz_tol_ppm * $mz / 1000000;
            $mz_floor   = floor($mz);
            $mz_minus1  = floor($mz - $ppm_offset);
            $mz_plus1   = floor($mz + $ppm_offset);
            
            $mz_row_bins_hash{$mz_minus1}{$row} = 1;
            $mz_row_bins_hash{$mz_floor}{$row}  = 1;
            $mz_row_bins_hash{$mz_plus1}{$row}  = 1;

            if (1 || $annotation_adduct_flag == 0)
            {
                $mz_minus3 = floor($mz - 2.014552 - $ppm_offset);
                $mz_minus2 = floor($mz - 2.014552 + $ppm_offset);
                $mz_plus2  = floor($mz + 2.014552 - $ppm_offset);
                $mz_plus3  = floor($mz + 2.014552 + $ppm_offset);

                # include -2.014552 if pos
                if ($pos_flag)
                {
                    $mz_row_bins_hash{$mz_minus3}{$row} = 1;
                    $mz_row_bins_hash{$mz_minus2}{$row} = 1;
                }
                # include +2.014552 if neg
                if ($neg_flag)
                {
                    $mz_row_bins_hash{$mz_plus2}{$row}  = 1;
                    $mz_row_bins_hash{$mz_plus3}{$row}  = 1;
                }

                # if neither, include both
                if ($pos_flag == 0 && $neg_flag == 0)
                {
                    $mz_row_bins_hash{$mz_minus3}{$row} = 1;
                    $mz_row_bins_hash{$mz_minus2}{$row} = 1;
                    $mz_row_bins_hash{$mz_plus2}{$row}  = 1;
                    $mz_row_bins_hash{$mz_plus3}{$row}  = 1;
                }
            }

            $row++;
        }
    }
}
#$row_count = $row;


# all rows, to be used when looking for matches regardless of m/z
@row_all_array = sort cmp_rows keys %annotation_hash;


# bogus row for unmapped hits
$annotation_hash{$bad_row_id}{mz}       = '';
$annotation_hash{$bad_row_id}{formula}  = '';
$annotation_hash{$bad_row_id}{name}     = '';
$annotation_hash{$bad_row_id}{kegg}     = '';
$annotation_hash{$bad_row_id}{kegg_map} = '';
$annotation_hash{$bad_row_id}{hmdb}     = '';
$annotation_hash{$bad_row_id}{pubchem}  = '';
$annotation_hash{$bad_row_id}{rt}       = '';


# read in data file

# skip down to first line that has anything on it
# lipidomics data has this issue sometimes
while($line=<DATA>)
{
    # skip comment lines
    if ($line =~ /^#/)
    {
        next;
    }

    # this line isn't purely whitespace, assume it is the header line
    if ($line =~ /\S/)
    {
        last;
    }
}


# data header line
$line =~ s/[\r\n]+//g;
$line =~ s/\"//g;
@array = split /\t/, $line;    # skip empty headers at and
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    $array[$i] =~ s/\s+/ /g;

    # clean up sample names
    $array[$i] =~ s/^_+//;
    $array[$i] =~ s/_+$//;
    
    $field = $array[$i];
    $data_header_col_hash{$field} = $i;
    $data_header_col_array[$i] = $field;
}


# Excel gets very upset if first field is "ID", change it
if ($data_header_col_array[0] =~ /^id$/i)
{
    $data_header_col_array[0] = 'Index';
    $data_header_col_hash{'Index'} = 0;
}


$elmaven_isotope_flag = 0;
if (defined($data_header_col_hash{'isotopeLabel'}) &&
    defined($data_header_col_hash{'parent'}))
{
    $elmaven_isotope_flag = 1;
}


# important data column headers
for ($i = 0; $i < @array; $i++)
{
    $header = $array[$i];

    if (!defined($data_mz_col) &&
        $header =~ /m\/*z/i)
    {
        $data_mz_col = $i;
    }
    elsif (!defined($data_name_col) &&
           $header =~ /(identity|compound)/i)
    {
        $data_name_col = $i;
        
        if ($header =~ /main ID/i)
        {
            printf STDERR "WARNING -- (main ID) used instead of (all ids), will miss hits\n";
        }
    }
    elsif (!defined($data_rt_col) &&
           $header =~ /retention time/i)
    {
        $data_rt_col = $i;
    }
    # heavy labeled flag, added by the Moffitt pipeline
    # the heavy labeled detection is best left to other scripts,
    # since it can get complicated, especially for El-MAVEN
    elsif (!defined($data_heavy_col) &&
           $header =~ /Heavy-labeled flag/i)
    {
        $data_heavy_col = $i;
    }
    elsif (!defined($data_pos_neg_col) &&
           $header =~ /pos.*neg/i)
    {
        $data_pos_neg_col = $i;
    }
    elsif (!defined($data_formula_col) &&
           $header =~ /formula/i)
    {
        $data_formula_col = $i;
    }
}

# use parent m/z if it is El-MAVEN isotope data
if ($elmaven_isotope_flag)
{
    $data_mz_col = $data_header_col_hash{'parent'};
}

# lipidomics m/z
if (!defined($data_mz_col))
{
    $data_mz_col = $data_header_col_hash{'CalcMz'};
}
if (!defined($data_name_col))
{
    $data_name_col = $data_header_col_hash{'LipidIon'};
    
    if (defined($data_name_col))
    {
        $lipidomics_flag = 1;
    }
}


if (!defined($data_mz_col))
{
    printf STDERR "ABORT -- m/z column not found in data file %s\n",
                 $data_filename;
    exit(1);
}
if (!defined($data_name_col))
{
    printf STDERR "ABORT -- name/identity column not found in data file %s\n",
                 $data_filename;
    exit(1);
}


# output new header line
# insert new annotation immediately after identity/name column
$tab_flag = 0;
for ($i = 0; $i < @data_header_col_array; $i++)
{
    $header = $data_header_col_array[$i];
    
    if ($tab_flag)
    {
        print "\t";
    }
    
    printf "$header";
    
    if ($i == $data_name_col)
    {
        printf "\t%s", 'Identity Mapped';
        printf "\t%s", 'Mapping Type';
        printf "\t%s", 'FormulaMapped';
        printf "\t%s", 'KEGG';
        printf "\t%s", 'HMDB';
        printf "\t%s", 'PubChem';
    }
    
    $tab_flag = 1;
}
printf "\n";


# read in and output annotated data
$row_data = 0;
while(defined($line=<DATA>))
{
    $line =~ s/[\r\n]+//g;
    $line =~ s/\"//g;

    @array = split /\t/, $line, -1;    # don't skip empty fields at and

    # clean up fields
    for ($col = 0; $col < @array; $col++)
    {
        $array[$col] =~ s/^\s+//;
        $array[$col] =~ s/\s+$//;
        $array[$col] =~ s/\s+/ /g;
        
        $array[$col] = reformat_sci($array[$col]);
    }

    # assume first field is row identifier for error reporting
    $first_field = $array[0];

    $mz   = $array[$data_mz_col];
    $name = $array[$data_name_col];
    
    # retention time sanity checks, if available
    $rt_data = '';
    if (defined($data_rt_col))
    {
        $rt_data = $array[$data_rt_col];
    }
    
    # formula, for sanity checks, if available
    $fsanity = '';
    if (defined($data_formula_col))
    {
        $formula = $array[$data_formula_col];

        if ($formula =~ /[A-Za-z0-9]/)
        {
            $fsanity = conform_formula($formula);
        
            # remove H's, since protonation state is questionable
            $fsanity =~ s/(D|H)[0-9]*//g;
        }
    }
    
    $pos_neg = '';
    if (defined($data_pos_neg_col))
    {
        $pos_flag = 0;
        $neg_flag = 0;
        if ($array[$data_pos_neg_col] =~ /pos/i)
        {
            $pos_flag = 1;
        }
        if ($array[$data_pos_neg_col] =~ /neg/i)
        {
            $neg_flag = 1;
        }
        
        if ($pos_flag == 1 && $neg_flag == 0)
        {
            $pos_neg = 'pos';
        }
        elsif ($pos_flag == 0 && $neg_flag == 1)
        {
            $pos_neg = 'neg';
        }
    }


    # matched (mapped) rows
    %matched_row_hash  = ();
    %matched_row_type_hash = ();
    @matched_row_array = ();
    $num_matches       = 0;
    
    # keep track of originally incorrectly identified metabolites
    %bad_mz_name_hash  = ();
    %bad_mz_row_hash   = ();
    $has_bad_mz_flag   = 0;
    $has_good_mz_flag  = 0;


    # scan for match(es) in annotation database
    # skip heavy labeled metabolites
    if (is_number($mz) && $name =~ /\S/ &&
        !is_heavy_labeled($name))
    {
        %candidate_row_mz_hash         = ();
        %candidate_row_pos_neg_hash    = ();
        
        # combine the two, for synonym matching
        %candidate_row_mz_pos_neg_hash = ();

        $mz_floor = floor($mz);
        
        # scan annotations for m/z within tolerance
        #
        # also include +/- 2.014552 (pos - neg difference)
        #  already pre-binned earlier during +/- 3 binning
        @row_array = sort cmp_rows keys %{$mz_row_bins_hash{$mz_floor}};
        foreach $row (@row_array)
        {
            # sanity check the formulas
            # skip rows that have wrong numbers of non-hydrogens
            #
            $fsanity_db = $annotation_hash{$row}{fsanity};
            if ($fsanity ne '' && $fsanity_db ne '' &&
                $fsanity ne $fsanity_db)
            {
                next;
            }
        
            $pos_neg_db = $annotation_hash{$row}{pos_neg};
            if (!defined($pos_neg_db))
            {
                $pos_neg_db = '';
            }

            $mz_db  = $annotation_hash{$row}{mz};
            $ppm    = 1000000 * abs(( $mz_db             - $mz) / $mz_db);

            # check only given m/z
            # don't check any rows for the wrong ionization mode
            if ($ppm <= $mz_tol_ppm &&
                !($pos_neg_db ne '' &&
                  $pos_neg ne '' &&
                  $pos_neg ne $pos_neg_db))
            {
                $candidate_row_mz_hash{$row}         = 1;
                $candidate_row_mz_pos_neg_hash{$row} = 1;
            }
            
            # also check +/- 2.014552, in case we only have one of pos/neg
            if (1 || $annotation_adduct_flag == 0)
            {
                $ppm_m2 = 1000000 * abs((($mz_db - 2.014552) - $mz) /
                                        ( $mz_db - 2.014552));
                $ppm_p2 = 1000000 * abs((($mz_db + 2.014552) - $mz) /
                                        ( $mz_db + 2.014552));

                if ($ppm_m2 <= $mz_tol_ppm)
                {
                    if ($pos_neg ne 'pos' && $pos_neg_db ne 'neg')
                    {
                        $candidate_row_pos_neg_hash{$row}    = 1;
                        $candidate_row_mz_pos_neg_hash{$row} = 1;
                    }
                }
                if ($ppm_p2 <= $mz_tol_ppm)
                {
                    if ($pos_neg ne 'neg' && $pos_neg_db ne 'pos')
                    {
                        $candidate_row_pos_neg_hash{$row}    = 1;
                        $candidate_row_mz_pos_neg_hash{$row} = 1;
                    }
                }
            }
        }
        
        ## original and +/- ~2 m/z row arrays
        #@row_mz_array         =
        #    sort cmp_rows keys %candidate_row_mz_hash;
        #@row_pos_neg_array    =
        #    sort cmp_rows keys %candidate_row_pos_neg_hash;
        @row_mz_pos_neg_array =
            sort cmp_rows keys %candidate_row_mz_pos_neg_hash;
        
        # support multiple ;-delimited names per data row
        $name_delim = $name;
        $name_delim = bless_delimiter_bar_metabolomics($name_delim);
        @name_array = split /\|/, $name_delim;

        # checks in less confident order:
        #
        #   1:  conformed name,     original + pos/neg m/z
        #   2:  conformed synonyms, original + pos/neg m/z
        #   3:  conformed name,     all rows
        #   4:  conformed synonyms, all rows
        #   5:  fuzzy,              original + pos/neg m/z
        #   9X: unmatched

        foreach $name_oc (@name_array)
        {
            $match_flag = 0;
            $match_type = '9X';    # 9: unmatched

            $name_lc           = lc $name_oc;
            $name_lc_conformed = conform_name($name_lc);
        
            # 1A/1B: conformed name, original + pos/neg m/z
            if ($match_flag == 0)
            {
                %candidate_name_hash = ();
                %temp_row_score_hash = ();
                @candidate_row_array = ();

                # score each row
                foreach $row (@row_mz_pos_neg_array)
                {
                    $name_db           = $annotation_hash{$row}{name};
                    $name_db_lc        = lc $name_db;
                    $name_db_conformed = $conformed_name_hash{$name_db};

                    if ($name_lc_conformed eq $name_db_conformed)
                    {
                        $frac_id = 0;
                        $score =
                            score_substring_mismatch($name_lc,
                                                     $name_db_lc,
                                                     $align_method,
                                                     \$frac_id);

                        $candidate_name_hash{$name_db_lc} = 1;

                        $temp_row_score_hash{$row} = $score;
                    }
                }
                
                @candidate_row_array =
                    sort cmp_conformed_rows keys %temp_row_score_hash;
                @temp_name_array     = keys %candidate_name_hash;

                # only 1 good candidate in database, should be right?
                if (@temp_name_array == 1)
                {
                    $match_type = '1A';
                }
                # less certain, since there are multiple candidates
                else
                {
                    $match_type = '1B';
                }

                # keep highest score (and ties) as matches
                if (@candidate_row_array)
                {
                    $best_row      = $candidate_row_array[0];
                    $best_score    = $temp_row_score_hash{$best_row};
                    $best_count    = $annotation_hash{$best_row}{col_count};
                    $best_has_kegg = $annotation_hash{$best_row}{kegg};
                    if (defined($best_has_kegg) && $best_has_kegg ne '')
                    {
                        $best_has_kegg = 1;
                    }
                    else
                    {
                        $best_has_kegg = 0;
                    }
                    # kegg id is part of a pathway map
                    if ($best_has_kegg &&
                        $annotation_hash{$best_row}{kegg_map} ne '')
                    {
                        $best_has_kegg = 2;
                    }
                }
                foreach $row (@candidate_row_array)
                {
                    $score    = $temp_row_score_hash{$row};
                    $count    = $annotation_hash{$row}{col_count};
                    $has_kegg = $annotation_hash{$row}{kegg};
                    if (defined($has_kegg) && $has_kegg ne '')
                    {
                        $has_kegg = 1;
                    }
                    else
                    {
                        $has_kegg = 0;
                    }
                    # kegg id is part of a pathway map
                    if ($has_kegg &&
                        $annotation_hash{$row}{kegg_map} ne '')
                    {
                        $has_kegg = 2;
                    }
                    
                    if ($has_kegg == $best_has_kegg &&
                        $count    == $best_count &&
                        $score    == $best_score)
                    {
                        $match_flag = 1;
                        
                        # store each row only once per row
                        if (!defined($matched_row_hash{$row}))
                        {
                            $matched_row_array[$num_matches]  = $row;
                            $matched_type_array[$num_matches] = $match_type;
                            $num_matches++;
                        }
                        $matched_row_hash{$row} = 1;
                        $matched_row_type_hash{$row}{$match_type} = 1;
                    }
                    # no longer tied
                    else
                    {
                        last;
                    }
                }
            }

            # 2A/2B: conformed,   original + pos/neg m/z, synonyms
            #
            # we need to limit m/z, otherwise weird synonyms can match:
            #
            # example:  D-lactose  C12H22O11
            #   HMDB0035312  Hebevinoside I  C44H72O13
            #
            # example:  Vitamin k2  C31H40O2
            #   HMDB0030017  Menatetrenone  C31H40O2
            #   HMDB0030020  Withanolide    C28H38O5
            #   HMDB0001892  Menadione      C11H8O2
            #
            # However, we now miss Glycochenodeoxycholate
            #   HMDB0000637  Chenodeoxycholic acid glycine conjugate
            #
            if ($match_flag == 0)
            {
                %candidate_name_hash = ();
                %temp_row_score_hash = ();
                @candidate_row_array = ();

                # score each row
                foreach $row (@row_mz_pos_neg_array)
                {
                    @synonym_array     = @{$annotation_hash{$row}{syn_arr}};
                    $name_db           = $annotation_hash{$row}{name};
                    $name_db_lc        = lc $name_db;
                    $name_db_conformed = $conformed_name_hash{$name_db};
                
                    foreach $synonym_db_conformed (@synonym_array)
                    {
                        if ($name_lc_conformed eq $synonym_db_conformed)
                        {
                            $frac_id = 0;
                            $score =
                                score_substring_mismatch($name_lc_conformed,
                                                         $name_db_conformed,
                                                         $align_method,
                                                         \$frac_id);

                            $candidate_name_hash{$name_db_lc} = 1;

                            $temp_row_score_hash{$row} = $score;
                            
                            last;
                        }
                    }
                }
                
                @candidate_row_array =
                    sort cmp_conformed_rows keys %temp_row_score_hash;
                @temp_name_array     = keys %candidate_name_hash;

                # only 1 good candidate in database, should be right?
                if (@temp_name_array == 1)
                {
                    $match_type = '2A';
                }
                # less certain, since there are multiple candidates
                else
                {
                    $match_type = '2B';
                }

                # keep highest score (and ties) as matches
                if (@candidate_row_array)
                {
                    $best_row      = $candidate_row_array[0];
                    $best_score    = $temp_row_score_hash{$best_row};
                    $best_count    = $annotation_hash{$best_row}{col_count};
                    $best_has_kegg = $annotation_hash{$best_row}{kegg};
                    if (defined($best_has_kegg) && $best_has_kegg ne '')
                    {
                        $best_has_kegg = 1;
                    }
                    else
                    {
                        $best_has_kegg = 0;
                    }
                    # kegg id is part of a pathway map
                    if ($best_has_kegg &&
                        $annotation_hash{$best_row}{kegg_map} ne '')
                    {
                        $best_has_kegg = 2;
                    }
                }
                foreach $row (@candidate_row_array)
                {
                    $score    = $temp_row_score_hash{$row};
                    $count    = $annotation_hash{$row}{col_count};
                    $has_kegg = $annotation_hash{$row}{kegg};
                    if (defined($has_kegg) && $has_kegg ne '')
                    {
                        $has_kegg = 1;
                    }
                    else
                    {
                        $has_kegg = 0;
                    }
                    # kegg id is part of a pathway map
                    if ($has_kegg &&
                        $annotation_hash{$row}{kegg_map} ne '')
                    {
                        $has_kegg = 2;
                    }
                    
                    if ($has_kegg == $best_has_kegg &&
                        $count    == $best_count &&
                        $score    == $best_score)
                    {
                        $match_flag = 1;
                        
                        # store each row only once per row
                        if (!defined($matched_row_hash{$row}))
                        {
                            $matched_row_array[$num_matches]  = $row;
                            $matched_type_array[$num_matches] = $match_type;
                            $num_matches++;
                        }
                        $matched_row_hash{$row} = 1;
                        $matched_row_type_hash{$row}{$match_type} = 1;
                    }
                    # no longer tied
                    else
                    {
                        last;
                    }
                }
            }

            # 3A/3B: conformed name, all rows
            if ($match_flag == 0)
            {
                %candidate_name_hash = ();
                %temp_row_score_hash = ();
                @candidate_row_array = ();

                # score each row
                foreach $row (@row_all_array)
                {
                    # sanity check the formulas
                    # skip rows that have wrong numbers of non-hydrogens
                    #
                    $fsanity_db = $annotation_hash{$row}{fsanity};
                    if ($fsanity ne '' && $fsanity_db ne '' &&
                        $fsanity ne $fsanity_db)
                    {
                        next;
                    }

                    $name_db           = $annotation_hash{$row}{name};
                    $name_db_lc        = lc $name_db;
                    $name_db_conformed = $conformed_name_hash{$name_db};
                
                    if ($name_lc_conformed eq $name_db_conformed)
                    {
                        $frac_id = 0;
                        $score =
                            score_substring_mismatch($name_lc,
                                                     $name_db_lc,
                                                     $align_method,
                                                     \$frac_id);

                        $candidate_name_hash{$name_db_lc} = 1;

                        $temp_row_score_hash{$row} = $score;
                    }
                }
                
                @candidate_row_array =
                    sort cmp_conformed_rows keys %temp_row_score_hash;
                @temp_name_array     = keys %candidate_name_hash;

                # only 1 good candidate in database, should be right?
                if (@temp_name_array == 1)
                {
                    $match_type = '3A';
                }
                # less certain, since there are multiple candidates
                else
                {
                    $match_type = '3B';
                }

                # keep highest score (and ties) as matches
                if (@candidate_row_array)
                {
                    $best_row      = $candidate_row_array[0];
                    $best_score    = $temp_row_score_hash{$best_row};
                    $best_count    = $annotation_hash{$best_row}{col_count};
                    $best_has_kegg = $annotation_hash{$best_row}{kegg};
                    if (defined($best_has_kegg) && $best_has_kegg ne '')
                    {
                        $best_has_kegg = 1;
                    }
                    else
                    {
                        $best_has_kegg = 0;
                    }
                    # kegg id is part of a pathway map
                    if ($best_has_kegg &&
                        $annotation_hash{$best_row}{kegg_map} ne '')
                    {
                        $best_has_kegg = 2;
                    }
                }
                foreach $row (@candidate_row_array)
                {
                    $score    = $temp_row_score_hash{$row};
                    $count    = $annotation_hash{$row}{col_count};
                    $has_kegg = $annotation_hash{$row}{kegg};
                    if (defined($has_kegg) && $has_kegg ne '')
                    {
                        $has_kegg = 1;
                    }
                    else
                    {
                        $has_kegg = 0;
                    }
                    # kegg id is part of a pathway map
                    if ($has_kegg &&
                        $annotation_hash{$row}{kegg_map} ne '')
                    {
                        $has_kegg = 2;
                    }
                    
                    if ($has_kegg == $best_has_kegg &&
                        $count    == $best_count &&
                        $score    == $best_score)
                    {
                        $match_flag = 1;
                        
                        # store each row only once per row
                        if (!defined($matched_row_hash{$row}))
                        {
                            $matched_row_array[$num_matches]  = $row;
                            $matched_type_array[$num_matches] = $match_type;
                            $num_matches++;

                            $ppm = 1000000 *
                                   abs(($annotation_hash{$row}{mz} - $mz)
                                        / $annotation_hash{$row}{mz});

                            if ($annotation_hash{$row}{mz} < 99999)
                            {
                                printf STDERR "WARNING -- mz differ:  %s  %f  %f  %.1f  %s  %s\n",
                                    $first_field,
                                    $mz,
                                    $annotation_hash{$row}{mz},
                                    $ppm,
                                    $name_delim,
                                    $annotation_hash{$row}{name},
                            }
                        }
                        $matched_row_hash{$row} = 1;
                        $matched_row_type_hash{$row}{$match_type} = 1;
                    }
                    # no longer tied
                    else
                    {
                        last;
                    }
                }
            }

            # 4A/4B: conformed, all rows, synonyms
            if ($match_flag == 0)
            {
                %candidate_name_hash = ();
                %temp_row_score_hash = ();
                @candidate_row_array = ();

                # score each row
                foreach $row (@row_all_array)
                {
                    # sanity check the formulas
                    # skip rows that have wrong numbers of non-hydrogens
                    #
                    $fsanity_db = $annotation_hash{$row}{fsanity};
                    if ($fsanity ne '' && $fsanity_db ne '' &&
                        $fsanity ne $fsanity_db)
                    {
                        next;
                    }

                    @synonym_array     = @{$annotation_hash{$row}{syn_arr}};
                    $name_db           = $annotation_hash{$row}{name};
                    $name_db_lc        = lc $name_db;
                    $name_db_conformed = $conformed_name_hash{$name_db};
                
                    foreach $synonym_db_conformed (@synonym_array)
                    {
                        if ($name_lc_conformed eq $synonym_db_conformed)
                        {
                            $frac_id = 0;
                            $score =
                                score_substring_mismatch($name_lc_conformed,
                                                         $name_db_conformed,
                                                         $align_method,
                                                         \$frac_id);

                            $candidate_name_hash{$name_db_lc} = 1;

                            $temp_row_score_hash{$row} = $score;
                            
                            last;
                        }
                    }
                }
                
                @candidate_row_array =
                    sort cmp_conformed_rows keys %temp_row_score_hash;
                @temp_name_array     = keys %candidate_name_hash;

                # only 1 good candidate in database, should be right?
                if (@temp_name_array == 1)
                {
                    $match_type = '4A';
                }
                # less certain, since there are multiple candidates
                else
                {
                    $match_type = '4B';
                }

                # keep highest score (and ties) as matches
                if (@candidate_row_array)
                {
                    $best_row      = $candidate_row_array[0];
                    $best_score    = $temp_row_score_hash{$best_row};
                    $best_count    = $annotation_hash{$best_row}{col_count};
                    $best_has_kegg = $annotation_hash{$best_row}{kegg};
                    if (defined($best_has_kegg) && $best_has_kegg ne '')
                    {
                        $best_has_kegg = 1;
                    }
                    else
                    {
                        $best_has_kegg = 0;
                    }
                    # kegg id is part of a pathway map
                    if ($best_has_kegg &&
                        $annotation_hash{$best_row}{kegg_map} ne '')
                    {
                        $best_has_kegg = 2;
                    }
                }
                foreach $row (@candidate_row_array)
                {
                    $score    = $temp_row_score_hash{$row};
                    $count    = $annotation_hash{$row}{col_count};
                    $has_kegg = $annotation_hash{$row}{kegg};
                    if (defined($has_kegg) && $has_kegg ne '')
                    {
                        $has_kegg = 1;
                    }
                    else
                    {
                        $has_kegg = 0;
                    }
                    # kegg id is part of a pathway map
                    if ($has_kegg &&
                        $annotation_hash{$row}{kegg_map} ne '')
                    {
                        $has_kegg = 2;
                    }
                    
                    if ($has_kegg == $best_has_kegg &&
                        $count    == $best_count &&
                        $score    == $best_score)
                    {
                        $match_flag = 1;
                        
                        # store each row only once per row
                        if (!defined($matched_row_hash{$row}))
                        {
                            $matched_row_array[$num_matches]  = $row;
                            $matched_type_array[$num_matches] = $match_type;
                            $num_matches++;

                            $ppm = 1000000 *
                                   abs(($annotation_hash{$row}{mz} - $mz)
                                        / $annotation_hash{$row}{mz});

                            if ($annotation_hash{$row}{mz} < 99999)
                            {
                                printf STDERR "WARNING -- mz differ:  %s  %f  %f  %.1f  %s  %s\n",
                                    $first_field,
                                    $mz,
                                    $annotation_hash{$row}{mz},
                                    $ppm,
                                    $name_delim,
                                    $annotation_hash{$row}{name},
                            }
                        }
                        $matched_row_hash{$row} = 1;
                        $matched_row_type_hash{$row}{$match_type} = 1;
                    }
                    # no longer tied
                    else
                    {
                        last;
                    }
                }
            }


            # 5A/5B: fuzzy, original m/z + pos/neg, name + synonyms
            if ($match_flag == 0)
            {
                %candidate_name_hash = ();
                %temp_row_score_hash = ();
                $best_score          = 0;
                #$name_lc_conformed   = conform_name($name_lc);

                # score each row
                foreach $row (@row_mz_pos_neg_array)
                {
                    $name_db           = $annotation_hash{$row}{name};
                    $name_db_lc        = lc $name_db;
                    $name_db_conformed = $conformed_name_hash{$name_db};

                    $frac_id = 0;
                    $score = score_substring_mismatch($name_lc_conformed,
                                                      $name_db_conformed,
                                                      $align_method,
                                                      \$frac_id);
                    #if ($frac_id < 0.5) { $score = 0; }
                    #if ($score < 0.5)   { $score = 0; }

                    if ($score > 0)
                    {
                        $candidate_name_hash{$name_db_lc} = 1;

                        if ($score >= $best_score)
                        {
                            $temp_row_score_hash{$row} = $score;
                            $best_score = $score;

                            if (0 && $name_lc =~ /lactose/ && $name_db_lc =~ /lactose/)
                            {
                                printf STDERR "LACTOSE  %s  %f  %s  %s  %s\n",
                                              $name_lc, $score,
                                              $name_lc_conformed,
                                              $name_db_conformed,
                                              $name_db_conformed;
                            }
                        }
                    }

                    if (1)
                    {
                      @synonym_array = @{$annotation_hash{$row}{syn_arr}};
                      foreach $synonym_db_conformed (@synonym_array)
                      {
                        $frac_id = 0;
                        $score =
                            score_substring_mismatch($name_lc_conformed,
                                                     $synonym_db_conformed,
                                                     $align_method,
                                                     \$frac_id);

                        if ($score > 0)
                        {
                            $candidate_name_hash{$name_db_lc} = 1;

                            if ($score >= $best_score)
                            {
                                $temp_row_score_hash{$row} = $score;
                                $best_score = $score;

                                if (0 && $name_lc =~ /lactose/ && $name_db_lc =~ /lactose/)
                                {
                                    printf STDERR "LACTOSE  %s  %f  %s  %s  %s\n",
                                                  $name_lc, $score,
                                                  $name_lc_conformed,
                                                  $name_db_conformed,
                                                  $synonym_db_conformed;
                                }
                            }
                        }
                      }
                    }
                }

                @candidate_row_array =
                    sort cmp_fuzzy_rows keys %temp_row_score_hash;
                @temp_name_array     = keys %candidate_name_hash;

                # only 1 good candidate in database, should be right?
                if (@temp_name_array == 1)
                {
                    $match_type = '5A';
                }
                # less certain, since there are multiple candidates
                else
                {
                    $match_type = '5B';
                }

                # keep highest score (and ties) as matches
                if (@candidate_row_array)
                {
                    $best_row      = $candidate_row_array[0];
                    $best_score    = $temp_row_score_hash{$best_row};
                    $best_count    = $annotation_hash{$best_row}{col_count};
                    $best_has_kegg = $annotation_hash{$best_row}{kegg};
                    if (defined($best_has_kegg) && $best_has_kegg ne '')
                    {
                        $best_has_kegg = 1;
                    }
                    else
                    {
                        $best_has_kegg = 0;
                    }
                    # kegg id is part of a pathway map
                    if ($best_has_kegg &&
                        $annotation_hash{$best_row}{kegg_map} ne '')
                    {
                        $best_has_kegg = 2;
                    }
                }
                foreach $row (@candidate_row_array)
                {
                    $score    = $temp_row_score_hash{$row};
                    $count    = $annotation_hash{$row}{col_count};
                    $has_kegg = $annotation_hash{$row}{kegg};
                    if (defined($has_kegg) && $has_kegg ne '')
                    {
                        $has_kegg = 1;
                    }
                    else
                    {
                        $has_kegg = 0;
                    }
                    # kegg id is part of a pathway map
                    if ($has_kegg &&
                        $annotation_hash{$row}{kegg_map} ne '')
                    {
                        $has_kegg = 2;
                    }
                    
                    if ($has_kegg == $best_has_kegg &&
                        $count    == $best_count &&
                        $score    == $best_score)
                    {
                        $match_flag = 1;
                        
                        # store each row only once per row
                        if (!defined($matched_row_hash{$row}))
                        {
                            $matched_row_array[$num_matches]  = $row;
                            $matched_type_array[$num_matches] = $match_type;
                            $num_matches++;
                        }
                        $matched_row_hash{$row} = 1;
                        $matched_row_type_hash{$row}{$match_type} = 1;
                    }
                    # no longer tied
                    else
                    {
                        last;
                    }
                }
            }

            # still no hits
            if ($match_flag == 0)
            {
                $match_type = '9X';

                # we don't currently have any lipidomics identifier mappings
                # so, for now, don't output any warnings
                if ($lipidomics_flag == 0)
                {
#                    printf STDERR "UNMATCHED   %s   %s   %s   %s\n",
#                                  $first_field, $mz, $name, $name_oc;
                }

                $matched_row_array[$num_matches]  = $bad_row_id;
                $matched_type_array[$num_matches] = $match_type;
                $num_matches++;

                $matched_row_hash{$bad_row_id} = 1;
                $matched_row_type_hash{$bad_row_id}{$match_type} = 1;

                #$bad_mz_name_hash{$name_oc} = 'unmatched';
                #$has_bad_mz_flag            = 1;
            }
            

            # retention time sanity checks
            if ($num_matches && is_number($rt_data))
            {
                %rt_ok_hash = ();

                foreach $row (@matched_row_array)
                {
                    $rt_db   = $annotation_hash{$row}{rt};

                    if (is_number($rt_db))
                    {
                        $delta = abs($rt_db - $rt_data);
                    
                        if ($delta <= $rt_tol)
                        {
                            $name_db = $annotation_hash{$row}{name};

                            $rt_ok_hash{$name_db} = 1;
                        }
                    }
                }
                
                foreach $row (@matched_row_array)
                {
                    $rt_db = $annotation_hash{$row}{rt};

                    if (is_number($rt_db))
                    {
                        $delta   = abs($rt_db - $rt_data);
                        $name_db = $annotation_hash{$row}{name};
                    
                        if ($delta > $rt_tol &&
                            !defined($rt_ok_hash{$name_db}))
                        {
                          if (defined($bad_mz_row_hash{$row}))
                          {
                            printf STDERR "WARNING -- rt differ:\t%s\t%f\t%f\t%s\t%s\t%f\t%f\n",
                                $first_field,
                                $mz,
                                $annotation_hash{$row}{mz},
                                $name_delim,
                                $annotation_hash{$row}{name},
                                $rt_data, $rt_db;
                          }
                        }
                    }
                }
            }
            

            # at least one good m/z match
            if ($num_matches && !($match_type =~ /[34]/))
            {
                $has_good_mz_flag = 1;
            }
            
            # printf STDERR "%s  %s\n", $name, @matched_type_array;
        }
    }
    

    # remove duplicate matches (usually from fuzzy match pos/neg offsets)
    if ($num_matches > 1)
    {
        @new_matched_row_array  = ();
        @new_matched_type_array = ();
        %seen_match_hash = ();
        $j = 0;
        for ($i = 0; $i < $num_matches; $i++)
        {
            $row     = $matched_row_array[$i];

            $name_db = '';
            $formula = '';
            $kegg    = '';
            $hmdb    = '';
            $pubchem = '';

            if (defined($annotation_hash{$row}))
            {
                $name_db = $annotation_hash{$row}{name};
                $formula = $annotation_hash{$row}{formula};
                $kegg    = $annotation_hash{$row}{kegg};
                $hmdb    = $annotation_hash{$row}{hmdb};
                $pubchem = $annotation_hash{$row}{pubchem};

                if (!defined($name_db)) { $name_db    = ''; }
                if (!defined($formula)) { $formula_db = ''; }
                if (!defined($kegg))    { $kegg_db    = ''; }
                if (!defined($hmdb))    { $hmdb_db    = ''; }
                if (!defined($pubchem)) { $pubchem_db = ''; }
            }


            $temp_id = sprintf "%s|%s|%s|%s|%s",
                $name_db, $formula, $kegg, $hmdb, $pubchem;


            # store only matches we haven't encountered yet
            if ($row == $bad_row_id || !defined($seen_match_hash{$temp_id}))
            {
                $new_matched_row_array[$j]  = $row;
                $new_matched_type_array[$j] = $matched_type_array[$i];

                $j++;
            }
            
            $seen_match_hash{$temp_id} = 1;
        }
        @matched_row_array  = @new_matched_row_array;
        @matched_type_array = @new_matched_type_array;
        $num_matches        = $j;
    }
    
    
    # set match type to worst observed per row
    if ($num_matches >= 1)
    {
        for ($i = 0; $i < $num_matches; $i++)
        {
            $row = $matched_row_array[$i];
            
            @temp_array = sort {$b cmp $a}
                          keys %{$matched_row_type_hash{$row}};

            if (@temp_array > 1)
            {
                $matched_type_array[$i] = $temp_array[0];
            }
        }
    }


    # remove originally bad MZMine identifications
    $name_corrected = $name;
    if ($has_good_mz_flag && $has_bad_mz_flag)
    {
        # remove hit from original name field
        @name_array_new = ();
        $j = 0;
        for ($i = 0; $i < @name_array; $i++)
        {
            $name_oc = $name_array[$i];
            
            if (!defined($bad_mz_name_hash{$name_oc}))
            {
                $name_array_new[$j++] = $name_oc;
            }
        }
        $name_corrected = join ';', @name_array_new;

        printf STDERR "CORRECTION -- remove poor ids:\t%s\t%s --> %s\n",
            $first_field, $name, $name_corrected;
        
        # convert clean | delimited text back to ; delimited
        $name_corrected =~ s/\|/\;/g;
        
        # overwrite original field with new field
        $array[$data_name_col] = $name_corrected;



        # remove matches to removed identification
        @new_matched_row_array  = ();
        @new_matched_type_array = ();
        $j = 0;
        for ($i = 0; $i < $num_matches; $i++)
        {
            $row = $matched_row_array[$i];

            # store only good matches
            if (!defined($bad_mz_row_hash{$row}))
            {
                $new_matched_row_array[$j]  = $row;
                $new_matched_type_array[$j] = $matched_type_array[$i];

                $j++;
            }
        }
        @matched_row_array  = @new_matched_row_array;
        @matched_type_array = @new_matched_type_array;
        $num_matches        = $j;
    }
    
    
    # output new line
    # insert new annotation immediately after identity/name column
    $tab_flag = 0;
    for ($col = 0; $col < @array; $col++)
    {
        $value = $array[$col];
        
        if ($tab_flag)
        {
            print "\t";
        }
        
        print "$value";
        

        # insert new annotation columns
        if ($col == $data_name_col)
        {
            $name_db_str    = '';
            $match_type_str = '';
            $formula_str    = '';
            $kegg_str       = '';
            $hmdb_str       = '';
            $pubchem_str    = '';

            #concatenate multiple annotation entries
            for ($i = 0; $i < @matched_row_array; $i++)
            {
                $row = $matched_row_array[$i];

                $name_db    = $annotation_hash{$row}{name};
                $match_type = $matched_type_array[$i];
                $formula    = $annotation_hash{$row}{formula};
                $kegg       = $annotation_hash{$row}{kegg};
                $hmdb       = $annotation_hash{$row}{hmdb};
                $pubchem    = $annotation_hash{$row}{pubchem};
                
                if (!defined($formula)) { $formula_db = ''; }
                if (!defined($kegg))    { $kegg_db    = ''; }
                if (!defined($hmdb))    { $hmdb_db    = ''; }
                if (!defined($pubchem)) { $pubchem_db = ''; }
                
                $formula_new = conform_formula($formula);
                if ($formula_new ne $formula)
                {
                    printf STDERR "Reordering formula %s --> %s\n",
                        $formula, $formula_new;
                }
                $formula = $formula_new;
            
                # any matches after the first one
                if ($i)
                {
                    $name_db_str    .= '|';
                    $match_type_str .= '|';
                    $formula_str    .= '|';
                    $kegg_str       .= '|';
                    $hmdb_str       .= '|';
                    $pubchem_str    .= '|';
                }

                $name_db_str    .= $name_db;
                $match_type_str .= $match_type;
                $formula_str    .= $formula;
                $kegg_str       .= $kegg;
                $hmdb_str       .= $hmdb;
                $pubchem_str    .= $pubchem;
            }
            
            # blank multiple match fields without any mappings
            if (!($name_db_str    =~ /[^|]/)) { $name_db_str    = ''; }
            if (!($match_type_str =~ /[^|]/)) { $match_type_str = ''; }
            if (!($formula_str    =~ /[^|]/)) { $formula_str    = ''; }
            if (!($kegg_str       =~ /[^|]/)) { $kegg_str       = ''; }
            if (!($hmdb_str       =~ /[^|]/)) { $hmdb_str       = ''; }
            if (!($pubchem_str    =~ /[^|]/)) { $pubchem_str    = ''; }


            # no matches found for this row
            # report match type of 9X
            if ($name =~ /\S/ && $match_type_str eq '')
            {
                $match_type_str = '9X';
            }


            printf "\t%s", $name_db_str;
            printf "\t%s", $match_type_str;
            printf "\t%s", $formula_str;
            printf "\t%s", $kegg_str;
            printf "\t%s", $hmdb_str;
            printf "\t%s", $pubchem_str;
        }
        
        $tab_flag = 1;
    }
    printf "\n";
    
    $row_data++;
}
