#!/usr/bin/perl -w


# TODO -- deal with multiple matches to the same metabolite per row,
#         set match type to worst of the multiple matches for QC purposes


# set lib search path to directory the script is run from
use File::Basename;
use lib dirname (__FILE__);

use Scalar::Util qw(looks_like_number);
use POSIX;
use align_text;    # text string alignment module

$mz_tol_ppm = 10;    # 10 ppm
$rt_tol     = 1.0;   # minutes

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
        $temp_array[$i] =~ tr/\//\|/;
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

# L-aminoacid     : aminoacid        ex. L-Alanine     --> Alanine
# "anoic"         : "yric"           ex. Butanoic acid --> Butyric acid
# "anoate\b"      : "yrate"          ex. Butanaoate    --> Butyrate
# "ic acid\b"     : "ate"            ex. glutamic acid --> glutamate
# "monosomething" : "something"      ex. monophosphate --> phosphate

# strip ()
# strip spaces
# maybe strip -

sub preprocess_name
{
    my $name       = $_[0];
    my $name_len;
    my $half_name_len;
    my $name_half1 = '';
    my $name_half2 = '';
    
    # lowercase everything
    $name = lc $name;

    # replace Greek letters at word boundaries
    $name =~ s/\balpha\b/a/g;
    $name =~ s/\bbeta\b/b/g;
    $name =~ s/\bgamma\b/g/g;
    $name =~ s/\bdelta\b/d/g;
    $name =~ s/\bepsilon\b/e/g;
    $name =~ s/\bzeta\b/z/g;
    $name =~ s/\beta\b/h/g;
    
    # strip the L- from L-aminoacids
    # only on word boundary, so we don't strip DL-aminoacid
    $name =~ s/\bl-(alanine)/$1/g;
    $name =~ s/\bl-(arginine)/$1/g;
    $name =~ s/\bl-(asparagine)/$1/g;
    $name =~ s/\bl-(aspartic acid)/$1/g;
    $name =~ s/\bl-(cysteine)/$1/g;
    $name =~ s/\bl-(glutamine)/$1/g;
    $name =~ s/\bl-(glutamic acid)/$1/g;
    # $name =~ s/\bl-(glycine)/$1/g;           # L- is never used !!
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
    
    # conform acids
    $name =~ s/anoic/yric/g;          # Butanoic acid --> Butyric acid
    $name =~ s/anoate\b/yrate/g;      # Butanoate     --> Butyrate
    $name =~ s/ic acid\b/ate/g;       # Glutamic acid --> Glutamate
    
    # mono is redundant and often left out in synonyms
    $name =~ s/mono(\S)/$1/g;
    
    # condense everything that isn't a letter, number, comma, or plus/minus
    # except when between two numbers
    #
    # protect -) as in (+/-) or (-) using capital letters
    $name =~ s/\(-|-\)/MINUS/g;       # protect minus signs
    $name =~ s/[^A-Za-z0-9,+]/-/g;    # convert to hyphens
    $name =~ s/-+/-/g;                # condense multiple hyphens in a row
    $name =~ s/(^-|-$)//g;            # strip leading/trailing hyphens
    $name =~ s/([,+])-/$1/g;          # strip hyphens next to comma or plus
    $name =~ s/-([,+])/$1/g;          # strip hyphens next to comma or plus
    $name =~ s/([a-z])-/$1/g;         # strip hyphens next to letters
    $name =~ s/-([a-z])/$1/g;         # strip hyphens near to letters
    $name =~ s/,+/,/g;                # condense multiple commas in a row
    $name =~ s/\++/\+/g;              # condense multiple pluses in a row
    
    # de-protect and condense minus signs
    #
    # hypothetical example: 1-(-)-galacturonate   -->   1-galacturonate
    #                       1-(+)-galacturonate   -->   1+galacturonate
    #                       1-galacturonate       -->   1galacturonate
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
        printf STDERR "ABORT -- unknown option %s\n", $field;
        $syntax_error_flag = 1;
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

    # $field = $array[$i];
    # $annotation_header_col_hash{$field} = $i;
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
    elsif (!defined($annotation_kegg_col) &&
           $header =~ /^kegg/i)
    {
        $annotation_kegg_col = $i;
    }
    elsif (!defined($annotation_hmdb_col) &&
           $header =~ /^hmdb/i)
    {
        $annotation_hmdb_col = $i;
    }
    elsif (!defined($annotation_pubchem_col) &&
           $header =~ /^pubchem/i)
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
if (!defined($annotation_kegg_col))
{
    printf STDERR "ABORT -- KEGG ID column not found in annotation file %s\n",
                 $annotation_filename;
    exit(1);
}
if (!defined($annotation_hmdb_col))
{
    printf STDERR "ABORT -- HMDB column not found in annotation file %s\n",
                 $annotation_filename;
    exit(1);
}
if (!defined($annotation_pubchem_col))
{
    printf STDERR "ABORT -- PubChem column not found in annotation file %s\n",
                 $annotation_filename;
    exit(1);
}



# read in the annotation
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

    $mz      = $array[$annotation_mz_col];
    $formula = $array[$annotation_formula_col];
    $name    = $array[$annotation_name_col];
    $kegg    = $array[$annotation_kegg_col];
    $hmdb    = $array[$annotation_hmdb_col];
    $pubchem = $array[$annotation_pubchem_col];
    
    # retention time sanity checks, if available
    if (defined($annotation_rt_col))
    {
        $rt  = $array[$annotation_rt_col];
    }
    else
    {
        $rt  = '';
    }
    
    if (!defined($mz))      { $mz      = ''; }
    if (!defined($formula)) { $formula = ''; }
    if (!defined($name))    { $name    = ''; }
    if (!defined($kegg))    { $kegg    = ''; }
    if (!defined($hmdb))    { $hmdb    = ''; }
    if (!defined($pubchem)) { $pubchem = ''; }
    if (!defined($pubchem)) { $rt      = ''; }
    
    # skip bad rows
    if (!is_number($mz))  { next; }
    if (!($name =~ /\S/)) { next; }
    
    # keep one row per combination of identifiers, ignoring pos/neg, rt, etc.
    #
    # no longer necessary, since hits are now de-duped during printing
    if (0)
    {
        # unique identifier for mapping
        $unique_id  =       $mz;
        $unique_id .= '~' . $formula;
        $unique_id .= '~' . $name;
        $unique_id .= '~' . $kegg;
        $unique_id .= '~' . $hmdb;
        $unique_id .= '~' . $pubchem;
        
        # only keep one copy of each mapping
        # "duplicates" can happen due to peaks at multiple retention times
        if (defined($unique_id_hash{$unique_id}))
        {
            next;
        }
        $unique_id_hash{$unique_id} = 1;
    }

    $annotation_hash{$row}{mz}      = $mz;
    $annotation_hash{$row}{formula} = $formula;
    $annotation_hash{$row}{name}    = $name;
    $annotation_hash{$row}{kegg}    = $kegg;
    $annotation_hash{$row}{hmdb}    = $hmdb;
    $annotation_hash{$row}{pubchem} = $pubchem;
    $annotation_hash{$row}{rt}      = $rt;
    
    # bin m/z for rapid scanning later
    # largest m/z we see is 900
    # +/-1 is ~1000 ppm, which is *way* over our 10 ppm tolerance
    # bins are too large, but are very simple to code up
    #
    # also factor in +/- 2.014552 (pos - neg) difference
    #
    $mz_floor = floor($mz);
    $mz_minus3 = $mz_floor - 3;
    $mz_minus2 = $mz_floor - 2;
    $mz_minus1 = $mz_floor - 1;
    $mz_plus1  = $mz_floor + 1;
    $mz_plus2  = $mz_floor + 2;
    $mz_plus3  = $mz_floor + 3;

    $mz_row_bins_hash{$mz_minus1}{$row} = 1;
    $mz_row_bins_hash{$mz_floor}{$row}  = 1;
    $mz_row_bins_hash{$mz_plus1}{$row}  = 1;

    if (defined($annotation_pos_neg_col))
    {
        $pos_neg = $array[$annotation_pos_neg_col];
        
        $pos_flag = 0;
        $neg_flag = 0;
        
        if ($pos_neg =~ /pos/i) { $pos_flag = 1; }
        if ($pos_neg =~ /neg/i) { $neg_flag = 1; }
        
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
    # include all +/- 2.014552, without knowing which are pos/neg
    else
    {
        $mz_row_bins_hash{$mz_minus3}{$row} = 1;
        $mz_row_bins_hash{$mz_minus2}{$row} = 1;
        $mz_row_bins_hash{$mz_plus2}{$row}  = 1;
        $mz_row_bins_hash{$mz_plus3}{$row}  = 1;
    }

    $mz_row_bins_hash{$mz_minus3}{$row} = 1;
    $mz_row_bins_hash{$mz_minus2}{$row} = 1;
    $mz_row_bins_hash{$mz_plus2}{$row}  = 1;
    $mz_row_bins_hash{$mz_plus3}{$row}  = 1;

    $row++;
}
#$row_count = $row;


# all rows, to be used when looking for matches regardless of m/z
@row_all_array = sort {$a <=> $b} keys %annotation_hash;


# read in data file

# skip down to first line that has anything on it
# lipidomics data has this issues sometimes
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
        printf "\t%s", 'Formula';
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


    # matched (mapped) rows
    %matched_row_hash  = ();
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
        %candidate_row_mz_hash      = ();
        %candidate_row_pos_neg_hash = ();

        $mz_floor = floor($mz);
        
        # scan annotations for m/z within tolerance
        #
        # also include +/- 2.014552 (pos - neg difference)
        #  already pre-binned earlier during +/- 3 binning
        @row_array = sort {$a<=>$b} keys %{$mz_row_bins_hash{$mz_floor}};
        foreach $row (@row_array)
        {
            $mz_db = $annotation_hash{$row}{mz};
            
            $ppm    = 1000000 * abs(( $mz_db             - $mz) / $mz_db);
            $ppm_m2 = 1000000 * abs((($mz_db - 2.014552) - $mz) /
                                    ( $mz_db - 2.014552));
            $ppm_p2 = 1000000 * abs((($mz_db + 2.014552) - $mz) /
                                    ( $mz_db + 2.014552));

            # check only given m/z
            if ($ppm <= $mz_tol_ppm)
            {
                $candidate_row_mz_hash{$row}      = 1;
                $candidate_row_pos_neg_hash{$row} = 1;
            }
            
            # also check +/- 2.014552, in case we only have one of pos/neg
            if ($ppm_m2 <= $mz_tol_ppm || $ppm_p2 <= $mz_tol_ppm)
            {
                $candidate_row_pos_neg_hash{$row} = 1;
            }
        }
        
        # original and +/- ~2 m/z row arrays
        @row_mz_array      = sort {$a<=>$b} keys %candidate_row_mz_hash;
        @row_pos_neg_array = sort {$a<=>$b} keys %candidate_row_pos_neg_hash;
        
        # support multiple ;-delimited names per data row
        $name_delim = $name;
        $name_delim = bless_delimiter_bar_metabolomics($name_delim);
        @name_array = split /\|/, $name_delim;

        # checks in less confident order:
        #
        #   1A: exact match, original m/z
        #   1B: conformed,   original m/z
        #   2A: exact match, pos/neg  m/z
        #   2B: conformed,   pos/neg  m/z
        #   3A: exact_match, all rows, regardless of m/z
        #   3B: conformed,   all rows, regardless of m/z
        #   4A: fuzzy,       original m/z, one candidate
        #   4B: fuzzy,       original m/z, >= 2 candidates
        #   5A: fuzzy,       pos/neg  m/z, one candidate
        #   5B: fuzzy,       pos/neg  m/z, >= 2 candidates
        #   9X: unmatched
        
        foreach $name_oc (@name_array)
        {
            $name_lc    = lc $name_oc;
        
            $match_flag = 0;
            $match_type = '9x';    # 9: unmatched
        
            # 1A: exact match, original m/z
            foreach $row (@row_mz_array)
            {
                $name_db_lc = lc $annotation_hash{$row}{name};
            
                if ($name_lc eq $name_db_lc)
                {
                    $match_flag = 1;
                    $match_type = '1A';
                
                    # store each row only once per row
                    if (!defined($matched_row_hash{$row}))
                    {
                        $matched_row_array[$num_matches]  = $row;
                        $matched_type_array[$num_matches] = $match_type;
                        $num_matches++;
                    }
                    $matched_row_hash{$row} = 1;
                }
            }

            # 1B: conformed,   original m/z
            if ($match_flag == 0)
            {
                $name_lc_conformed = preprocess_name($name_lc);
            
                foreach $row (@row_mz_array)
                {
                    $name_db_conformed =
                      preprocess_name($annotation_hash{$row}{name});
                
                    if ($name_lc_conformed eq $name_db_conformed)
                    {
                        $match_flag = 1;
                        $match_type = '1B';
                    
                        # store each row only once per row
                        if (!defined($matched_row_hash{$row}))
                        {
                            $matched_row_array[$num_matches]  = $row;
                            $matched_type_array[$num_matches] = $match_type;
                            $num_matches++;
                        }
                        $matched_row_hash{$row} = 1;
                    }
                }
            }

            # 2A: exact match, pos/neg  m/z
            if ($match_flag == 0)
            {
                foreach $row (@row_pos_neg_array)
                {
                    $name_db_lc = lc $annotation_hash{$row}{name};
            
                    if ($name_lc eq $name_db_lc)
                    {
                        $match_flag = 1;
                        $match_type = '2A';
                
                        # store each row only once per row
                        if (!defined($matched_row_hash{$row}))
                        {
                            $matched_row_array[$num_matches]  = $row;
                            $matched_type_array[$num_matches] = $match_type;
                            $num_matches++;
                        }
                        $matched_row_hash{$row} = 1;
                    }
                }
            }

            # 2B: conformed,   pos/neg  m/z
            if ($match_flag == 0)
            {
                $name_lc_conformed = preprocess_name($name_lc);
            
                foreach $row (@row_pos_neg_array)
                {
                    $name_db_conformed =
                      preprocess_name($annotation_hash{$row}{name});
                
                    if ($name_lc_conformed eq $name_db_conformed)
                    {
                        $match_flag = 1;
                        $match_type = '2B';
                    
                        # store each row only once per row
                        if (!defined($matched_row_hash{$row}))
                        {
                            $matched_row_array[$num_matches]  = $row;
                            $matched_type_array[$num_matches] = $match_type;
                            $num_matches++;
                        }
                        $matched_row_hash{$row} = 1;
                    }
                }
            }

            # 3A: exact_match, all rows, regardless of m/z
            if ($match_flag == 0)
            {
                foreach $row (@row_all_array)
                {
                    $name_db_lc = lc $annotation_hash{$row}{name};
            
                    if ($name_lc eq $name_db_lc)
                    {
                        $match_flag = 1;
                        $match_type = '3A';
                        
                        $bad_mz_name_hash{$name_oc} = $row;
                        $bad_mz_row_hash{$row}      = $name_oc;
                        $has_bad_mz_flag            = 1;
                
                        # store each row only once per row
                        if (!defined($matched_row_hash{$row}))
                        {
                            $matched_row_array[$num_matches]  = $row;
                            $matched_type_array[$num_matches] = $match_type;
                            $num_matches++;

                            printf STDERR "WARNING -- mz differ:\t%s\t%f\t%f\t%s\t%s\n",
                                $first_field,
                                $mz,
                                $annotation_hash{$row}{mz},
                                $name_delim,
                                $annotation_hash{$row}{name},
                        }
                        $matched_row_hash{$row} = 1;
                    }
                }
            }

            # 3B: conformed,   all rows, regardless of m/z
            if ($match_flag == 0)
            {
                $name_lc_conformed = preprocess_name($name_lc);
            
                foreach $row (@row_all_array)
                {
                    $name_db_conformed =
                      preprocess_name($annotation_hash{$row}{name});
                
                    if ($name_lc_conformed eq $name_db_conformed)
                    {
                        $match_flag = 1;
                        $match_type = '3B';

                        $bad_mz_name_hash{$name_oc} = $row;
                        $bad_mz_row_hash{$row}      = $name_oc;
                        $has_bad_mz_flag            = 1;
                    
                        # store each row only once per row
                        if (!defined($matched_row_hash{$row}))
                        {
                            $matched_row_array[$num_matches]  = $row;
                            $matched_type_array[$num_matches] = $match_type;
                            $num_matches++;

                            printf STDERR "WARNING -- mz differ:\t%s\t%f\t%f\t%s\t%s\n",
                                $first_field,
                                $mz,
                                $annotation_hash{$row}{mz},
                                $name_delim,
                                $annotation_hash{$row}{name},
                        }
                        $matched_row_hash{$row} = 1;
                    }
                }
            }

            # 4A/4B: fuzzy,       original m/z
            if ($match_flag == 0)
            {
                %candidate_name_hash = ();
                %temp_row_score_hash = ();
                $best_score          = 0;
                $name_lc_conformed   = preprocess_name($name_lc);

                # score each row
                foreach $row (@row_mz_array)
                {
                    $name_db_conformed =
                        preprocess_name($annotation_hash{$row}{name});

                    $frac_id = 0;
                    $score = score_substring_mismatch($name_lc_conformed,
                                                      $name_db_conformed,
                                                      'glocal',
                                                      \$frac_id);
                    if ($frac_id < 0.5) { $score = 0; }
                    $temp_row_score_hash{$row} = $score;
                    
                    if ($score > $best_score)
                    {
                        $best_score = $score;
                    }
                }

                # build list of non- zero-scoring candidate names
                foreach $row (@row_mz_array)
                {
                    $score = $temp_row_score_hash{$row};
                    
                    if ($score > 0)
                    {
                        $name_db_conformed =
                            preprocess_name($annotation_hash{$row}{name});

                        $candidate_name_hash{$name_db_conformed} = 1;
                    }
                }
                @temp_name_array = keys %candidate_name_hash;

                # keep highest score (and ties) as matches
                foreach $row (@row_mz_array)
                {
                    $score = $temp_row_score_hash{$row};
                    
                    if ($score > 0 && $score == $best_score)
                    {
                        $match_flag = 1;
                        
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
                    
                        # store each row only once per row
                        if (!defined($matched_row_hash{$row}))
                        {
                            $matched_row_array[$num_matches]  = $row;
                            $matched_type_array[$num_matches] = $match_type;
                            $num_matches++;
                        }
                        $matched_row_hash{$row} = 1;

                        if (0)
                        {
                            $mz_db = $annotation_hash{$row}{mz};
                            $name_db_conformed =
                                preprocess_name($annotation_hash{$row}{name});

                            printf STDERR "%s   %s   %s   %s   %s   %0.4f   %s|%s\n",
                                $first_field,
                                $row_data,
                                $match_type,
                                $mz, $mz_db,
                                $score,
                                $name_lc_conformed, $name_db_conformed;
                        }
                    }
                }
            }

            # 5A/5B: fuzzy,       pos/neg  m/z
            # currently too prone to mis-mappings, disable it
            if (0 && $match_flag == 0)
            {
                %candidate_name_hash = ();
                %temp_row_score_hash = ();
                $best_score          = 0;
                $name_lc_conformed   = preprocess_name($name_lc);

                # score each row
                foreach $row (@row_pos_neg_array)
                {
                    $name_db_conformed =
                        preprocess_name($annotation_hash{$row}{name});

                    $frac_id = 0;
                    $score = score_substring_mismatch($name_lc_conformed,
                                                      $name_db_conformed,
                                                      'glocal',
                                                      \$frac_id);
                    if ($frac_id < 0.5) { $score = 0; }
                    $temp_row_score_hash{$row} = $score;
                    
                    if ($score > $best_score)
                    {
                        $best_score = $score;
                    }
                }

                # build list of non- zero-scoring candidate names
                foreach $row (@row_pos_neg_array)
                {
                    $score = $temp_row_score_hash{$row};
                    
                    if ($score > 0)
                    {
                        $name_db_conformed =
                            preprocess_name($annotation_hash{$row}{name});

                        $candidate_name_hash{$name_db_conformed} = 1;
                    }
                }
                @temp_name_array = keys %candidate_name_hash;

                # keep highest score (and ties) as matches
                foreach $row (@row_pos_neg_array)
                {
                    $score = $temp_row_score_hash{$row};
                    
                    if ($score > 0 && $score == $best_score)
                    {
                        $match_flag = 1;
                        
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
                    
                        # store each row only once per row
                        if (!defined($matched_row_hash{$row}))
                        {
                            $matched_row_array[$num_matches]  = $row;
                            $matched_type_array[$num_matches] = $match_type;
                            $num_matches++;
                        }
                        $matched_row_hash{$row} = 1;

                        if (0)
                        {
                            $mz_db = $annotation_hash{$row}{mz};
                            $name_db_conformed =
                                preprocess_name($annotation_hash{$row}{name});

                            printf STDERR "%s   %s   %s   %s   %s   %0.4f   %s|%s\n",
                                $first_field,
                                $row_data,
                                $match_type,
                                $mz, $mz_db,
                                $score,
                                $name_lc_conformed, $name_db_conformed;
                        }
                    }
                }
            }

            # still no hits
            if ($match_flag == 0)
            {
                $match_type = '9X';
                printf STDERR "UNMATCHED   %s   %s   %s   %s\n",
                              $first_field, $mz, $name, $name_oc;

                $bad_mz_name_hash{$name_oc} = 'unmatched';
                $has_bad_mz_flag            = 1;
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
            if ($num_matches && !($match_type =~ /3/))
            {
                $has_good_mz_flag = 1;
            }
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
            $row        = $matched_row_array[$i];

            $name_db    = $annotation_hash{$row}{name};
            $formula    = $annotation_hash{$row}{formula};
            $kegg       = $annotation_hash{$row}{kegg};
            $hmdb       = $annotation_hash{$row}{hmdb};
            $pubchem    = $annotation_hash{$row}{pubchem};

            if (!defined($formula)) { $formula_db = ''; }
            if (!defined($kegg))    { $kegg_db    = ''; }
            if (!defined($hmdb))    { $hmdb_db    = ''; }
            if (!defined($pubchem)) { $pubchem_db = ''; }
            
            $temp_id = sprintf "%s|%s|%s|%s|%s",
                $name_db, $formula, $kegg, $hmdb, $pubchem;

            # store only matches we haven't encountered yet
            if (!defined($seen_match_hash{$temp_id}))
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
            if ($name_db_str    eq '|') { $name_db_str    = ''; }
            if ($match_type_str eq '|') { $match_type_str = ''; }
            if ($formula_str    eq '|') { $formula_str    = ''; }
            if ($kegg_str       eq '|') { $kegg_str       = ''; }
            if ($hmdb_str       eq '|') { $hmdb_str       = ''; }
            if ($pubchem_str    eq '|') { $pubchem_str    = ''; }


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
