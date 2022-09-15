#!/usr/bin/perl -w

# 2022-09-15:  add ParentFormula column for LipidSearch adducts
# 2022-09-14:  conform formulas containing spaces (LipidSearch)
# 2022-08-11:  improve sample blank detection
# 2022-06-07:  less stringent heavy label string matching
# 2022-05-09:  added BOC to heavy labeled detection
# 2022-02-23:  extend sample renaming to non- height/area columns
# 2022-02-23:  use #[old_name]:new_name.raw headers to rename samples
# 2021-11-30:  better lipidomics support
# 2021-11-30:  clean up lipidomics header code
# 2021-11-29:  support Area[sample name] Height[sample name] format
# 2021-08-19:  change default back to leaving heavy unscaled
# 2021-08-19:  --scale-heavy --no-scale-heavy to --heavy-tracer --heavy-spikein
# 2021-08-18:  add new --scale-heavy and --no-scale-heavy flags
# 2021-08-11:  improve sample blank detection
# 2021-08-11:  change low signal value warning messages
# 2021-08-10:  more lipidomics molecule name header aliases
# 2021-08-09:  fix sample name terminal period, underscore, space typos
# 2021-08-09:  automatically fix ..mzXML (double-dot) typos
# 2021-07-26:  issue separate floor warnings for blank and non-blank samples
# 2021-07-23:  raise noise floor to 1000, report correct # floored peaks
# 2021-07-23:  raise noise floor to 50
# 2021-07-21:  zero out any data with a value < 10 as bad data
# 2021-07-14:  conform with buest guess when no standard height/peak found
# 2021-07-12:  more lipidomics support
# 2021-07-02:  issue a warning for "main ID" files, only one hit reported
# 2021-06-22:  fallback to 'row identity (main ID)' if (all IDs) not found
# 2021-06-16:  improved is_heavy_labeled() function
# 2021-05-28:  support non-abundance columns inserted at end of file,
#              output right-aligned contiguous block of sample data
# 2021-04-14:  begin adding support for lipidomics
# 2021-04-05:  El-MAVEN groupID is not unique, create missing identifiers
# 2021-04-05:  add heavy label detection to El-MAVEN
# 2021-01-06:  add --discard-unidentified and --discard-heavy flags
#              change some output column header names and capitalization
# 2020-09-24:  add more support for El-MAVEN
# 2020-09-10:  add check for Excel corruption of values >= 1E7
# 2020-09-08:  rename Spikein Flag header to Potential Spikein Flag
# 2020-08-26:  default to deleting "one-hit wonders" (single pre gap-fill peak)
# 2020-08-26:  attempt to deal with more pos/neg sample name typos
# 2020-08-26:  add Number of non-zero peaks, Percent pre-gap peaks,
#              Percent non-zero peaks columns
# 2020-07-24:  add Non-Spikein Identified Flag column to output
#              begin adding support for El-MAVEN
# 2020-07-13:  add additional command line options for output file names
#              more header conforming (strip .cdf extension)
#              change "assigned" to "identified"
#              no longer output assigned/identified file
#
# 2020-05-19: fix more header conforming
#
# 2020-03-20: detect and insert missing row ID column
#             case insensitive search now for .mzXML, due to case issues
#             strip .mzXML from sample names, as they aren't always there

use Scalar::Util qw(looks_like_number);
use POSIX;

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
# I don't currently check to see if the given 1- or 2- letter elements are
# valid known elements or not.  I could, but that would require a good bit
# more work than I have time for at the moment.  I'm not *quite* that paranoid
# about the formulas just yet, although part of me still worries about it...
#
sub conform_formula
{
    my $formula_orig = $_[0];
    my $adduct       = $_[1];
    my $formula;
    my $formula_new = '';
    my @match_array;
    my @element_array;
    my @ion_array;    # each +/- group of elements in the adduct
    my $ion;
    my $match;
    my $heavy;
    my $element;
    my $heavy_plus_element;
    my $count;
    my %count_hash = ();
    my %count_adduct_hash = ();
    my $has_carbon_flag = 0;
    my $i;
    
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
    
    
    # deal with adduct
    if (defined($adduct))
    {
        @ion_array = split /([+-])/, $adduct;
        
        for ($i = 1; $i < @ion_array - 1; $i += 2)
        {
            $add_sub = $ion_array[$i];
            $ion     = $ion_array[$i+1];

            @match_array = $ion =~ m/(?:\[[0-9]+\])*[A-Z][a-z]*(?:[0-9]+)*/g;
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

                # (-) adduct loses a C, so add back into to parent
                if ($add_sub eq '-' && $element eq 'C')
                {
                    $has_carbon_flag = 1;
                }

                $heavy_plus_element = $heavy . $element;

                if (length $element > 2)
                {
                    printf STDERR "WARNING -- error in adduct %s %s\n",
                        $adduct, $ion;

                    return $formula_orig;
                }
                
                if ($add_sub eq '-')
                {
                    $count_hash{$heavy_plus_element} += $count;
                }
                elsif ($add_sub eq '+')
                {
                    $count_hash{$heavy_plus_element} -= $count;
                }
            }
        }
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
                
                # skip elements we've subtracted away with the adduct
                if ($count <= 0)
                {
                    next;
                }

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

                # skip elements we've subtracted away with the adduct
                if ($count <= 0)
                {
                    next;
                }

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

        # skip elements we've subtracted away with the adduct
        if ($count <= 0)
        {
            next;
        }
        
        # replace 1 count with blank
        if ($count == 1)
        {
            $count = '';
        }
        
        $formula_new .= $heavy_plus_element . $count;
    }
    
    return $formula_new;
}


$keep_single_pregap_flag   = 0;
$discard_unidentified_flag = 0;
$discard_heavy_flag        = 0;
$scale_heavy_flag          = 0;
$seen_heavy_flag           = 0;

$syntax_error_flag         = 0;
$num_files                 = 0;

$floor_cutoff              = 1000;    # absurdly low abundance, floor to zero

for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        if ($field =~ /^--keep-single-pregap$/)
        {
            $keep_single_pregap_flag = 1;
        }
        elsif ($field =~ /^--discard-single-pregap$/)
        {
            $keep_single_pregap_flag = 0;
        }
        elsif ($field =~ /^--discard-unidentified$/)
        {
            $discard_unidentified_flag = 1;
        }
        elsif ($field =~ /^--discard-heavy$/)
        {
            $discard_heavy_flag = 1;
        }
        elsif ($field =~ /^--heavy-tracer$/)
        {
            $scale_heavy_flag = 1;
        }
        elsif ($field =~ /^--heavy-spikein$/)
        {
            $scale_heavy_flag = 0;
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
            $filename = $field;
            $num_files++;
        }
        elsif ($num_files == 1)
        {
            $output_exclusions_filename = $field;
            $num_files++;
        }
        elsif ($num_files == 2)
        {
            $output_spikeins_filename = $field;
            $num_files++;
        }
    }
}


if (!defined($filename) || $syntax_error_flag)
{
    printf STDERR "Usage: strip_mzmine_columns.pl [options] mzmine_tab_delimited.txt [unidentified.txt spikeins.txt]\n";
    printf STDERR "\n";
    printf STDERR "Options:\n";
    printf STDERR "    --heavy-spikein            treat heavy rows as spike-ins\n";
    printf STDERR "    --heavy-tracer             treat heavy rows as biological\n";
    printf STDERR "    --ppm N                    override default m/z PPM tolerance\n";
    printf STDERR "\n";
    printf STDERR "    --discard-heavy            discard heavy labeled rows\n";
    printf STDERR "    --discard-unidentified     discard unidentified rows\n";
    printf STDERR "\n";
    printf STDERR "  options which use the MZmine \"row number of detected peaks\" column:\n";
    printf STDERR "    --discard-single-pregap    discard pre gap-filled single-hit rows (default)\n";
    printf STDERR "    --keep-single-pregap       keep pre gap-filled single-hit rows\n";
    exit(1);
}


if (!defined($output_exclusions_filename))
{
    $output_exclusions_filename = 'metabolomics_auto_unidentified.txt';
}
if (!defined($output_spikeins_filename))
{
    $output_spikeins_filename     = 'metabolomics_auto_spikeins.txt';
}


open INFILE, "$filename" or die "can't open $filename\n";


# skip down to first line that has anything on it
# lipidomics data has this issues sometimes
while($line=<INFILE>)
{
    # skip comment lines
    # comment may contain sample renaming information, so check for it
    if ($line =~ /^#/)
    {
        if ($line =~ /^#\[([^\]]+)\]:(.*?)\.raw\s*$/)
        {
            $sample_rename_hash{$1} = $2;
            printf STDERR "Renaming:\t%s\t%s\n",
                $1, $2;
        }
    
        next;
    }

    # this line isn't purely whitespace, assume it is the header line
    if ($line =~ /\S/)
    {
        last;
    }
}


# header line
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
    $header_col_hash{$field} = $i;
    $header_col_array[$i] = $field;
}


# Excel gets very upset if first field is "ID", change it
if ($header_col_array[0] =~ /^id$/i)
{
    $header_col_array[0] = 'Index';
    $header_col_hash{'Index'} = 0;
}


# conform sample columns
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];

    # strip mzXML from sample names
    $field =~ s/\.mzX?ML( Peak \S+)$/$1/i;
    $field =~ s/\.cdf( Peak \S+)$/$1/i;

    if ($field =~ / Peak (\S+)$/i)
    {

        $second_word = $1;
        $second_word = lc $second_word;
        $field =~ s/[ ._]+Peak (\S+)$/ Peak $second_word/i;
    }
    $field =~ s/[ ._]+$//;


    # rename samples from comment lines in lipidomics data
    if ($field =~ /\[([^]]+)\]/)
    {
        $new_name = $sample_rename_hash{$1};
        if (defined($new_name))
        {
            $field =~ s/\[[^]]+\]/\[$new_name\]/g;
        }
    }

    $header_col_array[$col] = $field;
}
$num_header_cols = @header_col_array;


# scan for peak height and peak area
$peak_height_flag = 0;
$peak_area_flag   = 0;

for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];

    if ($field =~ / Peak height$/i)
    {
        $peak_height_flag = 1;
        next;
    }
    if ($field =~ / Peak area$/i)
    {
        $peak_area_flag = 1;
        next;
    }
    
    # lipidomics
    if ($field =~ /^Area[,[]/i)
    {
         $peak_area_flag = 1;
         next;
    }
    if ($field =~ /^Height[,[]/i)
    {
         $peak_height_flag = 1;
         next;
    }
}


# conform the Height[] and Area[] columns so the rest of the pipeline is happy
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];
    
    if ($field =~ /^Height\[([^]]+)\]$/)
    {
        $header_col_array[$col] = $1 . ' Peak height';
    }
    if ($field =~ /^Area\[([^]]+)\]$/)
    {
        $header_col_array[$col] = $1 . ' Peak area';
    }
}


# uh oh, no data columns found
# conform them so that Area or Height start at the beginning or end
if ($peak_height_flag == 0 && $peak_area_flag == 0)
{
    printf STDERR "WARNING -- non-standard Height/Area nomenclature, conforming with best guess\n";

    for ($col = 0; $col < @header_col_array; $col++)
    {
        $field = $header_col_array[$col];

        if ($field =~ /([^A-Za-z0-9]*Height[^A-Za-z0-9]*)/i)
        {
            $height_area_str = $1;

            # conform the sample header
            # string may contain () or [], so escape it with \Q \E
            $field =~ s/\Q$height_area_str\E/_/;
            
            # deal with inserted _ at beginning/end
            if ($field =~ /^(_+)/)
            {
                $count_underscores_new  = length $1;
                $count_underscores_orig = 0;
                
                if ($header_col_array[$col] =~ /^(_+)/)
                {
                    $count_underscores_orig = length $1;
                }
                
                if ($count_underscores_new != $count_underscores_orig)
                {
                    $field =~ s/^_//;
                }
            }
            if ($field =~ /(_+)$/)
            {
                $count_underscores_new  = length $1;
                $count_underscores_orig = 0;
                
                if ($header_col_array[$col] =~ /(_+)$/)
                {
                    $count_underscores_orig = length $1;
                }
                
                if ($count_underscores_new != $count_underscores_orig)
                {
                    $field =~ s/_$//;
                }
            }

            # deal with dangling ( or ]
            if ($field =~ /\]$/ && !($field =~ /\[/))
            {
                $field =~ s/\]$//;
            }
            if ($field =~ /\)$/ && !($field =~ /\(/))
            {
                $field =~ s/\)$//;
            }
            
            # rename samples from comment lines
            $new_name = $sample_rename_hash{$field};
            if (defined($new_name))
            {
                $field = $new_name;
            }
            
            $field .= ' Peak height';

            $header_col_array[$col] = $field;
            $peak_height_flag = 1;
        }
        if ($field =~ /([^A-Za-z0-9]*Area[^A-Za-z0-9]*)/i)
        {
            $height_area_str = $1;
        
            # conform the sample header
            # string may contain () or [], so escape it with \Q \E
            $field =~ s/\Q$height_area_str\E/_/;

            # deal with inserted _ at beginning/end
            if ($field =~ /^(_+)/)
            {
                $count_underscores_new  = length $1;
                $count_underscores_orig = 0;
                
                if ($header_col_array[$col] =~ /^(_+)/)
                {
                    $count_underscores_orig = length $1;
                }
                
                if ($count_underscores_new != $count_underscores_orig)
                {
                    $field =~ s/^_//;
                }
            }
            if ($field =~ /(_+)$/)
            {
                $count_underscores_new  = length $1;
                $count_underscores_orig = 0;
                
                if ($header_col_array[$col] =~ /(_+)$/)
                {
                    $count_underscores_orig = length $1;
                }
                
                if ($count_underscores_new != $count_underscores_orig)
                {
                    $field =~ s/_$//;
                }
            }

            # deal with dangling ( or ]
            if ($field =~ /\]$/ && !($field =~ /\[/))
            {
                $field =~ s/\]$//;
            }
            if ($field =~ /\)$/ && !($field =~ /\(/))
            {
                $field =~ s/\)$//;
            }

            # rename samples from comment lines
            $new_name = $sample_rename_hash{$field};
            if (defined($new_name))
            {
                $field = $new_name;
            }

            $field .= ' Peak area';

            $header_col_array[$col] = $field;
            $peak_area_flag = 1;
        }
    }
}


# we're going to use the compound names to identify spikeins
# and signal/background peaks
$rowid_col = $header_col_hash{'row ID'};
$name_col  = $header_col_hash{'row identity (all IDs)'};


if (!defined($name_col))
{
    $name_col = $header_col_hash{'row identity (main ID)'};
    
    if (defined($name_col))
    {
        printf STDERR "WARNING -- (main ID) name used, will miss multiple hits per row\n";
    }
}
if (!defined($name_col))
{
    $name_col = $header_col_hash{'compound'};
}
if (!defined($name_col))
{
    $name_col = $header_col_hash{'compoundId'};
}

# lipidomics
$lipidsearch_flag = 0;
if (!defined($name_col))
{
    $name_col = $header_col_hash{'LipidIon'};
    $lipidsearch_flag = 1;
}
if (!defined($name_col))
{
    $name_col = $header_col_hash{'LipidGroup'};
    $lipidsearch_flag = 1;
}
if (!defined($name_col))
{
    $name_col = $header_col_hash{'IonFormula'};
    $lipidsearch_flag = 1;
}
if (!defined($name_col))
{
    $name_col = $header_col_hash{'LipidMolec'};
    $lipidsearch_flag = 1;
}


# sometimes use formula for heavy label scanning instead of name
if ($lipidsearch_flag)
{
    $formula_col = $header_col_hash{IonFormula};
}
# scan for formula col
if (!defined($formula_col))
{
    for ($i = 0; $i < @header_col_array; $i++)
    {
        $header = $header_col_array[$i];

        if (!defined($formula_col) && $header =~ /formula/i)
        {
            $formula_col = $i;
        }
    }
}


# flag columns to remove, as they clutter up the file
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];
    
    if ($field =~ /\.mzX?ML[^.]+$/i ||
        $field =~ /\.cdf[^.]+$/i ||
        $field =~ / Peak \S+$/ ||
        $field =~ /row identity/ ||
        $field =~ /row comment/ ||
        $field =~ /^Area[,[]/i)
    {
        # always keep peak height
        if ($field =~ / Peak height$/i ||
            $field =~ /^Height[,[]/i)
        {
            next;
        }

        # only keep peak area if no peak height
        if ($peak_height_flag == 0 &&
            ($field =~ / Peak area$/i ||
             $field =~ /^Area[,[]/i))
        {
            next;
        }
        
        # only keep the chosen row identify column
        if ($col eq $name_col)
        {
            next;
        }
 
        $col_to_remove_hash{$col} = 1;
    }
    
    # mzMine exports a blank column at the end of every row,
    # which can majorly screw up other software
    # just remove ALL blank column headers
    if ($field eq '')
    {
        $col_to_remove_hash{$col} = 1;
    }
}


# heavy/light labeled El-MAVEN data
$isotope_col = $header_col_hash{isotopeLabel};


$first_abundance_col = 9E99;
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];

    if ($field =~ / Peak height$/i ||
        $field =~ / Peak area$/i ||
        $field =~ /^Area[,[]/i ||
        $field =~ /^Height[,[]/i)
    {
        if (!defined($col_to_remove_hash{$col}))
        {
            $sample_col_hash{$col} = 1;
        
            if ($col < $first_abundance_col)
            {
                $first_abundance_col = $col;
            }
        }
    }
}

# none found, check for pos/neg
if ($first_abundance_col == 9E99)
{
  for ($col = 0; $col < @header_col_array; $col++)
  {
    $field = $header_col_array[$col];

    if ($field =~ /(^|[^A-Za-z0-9]+)(pos|neg)([^A-Za-z0-9]+|$)/i ||
        $field =~ /(pos|neg)$/i)
    {
        if (!defined($col_to_remove_hash{$col}))
        {
            $sample_col_hash{$col} = 1;

            if ($col < $first_abundance_col)
            {
                $first_abundance_col = $col;
            }
        }
    }
  }
}

# uh oh, column headers don't have anything obviously abundance-looking
# use the known-last metadata column to denote where samples start
if ($first_abundance_col == 9E99)
{
    $col = $header_col_hash{'Non-heavy identified flag'};
    if (defined($col))
    {
        printf STDERR "WARNING -- using pipeline known-last metadata column to locate data columns\n";
        $first_abundance_col = $col + 1;

        for ($col = 0; $col < $first_abundance_col; $col++)
        {
            $field = $header_col_array[$col];

            #if ($field =~ /\S/)
            #{
            #    $metadata_col_hash{$col} = 1;
            #}
        }
        for ($col = $first_abundance_col; $col < $num_header_cols; $col++)
        {
            $field = $header_col_array[$col];

            if ($field =~ /\S/)
            {
                $sample_col_hash{$col} = 1;
            }
        }
    }
}


printf STDERR "First abundance col: %s\n", $first_abundance_col;


@sample_col_array = sort {$a<=>$b} keys %sample_col_hash;
$n_total = @sample_col_array;


# scan for number of detected peaks prior to gap filling
$has_pregap_flag = 0;
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];
    
    if ($field =~ /row number of detected peaks/i)
    {
        $has_pregap_flag = 1;
        
        $pregap_col = $col;
    }
}

# it might be an El-Maven file
if ($has_pregap_flag == 0)
{
    for ($col = 0; $col < @header_col_array; $col++)
    {
        $field = $header_col_array[$col];
    
        if ($field =~ /goodPeakCount/i)
        {
            $has_pregap_flag = 1;
        
            $pregap_col = $col;
        }
    }
}


if ($has_pregap_flag == 0)
{
    print STDERR "WARNING -- \"one-hit wonder\" detection disabled;\n";
    print STDERR "           'row number of detected peaks' or 'goodPeakCount' columns not found\n";
}



if (!defined($name_col))
{
    printf STDERR "ABORT -- can't find 'row identity (all IDs)' column\n";
    exit(1);
}


if (!defined($rowid_col))
{
    printf STDERR "No row identifier column found, creating arbitrary identifiers\n";
}


# print missing rowid
$print_flag = 0;
if (!defined($rowid_col))
{
    print "row ID";
    $print_flag = 1;
}

# print header line, non-sample data first
for ($col = 0; $col < @header_col_array; $col++)
{
    if (defined($col_to_remove_hash{$col}))
    {
        next;
    }
    if (defined($sample_col_hash{$col}))
    {
        next;
    }
    
    if ($print_flag == 1)
    {
        print "\t";
    }
    else
    {
        $print_flag = 1;
    }

    print $header_col_array[$col];

    # insert new ParentFormula column
    if ($lipidsearch_flag && defined($formula_col) &&
        $col == $formula_col)
    {
        print "\tParentFormula";
    }
}


# insert new headers
if ($print_flag == 1)
{
    print "\t";
}
else
{
    $print_flag = 1;
}
print "Number of non-zero peaks";
if ($has_pregap_flag)
{
    print "\tPercent pre-gap peaks";
}
print "\tPercent non-zero peaks";
print "\tHeavy-labeled flag";
print "\tIdentified flag";
print "\tNon-heavy identified flag";


# print sample headers
for ($col = 0; $col < @header_col_array; $col++)
{
    if (defined($col_to_remove_hash{$col}))
    {
        next;
    }
    if (!defined($sample_col_hash{$col}))
    {
        next;
    }
    
    if ($print_flag == 1)
    {
        print "\t";
    }
    else
    {
        $print_flag = 1;
    }

    print $header_col_array[$col];
}
print "\n";



open OUTFILE_EXCLUSIONS, ">$output_exclusions_filename" or die "can't open file $output_exclusions_filename for writing\n";
open OUTFILE_SPIKEINS,   ">$output_spikeins_filename"   or die "can't open file $output_spikeins_filename for writing\n";


$num_excel_large = 0;	# number of potential Excel-corruptible values
$sum_excel_large = 0;	# sum of number of non-zero last 5 digits

$row_count               = 0;
$too_low_blank_count     = 0;
$too_low_non_blank_count = 0;

while(defined($line=<INFILE>))
{
    $row_count++;

    $line =~ s/[\r\n]+//g;
    $line =~ s/\"//g;

    @array = split /\t/, $line, -1;    # don't skip empty fields at and

    $print_flag = 0;

    # clean up fields
    for ($col = 0; $col < @array; $col++)
    {
        $array[$col] =~ s/^\s+//;
        $array[$col] =~ s/\s+$//;
        $array[$col] =~ s/\s+/ /g;
        
        $array[$col] = reformat_sci($array[$col]);
    }
    

    if (defined($rowid_col))
    {
        $rowid = $array[$rowid_col];
    }
    else
    {
        $rowid = $row_count;
    }


    $name    = $array[$name_col];
    $formula = '';
    
    if (defined($formula_col))
    {
        $formula = $array[$formula_col];
        $formula = conform_formula($formula);

        $array[$formula_col] = $formula;
    }
    
    # remove adduct from formula to yield parent formula
    $formula_parent = '';
    if (defined($formula_col) && $name =~ /([+-][A-Za-z0-9+-]+)$/)
    {
        $adduct = $1;

        $formula_parent = conform_formula($formula, $adduct);
    }


    # flag heavy labeled rows (including spike-ins)
    #
    # probably MZMine data if there is no isotopeLabel column
    #  could also be El-MAVEN without heavy labeling enabled,
    #  in which case is_heavy_labeled() may always return false
    $heavy_flag = 0;
    
    # use formula column on lipidomics data
    if ($lipidsearch_flag && defined($formula_col))
    {
        if (is_heavy_labeled($formula))
        {
            $heavy_flag = 1;
        
            if ($scale_heavy_flag == 0)
            {
                print OUTFILE_SPIKEINS "$rowid\n";
            }
        }
    }
    # use name column on MZmine data, or El-MAVEN with no isotope column
    elsif (!defined($isotope_col))
    {
        if (is_heavy_labeled($name))
        {
            $heavy_flag = 1;
        
            if ($scale_heavy_flag == 0)
            {
                print OUTFILE_SPIKEINS "$rowid\n";
            }
        }
    }
    # use isotope column on El-MAVEN data
    else
    {
        $isotope = $array[$isotope_col];
        
        # C12 PARENT
        # C13-label-#
        #
        # I've only seen C13 experiments so far, so, rather than
        # make a list of all heavy labels we may see, I'm going to key off
        # of "PARENT" and "label" instead, to hopefully future-proof it
        if ($isotope =~ /PARENT/i)
        {
            $heavy_flag = 0;
        }
        elsif ($isotope =~ /label/i)
        {
            $heavy_flag = 1;

            if ($scale_heavy_flag == 0)
            {
                print OUTFILE_SPIKEINS "$rowid\n";
            }
        }
    }

    
    # store identified/unidentified
    # consider it identified if it has two letters in a row
    # this should result in treating purely chemical formulas as unidentified
    #
    # include heavy labeled metabolites in the unidentified list, since
    # we're going to use that for iron training exclusions
    $identified_flag = 0;
    if (($name =~ /[A-Za-z][A-Za-z]/ ||
         $name =~ /[A-Za-z][0-9]/) &&
        !($name =~ /\d m\/z adduct of \d/) &&
        !($name =~ /Complex of [0-9.]+ and \d+/))
    {
        $identified_flag = 1;
        
        if ($heavy_flag)
        {
            print OUTFILE_EXCLUSIONS "$rowid\n";
        }
    }
    else
    {
        print OUTFILE_EXCLUSIONS "$rowid\n";
    }
    
    $nsid_flag = 0;
    if ($heavy_flag == 0 && $identified_flag == 1)
    {
        $nsid_flag = 1;
    }


    # count number of non-zero samples
    $n = 0;
    foreach $col (@sample_col_array)
    {
        $value = $array[$col];
        
        if (is_number($value) && $value != 0)
        {
            $n++;
        }
    }
    
    $n_percent = 0;
    if ($n_total)
    {
        $n_percent = 100 * $n / $n_total;
    }

    $n_pregap = 0;
    $n_percent_pregap = 0;    
    if ($has_pregap_flag && $n_total)
    {
        $n_pregap = $array[$pregap_col];
        
        if (!is_number($n_pregap))
        {
            $n_pregap = 0;
        }
    
        $n_percent_pregap = 100 * $n_pregap / $n_total;
    }
    
    # skip "one-hit wonders"
    if ($has_pregap_flag && $keep_single_pregap_flag == 0)
    {
        # must use == 1 instead of <= 1, since 0 means manually assigned
        # should I change 0 into the number of manually picked peaks?
        # or leave 0 as 0, as I am currently doing?
        if ($n_pregap == 1)
        {
            next;
        }
    }
    
    # discard unidentified rows
    if ($discard_unidentified_flag && $identified_flag == 0)
    {
        next;
    }
    # discard heavy labeled rows
    if ($discard_heavy_flag && $heavy_flag)
    {
        next;
    }
    
    
    # set seen heavy flag if any heavy metabolites made it through
    if ($heavy_flag)
    {
        $seen_heavy_flag = 1;
    }


    # print missing rowid
    if (!defined($rowid_col))
    {
        print $rowid;
        $print_flag = 1;
    }

    # print non-sample data
    for ($col = 0; $col < $num_header_cols; $col++)
    {
        if (defined($col_to_remove_hash{$col}))
        {
            next;
        }
        if (defined($sample_col_hash{$col}))
        {
            next;
        }

        $field = $array[$col];
        if (!defined($field))
        {
            $field = '';
        }
    
        if ($print_flag == 1)
        {
            print "\t";
        }
        else
        {
            $print_flag = 1;
        }
        
        print $field;

        # insert new ParentFormula column
        if ($lipidsearch_flag && defined($formula_col) &&
            $col == $formula_col)
        {
            print "\t$formula_parent";
        }
    }


    # insert new fields
    if ($print_flag == 1)
    {
        print "\t";
    }
    else
    {
        $print_flag = 1;
    }
    print "$n";
    if ($has_pregap_flag)
    {
        print "\t$n_percent_pregap";
    }
    print "\t$n_percent";
    print "\t$heavy_flag";
    print "\t$identified_flag";
    print "\t$nsid_flag";


    # print sample data
    for ($col = 0; $col < $num_header_cols; $col++)
    {
        if (defined($col_to_remove_hash{$col}))
        {
            next;
        }
        if (!defined($sample_col_hash{$col}))
        {
            next;
        }

        $field = $array[$col];
        if (!defined($field))
        {
            $field = '';
        }
        
        # check for bad data
        # mass spec data is generally at least in the 1000's
        # a value < 10 is almost certainly garbage
        if (is_number($field) && $field > 0 && $field < $floor_cutoff)
        {
            $field = '0';
        
            $sample = $header_col_array[$col];

            # don't warn about blank samples

            # if ($sample =~ /processing_bla?n?k\d*([^A-Za-z0-9]|$)/i ||
            #     $sample =~ /(^|[^A-Za-z0-9])blank\d*([^A-Za-z0-9]|$)/i)
            if ($sample =~ /Bla?n?k(\b|[A-Z_])/ ||
                $sample =~ /(\b|_)bla?n?k(\b|_)/i)
            {
                $too_low_blank_count++;
            }
            else
            {
                $too_low_non_blank_count++;
            }
        }
    
        if ($print_flag == 1)
        {
            print "\t";
        }
        else
        {
            $print_flag = 1;
        }
        
        print $field;
    }
    print "\n";


    # scan samples for number of significant digits
    foreach $col (@sample_col_array)
    {
        $value = $array[$col];
        if (!defined($value))
        {
            next;
        }
        
        # largeish numbers succeptible to Excel corruption
        if (is_number($value) && $value != 0 &&
            $value >= 1E7)
        {
            if ($value =~ /([0-9]{5})$/)
            {
                $last_five_digits = $1;

                $num_excel_large += 1;

                @temp_array = ($last_five_digits =~ m/[1-9]/g);
                $count = @temp_array;
                
                $sum_excel_large += $count;
            }
        }
    }
}


if ($too_low_blank_count)
{
    printf STDERR "NOPROBLEM -- %d value(s) < %d were floored to zero in blank samples\n",
        $too_low_blank_count, $floor_cutoff;
}
if ($too_low_non_blank_count)
{
    printf STDERR "LOWSIGNAL -- %d value(s) < %d were floored to zero in non-blank samples\n",
        $too_low_non_blank_count, $floor_cutoff;
}


$avg_excel_large = 9E99;
if ($num_excel_large)
{
    $avg_excel_large = $sum_excel_large / $num_excel_large;
}

if ($avg_excel_large < 3.1)
{
    printf STDERR "WARNING -- !!! Excel-corrupted values >= 1E7 detected !!!\n";
}

if ($seen_heavy_flag && $scale_heavy_flag == 0)
{
    printf STDERR "CAUTION -- treating heavy rows as spike-ins: leave unnormalized\n";
    printf STDERR "           use --heavy-tracer in isotope tracer experiments\n";
}
if ($seen_heavy_flag && $scale_heavy_flag == 1)
{
    printf STDERR "CAUTION -- treating heavy rows as biological: normalize them\n";
    printf STDERR "           use --heavy-spikein if they are all spike-ins\n";
}
