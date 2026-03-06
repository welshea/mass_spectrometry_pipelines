#!/usr/bin/perl -w


# set lib search path to directory the script is run from
use File::Basename;
use lib dirname (__FILE__);

use Scalar::Util qw(looks_like_number);
use POSIX;

# force lines to be flushed
# otherwise, STDOUT and STDERR can blend together on same output line...
select(STDERR); $| = 1;
select(STDOUT); $| = 1;


# which elements to add/subtract for each supported adduct
$adduct_hash{'[M+H]+'}{H}      = +1;
#$adduct_hash{'[M+H-2H2O]+'}{H} = -3;
#$adduct_hash{'[M+H-2H2O]+'}{O} = -2;
#$adduct_hash{'[M+NH4]+'}{H}    = +4;   # lose 4 from water, gain one from +H
#$adduct_hash{'[M+NH4]+'}{N}    = +1;
$adduct_hash{'[M-H]-'}{H}      = -1;


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
sub modify_formula
{
    my $formula_orig = $_[0];
    my $label        = $_[1];    # labeled element
    my $adduct       = $_[2];    # adduct to add to the formula
    my $formula;
    my $formula_new = '';
    my @match_array;
    my @element_array;
    my @adduct_ele_array;
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

        if (length $element > 2)
        {
            printf STDERR "WARNING -- error in formula:\t%s\n", $formula_orig;
            return $formula_orig;
        }


        if (defined($count_hash{$heavy_plus_element}))
        {
            printf STDERR "WARNING -- duplicate elements in formula:\t%s\n", $formula_orig;
            
            # add to existing counts
            $count_hash{$heavy_plus_element} += $count;
        }
        else
        {
            $count_hash{$heavy_plus_element} = $count;
        }
    }
    
    # modify count hash as appropriate for adduct
    if (defined($adduct) && $adduct =~ /[A-Za-z]/)
    {
        # I am lazy and don't want to write complicated generic adduct code;
        # only adducts defined in our adduct hash are supported for now
        if (defined($adduct_hash{$adduct}))
        {
            @adduct_ele_array = keys %{$adduct_hash{$adduct}};
            
            foreach $element (@adduct_ele_array)
            {
                $offset = $adduct_hash{$adduct}{$element};
                
                $count_hash{$element} += $offset;
                
                # zero counts left for this element
                if ($count_hash{$element} == 0)
                {
                    delete $count_hash{$element};
                }
                elsif ($count_hash{$element} < 0)
                {
                    printf STDERR "WARNING -- impossible adduct, edit data and re-run:\t|%s|\t|%s|\n",
                        $adduct, $formula_orig
                }
            }
        }
        else
        {
            printf STDERR "WARNING -- unsupported adduct, formula left as-is:\t%s\t%s\n",
                $adduct, $formula_orig;
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

                ## IsoCorrectoR requires a number after every element
                ## replace 1 count with blank
                #if ($count == 1)
                #{
                #    $count = '';
                #}

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

                ## IsoCorrectoR requires a number after every element
                ## replace 1 count with blank
                #if ($count == 1)
                #{
                #    $count = '';
                #}

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
        
        ## IsoCorrectoR requires a number after every element
        ## replace 1 count with blank
        #if ($count == 1)
        #{
        #    $count = '';
        #}
        
        $formula_new .= $heavy_plus_element . $count;
    }
    
    # IsoCorrectoR:
    #   add LabX# to the end, X is the label element, # is max label count
    $count        = $count_hash{$label};
    $formula_new .= Lab . $label . $count;
    
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
    if ($string =~ /\([^()]*[)-]D[0-9]+\b[^()]*\)/) { return 1; }
    if ($string =~ /\([^()]*\bBOC\b[^()]*\)/)       { return 1; }

    if ($string =~ /\b13C[0-9]+\b/) { return 1; }
    if ($string =~ /\b14N[0-9]+\b/) { return 1; }
    if ($string =~ /[)-]D[0-9]+\b/) { return 1; }
    if ($string =~ /\bBOC\b/)       { return 1; }
    
    return 0;
}


sub output_molecule_file
{
    my @name_array;
    my @formula_array;
    my $name;
    my $label;
    my $adduct;
    my $formula;
    my $formula_new;
    my $formula_str;
    my $name_escaped;

    @name_mol_array = sort { cmp_args_alphanumeric($a, $b) }
        keys %name_mol_hash;
    
    foreach $name_mol (@name_mol_array)
    {
        $formula = $name_mol_hash{$name_mol}{formula};
        $adduct  = $name_mol_hash{$name_mol}{adduct};
        
        @label_array = sort { $a cmp $b }
            keys %{$name_mol_hash{$name_mol}{element_hash}};
        
        ## for now, we only support a single label element
        $label   = $label_array[0];
        
        $formula_new = modify_formula($formula, $label, $adduct);

        $name_mol_escaped = $name_mol;
        if ($name_mol =~ /,/ || $name_mol =~ /\"/)
        {
            $name_mol_escaped =~ s/\"/\"\"/g;
            $name_mol_escaped = '"' . $name_mol_escaped . '"';
        }

        $formula_str = $formula_new;
        if ($formula =~ /,/ || $formula =~ /\"/)
        {
            $formula_str =~ s/\"/\"\"/g;
            $formula_str = '"' . $formula_str . '"';
        }
        
        print MOL_FILE "$name_mol_escaped,$formula_str,\n";
    }
}


# begin main()

$num_files = 0;


$output_prefix = 'isocorrector';
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
            $data_filename = $field;
            $num_files++;
        }
        elsif ($num_files == 1)
        {
            $output_prefix = $field;
            $num_files++;
        }
    }
}



if (!defined($data_filename) || $syntax_error_flag)
{
    $program_name = basename($0);

    print STDERR "Usage: $program_name merged_pipeline_log2_output.txt [output_file_prefix]\n";
    exit(1);
}


open DATA,       "$data_filename"       or die "can't open $data_filename\n";


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
@array = split /\t/, $line;    # skip empty headers at and
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    $array[$i] =~ s/\s+/ /g;
    $array[$i] =~ s/^\"(.*)\"$/$1/;
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


# formula column isn't the first encountered formula column
# FormulaMapped comes before it,
# so check for it first, so that the correct formula column is used
$data_formula_col = $data_header_col_hash{'formula'};
if (!defined($data_formula_col))
{
    $data_formula_col = $data_header_col_hash{'Formula'};
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
    if (!defined($data_isotope_col) &&
        $header =~ /isotopeLabel/i)
    {
        $data_isotope_col = $i;
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
    # fallback to row ID column if pos/neg column not present
    # identifiers will, after pipeline processing, start with pos_ / neg_
    elsif (!defined($data_pos_neg_col) &&
           $header =~ /^row ID$/i)
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
$lipidomics_flag = 0;
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

if (!defined($data_name_col))
{
    $data_name_col = $data_header_col_hash{'Name'};
}
if (!defined($data_name_col))
{
    $data_name_col = $data_header_col_hash{'name'};
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


$data_rowid_col     = $data_header_col_hash{'row ID'};
$data_adduct_col    = $data_header_col_hash{'adductName'};
$data_nhid_col      = $data_header_col_hash{'Non-heavy identified flag'};
$data_bad_label_col = $data_header_col_hash{'label'};
$data_metagroup_col = $data_header_col_hash{'metaGroupId'};
$data_parent_col    = $data_header_col_hash{'parent'};


if (!defined($data_rowid_col))
{
    printf STDERR "ABORT -- row ID column not found in data file %s\n",
                 $data_filename;
    exit(1);
}


# determine where sample columns start
$data_sample_start_col = 1;
if (defined($data_mz_col) && $data_mz_col > $data_sample_start_col)
{
    $data_sample_start_col = $data_mz_col;
}
if (defined($data_name_col) && $data_name_col > $data_sample_start_col)
{
    $data_sample_start_col = $data_name_col;
}
if (defined($data_rt_col) && $data_rt_col > $data_sample_start_col)
{
    $data_sample_start_col = $data_rt_col;
}
if (defined($data_heavy_col) && $data_heavy_col > $data_sample_start_col)
{
    $data_sample_start_col = $data_heavy_col;
}
if (defined($data_pos_neg_col) && $data_pos_neg_col > $data_sample_start_col)
{
    $data_sample_start_col = $data_pos_neg_col;
}
if (defined($data_formula_col) && $data_formula_col > $data_sample_start_col)
{
    $data_sample_start_col = $data_formula_col;
}
if (defined($data_nhid_col))
{
    $data_sample_start_col = $data_nhid_col;
}
$data_sample_start_col += 1;

printf STDERR "First sample column:\t%d\t%s\n",
    $data_sample_start_col, $data_header_col_array[$data_sample_start_col];


$iso_mol_filename  = sprintf "%s_MoleculeFile.csv",    $output_prefix;
$iso_meas_filename = sprintf "%s_MeasurementFile.csv", $output_prefix;

open MOL_FILE,  ">$iso_mol_filename"  or die "ABORT -- cannot open $iso_mol_filename for writing\n";
open MEAS_FILE, ">$iso_meas_filename" or die "ABORT -- cannot open $iso_meas_filename for writing\n";


# write headers to output CSV files
print MOL_FILE    "Molecule";
print MOL_FILE    ",MS ion or MS/MS product ion";
print MOL_FILE    ",MS/MS neutral loss";
print MOL_FILE    "\n";

print MEAS_FILE "Measurements/Samples";
for ($i = $data_sample_start_col; $i < @data_header_col_array; $i++)
{
    $header = $data_header_col_array[$i];
    
    # escape for CSV
    if ($header =~ /,/ || $header =~ /\"/)
    {
        $header =~ s/\"/\"\"/g;
        $header = '"' . $header . '"';
    }
    
    print MEAS_FILE ",$header";
}
print MEAS_FILE "\n";


# read in data
$row_data          = 0;
%seen_isotope_hash = ();    # have we seen this isotope in the current group
$group_id          = 'first_group';
$name_old          = '__old_name__';
$metagroup_old     = '__old_metagroup__';
$parent_old        = '__old_parent__';

while(defined($line=<DATA>))
{
    $line =~ s/[\r\n]+//g;
    @array = split /\t/, $line, -1;    # don't skip empty fields at and

    # clean up fields
    for ($col = 0; $col < @array; $col++)
    {
        $array[$col] =~ s/^\s+//;
        $array[$col] =~ s/\s+$//;
        $array[$col] =~ s/\s+/ /g;
        $array[$col] =~ s/^\"(.*)\"$/$1/;
        $array[$col] =~ s/^\s+//;
        $array[$col] =~ s/\s+$//;
        $array[$col] =~ s/\s+/ /g;
        
        $array[$col] = reformat_sci($array[$col]);
    }

    $rowid     = $array[$data_rowid_col];
    $name      = $array[$data_name_col];
    $isotope   = $array[$data_isotope_col];

    $adduct    = '';
    $metagroup = '';
    $parent    = '';

    if (defined($data_adduct_col))
    {
        $adduct    = $array[$data_adduct_col];
    }
    if (defined($data_metagroup_col))
    {
        $metagroup = $array[$data_metagroup_col];
    }
    if (defined($data_parent_col))
    {
        $parent    = $array[$data_parent_col];
    }


    $pos_neg = 'unk';
    if (defined($data_pos_neg_col))
    {
        $pos_flag = 0;
        $neg_flag = 0;
        if ($array[$data_pos_neg_col] =~ /^pos/i)
        {
            $pos_flag = 1;
        }
        if ($array[$data_pos_neg_col] =~ /^neg/i)
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
    
    
    # assign current group id, store some info for it
    if ($adduct =~ /[A-Z]/ ||
        $name      ne $name_old ||
        $metagroup ne $metagroup_old ||
        $parent    ne $parent_old ||
        defined($seen_isotope_hash{$isotope}))
    {
        $group_id = $rowid;
        
        # guess adduct if it is missing
        if (!($adduct =~ /[A-Z]/))
        {
            if ($pos_neg eq 'pos')
            {
                $adduct = '[M+H]+';
            }
            if ($pos_neg eq 'neg')
            {
                $adduct = '[M-H]-';
            }
        }
        
        $group_hash{$group_id}{adduct} = $adduct;
        $group_hash{$group_id}{name}   = $name;
        
        %seen_isotope_hash = ();
    }
    
    $seen_isotope_hash{$isotope} = 1;
    $name_old                    = $name;
    $metagroup_old               = $metagroup;
    $parent_old                  = $parent;
    
    # rows flagged as bad labels
    $bad_label = '';
    if (defined($data_bad_label_col))
    {
        $bad_label = $array[$data_bad_label_col];
    }
    
    # Leave bad peaks in for now.
    # Most (~90%) of them are noise-level low-abundance peaks.
    # The rest appear to be messy peaks that smear into each other.
    # Either way, I think having the less-exactly-quantified abundances
    # should help the isotopolog abundance corrections more than
    # leaving them out.  IsoCorrectoR can handle missing data,
    # but my gut feeling is that noisy labeled peaks are still better than
    # entirely missing labeled peaks.
    #
    #if ($bad_label =~ /b/i)
    #{
    #    next;
    #}
    
    
    # retention time sanity checks, if available
    $rt_data = '';
    if (defined($data_rt_col))
    {
        $rt_data = $array[$data_rt_col];
    }
    
    # formula, for sanity checks, if available
    if (defined($data_formula_col))
    {
        $formula = $array[$data_formula_col];
    }
    
    

    # Citicoline had an incorrect formula in our database
    # causing the peak assignments to be bogus
    # discard the row entirely
    if ($formula =~ /C14H27N4O11P2/i &&
        ($name =~ /\bCiticoline\b/i ||
         $name =~ /\bCytidine 5'-Diphosphocholine\b/i))
    {
        next;
    }


    # some ions have different adducts
    # this should all turn into +/- H due to the formulas already having
    # removed the 2H2O
    if (0 && $pos_neg eq 'pos')
    {
        if ($formula eq 'C24H36O2' &&
            $name =~ /Chenodeoxycholic acid (CDCA)/i)
        {
            next;
        }
        elsif ($formula eq 'C24H36O2' &&
               $name =~ /Deoxycholate (DCA)/i)
        {
            next;
        }
        elsif ($formula eq 'C24H43O5N' &&
               $name =~ /Cholate (DCA)/i)
        {
            next;
        }
        elsif ($formula eq 'C26H48N2O6S' &&
               $name =~ /Taurochenodeoxycholic acid (TCDCA)/i)
        {
            next;
        }
        elsif ($formula eq 'C26H48N2O7S' &&
               $name =~ /Taurocholic acid (TCA)/i)
        {
            next;
        }
    }
    
    
    # convert isotopes into IsoCorrectoR format
    # TODO -- only supports single element for now
    $element   = '';
    $ele_count = 0;
    if ($isotope =~ /([A-Za-z]+)[0-9]+\s*PARENT/i)
    {
        $element   = $1;
        $ele_count = 0;
    }
    if ($isotope =~ /([A-Za-z]+)[0-9]+-label-([0-9]+)/i)
    {
        $element   = $1;
        $ele_count = $2;
    }
    
    if ($element eq '')
    {
        printf STDERR "WARNING -- unsupported label:\t%s\t%s\n",
            $isotope, $name;

        next;
    }

    # uniquify metabolite group using row identifier
    $name_mol =  sprintf "%s %s",   $name, $group_id;
    
    # Reformat name to what IsoCorrectoR requires.
    # IsoCorrectoR cannot handle + signs or ()[]{} in molecule names.
    # !!! It tries to interpret molecule names as regular expressions !!!
    # Strip underscores just in case that might cause problems with parsing.
    # Also, out of parnoia, deal with more special characters that are
    # unlikely to ever occur.
    #
    # WARNING -- Backslashes are not currently dealt with, but should not
    # be a problem, as they are unlikely to occur in real data?
    #$name_mol =~ s/\+/\-plus\-/g;   # replace reserved + with the word plus
    #$name_mol =~ s/[\(\)\[\]\{\}]/#/g;   # use hash mark for various brackets
    $name_mol =~ s/_//g;            # paranoia over _Tracer# parser
    #$name_mol =~ s/[\.\|]/;/g;      # replace other reserved characters
    #$name_mol =~ s/[\*\?\^\$]//g;   # strip remaining reserved characters
    #$name_mol =~ s/\-+/\-/g;        # there could be sequential dashes now
    $name_new =  sprintf "%s_%s%s", $name_mol, $element, $ele_count;


    # escape for CSV
    $name_new_escaped = $name_new;
    if ($name_new =~ /,/ || $name_new =~ /\"/)
    {
        $name_new_escaped =~ s/\"/\"\"/g;
        $name_new_escaped = '"' . $name_new_escaped . '"';
    }

    # write to IsoCorrectoR Measurements file
    print MEAS_FILE "$name_new_escaped";
    for ($i = $data_sample_start_col; $i < @data_header_col_array; $i++)
    {
        $field = $array[$i];
        
        # input is log2 abundances, unlog them
        if (is_number($field))
        {
            $field = 2**$field;
        }

        # escape for CSV
        if ($field =~ /,/ || $field =~ /\"/)
        {
            $field =~ s/\"/\"\"/g;
            $field = '"' . $field . '"';
        }

        print MEAS_FILE ",$field";
    }
    print MEAS_FILE "\n";


    # store formula and adduct for molecule file later
    $name_mol_hash{$name_mol}{formula} = $formula;
    $name_mol_hash{$name_mol}{adduct}  = $adduct;
    $name_mol_hash{$name_mol}{element_hash}{$element} = 1;

    # store mappings for output mapping table later
    $mapping_hash{$name_new} = $rowid;

    $row_data++;
}

output_molecule_file();


# output mapping table
$mapping_filename = sprintf "%s_id_mapping_table.txt", $output_prefix;
open MAPPING_FILE,  ">$mapping_filename"  or die "ABORT -- cannot open $mapping_filename for writing\n";

@name_new_array = sort { cmp_args_alphanumeric($a, $b) }
                  keys %mapping_hash;

printf MAPPING_FILE "%s\t%s\n", 'row ID', 'NameIsoCorrectoR';
foreach $name_new (@name_new_array)
{
    $row_id = $mapping_hash{$name_new};
    
    printf MAPPING_FILE "%s\t%s\n", $row_id, $name_new;
}
close MAPPING_FILE;

close DATA;
close MOL_FILE;
close MEAS_FILE;
