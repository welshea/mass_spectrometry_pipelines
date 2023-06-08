#!/usr/bin/perl -w

# 2023-06-08:  remove Formula and SMILES columns from LipidMatch output
# 2023-06-08:  add more LipidMatch support
# 2023-05-26:  begin adding LipidMatch support
# 2023-05-26:  fix %neg_col_hash = () typo, don't think it hurt anything
# 2023-02-13:  warn if tScore headers are detected
# 2022-09-15:  handle more mixed +/- adducts
# 2022-08-04:  default to removing unobservable adducts
# 2022-08-04:  change most ion related variable names to adduct related
# 2022-08-03:  add additional known adduct charges
# 2022-03-01:  infer unknown adduct charge from known adducts
# 2022-02-28:  major changes to support LipidSearch mergeResult.txt
# 2021-11-30:  initial release


use Scalar::Util qw(looks_like_number);
use POSIX;


# known adduct charges
$known_adduct_charge_hash{'+H'}     = 'pos';
$known_adduct_charge_hash{'+NH4'}   = 'pos';
$known_adduct_charge_hash{'+Na'}    = 'pos';
$known_adduct_charge_hash{'+K'}     = 'pos';
$known_adduct_charge_hash{'+H+ACN'} = 'pos';
$known_adduct_charge_hash{'+H-H2O'} = 'pos';
$known_adduct_charge_hash{'+H-NH3'} = 'pos';

$known_adduct_charge_hash{'-H'}      = 'neg';
$known_adduct_charge_hash{'-2H'}     = 'neg';   # -2 charge
$known_adduct_charge_hash{'-CH3'}    = 'neg';
$known_adduct_charge_hash{'+Cl'}     = 'neg';
$known_adduct_charge_hash{'-H-NH3'}  = 'neg';
$known_adduct_charge_hash{'-H-CO2'}  = 'neg';
$known_adduct_charge_hash{'+HCOO'}   = 'neg';
$known_adduct_charge_hash{'+C2H3O2'} = 'neg';
$known_adduct_charge_hash{'+CH3COO'} = 'neg';


# keep only these adducts, filter out the rest
# any others usually indicate forgetting to deselect them in LipidSearch
$adduct_observable_hash{'+H'}      = 'pos';
$adduct_observable_hash{'+NH4'}    = 'pos';
$adduct_observable_hash{'+Na'}     = 'pos';
$adduct_observable_hash{'+H-H2O'}  = 'pos';
$adduct_observable_hash{'-H'}      = 'neg';
$adduct_observable_hash{'+HCOO'}   = 'neg';
$adduct_observable_hash{'-2H'}     = 'neg';
$adduct_observable_hash{'+C2H3O2'} = 'neg';
$adduct_observable_hash{'+CH3COO'} = 'neg';

# LipidSearch headers
$keep_header_hash{'Rej.'} = 1;
$keep_header_hash{'LipidIon'} = 1;
$keep_header_hash{'LipidGroup'} = 1;
$keep_header_hash{'Class'} = 1;
$keep_header_hash{'FattyAcid'} = 1;
$keep_header_hash{'FA1'} = 1;
$keep_header_hash{'FA2'} = 1;
$keep_header_hash{'FA3'} = 1;
$keep_header_hash{'FA4'} = 1;
$keep_header_hash{'CalcMz'} = 1;
$keep_header_hash{'IonFormula'} = 1;

# LipidMatch headers
$keep_header_hash{'Score'} = 1;
$keep_header_hash{'SeriesType_Identifier'} = 1;
$keep_header_hash{'Name_or_Class'} = 1;
#$keep_header_hash{'Formula'} = 1;       # useless and wrong, always CH2
#$keep_header_hash{'SMILES'} = 1;        # useless and wrong, always CC
$keep_header_hash{'m/z'} = 1;
$keep_header_hash{'Retention Time'} = 1;
$keep_header_hash{'Adduct'} = 1;
$keep_header_hash{'Unique'} = 1;
$keep_header_hash{'Score_Description'} = 1;
$keep_header_hash{'Needs_Validation'} = 1;
$keep_header_hash{'row.ID'} = 1;
$keep_header_hash{'row.number.of.detected.peaks'} = 1;
$keep_header_hash{'row ID'} = 1;
$keep_header_hash{'row number of detected peaks'} = 1;

# more LipidMatch headers
$keep_header_hash{'Class_At_Max_Intensity'} = 1;
$keep_header_hash{'Adduct_At_Max_Intensity'} = 1;
$keep_header_hash{'Only_One_Class'} = 1;
$keep_header_hash{'Num_Frags'} = 1;
$keep_header_hash{'Molecular'} = 1;
$keep_header_hash{'Adducts.Confirmed.by.MS.MS'} = 1;
$keep_header_hash{'nominal mass'} = 1;
$keep_header_hash{'mass defect'} = 1;
$keep_header_hash{'[CH2]n_Kendrick_mass'} = 1;
$keep_header_hash{'[CH2]n_nominal_Kendrick_mass'} = 1;
$keep_header_hash{'[CH2]n_Kendrick_mass_defect'} = 1;
$keep_header_hash{'Number_of_F_Fragments'} = 1;


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


sub csv2tsv_not_excel
{
    my $line = $_[0];
    my $i;
    my $n;
    my @temp_array;

    # placeholder strings unlikely to ever be encountered normally
    #    my $tab = '___TaBChaR___';
    #    my $dq  = '___DqUOteS___';
    #
    # https://stackoverflow.com/questions/8695118/whats-the-file-group-record-unit-separator-control-characters-and-its-usage#:~:text=30%20%E2%80%93%20RS%20%E2%80%93%20Record%20separator%20Within,units%20in%20the%20ASCII%20definition.
    # interestingly, Microsoft Word appears to use 1E and 1F as
    #  non-breaking and optional hyphen
    #
    # http://jkorpela.fi/chars/c0.html
    #
    # I'm going to use \x1A (substitute) for tab, since it is
    #  "used in the place of a character that has been found to be invalid
    #   or in error. SUB is intended to be introduced by automatic means",
    #   and that is exactly how this function uses it.
    #
    # I'll use \x1D for "", since Word may use 1E and 1F internally,
    #  and who knows if they may ever accidentally show up in exported files,
    #  plus "group separator" seems somewhat appropriate, given that regular
    #  double-quotes are used for "grouping".
    #
    my $tab = "\x1A";    # (substitute)      need single char for split regex
    my $dq  = "\x1D";    # (group separator) need single char for split regex
    
    # remove null characters, since I've only seen them randomly introduced
    #  during NTFS glitches; "Null characters may be inserted into or removed
    #  from a stream of data without affecting the information content of that
    #  stream."
    $line =~ s/\x00//g;
    
    # remove UTF8 byte order mark, since it corrupts the first field
    # also remove some weird corruption of the UTF8 byte order mark (??)
    #
    # $line =~ s/^\xEF\xBB\xBF//;      # actual BOM
    # $line =~ s/^\xEF\x3E\x3E\xBF//;  # corrupted BOM I have seen in the wild
    $line =~ s/^(\xEF\xBB\xBF|\xEF\x3E\x3E\xBF)//;

    # replace any (incredibly unlikely) instances of $dq with $tab
    $line =~ s/$dq/$tab/g;
    
    # replace embedded tabs with placeholder string, to deal with better later
    $line =~ s/\t/$tab/g;
    
    # HACK -- handle internal \r and \n the same way we handle tabs
    $line =~ s/[\r\n]+(?!$)/$tab/g;
    
    # further escape ""
    $line =~ s/""/$dq/g;

    # only apply slow tab expansion to lines still containing quotes
    if ($line =~ /"/)
    {
        # convert commas only if they are not within double quotes
        # incrementing $i within array access for minor speed increase
        #   requires initializing things to -2
        @temp_array = split /((?<![^, $tab$dq])"[^\t"]+"(?![^, $tab$dq\r\n]))/, $line;
        $n = @temp_array - 2;
        for ($i = -2; $i < $n;)
        {
            $temp_array[$i += 2] =~ tr/,/\t/;
        }
        $line = join '', @temp_array;
        
        # slightly faster than split loop on rows with many quoted fields,
        #  but *much* slower on lines containing very few quoted fields
        # use split loop instead
        #
        # /e to evaluate code to handle different capture cases correctly
        #$line =~ s/(,?)((?<![^, $tab$dq])"[^\t"]+"(?![^, $tab$dq\r\n]))|(,)/defined($3) ? "\t" : ((defined($1) && $1 ne '') ? "\t$2" : $2)/ge;
    }
    else
    {
        $line =~ tr/,/\t/;
    }

    # unescape ""
    $line =~ s/$dq/""/g;

    # finish dealing with embedded tabs
    # remove tabs entirely, preserving surrounding whitespace
    $line =~ s/(\s|^)($tab)+/$1/g;
    $line =~ s/($tab)+(\s|$)/$2/g;
    # replace remaining tabs with spaces so text doesn't abutt together
    $line =~ s/($tab)+/ /g;

    # Special case "" in a field by itself.
    #
    # This generally results from lazily-coded csv writers that enclose
    #  every single field in "", even empty fields, whether they need them
    #  or not.
    #
    # \K requires Perl >= v5.10.0 (2007-12-18)
    #   (?: *)\K is faster than replacing ( *) with $1
    $line =~ s/(?:(?<=\t)|^)(?: *)\K""( *)(?=\t)/$1/g;  # start|tabs "" tabs
    $line =~    s/(?<=\t)(?: *)\K""( *)(?=\t|$)/$1/g;   # tabs  "" tabs|end
    $line =~          s/^(?: *)\K""( *)$/$1/g;          # start "" end

    # strip enclosing double-quotes, preserve leading/trailing spaces
    #
    # \K requires Perl >= v5.10.0 (2007-12-18)
    #   (?: *)\K is faster than replacing ( *) with $1
    #
    #$line =~ s/(?:(?<=\t)|^)(?: *)\K"([^\t]+)"( *)(?=\t|[\r\n]*$)/$1$2/g;

    # remove enclosing spaces, to support space-justified quoted fields
    #   ( *) might be faster without () grouping, but left in for clarity
    $line =~ s/(?:(?<=\t)|^)( *)"([^\t]+)"( *)(?=\t|[\r\n]*$)/$2/g;

    # unescape escaped double-quotes
    $line =~ s/""/"/g;

    
    return $line;
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



%keep_col_hash = ();
$num_files     = 0;
$filename_pos  = 'lipidomics_split_pos.txt';
$filename_neg  = 'lipidomics_split_neg.txt';


$all_adducts_flag     = 0;    # default to removing adducts we don't want
$tscore_detected_flag = 0;    # only present in incorrect exports?
$lipidmatch_flag      = 0;    # is this a LipidMatch CombinedIDed_FIN file?


for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        if ($field eq '--all-adducts')
        {
            $all_adducts_flag = 1;
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
            $filename_pos = $field;
            $num_files++;
        }
        elsif ($num_files == 2)
        {
            $filename_neg = $field;
            $num_files++;
        }
    }
}


if (!defined($filename) || $syntax_error_flag)
{
    printf STDERR "Usage: split_lipidomics_into_pos_neg.pl [options] tab_delimited.txt [outfile_pos outfile_neg]\n";
    printf STDERR "\n";
    printf STDERR "  Options:\n";
    printf STDERR "    --all-adducts     keep all adducts, including unobservable\n";
    exit(1);
}


open INFILE,      "$filename"      or die "can't open $filename\n";
open OUTFILE_POS, ">$filename_pos" or die "can't open output $filename_pos\n";
open OUTFILE_NEG, ">$filename_neg" or die "can't open output $filename_neg\n";


# use filename extension to determine whether it is tsv or csv
$delimiter_type = 'tsv';
if ($filename =~ /\.csv$/i)
{
    $delimiter_type = 'csv';
}


# skip down to first line that has anything on it
# lipidomics data has this issues sometimes
%sample_rename_hash = ();
while($line=<INFILE>)
{
    $line =~ s/[\r\n]+$//;

    # skip comment lines
    # comment may contain sample renaming information, so check for it
    if ($line =~ /^#/)
    {
        if ($line =~ /^#\[([^\]]+)\]:(.*?)\.raw\s*$/)
        {
            $sample_old = $1;

            ## comment out for now, we'll rename later in the strip script
            #$sample_new = $2;
            #$sample_rename_hash{$sample_old} = $sample_new;
            #printf STDERR "Renaming:\t%s\t%s\n", $sample_old, $sample_new;

            # strip all the extra blank tabs/whitespace from end of line
            $line =~ s/\s+$//;
            
            $rename_sample_line_hash{$sample_old} = $line;
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
if ($delimiter_type eq 'csv') { $line = csv2tsv_not_excel($line); }
$line =~ s/[\r\n]+//g;
$line =~ s/\"//g;
@array = split /\t/, $line;    # skip empty headers at and

# clean up fields, detect file type
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    $array[$i] =~ s/\s+/ /g;

    # clean up sample names
    $array[$i] =~ s/^_+//;
    $array[$i] =~ s/_+$//;

    $field = $array[$i];

    # LipidMatch CombinedIDed_FIN.csv
    if ($field =~ /SeriesType_Identifier/)
    {
        $lipidmatch_flag = 1;
    }
}


for ($i = 0; $i < @array; $i++)
{
    # handle LipidMatch strangeness
    if ($lipidmatch_flag)
    {
        # strip pos/neg from sample names,
        # since only one is given, but not the other
        #
        # FIXME -- does not currently handle p[ and n[
        # FIXME -- currently not sufficiently pos/neg paranoid

        $field = $array[$i];
        $field_orig = $field;
        
        # strip pos from sample name
        if (!($field =~ s/([^A-Za-z0-9])pos[0-9]{0,1}(?:\]*)$/$1/i ||
              $field =~ s/(^|[^\]\)\}A-Za-z0-9]+)pos([^\[\(\{A-Za-z0-9]+|$)/$2/i))
        {
            $field =~ s/^pos([^A-Za-z0-9])/$1/i;
        }

        # strip neg from sample name if no pos was stripped
        if ($field eq $field_orig)
        {
            if (!($field =~ s/([^A-Za-z0-9])neg[0-9]{0,1}(?:\]*)$/$1/i ||
                 $field =~ s/(^|[^\]\)\}A-Za-z0-9]+)neg([^\[\(\{A-Za-z0-9]+|$)/$2/i))
            {
                $field =~ s/^neg([^A-Za-z0-9])/$1/i;
            }
        }

        # clean up underscores, etc.
        $field =~ s/[_ ]+/_/g;
        $field =~ s/\-+/\-/g;
        $field =~ s/^[_ -]//;
        $field =~ s/[_ -]$//;
        
        $array[$i] = $field;


        # replace . with spaces
        $array[$i] =~ s/\./ /g;
    }
    

    $field = $array[$i];
    
    $header_col_hash{$field} = $i;
    $header_col_array[$i] = $field;

    # flag suspicious columns that may indicate wrong export options
    if ($field =~ /^tScore/)
    {
       $tscore_detected_flag = 1;
    }
    # LipidMatch CombinedIDed_FIN.csv
    if ($field =~ /SeriesType_Identifier/)
    {
        $lipidmatch_flag = 1;
    }
}


# Excel gets very upset if first field is "ID", change it
if ($header_col_array[0] =~ /^id$/i)
{
    $header_col_array[0] = 'Index';
    $header_col_hash{'Index'} = 0;
}


# rename samples from comment lines in lipidomics data
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];

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


# conform sample columns
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];

    # strip mzXML from sample names
    $field =~ s/[. ]mzX?ML( Peak \S+)$/$1/i;
    $field =~ s/[. ]cdf( Peak \S+)$/$1/i;

    if ($field =~ / Peak (\S+)$/i)
    {

        $second_word = $1;
        $second_word = lc $second_word;
        $field =~ s/[ ._]+Peak (\S+)$/ Peak $second_word/i;
    }
    $field =~ s/[ ._]+$//;

    $header_col_array[$col] = $field;
}
$num_header_cols = @header_col_array;


# scan for peak height and peak area
$peak_height_flag = 0;
$peak_area_flag   = 0;

for ($col = 0; $col < $num_header_cols; $col++)
{
    $field = $header_col_array[$col];

    if ($field =~ / Peak height$/i)
    {
        $peak_height_flag = 1;
        $keep_col_hash{$col} = 1;

        next;
    }
    if ($field =~ / Peak area$/i)
    {
        $peak_area_flag = 1;
        $keep_col_hash{$col} = 1;

        next;
    }
    
    # lipidomics
    if ($field =~ /^Area[,[]/i)
    {
        $peak_area_flag = 1;
        $keep_col_hash{$col} = 1;

        next;
    }
    if ($field =~ /^Height[,[]/i)
    {
        $peak_height_flag = 1;
        $keep_col_hash{$col} = 1;

        next;
    }
    if ($field =~ /^Rt[,[]/i)
    {
        $keep_col_hash{$col} = 1;
    }
    if (defined($keep_header_hash{$field}))
    {
        $keep_col_hash{$col} = 1;
    }
}


# uh oh, no data columns found
# conform them so that Area or Height start at the beginning or end
if ($peak_height_flag == 0 && $peak_area_flag == 0)
{
    printf STDERR "WARNING -- non-standard Height/Area nomenclature, conforming with best guess\n";

    for ($col = 0; $col < $num_header_cols; $col++)
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

            $field .= ' Peak area';

            $header_col_array[$col] = $field;
            $peak_area_flag = 1;
        }
    }
}


$name_col = $header_col_hash{'LipidIon'};

# LipidMatch
if (!defined($name_col))
{
    $name_col = $header_col_hash{'Molecular'};
}
$adduct_col = $header_col_hash{'Adduct'};


$first_abundance_col = 9E99;
for ($col = 0; $col < $num_header_cols; $col++)
{
    $field = $header_col_array[$col];

    if ($field =~ / Peak height$/i ||
        $field =~ / Peak area$/i ||
        $field =~ /^Area[,[]/i ||
        $field =~ /^Height[,[]/i)
    {
        $sample_col_hash{$col} = 1;
        $keep_col_hash{$col} = 1;

        if ($col < $first_abundance_col)
        {
            $first_abundance_col = $col;
        }
    }
    else
    {
        if ($field =~ /\S/)
        {
#            $metadata_col_hash{$col} = 1;
        }
    }
}


printf STDERR "First abundance col: %s\n", $first_abundance_col;


@sample_col_array = sort {$a<=>$b} keys %sample_col_hash;
#$num_sample_cols  = @sample_col_array;


if (!defined($name_col))
{
    printf STDERR "ABORT -- can't find 'LipidIon' column\n";
    exit(1);
}


#$rowid_col = $header_col_hash{'row ID'};
#if (!defined($rowid_col))
#{
#    printf STDERR "No row identifier column found, creating arbitrary identifiers\n";
#}


#@sample_col_array   = sort {$a <=> $b} keys %sample_col_hash;
#@metadata_col_array = sort {$a <=> $b} keys %metadata_col_hash;


# read in file, figure out which rows/columns are pos/neg later
$row = 0;
while(defined($line=<INFILE>))
{
    if ($delimiter_type eq 'csv') { $line = csv2tsv_not_excel($line); }

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
    
    # store data
    for ($col = 0; $col < $num_header_cols; $col++)
    {
        $field = $array[$col];
        
        if (!defined($field))
        {
            $field = '';
        }
        
        # set missing or non-numeric abundance data to zero
        if (defined($sample_col_hash{$col}) &&
            !is_number($field))
        {
            $field = '0';
        }
        
        $row_data_array[$row][$col] = $field;
    }

    $row++;
}
$num_rows = $row;


foreach $col (@sample_col_array)
{
    $col_counts_hash{$col}{pos} = 0;
    $col_counts_hash{$col}{neg} = 0;
}


# flag pos/neg rows
# count pos/neg for each sample
for ($row = 0; $row < $num_rows; $row++)
{
    $lipid_ion = $row_data_array[$row][$name_col];

    $pos_neg    = '';
    $adduct_str = '';


    if (defined($adduct_col))
    {
        $adduct_str = $row_data_array[$row][$adduct_col];
    }
    elsif ($lipid_ion =~ /([+-][A-Za-z0-9+-]+)$/)
    {
        $adduct_str = $1;
    }
    
    # strip []charge from adducts, since LipidSearch doesn't have them
    $adduct_str_orig = $adduct_str;
    $adduct_str      =~ s/\[M(.*)\][0-9]*[+-]/$1/;


    # skip adducts that shouldn't be observed
    if ($all_adducts_flag == 0 &&
        !defined($adduct_observable_hash{$adduct_str}))
    {
        $filtered_adduct_hash{$adduct_str} = 1;

        next;
    }

    # look up known ion charges
    if (defined($known_adduct_charge_hash{$adduct_str}))
    {
        $pos_neg = $known_adduct_charge_hash{$adduct_str};
        $pos_neg_row_hash{$row} = $pos_neg;
    }
    
    # unknown adduct charge, see if the charge is in the file
    if ($pos_neg eq '')
    {
        $charge = $adduct_str_orig;
        $charge =~ s/\[M(.*)\][0-9]*([+-])/$2/;
        
        if ($charge eq '+')
        {
            $pos_neg                = 'pos';
            $pos_neg_row_hash{$row} = $pos_neg;
        }
        if ($charge eq '-')
        {
            $pos_neg                = 'neg';
            $pos_neg_row_hash{$row} = $pos_neg;
        }
    }
    
    # count pos/neg for each sample
    # so that we can check for un-merged pos/neg columns later
    if ($pos_neg eq 'pos' || $pos_neg eq 'neg')
    {
        foreach $col (@sample_col_array)
        {
            $field = $row_data_array[$row][$col];
            
            if (is_number($field) && $field != 0)
            {
                $col_counts_hash{$col}{$pos_neg} += 1;
            }
        }
    }
    else
    {
        $unknown_row_charge_hash{$row} = $adduct_str;
    }
}


%pos_col_hash = ();
%neg_col_hash = ();
foreach $col (@sample_col_array)
{
    $num_pos = $col_counts_hash{$col}{pos};
    $num_neg = $col_counts_hash{$col}{neg};
    
    if ($num_pos && $num_neg == 0)
    {
        $pos_col_hash{$col} = 1;
    }
    if ($num_neg && $num_pos == 0)
    {
        $neg_col_hash{$col} = 1;
    }
}


# count pos/neg samples for all unknown adduct charges
@unknown_row_charge_array = sort {$a<=>$b} keys %unknown_row_charge_hash;
%unknown_adduct_hash = ();
foreach $row (@unknown_row_charge_array)
{
    $adduct_str = $unknown_row_charge_hash{$row};
    
    if (!defined($unknown_adduct_hash{$adduct_str}))
    {
        $unknown_adduct_hash{$adduct_str}{pos} = 0;
        $unknown_adduct_hash{$adduct_str}{neg} = 0;
    }

    foreach $col (@sample_col_array)
    {
        $field = $row_data_array[$row][$col];
            
        if (is_number($field) && $field != 0)
        {
            if (defined($pos_col_hash{$col}))
            {
                $unknown_adduct_hash{$adduct_str}{pos} += 1;
            }
            elsif (defined($neg_col_hash{$col}))
            {
                $unknown_adduct_hash{$adduct_str}{neg} += 1;
            }
        }
    }
}

# infer unknown adduct charges from pos/neg sample patterns
foreach $adduct_str (sort keys %unknown_adduct_hash)
{
    $count_pos = $unknown_adduct_hash{$adduct_str}{pos};
    $count_neg = $unknown_adduct_hash{$adduct_str}{neg};
    
    if ($count_pos && $count_neg == 0)
    {
        $known_adduct_charge_hash{$adduct_str} = 'pos';

        printf STDERR "Inferring unknown adduct charge:\t%s\t%s\n",
            $adduct_str, 'pos';
    }
    elsif ($count_neg && $count_pos == 0)
    {
        $known_adduct_charge_hash{$adduct_str} = 'neg';

        printf STDERR "Inferring unknown adduct charge:\t%s\t%s\n",
            $adduct_str, 'neg';
    }
    else
    {
        printf STDERR "WARNING -- skipping misbehaving unknown adduct:\t%s\n",
            $adduct_str;
    }
}

# assign the newly inferred charges to remaining rows
foreach $row (@unknown_row_charge_array)
{
    $adduct_str = $unknown_row_charge_hash{$row};
    
    $pos_neg = $known_adduct_charge_hash{$adduct_str};
    
    if (defined($pos_neg))
    {
        $pos_neg_row_hash{$row} = $pos_neg;
    }
}


@pos_col_array = sort {$a <=> $b} keys %pos_col_hash;
@neg_col_array = sort {$a <=> $b} keys %neg_col_hash;

#$num_pos_cols = @pos_col_array;
#$num_neg_cols = @neg_col_array;


# scan through positive and negative columns
# extract out lipidomics sample names so that we can pull in matching Rt
foreach $col (@pos_col_array)
{
    $header = $header_col_array[$col];
    
    if ($header =~ /^[^[]+[,[]([^[]+)\]/)
    {
        $sample = $1;
        
        $pos_sample_hash{$sample} = 1;
    }
}
foreach $col (@neg_col_array)
{
    $header = $header_col_array[$col];
    
    if ($header =~ /^[^[]+[,[]([^[]+)\]/)
    {
        $sample = $1;
        
        $neg_sample_hash{$sample} = 1;
    }
}
for ($col = 0; $col < @header_col_array; $col++)
{
    if (!defined($keep_col_hash{$col})) { next; }

    $header = $header_col_array[$col];

    if ($header =~ /^[^[]+[,[]([^[]+)\]/)
    {
        $sample = $1;
        
        if (defined($pos_sample_hash{$sample}))
        {
            $pos_col_hash{$col} = 1;
        }
        if (defined($neg_sample_hash{$sample}))
        {
            $neg_col_hash{$col} = 1;
        }
    }
}


# print sample rename comment lines if we haven't renamed anything
@temp_array = sort keys %sample_rename_hash;
if (@temp_array == 0)
{
    $temp_pos_flag = 0;
    $temp_neg_flag = 0;

    foreach $sample (sort keys %pos_sample_hash)
    {
        $line = $rename_sample_line_hash{$sample};
        if (defined($line))
        {
            print OUTFILE_POS "$line\n";
            
            $temp_pos_flag = 1;
        }
    }
    foreach $sample (sort keys %neg_sample_hash)
    {
        $line = $rename_sample_line_hash{$sample};
        if (defined($line))
        {
            print OUTFILE_NEG "$line\n";

            $temp_neg_flag = 1;
        }
    }
    
    # print blank lines to separate comments from true header
    if ($temp_pos_flag)
    {
        print OUTFILE_POS "\n";
    }
    if ($temp_neg_flag)
    {
        print OUTFILE_NEG "\n";
    }
}


# print new header lines
$header_pos = '';
$header_neg = '';
for ($col = 0; $col < @header_col_array; $col++)
{
    if (!defined($keep_col_hash{$col})) { next; }

    $header = $header_col_array[$col];

    if (defined($keep_header_hash{$header}) ||
        defined($pos_col_hash{$col}) ||
        ($lipidmatch_flag && defined($keep_col_hash{$col})))
    {
        if ($header_pos ne '') { $header_pos .= "\t"; }
        $header_pos .= $header;
    }

    if (defined($keep_header_hash{$header}) ||
        defined($neg_col_hash{$col}) ||
        ($lipidmatch_flag && defined($keep_col_hash{$col})))
    {
        if ($header_neg ne '') { $header_neg .= "\t"; }
        $header_neg .= $header;
    }
}
print OUTFILE_POS "$header_pos\n";
print OUTFILE_NEG "$header_neg\n";



# output data, split into pos and neg files
$line_pos = '';
$line_neg = '';
for ($row = 0; $row < $num_rows; $row++)
{
    $line_pos = '';
    $line_neg = '';
    for ($col = 0; $col < @header_col_array; $col++)
    {
        if (!defined($keep_col_hash{$col})) { next; }

        $header = $header_col_array[$col];
        $field  = $row_data_array[$row][$col];
    
        if (defined($keep_header_hash{$header}) ||
            defined($pos_col_hash{$col}) ||
            ($lipidmatch_flag && defined($keep_col_hash{$col})))
        {
            if ($line_pos ne '') { $line_pos .= "\t"; }
            $line_pos .= $field;
        }

        if (defined($keep_header_hash{$header}) ||
            defined($neg_col_hash{$col}) ||
            ($lipidmatch_flag && defined($keep_col_hash{$col})))
        {
            if ($line_neg ne '') { $line_neg .= "\t"; }
            $line_neg .= $field;
        }
    }
    

    $lipid_ion  = $row_data_array[$row][$name_col];
    $adduct_str = '';


    if (defined($adduct_col))
    {
        $adduct_str = $row_data_array[$row][$adduct_col];
    }
    elsif ($lipid_ion =~ /([+-][A-Za-z0-9+-]+)$/)
    {
        $adduct_str = $1;
    }

    # strip []charge from adducts, since LipidSearch doesn't have them
    $adduct_str =~ s/\[M(.*)\][0-9]*[+-]/$1/;


    # skip adducts that shouldn't be observed
    if ($all_adducts_flag == 0 &&
        !defined($adduct_observable_hash{$adduct_str}))
    {
        $filtered_adduct_hash{$adduct_str} = 1;

        next;
    }

    
    # output positive ions to positive file
    # output negative ions to negative file
    # skip uncertain charge lines
    #
    $pos_neg = $pos_neg_row_hash{$row};
    if (!defined($pos_neg))
    {
        $pos_neg = '';
    }
    if ($pos_neg eq 'pos')
    {
        print OUTFILE_POS "$line_pos\n";
    }
    elsif ($pos_neg eq 'neg')
    {
        print OUTFILE_NEG "$line_neg\n";
    }
    else
    {
        printf STDERR "WARNING -- skipping unknown adduct row:\t%s\n", $lipid_ion;
    }
}

close OUTFILE_POS;
close OUTFILE_NEG;


@filtered_adduct_array = sort keys %filtered_adduct_hash;
foreach $adduct (@filtered_adduct_array)
{
    printf STDERR "WARNING -- removed %s adducts, not in observeable adduct list\n",
        $adduct;
}

if ($tscore_detected_flag)
{
    printf STDERR "WARNING -- tScore column headers detected; check LipidSearch export options\n";
}
