#!/usr/bin/perl -w

# 2022-08-05:  flag main ions from LipidSearch summary file
# 2022-08-05:  denote non-canonical main ions
# 2022-08-03:  rename SuperClass to Class
# 2022-03-10:  annotate lipid class with more descriptive categories
# 2022-03-10:  edit usage statement to use current program name
# 2022-03-10:  surround adduct column in [] to keep Excel happy
# 2022-03-10:  change singleton/primary,secondary to main,blank


# run this on already normalized, logged data
# assume we've only kept height or area, not both
# expects LipidSearch formatted data


use Scalar::Util qw(looks_like_number);
use POSIX;
use File::Spec;
use File::Basename;

# this correctly doesn't follow symlinked scripts, so you get the right path
# various other methods return the destination linked path instead,
#  which we don't want here
$script_path   = dirname(File::Spec->rel2abs(__FILE__));


# number of standard deviations outwards from mean
#
# Kissel #3645 PA(18:1_18:1)-H
#  2sd isn't enough to count them as the same
$num_sd = 3;


# for each lipid class, we should see the following as a main ion
# if we don't, something may be wrong

# MainGrade A, B, C
$optimal_main_hash{'Cer'}{'+H'}        = 'ABC';
$optimal_main_hash{'Cer'}{'+HCOO'}     = 'ABC';
$optimal_main_hash{'ChE'}{'+NH4'}      = 'ABC';
$optimal_main_hash{'D7ChE'}{'+NH4'}    = 'ABC';
$optimal_main_hash{'FA'}{'-H'}         = 'ABC';
$optimal_main_hash{'Hex1Cer'}{'+H'}    = 'ABC';
$optimal_main_hash{'Hex1Cer'}{'+HCOO'} = 'ABC';
$optimal_main_hash{'Hex2Cer'}{'+H'}    = 'ABC';
$optimal_main_hash{'Hex2Cer'}{'+HCOO'} = 'ABC';
$optimal_main_hash{'Hex3Cer'}{'+H'}    = 'ABC';
$optimal_main_hash{'Hex3Cer'}{'+HCOO'} = 'ABC';
$optimal_main_hash{'LSM'}{'+H'}        = 'ABC';
$optimal_main_hash{'MG'}{'+H'}         = 'ABC';
$optimal_main_hash{'SM'}{'+H'}         = 'ABC';
$optimal_main_hash{'SPH'}{'+H'}        = 'ABC';
$optimal_main_hash{'SPH'}{'+H-H2O'}    = 'ABC';

# MainGrade A, B
$optimal_main_hash{'AcCa'}{'+H'}       = 'AB';
$optimal_main_hash{'CL'}{'-H'}         = 'AB';
$optimal_main_hash{'Co'}{'+H'}         = 'AB';
$optimal_main_hash{'Co'}{'+NH4'}       = 'AB';
$optimal_main_hash{'DG'}{'+NH4'}       = 'AB';
$optimal_main_hash{'DG'}{'+Na'}        = 'AB';
$optimal_main_hash{'DLCL'}{'-H'}       = 'AB';
$optimal_main_hash{'LPA'}{'-H'}        = 'AB';
$optimal_main_hash{'LPC'}{'+H'}        = 'AB';
$optimal_main_hash{'LPC'}{'+HCOO'}     = 'AB';
$optimal_main_hash{'LPE'}{'+H'}        = 'AB';
$optimal_main_hash{'LPE'}{'-H'}        = 'AB';
$optimal_main_hash{'LPG'}{'+NH4'}      = 'AB';
$optimal_main_hash{'LPG'}{'-H'}        = 'AB';
$optimal_main_hash{'LPI'}{'+NH4'}      = 'AB';
$optimal_main_hash{'LPI'}{'-H'}        = 'AB';
$optimal_main_hash{'LPS'}{'+H'}        = 'AB';
$optimal_main_hash{'LPS'}{'-H'}        = 'AB';
$optimal_main_hash{'MLCL'}{'-H'}       = 'AB';
$optimal_main_hash{'PA'}{'-H'}         = 'AB';
$optimal_main_hash{'PC'}{'+H'}         = 'AB';
$optimal_main_hash{'PC'}{'+HCOO'}      = 'AB';
$optimal_main_hash{'PE'}{'+H'}         = 'AB';
$optimal_main_hash{'PE'}{'-H'}         = 'AB';
$optimal_main_hash{'PG'}{'+NH4'}       = 'AB';
$optimal_main_hash{'PG'}{'-H'}         = 'AB';
$optimal_main_hash{'PI'}{'+NH4'}       = 'AB';
$optimal_main_hash{'PI'}{'-H'}         = 'AB';
$optimal_main_hash{'PS'}{'+H'}         = 'AB';
$optimal_main_hash{'PS'}{'-H'}         = 'AB';
$optimal_main_hash{'TG'}{'+NH4'}       = 'AB';



sub cmp_row_rt
{
    my $rt_a     = $row_rt_avg_array[$a];
    my $rt_b     = $row_rt_avg_array[$b];
    my $adduct_a = $row_adduct_array[$a];
    my $adduct_b = $row_adduct_array[$b];
    my $ion_a    = $row_fattyacid_ion_array[$a];
    my $ion_b    = $row_fattyacid_ion_array[$b];

    my $sd_a     = $row_rt_sd_array[$a];
    my $sd_b     = $row_rt_sd_array[$b];
    my $ub_a     = $rt_a + $num_sd * $sd_a;
    my $ub_b     = $rt_b + $num_sd * $sd_b;
    my $lb_a     = $rt_a - $num_sd * $sd_a;
    my $lb_b     = $rt_b - $num_sd * $sd_b;

    if ($rt_a < $rt_b) { return -1; }
    if ($rt_a > $rt_b) { return  1; }

    if ($ub_a < $ub_b) { return -1; }
    if ($ub_a > $ub_b) { return  1; }

    if ($lb_a < $lb_b) { return -1; }
    if ($lb_a > $lb_b) { return  1; }
    
    
    if ($adduct_a lt $adduct_b) { return -1; }
    if ($adduct_a gt $adduct_b) { return  1; }

    if ($ion_a lt $ion_b) { return -1; }
    if ($ion_a gt $ion_b) { return  1; }
    
    return $a<=>$b;
}


sub cmp_row_main
{
    my $rt_a       = $row_rt_avg_array[$a];
    my $rt_b       = $row_rt_avg_array[$b];
    my $adduct_a   = $row_adduct_array[$a];
    my $adduct_b   = $row_adduct_array[$b];
    my $ion_a      = $row_fattyacid_ion_array[$a];
    my $ion_b      = $row_fattyacid_ion_array[$b];
    my $data_avg_a = $row_data_avg_array[$a];
    my $data_avg_b = $row_data_avg_array[$b];
    my $data_sum_a = $row_data_sum_array[$a];
    my $data_sum_b = $row_data_sum_array[$b];
    my $count_a    = $row_data_count_array[$a];
    my $count_b    = $row_data_count_array[$b];

    # rank by summed abundance
    # takes both mean abundance and #present into account
    if ($data_sum_a > $data_sum_b) { return -1; }
    if ($data_sum_a < $data_sum_b) { return  1; }

    # prefer rows with less missing data
    if ($count_a > $count_b) { return -1; }
    if ($count_a < $count_b) { return  1; }
    
    # if we're still equal by now, the means are probably equal anyways...
    if ($data_avg_a > $data_avg_b) { return -1; }
    if ($data_avg_a < $data_avg_b) { return  1; }

    if ($rt_a < $rt_b) { return -1; }
    if ($rt_a > $rt_b) { return  1; }
    
    if ($adduct_a lt $adduct_b) { return -1; }
    if ($adduct_a gt $adduct_b) { return  1; }

    if ($ion_a lt $ion_b) { return -1; }
    if ($ion_a gt $ion_b) { return  1; }
    
    return $a<=>$b;
}


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


# for now, hard code the annotation file name
sub read_in_lipid_class_annotation
{
    $annotation_filename = 'lipid_class_annotation.txt';
    $full_path = $script_path . '/' . $annotation_filename;

    %header_col_hash  = ();
    @header_col_array = ();
    
    open ANNOTATION, "$full_path" or die "ABORT -- can't open $full_path\n";
    
    $line = <ANNOTATION>;
    $line =~ s/[\r\n]+$//;
    @array = split /\t/, $line;
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        
        $changed_flag = 1;
        while ($changed_flag)
        {
            $changed_flag = ($array[$i] =~ s/^\"(.*)\"$/$1/);
        }
        
        if ($array[$i] =~ /\S/)
        {
            $header_col_array[$i] = $array[$i];
            $header_col_hash{$array[$i]} = $i;
        }
    }
    
    $abbrev_col   = $header_col_hash{Abbreviation};
    $supclass_col = $header_col_hash{Class};
    $subclass_col = $header_col_hash{SubClass};
    
    # backwards compatability with older SuperClass terminology
    if (!defined($supclass_col))
    {
        $supclass_col = $header_col_hash{SuperClass};
    }
    
    if (!defined($abbrev_col))
    {
        printf "ABORT -- Abbreviation column not found in %s\n", $full_path;
        exit(2);
    }
    if (!defined($supclass_col))
    {
        printf "ABORT -- Class column not found in %s\n", $full_path;
        exit(2);
    }
    if (!defined($subclass_col))
    {
        printf "ABORT -- SubClass column not found in %s\n", $full_path;
        exit(2);
    }
    
    while(defined($line=<ANNOTATION>))
    {
        $line =~ s/[\r\n]+$//;
        @array = split /\t/, $line;
        for ($i = 0; $i < @array; $i++)
        {
            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;
        
            $changed_flag = 1;
            while ($changed_flag)
            {
                $changed_flag = ($array[$i] =~ s/^\"(.*)\"$/$1/);
            }
        }
        
        $abbrev   = $array[$abbrev_col];
        $supclass = $array[$supclass_col];
        $subclass = $array[$subclass_col];
        
        if ($abbrev =~ /\S/)
        {
            $class_annotation_hash{$abbrev}{super} = $supclass;
            $class_annotation_hash{$abbrev}{sub}   = $subclass;
        }
    }
    close ANNOTATION;
}


sub read_in_lipidsearch_summary_file
{
    my $summary_filename = $_[0];

    %header_col_hash  = ();
    @header_col_array = ();
    
    open SUMMARY, "$summary_filename" or die "ABORT -- can't open $summary_filename\n";
    
    # continue reading until we get to the real header line
    $line = <SUMMARY>;
    while (!($line =~ /MainIon/ && $line =~ /BaseRt/))
    {
        $line = <SUMMARY>;
    }

    $line =~ s/[\r\n]+$//;
    
    @array = split /\t/, $line;
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        
        $changed_flag = 1;
        while ($changed_flag)
        {
            $changed_flag = ($array[$i] =~ s/^\"(.*)\"$/$1/);
        }
        
        if ($array[$i] =~ /\S/)
        {
            $header_col_array[$i] = $array[$i];
            $header_col_hash{$array[$i]} = $i;
        }
    }
    
    $lipid_molec_col = $header_col_hash{'LipidMolec'};
    $base_rt_col     = $header_col_hash{'BaseRt'};
    $main_ion_col    = $header_col_hash{'MainIon'};
    
    if (!defined($lipid_molec_col))
    {
        printf STDERR "ABORT -- LipidMolec column not found in summary file %s\n",
            $summary_filename;
        die(2);
    }
    if (!defined($base_rt_col))
    {
        printf STDERR "ABORT -- BaseRt column not found in summary file %s\n",
            $summary_filename;
        die(2);
    }
    if (!defined($main_ion_col))
    {
        printf STDERR "ABORT -- MainIon column not found in summary file %s\n",
            $summary_filename;
        die(2);
    }
    
    while(defined($line=<SUMMARY>))
    {
        $line =~ s/[\r\n]+$//;
        @array = split /\t/, $line;
        for ($i = 0; $i < @array; $i++)
        {
            $array[$i] =~ s/^\s+//;
            $array[$i] =~ s/\s+$//;
        
            $changed_flag = 1;
            while ($changed_flag)
            {
                $changed_flag = ($array[$i] =~ s/^\"(.*)\"$/$1/);
            }
        }
        
        $lipid_molec = $array[$lipid_molec_col];
        $base_rt     = $array[$base_rt_col];
        $main_ion    = $array[$main_ion_col];
        
        if ($lipid_molec =~ /[A-Za-z0-9]/ &&
            $main_ion    =~ /[A-Za-z]/ &&
            is_number($base_rt))
        {
            $main_id = $base_rt . $main_ion;
        
            $lipidsearch_summary_hash{$lipid_molec}{$main_id}{adduct} =
                $main_ion;
            $lipidsearch_summary_hash{$lipid_molec}{$main_id}{rt} =
                $base_rt;
            
            # read in a valid LipidSearch summary main ion data row
            $valid_summary_flag = 1;
        }
    }
    close SUMMARY;
}



$num_files          = 0;
$valid_summary_flag = 0;

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
            $filename = $field;
            $num_files++;
        }
        elsif ($num_files == 1)
        {
            $summary_filename = $field;
            $num_files++;
        }
    }
}


if (!defined($filename) || $syntax_error_flag)
{
    $program_name = basename($0);

    printf STDERR "Usage: $program_name iron_log2_merged.txt [lipidsearch_summary.txt]\n";
    exit(1);
}


# read in the lipid class annotation file, hardcoded file name
read_in_lipid_class_annotation();

# read in LipidSearch summary file
if (defined($summary_filename))
{
    read_in_lipidsearch_summary_file($summary_filename);
}

open INFILE,      "$filename"      or die "ABORT -- can't open $filename\n";


%header_col_hash  = ();
@header_col_array = ();

# skip down to first line that has anything on it
# lipidomics data has this issues sometimes
while($line=<INFILE>)
{
    $line =~ s/[\r\n]+$//;

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


$lipidion_col  = $header_col_hash{'LipidIon'};
$class_col     = $header_col_hash{'Class'};
$fattyacid_col = $header_col_hash{'FattyAcid'};

if (!defined($lipidion_col))
{
    printf STDERR "ABORT -- can't find 'LipidIon' column\n";
    exit(1);
}
if (!defined($class_col))
{
    printf STDERR "ABORT -- can't find 'Class' column\n";
    exit(1);
}
if (!defined($fattyacid_col))
{
    printf STDERR "ABORT -- can't find 'FattyAcid' column\n";
    exit(1);
}


# find first Rt column
# we're going to insert our new columns right before this
$first_rt_col = -42;
for ($col = 0; $col < @header_col_array; $col++)
{
    $header = $header_col_array[$col];
    
    if ($header =~ /^Rt/)
    {
        $first_rt_col = $col;
        last;
    }
}
if ($first_rt_col == -42)
{
    printf STDERR "ABORT -- no Rt columns found\n";
    exit(2);
}
printf STDERR "First Rt col: %d\n", $first_rt_col;

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
    if ($field =~ /^Rt[,[]/i)
    {
        $rt_col_hash{$col} = 1;
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

        if ($col < $first_abundance_col)
        {
            $first_abundance_col = $col;
        }
    }
    else
    {
        if ($field =~ /\S/)
        {
            $metadata_col_hash{$col} = 1;
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

            if ($field =~ /\S/)
            {
                $metadata_col_hash{$col} = 1;
            }
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
printf STDERR "First Rt col: %d\n", $first_rt_col;


@sample_col_array = sort {$a<=>$b} keys %sample_col_hash;
@rt_col_array     = sort {$a<=>$b} keys %rt_col_hash;


# read in file
$row = 0;
while(defined($line=<INFILE>))
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
    
    # store data
    for ($col = 0; $col < $num_header_cols; $col++)
    {
        $field = $array[$col];
        
        if (!defined($field))
        {
            $field = '';
        }
        
        # set missing or non-numeric abundance data to blank
        # assume that zero data is missing as well, even if input is logged
        if (defined($sample_col_hash{$col}) &&
            !is_number($field))
        {
            $field = '';
        }
        
        if (is_number($field) && $field == 0.0)
        {
            $field = '';
        }
        
        $row_data_array[$row][$col] = $field;
    }

    $row++;
}
$num_rows = $row;


# calculate row average abundance
$row = 0;
for ($row = 0; $row < $num_rows; $row++)
{
    $avg = 0;
    $n   = 0;

    foreach $col (@sample_col_array)
    {
        $value = $row_data_array[$row][$col];

        if (is_number($value) && $value != 0.0)
        {
            $avg += $value;
            $n++;
        }
    }

    $row_data_sum_array[$row]   = $avg;

    if ($n)
    {
        $avg /= $n;
    }
    
    $row_data_avg_array[$row]   = $avg;
    $row_data_count_array[$row] = $n;
}


# calculate row average Rt and standard deviation of Rt
$row = 0;
for ($row = 0; $row < $num_rows; $row++)
{
    $avg = 0;
    $n   = 0;
    $sd  = 0;

    foreach $col (@rt_col_array)
    {
        $value = $row_data_array[$row][$col];
        
        if (is_number($value) && $value != 0.0)
        {
            $avg += $value;
            $n++;
        }
    }
    if ($n)
    {
        $avg /= $n;
    }

    foreach $col (@rt_col_array)
    {
        $value = $row_data_array[$row][$col];
        
        if (is_number($value) && $value != 0.0)
        {
            $diff  = $avg - $value;
            $sd   += $diff * $diff;
        }
    }
    if ($n >= 2)
    {
        # sample standard deviation, (n - 1)
        #$sd = sqrt($sd / ($n - 1));

        # population standard deviation (n) gives better groupings here,
        # as it usually does in this sort of use case
        $sd = sqrt($sd / $n);
    }
    
    $row_rt_avg_array[$row] = $avg;
    $row_rt_sd_array[$row]  = $sd;
}


# store ions/adducts for each row and lipid
for ($row = 0; $row < $num_rows; $row++)
{
    $lipid_ion  = $row_data_array[$row][$lipidion_col];
    $fattyacid  = $row_data_array[$row][$fattyacid_col];
    $class      = $row_data_array[$row][$class_col];
    
    $adduct = '';
    if ($lipid_ion =~ /([+-][A-Za-z0-9-]+)$/)
    {
        $adduct = $1;
    }
    elsif ($lipid_ion =~ /([+-][A-Za-z0-9-]+)$/)
    {
        $adduct = $1;
    }
    

    # LipidIon can have more general nomenclature, which FattyAcid then
    # assigns to a more specific lipid using assignments from row with
    # similar Rt to narrow it down.  We need to reconstruct the true
    # lipid ion from the class, fattyacid, and adduct.
    #
    $fattyacid_base = $class . $fattyacid;
    $fattyacid_ion  = $class . $fattyacid . $adduct;


    if ($adduct eq '')
    {
        $adduct = 'missing';
    }
    
    $row_fattyacid_ion_array[$row]                 = $fattyacid_ion;
    $row_adduct_array[$row]                        = $adduct;
    $fattyacid_hash{$fattyacid_base}{$row}{adduct} = $adduct;
}


@fattyacid_base_array = sort keys %fattyacid_hash;

foreach $fattyacid_base (@fattyacid_base_array)
{
    @adduct_row_array = sort cmp_row_rt keys %{$fattyacid_hash{$fattyacid_base}};
    $num_adducts = @adduct_row_array;
    
    %temp_adduct_hash    = ();
    %temp_group_row_hash = ();
    $group               = 0;
    $rt_upper_max        = -42;
    $rt_upper_min        = -42;

    # assign each adduct within a lipid to an Rt group
    for ($i = 0; $i < $num_adducts; $i++)
    {
        $row      = $adduct_row_array[$i];
        $adduct   = $fattyacid_hash{$fattyacid_base}{$row}{adduct};
        $rt_avg   = $row_rt_avg_array[$row];
        $rt_sd    = $row_rt_sd_array[$row];
        # $data_avg = $row_data_avg_array[$row];
        # $data_sum = $row_data_sum_array[$row];
        
        $rt_min = $rt_avg - $num_sd * $rt_sd;
        $rt_max = $rt_avg + $num_sd * $rt_sd;
        
        # new group
        # test its lower bound against the group lowest upper bound so far
        # we test against the lowest upper bound, rather than highest,
        #  in order to minimize potential chaining into a too-large group
        #
        # Kissel #3645 PG(18:1_18:1)+NH4 is a difficult case
        #
        if ($rt_min <= $rt_upper_min)
        {
            # check if overlaps with next row in new group
            if ($i && $i < $num_adducts - 1)
            {
                $row_next    = $adduct_row_array[$i+1];
                $rt_avg_next = $row_rt_avg_array[$row_next];
                $rt_sd_next  = $row_rt_sd_array[$row_next];
                $rt_min_next = $rt_avg_next - $num_sd * $rt_sd_next;
                
                if ($rt_max      >= $rt_min_next &&
                    $rt_min_next >  $rt_upper_min)
                {
                    $row_prev    = $adduct_row_array[$i-1];
                    $rt_avg_prev = $row_rt_avg_array[$row_prev];
                    $dist_prev   = $rt_avg - $rt_avg_prev;
                    $dist_next   = $rt_avg_next - $rt_avg;

                    # go with the closer average
                    $choice_str = '';
                    $choice_cmp = '';
                    
                    # include current row in previous group
                    if ($dist_prev < $dist_next)
                    {
                        $choice_str = 'prev';
                        $choice_cmp = '<';

                        # we might use this for sanity checking, don't remove code yet
                        $temp_adduct_hash{$adduct} = 1;

                        # update bounds of group
                        if ($rt_max > $rt_upper_max)
                        {
                            $rt_upper_max = $rt_max;
                        }
                        if ($rt_max < $rt_upper_min)
                        {
                            $rt_upper_min = $rt_max;
                        }
                    }
                    # new group
                    else
                    {
                        $choice_str = 'next';
                        $choice_cmp = '>';
                    
                        %temp_adduct_hash = ();
                        $group++;

                        $rt_upper_max = $rt_max;
                        $rt_upper_min = $rt_max;
                    }

                    printf STDERR "RT_BETWEEN_GROUPS %s %s %.4f   add to %s: %.4f %s %.4f\n",
                        $fattyacid_base, $adduct, $rt_avg,
                        $choice_str, $dist_prev, $choice_cmp, $dist_next;
                }
            }
            # include current row in previous group
            else
            {
                # we might use this for sanity checking, don't remove code yet
                $temp_adduct_hash{$adduct} = 1;
        
                # update bounds of group
                if ($rt_max > $rt_upper_max)
                {
                    $rt_upper_max = $rt_max;
                }
                if ($rt_max < $rt_upper_min)
                {
                    $rt_upper_min = $rt_max;
                }
            }
        }
        # new group
        else
        {
            %temp_adduct_hash = ();
            $group++;
            
            $rt_upper_max = $rt_max;
            $rt_upper_min = $rt_max;
        }

        $fattyacid_hash{$fattyacid_base}{$row}{group} = $group;
        $temp_group_row_hash{$group}{$row}    = 1;
    }
    
    #$num_groups = $group;


    # identify the "main" ion within each group
    foreach $group (sort {$a<=>$b} keys %temp_group_row_hash)
    {
        @group_row_array =
            sort cmp_row_main keys %{$temp_group_row_hash{$group}};

        $row_main = $group_row_array[0];
        $fattyacid_hash{$fattyacid_base}{$row_main}{main} = 'main';
        ##$fattyacid_hash{$fattyacid_base}{$row_main}{main} = 'primary';

        #if (@group_row_array > 1)
        #{
        #    $row_main = $group_row_array[1];
        #    $fattyacid_hash{$fattyacid_base}{$row_main}{main} = 'secondary';
        #}

        #if (@group_row_array == 1)
        #{
        #    $row_main = $group_row_array[0];
        #    $fattyacid_hash{$fattyacid_base}{$row_main}{main} = 'main';
        #    ##$fattyacid_hash{$fattyacid_base}{$row_main}{main} = 'singleton';
        #}
    }
    
    # identify the main ions from LipidSearch
    if (defined($lipidsearch_summary_hash{$fattyacid_base}))
    {
        @main_id_array = sort keys %{$lipidsearch_summary_hash{$fattyacid_base}};
        
        foreach $main_id (@main_id_array)
        {
            $adduct_ls = $lipidsearch_summary_hash{$fattyacid_base}{$main_id}{adduct};
            $rt_ls     = $lipidsearch_summary_hash{$fattyacid_base}{$main_id}{rt};
            
            $best_row     = '';
            $best_rt_diff = 9E99;
            for ($i = 0; $i < $num_adducts; $i++)
            {
                $row      = $adduct_row_array[$i];
                $adduct   = $fattyacid_hash{$fattyacid_base}{$row}{adduct};
                $rt_avg   = $row_rt_avg_array[$row];
                
                if ($adduct eq $adduct_ls)
                {
                    $diff = abs($rt_avg - $rt_ls);
                    if ($diff < $best_rt_diff)
                    {
                        $best_row     = $row;
                        $best_rt_diff = $diff;
                    }
                }
            }
            
            # flag row as identified as main ion by LipidSearch
            if ($best_row ne '')
            {
                $fattyacid_hash{$fattyacid_base}{$best_row}{main_ls} = 'main';
            }
        }
    }
}


# scan main ions for potential problems
foreach $fattyacid_base (sort keys %fattyacid_hash)
{
    $lipid_class = $fattyacid_base;
    $lipid_class =~ s/\(.*//;
    
    foreach $row (sort keys %{$fattyacid_hash{$fattyacid_base}})
    {
        $main_str = $fattyacid_hash{$fattyacid_base}{$row}{main};
        
        if (defined($main_str) &&
            ($main_str eq 'main' ||
             $main_str eq 'primary' || $main_str eq 'singleton'))
        {
            $adduct = $row_adduct_array[$row];

            # only check lipids we have prior knowledge for
            if (defined($optimal_main_hash{$lipid_class}))
            {
                # main ion is not one generally observed to be highest
                if (!defined($optimal_main_hash{$lipid_class}{$adduct}))
                {
                    # overwrite with non-canonical
                    $fattyacid_hash{$fattyacid_base}{$row}{main} =
                        'main, non-canonical';
                
                    $fattyacid_ion = $row_fattyacid_ion_array[$row];
                
                    printf STDERR
                        "WARNING -- non-canonical main ion:  %s\t%s\t%s\t%s\n",
                        $fattyacid_ion, $lipid_class, $adduct, $main_str;
                }
            }
        }
    }
}


# print headers before first Rt column
for ($col = 0; $col < $first_rt_col; $col++)
{
    $header = $header_col_array[$col];

    if ($col)
    {
        printf "\t";
    }
    
    printf "%s", $header;
}

# insert our new columns
printf "\t%s", 'Lipid';
printf "\t%s", 'Class';
printf "\t%s", 'SubClass';
printf "\t%s", 'Adduct';
printf "\t%s", 'IsomerBBSR';
printf "\t%s", 'MainIonBBSR';
if ($valid_summary_flag)
{
    printf "\t%s", 'MainIonLipidSearch';
}
printf "\t%s", 'MeanRt';
printf "\t%s", 'MeanAbundance';
printf "\t%s", 'SumAbundance';

# print the rest of the headers
for ($col = $first_rt_col; $col < @header_col_array; $col++)
{
    $header = $header_col_array[$col];

    if ($col)
    {
        printf "\t";
    }
    
    printf "%s", $header;
}
printf "\n";




foreach $fattyacid_base (@fattyacid_base_array)
{
    @adduct_row_array = sort cmp_row_rt keys %{$fattyacid_hash{$fattyacid_base}};
    $num_adducts = @adduct_row_array;
    
    ## for debugging, skip one-row lipids, since those are their own main
    #if ($num_adducts <= 1)
    #{
    #    next;
    #}

    for ($i = 0; $i < $num_adducts; $i++)
    {
        $row      = $adduct_row_array[$i];
        $adduct   = $fattyacid_hash{$fattyacid_base}{$row}{adduct};
        $rt_avg   = $row_rt_avg_array[$row];
        $rt_sd    = $row_rt_sd_array[$row];
        $data_avg = $row_data_avg_array[$row];
        $data_sum = $row_data_sum_array[$row];
        
        $rt_min = $rt_avg - $num_sd * $rt_sd;
        $rt_max = $rt_avg + $num_sd * $rt_sd;
        
        $group = $fattyacid_hash{$fattyacid_base}{$row}{group};

        $main  = $fattyacid_hash{$fattyacid_base}{$row}{main};
        if (!defined($main))
        {
            $main = '';
        }

        $main_ls = $fattyacid_hash{$fattyacid_base}{$row}{main_ls};
        if (!defined($main_ls))
        {
            $main_ls = '';
        }

        $lipid_class = $fattyacid_base;
        $lipid_class =~ s/\(.*//;
        $supclass = $class_annotation_hash{$lipid_class}{super};
        $subclass = $class_annotation_hash{$lipid_class}{sub};
        if (!defined($supclass) && !defined($subclass))
        {
            $supclass = '';
            $subclass = '';
            
            $missing_class_annotation_hash{$lipid_class} = 1;
        }
        elsif (!defined($supclass))
        {
            $supclass = '';
        }
        elsif (!defined($subclass))
        {
            $subclass = '';
        }
        
        # print columns before first Rt column
        for ($col = 0; $col < $first_rt_col; $col++)
        {
            $value = $row_data_array[$row][$col];
            if (!defined($value))
            {
                $value = '';
            }
        
            if ($col)
            {
                printf "\t";
            }
            
            printf "%s", $value;
        }
        
        # insert our new columns
        printf "\t%s",     $fattyacid_base;
        printf "\t%s",     $supclass;
        printf "\t%s",     $subclass;
        printf "\t\[%s\]", $adduct;
        printf "\t%s",     $group;
        printf "\t%s",     $main;
        if ($valid_summary_flag)
        {
            printf "\t%s", $main_ls;
        }
        printf "\t%s",     $rt_avg;
        printf "\t%s",     $data_avg;
        printf "\t%s",     $data_sum;

        # print the rest of the columns
        for ($col = $first_rt_col; $col < @header_col_array; $col++)
        {
            $value = $row_data_array[$row][$col];
            if (!defined($value))
            {
                $value = '';
            }
        
            if ($col)
            {
                printf "\t";
            }
            
            printf "%s", $value;
        }
        printf "\n";
    }
}


foreach $lipid_class (sort keys %missing_class_annotation_hash)
{
    printf STDERR "WARNING -- missing lipid class annotation:  %s\n",
        $lipid_class;
}
