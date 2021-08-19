#!/usr/bin/perl -w

# 2021-08-19:  update csv2tsv_not_excel() function
# 2021-03-05:  expand usage help message
# 2020-08-05:  more robust csb2tab_not_excel() function
# 2020-07-31:  faster and more robust csb2tab_not_excel() function
# 2020-07-28:  optimize csv2tab_not_excel() function
# 2020-07-24:  fix csv2tab, was condensing multiple tabs in front of quotes


use Scalar::Util qw(looks_like_number);
use POSIX;

# Generally only needs to run on POS ion data.
# The small number of NEG rows detected as OCT are either weakly correlated
#  with the other OCT rows, or not at all.  They also aren't convincing in
#  the m/z vs. retention time plot.

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


# use global array/hash
sub cmp_row_mz
{
    my $row_id_1;
    my $row_id_2;
    my $mz_1;
    my $mz_2;
    my $rt_1;
    my $rt_2;
    my $geomean_1;
    my $geomean_2;

    $mz_1 = $row_data[$a]{'mz'};
    $mz_2 = $row_data[$b]{'mz'};

    $rt_1 = $row_data[$a]{'rt'};
    $rt_2 = $row_data[$b]{'rt'};

    if ($mz_1 < $mz_2) { return -1; }
    if ($mz_1 > $mz_2) { return 1; }

    if ($rt_1 < $rt_2) { return -1; }
    if ($rt_1 > $rt_2) { return 1; }
    
    $geomean_1 = $mz_1 * $rt_1;
    $geomean_2 = $mz_2 * $rt_2;
    if ($geomean_1 < $geomean_2) { return -1; }
    if ($geomean_1 > $geomean_2) { return 1; }
    
    
    $row_id_1 = $row_data[$a]{'row_id'};
    $row_id_2 = $row_data[$b]{'row_id'};
    return $row_id_1 cmp $row_id_2;
}


# use global array/hash
sub cmp_row_rt
{
    my $row_id_1;
    my $row_id_2;
    my $mz_1;
    my $mz_2;
    my $rt_1;
    my $rt_2;
    my $geomean_1;
    my $geomean_2;

    $mz_1 = $row_data[$a]{'mz'};
    $mz_2 = $row_data[$b]{'mz'};

    $rt_1 = $row_data[$a]{'rt'};
    $rt_2 = $row_data[$b]{'rt'};

    if ($rt_1 < $rt_2) { return -1; }
    if ($rt_1 > $rt_2) { return 1; }

    if ($mz_1 < $mz_2) { return -1; }
    if ($mz_1 > $mz_2) { return 1; }
    
    $geomean_1 = $mz_1 * $rt_1;
    $geomean_2 = $mz_2 * $rt_2;
    if ($geomean_1 < $geomean_2) { return -1; }
    if ($geomean_1 > $geomean_2) { return 1; }
    
    
    $row_id_1 = $row_data[$a]{'row_id'};
    $row_id_2 = $row_data[$b]{'row_id'};
    return $row_id_1 cmp $row_id_2;
}


# This function was written to handle invalid CSV escaping in a consistent
# manner that will, hopefully, result in output fields that are closer to the
# intent of the malformed fields (guessing intent of malformed fields is,
# admittedly, fraught with peril).
#
# Due to the way that Excel handles quotes and using them to embed actual
# newlines (one line wrapped onto multiple lines), as well as how this script
# handles invalid CSV fields, the resulting tab-delimited text file may not
# load into Excel the same as if the .csv file had been loaded directly into
# Excel.  Thus the _not_excel part of the function name.  The resulting tsv
# file could, depending on how its quotes are arranged, result in strange
# Excel behavior.  The output should be fine for most other programs that read
# tab-delimited text, which is mainly what I wrote this for.  If you need
# something with better Excel compatability, I suggest using Python's
# csv.reader and csv.writer set to "dialect=csv.excel_tab".
#
#
# This implementation results in the following behavior:
#
#   1) One or more embedded tab characters (not escaped as \t) in a row are
#      either removed (at beginning/end of field or space on either side),
#      or replaced with a single space (if removing would otherwise merge
#      non-space text together).  Existing spaces are preserved, only
#      multiple unescaped tabs in a row are stripped/condensed.
#
#      While not exactly standards compliant (is there even a standard for
#      for how to handle embedded tabs??), I have observed this to be closest
#      to original intent in the wild, where embedded tabs are almost always
#      due to user input error (such as pressing the tab key in an attempt to
#      advance to the next field/line in a form, etc.) and improper validation
#      or sanitizing of the data prior to storing it in a database.
#
#      Internal \r and \n that are not part of the EOL are treated the same
#      way as we treat embedded tabs, since they usually arise from similar
#      issues.  This behavior can be disabled by commenting out a single
#      substitution line in the function below.
#
#      Feel free to modify the handling of embedded tabs at the bottom of the
#      function below if you desire different tab handling behavior.
#
#   2) Leading/trailing spaces are ignored for the purposes of determining
#      if a field is enclosed in double-quotes.
#
#   3) "" is treated as escaped double-quotes, regardless of the position
#      within the field.  This includes the beginning/end of the field (with
#      or without leading/trailing spaces).  Therefore, in order to be
#      considered enclosed, the beginning/end of the field (after ignoring
#      leading/trailing spaces) must have an odd number of double-quotes,
#      so that the outer double-quote on each end is not treated as escaped.
#
#        !! ""aa,bb"" is not considered to be enclosed in double-quotes !!
#
#   4) Fields consisting solely of "", with or without leading/trailing
#      spaces, are special-cased as empty text.  I have now observed this
#      empty field double-quoting behavior in the wild, where it was clearly
#      meant to indicate an empty field.
#
#   5) Enclosing double-quotes are stripped from the output fields, along
#      with any leading/trailing spaces surrounding them.  Unenclosed
#      leading/trailing spaces are left as-is.
#
#   6) After removing enclosing double-quotes, "" are unescaped back into ",
#      regardless of where they are located within the field.
#
#   7) *MOST*, but not all, white space is preserved.  See 1) and 5).
#
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
    my $dq  = "\x1D";    # (group separator) seed single char for split regex
    
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
        # /e to evaluates code to handle different capture cases correctly
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


sub debug_print_stuff
{
  printf "%s",   'Combo';
  printf "\t%s", 'row_id_1';
  printf "\t%s", 'row_id_2';
  printf "\t%s", 'count';
  printf "\t%s", 'rt_1';
  printf "\t%s", 'mz_1';
  printf "\t%s", 'rt_2';
  printf "\t%s", 'mz_2';
  printf "\t%s", 'rt_min';
  printf "\t%s", 'mz_min';
  printf "\t%s", 'error';
  printf "\t%s", 'delta rt';
  printf "\t%s", 'delta mz';
  printf "\t%s", 'delta rt norm';
  printf "\t%s", 'delta mz norm';
  printf "\t%s", 'rt_mz';
  printf "\t%s", 'rt_mz_norm';
  printf "\n";

  foreach $index_1 (@kept_array)
  {
    $row_id_1 = $row_data[$index_1]{'row_id'};
    $mz_1     = $row_data[$index_1]{'mz'};
    $rt_1     = $row_data[$index_1]{'rt'};

    @temp_array = sort cmp_row_mz keys %{$overlap_hash{$index_1}};
#    @temp_array = sort cmp_row_mz keys %{$filt_1_hash{$index_1}};
    $count = @temp_array + 1;		# add self back in
    
       printf "%s_%s", $row_id_1, $row_id_1;
       printf "\t%s", $row_id_1;
       printf "\t%s", $row_id_1;
       printf "\t%s", $count;
       printf "\t%s", $rt_1;
       printf "\t%s", $mz_1;
       printf "\t%s", $rt_1;
       printf "\t%s", $mz_1;
       printf "\t%s", $rt_1;
       printf "\t%s", $mz_1;
       printf "\n";
   
    foreach $index_2 (@temp_array)
    {
       $row_id_2 = $row_data[$index_2]{'row_id'};
       $mz_2     = $row_data[$index_2]{'mz'};
       $rt_2     = $row_data[$index_2]{'rt'};

       if (!defined($kept_hash{$index_2}))
       {
           $kept_hash{$index_2} = 0;
       }
       $kept_hash{$index_2} += 1;

       if ($row_id_2 eq $row_id_1)
       {
           next;
       }
       
       $delta_rt = $rt_2 - $rt_1;
       $delta_mz = $mz_2 - $mz_1;
       
       $ratio = $delta_rt;
       if ($delta_mz)
       {
           $ratio = $delta_rt / $delta_mz;
       }

       $rt_min = $rt_1;
       if ($rt_2 < $rt_min) { $rt_min = $rt_2; }
       $mz_min = $mz_1;
       if ($mz_2 < $mz_min) { $mz_min = $mz_2; }

       $division = abs($delta_mz / $ref_mz);
       $rounded  = floor($division + 0.5);
       $error    = $rounded * $ref_mz - abs($delta_mz);
       
       printf "%s_%s", $row_id_1, $row_id_2;
       printf "\t%s", $row_id_1;
       printf "\t%s", $row_id_2;
       printf "\t%s", $count;
       printf "\t%s", $rt_1;
       printf "\t%s", $mz_1;
       printf "\t%s", $rt_2;
       printf "\t%s", $mz_2;
       printf "\t%s", $rt_min;
       printf "\t%s", $mz_min;
       printf "\t%s", $error;
       printf "\t%s", $delta_rt;
       printf "\t%s", $delta_mz;
       printf "\t%s", ($delta_rt / $rt_min);
       printf "\t%s", ($delta_mz / ($mz_min));
       printf "\t%s", $ratio;
       printf "\t%s", $ratio * $mz_min / $rt_min;
       printf "\n";
    }
  }
}


$ref_mz              = 22.0131;
#$ref_mz *= 2;
#$ref_err             = 0.0015;	# super lax
#$ref_err             = 0.0012;	# loses pos_02665
$ref_err             = 0.0009;	# 0.0008 loses pos_05796

$rt_mz_lb = -0.003;		# -0.001255 is lowest observed in ~12 rt region
$rt_mz_ub =  0.01;
$rt_mz_norm_lb = -0.1;		# ~ -0.046 for trends ~12 rt
$rt_mz_norm_ub =  0.4;		# 0.3 loses a few minor ones, 0.31 is ok


#$mz_min_cutoff       = 4 * 22.0131 - $ref_err;	# keeps them all
#$mz_1_min_cutoff      = 11 * 22.0131 - $ref_err;
#$mz_1_min_cutoff      = 140;
#$mz_2_min_cutoff      = 240;	# loses a few somewhat questionable ones
#$mz_2_min_cutoff      = 6 * 44.0262;

#$mz_1_min_cutoff      = 7 * 22.0131 - $ref_err;
#$mz_2_min_cutoff      = 7 * 22.0131 - $ref_err;

# recommend minimum of 6x, any smaller and not-so-great looking groups appear
$mz_1_min_cutoff      = 11 * 22.0131 - $ref_err; # ~240
$mz_2_min_cutoff      = 11 * 22.0131 - $ref_err; # ~240

#$rt_1_cutoff          = 5;
#$rt_2_cutoff          = 5;
#$rt_1_cutoff          = 4.5;
#$rt_2_cutoff          = 4.5;

#$win_inc   = 0.1;
#$win_width = 1.0;

# override cutoffs
#$zero_rt_mz_norm_tol  = 0.015;    # 0.01 has some false positive cell lines
$zero_rt_mz_norm_tol  = 0.01;     # 0.011 loses some <240 m/z, 0.005 too lax
$rt_1_cutoff          = 0;
$rt_2_cutoff          = 0;
$rt_mz_lb             = -9E99;
$rt_mz_ub             =  9E99;
#$delta_mz_norm_cutoff = 2.0;
$delta_mz_norm_cutoff  = 9E99;
$delta_rt_norm_cutoff = 0.5;
#$rt_mz_norm_lb        = -0.25;
#$rt_mz_norm_ub        =  0.50;
$rt_mz_norm_lb        = -0.1;
$rt_mz_norm_ub        =  0.4;    # 0.35 appears safe?, 0.03 is NOT sfae
$nearest_delta_mz_cutoff = 3 * 22.0131 + $ref_err;
$mz_1_min_cutoff      = 11 * 22.0131 - $ref_err;    # < 11 not safe
$mz_2_min_cutoff      = 11 * 22.0131 - $ref_err;    # < 11 not safe

$filename = shift;
$outfile_oct_name = shift;
$outfile_metrics_name = shift;

if (!defined($filename) || $filename =~ /^-/)
{
    print "Usage:\n";
    print "\n";
    print "  scan_oct.pl mzmine_pos.csv [output_oct_rows.txt [output_oct_metrics.txt]]\n";
    print "\n";
    print "  Both CSV and tab-delimited input is accepted.\n";
    print "  (+) ion mode files yield best results, since OCT doesn't show up in (-).\n";
    print "  Any (-) ion mode rows detected as OCT are likely false positives.\n";
    print "\n";
    print "  If output file names are not specified, default file names will be used:\n";
    print "    oct_detected_rows.txt oct_sample_metrics.txt\n";
    print "\n";
    print "  Output files contain rows identified as potential OCT, and OCT abundance-\n";
    print "  related metrics, respectively.\n";

    exit(1);
}


if (!defined($outfile_oct_name))
{
    $outfile_oct_name = 'oct_detected_rows.txt';
}
if (!defined($outfile_metrics_name))
{
    $outfile_metrics_name = 'oct_sample_metrics.txt';
}

open INFILE, "$filename" or die "can't open $filename\n";

# header line
$line = <INFILE>;
$line =~ s/[\r\n]+//g;

# guess whether tsv or csv
# assume that it won't have any tabs at all if it is a csv file
$csv_flag     = 0;
$count_tabs   = 0;
$count_commas = 0;
@matches      = ($line =~ /\t/g);
$count_tabs   = @matches;
@matches      = ($line =~ /\,/g);
$count_commas = @matches;
if ($count_commas && $count_tabs == 0) { $csv_flag = 1; }

if ($csv_flag) { $line = csv2tsv_not_excel($line); }
$line =~ s/\"//g;

@array = split /\t/, $line;
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
    $header_col_array[$i] = $array[$i];
}

$row_id_col = $header_col_hash{'row ID'};
$mz_col     = $header_col_hash{'row m/z'};
$rt_col     = $header_col_hash{'row retention time'};
$name_col   = $header_col_hash{'row identity (all IDs)'};

if (!defined($row_id_col))
{
    print STDERR "WARNING -- 'row ID' column not found\n";
}
if (!defined($mz_col))
{
    print STDERR "ABORT -- 'row m/z' column not found\n";
    exit(1);
}
if (!defined($rt_col))
{
    print STDERR "ABORT -- 'row retention time' column not found\n";
    exit(1);
}


# conform sample columns
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];

    # strip mzXML from sample names
    $field =~ s/\.mzX?ML( Peak \S+)$/$1/i;

    if ($field =~ / Peak (\S+)$/i)
    {

        $second_word = $1;
        $second_word = lc $second_word;
        $field =~ s/ Peak (\S+)$/ Peak $second_word/i;
    }

    $header_col_array[$col] = $field;
}


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
}


# flag columns to remove, as they clutter up the file
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];
    
    if ($field =~ /\.mzX?ML[^.]+$/i ||
        $field =~ / Peak \S+$/)
    {
        # always keep peak height
        if ($field =~ / Peak height$/i)
        {
            next;
        }
        
        # only keep peak area if no peak height
        if ($peak_height_flag == 0 &&
            $field =~ / Peak area$/i)
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


$first_abundance_col = 9E99;
for ($col = 0; $col < @header_col_array; $col++)
{
    $field = $header_col_array[$col];

    if ($field =~ / Peak height$/i ||
        $field =~ / Peak area$/i)
    {
        if (!defined($col_to_remove_hash{$col}) &&
            $col < $first_abundance_col)
        {
            $first_abundance_col = $col;
        }
    }
}


# no sample columns found, take a guess
# assume mzmine format
if ($first_abundance_col == 9E99)
{
    printf STDERR "WARNING -- no Peak height or Peak area columns found\n";
    printf STDERR "WARNING == guessing first sample column\n";

    if (defined($name_col))
    {
        $first_abundance_col = $name_col;
    }
    else
    {
        $first_abundance_col = 0;
    }

    if ($row_id_col > $first_abundance_col)
    {
        $first_abundance_col = $row_id_col;
    }
    if ($mz_col > $first_abundance_col)
    {
        $first_abundance_col = $mz_col;
    }
    if ($rt_col > $first_abundance_col)
    {
        $first_abundance_col = $rt_col;
    }
    $first_abundance_col += 1;

}

#printf STDERR "First sample column: %s\n",
#    $header_col_array[$first_abundance_col];


$row = 0;
$max_col = 0;
while(defined($line=<INFILE>))
{
    $line =~ s/[\r\n]+//g;
    if ($csv_flag) { $line = csv2tsv_not_excel($line); }
    $line =~ s/\"//g;

    @array = split /\t/, $line;

    # clean up fields
    for ($col = 0; $col < @array; $col++)
    {
        $array[$col] =~ s/^\s+//;
        $array[$col] =~ s/\s+$//;
        $array[$col] =~ s/\s+/ /g;
    }

    # no row id column found, use the row # as the identifier
    if (!defined($row_id_col))
    {
        $row_id = $row;
    }
    # use the provided identifer
    else
    {
        $row_id = $array[$row_id_col];
    }

    $mz     = $array[$mz_col];
    $rt     = $array[$rt_col];

    # mz too low, don't even read it in
    # NOTE -- we need to store it all for the histogram
#    if ($mz < $mz_2_min_cutoff)
#    {
#        next;
#    }
    
    if ($row_id =~ /[A-Za-z0-9]/ &&
        is_number($mz) && is_number($rt))
    {
        $row_data[$row]{'row_id'} = $row_id;
        $row_data[$row]{'mz'}     = $mz;
        $row_data[$row]{'rt'}     = $rt;
        
        $sorted_rows[$row] = $row;
        
        # store line into data matrix, to extract samples from later
        for ($i = 0; $i < @array; $i++)
        {
            if ($array[$i] =~ /\S/ && is_number($array[$i]))
            {
                $data_matrix[$row][$i] = $array[$i];
            }
            
            if ($i > $max_col)
            {
                $max_col = $i;
            }
        }
        
        $row++;
    }
}
$num_rows = $row;



# identify dense regions at beginning and end to exclude
@sorted_rows = sort cmp_row_rt @sorted_rows;

$min_rt    = $row_data[$sorted_rows[0]]{'rt'};
$max_rt    = $row_data[$sorted_rows[$num_rows - 1]]{'rt'};
$range_rt  = $max_rt - $min_rt;

$win_width = $range_rt / 20.0;
$win_inc   = $win_width / 10.0;

$i_start    = 0;
$i_end      = 0;
$rt_start   = $min_rt;
#$div        = floor($rt_start / $win_inc);
#$rt_start   = $div * $win_inc;
$rt_end     = $rt_start + $win_width;
if ($rt_end > $max_rt)
{
    $rt_end = $max_rt;
}
$rt_center  = 0.5 * ($rt_start + $rt_end);

$num_histo_win = 0;
$histo_sum     = 0;
while ($rt_center < $max_rt)
{
    # advance $i_start to within $rt_start
    while ($i_start < $num_rows &&
           $row_data[$sorted_rows[$i_start]]{'rt'} < $rt_start - 1E-5)
    {
        $i_start++;
    }
    # advance $i_end to past $rt_end
    while ($i_end < $num_rows &&
           $row_data[$sorted_rows[$i_end]]{'rt'} < $rt_end + 1E-5)
    {
        $i_end++;
    }
    
    # we've gone off the end of the array, so we're finished scanning
    if ($i_start >= $num_rows)
    {
        last;
    }

    # count number of points within the rt range
    $count = 0;
    for ($j = $i_start; $j < $i_end; $j++)
    {
        $mz    = $row_data[$sorted_rows[$j]]{'mz'};
        $rt    = $row_data[$sorted_rows[$j]]{'rt'};

        if ($rt > $rt_start - 1E-5 &&
            $rt < $rt_end   + 1E-5)
        {
            $count++;
        }
    }

    $histo_data[$num_histo_win]{'count'} = $count;
    $histo_data[$num_histo_win]{'rt'}    = $rt_center;
    $histo_sum += $count;

    $temp_end = $i_end;
    if ($i_end >= $num_rows) { $temp_end = $num_rows - 1; }

#    printf STDERR "HISTO:\t%d\t%f\t%f\t%d\t%d\t%f\t%f\t%d\n",
#        $num_histo_win, $rt_start, $rt_end,
#        $i_start, $i_end,
#        $row_data[$sorted_rows[$i_start]]{'rt'},
#        $row_data[$sorted_rows[$temp_end]]{'rt'},
#        $count;

#    printf STDERR "HISTO:\t%d\t%f\t%d\n",
#        $num_histo_win, $rt_center, $count;

    # we've reached our final window
    if ($i_end >= $num_rows)
    {
        last;
    }

    $rt_start  += $win_inc;
    $rt_end    += $win_inc;
    $rt_center  = 0.5 * ($rt_start + $rt_end);

    $num_histo_win++;
}


if (0)
{
    # find median value of histogram counts
    @temp_array = ();
    for ($i = 0; $i < $num_histo_win; $i++)
    {
        $temp_array[$i] = $histo_data[$i]{'count'};
    }
    @temp_array = sort {$a<=>$b} @temp_array;
    $n_half = floor(0.5 * $num_histo_win);
    if ($num_histo_win % 2)
    {
        $histo_median = $temp_array[$n_half];
    }
    else
    {
        $histo_median = 0.5 * ($temp_array[$n_half - 1] +
                               $temp_array[$n_half]);
    }

#    $histo_noise = 0.5 * $histo_median;
#    printf STDERR "NOISE:\t%f\n", $histo_noise;
}

# find initial peak max
$histo_peak_max_bin = 0;
for ($j = 0; $j < $num_histo_win; $j++)
{
    $count = $histo_data[$j]{'count'};
    if ($count < $histo_data[$histo_peak_max_bin]{'count'})
    {
        last;
    }

    $histo_peak_max_bin = $j;
}
# walk down slope
$histo_peak_end = $histo_peak_max_bin;
for ($j = $histo_peak_max_bin; $j < $num_histo_win; $j++)
{
    $count = $histo_data[$j]{'count'};
    if ($count > $histo_data[$histo_peak_end]{'count'})
    {
        last;
    }
    
    $histo_peak_end = $j;
}
$rt_lb = $histo_data[$histo_peak_end]{'rt'};


#printf STDERR "RT_LB:\t%f\n", $rt_lb;


@sorted_rows = sort cmp_row_mz @sorted_rows;

%filt_1_hash = ();

for ($i = 0; $i < $num_rows; $i++)
{
    $index_1  = $sorted_rows[$i];

    $row_id_1 = $row_data[$index_1]{'row_id'};
    $mz_1     = $row_data[$index_1]{'mz'};
    $rt_1     = $row_data[$index_1]{'rt'};

    if ($rt_1 < $rt_1_cutoff)
    {
        next;
    }

    if ($mz_1 < $mz_1_min_cutoff)
    {
        next;
    }

    if ($rt_1 < $rt_lb)
    {
        next;
    }

    # *** TODO ***
    # store all deltas
    # sort on delta_rt
    # when mz is tied within error bounds, keep smallest delta_rt ??
    # *** doesn't work with pos_03692 relative to initial point, jumps back
    #
    # will need to check for distances relative to neighbors

    for ($j = 0; $j < $num_rows; $j++)
    {
        if ($i == $j)
        {
            next;
        }
    
        $index_2  = $sorted_rows[$j];

#        $row_id_2 = $row_data[$index_2]{'row_id'};
        $mz_2     = $row_data[$index_2]{'mz'};
        $rt_2     = $row_data[$index_2]{'rt'};

        if ($mz_2 < $mz_2_min_cutoff)
        {
            next;
        }

        if ($rt_2 < $rt_2_cutoff)
        {
            next;
        }

        if ($rt_2 < $rt_lb)
        {
            next;
        }

        $delta_mz    = $mz_2 - $mz_1;
        
        $division    = abs($delta_mz / $ref_mz);
        $rounded     = floor($division + 0.5);
        $error       = $rounded * $ref_mz - abs($delta_mz);
        
        $delta_rt    = $rt_2 - $rt_1;

        $rt_min = $rt_1;
        if ($rt_2 < $rt_min) { $rt_min = $rt_2; }
        $mz_min = $mz_1;
        if ($mz_2 < $mz_min) { $mz_min = $mz_2; }
        
        $rt_mz = $delta_rt;
        if ($delta_mz)
        {
            $rt_mz = $delta_rt / $delta_mz;
        }

        $delta_rt_norm = $delta_rt / $rt_min;
        $delta_mz_norm = $delta_mz / $mz_min;
        $rt_mz_norm = $rt_mz * $mz_min / $rt_min;

        # skip m/z too near the starting point
        if ($index_1 != $index_2 &&
            abs($delta_mz) <= $ref_err)
        {
            next;
        }
        
        # skip m/z that don't differ by a multiple of our target m/z
        if (abs($error) > $ref_err)
        {
            next;
        }

        $pairwise_hash{$index_1}{$index_2}{'rt_mz_norm'}    = $rt_mz_norm;
#       $pairwise_hash{$index_1}{$index_2}{'err'}           = $error;
#       $pairwise_hash{$index_1}{$index_2}{'rt_mz'}         = $rt_mz;
#       $pairwise_hash{$index_1}{$index_2}{'delta_rt_norm'} = $delta_rt_norm;
        
        if ($rt_mz >= $rt_mz_lb &&
            $rt_mz <= $rt_mz_ub &&
            $rt_mz_norm >= $rt_mz_norm_lb &&
            $rt_mz_norm <= $rt_mz_norm_ub &&
            abs($rt_mz_norm) > $zero_rt_mz_norm_tol &&
            abs($delta_rt_norm) <= $delta_rt_norm_cutoff &&
            abs($delta_mz_norm) <= $delta_mz_norm_cutoff)
        {
            $filt_1_hash{$index_1}{$index_2}{'delta_mz'} = $delta_mz;
        }
    }
}


%overlap_hash = ();
@filt_1_index_array = sort cmp_row_mz keys %filt_1_hash;

foreach $index_1 (@filt_1_index_array)
{
    @index_2_array = sort cmp_row_mz keys %{$filt_1_hash{$index_1}};

    # skip small groups
    if (@index_2_array < 3)
    {
        next;
    }
    
    foreach $index_2 (@index_2_array)
    {
        %temp_overlap_hash = ();

        $rt_mz_norm_1_2 = $pairwise_hash{$index_1}{$index_2}{'rt_mz_norm'};

        @index_3_array = sort cmp_row_mz keys %{$filt_1_hash{$index_2}};

        # skip small groups
        if (@index_3_array < 3)
        {
            next;
        }
        
        foreach $index_3 (@index_3_array)
        {
            if ($index_3 eq $index_1)
            {
                next;
            }

            # skip m/z too near each other
            $mz_1 = $row_data[$index_1]{'mz'};
            $mz_3 = $row_data[$index_3]{'mz'};
            $delta_mz = $mz_3 - $mz_1;
            if (abs($delta_mz) <= $ref_err)
            {
                next;
            }

            $division    = abs($delta_mz / $ref_mz);
            $rounded     = floor($division + 0.5);
            $error       = $rounded * $ref_mz - abs($delta_mz);

            # make sure the delta_mz is within tolerance
            if (abs($error) > $ref_err)
            {
                next;
            }
            
            # make sure ratios are in same direction
            $rt_mz_norm_1_3 = $pairwise_hash{$index_1}{$index_3}{'rt_mz_norm'};
            $rt_mz_norm_2_3 = $pairwise_hash{$index_2}{$index_3}{'rt_mz_norm'};

            if ($rt_mz_norm_1_2 * $rt_mz_norm_1_3 < 0 ||
                $rt_mz_norm_1_2 * $rt_mz_norm_2_3 < 0)
            {
#                # let slopes near zero through
#                if (abs($rt_mz_norm_1_2) > abs($rt_mz_norm_lb) ||
#                    abs($rt_mz_norm_1_3) > abs($rt_mz_norm_lb) ||
#                    abs($rt_mz_norm_2_3) > abs($rt_mz_norm_lb))
#                {
                    next;
#                }
            }

            if (defined($filt_1_hash{$index_1}{$index_3}))
            {
                $temp_overlap_hash{$index_3} = 1;
            }
        }

        @temp_overlap_array = keys %temp_overlap_hash;
        if (@temp_overlap_array >= 2)
        {
            $overlap_hash{$index_1}{$index_2} = 1;

            foreach $index_3 (@temp_overlap_array)
            {
                $overlap_hash{$index_1}{$index_3} = 1;
            }
        }
    }
}


@overlap_array = sort cmp_row_mz keys %overlap_hash;
#@overlap_array = sort cmp_row_mz keys %filt_1_hash;


# find closest overlaps in delta_mz
%best_delta_mz_hash = ();
foreach $index_1 (@overlap_array)
{
    @temp_array = sort cmp_row_mz keys %{$overlap_hash{$index_1}};
    
    foreach $index_2 (@temp_array)
    {
        $delta_mz = abs($filt_1_hash{$index_1}{$index_2}{'delta_mz'});
        
        if (!defined($best_delta_mz_hash{$index_1}))
        {
            $best_delta_mz_hash{$index_1} = $delta_mz;
        }
        if (!defined($best_delta_mz_hash{$index_2}))
        {
            $best_delta_mz_hash{$index_2} = $delta_mz;
        }
        
        if ($delta_mz < $best_delta_mz_hash{$index_1})
        {
            $best_delta_mz_hash{$index_1} = $delta_mz;
        }
        if ($delta_mz < $best_delta_mz_hash{$index_2})
        {
            $best_delta_mz_hash{$index_2} = $delta_mz;
        }
    }
}
# remove distant overlaps
foreach $index_1 (@overlap_array)
{
    if ($best_delta_mz_hash{$index_1} > $nearest_delta_mz_cutoff)
    {
        $row_id_1 = $row_data[$index_1]{'row_id'};
        $mz_1     = $row_data[$index_1]{'mz'};
        $rt_1     = $row_data[$index_1]{'rt'};

#        printf STDERR "REMOVED\t%s\t%f\t%f\t%f\n",
#            $row_id_1, $rt_1, $mz_1,
#            $best_delta_mz_hash{$index_1};

        delete $overlap_hash{$index_1};

        next;
    }

    @temp_array = sort cmp_row_mz keys %{$overlap_hash{$index_1}};
    
    foreach $index_2 (@temp_array)
    {
        if ($best_delta_mz_hash{$index_2} > $nearest_delta_mz_cutoff)
        {
            delete $overlap_hash{$index_1}{$index_2};
        }
    }
}
@overlap_array = sort cmp_row_mz keys %overlap_hash;


# filter slopes
# positive slopes are on the left side of the m/z vs. rt plot
# negative slopes are on the right side of the m/z vs. rt plot
%pos_slopes_rt_array = ();
%neg_slopes_rt_array = ();
$count_pos_slopes_rt = 0;
$count_neg_slopes_rt = 0;
foreach $index_1 (@overlap_array)
{
    @temp_array = sort cmp_row_mz keys %{$overlap_hash{$index_1}};
    
    foreach $index_2 (@temp_array)
    {
        $rt_2 = $row_data[$index_2]{'rt'};
        $rt_mz_norm_1_2 = $pairwise_hash{$index_1}{$index_2}{'rt_mz_norm'};
        
        if ($rt_mz_norm_1_2 > 0)
        {
            $pos_slopes_rt_array[$count_pos_slopes_rt++] = $rt_2;
        }
        if ($rt_mz_norm_1_2 < 0)
        {
            $neg_slopes_rt_array[$count_neg_slopes_rt++] = $rt_2;
        }
    }
}

$avg_pos_slopes_rt = 0;
$avg_neg_slopes_rt = 0;
$sd_pos_slopes_rt  = 0;
$sd_neg_slopes_rt  = 0;
foreach $rt (@pos_slopes_rt_array)
{
    $avg_pos_slopes_rt += $rt;
}
if (@pos_slopes_rt_array)
{
    $avg_pos_slopes_rt /= @pos_slopes_rt_array;
}
foreach $rt (@pos_slopes_rt_array)
{
    $diff = $avg_pos_slopes_rt - $rt;
    $sd_pos_slopes_rt += $diff * $diff;
}
if (@pos_slopes_rt_array)
{
    $sd_pos_slopes_rt = sqrt($sd_pos_slopes_rt / @pos_slopes_rt_array);
}
foreach $rt (@neg_slopes_rt_array)
{
    $avg_neg_slopes_rt += $rt;
}
if (@neg_slopes_rt_array)
{
    $avg_neg_slopes_rt /= @neg_slopes_rt_array;
}
foreach $rt (@neg_slopes_rt_array)
{
    $diff = $avg_neg_slopes_rt - $rt;
    $sd_neg_slopes_rt += $diff * $diff;
}
if (@neg_slopes_rt_array)
{
    $sd_neg_slopes_rt = sqrt($sd_neg_slopes_rt / @neg_slopes_rt_array);
}

printf STDERR "Pos slopes rt:\t%f\t%f\t%f\n",
    $avg_pos_slopes_rt - 2 * $sd_pos_slopes_rt,
    $avg_pos_slopes_rt,
    $avg_pos_slopes_rt + 2 * $sd_pos_slopes_rt,
printf STDERR "Neg slopes rt:\t%f\t%f\t%f\n",
    $avg_neg_slopes_rt - 2 * $sd_neg_slopes_rt,
    $avg_neg_slopes_rt,
    $avg_neg_slopes_rt + 2 * $sd_neg_slopes_rt,
$pos_slopes_rt_ub = $avg_pos_slopes_rt + 2 * $sd_pos_slopes_rt;
$neg_slopes_rt_lb = $avg_neg_slopes_rt - 2 * $sd_neg_slopes_rt;
$rt_mirror_point = 0.5 * ($pos_slopes_rt_ub + $neg_slopes_rt_lb);
printf STDERR "Pos/Neg rt mirror point:\t%f\n", $rt_mirror_point;


if ($count_pos_slopes_rt && $count_neg_slopes_rt &&
    $pos_slopes_rt_ub <= $neg_slopes_rt_lb)
{
  # remove slopes on the wrong side of the cutoff
  foreach $index_1 (@overlap_array)
  {
    @temp_array = sort cmp_row_mz keys %{$overlap_hash{$index_1}};
    
    foreach $index_2 (@temp_array)
    {
        $rt_2 = $row_data[$index_2]{'rt'};
        $rt_mz_norm_1_2 = $pairwise_hash{$index_1}{$index_2}{'rt_mz_norm'};
        
        if ($rt_mz_norm_1_2 > 0 && $rt_2 >= $rt_mirror_point)
        {
            delete $overlap_hash{$index_1}{$index_2};
        }
        elsif ($rt_mz_norm_1_2 < 0 && $rt_2 <= $rt_mirror_point)
        {
            delete $overlap_hash{$index_1}{$index_2};
        }
    }
  }
  @overlap_array = sort cmp_row_mz keys %overlap_hash;
}


# remove small clusters, make new filtered list
%kept_hash = ();
foreach $index_1 (@overlap_array)
{
    $row_id_1 = $row_data[$index_1]{'row_id'};
    $mz_1     = $row_data[$index_1]{'mz'};
    $rt_1     = $row_data[$index_1]{'rt'};

    @temp_array = sort cmp_row_mz keys %{$overlap_hash{$index_1}};
#    @temp_array = sort cmp_row_mz keys %{$filt_1_hash{$index_1}};
    $count = @temp_array + 1;		# add self back in
    
    # skip small clusters
    if ($count < 4)
    {
        next;
    }

    # count number of times we identified points within a cluster
    foreach $index_2 (@temp_array)
    {
       if (!defined($kept_hash{$index_1}))
       {
           $kept_hash{$index_1} = 0;
       }
       $kept_hash{$index_1} += 1;

       if (!defined($kept_hash{$index_2}))
       {
           $kept_hash{$index_2} = 0;
       }
       $kept_hash{$index_2} += 1;
    }
}

@kept_array = sort cmp_row_rt keys %kept_hash;


#debug_print_stuff();
#die;


open OUTFILE, ">$outfile_oct_name" or die "ABORT -- can't open $outfile_oct_name for writing\n";

# the small cutoff isn't very convincing... but filters out the weakest stuff
# if anything slips through, I'll want to know so I can tweak the algorithm
$num_detected = 0;
if (@kept_array >= 16)
{
    $num_detected = @kept_array;
    printf STDERR "Finished: %d potential OCT detected\n", $num_detected;

    printf OUTFILE "%s\t%s\t%s\t%s\n",
        'Row_ID', 'Row_ID_merged',
        'Retention_Time', 'm/z';

    foreach $index_1 (@kept_array)
    {
        $row_id_1 = $row_data[$index_1]{'row_id'};
        $mz_1     = $row_data[$index_1]{'mz'};
        $rt_1     = $row_data[$index_1]{'rt'};

        $row_id_pos_neg = $row_id_1;
        if ($row_id_1 =~ /^[0-9]+$/)
        {
            $row_id_pos_neg = sprintf "%05d", $row_id_1;
        }
        $pos_flag = 0;
        $neg_flag = 0;
        if ($filename =~ /pos/i)
        {
            $pos_flag = 1;
        }
        if ($filename =~ /neg/i)
        {
            $neg_flag = 1;
        }
        if (($pos_flag || $neg_flag) &&
            !($pos_flag && $neg_flag))
        {
            if ($pos_flag)
            {
                $row_id_pos_neg = sprintf "%s_%s",
                    'pos', $row_id_pos_neg;
            }
            elsif ($neg_flag)
            {
                $row_id_pos_neg = sprintf "%s_%s",
                    'neg', $row_id_pos_neg;
            }
        }
        
        printf OUTFILE "%s\t%s\t%f\t%f\n",
            $row_id_1, $row_id_pos_neg,
            $rt_1, $mz_1;
    }
}
else
{
    printf STDERR "Finished: 0 potential OCT detected\n";
    exit(0);
}
close OUTFILE;


# scan over sample abundances, determine if it is probably logged or not
$count_points_all   = 0;
$count_points_large = 0;
$count_points_small = 0;
for ($row = 0; $row < $num_rows; $row++)
{
    for ($col = $first_abundance_col; $col < $max_col; $col++)
    {
        if (defined($col_to_remove_hash{$col}))
        {
            next;
        }

        $value = $data_matrix[$row][$col];
        
        if (!defined($value) || !($value =~ /\S/))
        {
            $value = '';
        }
        
        if ($value ne '')
        {
            $count_points_all++;
            
            if ($value >= 100)
            {
                $count_points_large++;
            }
            elsif ($value > 0)
            {
                $count_points_small++;
            }
        }
    }
}

#printf STDERR "%d\t%d\t%d\n",
#    $count_points_all, $count_points_large, $count_points_small;

$already_log_flag = 0;
if ($count_points_small > $count_points_large)
{
    $already_log_flag = 1;
}


# take the log2 if we haven't already, delete missing data
$one_over_log2 = 1.0 / log(2.0);
for ($row = 0; $row < $num_rows; $row++)
{
    for ($col = $first_abundance_col; $col < $max_col; $col++)
    {
        if (defined($col_to_remove_hash{$col}))
        {
            next;
        }

        $value = $data_matrix[$row][$col];
        
        if (!defined($value) || !($value =~ /\S/))
        {
            $value = '';
        }
        
        if ($value ne '')
        {
            if ($already_log_flag == 0)
            {
                if ($value > 0)
                {
                    $data_matrix[$row][$col] =
                        log($value) * $one_over_log2;
                }
                else
                {
                    delete $data_matrix[$row][$col];
                }
            }
        }
    }
}


# sum log abundances
@col_sum_all_array     = ();
@col_sum_oct_array     = ();
@col_sum_non_oct_array = ();
for ($col = $first_abundance_col; $col <= $max_col; $col++)
{
    if (defined($col_to_remove_hash{$col}))
    {
        next;
    }

    $col_sum_all_array[$col]{sum}       = 0;
    $col_sum_oct_array[$col]{sum}       = 0;
    $col_sum_non_oct_array[$col]{sum}   = 0;
    $col_sum_all_array[$col]{count}     = 0;
    $col_sum_oct_array[$col]{count}     = 0;
    $col_sum_non_oct_array[$col]{count} = 0;

    for ($row = 0; $row < $num_rows; $row++)
    {
        $value = $data_matrix[$row][$col];
        
        if (!defined($value) || !($value =~ /\S/))
        {
            $value = '';
        }

        # do not sum values <= 0
        # zeroes due to unlogged values of 1 are unlikely
        if ($value ne '' && $value > 0)
        {
            $col_sum_all_array[$col]{sum}   += $value;
            $col_sum_all_array[$col]{count} += 1;
            
            if (defined($kept_hash{$row}))
            {
                $col_sum_oct_array[$col]{sum}   += $value;
                $col_sum_oct_array[$col]{count} += 1;
            }
            else
            {
                $col_sum_non_oct_array[$col]{sum}   += $value;
                $col_sum_non_oct_array[$col]{count} += 1;
            }
        }
    }
}


open OUTFILE, ">$outfile_metrics_name" or die "ABORT -- can't open $outfile_metrics_name for writing\n";
print OUTFILE "Sample\t%OCT_vs_non-OCT\n";

# print OCT percentages
for ($col = $first_abundance_col; $col <= $max_col; $col++)
{
    if (defined($col_to_remove_hash{$col}))
    {
        next;
    }

    $count_all     = $col_sum_all_array[$col]{count};
    $count_oct     = $col_sum_oct_array[$col]{count};
    $count_non_oct = $col_sum_non_oct_array[$col]{count};
    $sum_all       = $col_sum_all_array[$col]{sum};
    $sum_oct       = $col_sum_oct_array[$col]{sum};
    $sum_non_oct   = $col_sum_non_oct_array[$col]{sum};
    $avg_all       = 0;
    $avg_oct       = 0;
    $avg_non_oct   = 0;
    
    if ($count_all)
    {
        $avg_all = $sum_all / $count_all;
    }
    if ($count_oct)
    {
        $avg_oct = $sum_oct / $count_oct;
    }
    if ($count_non_oct)
    {
        $avg_non_oct = $sum_non_oct / $count_non_oct;
    }
    
    $signal_fraction = 0.0;
    if ($sum_non_oct)
    {
        $signal_fraction = $sum_oct / $sum_non_oct;
#        $signal_fraction = $sum_oct / $sum_all;
    }
    elsif ($sum_oct)
    {
        $signal_fraction = 1.0;
    }
    $metric = 100 * $signal_fraction;

    # assume data is log2
#    $log2_ratio = $avg_oct - $avg_non_oct;
#    $metric = $log2_ratio;
    
    
    $sample_name = $header_col_array[$col];
    
#    printf STDERR "Sample:\t%s\t%f\n", $sample_name, $metric;
    printf OUTFILE "%s\t%f\n", $sample_name, $metric;
}
close OUTFILE;
