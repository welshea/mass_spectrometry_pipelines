#!/usr/bin/perl -w

# 2021-08-19  improve/update documentation of csv2tsv_not_excel() function
# 2021-08-16  strip leading/trailing spaces from enclosed fields,
#             instead of the previous intentional preservation behavior
# 2021-08-16  bugfix removal of enclosing double-quotes with trailing spaces
# 2021-06-21  special case "" in a field by itself
# 2020-08-05  modify split regex to handle a few more invalid edge cases
# 2020-08-04  modify split regex to handle a few more invalid edge cases
# 2020-07-31  more improvements to robustness and optimizations
# 2020-07-28  speed up by replacing more quick and easy cases first
# 2020-07-24  bugfix, was condensing multiple tabs ahead of quoted field


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
# The following requirements of RFC 4180 are intentionally ignored:
#
# 1) "Fields containing line breaks (CRLF), double quotes, and commas
#     should be enclosed in double-quotes."
# 2) "Each record is located on a separate line, delimited by a line
#     break (CRLF)."
# 3) "Each line should contain the same number of fields throughout the file."
# 4) "The last field in the record must not be followed by a comma."
#
# I'm going to assume that all four of these are likely defined in the RFC
# due to the allowance of embedded newlines in 1).  The common implementation
# of 1) appears to assumes MS-DOS CRLF end-of-line (EOL) characters, rather
# than unix LF or a mix-and-match of some lines ending in CRLF and others in
# LF (frequently occurs aftering editing a CRLF file in Unix).  For example,
# Excel appears to assume that LF will be used to denote an embedded newline,
# with CRLF indicating a "real" newline.  This, as well as 2), is problematic
# in Unix software, which generally assumes LF is the EOL character, since
# lines containing embedded newlines are then split across multiple lines of
# text.  This is where 3) and 4) come in.  If we know that all lines contain
# the same number of fields, whenever we encounter a line that is too short,
# and its last field begins with a " that doesn't have a matching
# right-enclosing ", we can assume it is line-wrapped, read in the next line,
# and append it to the end of the original line, continuing to do so until we
# reach our fixed total number of fields per line.  Due to splitting of
# embedded newlines into multiple lines, escaped commas within embedded quotes
# could be left dangling at the end of a wrapped line instead of indicating a
# new field, so 4) is probably there to better handle these resulting
# problematic line-wrapped lines.  By requiring 3) and 4), we can sanity check
# the un-wrapped lines as well.  I wouldn't be surprised if 3) and 4) also
# make the embedded newline-aware parsing simpler?
#
# This tool is intended to be used within the Unix paradigm of line-based text
# processing tools, in which each line is generally independent of all other
# lines.  Allowing embedded newlines in 1) would break grep, sed, sort, wc,
# etc., and probably most Unix-based tools in general.  Wrapped lines could be
# reordered or filtered so that they are no longer sequential with their
# original line-wrapped partner lines or ordered for correct unwrapping.
# Thus, we do not support embedded newlines in CSV files.  If you wish to
# handle them, you will need to write a separate function to undo the line
# wrap and escape the embedded newlines as something such as \r or \n or \r\n
# prior to passing the results to this function.
#
# CRLF in 2) is not enforced, since we expect to receive data from Unix text
# files, which end in LF instead of CRLF.  Since 2), 3), and 4) are no longer
# needed to support 1), we have no reason to enforce them.  It is not uncommon
# for real data to have lines with variable numbers of fields.  After
# disallowing embedded newlines and ignoring 1, 2), and 3), we see no reason
# why terminal empty fields should be encoded as ,"" to satisfy 4) instead of
# simply using a terminal comma.  Excel also violates 4).  Thus, 4) is ignored
# as well, since it is frequently violated in real data.
#
# Additionally, the following requirement of RFC 4180 is loosened so as to
# attempt to accept invalid syntax wherever valid syntax could not be parsed.
# Excel also selectively loosens this requirement when reading in files:
#   
#   "If fields are not enclosed with double quotes, then double quotes may not
#    appear inside the fields."
#
#
# The specific invalid edge cases I tested for involve quotes inside fields
# that *aren't* enclosed in double-quotes, as well as leading/trailing spaces
# before or after enclosing double-quotes.
#
#
#   RFC 4180 on double-quotes:
#
#     "If double-quotes are used to enclose fields, then a double-quote
#      appearing inside a field must be escaped by preceding it with
#      another double quote.
#
#     "If fields are not enclosed with double quotes, then double quotes
#      may not appear inside the fields."
#
#
# This would appear to rule out treating fields surrounded in double-quotes
# as "enclosed" if there are any leading/trailing spaces, or if there are
# an even number of leading or trailing double-quotes.  However, this all
# depends on the defintion of "enclosed".
#
#
#   From Wikipedia's discussion on leading/trailing spaces:
#
#     "According to RFC 4180, spaces outside quotes in a field are not
#      allowed; however, the RFC also says that "Spaces are considered
#      part of a field and should not be ignored." and "Implementors
#      should 'be conservative in what you do, be liberal in what you
#      accept from others' (RFC 793, section 2.10) when processing CSV
#      files.""
#
#
# I choose to ignore leading/trailing spaces for the purposes of determing
# whether a field is enclosed within double-quotes or not.  I shall peserve
# them if the field is not enclosed ('spaces are considered part of a field
# and should not be ignored' and 'be conservative in what you do'), but remove
# them if the field is considered enclosed ('be liberal in what you accept
# from others').  Detecting and removing spaces that surround enclosing
# double-quotes allows for using spaces to visually align columns of variable
# width fields when displayed on a simple terminal or text editor that would,
# otherwise, not be vertically aligned.  I believe that this strikes a useful
# balance between simply preserving or removing all leading/trailing spaces.
# This attempts to follow the principle of "least surprise" to the user.
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
#
# A relatively simple example: ""word""
#
# If we consider it to be enclosed within double-quotes, then the remaining
# single double-quotes are invalid, but if it is not considered enclosed, then
# the original "" are also invalid.  How best should we deal with it?
#
# Various tools result in different results:
#
#   word""    python and Excel, this looks not-right to me
#   "wor      many conversion websites, truncation is just plain wrong
#   "word"    this function, which I think handles it the most reasonably
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


$filename = shift;

open FILE, "$filename" or die "ABORT -- can't open file $filename\n";

while(defined($line=<FILE>))
{
    $line_new = csv2tsv_not_excel($line);
    
    print "$line_new";
}
close FILE;
