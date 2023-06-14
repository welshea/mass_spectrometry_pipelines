#!/usr/bin/perl -w

# Be sure the annotation file hasn't gone through Excel,
# since Excel will have done some problematic escapes of ", etc.


# 2023-06-14:  initial release


sub fix_some_parentheses
{
    my @parts_array;
    my $name_orig = $_[0];
    my $name_new;
    my $class;
    my $not_class;
    my $count_open;
    my $count_open2;
    my $count_close;
    my $count_close2;
    my $unmatched_flag;
    my $ok_flag;
    my $i;

    # example entry with extra ) at end
    #    MGDG(18:4;O/16:4;O))
    #
    # lots of entries missing terminal )
    #    GlcCer(d16:1(4E)(15Me)/20:0(2OH)
    #
    if ($name_orig =~ /^([A-Za-z\(][-A-Za-z0-9 \(\)]*[A-Za-z0-9])(\([^\(\)]+:.*[_\/]*){1,3}$/)
    {
        $class     = $1;
        $not_class = $2;
        
        $count_open  = $class =~ tr/\(//;    # empty // means no replacement
        $count_close = $class =~ tr/\)//;

        # unmatched parantheses in class
        $unmatched_flag = 0;
        if ($count_open != $count_close)
        {
            $unmatched_flag = 1;
        
            # LMSP0601CC01: NAME = KDN)GD1a(d18:1/16:0)
            # LMSP0601CB01: NAME = Neu5Ac)GD1a(d18:1/16:0)
            #
            # based on synonyms and other LIMDs, it looks like these
            # should be "(KDN,Neu5Ac)GD1a"
            #
            if ($count_close && $count_open == 0 &&
                $class =~ /^(KDN|Neu5Ac)\)GD1a/)
            {
                $class =~ s/^(KDN|Neu5Ac)/\(KDN,Neu5Ac/;
                
                $message = 'missing either KDN or Neu5Ac in (KDN,Neu5Ac)GD1a';
            }
        }


        $count_open  = $not_class =~ tr/\(//;    # empty // means no replacement
        $count_close = $not_class =~ tr/\)//;

        # unmatched parantheses in the rest of the name
        if ($count_open != $count_close)
        {
            $unmatched_flag = 1;

            # more ( than )
            if ($count_open > $count_close)
            {
                if ($count_open - $count_close == 1)
                {
                    # missing terminal )
                    #    PC(26:2/26:2
                    #    LPC(18:1
                    #    Cer(20:0
                    #    Cer(22:0
                    #      ...
                    #    Cer(26:0
                    if (!($not_class =~ /\)$/))
                    {
                        $not_class .= ')';

                        $message = 'missing terminal ) entirely';
                    }
                
                    # insert missing ) before closing ]
                    #    DGDG(18:2(9Z,12Z)/18:2(9Z,12Z)[15(R)OH-18:2(9Z,12Z])
                    #
                    # I'm not sure this is the correct fix?
                    #
                    elsif ($not_class =~ /(\([^\(\)\[\]]+)\]/)
                    {
                        $not_class =~ s/(\([^\(\)\[\]]+)\]/$1\)\]/;

                        $message = 'missing ) before ] in (x,y]';
                    }
                    
                    # delete extra ( in 3((
                    #     GlcADG(18:3((9Z,12Z,15Z)/16:0)
                    elsif ($not_class =~ /\((\([^\(\)]+\)[_\/:])/)
                    {
                        $not_class =~ s/\((\([^\(\)]+\)[_\/:])/$1/;

                        $message = 'extra ( in x((...)';
                    }

                    # each part is OK, just need a closing ) added
                    elsif ($not_class =~ /^\(/)
                    {
                        @part_array    = split /[:_\/]/, $not_class;
                        $part_array[0] =~ s/^\(//;
                        
                        $ok_flag = 1;
                        for ($i = 0; $i < @part_array; $i++)
                        {
                            $count_open2  = $part_array[$i] =~ tr/\(//;
                            $count_close2 = $part_array[$i] =~ tr/\)//;
                            
                            if ($count_open2 != $count_close2)
                            {
                                $ok_flag = 0;
                                last;
                            }
                        }
                        
                        if ($ok_flag)
                        {
                            $not_class .= ')';

                            $message = 'missing additional terminal )';
                        }
                    }
                }
            }
            # more ) than (
            elsif ($count_close - $count_open == 1)
            {
                # insert ( after :1
                #    DGMG(16:1(3E)3Me,7ME,11Me,15Me)/0:0)
                #    this also has a capitalization typo in 7ME
                #
                if ($not_class =~ /:([^_:\(]+)([^_:]+,)+/)
                {
                    # fix the parentheses typo
                    $not_class =~ s/:([^_:\(]+)([^_:]+,)+/:$1\($2/;
                    
                    # fix ME capitalization typo while we're at it
                    $not_class =~ s/\b([0-9])+ME\b/$1Me/g;

                    $message = 'missing opening ( after :x';
                }

                # (t18:0) shouldn't have () around it
                #    Cer(t18:0)/24:0(2OH[R]))
                #    PI-Cer(t18:0)/24:0(2OH[S]))
                elsif ($not_class =~ /^(\([mdt][0-9]+:[0-9]+)\)/)
                {
                    $not_class =~ s/^(\([mdt][0-9]+:[0-9]+)\)/$1/;

                    $message = ') in (...) that should not be enclosed';
                }

                # strip extra ) at end
                #    MGDG(18:4;O/16:4;O))
                #    MG(21:1(5Z)20Me))
                elsif ($not_class =~
                       /^(\([^\(\)]+(\([^\(\)]+\))*[^\(\)]*\))\)$/)
                {
                    $not_class =~
                        s/^(\([^\(\)]+(\([^\(\)]+\))*[^\(\)]*\))\)$/$1/;

                    $message = 'extra terminal )';
                }
            }
        }
        
        # we had unmatched parentheses somewhere
        if ($unmatched_flag)
        {
            $name_new = $class . $not_class;
            
            # fixed it with a rule, return corrected string
            if ($name_new ne $name_orig)
            {
                return $name_new;
            }
            # did not match current rules, print a warning
            else
            {
                $message = 'UNCORRECTED';
            }
        }
    }
    
    return $name_orig;
}


# force lines to be flushed
# otherwise, STDOUT and STDERR can blend together on same output line...
select(STDERR); $| = 1;
select(STDOUT); $| = 1;


$lmaps_header_names_col_hash{'NAME'}            = -1;
$lmaps_header_names_col_hash{'SYSTEMATIC_NAME'} = -1;
$lmaps_header_names_col_hash{'ABBREVIATION'}    = -1;
$lmaps_header_names_col_hash{'SYNONYMS'}        = -1;


$lipidmaps_filename = shift;


open LIPIDMAPS, "$lipidmaps_filename" or die "can't open Lipid Maps file $lipidmaps_filename\n";


# read in lipidmaps header
$line = <LIPIDMAPS>;
$line =~ s/[\r\n]+//g;
@array = split /\t/, $line;
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    
    if ($array[$i] =~ /[A-Za-z0-9]/)
    {
        $lmaps_header_col_hash{$array[$i]} = $i;
        $lmaps_header_col_array[$i]        = $array[$i];
    }
}

$lmaps_id_col      = $lmaps_header_col_hash{'LM_ID'};

if (!defined($lmaps_id_col))
{
    printf STDERR "ABORT -- LM_ID column not found in file %s\n",
        $lipidmaps_filename;

    exit(2);
}


# update columns for header names hash
@lmaps_header_names_array     = sort keys %lmaps_header_names_col_hash;
@lmaps_header_names_col_array = ();
for ($i = 0, $j = 0; $i < @lmaps_header_names_array; $i++)
{
    $header = $lmaps_header_names_array[$i];
    $col    = $lmaps_header_col_hash{$header};

    if (defined($col))
    {
        $lmaps_header_names_col_hash{$header} = $col;
        $lmaps_header_names_col_array[$j++]   = $col;
    }
    else
    {
        delete $lmaps_header_names_col_hash{$header};
    }
}
@lmaps_header_names_array = sort keys %lmaps_header_names_col_hash;


if (@lmaps_header_names_array == 0)
{
    print "ABORT -- no name related columns found in file %s\n",
        $lipidmaps_filename;
    exit(2);
}



# print header line
print "$line\n";

printf STDERR "%s\t%s\t%s\t%s\t%s\n",
    'LM_ID', 'Field', 'TypeOfError', 'Original', 'Corrected';

# print the rest of the file, correcting parentheses as we go
while(defined($line=<LIPIDMAPS>))
{
    $line =~ s/[\r\n]+//g;
    
    @array = split /\t/, $line;
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
    }
    
    $lm_id   = $array[$lmaps_id_col];
    
    if (!defined($lm_id) || !($lm_id =~ /[A-Za-z0-9]/))
    {
        next;
    }


    # detect and correct for unmatched parentheses
    foreach $col (@lmaps_header_names_col_array)
    {
        $name_str = $array[$col];

        if (defined($name_str) && $name_str =~ /[A-Za-z0-9]/)
        {
            @split_array = split /\s*\|+\s*/, $name_str;
            
            for ($i = 0; $i < @split_array; $i++)
            {
                $name_orig = $split_array[$i];
                $message   = '';

                $name_new  = fix_some_parentheses($name_orig);
                
                if ($message ne '')
                {
                    printf STDERR "%s\t%s\t%s\t%s\t%s\n",
                        $lm_id, $lmaps_header_col_array[$col],
                        $message,
                        $name_orig, $name_new;
                
                    $split_array[$i] = $name_new;
                }
            }
        }
        
        $name_str_new = join " | ", $name_str;
        
        $array[$col] = $name_str_new;
    }
    
    $line_new = join "\t", @array;
    
    print STDOUT "$line_new\n";
}
close LIPIDMAPS;
