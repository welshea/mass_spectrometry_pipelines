#!/usr/bin/perl -w


# 2022-10-10  initial commit to git


# Be sure the annotation file hasn't gone through Excel,
# since Excel will have done some problematic escapes of ", etc.

# For now, require the data input file to have gone all the way through
# the BBSR pipeline, so that the ParentForumla and Lipid columns are added.
# I may add support for re-deriving them if they are missing, eventually,
# from the FattyAcid (lipid) and Class (lipid), and LipidIon (adduct)
# columns.


use File::Basename;


$lmaps_header_names_col_hash{'NAME'}            = -1;
$lmaps_header_names_col_hash{'SYSTEMATIC_NAME'} = -1;
$lmaps_header_names_col_hash{'ABBREVIATION'}    = -1;
$lmaps_header_names_col_hash{'SYNONYMS'}        = -1;


$lipidmaps_filename = shift;
$data_filename = shift;


if (!defined($lipidmaps_filename) || !defined($data_filename))
{
    $program_name = basename($0);
    
    print "Usage: $program_name lipidmaps_sdf_parsed.txt data.txt\n";
    exit(1);
}



open LIPIDMAPS, "$lipidmaps_filename" or die "can't open Lipid Maps file $lipidmaps_filename\n";
open DATA,      "$data_filename"      or die "can't open data file $lipidmaps_filename\n";


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
$lmaps_formula_col = $lmaps_header_col_hash{'FORMULA'};
$lmaps_abbrev_col  = $lmaps_header_col_hash{'ABBREVIATION'};

if (!defined($lmaps_id_col))
{
    printf STDERR "ABORT -- LM_ID column not found in file %s\n",
        $lmaps_filename;

    exit(2);
}
if (!defined($lmaps_formula_col))
{
    printf STDERR "ABORT -- FORMULA column not found in file %s\n",
        $lmaps_filename;

    exit(2);
}
if (!defined($lmaps_abbrev_col))
{
    printf STDERR "ABORT -- ABBREVIATION column not found in file %s\n",
        $lmaps_filename;

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


# read in annotation
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
    $formula = $array[$lmaps_formula_col];
    $abbrev  = $array[$lmaps_abbrev_col];
    
    if (!defined($lm_id) || !($lm_id =~ /[A-Za-z0-9]/))
    {
        next;
    }

    # store formula --> LM_ID lookup table    
    if (defined($formula) && $formula =~ /[A-Za-z0-9]/)
    {
        $formula_lmid_hash{$formula}{$lm_id} = 1;
    }


    # store name --> LM_ID lookup table    
    foreach $col (@lmaps_header_names_col_array)
    {
        $name_str = $array[$col];

        if (defined($name_str) && $name_str =~ /[A-Za-z0-9]/)
        {
            @split_array = split /\s*\|+\s*/, $name_str;
            foreach $name (@split_array)
            {
                $name_lc = lc $name;

                # strip inner (...) lists with numbers/letters
                # strip (O-, (P-
                # LipidSearch doesn't go into this much detail
                if (1)
                {
                    if ($name_lc =~ m/^([^(]+\([^(]+:)(.*)\)$/)
                    {
                        $first = $1;
                        $inner = $2;

                        # strip inner (...)
                        if ($inner =~ /\(/)
                        {
                            $inner =~ s/\(([0-9]+[a-z],*)+\)//g;
                            $name_lc = $first . $inner . ')';
                        }
                        
                        # strip (O-, (P-
                        $name_lc =~ s/\([OP]-/\(/ig;
                    }
                }
                
                if ($name_lc =~ /[a-z0-9]/)
                {
                    $name_lc_lmid_hash{$name_lc}{$lm_id} = 1;
                }
                
                # abbreviate coenzymes
                if ($name_lc =~ /^coenzyme\s+(.*)/)
                {
                    $string = $1;
                    $abbrev = 'co(' . $string . ')';
                    
                    #printf STDERR "FOOBAR\t%s\t%s\n", $name_lc, $abbrev;

                    $name_lc_lmid_hash{$abbrev}{$lm_id} = 1;
                }
            }
        }
    }
    
    # store abbreviation lookup hash
    # conform abbreviation to what we expect from LipidSearch
    if ($abbrev =~ /[A-Za-z0-9]/)
    {
        $abbrev =~ s/^(\S+\s+)[OP]-/$1/;

        $abbrev =~ s/ /\(/;
        $abbrev =~ s/;.*//;
        $abbrev .= ')';
        
        $abbrev_lmid_hash{$abbrev}{$lm_id} = 1;
    }

    # store all annotation fields by LM_ID
    for ($col = 0; $col < @lmaps_header_col_array; $col++)
    {
        $field = $array[$col];
        if (!defined($field) || !($field =~ /\S/))
        {
            next;
        }
        
        $lmid_col_annotation_hash{$lm_id}{$col} = $field;
    }
}
close LIPIDMAPS;

# array of all names, for subtring matches later
@lmaps_name_lc_all_array = sort keys %name_lc_lmid_hash;



# read in data header
$line = <DATA>;
$line =~ s/[\r\n]+//g;
@array = split /\t/, $line;
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    
    if ($array[$i] =~ /[A-Za-z0-9]/)
    {
        $data_header_col_hash{$array[$i]} = $i;
        $data_header_col_array[$i]        = $array[$i];
    }
}

$data_name_col    = $data_header_col_hash{'Lipid'};
$data_formula_col = $data_header_col_hash{'ParentFormula'};
$data_lgroup_col  = $data_header_col_hash{'LipidGroup'};

if (!defined($data_name_col))
{
    printf STDERR "ABORT -- Lipid column not found in file %s\n",
        $data_filename;

    exit(2);
}
if (!defined($data_formula_col))
{
    printf STDERR "ABORT -- ParentFormula column not found in file %s\n",
        $data_filename;

    exit(2);
}
if (!defined($data_lgroup_col))
{
    printf STDERR "ABORT -- LipidGroup column not found in file %s\n",
        $data_filename;

    exit(2);
}


# match LM_IDs to data
$row = 0;
while(defined($line=<DATA>))
{
    $line =~ s/[\r\n]+//g;
    
    @array = split /\t/, $line;
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
    }
    
    $line_new = join "\t", @array;
    $data_row_line_array[$row] = $line_new;
    
    $name    = $array[$data_name_col];
    $name_lc = lc $name;
    $formula = $array[$data_formula_col];
    $lgroup  = $array[$data_lgroup_col];

    # remove any ions that may be there
    $name_lc =~ s/[+-][A-Za-z0-9+-]+//g;

    # replace underscores with slashes for name mapping
    $name_lc_orig  = $name_lc;
    $name_lc_slash = $name_lc;
    $name_lc_slash =~ tr/_/\//;
    
    %lmid_hits_hash = ();
    $num_hits       = 0;
    $match_type     = '';

    # name lookup, original, generally with underscores
    if (defined($name_lc_lmid_hash{$name_lc_orig}) &&
        defined($formula_lmid_hash{$formula}))
    {
        @temp_array = keys %{$name_lc_lmid_hash{$name_lc_orig}};
        foreach $lm_id (@temp_array)
        {
            if (defined($formula_lmid_hash{$formula}{$lm_id}))
            {
                $match_type             = '01_strict';
                $lmid_hits_hash{$lm_id} = $match_type;
                $num_hits++;

                #printf STDERR "%s  %s  %s  %s  %s\n",
                #    $row, $name, $lm_id, $name_lc_orig, $match_type;
            }
        }
    }

    # name lookup, replace underscores with slashes
    if ($num_hits == 0)
    {
        if (defined($name_lc_lmid_hash{$name_lc_slash}) &&
            defined($formula_lmid_hash{$formula}))
        {
            @temp_array = keys %{$name_lc_lmid_hash{$name_lc_slash}};
            foreach $lm_id (@temp_array)
            {
                if (!defined($lmid_hits_hash{$lm_id}) &&
                    defined($formula_lmid_hash{$formula}{$lm_id}))
                {
                    $match_type             = '02_relaxed';
                    $lmid_hits_hash{$lm_id} = $match_type;
                    $num_hits++;

                    #printf STDERR "%s  %s  %s  %s  %s\n",
                    #    $row, $name, $lm_id, $name_lc_slash, $match_type;
                }
            }
        }
    }
    
    # LipidGroup / ABBREVIATION lookup
    # we'll use the lipid group conforming for formula sanity checking later
    if ($num_hits == 0)
    {
        # conform LipidGroup to what we expect from conformed abbrev
        $lgroup_conformed = $lgroup;

        # strip lowercase letters from in front/behind numbers
        if ($lgroup =~ m/^([^(]+\()(.*)\)/)
        {
            $first = $1;
            $inner = $2;

            $inner =~ s/[a-z]*([0-9]+)[a-z]*/$1/g;
            
            $lgroup_conformed = $first . $inner . ')';
        }
        
        # remove ions
        $lgroup_conformed =~ s/[+-][A-Za-z0-9+-]+//g;
        
        # conform abbreviations
        $lgroup_conformed =~ s/^AcCa\(/CAR\(/;    # Acyl carnitine
        $lgroup_conformed =~ s/^ChE\(/CE\(/;      # Cholesterol esters
        $lgroup_conformed =~ s/^Hex1(?![0-9])/Hex/;      # 1 is left off
        $lgroup_conformed =~ s/^SPH\(/SPB\(/;     # Sphingoid bases
        
        # check abbreviations
        if (defined($abbrev_lmid_hash{$lgroup_conformed}) &&
            defined($formula_lmid_hash{$formula}))
        {
            @temp_array = keys %{$abbrev_lmid_hash{$lgroup_conformed}};
            foreach $lm_id (@temp_array)
            {
                # must match formula as well
                if (!defined($lmid_hits_hash{$lm_id}) &&
                    defined($formula_lmid_hash{$formula}{$lm_id}))
                {
                    $match_type             = '03_lipidgroup';
                    $lmid_hits_hash{$lm_id} = $match_type;
                    $num_hits++;
                }
            }
        }

        # also check synonyms, names, etc.
        $lgroup_conformed_lc = lc $lgroup_conformed;
        if (defined($name_lc_lmid_hash{$lgroup_conformed_lc}) &&
            defined($formula_lmid_hash{$formula}))
        {
            @temp_array = keys %{$name_lc_lmid_hash{$lgroup_conformed_lc}};
            foreach $lm_id (@temp_array)
            {
                # must match formula as well
                if (!defined($lmid_hits_hash{$lm_id}) &&
                    defined($formula_lmid_hash{$formula}{$lm_id}))
                {
                    $match_type             = '03_lipidgroup';
                    $lmid_hits_hash{$lm_id} = $match_type;
                    $num_hits++;
                }
            }
        }
    }
    
    # formula lookup, must share class with abbreviation
    if (1 && $num_hits == 0)
    {
        if (defined($formula_lmid_hash{$formula}))
        {
            @temp_array = keys %{$formula_lmid_hash{$formula}};
            foreach $lm_id (@temp_array)
            {
                $abbrev = $lmid_col_annotation_hash{$lm_id}{$lmaps_abbrev_col};
                if (!defined($abbrev))
                {
                    $abbrev = '';
                }
                $abbrev =~ s/^(\S+\s+)[OP]-/$1/;
                $abbrev =~ s/ /\(/;
                $abbrev =~ s/;.*//;
                $abbrev .= ')';
                
                if ($abbrev =~ /(^[^)]+)\(/)
                {
                    $lmaps_class = $1;
                    
                    if ($lgroup_conformed =~ /(^[^)]+)\(/)
                    {
                        $lgroup_class = $1;

                        if ($lmaps_class eq $lgroup_class &&
                            !defined($lmid_hits_hash{$lm_id}))
                        {
                            $match_type             = '04_lipidclass';
                            $lmid_hits_hash{$lm_id} = $match_type;
                            $num_hits++;
                        }
                    }
                }
            }
        }
    }

    # formula lookup, no sanity checking
    # these look like poor matches to me, disable them
    if (0 && $num_hits == 0)
    {
        if (defined($formula_lmid_hash{$formula}))
        {
            @temp_array = keys %{$formula_lmid_hash{$formula}};
            foreach $lm_id (@temp_array)
            {
                if (!defined($lmid_hits_hash{$lm_id}))
                {
                    $match_type             = '99_formula';
                    $lmid_hits_hash{$lm_id} = $match_type;
                    $num_hits++;
                }
            }
        }
    }

    # substring matching, if no name matches found
    #
    # I think this is a bad idea, disable it
    #
    # this doesn't appear to be needed, currently adds zero new hits
    if (0 && $num_hits == 0)
    {
       foreach $name_lc_lmaps (@lmaps_name_lc_all_array)
       {
           if ($name_lc_lmaps =~ /\Q$name_lc_orig\E/)
           {
               @temp_array = keys %{$name_lc_lmid_hash{$name_lc_lmaps}};
               foreach $lm_id (@temp_array)
               {
                   if (!defined($lmid_hits_hash{$lm_id}))
                   {
                       $match_type             = 'substring';
                       $lmid_hits_hash{$lm_id} = $match_type;
                       $num_hits++;

                       #printf STDERR "%s  %s  %s  %s  %s\n",
                       #    $row, $name, $lm_id, $name_lc_lmaps, $match_type;
                   }
               }
           }
       }
    }

    # substring matching, if no name matches found
    #
    # I think this is a bad idea, disable it
    if (0 && $num_hits == 0)
    {
       foreach $name_lc_lmaps (@lmaps_name_lc_all_array)
       {
           if ($name_lc_lmaps =~ /\Q$name_lc_slash\E/)
           {
               @temp_array = keys %{$name_lc_lmid_hash{$name_lc_lmaps}};
               foreach $lm_id (@temp_array)
               {
                   if (!defined($lmid_hits_hash{$lm_id}))
                   {
                       $match_type = 'substring_relaxed';
                       $lmid_hits_hash{$lm_id} = $match_type;
                       $num_hits++;

                       #printf STDERR "%s  %s  %s  %s  %s\n",
                       #    $row, $name, $lm_id, $name_lc_lmaps, $match_type;
                   }
               }
           }
       }
    }
    
    if ($num_hits)
    {
        foreach $lm_id (sort keys %lmid_hits_hash)
        {
            $row_lmid_hash{$row}{$lm_id} = $lmid_hits_hash{$lm_id};
        }
    }

    $row++;
}
$num_data_rows = $row;
close DATA;


#foreach $row (sort {$a<=>$b} keys %row_lmid_hash)
#{
#    @row_lmid_array = sort keys %{$row_lmid_hash{$row}};
#
#    $match_type = $row_lmid_hash{$row}{$row_lmid_array[0]};
#    
#    printf "%s\t%s\n", $row, $match_type;
#}


$insertion_col = $data_header_col_hash{'SubClass'};
if (!defined($insertion_col))
{
    $insertion_col = $data_name_col;
}

# print new header line
# insert LipidMaps annotation after Lipid column
for ($i = 0; $i < @data_header_col_array; $i++)
{
    if ($i)
    {
        print "\t";
    }
    
    print $data_header_col_array[$i];
    
    if ($i == $insertion_col)
    {
        print "\tMatchTypeLipidMaps";
    
        for ($col = 0; $col < @lmaps_header_col_array; $col++)
        {
            $header = $lmaps_header_col_array[$col];
            print "\t$header";
        }
    }
}
print "\n";


for ($row = 0; $row < $num_data_rows; $row++)
{
    $line = $data_row_line_array[$row];
    
    @array = split /\t/, $line;
    
    # insert LipidMaps annotation after Lipid column
    for ($i = 0; $i < @data_header_col_array; $i++)
    {
        if ($i)
        {
            print "\t";
        }
        
        $value = $array[$i];
        print "$value";
        
        
        if ($i == $insertion_col)
        {
            if (defined($row_lmid_hash{$row}))
            {
                @row_lmid_array = sort keys %{$row_lmid_hash{$row}};
                $match_type     = $row_lmid_hash{$row}{$row_lmid_array[0]};
                
                print "\t$match_type";
                
                for ($col = 0; $col < @lmaps_header_col_array; $col++)
                {
                    %temp_seen_hash = ();
                    @temp_array     = ();
                    $k              = 0;
                    foreach $lm_id (@row_lmid_array)
                    {
                        $value_str = $lmid_col_annotation_hash{$lm_id}{$col};
                        if (!defined($value_str))
                        {
                            $value_str = '';
                        }

                        @sub_array = split /\s*\|+\s*/, $value_str;
                        
                        foreach $value (@sub_array)
                        {
                            if (!($value =~ /\S/))
                            {
                                next;
                            }

                            if (!defined($temp_seen_hash{$value}))
                            {
                                $temp_array[$k++] = $value;
                            }

                            $temp_seen_hash{$value} = 1;
                        }
                    }
                    
                    $merged_str = join " | ", @temp_array;

                    print "\t$merged_str";
                }
            }
            else
            {
                print "\t";
            
                for ($col = 0; $col < @lmaps_header_col_array; $col++)
                {
                    print "\t";
                }
            }

        }
    }
    print "\n";
}
