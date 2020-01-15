#!/usr/bin/perl -w

# TODO --
#  add trypsin exclusion support
#  add unique + razor maxquant group inclusion/exclusion criteria

# can't pipe input, since we need to read through the file twice

sub bless_delimiter_bar_no_space
{
    my $text = $_[0];

    $text =~ s/\;/\|/g;
    $text =~ s/\/\//\|/g;
    $text =~ s/,/\|/g;
#    $text =~ s/\s+/\|/g;
    $text =~ s/\|+/\|/g;
    $text =~ s/^\|//;
    $text =~ s/\|$//;
    
    return $text;
}

sub bless_delimiter_bar
{
    my $text = $_[0];

    $text =~ s/\;/\|/g;
    $text =~ s/\/\//\|/g;
    $text =~ s/,/\|/g;
    $text =~ s/\s+/\|/g;
    $text =~ s/\|+/\|/g;
    $text =~ s/^\|//;
    $text =~ s/\|$//;
    
    return $text;
}

sub bless_delimiter_space
{
    my $text = $_[0];

    $text =~ s/\;/\|/g;
    $text =~ s/,/\|/g;
    $text =~ s/\s+/\|/g;
    $text =~ s/\|+/\|/g;
    $text =~ s/^\|//;
    $text =~ s/\|$//;
    $text =~ s/\|/ /g;
    
    return $text;
}


# read in command line arguments
$num_files = 0;
$syntax_error_flag = 0;
$use_maxquant_unique_flag = 0;
$use_maxquant_razor_flag  = 0;

for ($i = 0; $i < @ARGV; $i++)
{
    $field = $ARGV[$i];

    if ($field =~ /^--/)
    {
        if ($field eq '--regrouped' ||
            $field eq '--regroup')
        {
            $use_maxquant_razor_flag  = 0;
            $use_maxquant_unique_flag = 0;
        }
        elsif ($field eq '--maxquant-razor')
        {
            $use_maxquant_razor_flag = 1;
        }
        elsif ($field eq '--maxquant-unique')
        {
            $use_maxquant_unique_flag = 1;
        }
        else
        {
            printf "ABORT -- unknown option %s\n", $field;
            $syntax_error_flag = 1;
        }
    }
    else
    {
        if ($num_files == 0)
        {
            $peptides_filename = $field;
            $num_files++;
        }
        elsif ($num_files == 1)
        {
            $summary_filename = $field;
            $num_files++;
        }
    }
}

if ($syntax_error_flag)
{
    exit(1);
}



open SUMMARY, "$summary_filename" or die "can't open file $summary_filename\n";

$line = <SUMMARY>;
$line =~ s/[\r\n]+//g;
@array = split /\t/, $line;
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    $array[$i] =~ s/\"//g;
    
    $summary_header_col_hash{$array[$i]} = $i;
}

$summary_group_col = $summary_header_col_hash{'ProteinGroup'};

if (!defined($summary_group_col))
{
    print STDERR "ABORT -- missing summary ProteinGroup column\n";
    exit(1);
}


%group_kept_hash = ();
while(defined($line=<SUMMARY>))
{
    $line =~ s/[\r\n]+//g;

    @array = split /\t/, $line;
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\"//g;
    }

    $group = $array[$summary_group_col];
    
    $group_kept_hash{$group} = 1;
}



open PEPTIDES, "$peptides_filename" or die "can't open file $peptides_filename\n";

$line = <PEPTIDES>;
$line =~ s/[\r\n]+//g;
@array = split /\t/, $line;
for ($i = 0; $i < @array; $i++)
{
    $array[$i] =~ s/^\s+//;
    $array[$i] =~ s/\s+$//;
    $array[$i] =~ s/\"//g;
    
    $header_col_hash{$array[$i]} = $i;
}
$header_line = join "\t", @array;


# new headers from new maxquant pipeline
$sequence_col        = $header_col_hash{'Sequence'};

if ($use_maxquant_unique_flag)
{
    $group_col       = $header_col_hash{'Group_Maxquant'};
}
elsif ($use_maxquant_razor_flag)
{
    $group_col       = $header_col_hash{'Group_Maxquant_Razor'};
}
else
{
    $group_col       = $header_col_hash{'Group_Regrouped'};
}

if (!defined($sequence_col))
{
    print STDERR "ABORT -- missing sequence column\n";
    exit(1);
}
if (!defined($group_col))
{
    print STDERR "ABORT -- missing group column\n";
    exit(1);
}


print "$header_line\n";

while(defined($line=<PEPTIDES>))
{
    $line =~ s/[\r\n]+//g;

    @array = split /\t/, $line;
    for ($i = 0; $i < @array; $i++)
    {
        $array[$i] =~ s/^\s+//;
        $array[$i] =~ s/\s+$//;
        $array[$i] =~ s/\"//g;
    }
    $line_new = join "\t", @array;

    $group = $array[$group_col];
    
    if (defined($group_kept_hash{$group}))
    {
        print "$line_new\n";
    }
}
