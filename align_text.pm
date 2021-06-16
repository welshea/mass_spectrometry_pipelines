#!/usr/bin/perl -w

# 2021-06-10:  rewrite large portions to support affine gaps correctly
# 2021-06-03:  enable debug output via appending _debug to $type
# 2014-02-03:  last time I worked on it before 2021

#use strict;
use List::Util qw[min max];

# types of alignment:
#   global      global alignment of both sequences; Needleman-Wunch
#   local       local alignment; Smith-Waterman
#   overlap     best overlap, mostly anchor at least one end of one or both
#   glocal      best overlap, no anchor
#
#   Usually, glocal and local will yield the same alignment.
#   I needed to specifically craft an example for which they are different.
#
#   Glocal can sometimes yield longer alignments than local due to not
#   ending the traceback when it encounters a non-positive score, if the
#   scores of the blocks on either side of the unmatched region are high
#   enough.

# test strings, yield different alignments for the different methods
#
#   match, mis-match, gap-open, gap-extend = [25, 11, 12, 9]
#
#   xxxzdaaaaaaa0000000000bbbbbyyycz  qdccaaaa4444bbbbbzzzzcqqqqq
#
#   global:    raw score = 7
#     XXXZDAAAAAAA0000000000BBBBB~YYYC~~~~Z
#     ~~~QD~CCAAAA~~~~~~4444BBBBBZZZZCQQQQQ
#
#   overlap:   raw score = 93
#     ZDAAAAAAA0000000000BBBBBYYYCZ
#     QD~CCAAAA~~~~~~4444BBBBB~~~~Z
#
#   glocal:    raw score = 124
#     AAAA0000000000BBBBB
#     AAAA~~~~~~4444BBBBB
#
#   local:     raw score = 125
#     BBBBB
#     BBBBB
#
# -------------------------------------------------------------------------
#
#   cccccaaaaaddddd  eeeeeaaaaafffff
#
#   global:     raw score = 15
#     CCCCCAAAAADDDDD
#     EEEEEAAAAAFFFFF
#
#   overlap:    raw score = 29    (overhanging sequences shown for clarity)
#     CCCCC~~~~AAAAADDDDD
#         EEEEEAAAAA~~~~~DDDDD
#
#   local
#   glocal:     raw score = 125   (overhanging sequences omitted)
#     AAAAA
#     AAAAA
#
sub score_substring_mismatch
{
    my $string1            = $_[0];
    my $string2            = $_[1];
    my $type_orig          = $_[2];
    my $return_frac_id_ptr = $_[3];    # passed by reference
    my $match_score        = $_[4];
    my $mismatch_penalty   = $_[5];
    my $gap_open_penalty   = $_[6];
    my $gap_extend_penalty = $_[7];
    my @char_array1        = ();
    my @char_array2        = ();
    my @matrix             = ();
    my $len1               = length $string1;
    my $len2               = length $string2;
    my $min_len;
    my $max_len;
    my $row;
    my $col;
    my $pointer;
    my $score;
    my $score_best;
    my $score_diag;
    my $score_up;
    my $score_left;
    my $tb_row;
    my $tb_col;
    my $best_tb_row        = 0;
    my $best_tb_col        = 0;
    my $best_tb_score      = -9E99;
    my $was_positive;
    my $is_positive;
    my $seq_align1;
    my $seq_align2;

    my $num_match          = 0;
    my $num_mismatch       = 0;
    my $num_gap1           = 0;
    my $num_gap2           = 0;
    my $left_overhang      = 0;
    my $right_overhang     = 0;
    my $align_length       = 0;
    my $my_score           = 0;
    my $min_padded_len     = 0;
    my $max_padded_len     = 0;

    my $first_match_pos    =  9E99;
    my $last_match_pos     = -9E99;
    my $first_match_row    = 0;
    my $first_match_col    = 0;
    my $last_match_row     = 0;
    my $last_match_col     = 0;
    my $pos;

    # HACK -- enable debug output by appending _debug to type
    $type = $type_orig;
    $type =~ s/_debug$//;

    # default to glocal
    if ($type ne 'global' && $type ne 'local' && $type ne 'overlap' &&
        $type ne 'glocal')
    {
        printf STDERR "ERROR - no method for alignment of type %s\n", $type;
        return '';
    }
    
    # |open + extend| *MUST* be > |mismatch|
    #   otherwise, unexpected things may happen
    #
    # [25, 12, 11, 10]  seems to work well ??
    # [25, 11, 12,  9]  an alternative to discourage gap opening
    #
    # real test case that glocal needs to handle well:
    #   cysgly  lcysteinylglycine
    #
    # contrived variant, works with [25, 11, 12, 9]
    #   cysXXgly  lcysteinylgXlycine
    #
    if (!defined($match_score))
    {
        $match_score        =  25;    # > 2x largest penalty
    }
    if (!defined($mismatch_penalty))
    {
        $mismatch_penalty   =  11;
    }
    if (!defined($gap_open_penalty))
    {
        $gap_open_penalty   =  12;    # |open + extend| < 2x |mismatch|
    }
    if (!defined($gap_extend_penalty))
    {
        $gap_extend_penalty =   9;    # |open + extend| < 2x |mismatch|
    }

    $min_len = $len1;
    if ($len2 < $min_len) { $min_len = $len2 };
    
    $max_len = $len1;
    if ($len2 > $max_len) { $max_len = $len2 };

    @char_array1 = split //, uc $string1, -1;
    @char_array2 = split //, uc $string2, -1;
    
    # allocate max dimensions
    $matrix[$len1][$len2]{score_best} = 0;    # best of the 3 scores
    $matrix[$len1][$len2]{score_diag} = 0;    # from previous aligned seq1:seq2
    $matrix[$len1][$len2]{score_up}   = 0;    # from previous gap in seq1
    $matrix[$len1][$len2]{score_left} = 0;    # from previous gap in seq2
    
    # fill 0,0
    # the full first row and first col are initialized further down
    $matrix[0][0]{score_best}   = 0;
    $matrix[0][0]{score_diag}   = 0;
    $matrix[0][0]{score_up}     = 0;
    $matrix[0][0]{score_left}   = 0;
    $matrix[0][0]{is_positive}  = 0;
    $matrix[0][0]{was_positive} = 0;


    # initialize first rows and cols

    $row = 0;
    for ($col = 1; $col <= $len2; $col++)
    {
        $pointer                  = \%{$matrix[$row][$col]};

        # traceback
        $$pointer{tb_row}         = $row;
        $$pointer{tb_col}         = $col - 1;
        $$pointer{is_positive}    = 0;
        $$pointer{was_positive}   = 0;
        
        # local alignments, set first rows/cols to zero
        if ($type ne 'global')
        {
            $$pointer{score_best} = 0;
        }
        # global, initialize gaps
        else
        {
            $score_best = $matrix[$row][$col-1]{score_best};

            # gap open
            if ($row == 1)
            {
                $$pointer{score_best} = $score_best - $gap_open_penalty;
            }
            # gap extend
            else
            {
                $$pointer{score_best} = $score_best - $gap_extend_penalty;
            }
        }

        # there are no previous paths, so all paths are the same value
        $$pointer{score_diag} = $$pointer{score_best};
        $$pointer{score_up}   = $$pointer{score_best};
        $$pointer{score_left} = $$pointer{score_best};
    }

    $col = 0;
    for ($row = 1; $row <= $len1; $row++)
    {
        $pointer                 = \%{$matrix[$row][$col]};

        # traceback
        $$pointer{tb_row}        = $row - 1;
        $$pointer{tb_col}        = $col;
        $$pointer{is_positive}   = 0;
        $$pointer{was_positive}  = 0;

        # local alignments, set first rows/cols to zero
        if ($type ne 'global')
        {
            $$pointer{score_best} = 0;
        }
        # global, initialize gaps
        else
        {
            $score_best = $matrix[$row-1][$col]{score_best};

            # gap open
            if ($row == 1)
            {
                $$pointer{score_best} = $score_best - $gap_open_penalty;
            }
            # gap extend
            else
            {
                $$pointer{score_best} = $score_best - $gap_extend_penalty;
            }
        }

        # there are no previous paths, so all paths are the same value
        $$pointer{score_diag} = $$pointer{score_best};
        $$pointer{score_up}   = $$pointer{score_best};
        $$pointer{score_left} = $$pointer{score_best};
    }
    

    # fill the rest of the matrix
    for ($row = 1; $row <= $len1; $row++)
    {
        for ($col = 1; $col <= $len2; $col++)
        {
            $pointer          = \%{$matrix[$row][$col]};



            # diag: [row - 1][col - 1]
            $tb_row = $row - 1;
            $tb_col = $col - 1;

            # score if position is aligned instead of gapped
            $is_positive      = 0;
            if ($char_array1[$row-1] eq $char_array2[$col-1])
            {
                $score        = $match_score;
                $is_positive  = 1;
            }
            else
            {
                $score        = -$mismatch_penalty;
            }

            # check to see if we have ever had a positive match
            $was_positive = 0;
            if ($is_positive)
            {
                $was_positive = 1;
            }
            # else propagate the flag
            else
            {
                $was_positive = $matrix[$tb_row][$tb_col]{was_positive};
            }

            $score_diag = $matrix[$tb_row][$tb_col]{score_diag} + $score;
            $score_up   = $matrix[$tb_row][$tb_col]{score_up}   + $score;
            $score_left = $matrix[$tb_row][$tb_col]{score_left} + $score;

            # glocal doesn't allow negative scores before first positive hit
            # local  doesn't allow negative scores anywhere
            if ($type eq 'local' ||
                ($type eq 'glocal' && $was_positive == 0))
            {
                if ($score_diag < 0) { $score_diag = 0; }
                if ($score_up   < 0) { $score_up   = 0; }
                if ($score_left < 0) { $score_left = 0; }
            }
            
            # maximum of potential states
            $score_best = $score_diag;
            if ($score_up > $score_best)
            {
                $score_best = $score_up;
            }
            if ($score_left > $score_best)
            {
                $score_best = $score_left;
            }
            $$pointer{score_diag} = $score_best;

            # printf "%d\t%d\t%d\t%d\t%f\n",
            #    $row, $col, $tb_row, $tb_col, $score_best;

            # initialize best of the 3 paths
            $$pointer{score_best}   = $score_best;
            $$pointer{is_positive}  = $is_positive;
            $$pointer{was_positive} = $was_positive;
            $$pointer{tb_col}       = $tb_col;
            $$pointer{tb_row}       = $tb_row;



            # left: [row][col-1]
            $tb_row = $row;
            $tb_col = $col - 1;

            # check to see if we have ever had a positive match before
            $was_positive   = $matrix[$tb_row][$tb_col]{was_positive};

            $score_diag     = $matrix[$tb_row][$tb_col]{score_diag} -
                              $gap_open_penalty;
            if ($tb_col == 0)
            {
                # col 0 can't have come from a gap
                $score_left = $matrix[$tb_row][$tb_col]{score_left} -
                              $gap_open_penalty;
            }
            else
            {
                # gap extension
                $score_left = $matrix[$tb_row][$tb_col]{score_left} -
                              $gap_extend_penalty;
            }

            # glocal doesn't allow negative scores before first positive hit
            # local  doesn't allow negative scores anywhere
            if ($type eq 'local' ||
                ($type eq 'glocal' && $was_positive == 0))
            {
                if ($score_diag < 0) { $score_diag = 0; }
                if ($score_left < 0) { $score_left = 0; }
            }
            
            # maximum of potential states
            $score_best = $score_diag;
            if ($score_left > $score_best)
            {
                $score_best = $score_left;
            }
            $$pointer{score_left} = $score_best;

            # printf "%d\t%d\t%d\t%d\t%f\n",
            #     $row, $col, $tb_row, $tb_col, $score_best;
            
            # store best of the 3 paths
            if ($score_best > $$pointer{score_best})
            {
                $$pointer{score_best}   = $score_best;
                $$pointer{is_positive}  = 0;              # it's a gap
                $$pointer{was_positive} = $was_positive;
                $$pointer{tb_col}       = $tb_col;
                $$pointer{tb_row}       = $tb_row;
            }


            
            # up: [row-1][col]
            $tb_row = $row - 1;
            $tb_col = $col;

            # check to see if we have ever had a positive match before
            $was_positive = $matrix[$tb_row][$tb_col]{was_positive};

            $score_diag   = $matrix[$tb_row][$tb_col]{score_diag} -
                            $gap_open_penalty;
            if ($tb_row == 0)
            {
                # row 0 can't have come from a gap
                $score_up = $matrix[$tb_row][$tb_col]{score_up} -
                            $gap_open_penalty;
            }
            else
            {
                # gap extension
                $score_up = $matrix[$tb_row][$tb_col]{score_up} -
                            $gap_extend_penalty;
            }

            # glocal doesn't allow negative scores before first positive hit
            # local  doesn't allow negative scores anywhere
            if ($type eq 'local' ||
                ($type eq 'glocal' && $was_positive == 0))
            {
                if ($score_diag < 0) { $score_diag = 0; }
                if ($score_up   < 0) { $score_up   = 0; }
            }
            
            # maximum of potential states
            $score_best = $score_diag;
            if ($score_up > $score_best)
            {
                $score_best = $score_up;
            }
            $$pointer{score_up} = $score_best;

            # printf "%d\t%d\t%d\t%d\t%f\n",
            #     $row, $col, $tb_row, $tb_col, $score_best;
            
            # store best of the 3 paths
            if ($score_best > $$pointer{score_best})
            {
                $$pointer{score_best}   = $score_best;
                $$pointer{is_positive}  = 0;              # it's a gap
                $$pointer{was_positive} = $was_positive;
                $$pointer{tb_col}       = $tb_col;
                $$pointer{tb_row}       = $tb_row;
            }



            # keep track of best non-global alignments

            # maximize local score
            $score_best = $$pointer{score_best};
            if ($type eq 'local')
            {
                if ($score_best > $best_tb_score)
                {
                    $best_tb_score = $score_best;
                    $best_tb_row   = $row;
                    $best_tb_col   = $col;
                }
            }
            # anchor at least the far end of either sequence
            # gaps may still wind up causing it to be not-quite anchored
            elsif ($type eq 'overlap')
            {
                # best score, must be on far edges
                if ($score_best >= $best_tb_score &&
                    ($row == $len1 || $col == $len2))
                {
                    $best_tb_score = $score_best;
                    $best_tb_row   = $row;
                    $best_tb_col   = $col;
                }
            }
            # best overlap, whether anchored or not
            # alignment must end on a match
            elsif ($type eq 'glocal')
            {
                # best score, nearest to either end for ties
                if ($score_best >= $best_tb_score &&
                    $$pointer{is_positive})
                {
                    $best_tb_score = $score_best;
                    $best_tb_row   = $row;
                    $best_tb_col   = $col;
                }
            }
            
            # printf "%d\t%d\t\t\t%f\n", $row, $col, $$pointer{score_best};
        }
    }


    
    # traceback

    # always start in the far corner
    if ($type eq 'global')
    {
        $row        = $len1;
        $col        = $len2;
    }
    # start from the highest scoring position
    else
    {
        $row        = $best_tb_row;
        $col        = $best_tb_col;
    }
    
    $seq_align1     = '';
    $seq_align2     = '';
    $score          = $matrix[$row][$col]{score_best};
    $tb_row         = $matrix[$row][$col]{tb_row};
    $tb_col         = $matrix[$row][$col]{tb_col};
    $right_overhang = max($len1 - $row, $len2 - $col);


    # trace back through the matrix to assemble the alignment
    $pos = 0;
    while (defined($tb_row) && defined($tb_col))
    {
        # local traceback ends once alignment score is zero
        if ($type eq 'local' && $matrix[$row][$col]{score_best} <= 0)
        {
            last;
        }

        # non-global traceback ends when first row/col is reached
        if ($type ne 'global' && ($row == 0 || $col == 0))
        {
            last;
        }

        # glocal traceback ends at the last (first in sequence) match
        if ($type eq 'glocal' && $matrix[$row][$col]{was_positive} == 0)
        {
            last;
        }

        # printf "SCORE %s\n", $matrix[$row][$col]{score_best};
       
        # gap in seq1
        if ($row == $tb_row)
        {
          $seq_align1 .= '~';
          $seq_align2 .= $char_array2[$col-1];
          $num_gap1++;
        }
        # gap in seq2
        elsif ($col == $tb_col)
        {
          $seq_align1 .= $char_array1[$row-1];
          $seq_align2 .= '~';
          $num_gap2++;
        }
        # match state
        else
        {
          $seq_align1 .= $char_array1[$row-1];
          $seq_align2 .= $char_array2[$col-1];
          
          if ($char_array1[$row-1] eq $char_array2[$col-1])
          {
              $num_match++;
              
              # store first and last match positions for later
              if ($pos < $first_match_pos)
              {
                  $first_match_pos = $pos;
                  $first_match_row = $row;
                  $first_match_col = $col;
              }
              if ($pos > $last_match_pos)
              {
                  $last_match_pos  = $pos;
                  $last_match_row  = $row;
                  $last_match_col  = $col;
              }
          }
          else
          {
              $num_mismatch++;
          }
        }
        $pos--;
        
        $row = $tb_row;
        $col = $tb_col;

        $tb_row = $matrix[$row][$col]{tb_row};
        $tb_col = $matrix[$row][$col]{tb_col};
    }

    $align_length     = length $seq_align1;
    $first_match_pos += $align_length;
    $last_match_pos  += $align_length;
    
    if ($first_match_pos == 9E99 || $last_match_pos == -9E99)
    {
        $first_match_pos = $align_length;
        $last_match_pos  = $align_length;
    }
    
    $left_overhang = max($row, $col);
    
    # alignment was assembled in reverse, so un-reverse the alignment
    $seq_align1 = reverse $seq_align1;
    $seq_align2 = reverse $seq_align2;
    
    $min_padded_len = $len1 + $num_gap1;
    if ($len2 + $num_gap2 < $min_padded_len)
    {
        $min_padded_len = $len2 + $num_gap2;
    }

    $max_padded_len = $len1 + $num_gap1;
    if ($len2 + $num_gap2 > $max_padded_len)
    {
        $max_padded_len = $len2 + $num_gap2;
    }
    
    # adjust metrics so as to treat global/overlap more like local
    $num_mismatch    = $last_match_pos - $first_match_pos + 1 - $num_match -
                       max($num_gap1, $num_gap2);
    $left_overhang   = max($first_match_row, $first_match_col) - 1;
    $right_overhang  = max($len1 - $last_match_row, $len2 - $last_match_col);
    $align_length    = $last_match_pos - $first_match_pos + 1;
    

    # match score and penalties optimized for use with glocal alignments
    #
    # penalize less-end-anchored overhang        in numerator
    # penalize alignment length + sqrt(overhang) in denominator
    #
    # glocal:
    #   AAAAAcccccddddd             0.612574
    #   AAAAAeeeeefffff
    #
    #   cccccAAAAAddddd             0.338623
    #   eeeeeAAAAAfffff
    #
    #        AAAAAcccccddddd        0.311500
    #   eeeeeAAAAAfffff
    #
    #             AAAAAcccccddddd   0.194014
    #   eeeeefffffAAAAA
    #
    #
    # potentially problematic examples:
    #   ABCD1         ABCD1-DEFGHI2
    #   cysgly        lcysteinylglycine
    #   fumarate13c4  succinated4
    #   fumarate      succinate   (check for < 50% coverage after returning)
    #
    $my_score = ($num_match - min($left_overhang, $right_overhang)) /
                ($align_length + sqrt($left_overhang + $right_overhang));


    # scores < 0 aren't good for comparing, since they can behave oddly
    if ($my_score < 0)
    {
        $my_score = 0;
    }
    
    if ($type_orig =~ /_debug$/)
    {
        printf STDERR "%s\n",             $seq_align1;
        printf STDERR "%s\n",             $seq_align2;
        printf STDERR "M:\t%d\n",         $num_match;
        printf STDERR "MM:\t%d\n",        $num_mismatch;
        printf STDERR "Gap1:\t%d\n",      $num_gap1;
        printf STDERR "Gap2:\t%d\n",      $num_gap2;
        printf STDERR "AlignLen:\t%s\n",  $align_length;
        printf STDERR "MinPadLen:\t%d\n", $min_padded_len;
        printf STDERR "MaxPadLen:\t%d\n", $max_padded_len;
        printf STDERR "FirstMatch:\t%s\n", $first_match_pos;
        printf STDERR "LastMatch:\t%s\n",  $last_match_pos;
        printf STDERR "LeftOH:\t%d\n",    $left_overhang;
        printf STDERR "RightOH:\t%d\n",   $right_overhang;
        printf STDERR "Score:\t%s\n",     $score;
        printf STDERR "MyScore:\t%s\n",   $my_score;
        printf STDERR "%s\n",             $seq_align1;
        printf STDERR "%s\n",             $seq_align2;
    }
    
    ${$return_frac_id_ptr} = $num_match / $min_len;
    
    return $my_score;
}
1;