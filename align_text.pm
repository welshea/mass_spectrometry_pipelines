#!/usr/bin/perl -w

# 2021-10-26:  restructure affine gap tracebacks to be more correct
# 2021-10-25:  more affine gap traceback buxfixes
# 2021-10-22:  bugfix affine gap traceback, should only affect overlap
# 2021-10-21:  bugfix global matrix column initialization typo
# 2021-10-21:  rename glocal to gmiddle to avoid literature confusion
# 2021-10-21:  more consistent tie-breaking behavior between modes
# 2021-10-21:  update post-alignment scoring examples in comments
# 2021-06-14:  further tweaks to post-alignment scoring and comments
# 2021-06-10:  rewrite large portions to support affine gaps correctly
# 2021-06-03:  enable debug output via appending _debug to $type
# 2014-02-03:  last time I worked on it before 2021

#use strict;
use List::Util qw[min max];

use POSIX;
$BOGUS_SCORE = -DBL_MAX();


# types of alignment:
#   global      global alignment of both sequences; Needleman-Wunch
#   local       local alignment; Smith-Waterman
#   overlap     best overlap, anchor at least one end of one or both
#   gmiddle     no gap penalties at begin *and* end of *both* strings
#
#   I had previously used the term "glocal" to refer to a hybrid global-local
#   alignment.  It appears that many people use the term glocal, or
#   semi-global, as another name for overlap alignment.  Overlap alignments
#   anchor at least one end of the alignment, whereas gmiddle does not anchor
#   *any* ends of the sequences.  I now refer to it as gmiddle to avoid
#   confusion with common usage of the term glocal.  gmiddle -- find the best
#   alignment within the middle of the sequences, ignoring poorly aligned
#   ends entirely.
#
#   Global:
#     Initialize first row and col with affine gap penalties
#     Traceback begins at m,n
#     Traceback ends   at 0,0
#
#   Overlap:
#     Initialize first row and column to zero
#     Traceback begins at highest score within last row and column
#     Traceback ends at either first row/col
#
#   Local:
#     Initialize first row and column to zero
#     Negative scores are floored to zero as matrix is filled
#     Traceback begins at highest score within entire matrix
#     Traceback ends at first encountered zero score
#
#   Gmiddle:
#     Initialize first row and column to zero
#     Negative scores are floored to zero as matrix is filled, but only until
#      the first positive match along the path
#     Traceback begins at highest score within entire matrix
#     Traceback ends at aforementioned first positive match along the path
#
#   Gmiddle can be thought of as a global alignment that completely ignores
#   poorly aligned end regions (a "global" alignment of the "middle" regions),
#   or as a local alignment that can pass through negatively scoring regions
#   on its way to finding a longer, higher scoring alignment.
#
#   Usually, gmiddle and local will yield the same alignment.
#   I needed to specifically craft an example for which they are different.
#
#   Gmiddle can sometimes yield longer alignments than local due to not
#   ending the traceback when it encounters a non-positive score, if the
#   scores of the blocks on either side of the unmatched region are high
#   enough.
#
# test strings, yield different alignments for the different methods
#
#   match, mis-match, gap-open, gap-extend = [25, 11, 12, 9]
#
#
#   22AA55555CCC  333AABCCC           AA55555CCC22  AABCCC333
#
#   global:    raw score = 41         global:    raw score = 41
#     ~22AA55555CCC                        AA55555CCC~22
#     333AA~~~~BCCC                        AA~~~~BCCC333
#
#   overlap:   raw score = 54         overlap:   raw score = 54
#      22AA55555CCC                        AA55555CCC22
#      ~~AA~~~~BCCC                        AA~~~~BCCC~~
#
#   local:     raw score = 75         local:     raw score = 75
#               CCC                               CCC
#               CCC                               CCC
#
#   gmiddle:   raw score = 75         gmiddle:   raw score = 75
#        AA55555CCC                        AA55555CCC
#        AA~~~~BCCC                        AA~~~~BCCC
#
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
    my $num_positive;
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

    # default to gmiddle
    if ($type ne 'global' && $type ne 'local' && $type ne 'overlap' &&
        $type ne 'gmiddle')
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
    # real test case that gmiddle needs to handle well:
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
    $matrix[$len1][$len2]{score_best}  = 0;    # best of the 3 scores
    $matrix[$len1][$len2]{diag}{score} = 0;    # from previous aligned seq1:seq2
    $matrix[$len1][$len2]{left}{score} = 0;    # from previous gap in seq2
    $matrix[$len1][$len2]{up}{score}   = 0;    # from previous gap in seq1
    
    # fill 0,0
    # the full first row and first col are initialized further down
    $matrix[0][0]{diag}{score}        = 0;
    $matrix[0][0]{diag}{state_best}   = 'origin';
    $matrix[0][0]{diag}{is_positive}  = 0;
    $matrix[0][0]{diag}{num_positive} = 0;
    $matrix[0][0]{left}{score}        = 0;
    $matrix[0][0]{left}{state_best}   = 'origin';
    $matrix[0][0]{left}{is_positive}  = 0;
    $matrix[0][0]{left}{num_positive} = 0;
    $matrix[0][0]{up}{score}          = 0;
    $matrix[0][0]{up}{state_best}     = 'origin';
    $matrix[0][0]{up}{is_positive}    = 0;
    $matrix[0][0]{up}{num_positive}   = 0;
    $matrix[0][0]{score_best}         = 0;
    $matrix[0][0]{state_best}         = 'origin';
    $matrix[0][0]{num_positive}       = 0;


    # initialize first rows and cols

    $row = 0;
    for ($col = 1; $col <= $len2; $col++)
    {
        $pointer                  = \%{$matrix[$row][$col]};

        # traceback
        $$pointer{diag}{state_best}   = 'left';
        $$pointer{diag}{is_positive}  = 0;
        $$pointer{diag}{num_positive} = 0;
        $$pointer{left}{state_best}   = 'left';
        $$pointer{left}{is_positive}  = 0;
        $$pointer{left}{num_positive} = 0;
        $$pointer{up}{state_best}     = 'left';
        $$pointer{up}{is_positive}    = 0;
        $$pointer{up}{num_positive}   = 0;
        $$pointer{score_best}         = 0;
        $$pointer{state_best}         = 'left';
        $$pointer{num_positive}       = 0;
        
        # global, initialize gaps
        if ($type eq 'global')
        {
            $score_best = $matrix[$row][$col-1]{score_best};

            # gap open
            if ($col == 1)
            {
                $$pointer{score_best} = $score_best - $gap_open_penalty;
            }
            # gap extend
            else
            {
                $$pointer{score_best} = $score_best - $gap_extend_penalty;
            }
        }

        $$pointer{diag}{score} = $$pointer{score_best};
        $$pointer{left}{score} = $$pointer{score_best};
        $$pointer{up}{score}   = $$pointer{score_best};
        
        # bogify impossible paths through global alignment
        # only needed for global, messes up non-global
        # non-global is all zeroes and can start from "illegal" transitions
        if ($type eq 'global')
        {
            $$pointer{diag}{score} = $BOGUS_SCORE;
            $$pointer{up}{score}   = $BOGUS_SCORE;
        }
    }

    $col = 0;
    for ($row = 1; $row <= $len1; $row++)
    {
        $pointer                 = \%{$matrix[$row][$col]};

        # traceback
        $$pointer{diag}{state_best}   = 'up';
        $$pointer{diag}{is_positive}  = 0;
        $$pointer{diag}{num_positive} = 0;
        $$pointer{left}{state_best}   = 'up';
        $$pointer{left}{is_positive}  = 0;
        $$pointer{left}{num_positive} = 0;
        $$pointer{up}{state_best}     = 'up';
        $$pointer{up}{is_positive}    = 0;
        $$pointer{up}{num_positive}   = 0;
        $$pointer{score_best}         = 0;
        $$pointer{state_best}         = 'up';
        $$pointer{num_positive}       = 0;

        # global, initialize gaps
        if ($type eq 'global')
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
        
        $$pointer{diag}{score} = $$pointer{score_best};
        $$pointer{left}{score} = $$pointer{score_best};
        $$pointer{up}{score}   = $$pointer{score_best};

        # bogify impossible paths through global alignment
        # only needed for global, messes up non-global
        # non-global is all zeroes and can start from "illegal" transitions
        if ($type eq 'global')
        {
            $$pointer{diag}{score} = $BOGUS_SCORE;
            $$pointer{left}{score} = $BOGUS_SCORE;
        }
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
            $pointer_tb       = \%{$matrix[$tb_row][$tb_col]};

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

            $score_diag = $$pointer_tb{diag}{score} + $score;
            $score_left = $$pointer_tb{left}{score} + $score;
            $score_up   = $$pointer_tb{up}{score}   + $score;

            
            # maximum of potential states
            # break ties on num_positive
            $num_positive = $$pointer_tb{diag}{num_positive};
            $score_best   = $score_diag;
            $state_best   = 'diag';
            if ($score_left >= $score_best)
            {
                if ($score_left > $score_best ||
                    $$pointer_tb{left}{num_positive} > $num_positive)
                {
                    $num_positive = $$pointer_tb{left}{num_positive};
                    $score_best   = $score_left;
                    $state_best   = 'left';
                }
            }
            if ($score_up >= $score_best)
            {
                if ($score_up > $score_best ||
                    $$pointer_tb{up}{num_positive} > $num_positive)
                {
                    $num_positive = $$pointer_tb{up}{num_positive};
                    $score_best   = $score_up;
                    $state_best   = 'up';
                }
            }
            $$pointer{diag}{score}      = $score_best;
            $$pointer{diag}{state_best} = $state_best;


            # increment / propagate num positive matches so far
            $$pointer{diag}{num_positive} =
               $$pointer_tb{$state_best}{num_positive} + $is_positive;


            # local doesn't allow negative scores anywhere
            if ($type eq 'local')
            {
                if ($$pointer{diag}{score} < 0)
                {
                    $$pointer{diag}{score} = 0;
                }
            }
            # gmiddle doesn't allow negative scores before first positive hit
            if ($type eq 'gmiddle' && $$pointer{diag}{num_positive} == 0)
            {
                if ($$pointer{diag}{score} < 0)
                {
                    $$pointer{diag}{score} = 0;
                }
            }



            # left: [row][col-1]
            $tb_row = $row;
            $tb_col = $col - 1;
            $pointer_tb       = \%{$matrix[$tb_row][$tb_col]};

            $score_diag     = $$pointer_tb{diag}{score} -
                              $gap_open_penalty;
            if ($tb_col == 0)
            {
                # col 0 can't have come from a gap
                $score_left = $$pointer_tb{left}{score} -
                              $gap_open_penalty;
            }
            else
            {
                # gap extension
                $score_left = $$pointer_tb{left}{score} -
                              $gap_extend_penalty;
            }


            # maximum of potential states
            # break ties on num_positive
            $num_positive = $$pointer_tb{diag}{num_positive};
            $score_best   = $score_diag;
            $state_best   = 'diag';
            if ($score_left >= $score_best)
            {
                if ($score_left > $score_best ||
                    $$pointer_tb{left}{num_positive} > $num_positive)
                {
                    $num_positive = $$pointer_tb{left}{num_positive};
                    $score_best   = $score_left;
                    $state_best   = 'left';
                }
            }
            $$pointer{left}{score}      = $score_best;
            $$pointer{left}{state_best} = $state_best;


            # increment / propagate num positive matches so far
            $$pointer{left}{num_positive} =
               $$pointer_tb{$state_best}{num_positive} + $is_positive;


            # local doesn't allow negative scores anywhere
            if ($type eq 'local')
            {
                if ($$pointer{left}{score} < 0)
                {
                    $$pointer{left}{score} = 0;
                }
            }
            # gmiddle doesn't allow negative scores before first positive hit
            if ($type eq 'gmiddle' && $$pointer{left}{num_positive} == 0)
            {
                if ($$pointer{left}{score} < 0)
                {
                    $$pointer{left}{score} = 0;
                }
            }


            
            # up: [row-1][col]
            $tb_row = $row - 1;
            $tb_col = $col;
            $pointer_tb       = \%{$matrix[$tb_row][$tb_col]};

            $score_diag   = $$pointer_tb{diag}{score} -
                            $gap_open_penalty;
            if ($tb_row == 0)
            {
                # row 0 can't have come from a gap
                $score_up = $$pointer_tb{up}{score} -
                            $gap_open_penalty;
            }
            else
            {
                # gap extension
                $score_up = $$pointer_tb{up}{score} -
                            $gap_extend_penalty;
            }


            # maximum of potential states
            # break ties on num_positive
            $num_positive = $$pointer_tb{diag}{num_positive};
            $score_best   = $score_diag;
            $state_best   = 'diag';
            if ($score_up >= $score_best)
            {
                if ($score_up > $score_best ||
                    $$pointer_tb{up}{num_positive} > $num_positive)
                {
                    $num_positive = $$pointer_tb{up}{num_positive};
                    $score_best   = $score_up;
                    $state_best   = 'up';
                }
            }
            $$pointer{up}{score}      = $score_best;
            $$pointer{up}{state_best} = $state_best;


            # increment / propagate num positive matches so far
            $$pointer{up}{num_positive} =
               $$pointer_tb{$state_best}{num_positive} + $is_positive;


            # local doesn't allow negative scores anywhere
            if ($type eq 'local')
            {
                if ($$pointer{diag}{score} < 0)
                {
                    $$pointer{diag}{score} = 0;
                }
            }
            # gmiddle doesn't allow negative scores before first positive hit
            if ($type eq 'gmiddle' && $$pointer{up}{num_positive} == 0)
            {
                if ($$pointer{up}{score} < 0)
                {
                    $$pointer{up}{score} = 0;
                }
            }



            # store best of the current 3 states
            # break ties on num_positive
            $$pointer{num_positive} = $$pointer{diag}{num_positive};
            $$pointer{score_best}   = $$pointer{diag}{score};
            $$pointer{state_best}   = 'diag';
            if ($$pointer{left}{score} >= $$pointer{score_best})
            {
                if ($$pointer{left}{score}        > $$pointer{score_best} ||
                    $$pointer{left}{num_positive} > $$pointer{num_positive})
                {
                    $$pointer{num_positive} = $$pointer{left}{num_positive};
                    $$pointer{score_best}   = $$pointer{left}{score};
                    $$pointer{state_best}   = 'left';
                }
            }
            if ($$pointer{up}{score} >= $$pointer{score_best})
            {
                if ($$pointer{up}{score}        > $$pointer{score_best} ||
                    $$pointer{up}{num_positive} > $$pointer{num_positive})
                {
                    $$pointer{num_positive} = $$pointer{up}{num_positive};
                    $$pointer{score_best}   = $$pointer{up}{score};
                    $$pointer{state_best}   = 'up';
                }
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
                if ($score_best > $best_tb_score &&
                    ($row == $len1 || $col == $len2))
                {
                    $best_tb_score = $score_best;
                    $best_tb_row   = $row;
                    $best_tb_col   = $col;
                }
            }
            # best middle alignment, whether anchored or not
            elsif ($type eq 'gmiddle')
            {
                # best score
                if ($score_best > $best_tb_score)
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
    $state_best     = $matrix[$row][$col]{state_best};
    $right_overhang = max($len1 - $row, $len2 - $col);


    # trace back through the matrix to assemble the alignment
    $pos = 0;
    while ($row || $col)
    {
        # local traceback ends once alignment score is zero
        if ($type eq 'local' && $matrix[$row][$col]{$state_best}{score} <= 0)
        {
            last;
        }

        # non-global traceback ends when first row/col is reached
        if ($type ne 'global' &&
            ($row == 0 || $col == 0))
        {
            last;
        }

        # gmiddle traceback ends at the last (first in sequence) match
        if ($type eq 'gmiddle' &&
            $matrix[$row][$col]{$state_best}{num_positive} == 0)
        {
            last;
        }
        
        $tb_row = $row;
        $tb_col = $col;

        # commented out bounding sanity checks
        # current code shouldn't be taking any impossible paths...
        if ($state_best eq 'diag')
        {
            #if ($row)
            #{
                $tb_row = $row - 1;
            #}

            #if ($col)
            #{
                $tb_col = $col - 1;
            #}
        }
        elsif ($state_best eq 'left')
        {
            #if ($col)
            #{
                $tb_col = $col - 1;
            #}
        }
        elsif ($state_best eq 'up')
        {
            #if ($row)
            #{
                $tb_row = $row - 1;
            #}
        }
        

        # gap in seq1
        if ($row == $tb_row)
        {
          $seq_align1 .= '~';
          $seq_align2 .= $char_array2[$tb_col];
          $num_gap1++;
        }
        # gap in seq2
        elsif ($col == $tb_col)
        {
          $seq_align1 .= $char_array1[$tb_row];
          $seq_align2 .= '~';
          $num_gap2++;
        }
        # match state
        else
        {
          $seq_align1 .= $char_array1[$tb_row];
          $seq_align2 .= $char_array2[$tb_col];
          
          if ($char_array1[$tb_row] eq $char_array2[$tb_col])
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


        #printf "DEBUG\t%d\t%d\t%s\t%d\t%f\t%d\t%d\n",
        #    $row, $col, $state_best,
        #    $matrix[$row][$col]{$state_best}{num_positive},
        #    $matrix[$row][$col]{score_best},
        #    $tb_row, $tb_col;

        $state_best = $matrix[$row][$col]{$state_best}{state_best};
        
        $row = $tb_row;
        $col = $tb_col;
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
    $left_overhang   = max($first_match_row, $first_match_col) - 1;
    $right_overhang  = max($len1 - $last_match_row, $len2 - $last_match_col);
    $align_length    = $last_match_pos - $first_match_pos + 1;
    

    # match score and penalties optimized for use with gmiddle alignments
    #
    # penalize less-end-anchored overhang        in numerator
    # penalize alignment length + sqrt(overhang) in denominator
    #
    # gmiddle:
    #   
    #       AAAAAccee       0.714286
    #       AAAAAggjj
    #
    #     ccAAAAAee         0.428571
    #     ggAAAAAjj
    #
    #       AAAAAccee       0.402712
    #     ggAAAAAjj
    #
    #       AAAAAccee       0.127740
    #   ggjjAAAAA
    #
    #
    # potentially problematic examples:
    #   ABCD1         ABCD1-DEFGHI2       0.653958
    #   cysgly        lcysteinylglycine   0.351221
    #   fumarate13c4  succinated4         0.0  (global/overlap align too much)
    #   fumarate      succinate           0.550510
    #     [check for < 50% coverage to zero it after returning]
    #
    # I may want to consider zeroing out the score at < 50% coverage here,
    # rather than relying on the user to do so in the calling function?
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
