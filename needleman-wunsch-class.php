<?php
/*
 * Copyright (c) 2011 Andrew E. Bruno <aeb@qnot.org> 
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * This class implements the Needleman-Wunsch global alignment algorithm and
 * was created for educational purposes to demonstrate how the algorithm works.
 * It is not intended for use in real sequence alignment or use with very large
 * sequences. It computes the alignment scoring table and provides methods to
 * display the alignment table and optimal global alignment in both HTML
 * and ASCII format. 
 */
class NeedlemanWunsch {
    private static $arrow_up = '&#8593;';
    private static $arrow_left = '&#8592;';
    private static $arrow_nw = '&#8598;';
    private $match_score = 1;
    private $mismatch_score = 0;
    private $gap_penalty = -1;
    private $matrix = array();
    private $optimal_alignment = array();
    /**
     * Constructor
     */
    public function __construct($str1, $str2) {
			$this->str1 = $str1;
			$this->str2 = $str2;
			$this->blosum85 = Array(
				"A"=> Array("A"=>5, "R"=>2, "N"=>2, "D"=>2, "C"=>1, "Q"=>1, "E"=>1, "G"=>0, "H"=>2, "I"=>2, "L"=>2, "K"=>1, "M"=>2, "F"=>3, "P"=>1, "S"=>1, "T"=>0, "W"=>3, "Y"=>3, "V"=>1, "B"=>2, "Z"=>1, "X"=>1, "*"=>6),
				"R"=> Array("A"=>2, "R"=>6, "N"=>1, "D"=>2, "C"=>4, "Q"=>1, "E"=>1, "G"=>3, "H"=>0, "I"=>4, "L"=>3, "K"=>2, "M"=>2, "F"=>4, "P"=>2, "S"=>1, "T"=>2, "W"=>4, "Y"=>3, "V"=>3, "B"=>2, "Z"=>0, "X"=>2, "*"=>6),
				"N"=> Array("A"=>2, "R"=>1, "N"=>7, "D"=>1, "C"=>4, "Q"=>0, "E"=>1, "G"=>1, "H"=>0, "I"=>4, "L"=>4, "K"=>0, "M"=>3, "F"=>4, "P"=>3, "S"=>0, "T"=>0, "W"=>5, "Y"=>3, "V"=>4, "B"=>4, "Z"=>1, "X"=>2, "*"=>6),
				"D"=> Array("A"=>2, "R"=>2, "N"=>1, "D"=>7, "C"=>5, "Q"=>1, "E"=>1, "G"=>2, "H"=>2, "I"=>5, "L"=>5, "K"=>1, "M"=>4, "F"=>4, "P"=>2, "S"=>1, "T"=>2, "W"=>6, "Y"=>4, "V"=>4, "B"=>4, "Z"=>1, "X"=>2, "*"=>6),
				"C"=> Array("A"=>1, "R"=>4, "N"=>4, "D"=>5, "C"=>9, "Q"=>4, "E"=>5, "G"=>4, "H"=>5, "I"=>2, "L"=>2, "K"=>4, "M"=>2, "F"=>3, "P"=>4, "S"=>2, "T"=>2, "W"=>4, "Y"=>3, "V"=>1, "B"=>4, "Z"=>5, "X"=>3, "*"=>6),
				"Q"=> Array("A"=>1, "R"=>1, "N"=>0, "D"=>1, "C"=>4, "Q"=>6, "E"=>2, "G"=>3, "H"=>1, "I"=>4, "L"=>3, "K"=>1, "M"=>0, "F"=>4, "P"=>2, "S"=>1, "T"=>1, "W"=>3, "Y"=>2, "V"=>3, "B"=>1, "Z"=>4, "X"=>1, "*"=>6),
				"E"=> Array("A"=>1, "R"=>1, "N"=>1, "D"=>1, "C"=>5, "Q"=>2, "E"=>6, "G"=>3, "H"=>1, "I"=>4, "L"=>4, "K"=>0, "M"=>3, "F"=>4, "P"=>2, "S"=>1, "T"=>1, "W"=>4, "Y"=>4, "V"=>3, "B"=>0, "Z"=>4, "X"=>1, "*"=>6),
				"G"=> Array("A"=>0, "R"=>3, "N"=>1, "D"=>2, "C"=>4, "Q"=>3, "E"=>3, "G"=>6, "H"=>3, "I"=>5, "L"=>5, "K"=>2, "M"=>4, "F"=>4, "P"=>3, "S"=>1, "T"=>2, "W"=>4, "Y"=>5, "V"=>4, "B"=>1, "Z"=>3, "X"=>2, "*"=>6),
				"H"=> Array("A"=>2, "R"=>0, "N"=>0, "D"=>2, "C"=>5, "Q"=>1, "E"=>1, "G"=>3, "H"=>8, "I"=>4, "L"=>3, "K"=>1, "M"=>3, "F"=>2, "P"=>3, "S"=>1, "T"=>2, "W"=>3, "Y"=>2, "V"=>4, "B"=>1, "Z"=>0, "X"=>2, "*"=>6),
				"I"=> Array("A"=>2, "R"=>4, "N"=>4, "D"=>5, "C"=>2, "Q"=>4, "E"=>4, "G"=>5, "H"=>4, "I"=>5, "L"=>1, "K"=>3, "M"=>1, "F"=>1, "P"=>4, "S"=>3, "T"=>1, "W"=>3, "Y"=>2, "V"=>3, "B"=>5, "Z"=>4, "X"=>2, "*"=>6),
				"L"=> Array("A"=>2, "R"=>3, "N"=>4, "D"=>5, "C"=>2, "Q"=>3, "E"=>4, "G"=>5, "H"=>3, "I"=>1, "L"=>4, "K"=>3, "M"=>2, "F"=>0, "P"=>4, "S"=>3, "T"=>2, "W"=>3, "Y"=>2, "V"=>0, "B"=>5, "Z"=>4, "X"=>2, "*"=>6),
				"K"=> Array("A"=>1, "R"=>2, "N"=>0, "D"=>1, "C"=>4, "Q"=>1, "E"=>0, "G"=>2, "H"=>1, "I"=>3, "L"=>3, "K"=>6, "M"=>2, "F"=>4, "P"=>2, "S"=>1, "T"=>1, "W"=>5, "Y"=>3, "V"=>3, "B"=>1, "Z"=>1, "X"=>1, "*"=>6),
				"M"=> Array("A"=>2, "R"=>2, "N"=>3, "D"=>4, "C"=>2, "Q"=>0, "E"=>3, "G"=>4, "H"=>3, "I"=>1, "L"=>2, "K"=>2, "M"=>7, "F"=>1, "P"=>3, "S"=>2, "T"=>1, "W"=>2, "Y"=>2, "V"=>0, "B"=>4, "Z"=>2, "X"=>1, "*"=>6),
				"F"=> Array("A"=>3, "R"=>4, "N"=>4, "D"=>4, "C"=>3, "Q"=>4, "E"=>4, "G"=>4, "H"=>2, "I"=>1, "L"=>0, "K"=>4, "M"=>1, "F"=>7, "P"=>4, "S"=>3, "T"=>3, "W"=>0, "Y"=>3, "V"=>1, "B"=>4, "Z"=>4, "X"=>2, "*"=>6),
				"P"=> Array("A"=>1, "R"=>2, "N"=>3, "D"=>2, "C"=>4, "Q"=>2, "E"=>2, "G"=>3, "H"=>3, "I"=>4, "L"=>4, "K"=>2, "M"=>3, "F"=>4, "P"=>8, "S"=>1, "T"=>2, "W"=>5, "Y"=>4, "V"=>3, "B"=>3, "Z"=>2, "X"=>2, "*"=>6),
				"S"=> Array("A"=>1, "R"=>1, "N"=>0, "D"=>1, "C"=>2, "Q"=>1, "E"=>1, "G"=>1, "H"=>1, "I"=>3, "L"=>3, "K"=>1, "M"=>2, "F"=>3, "P"=>1, "S"=>5, "T"=>1, "W"=>4, "Y"=>2, "V"=>2, "B"=>0, "Z"=>1, "X"=>1, "*"=>6),
				"T"=> Array("A"=>0, "R"=>2, "N"=>0, "D"=>2, "C"=>2, "Q"=>1, "E"=>1, "G"=>2, "H"=>2, "I"=>1, "L"=>2, "K"=>1, "M"=>1, "F"=>3, "P"=>2, "S"=>1, "T"=>5, "W"=>4, "Y"=>2, "V"=>0, "B"=>1, "Z"=>1, "X"=>1, "*"=>6),
				"W"=> Array("A"=>3, "R"=>4, "N"=>5, "D"=>6, "C"=>4, "Q"=>3, "E"=>4, "G"=>4, "H"=>3, "I"=>3, "L"=>3, "K"=>5, "M"=>2, "F"=>0, "P"=>5, "S"=>4, "T"=>4, "W"=>11, "Y"=>2, "V"=>3, "B"=>5, "Z"=>4, "X"=>3, "*"=>6),
				"Y"=> Array("A"=>3, "R"=>3, "N"=>3, "D"=>4, "C"=>3, "Q"=>2, "E"=>4, "G"=>5, "H"=>2, "I"=>2, "L"=>2, "K"=>3, "M"=>2, "F"=>3, "P"=>4, "S"=>2, "T"=>2, "W"=>2, "Y"=>7, "V"=>2, "B"=>4, "Z"=>3, "X"=>2, "*"=>6),
				"V"=> Array("A"=>1, "R"=>3, "N"=>4, "D"=>4, "C"=>1, "Q"=>3, "E"=>3, "G"=>4, "H"=>4, "I"=>3, "L"=>0, "K"=>3, "M"=>0, "F"=>1, "P"=>3, "S"=>2, "T"=>0, "W"=>3, "Y"=>2, "V"=>5, "B"=>4, "Z"=>3, "X"=>1, "*"=>6),
				"B"=> Array("A"=>2, "R"=>2, "N"=>4, "D"=>4, "C"=>4, "Q"=>1, "E"=>0, "G"=>1, "H"=>1, "I"=>5, "L"=>5, "K"=>1, "M"=>4, "F"=>4, "P"=>3, "S"=>0, "T"=>1, "W"=>5, "Y"=>4, "V"=>4, "B"=>4, "Z"=>0, "X"=>2, "*"=>6),
				"Z"=> Array("A"=>1, "R"=>0, "N"=>1, "D"=>1, "C"=>5, "Q"=>4, "E"=>4, "G"=>3, "H"=>0, "I"=>4, "L"=>4, "K"=>1, "M"=>2, "F"=>4, "P"=>2, "S"=>1, "T"=>1, "W"=>4, "Y"=>3, "V"=>3, "B"=>0, "Z"=>4, "X"=>1, "*"=>6),
				"X"=> Array("A"=>1, "R"=>2, "N"=>2, "D"=>2, "C"=>3, "Q"=>1, "E"=>1, "G"=>2, "H"=>2, "I"=>2, "L"=>2, "K"=>1, "M"=>1, "F"=>2, "P"=>2, "S"=>1, "T"=>1, "W"=>3, "Y"=>2, "V"=>1, "B"=>2, "Z"=>1, "X"=>2, "*"=>6),
				"*"=> Array("A"=>6, "R"=>6, "N"=>6, "D"=>6, "C"=>6, "Q"=>6, "E"=>6, "G"=>6, "H"=>6, "I"=>6, "L"=>6, "K"=>6, "M"=>6, "F"=>6, "P"=>6, "S"=>6, "T"=>6, "W"=>6, "Y"=>6, "V"=>6, "B"=>6, "Z"=>6, "X"=>6, "*"=>1)
			);
			$this->blosum65 = Array(
				"A"=> Array("A"=>  4,"R"=> -1,"N"=> -2,"D"=> -2,"C"=>  0,"Q"=> -1,"E"=>-1,"G"=>0,"H"=>-2,"I"=>-1,"L"=>-2,"K"=>-1,"M"=>-1,"F"=>-2,"P"=>-1,"S"=>1,"T"=>0,"W"=>-3,"Y"=>-2,"V"=>0,"B"=>-2,"Z"=>-1,"X"=>-1,"*"=>-5 ),
				"R"=> Array("A"=> -1,"R"=>  6,"N"=>  0,"D"=> -2,"C"=> -4,"Q"=>  1,"E"=>  0,"G"=> -2,"H"=>  0,"I"=> -3,"L"=> -2,"K"=>  2,"M"=> -2,"F"=> -3,"P"=> -2,"S"=> -1,"T"=> -1,"W"=> -3,"Y"=> -2,"V"=> -3,"B"=> -1,"Z"=>  0,"X"=> -1,"*"=> -5 ),
				"N"=> Array("A"=> -2,"R"=>  0,"N"=>  6,"D"=>  1,"C"=> -3,"Q"=>  0,"E"=>  0,"G"=> -1,"H"=>  1,"I"=> -3,"L"=> -4,"K"=>  0,"M"=> -2,"F"=> -3,"P"=> -2,"S"=>  1,"T"=>  0,"W"=> -4,"Y"=> -2,"V"=> -3,"B"=>  3,"Z"=>  0,"X"=> -1,"*"=> -5 ),
				"D"=> Array("A"=> -2,"R"=> -2,"N"=>  1,"D"=>  6,"C"=> -4,"Q"=>  0,"E"=>  2,"G"=> -1,"H"=> -1,"I"=> -3,"L"=> -4,"K"=> -1,"M"=> -3,"F"=> -4,"P"=> -2,"S"=>  0,"T"=> -1,"W"=> -5,"Y"=> -3,"V"=> -3,"B"=>  4,"Z"=>  1,"X"=> -1,"*"=> -5 ),
				"C"=> Array("A"=>  0,"R"=> -4,"N"=> -3,"D"=> -4,"C"=>  9,"Q"=> -3,"E"=> -4,"G"=> -3,"H"=> -3,"I"=> -1,"L"=> -1,"K"=> -3,"M"=> -2,"F"=> -2,"P"=> -3,"S"=> -1,"T"=> -1,"W"=> -2,"Y"=> -2,"V"=> -1,"B"=> -3,"Z"=> -4,"X"=> -2,"*"=> -5 ),
				"Q"=> Array("A"=> -1,"R"=>  1,"N"=>  0,"D"=>  0,"C"=> -3,"Q"=>  6,"E"=>  2,"G"=> -2,"H"=>  1,"I"=> -3,"L"=> -2,"K"=>  1,"M"=>  0,"F"=> -3,"P"=> -1,"S"=>  0,"T"=> -1,"W"=> -2,"Y"=> -2,"V"=> -2,"B"=>  0,"Z"=>  3,"X"=> -1,"*"=> -5 ),
				"E"=> Array("A"=> -1,"R"=>  0,"N"=>  0,"D"=>  2,"C"=> -4,"Q"=>  2,"E"=>  5,"G"=> -2,"H"=>  0,"I"=> -3,"L"=> -3,"K"=>  1,"M"=> -2,"F"=> -3,"P"=> -1,"S"=>  0,"T"=> -1,"W"=> -3,"Y"=> -2,"V"=> -3,"B"=>  1,"Z"=>  4,"X"=> -1,"*"=> -5 ),
				"G"=> Array("A"=>  0,"R"=> -2,"N"=> -1,"D"=> -1,"C"=> -3,"Q"=> -2,"E"=> -2,"G"=>  6,"H"=> -2,"I"=> -4,"L"=> -4,"K"=> -2,"M"=> -3,"F"=> -3,"P"=> -2,"S"=>  0,"T"=> -2,"W"=> -3,"Y"=> -3,"V"=> -3,"B"=> -1,"Z"=> -2,"X"=> -2,"*"=> -5 ),
				"H"=> Array("A"=> -2,"R"=>  0,"N"=>  1,"D"=> -1,"C"=> -3,"Q"=>  1,"E"=>  0,"G"=> -2,"H"=>  8,"I"=> -3,"L"=> -3,"K"=> -1,"M"=> -2,"F"=> -1,"P"=> -2,"S"=> -1,"T"=> -2,"W"=> -2,"Y"=>  2,"V"=> -3,"B"=>  0,"Z"=>  0,"X"=> -1,"*"=> -5 ),
				"I"=> Array("A"=> -1,"R"=> -3,"N"=> -3,"D"=> -3,"C"=> -1,"Q"=> -3,"E"=> -3,"G"=> -4,"H"=> -3,"I"=>  4,"L"=>  2,"K"=> -3,"M"=>  1,"F"=>  0,"P"=> -3,"S"=> -2,"T"=> -1,"W"=> -2,"Y"=> -1,"V"=>  3,"B"=> -3,"Z"=> -3,"X"=> -1,"*"=> -5 ),
				"L"=> Array("A"=> -2,"R"=> -2,"N"=> -4,"D"=> -4,"C"=> -1,"Q"=> -2,"E"=> -3,"G"=> -4,"H"=> -3,"I"=>  2,"L"=>  4,"K"=> -3,"M"=>  2,"F"=>  0,"P"=> -3,"S"=> -3,"T"=> -1,"W"=> -2,"Y"=> -1,"V"=>  1,"B"=> -4,"Z"=> -3,"X"=> -1,"*"=> -5 ),
				"K"=> Array("A"=> -1,"R"=>  2,"N"=>  0,"D"=> -1,"C"=> -3,"Q"=>  1,"E"=>  1,"G"=> -2,"H"=> -1,"I"=> -3,"L"=> -3,"K"=>  5,"M"=> -2,"F"=> -3,"P"=> -1,"S"=>  0,"T"=> -1,"W"=> -3,"Y"=> -2,"V"=> -2,"B"=>  0,"Z"=>  1,"X"=> -1,"*"=> -5 ),
				"M"=> Array("A"=> -1,"R"=> -2,"N"=> -2,"D"=> -3,"C"=> -2,"Q"=>  0,"E"=> -2,"G"=> -3,"H"=> -2,"I"=>  1,"L"=>  2,"K"=> -2,"M"=>  6,"F"=>  0,"P"=> -3,"S"=> -2,"T"=> -1,"W"=> -2,"Y"=> -1,"V"=>  1,"B"=> -3,"Z"=> -2,"X"=> -1,"*"=> -5 ),
				"F"=> Array("A"=> -2,"R"=> -3,"N"=> -3,"D"=> -4,"C"=> -2,"Q"=> -3,"E"=> -3,"G"=> -3,"H"=> -1,"I"=>  0,"L"=>  0,"K"=> -3,"M"=>  0,"F"=>  6,"P"=> -4,"S"=> -2,"T"=> -2,"W"=>  1,"Y"=>  3,"V"=> -1,"B"=> -3,"Z"=> -3,"X"=> -2,"*"=> -5 ),
				"P"=> Array("A"=> -1,"R"=> -2,"N"=> -2,"D"=> -2,"C"=> -3,"Q"=> -1,"E"=> -1,"G"=> -2,"H"=> -2,"I"=> -3,"L"=> -3,"K"=> -1,"M"=> -3,"F"=> -4,"P"=>  8,"S"=> -1,"T"=> -1,"W"=> -4,"Y"=> -3,"V"=> -2,"B"=> -2,"Z"=> -1,"X"=> -2,"*"=> -5 ),
				"S"=> Array("A"=>  1,"R"=> -1,"N"=>  1,"D"=>  0,"C"=> -1,"Q"=>  0,"E"=>  0,"G"=>  0,"H"=> -1,"I"=> -2,"L"=> -3,"K"=>  0,"M"=> -2,"F"=> -2,"P"=> -1,"S"=>  4,"T"=>  1,"W"=> -3,"Y"=> -2,"V"=> -2,"B"=>  0,"Z"=>  0,"X"=> -1,"*"=> -5 ),
				"T"=> Array("A"=>  0,"R"=> -1,"N"=>  0,"D"=> -1,"C"=> -1,"Q"=> -1,"E"=> -1,"G"=> -2,"H"=> -2,"I"=> -1,"L"=> -1,"K"=> -1,"M"=> -1,"F"=> -2,"P"=> -1,"S"=>  1,"T"=>  5,"W"=> -3,"Y"=> -2,"V"=>  0,"B"=> -1,"Z"=> -1,"X"=> -1,"*"=> -5 ),
				"W"=> Array("A"=> -3,"R"=> -3,"N"=> -4,"D"=> -5,"C"=> -2,"Q"=> -2,"E"=> -3,"G"=> -3,"H"=> -2,"I"=> -2,"L"=> -2,"K"=> -3,"M"=> -2,"F"=>  1,"P"=> -4,"S"=> -3,"T"=> -3,"W"=> 10,"Y"=>  2,"V"=> -3,"B"=> -4,"Z"=> -3,"X"=> -2,"*"=> -5 ),
				"Y"=> Array("A"=> -2,"R"=> -2,"N"=> -2,"D"=> -3,"C"=> -2,"Q"=> -2,"E"=> -2,"G"=> -3,"H"=>  2,"I"=> -1,"L"=> -1,"K"=> -2,"M"=> -1,"F"=>  3,"P"=> -3,"S"=> -2,"T"=> -2,"W"=>  2,"Y"=>  7,"V"=> -1,"B"=> -3,"Z"=> -2,"X"=> -1,"*"=> -5 ),
				"V"=> Array("A"=>  0,"R"=> -3,"N"=> -3,"D"=> -3,"C"=> -1,"Q"=> -2,"E"=> -3,"G"=> -3,"H"=> -3,"I"=>  3,"L"=>  1,"K"=> -2,"M"=>  1,"F"=> -1,"P"=> -2,"S"=> -2,"T"=>  0,"W"=> -3,"Y"=> -1,"V"=>  4,"B"=> -3,"Z"=> -2,"X"=> -1,"*"=> -5 ),
				"B"=> Array("A"=> -2,"R"=> -1,"N"=>  3,"D"=>  4,"C"=> -3,"Q"=>  0,"E"=>  1,"G"=> -1,"H"=>  0,"I"=> -3,"L"=> -4,"K"=>  0,"M"=> -3,"F"=> -3,"P"=> -2,"S"=>  0,"T"=> -1,"W"=> -4,"Y"=> -3,"V"=> -3,"B"=>  4,"Z"=>  1,"X"=> -1,"*"=> -5 ),
				"Z"=> Array("A"=> -1,"R"=>  0,"N"=>  0,"D"=>  1,"C"=> -4,"Q"=>  3,"E"=>  4,"G"=> -2,"H"=>  0,"I"=> -3,"L"=> -3,"K"=>  1,"M"=> -2,"F"=> -3,"P"=> -1,"S"=>  0,"T"=> -1,"W"=> -3,"Y"=> -2,"V"=> -2,"B"=>  1,"Z"=>  4,"X"=> -1,"*"=> -5 ),
				"X"=> Array("A"=> -1,"R"=> -1,"N"=> -1,"D"=> -1,"C"=> -2,"Q"=> -1,"E"=> -1,"G"=> -2,"H"=> -1,"I"=> -1,"L"=> -1,"K"=> -1,"M"=> -1,"F"=> -2,"P"=> -2,"S"=> -1,"T"=> -1,"W"=> -2,"Y"=> -1,"V"=> -1,"B"=> -1,"Z"=> -1,"X"=> -1,"*"=> -5 ),
				"*"=> Array("A"=> -5,"R"=> -5,"N"=> -5,"D"=> -5,"C"=> -5,"Q"=> -5,"E"=> -5,"G"=> -5,"H"=> -5,"I"=> -5,"L"=> -5,"K"=> -5,"M"=> -5,"F"=> -5,"P"=> -5,"S"=> -5,"T"=> -5,"W"=> -5,"Y"=> -5,"V"=> -5,"B"=> -5,"Z"=> -5,"X"=> -5,"*"=>  1 )
				);
        // $this->match_score = $match_score;
        // $this->mismatch_score = $mismatch_score;
				$this->gap_penalty = -0.5;
        $this->gap_extension = -0.5;
				$this->compute();
    }
    /**
     * Computes the Needleman-Wunsch global alignment and returns a data structure
     * representing the alignment table.
     */
    public function compute() {
				$seq1 = $this->str1;
				$seq2 = $this->str2;
				$blosum65 = $this->blosum65;
        $this->init($seq1, $seq2);
        for($i = 1; $i < count($this->matrix); $i++) {
            for($j = 1; $j < count($this->matrix[$i]); $j++) {
                // $match_mismatch = ($seq1[$i-1] === $seq2[$j-1]) ? $this->match_score : $this->mismatch_score;
								$match_mismatch = $blosum65[$seq1[$i-1]][$seq2[$j-1]];
								$previous_up_pointer = $this->matrix[$i-1][$j]['pointer'];
								$previous_left_pointer = $this->matrix[$i][$j-1]['pointer'];
								if ($previous_up_pointer === self::$arrow_up) {
									$hgap_penalty = $this->gap_extension;
								} else {
									$hgap_penalty = $this->gap_penalty;
								}
								if ($previous_left_pointer === self::$arrow_left) {
									$lgap_penalty = $this->gap_extension;
								} else {
									$lgap_penalty = $this->gap_penalty;
								}
                $match = $this->matrix[$i-1][$j-1]['val'] + $match_mismatch;
								// $hgap = $this->matrix[$i-1][$j]['val'] + $this->gap_penalty;
                $hgap = $this->matrix[$i-1][$j]['val'] + $hgap_penalty;
								// $vgap = $this->matrix[$i][$j-1]['val'] + $this->gap_penalty;
                $vgap = $this->matrix[$i][$j-1]['val'] + $lgap_penalty;
                $max = max($match, $hgap, $vgap);
                $pointer = self::$arrow_nw;
                if($max === $hgap) {
                    $pointer = self::$arrow_up;
                } else if($max === $vgap) {
                    $pointer = self::$arrow_left;
                }
                $this->matrix[$i][$j]['pointer'] = $pointer;
                $this->matrix[$i][$j]['val'] = $max;
								$this->gap_penalty = -10;
            }
        }
        $i = count($this->matrix)-1;
        $j = count($this->matrix[0])-1;
        $this->optimal_alignment['seq1'] = array();
        $this->optimal_alignment['seq2'] = array();
        $this->optimal_alignment['aln'] = array();
        $this->optimal_alignment['score'] = $this->matrix[$i][$j]['val'];
        while($i !== 0 and $j !== 0) { //mudar logica de loop e de diminicao de i e j para conseguir incluir antes do primeiro match
				// while ($flag_trace === 1) {
            $base1 = $seq1[$i-1];
            $base2 = $seq2[$j-1];
            $this->matrix[$i][$j]['trace'] = true;
            $pointer = $this->matrix[$i][$j]['pointer'];
            if($pointer === self::$arrow_nw) {
                $i--;
                $j--;
                $this->optimal_alignment['seq1'][] = $base1;
                $this->optimal_alignment['seq2'][] = $base2;
                $this->optimal_alignment['aln'][] = ($base1 === $base2) ? '|' : ' ';
            } else if($pointer === self::$arrow_up) {
                $i--;
                $this->optimal_alignment['seq1'][] = $base1;
                $this->optimal_alignment['seq2'][] = '-';
                $this->optimal_alignment['aln'][] = ' ';
            } else if($pointer === self::$arrow_left) {
                $j--;
                $this->optimal_alignment['seq1'][] = '-';
                $this->optimal_alignment['seq2'][] = $base2;
                $this->optimal_alignment['aln'][] = ' ';
            } else {
                die("Invalid pointer: $i,$j");
            }
        }
				while ($i >= 0) {
					$base1 = $seq1[$i-1];
					if ($base1 !== '') {
						$this->optimal_alignment['seq1'][] = $base1;
						$this->optimal_alignment['seq2'][] = '-';
					}
					$i--;
				}
				while ($j >= 0) {
					$base2 = $seq2[$j-1];
					if ($base2 !== '') {
						$this->optimal_alignment['seq1'][] = '-';
						$this->optimal_alignment['seq2'][] = $base2;
					}
					$j--;
				}
        foreach(array('seq1', 'seq2', 'aln') as $k) {
            $this->optimal_alignment[$k] = array_reverse($this->optimal_alignment[$k]);
        }
        return $this->matrix;
    }
    /**
     * Returns the optimal alignment data structure
     */
    public function getOptimalGlobalAlignment() {
        return $this->optimal_alignment;
    }

		private function init($seq1, $seq2) {
			$this->matrix = array();
			$this->optimal_alignment = array();
			for($i = 0; $i < strlen($seq1)+1; $i++) {
				for($j = 0; $j < strlen($seq2)+1; $j++) {
					$this->matrix[$i][$j] = array(
						'pointer' => null,
						'trace' => null,
						'val' => 0
					);
				}
			}
			for($i = 0; $i < strlen($seq1); $i++) {
				$this->matrix[$i+1][0]['val'] = ($i+1) * $this->gap_penalty;
			}
			for($j = 0; $j < strlen($seq2); $j++) {
				$this->matrix[0][$j+1]['val'] = ($j+1) * $this->gap_penalty;
			}
		}
}

?>
