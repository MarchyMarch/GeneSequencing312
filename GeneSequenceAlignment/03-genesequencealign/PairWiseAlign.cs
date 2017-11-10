using System;
using System.Collections.Generic;
using System.Text;

namespace GeneticsLab
{
    enum directions { LEFT, TOP, DIAGONAL};
    class PairWiseAlign
    {
        int MaxCharactersToAlign;
        int distance;

        public PairWiseAlign()
        {
            // Default is to align only 5000 characters in each sequence.
            this.MaxCharactersToAlign = 5000;
            this.distance = 3;
        }

        public PairWiseAlign(int len)
        {
            // Alternatively, we can use an different length; typically used with the banded option checked.
            this.MaxCharactersToAlign = len;
            this.distance = 3;
        }

        public int getMaxCharactersToAlign()
        {
            return this.MaxCharactersToAlign;
        }



        /// <summary>
        /// this is the function you implement.
        /// </summary>
        /// <param name="sequenceA">the first sequence</param>
        /// <param name="sequenceB">the second sequence, may have length not equal to the length of the first seq.</param>
        /// <param name="banded">true if alignment should be band limited.</param>
        /// <returns>the alignment score and the alignment (in a Result object) for sequenceA and sequenceB.  The calling function places the result in the dispay appropriately.
        /// 
        public ResultTable.Result Align_And_Extract(GeneSequence sequenceA, GeneSequence sequenceB, bool banded)
        {
            ResultTable.Result result = new ResultTable.Result();
            int score;                                                       // place your computed alignment score here
            string[] alignment = new string[2];                              // place your two computed alignments here

            // ***************************************************************************************
            score = 0;                                                
            alignment[0] = "";
            alignment[1] = "";
            // ***************************************************************************************

            if (!banded)
            {
                unrestrictedAlignmentAlgorithm(ref score, ref alignment, ref sequenceA, ref sequenceB);
            }
            else
            {
                bandedAlignmentAlgorithm(ref score, ref alignment, ref sequenceA, ref sequenceB);
            }

            result.Update(score,alignment[0],alignment[1]);                  // bundling your results into the right object type 
            return(result);
        }

        /// <summary>
        /// This method fills the rows with the start values associated with starting the Needleman/Wunsch algorithm.
        /// Time Complexity is O(n+m), it iterates through both sequences, a and b, and which is associated with n and m.
        /// Space Complexity is O(1), because it only modifies already allocated space.
        /// </summary>
        /// <param name="values">The values of each cell</param>
        /// <param name="previous">A list of previous directions that for each cell</param>
        /// <param name="lengthA">The length of gene sequence A</param>
        /// <param name="lengthB">The length of gene sequence B</param>
        /// <param name="banded">Tells if the banded box is checked</param>
        private void fillStartCells(ref int[,] values, ref directions[,] previous, int lengthA, int lengthB, bool banded)
        {
            for(int col = 0; col < lengthB+1; col++)
            {
                if (banded && (col > distance)) break;

                values[0, col] = col * 5;
                previous[0, col] = directions.LEFT;
            }

            for(int row = 0; row < lengthA+1; row++)
            {
                if (banded && (row > distance)) break;

                values[row, 0] = row * 5;
                previous[row, 0] = directions.TOP;
            }
        }

        /// <summary>
        /// This method itertorates through the two gene sequences and finds the correct alignments
        /// The Time Complexity of this method is O(n) because it will iterate through the biggest gene sequence of size n
        /// The Space Complexity of this method is O(n) because a string of the same length of the longest gene sequence is
        ///     allocated
        /// </summary>
        /// <param name="alignment">An array of strings that is the optimal alignments of the two sequences</param>
        /// <param name="previous">An array that holds the previous square/acts as the guiding arrow</param>
        /// <param name="geneSequenceA">The first DNA sequence</param>
        /// <param name="geneSequenceB"> The second DNA sequence</param>
        /// <param name="lengthOfA">the length of geneSequenceA</param>
        /// <param name="lengthOfB">the length of geneSequenceB</param>
        private void createAlignments(ref string[] alignment, ref directions[,] previous, ref GeneSequence geneSequenceA, ref GeneSequence geneSequenceB, ref int lengthOfA, ref int lengthOfB)
        {
            int rowIterator = lengthOfA;
            int colIterator = lengthOfB;
            StringBuilder stringA = new StringBuilder();
            StringBuilder stringB = new StringBuilder();

            while(rowIterator != 0 || colIterator != 0)
            {
                // match or subtitute
                if (previous[rowIterator, colIterator] == directions.DIAGONAL)
                {
                    stringA.Insert(0, geneSequenceA.Sequence[rowIterator - 1]);
                    stringB.Insert(0, geneSequenceB.Sequence[colIterator - 1]);
                    rowIterator--;
                    colIterator--;
                }// insert
                else if (previous[rowIterator, colIterator] == directions.LEFT)
                {
                    stringA.Insert(0, '-');
                    stringB.Insert(0, geneSequenceB.Sequence[colIterator - 1]);
                    colIterator--;
                }// delete
                else
                {
                    stringA.Insert(0, geneSequenceA.Sequence[rowIterator - 1]);
                    stringB.Insert(0, '-');
                    rowIterator--;
                }
            }

            alignment[0] = stringA.ToString().Substring(0, Math.Min(lengthOfA, 100));
            alignment[1] = stringB.ToString().Substring(0, Math.Min(lengthOfB, 100));
        }

        // ********************************************************************************************
        // ************************** Unrestricted Alignment Algorithm ********************************
        // ********************************************************************************************

        /// <summary>
        /// This is the unrestricted algorithm. It tries to find the optimal alignment in the given DNA sequences
        /// The Time Complexity of this method is O(nm), n and m are the lengths of the two gene sequences.  It is 
        ///     nm because it iterates through both sequences.
        /// The Space Complexity of this method is O(nm), n and m are the lengths of the two gene sequences.  It is
        ///     nm because it makes a two demensional array of size n by m.
        /// </summary>
        /// <param name="score">The score of the cell for the given sequence</param>
        /// <param name="alignment">A string array of the optimal alignment of the two sequences</param>
        /// <param name="geneSequenceA">The first given DNA sequence</param>
        /// <param name="geneSequenceB">The second DNA sequence given</param>

        private void unrestrictedAlignmentAlgorithm(ref int score, ref string[] alignment, ref GeneSequence geneSequenceA, ref GeneSequence geneSequenceB)
        {
            int lengthOfA = Math.Min(geneSequenceA.Sequence.Length, MaxCharactersToAlign);
            int lengthOfB = Math.Min(geneSequenceB.Sequence.Length, MaxCharactersToAlign);

            int[,] values = new int[lengthOfA + 1, lengthOfB + 1];
            directions[,] previous = new directions[lengthOfA + 1, lengthOfB + 1];
            fillStartCells(ref values, ref previous, lengthOfA, lengthOfB, false);

            for(int row = 1; row < lengthOfA+1; row++)
            {
                for(int col = 1; col < lengthOfB+1; col++)
                {
                    int costOfTopDelete = values[row - 1, col] + 5;
                    int costOfLeftInsert = values[row, col - 1] + 5;
                    int costOfDiagonalMove = 1;
                    if (geneSequenceA.Sequence[row - 1] == geneSequenceB.Sequence[col - 1]) costOfDiagonalMove = -3;
                    int costOfDiagonal = values[row - 1, col - 1] + costOfDiagonalMove;
                    int minCost = Math.Min(costOfTopDelete, Math.Min(costOfLeftInsert, costOfDiagonal));

                    values[row, col] = minCost;

                    if(minCost == costOfDiagonal)
                    {
                        previous[row, col] = directions.DIAGONAL;
                    }
                    else if (minCost == costOfLeftInsert)
                    {
                        previous[row, col] = directions.LEFT;
                    }
                    else
                    {
                        previous[row, col] = directions.TOP;
                    }
                }
            }

            score = values[lengthOfA, lengthOfB];
            createAlignments(ref alignment, ref previous, ref geneSequenceA, ref geneSequenceB, ref lengthOfA, ref lengthOfB);
        }

        // ********************************************************************************************
        // ****************************** Banded Alignment Algorithm **********************************
        // ********************************************************************************************

        /// <summary>
        /// This is the banded alignment algorithm.
        /// The time complexity of the algorithm is O(n+m), where n and m are the sizes of the two gene sequences.  It is n + m
        ///     because it only iterates over a certain amount of cells.
        /// The Space Complexity is O(nm), n and m are the lengths of the two gene sequences.  It is
        ///     nm because it makes a two demensional array of size n by m.
        /// </summary>
        /// <param name="score">The score that is to be passed back to fill the cell</param>
        /// <param name="alignment">A string array that holds the two aligned strings or states there is no alignment</param>
        /// <param name="geneSequenceA">The first DNA sequence that is given</param>
        /// <param name="geneSequenceB"> The second DNA sequence that is given</param>

        private void bandedAlignmentAlgorithm(ref int score, ref string[] alignment, ref GeneSequence geneSequenceA, ref GeneSequence geneSequenceB)
        {
            int lengthOfA = Math.Min(geneSequenceA.Sequence.Length, MaxCharactersToAlign);
            int lengthOfB = Math.Min(geneSequenceB.Sequence.Length, MaxCharactersToAlign);

            int[,] values = new int[lengthOfA + 1, lengthOfB + 1];
            directions[,] previous = new directions[lengthOfA + 1, lengthOfB + 1];
            fillStartCells(ref values, ref previous, lengthOfA, lengthOfB, true);

            int colStart = 1;
            bool isAligned = false;
            int row = 1;
            int col = colStart;

            for(row = 1; row < lengthOfA+1; row++)
            {
                for(col = colStart; col < lengthOfB+1; col++)
                {
                    if ((distance + row) < col) break;

                    int costOfTopDelete = values[row - 1, col] + 5;
                    if ((distance + row) == col) costOfTopDelete = int.MaxValue;

                    int costOfLeftInsert = values[row, col - 1] + 5;
                    if ((distance + col) == row) costOfLeftInsert = int.MaxValue;

                    int costOfDiagonalMove = 1;
                    if (geneSequenceA.Sequence[row - 1] == geneSequenceB.Sequence[col - 1]) costOfDiagonalMove = -3;
                    int costOfDiagonal = values[row - 1, col - 1] + costOfDiagonalMove;

                    int minCost = Math.Min(costOfTopDelete, Math.Min(costOfLeftInsert, costOfDiagonal));
                    values[row, col] = minCost;

                    if(minCost == costOfDiagonal)
                    {
                        previous[row, col] = directions.DIAGONAL;
                    }
                    else if(minCost == costOfLeftInsert)
                    {
                        previous[row, col] = directions.LEFT;
                    }
                    else
                    {
                        previous[row, col] = directions.TOP;
                    }

                    if (col == lengthOfB && row == lengthOfA) isAligned = true;
                }

                if (row > distance) colStart++;
            }

            if(isAligned)
            {
                score = values[lengthOfA, lengthOfB];
                createAlignments(ref alignment, ref previous, ref geneSequenceA, ref geneSequenceB, ref lengthOfA, ref lengthOfB);
            }
            else
            {
                score = int.MaxValue;
                alignment[0] = "No Alignment Possible";
                alignment[1] = "No Alignment Possible";
            }
        }
    }
}
