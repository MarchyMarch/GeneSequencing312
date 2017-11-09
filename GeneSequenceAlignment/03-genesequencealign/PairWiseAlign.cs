using System;
using System.Collections.Generic;
using System.Text;

namespace GeneticsLab
{
    enum directions { LEFT, TOP, RIGHT};
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


            // ********* these are placeholder assignments that you'll replace with your code  *******
            score = 0;                                                
            alignment[0] = "";
            alignment[1] = "";
            // ***************************************************************************************
            

            result.Update(score,alignment[0],alignment[1]);                  // bundling your results into the right object type 
            return(result);
        }

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

        private void createAlignments(ref string[] alignment, ref directions[,] previous, ref GeneSequence geneSequenceA, ref GeneSequence geneSequenceB, ref int lengthOfA, ref int lengthOfB)
        {

        }

        // ********************************************************************************************
        // ************************** Unrestricted Alignment Algorithm ********************************
        // ********************************************************************************************

        private void unrestrictedAlignmentAlgorithm(ref int score, ref string[] alignment, ref GeneSequence geneSequenceA, ref GeneSequence geneSequenceB)
        {

        }

        // ********************************************************************************************
        // ****************************** Banded Alignment Algorithm **********************************
        // ********************************************************************************************

        private void bandedAlignmentAlgorithm(ref int score, ref string[] alignment, ref GeneSequence geneSequenceA, ref GeneSequence geneSequenceB)
        {

        }
    }
}
