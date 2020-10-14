using System;

namespace UMapx.Distribution
{
    /// <summary>
    /// Defines the Bayes probability class.
    /// </summary>
    [Serializable]
    public class Bayes
    {
        #region Private data
        private double[] Pp;
        private double Pa;
        private int N;
        #endregion

        #region Bayes components
        /// <summary>
        /// Initializes the Bayes probability class.
        /// </summary>
        /// <param name="stat">Array of statistical probabilities</param>
        /// <param name="prior">An array of a priori probabilities (before experiment)</param>
        public Bayes(double[] stat, double[] prior)
        {
            if (stat.Length != prior.Length)
                throw new Exception("Arrays must be of the same dimensions.");


            this.N = prior.Length;
            this.Pp = new double[N];
            this.Pa = 0;
            int i;

            for (i = 0; i < N; i++)
            {
                Pa += stat[i] * prior[i];
            }

            for (i = 0; i < N; i++)
            {
                Pp[i] = stat[i] * prior[i] / Pa;
            }
        }
        /// <summary>
        /// Returns the value of the total probability.
        /// </summary>
        public double General
        {
            get
            {
                return Pa;
            }
        }
        /// <summary>
        /// Returns an array of values of posterior probabilities (after the experiment).
        /// </summary>
        public double[] Probabilities
        {
            get
            {
                return Pp;
            }
        }
        #endregion
    }
}
