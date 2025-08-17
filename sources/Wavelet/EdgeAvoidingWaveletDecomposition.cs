using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the edge-avoiding wavelet decomposition.
    /// <remarks>
    /// More information can be found on the website:
    /// https://www.cs.huji.ac.il/w~raananf/projects/eaw/
    /// </remarks>
    /// </summary>
    [Serializable]
    public class EdgeAvoidingWaveletDecomposition : IPyramidTransform
    {
        #region Private data
        private float sigmaSpatial = 4f;
        private float sigmaRange = 0.1f;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the edge-avoiding wavelet decomposition.
        /// </summary>
        /// <param name="sigmaSpatial">Spatial smoothing factor</param>
        /// <param name="sigmaRange">Range smoothing factor</param>
        public EdgeAvoidingWaveletDecomposition(float sigmaSpatial = 4f, float sigmaRange = 0.1f)
        {
            SigmaSpatial = sigmaSpatial;
            SigmaRange = sigmaRange;
        }
        /// <summary>
        /// Gets or sets spatial smoothing factor.
        /// </summary>
        public float SigmaSpatial
        {
            get
            {
                return sigmaSpatial;
            }
            set
            {
                sigmaSpatial = value;
            }
        }
        /// <summary>
        /// Gets or sets range smoothing factor.
        /// </summary>
        public float SigmaRange
        {
            get
            {
                return sigmaRange;
            }
            set 
            { 
                sigmaRange = value; 
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Decomposes array into base and detail layers.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Array</returns>
        public float[][] Forward(float[] input)
        {
            var base1 = (float[])input.Clone();
            BilateralGridFilter.bilateralgridfilter(base1, sigmaSpatial, sigmaRange);
            var detail1 = input.Sub(base1);

            var base2 = (float[])base1.Clone();
            BilateralGridFilter.bilateralgridfilter(base2, sigmaSpatial, sigmaRange);
            var detail2 = base1.Sub(base2);

            return new float[][] { base2, detail2, detail1 };
        }
        /// <summary>
        /// Reconstructs array from base and detail layers.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[][] input)
        {
            return input[0].Add(input[1]).Add(input[2]);
        }
        /// <summary>
        /// Decomposes matrix into base and detail layers.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        public float[][,] Forward(float[,] input)
        {
            var base1 = (float[,])input.Clone();
            BilateralGridFilter.bilateralgridfilter(base1, sigmaSpatial, sigmaRange);
            var detail1 = input.Sub(base1);

            var base2 = (float[,])base1.Clone();
            BilateralGridFilter.bilateralgridfilter(base2, sigmaSpatial, sigmaRange);
            var detail2 = base1.Sub(base2);

            return new float[][,] { base2, detail2, detail1 };
        }
        /// <summary>
        /// Reconstructs matrix from base and detail layers.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[][,] input)
        {
            return input[0].Add(input[1]).Add(input[2]);
        }
        /// <summary>
        /// Decomposes array into base and detail layers.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Array</returns>
        public Complex32[][] Forward(Complex32[] input)
        {
            var base1 = (Complex32[])input.Clone();
            BilateralGridFilter.bilateralgridfilter(base1, sigmaSpatial, sigmaRange);
            var detail1 = input.Sub(base1);

            var base2 = (Complex32[])base1.Clone();
            BilateralGridFilter.bilateralgridfilter(base2, sigmaSpatial, sigmaRange);
            var detail2 = base1.Sub(base2);

            return new Complex32[][] { base2, detail2, detail1 };
        }
        /// <summary>
        /// Reconstructs array from base and detail layers.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[][] input)
        {
            return input[0].Add(input[1]).Add(input[2]);
        }
        /// <summary>
        /// Decomposes matrix into base and detail layers.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[][,] Forward(Complex32[,] input)
        {
            var base1 = (Complex32[,])input.Clone();
            BilateralGridFilter.bilateralgridfilter(base1, sigmaSpatial, sigmaRange);
            var detail1 = input.Sub(base1);

            var base2 = (Complex32[,])base1.Clone();
            BilateralGridFilter.bilateralgridfilter(base2, sigmaSpatial, sigmaRange);
            var detail2 = base1.Sub(base2);

            return new Complex32[][,] { base2, detail2, detail1 };
        }
        /// <summary>
        /// Reconstructs matrix from base and detail layers.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[][,] input)
        {
            return input[0].Add(input[1]).Add(input[2]);
        }
        #endregion
    }
}
