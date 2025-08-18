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
        private float sigmaSpatial;
        private float sigmaRange;
        private int levels;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the edge-avoiding wavelet decomposition.
        /// </summary>
        /// <param name="sigmaSpatial">Spatial smoothing factor</param>
        /// <param name="sigmaRange">Range smoothing factor</param>
        /// <param name="levels">Levels</param>
        public EdgeAvoidingWaveletDecomposition(float sigmaSpatial = 4f, float sigmaRange = 0.1f, int levels = 3)
        {
            SigmaSpatial = sigmaSpatial;
            SigmaRange = sigmaRange;
            Levels = levels;
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
        /// <summary>
        /// Gets or sets the number of pyramid levels (>=1), including the base level.
        /// </summary>
        public int Levels
        {
            get
            {
                return levels;
            }
            set
            {
                levels = (value < 1) ? 1 : value;
            }
        }
        #endregion

        #region Public voids
        /// <summary>
        /// Decomposes array into base and detail layers.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Array</returns>
        public float[][] Forward(float[] input)
        {
            var current = (float[])input.Clone();
            var details = new float[levels - 1][];

            for (int i = 0; i < levels - 1; i++)
            {
                var baseI = (float[])current.Clone();
                BilateralGridFilter.bilateralgridfilter(baseI, sigmaSpatial, sigmaRange);
                var detailI = current.Sub(baseI);
                details[i] = detailI;
                current = baseI;
            }

            var output = new float[levels][];
            output[0] = current;

            for (int i = 0; i < details.Length; i++)
                output[i + 1] = details[details.Length - 1 - i];

            return output;
        }
        /// <summary>
        /// Reconstructs array from base and detail layers.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[][] input)
        {
            var result = (float[])input[0].Clone();
            for (int i = 1; i < input.Length; i++)
                result = result.Add(input[i]);
            return result;
        }
        /// <summary>
        /// Decomposes matrix into base and detail layers.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        public float[][,] Forward(float[,] input)
        {
            var current = (float[,])input.Clone();
            var details = new float[levels - 1][,];

            for (int i = 0; i < levels - 1; i++)
            {
                var baseI = (float[,])current.Clone();
                BilateralGridFilter.bilateralgridfilter(baseI, sigmaSpatial, sigmaRange);
                var detailI = current.Sub(baseI);
                details[i] = detailI;
                current = baseI;
            }

            var output = new float[levels][,];
            output[0] = current;

            for (int i = 0; i < details.Length; i++)
                output[i + 1] = details[details.Length - 1 - i];

            return output;
        }
        /// <summary>
        /// Reconstructs matrix from base and detail layers.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[][,] input)
        {
            var result = (float[,])input[0].Clone();
            for (int i = 1; i < input.Length; i++)
                result = result.Add(input[i]);
            return result;
        }
        /// <summary>
        /// Decomposes array into base and detail layers.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Array</returns>
        public Complex32[][] Forward(Complex32[] input)
        {
            var current = (Complex32[])input.Clone();
            var details = new Complex32[levels - 1][];

            for (int i = 0; i < levels - 1; i++)
            {
                var baseI = (Complex32[])current.Clone();
                BilateralGridFilter.bilateralgridfilter(baseI, sigmaSpatial, sigmaRange);
                var detailI = current.Sub(baseI);
                details[i] = detailI;
                current = baseI;
            }

            var output = new Complex32[levels][];
            output[0] = current;

            for (int i = 0; i < details.Length; i++)
                output[i + 1] = details[details.Length - 1 - i];

            return output;
        }
        /// <summary>
        /// Reconstructs array from base and detail layers.
        /// </summary>
        /// <param name="input">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[][] input)
        {
            var result = (Complex32[])input[0].Clone();
            for (int i = 1; i < input.Length; i++)
                result = result.Add(input[i]);
            return result;
        }
        /// <summary>
        /// Decomposes matrix into base and detail layers.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[][,] Forward(Complex32[,] input)
        {
            var current = (Complex32[,])input.Clone();
            var details = new Complex32[levels - 1][,];

            for (int i = 0; i < levels - 1; i++)
            {
                var baseI = (Complex32[,])current.Clone();
                BilateralGridFilter.bilateralgridfilter(baseI, sigmaSpatial, sigmaRange);
                var detailI = current.Sub(baseI);
                details[i] = detailI;
                current = baseI;
            }

            var output = new Complex32[levels][,];
            output[0] = current;

            for (int i = 0; i < details.Length; i++)
                output[i + 1] = details[details.Length - 1 - i];

            return output;
        }
        /// <summary>
        /// Reconstructs matrix from base and detail layers.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[][,] input)
        {
            var result = (Complex32[,])input[0].Clone();
            for (int i = 1; i < input.Length; i++)
                result = result.Add(input[i]);
            return result;
        }
        #endregion
    }
}
