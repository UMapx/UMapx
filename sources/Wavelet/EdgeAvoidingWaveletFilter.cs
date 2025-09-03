using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the edge-avoiding wavelet filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://www.cs.huji.ac.il/w~raananf/projects/eaw/
    /// </remarks>
    /// </summary>
    [Serializable]
    public class EdgeAvoidingWaveletFilter : IFilter
    {
        #region Private data
        private EdgeAvoidingWaveletDecomposition waveletDecomposition;
        private float factor;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the edge-avoiding wavelet filter.
        /// </summary>
        /// <param name="waveletDecomposition">Edge-avoiding wavelet decomposition</param>
        /// <param name="factor">Factor [-1, 1]</param>
        public EdgeAvoidingWaveletFilter(EdgeAvoidingWaveletDecomposition waveletDecomposition, float factor = -1.0f)
        {
            WaveletDecomposition = waveletDecomposition;
            Factor = factor;
        }
        /// <summary>
        /// Gets or sets edge-avoiding wavelet decomposition.
        /// </summary>
        public EdgeAvoidingWaveletDecomposition WaveletDecomposition
        {
            get
            {
                return waveletDecomposition;
            }
            set
            {
                waveletDecomposition = value;
            }
        }
        /// <summary>
        /// Gets or sets the factor [-1, 1].
        /// </summary>
        public float Factor
        {
            get
            {
                return factor;
            }
            set
            {
                factor = value;
            }
        }
        #endregion

        #region Public apply voids
        /// <summary>
        /// Applies edge-avoiding wavelet filter.
        /// </summary>
        /// <param name="input">Input</param>
        public void Apply(float[] input)
        {
            var dcmp = waveletDecomposition.Forward(input);
            var levels = dcmp.Length;

            for (int i = 1; i < levels; i++)
            {
                dcmp[i] = dcmp[i].Mul(1 + factor);
            }

            var reconstruction = waveletDecomposition.Backward(dcmp);

            for (int i = 0; i < input.Length; i++)
                input[i] = reconstruction[i];
        }
        /// <summary>
        /// Applies edge-avoiding wavelet filter.
        /// </summary>
        /// <param name="input">Input</param>
        public void Apply(float[,] input)
        {
            var dcmp = waveletDecomposition.Forward(input);
            var levels = dcmp.Length;

            for (int i = 1; i < levels; i++)
            {
                dcmp[i] = dcmp[i].Mul(1 + factor);
            }

            var reconstruction = waveletDecomposition.Backward(dcmp);

            for (int i = 0; i < input.GetLength(0); i++)
                for (int j = 0; j < input.GetLength(1); j++)
                    input[i, j] = reconstruction[i, j];
        }
        /// <summary>
        /// Applies edge-avoiding wavelet filter.
        /// </summary>
        /// <param name="input">Input</param>
        public void Apply(ComplexF[] input)
        {
            var dcmp = waveletDecomposition.Forward(input);
            var levels = dcmp.Length;

            for (int i = 1; i < levels; i++)
            {
                dcmp[i] = dcmp[i].Mul(1 + factor);
            }

            var reconstruction = waveletDecomposition.Backward(dcmp);

            for (int i = 0; i < input.Length; i++)
                input[i] = reconstruction[i];
        }
        /// <summary>
        /// Applies edge-avoiding wavelet filter.
        /// </summary>
        /// <param name="input">Input</param>
        public void Apply(ComplexF[,] input)
        {
            var dcmp = waveletDecomposition.Forward(input);
            var levels = dcmp.Length;

            for (int i = 1; i < levels; i++)
            {
                dcmp[i] = dcmp[i].Mul(1 + factor);
            }

            var reconstruction = waveletDecomposition.Backward(dcmp);

            for (int i = 0; i < input.GetLength(0); i++)
                for (int j = 0; j < input.GetLength(1); j++)
                    input[i, j] = reconstruction[i, j];
        }
        #endregion
    }
}