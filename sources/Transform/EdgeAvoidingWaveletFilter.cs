using System;
using UMapx.Core;

namespace UMapx.Transform
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
        private float sigmaSpatial;
        private float sigmaRange;
        private float factor;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the edge-avoiding wavelet filter.
        /// </summary>
        /// <param name="sigmaSpatial">Spatial smoothing factor</param>
        /// <param name="sigmaRange">Range smoothing factor</param>
        /// <param name="factor">Factor [-1, 1]</param>
        public EdgeAvoidingWaveletFilter(float sigmaSpatial = 4f, float sigmaRange = 0.1f, float factor = -1.0f)
        {
            this.SigmaSpatial = sigmaSpatial;
            this.SigmaRange = sigmaRange;
            this.Factor = factor;
        }

        /// <summary>
        /// Gets or sets the spatial smoothing factor.
        /// </summary>
        public float SigmaSpatial
        {
            get
            {
                return this.sigmaSpatial;
            }
            set
            {
                this.sigmaSpatial = value;
            }
        }

        /// <summary>
        /// Gets or sets the range smoothing factor.
        /// </summary>
        public float SigmaRange
        {
            get
            {
                return this.sigmaRange;
            }
            set
            {
                this.sigmaRange = value;
            }
        }
        /// <summary>
        /// Gets or sets the factor [-1, 1].
        /// </summary>
        public float Factor
        {
            get
            {
                return this.factor;
            }
            set
            {
                this.factor = value;
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
            var (baseLevel, detail2, detail1) = Decompose(input, sigmaSpatial, sigmaRange);

            detail1 = Matrice.Mul(detail1, 1 + factor);
            detail2 = Matrice.Mul(detail2, 1 + factor);

            var reconstruction = Reconstruct(baseLevel, detail2, detail1);

            for (int i = 0; i < input.Length; i++)
                input[i] = reconstruction[i];
        }
        /// <summary>
        /// Applies edge-avoiding wavelet filter.
        /// </summary>
        /// <param name="input">Input</param>
        public void Apply(float[,] input)
        {
            var (baseLevel, detail2, detail1) = Decompose(input, sigmaSpatial, sigmaRange);

            detail1 = Matrice.Mul(detail1, 1 + factor);
            detail2 = Matrice.Mul(detail2, 1 + factor);

            var reconstruction = Reconstruct(baseLevel, detail2, detail1);

            for (int i = 0; i < input.GetLength(0); i++)
                for (int j = 0; j < input.GetLength(1); j++)
                    input[i, j] = reconstruction[i, j];
        }
        /// <summary>
        /// Applies edge-avoiding wavelet filter.
        /// </summary>
        /// <param name="input">Input</param>
        public void Apply(Complex32[] input)
        {
            var (baseLevel, detail2, detail1) = Decompose(input, sigmaSpatial, sigmaRange);

            detail1 = Matrice.Mul(detail1, 1 + factor);
            detail2 = Matrice.Mul(detail2, 1 + factor);

            var reconstruction = Reconstruct(baseLevel, detail2, detail1);

            for (int i = 0; i < input.Length; i++)
                input[i] = reconstruction[i];
        }
        /// <summary>
        /// Applies edge-avoiding wavelet filter.
        /// </summary>
        /// <param name="input">Input</param>
        public void Apply(Complex32[,] input)
        {
            var (baseLevel, detail2, detail1) = Decompose(input, sigmaSpatial, sigmaRange);

            detail1 = Matrice.Mul(detail1, 1 + factor);
            detail2 = Matrice.Mul(detail2, 1 + factor);

            var reconstruction = Reconstruct(baseLevel, detail2, detail1);

            for (int i = 0; i < input.GetLength(0); i++)
                for (int j = 0; j < input.GetLength(1); j++)
                    input[i, j] = reconstruction[i, j];
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Decomposes array into base and detail layers.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="sigmaSpatial">Spatial smoothing factor</param>
        /// <param name="sigmaRange">Range smoothing factor</param>
        /// <returns>Array</returns>
        internal static (float[], float[], float[]) Decompose(float[] input, float sigmaSpatial = 4f, float sigmaRange = 0.1f)
        {
            var base1 = (float[])input.Clone();
            BilateralGridFilter.bilateralgridfilter(base1, sigmaSpatial, sigmaRange);
            var detail1 = Matrice.Sub(input, base1);

            var base2 = (float[])base1.Clone();
            BilateralGridFilter.bilateralgridfilter(base2, sigmaSpatial, sigmaRange);
            var detail2 = Matrice.Sub(base1, base2);

            return (base2, detail2, detail1);
        }
        /// <summary>
        /// Reconstructs array from base and detail layers.
        /// </summary>
        /// <param name="baseLevel">Base level</param>
        /// <param name="detail2">Detail 1</param>
        /// <param name="detail1">Detail 2</param>
        /// <returns>Array</returns>
        internal static float[] Reconstruct(float[] baseLevel, float[] detail2, float[] detail1)
        {
            var level1 = Matrice.Add(baseLevel, detail2);
            return Matrice.Add(level1, detail1);
        }
        /// <summary>
        /// Decomposes array into base and detail layers.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="sigmaSpatial">Spatial smoothing factor</param>
        /// <param name="sigmaRange">Range smoothing factor</param>
        /// <returns>Array</returns>
        internal static (Complex32[], Complex32[], Complex32[]) Decompose(Complex32[] input, float sigmaSpatial = 4f, float sigmaRange = 0.1f)
        {
            var base1 = (Complex32[])input.Clone();
            BilateralGridFilter.bilateralgridfilter(base1, sigmaSpatial, sigmaRange);
            var detail1 = Matrice.Sub(input, base1);

            var base2 = (Complex32[])base1.Clone();
            BilateralGridFilter.bilateralgridfilter(base2, sigmaSpatial, sigmaRange);
            var detail2 = Matrice.Sub(base1, base2);

            return (base2, detail2, detail1);
        }
        /// <summary>
        /// Reconstructs array from base and detail layers.
        /// </summary>
        /// <param name="baseLevel">Base level</param>
        /// <param name="detail2">Detail 1</param>
        /// <param name="detail1">Detail 2</param>
        /// <returns>Array</returns>
        internal static Complex32[] Reconstruct(Complex32[] baseLevel, Complex32[] detail2, Complex32[] detail1)
        {
            var level1 = Matrice.Add(baseLevel, detail2);
            return Matrice.Add(level1, detail1);
        }

        /// <summary>
        /// Decomposes matrix into base and detail layers.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <param name="sigmaSpatial">Spatial smoothing factor</param>
        /// <param name="sigmaRange">Range smoothing factor</param>
        /// <returns>Matrix</returns>
        internal static (float[,], float[,], float[,]) Decompose(float[,] input, float sigmaSpatial = 4f, float sigmaRange = 0.1f)
        {
            var base1 = (float[,])input.Clone();
            BilateralGridFilter.bilateralgridfilter(base1, sigmaSpatial, sigmaRange);
            var detail1 = Matrice.Sub(input, base1);

            var base2 = (float[,])base1.Clone();
            BilateralGridFilter.bilateralgridfilter(base2, sigmaSpatial, sigmaRange);
            var detail2 = Matrice.Sub(base1, base2);

            return (base2, detail2, detail1);
        }
        /// <summary>
        /// Reconstructs matrix from base and detail layers.
        /// </summary>
        /// <param name="baseLevel">Base level</param>
        /// <param name="detail2">Detail 1</param>
        /// <param name="detail1">Detail 2</param>
        /// <returns>Matrix</returns>
        internal static float[,] Reconstruct(float[,] baseLevel, float[,] detail2, float[,] detail1)
        {
            var level1 = Matrice.Add(baseLevel, detail2);
            return Matrice.Add(level1, detail1);
        }
        /// <summary>
        /// Decomposes matrix into base and detail layers.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <param name="sigmaSpatial">Spatial smoothing factor</param>
        /// <param name="sigmaRange">Range smoothing factor</param>
        /// <returns>Matrix</returns>
        internal static (Complex32[,], Complex32[,], Complex32[,]) Decompose(Complex32[,] input, float sigmaSpatial = 4f, float sigmaRange = 0.1f)
        {
            var base1 = (Complex32[,])input.Clone();
            BilateralGridFilter.bilateralgridfilter(base1, sigmaSpatial, sigmaRange);
            var detail1 = Matrice.Sub(input, base1);

            var base2 = (Complex32[,])base1.Clone();
            BilateralGridFilter.bilateralgridfilter(base2, sigmaSpatial, sigmaRange);
            var detail2 = Matrice.Sub(base1, base2);

            return (base2, detail2, detail1);
        }
        /// <summary>
        /// Reconstructs matrix from base and detail layers.
        /// </summary>
        /// <param name="baseLevel">Base level</param>
        /// <param name="detail2">Detail 1</param>
        /// <param name="detail1">Detail 2</param>
        /// <returns>Matrix</returns>
        internal static Complex32[,] Reconstruct(Complex32[,] baseLevel, Complex32[,] detail2, Complex32[,] detail1)
        {
            var level1 = Matrice.Add(baseLevel, detail2);
            return Matrice.Add(level1, detail1);
        }
        #endregion
    }
}