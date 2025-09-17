using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Hartley transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_Hartley_transform
    /// </remarks>
    [Serializable]
    public class HartleyTransform : TransformBaseMatrixFloat, ITransform
    {
        #region Initialize
        /// <summary>
        /// Initializes the Hartley transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="spectrumType">Spectrum type</param>
        /// <param name="direction">Processing direction</param>
        public HartleyTransform(bool normalized = true, SpectrumType spectrumType = SpectrumType.Fourier, Direction direction = Direction.Vertical)
        {
            this.Normalized = normalized; 
            this.Direction = direction;
            this.SpectrumType = spectrumType;
        }
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        public bool Normalized { get; set; }
        /// <summary>
        /// Gets or sets spectrum type.
        /// </summary>
        public SpectrumType SpectrumType { get; set; }
        #endregion

        #region Hartley static components
        /// <summary>
        /// Implements the construction of the Hartley transform matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="spectrumType">Spectrum type</param>
        /// <returns>Matrix</returns>
        public static float[,] Matrix(int n, bool normalized = false, SpectrumType spectrumType = SpectrumType.Hartley)
        {
            if (spectrumType == SpectrumType.Fourier)
            {
                var F = FourierTransform.Matrix(n);
                return F.ToReal().Sub(F.ToImag());
            }

            int j, i;
            float[,] H = new float[n, n];
            float s = 2.0f * Maths.Pi / n;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = Special.Cas(s * i * j);
                }
            }

            if (normalized)
            {
                H = H.Div(Maths.Sqrt(n));
            }

            return H;
        }
        #endregion

        #region Hartley Transform
        /// <inheritdoc/>
        public override float[,] Matrix(int n, bool backward)
        {
            var U = Matrix(n, Normalized, SpectrumType);

            if (backward)
            {
                U = U.Transpose();
            }

            return U;
        }
        #endregion
    }
}
