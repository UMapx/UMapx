using System;
using System.Runtime.Serialization;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the discrete wavelet.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// <see href="https://en.wikipedia.org/wiki/Wavelet"/>.
    /// </remarks>
    [Serializable]
    public partial class WaveletPacket : ICloneable, ISerializable
    {
        #region Private data
        private float[] lp;        // Low-Pass filter
        private float[] hp;        // High-Pass filer
        private float[] ilp;       // Inverse Low-Pass filter
        private float[] ihp;       // Inverse High-Pass filter
        #endregion

        #region Wavelet components
        /// <summary>
        /// Initializes the discrete wavelet.
        /// </summary>
        /// <param name="lp">Scaling function of forward transform.</param>
        /// <param name="hp">Wavelet function of forward transform.</param>
        /// <param name="ilp">Scaling function of backward transform.</param>
        /// <param name="ihp">Wavelet function of backward transform.</param>
        public WaveletPacket(float[] lp, float[] hp, float[] ilp, float[] ihp)
        {
            this.lp = lp; this.hp = hp; this.ilp = ilp; this.ihp = ihp;
        }
        /// <summary>
        /// Gets or sets the scaling function of forward transform.
        /// </summary>
        public float[] LowPass
        {
            get
            {
                return this.lp;
            }
            set
            {
                this.lp = value;
            }
        }
        /// <summary>
        /// Gets or sets the wavelet function of forward transform.
        /// </summary>
        public float[] HighPass
        {
            get
            {
                return this.hp;
            }
            set
            {
                this.hp = value;
            }
        }
        /// <summary>
        /// Gets or sets the scaling function of backward transform.
        /// </summary>
        public float[] ILowPass
        {
            get
            {
                return ilp;
            }
            set
            {
                this.ilp = value;
            }
        }
        /// <summary>
        /// Gets or sets the wavelet function of backward transform.
        /// </summary>
        public float[] IHighPass
        {
            get
            {
                return ihp;
            }
            set
            {
                this.ihp = value;
            }
        }
        #endregion

        #region Public static voids
        /// <summary>
        /// Builds the analysis high-pass filter h1 from a real-valued analysis low-pass h0
        /// using the Conjugate-Quadrature-Filter (CQF) relation:
        ///     h1[n] = (-1)^n * h0[N-1-n].
        /// For a paraunitary (orthonormal) 2-channel bank this must be combined with a
        /// properly normalized low-pass (e.g., sum(h0)=√2 and H0(π)=0).
        /// </summary>
        /// <param name="v">Scaling function.</param>
        /// <returns>Wavelet function.</returns>
        public static float[] CQF(float[] v)
        {
            // High-pass by CQF:
            // h1[n] = (-1)^n * h0[(N-1-n) mod N]
            var N = v.Length;
            var h = new float[N];

            for (int i = 0; i < N; i++)
            {
                // reverse index (mod N)
                int r = N - 1 - i;
                float sign = ((i & 1) == 0) ? +1f : -1f;
                h[i] = sign * v[r];
            }

            return h;
        }
        /// <summary>
        /// Creates the discrete wavelet.
        /// </summary>
        /// <param name="scaling">Scaling function.</param>
        /// <returns>Discrete wavelet.</returns>
        public static WaveletPacket Create(float[] scaling)
        {
            float[] lp = scaling;
            float[] hp = CQF(lp);
            float[] ilp = Matrice.Flip(lp);
            float[] ihp = Matrice.Flip(hp);

            return new WaveletPacket(lp, hp, ilp, ihp);
        }
        /// <summary>
        /// Creates an orthogonal bank by reversing its analysis filters for synthesis.
        /// </summary>
        /// <param name="scaling">Scaling function.</param>
        /// <param name="wavelet">Wavelet function.</param>
        /// <returns>Discrete wavelet.</returns>
        /// <remarks>For a general biorthogonal bank, supply all four filters to the constructor.</remarks>
        public static WaveletPacket Create(float[] scaling, float[] wavelet)
        {
            float[] lp = scaling;
            float[] hp = wavelet;
            float[] ilp = Matrice.Flip(lp);
            float[] ihp = Matrice.Flip(hp);

            return new WaveletPacket(lp, hp, ilp, ihp);
        }

        /// <summary>
        /// Builds the dual synthesis filters of a real biorthogonal analysis pair.
        /// </summary>
        /// <param name="scaling">Analysis low-pass coefficients with nonzero DC gain.</param>
        /// <param name="wavelet">Analysis high-pass coefficients of a perfect-reconstruction pair.</param>
        /// <returns>A bank with equal, even filter lengths and unit reconstruction gain.</returns>
        /// <remarks>
        /// Padding preserves the analysis origin Length / 2 - 1. Alternating the
        /// signs of the opposite analysis filter cancels aliasing; reversing each
        /// filter independently is valid only for orthogonal banks.
        /// </remarks>
        private static WaveletPacket CreateBiorthogonal(float[] scaling, float[] wavelet)
        {
            int length = Math.Max(scaling.Length, wavelet.Length);
            length += length & 1;
            float[] low = new float[length], high = new float[length];
            Array.Copy(scaling, 0, low, length / 2 - scaling.Length / 2, scaling.Length);
            Array.Copy(wavelet, 0, high, length / 2 - wavelet.Length / 2, wavelet.Length);
            float[] inverseLow = new float[length], inverseHigh = new float[length];
            double lowSum = 0, inverseSum = 0;
            for (int i = 0; i < length; i++)
            {
                float sign = (i & 1) == 0 ? 1 : -1;
                inverseLow[i] = -sign * high[i];
                inverseHigh[i] = sign * low[i];
                lowSum += low[i];
                inverseSum += inverseLow[i];
            }
            // Downsampling halves the DC gain. The sign also depends on the
            // high-pass phase convention used by each coefficient table.
            double gain = lowSum * inverseSum / 2;
            for (int i = 0; i < length; i++)
            {
                inverseLow[i] = (float)(inverseLow[i] / gain);
                inverseHigh[i] = (float)(inverseHigh[i] / gain);
            }
            return new WaveletPacket(low, high, inverseLow, inverseHigh);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the discrete wavelet.
        /// </summary>
        /// <returns>Discrete wavelet.</returns>
        object ICloneable.Clone()
        {
            return new WaveletPacket(
                (float[])this.lp.Clone(),
                (float[])this.hp.Clone(),
                (float[])this.ilp.Clone(),
                (float[])this.ihp.Clone());
        }
        /// <summary>
        /// Creates a copy of the discrete wavelet.
        /// </summary>
        /// <returns>Discrete wavelet.</returns>
        public WaveletPacket Clone()
        {
            return new WaveletPacket(
                (float[])this.lp.Clone(),
                (float[])this.hp.Clone(),
                (float[])this.ilp.Clone(),
                (float[])this.ihp.Clone());
        }
        #endregion

        #region Serialization members
        /// <summary>
        /// Gets information about the object.
        /// </summary>
        /// <param name="info">Data needed for serialization and deserialization.</param>
        /// <param name="context">Source and destination of a given stream.</param>
        public void GetObjectData(SerializationInfo info, StreamingContext context)
        {
            info.AddValue("Scaling function of forward transform", this.lp);
            info.AddValue("Wavelet function of forward transform", this.hp);
            info.AddValue("Scaling function of backward transform", this.ilp);
            info.AddValue("Wavelet function of backward transform", this.ihp);
        }
        #endregion
    }
}
