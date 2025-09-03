using System;
using System.Runtime.Serialization;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the discrete wavelet.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Wavelet
    /// </remarks>
    /// </summary>
    [Serializable]
    public partial class WaveletPack : ICloneable, ISerializable
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
        /// <param name="lp">Scaling function of forward transform</param>
        /// <param name="hp">Wavelet function of forward transform</param>
        /// <param name="ilp">Scaling function of backward transform</param>
        /// <param name="ihp">Wavelet function of backward transform</param>
        public WaveletPack(float[] lp, float[] hp, float[] ilp, float[] ihp)
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
        /// <param name="v">Scaling function</param>
        /// <returns>Wavelet function</returns>
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
        /// <param name="scaling">Scaling function</param>
        /// <returns>Discrete wavelet</returns>
        public static WaveletPack Create(float[] scaling)
        {
            float[] lp = scaling;
            float[] hp = CQF(lp);
            float[] ilp = MatrixF.Flip(lp);
            float[] ihp = MatrixF.Flip(hp);

            return new WaveletPack(lp, hp, ilp, ihp);
        }
        /// <summary>
        /// Creates the discrete wavelet.
        /// </summary>
        /// <param name="scaling">Scaling function</param>
        /// <param name="wavelet">Wavelet function</param>
        /// <returns>Discrete wavelet</returns>
        public static WaveletPack Create(float[] scaling, float[] wavelet)
        {
            float[] lp = scaling;
            float[] hp = wavelet;
            float[] ilp = MatrixF.Flip(lp);
            float[] ihp = MatrixF.Flip(hp);

            return new WaveletPack(lp, hp, ilp, ihp);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Creates a copy of the discrete wavelet.
        /// </summary>
        /// <returns>Discrete wavelet</returns>
        object ICloneable.Clone()
        {
            return new WaveletPack(
                (float[])this.lp.Clone(),
                (float[])this.hp.Clone(),
                (float[])this.ilp.Clone(),
                (float[])this.ihp.Clone());
        }
        /// <summary>
        /// Creates a copy of the discrete wavelet.
        /// </summary>
        /// <returns>Discrete wavelet</returns>
        public WaveletPack Clone()
        {
            return new WaveletPack(
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
        /// <param name="info">Data needed for serialization and deserialization</param>
        /// <param name="context">Source and destination of a given stream</param>
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
