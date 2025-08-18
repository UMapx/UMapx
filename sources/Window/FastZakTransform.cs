using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the fast Zak transform.
    /// </summary>
    [Serializable]
    public class FastZakTransform : ZakTransform, IZakTransform
    {
        #region Initialize
        /// <summary>
        /// Initializes the fast Zak transform.
        /// </summary>
        /// <param name="m">Number of frequency shifts [4, N/2]</param>
        public FastZakTransform(int m) : base(m)
        {
            this.DFT = new FastFourierTransform(false, Direction.Vertical);
        }
        #endregion
    }
}
