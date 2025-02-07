using System;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the channel equalization filter.
    /// </summary>
    [Serializable]
    public class EqualizeChannel : IBitmapFilter
    {
        #region Private data
        private RGBA channel = RGBA.Red;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the channel equalization filter.
        /// </summary>
        /// <param name="channel">Channel of RGBA model</param>
        public EqualizeChannel(RGBA channel)
        {
            Channel = channel;
        }
        /// <summary>
        /// Initializes the channel equalization filter.
        /// </summary>
        public EqualizeChannel() { }
        /// <summary>
        /// Gets or sets the channel of the RGBA model.
        /// </summary>
        public RGBA Channel
        {
            get
            {
                return this.channel;
            }
            set
            {
                this.channel = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int c1 = (int)this.channel;

            for (int y = 0; y < bmData.Height; y++)
            {
                for (int x = 0; x < bmData.Width; x++, p += 4)
                {
                    p[2] = p[c1];
                    p[1] = p[c1];
                    p[0] = p[c1];
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            Apply(bmData);
            BitmapFormat.Unlock(Data, bmData);
        }
        #endregion
    }
}
