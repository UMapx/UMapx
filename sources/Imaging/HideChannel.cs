using System;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the channel hide filter.
    /// </summary>
    [Serializable]
    public class HideChannel : IBitmapFilter
    {
        #region Private data
        private RGBA channel = RGBA.Red;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the channel hide filter.
        /// </summary>
        /// <param name="channel">Channel of RGBA model</param>
        public HideChannel(RGBA channel)
        {
            Channel = channel;
        }
        /// <summary>
        /// Initializes the channel hide filter.
        /// </summary>
        public HideChannel() { }
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
                    p[c1] = 0;
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
