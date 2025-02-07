using System;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the channel show filter.
    /// </summary>
    [Serializable]
    public class ShowChannel : IBitmapFilter
    {
        #region Private data
        private RGBA channel = RGBA.Red;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the channel show filter.
        /// </summary>
        /// <param name="channel">Channel of RGBA model</param>
        public ShowChannel(RGBA channel)
        {
            Channel = channel;
        }
        /// <summary>
        /// Initializes the channel show filter.
        /// </summary>
        public ShowChannel() { }
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
            for (int y = 0; y < bmData.Height; y++)
            {
                for (int x = 0; x < bmData.Width; x++, p += 4)
                {
                    switch (channel)
                    {
                        case RGBA.Blue:
                            p[1] = 0; p[2] = 0;
                            break;

                        case RGBA.Green:
                            p[2] = 0; p[0] = 0;
                            break;

                        case RGBA.Red:
                            p[1] = 0; p[0] = 0;
                            break;

                        case RGBA.Alpha:
                            p[0] = p[1] = p[2] = 255;
                            break;
                    }
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
