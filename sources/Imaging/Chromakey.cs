using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the chromakey filter.
    /// </summary>
    [Serializable]
    public class Chromakey : IBitmapFilter
    {
        #region Class components
        /// <summary>
        /// Initializes the chromakey filter.
        /// </summary>
        /// <param name="channel">Channel</param>
        public Chromakey(RGBA channel = RGBA.Green)
        {
            Channel = channel;
        }
        /// <summary>
        /// Gets or sets channel.
        /// </summary>
        public RGBA Channel { get; set; }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // channel subtraction
            if (Channel == RGBA.Red)
            {
                Parallel.For(0, height, y =>
                {
                    int x, ystride, k;
                    ystride = y * stride;

                    for (x = 0; x < width; x++)
                    {
                        k = ystride + x * 4;
                        var luma = RGB.PAL(p[k + 2], p[k + 1], p[k]);
                        var diff = p[k + 2] - luma;

                        p[k] = p[k + 1] = p[k + 2] = Maths.Byte(127.5f + diff);
                    }
                });
            }
            else if (Channel == RGBA.Green)
            {
                Parallel.For(0, height, y =>
                {
                    int x, ystride, k;
                    ystride = y * stride;

                    for (x = 0; x < width; x++)
                    {
                        k = ystride + x * 4;
                        var luma = RGB.PAL(p[k + 2], p[k + 1], p[k]);
                        var diff = p[k + 1] - luma;

                        p[k] = p[k + 1] = p[k + 2] = Maths.Byte(127.5f + diff);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, y =>
                {
                    int x, ystride, k;
                    ystride = y * stride;

                    for (x = 0; x < width; x++)
                    {
                        k = ystride + x * 4;
                        var luma = RGB.PAL(p[k + 2], p[k + 1], p[k]);
                        var diff = p[k] - luma;

                        p[k] = p[k + 1] = p[k + 2] = Maths.Byte(127.5f + diff);
                    }
                });
            }

            // histogram processing
            int[] hist = Statistics.Histogram(bmData);
            int n = hist.Length;
            int z = n / 2;
            int[] lef = new int[z];
            int[] rig = new int[z];

            for (int i = 0; i < z; i++)
            {
                lef[i] = hist[i];
                rig[i] = hist[i + z];
            }

            _ = Statistics.Max(lef, out int index1);
            _ = Statistics.Max(rig, out int index2);
            index2 += z;

            int count = index2 - index1;
            int[] center = new int[count];

            for (int i = 0; i < count; i++)
            {
                center[i] = hist[i + index1];
            }

            _ = Statistics.Min(center, out int threshold);
            threshold += index1;

            // creating mask
            Parallel.For(0, height, y =>
            {
                int x, ystride, k;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;

                    p[k] = p[k + 1] = p[k + 2] = (p[k] >= threshold) ? (byte)0 : (byte)255;
                }
            });

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
