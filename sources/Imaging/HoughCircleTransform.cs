using System;
using System.Collections;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the Hough circle transform filter.
    /// </summary>
    [Serializable]
    public class HoughCircleTransform : IBitmapFilter
    {
        #region Private data
        private int radiusToDetect;

        private short[,] houghMap;
        private short maxMapIntensity = 0;

        private int width;
        private int height;

        private int localPeakRadius = 4;
        private short minCircleIntensity = 10;
        private ArrayList circles = new ArrayList();
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Hough circle transform filter.
        /// </summary>
        /// <param name="radiusToDetect">Radius</param>
        public HoughCircleTransform(int radiusToDetect)
        {
            this.radiusToDetect = radiusToDetect;
        }
        #endregion

        #region Class components
        /// <summary>
        /// Gets the minimum circle intensity.
        /// </summary>
        public short MinCircleIntensity
        {
            get { return minCircleIntensity; }
            set { minCircleIntensity = value; }
        }
        /// <summary>
        /// Gets or sets the radius search for the local peak value.
        /// </summary>
        public int LocalPeakRadius
        {
            get { return localPeakRadius; }
            set { localPeakRadius = Math.Max(1, Math.Min(10, value)); }
        }
        /// <summary>
        /// Gets the maximum line intensity.
        /// </summary>
        public short MaxIntensity
        {
            get { return maxMapIntensity; }
        }
        /// <summary>
        /// Gets the count of lines found.
        /// </summary>
        public int CirclesCount
        {
            get { return circles.Count; }
        }
        /// <summary>
        /// Returns an array of lines with absolute intensity.
        /// </summary>
        /// <param name="count">Count</param>
        /// <returns>Array</returns>
        public HoughCircle[] GetMostIntensiveCircles(int count)
        {
            // lines count
            int n = Math.Min(count, circles.Count);

            // result array
            HoughCircle[] dst = new HoughCircle[n];
            circles.CopyTo(0, dst, 0, n);

            return dst;
        }
        /// <summary>
        /// Returns an array of lines with relative intensity.
        /// </summary>
        /// <param name="minRelativeIntensity">Minimum relative intensity</param>
        /// <returns>Array</returns>
        public HoughCircle[] GetCirclesByRelativeIntensity(double minRelativeIntensity)
        {
            int count = 0, n = circles.Count;

            while ((count < n) && (((HoughCircle)circles[count]).RelativeIntensity >= minRelativeIntensity))
                count++;

            return GetMostIntensiveCircles(count);
        }
        #endregion

        #region Image processing
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            // lock source image
            BitmapData imageData = BitmapConverter.Lock8bpp(Data);

            try
            {
                // process the image
                Apply(imageData);
            }
            finally
            {
                // unlock image
                Data.UnlockBits(imageData);
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public void Apply(BitmapData bmData)
        {
            if (bmData.PixelFormat != PixelFormat.Format8bppIndexed)
            {
                throw new Exception("Image format should be: 8bppIndexed");
            }

            // get source image size
            width = bmData.Width;
            height = bmData.Height;

            int srcOffset = bmData.Stride - width;

            // allocate Hough map of the same size like image
            houghMap = new short[height, width];

            // do the job
            unsafe
            {
                byte* src = (byte*)bmData.Scan0.ToPointer();

                // for each row
                for (int y = 0; y < height; y++)
                {
                    // for each pixel
                    for (int x = 0; x < width; x++, src++)
                    {
                        if (*src != 0)
                        {
                            DrawHoughCircle(x, y);
                        }
                    }
                    src += srcOffset;
                }
            }

            // find max value in Hough map
            maxMapIntensity = 0;
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    if (houghMap[i, j] > maxMapIntensity)
                    {
                        maxMapIntensity = houghMap[i, j];
                    }
                }
            }

            CollectCircles();
        }
        /// <summary>
        /// Returns the Hough synogram.
        /// </summary>
        /// <returns>Bitmap</returns>
        public Bitmap ToBitmap()
        {
            // check if Hough transformation was made already
            if (houghMap == null)
            {
                throw new ApplicationException("Hough conversion was not made");
            }

            int width = houghMap.GetLength(1);
            int height = houghMap.GetLength(0);

            // create new image
            Bitmap image = new Bitmap(width, height, PixelFormat.Format8bppIndexed);

            // get palette
            ColorPalette cp = image.Palette;
            // init palette
            for (int i = 0; i < 256; i++)
            {
                cp.Entries[i] = Color.FromArgb(i, i, i);
            }
            // set palette back
            image.Palette = cp;

            // lock destination bitmap data
            BitmapData imageData = BitmapConverter.Lock8bpp(image);
            float scale = 255.0f / maxMapIntensity;

            // do the job
            unsafe
            {
                byte* dst = (byte*)imageData.Scan0.ToPointer();
                int y, x;

                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++, dst++)
                    {
                        *dst = Maths.Byte(scale * houghMap[y, x]);
                    }
                }
            }

            // unlock destination images
            image.UnlockBits(imageData);

            return image;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Collect circles with intesities greater or equal then specified.
        /// </summary>
        private void CollectCircles()
        {
            short intensity;
            bool foundGreater;

            // clean circles collection
            circles.Clear();
            int y, x, ty, tx, txMax, tyMax;

            // for each Y coordinate
            for (y = 0; y < height; y++)
            {
                // for each X coordinate
                for (x = 0; x < width; x++)
                {
                    // get current value
                    intensity = houghMap[y, x];

                    if (intensity < minCircleIntensity)
                        continue;

                    foundGreater = false;

                    // check neighboors
                    for (ty = y - localPeakRadius, tyMax = y + localPeakRadius; ty < tyMax; ty++)
                    {
                        // continue if the coordinate is out of map
                        if (ty < 0)
                            continue;
                        // break if it is not local maximum or coordinate is out of map
                        if ((foundGreater == true) || (ty >= height))
                            break;

                        for (tx = x - localPeakRadius, txMax = x + localPeakRadius; tx < txMax; tx++)
                        {
                            // continue or break if the coordinate is out of map
                            if (tx < 0)
                                continue;
                            if (tx >= width)
                                break;

                            // compare the neighboor with current value
                            if (houghMap[ty, tx] > intensity)
                            {
                                foundGreater = true;
                                break;
                            }
                        }
                    }

                    // was it local maximum ?
                    if (!foundGreater)
                    {
                        // we have local maximum
                        circles.Add(new HoughCircle(x, y, radiusToDetect, intensity, (double)intensity / maxMapIntensity));
                    }
                }
            }

            circles.Sort();
        }
        /// <summary>
        /// Draw Hough circle:
        /// http://www.cs.unc.edu/~mcmillan/comp136/Lecture7/circle.html
        /// TODO: more optimizations of circle drawing could be done.
        /// </summary>
        /// <param name="xCenter">Co. X</param>
        /// <param name="yCenter">Co. Y</param>
        private void DrawHoughCircle(int xCenter, int yCenter)
        {
            int x = 0;
            int y = radiusToDetect;
            int p = (5 - radiusToDetect * 4) / 4;

            SetHoughCirclePoints(xCenter, yCenter, x, y);

            while (x < y)
            {
                x++;
                if (p < 0)
                {
                    p += 2 * x + 1;
                }
                else
                {
                    y--;
                    p += 2 * (x - y) + 1;
                }
                SetHoughCirclePoints(xCenter, yCenter, x, y);
            }
        }
        /// <summary>
        /// Set circle points.
        /// </summary>
        /// <param name="cx">Cx</param>
        /// <param name="cy">Cy</param>
        /// <param name="x">Co. X</param>
        /// <param name="y">Co. Y</param>
        private void SetHoughCirclePoints(int cx, int cy, int x, int y)
        {
            if (x == 0)
            {
                SetHoughPoint(cx, cy + y);
                SetHoughPoint(cx, cy - y);
                SetHoughPoint(cx + y, cy);
                SetHoughPoint(cx - y, cy);
            }
            else if (x == y)
            {
                SetHoughPoint(cx + x, cy + y);
                SetHoughPoint(cx - x, cy + y);
                SetHoughPoint(cx + x, cy - y);
                SetHoughPoint(cx - x, cy - y);
            }
            else if (x < y)
            {
                SetHoughPoint(cx + x, cy + y);
                SetHoughPoint(cx - x, cy + y);
                SetHoughPoint(cx + x, cy - y);
                SetHoughPoint(cx - x, cy - y);
                SetHoughPoint(cx + y, cy + x);
                SetHoughPoint(cx - y, cy + x);
                SetHoughPoint(cx + y, cy - x);
                SetHoughPoint(cx - y, cy - x);
            }
        }
        /// <summary>
        /// Set point.
        /// </summary>
        /// <param name="x">Co. X</param>
        /// <param name="y">Co. Y</param>
        private void SetHoughPoint(int x, int y)
        {
            if ((x >= 0) && (y >= 0) && (x < width) && (y < height))
            {
                houghMap[y, x]++;
            }
        }
        #endregion
    }
}
