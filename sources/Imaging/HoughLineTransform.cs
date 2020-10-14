using System;
using System.Collections;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the Hough line transform filter.
    /// </summary>
    [Serializable]
    public class HoughLineTransform : IBitmapFilter
    {
        #region Private data
        private int stepsPerDegree;
        private int houghHeight;
        private double thetaStep;

        private double[] sinMap;
        private double[] cosMap;

        private short[,] houghMap;
        private short maxMapIntensity = 0;

        private int localPeakRadius = 4;
        private short minLineIntensity = 10;
        private ArrayList lines = new ArrayList();
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Hough line transform filter.
        /// </summary>
        public HoughLineTransform()
        {
            StepsPerDegree = 1;
        }
        #endregion

        #region Class components
        /// <summary>
        /// Gets or sets steps per degree.
        /// </summary>
        public int StepsPerDegree
        {
            get { return stepsPerDegree; }
            set
            {
                stepsPerDegree = Math.Max(1, Math.Min(10, value));
                houghHeight = 180 * stepsPerDegree;
                thetaStep = Math.PI / houghHeight;

                // precalculate Sine and Cosine values
                sinMap = new double[houghHeight];
                cosMap = new double[houghHeight];

                for (int i = 0; i < houghHeight; i++)
                {
                    sinMap[i] = Math.Sin(i * thetaStep);
                    cosMap[i] = Math.Cos(i * thetaStep);
                }
            }
        }
        /// <summary>
        /// Gets the minimum line intensity.
        /// </summary>
        public short MinLineIntensity
        {
            get { return minLineIntensity; }
            set { minLineIntensity = value; }
        }
        /// <summary>
        /// Gets the maximum line intensity.
        /// </summary>
        public short MaxIntensity
        {
            get { return maxMapIntensity; }
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
        /// Gets the count of lines found.
        /// </summary>
        public int LinesCount
        {
            get { return lines.Count; }
        }
        /// <summary>
        /// Returns an array of lines with absolute intensity.
        /// </summary>
        /// <param name="count">Count</param>
        /// <returns>Array</returns>
        public HoughLine[] GetMostIntensiveLines(int count)
        {
            // lines count
            int n = Math.Min(count, lines.Count);

            // result array
            HoughLine[] dst = new HoughLine[n];
            lines.CopyTo(0, dst, 0, n);

            return dst;
        }
        /// <summary>
        /// Returns an array of lines with relative intensity.
        /// </summary>
        /// <param name="minRelativeIntensity">Minimum relative intensity</param>
        /// <returns>Array</returns>
        public HoughLine[] GetLinesByRelativeIntensity(double minRelativeIntensity)
        {
            int count = 0, n = lines.Count;

            while ((count < n) && (((HoughLine)lines[count]).RelativeIntensity >= minRelativeIntensity))
                count++;

            return GetMostIntensiveLines(count);
        }
        #endregion

        #region Image processing
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

            int height = houghMap.GetLength(0);
            int width = houghMap.GetLength(1);

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
            int width = bmData.Width;
            int height = bmData.Height;
            Rectangle rect = new Rectangle(0, 0, width, height);
            int halfWidth = width / 2;
            int halfHeight = height / 2;

            // make sure the specified rectangle recides with the source image
            rect.Intersect(new Rectangle(0, 0, width, height));

            int startX = -halfWidth + rect.Left;
            int startY = -halfHeight + rect.Top;
            int stopX = width - halfWidth - (width - rect.Right);
            int stopY = height - halfHeight - (height - rect.Bottom);

            int offset = bmData.Stride - rect.Width;

            // calculate Hough map's width
            int halfHoughWidth = (int)Math.Sqrt(halfWidth * halfWidth + halfHeight * halfHeight);
            int houghWidth = halfHoughWidth * 2;

            houghMap = new short[houghHeight, houghWidth];

            // do the job
            unsafe
            {
                byte* src = (byte*)bmData.Scan0.ToPointer() + rect.Top * bmData.Stride + rect.Left;
                int y, x, theta, radius;

                // for each row
                for (y = startY; y < stopY; y++)
                {
                    // for each pixel
                    for (x = startX; x < stopX; x++, src++)
                    {
                        if (*src != 0)
                        {
                            // for each Theta value
                            for (theta = 0; theta < houghHeight; theta++)
                            {
                                radius = (int)Math.Round(cosMap[theta] * x - sinMap[theta] * y) + halfHoughWidth;

                                if ((radius < 0) || (radius >= houghWidth))
                                    continue;

                                houghMap[theta, radius]++;
                            }
                        }
                    }
                    src += offset;
                }
            }

            // find max value in Hough map
            maxMapIntensity = 0;
            int i, j;

            for (i = 0; i < houghHeight; i++)
            {
                for (j = 0; j < houghWidth; j++)
                {
                    if (houghMap[i, j] > maxMapIntensity)
                    {
                        maxMapIntensity = houghMap[i, j];
                    }
                }
            }

            CollectLines();
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Collect lines with intesities greater or equal then specified.
        /// </summary>
        private void CollectLines()
        {
            int maxTheta = houghMap.GetLength(0);
            int maxRadius = houghMap.GetLength(1);

            short intensity;
            bool foundGreater;

            int halfHoughWidth = maxRadius >> 1;
            int theta, radius, cycledTheta, cycledRadius, tt, tr, ttMax, trMax;

            // clean lines collection
            lines.Clear();

            // for each Theta value
            for (theta = 0; theta < maxTheta; theta++)
            {
                // for each Radius value
                for (radius = 0; radius < maxRadius; radius++)
                {
                    // get current value
                    intensity = houghMap[theta, radius];

                    if (intensity < minLineIntensity)
                        continue;

                    foundGreater = false;

                    // check neighboors
                    for (tt = theta - localPeakRadius, ttMax = theta + localPeakRadius; tt < ttMax; tt++)
                    {
                        // break if it is not local maximum
                        if (foundGreater == true)
                            break;

                        cycledTheta = tt;
                        cycledRadius = radius;

                        // check limits
                        if (cycledTheta < 0)
                        {
                            cycledTheta = maxTheta + cycledTheta;
                            cycledRadius = maxRadius - cycledRadius;
                        }
                        if (cycledTheta >= maxTheta)
                        {
                            cycledTheta -= maxTheta;
                            cycledRadius = maxRadius - cycledRadius;
                        }

                        for (tr = cycledRadius - localPeakRadius, trMax = cycledRadius + localPeakRadius; tr < trMax; tr++)
                        {
                            // skip out of map values
                            if (tr < 0)
                                continue;
                            if (tr >= maxRadius)
                                break;

                            // compare the neighboor with current value
                            if (houghMap[cycledTheta, tr] > intensity)
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
                        lines.Add(new HoughLine((double)theta / stepsPerDegree, (short)(radius - halfHoughWidth), intensity, (double)intensity / maxMapIntensity));
                    }
                }
            }

            lines.Sort();
        }
        #endregion
    }
}
