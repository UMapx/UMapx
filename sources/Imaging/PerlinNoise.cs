using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the Perlin noise.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Perlin_noise
    /// </remarks>
    /// </summary>
    public class PerlinNoise
    {
        #region Private data
        private double initFrequency = 1.0;
        private double initAmplitude = 1.0;
        private double persistence = 0.65;
        private int octaves = 4;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes the Perlin noise.   
        /// </summary>
        /// <param name="octaves">Octaves[1, 32]</param>
        /// <param name="persistence">Persistence</param>
        /// <param name="frequency">Frequency</param>
        /// <param name="amplitude">Amplitude</param>
        public PerlinNoise(int octaves = 4, double persistence = 0.65, double frequency = 1, double amplitude = 1)
        {
            Octaves = octaves; Persistence = persistence; Frequency = frequency; Amplitude = amplitude;
        }
        /// <summary>
        /// Gets or sets the frequency.
        /// </summary>
        public double Frequency
        {
            get { return initFrequency; }
            set { initFrequency = value; }
        }
        /// <summary>
        /// Gets or sets the amplitude value.
        /// </summary>
        public double Amplitude
        {
            get { return initAmplitude; }
            set { initAmplitude = value; }
        }
        /// <summary>
        /// Gets or sets the persistence value.
        /// </summary>
        public double Persistence
        {
            get { return persistence; }
            set { persistence = value; }
        }
        /// <summary>
        /// Gets or sets the number of octaves[1, 32].
        /// </summary>
        public int Octaves
        {
            get { return octaves; }
            set { octaves = System.Math.Max(1, System.Math.Min(32, value)); }
        }
        /// <summary>
        /// One-dimensional Perlin noise function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            double frequency = initFrequency;
            double amplitude = initAmplitude;
            double sum = 0;

            // octaves
            for (int i = 0; i < octaves; i++)
            {
                sum += PerlinNoise.SmoothedNoise(x * frequency) * amplitude;

                frequency *= 2;
                amplitude *= persistence;
            }
            return sum;
        }
        /// <summary>
        /// Two-dimensional Perlin noise function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="y">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public double Function2D(double x, double y)
        {
            double frequency = initFrequency;
            double amplitude = initAmplitude;
            double sum = 0;

            // octaves
            for (int i = 0; i < octaves; i++)
            {
                sum += PerlinNoise.SmoothedNoise(x * frequency, y * frequency) * amplitude;

                frequency *= 2;
                amplitude *= persistence;
            }
            return sum;
        }
        #endregion

        #region Private static voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private static double Noise(int x)
        {
            int n = (x << 13) ^ x;

            return (1.0 - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        private static double Noise(int x, int y)
        {
            int n = x + y * 57;
            n = (n << 13) ^ n;

            return (1.0 - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private static double SmoothedNoise(double x)
        {
            int xInt = (int)x;
            double xFrac = x - xInt;

            return PerlinNoise.CosineInterpolate(Noise(xInt), Noise(xInt + 1), xFrac);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        private static double SmoothedNoise(double x, double y)
        {
            // params
            int xInt = (int)x;
            int yInt = (int)y;
            double xFrac = x - xInt;
            double yFrac = y - yInt;

            // get four noise values
            double x0y0 = PerlinNoise.Noise(xInt, yInt);
            double x1y0 = PerlinNoise.Noise(xInt + 1, yInt);
            double x0y1 = PerlinNoise.Noise(xInt, yInt + 1);
            double x1y1 = PerlinNoise.Noise(xInt + 1, yInt + 1);

            // x interpolation
            double v1 = PerlinNoise.CosineInterpolate(x0y0, x1y0, xFrac);
            double v2 = PerlinNoise.CosineInterpolate(x0y1, x1y1, xFrac);

            // y interpolation
            return PerlinNoise.CosineInterpolate(v1, v2, yFrac);
        }
        ///
        private static double CosineInterpolate(double x1, double x2, double a)
        {
            double f = (1 - Math.Cos(a * Math.PI)) * 0.5;

            return x1 * (1 - f) + x2 * f;
        }
        #endregion
    }
}
