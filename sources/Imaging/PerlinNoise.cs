using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the Perlin noise.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Perlin_noise
    /// </remarks>
    public class PerlinNoise
    {
        #region Private data
        private float initFrequency = 1.0f;
        private float initAmplitude = 1.0f;
        private float persistence = 0.65f;
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
        public PerlinNoise(int octaves = 4, float persistence = 0.65f, float frequency = 1, float amplitude = 1)
        {
            Octaves = octaves; Persistence = persistence; Frequency = frequency; Amplitude = amplitude;
        }
        /// <summary>
        /// Gets or sets the frequency.
        /// </summary>
        public float Frequency
        {
            get { return initFrequency; }
            set { initFrequency = value; }
        }
        /// <summary>
        /// Gets or sets the amplitude value.
        /// </summary>
        public float Amplitude
        {
            get { return initAmplitude; }
            set { initAmplitude = value; }
        }
        /// <summary>
        /// Gets or sets the persistence value.
        /// </summary>
        public float Persistence
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
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Function(float x)
        {
            float frequency = initFrequency;
            float amplitude = initAmplitude;
            float sum = 0;

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
        /// <param name="x">Value</param>
        /// <param name="y">Value</param>
        /// <returns>Value</returns>
        public float Function2D(float x, float y)
        {
            float frequency = initFrequency;
            float amplitude = initAmplitude;
            float sum = 0;

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
        /// Generates a deterministic pseudo-random noise value for one dimension.
        /// </summary>
        /// <param name="x">Input coordinate</param>
        /// <returns>Noise value in the range [-1, 1]</returns>
        private static float Noise(int x)
        {
            int n = (x << 13) ^ x;

            return (1.0f - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0f);
        }
        /// <summary>
        /// Generates a deterministic pseudo-random noise value for two dimensions.
        /// </summary>
        /// <param name="x">X coordinate</param>
        /// <param name="y">Y coordinate</param>
        /// <returns>Noise value in the range [-1, 1]</returns>
        private static float Noise(int x, int y)
        {
            int n = x + y * 57;
            n = (n << 13) ^ n;

            return (1.0f - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0f);
        }
        /// <summary>
        /// Computes smoothed 1D noise using cosine interpolation.
        /// </summary>
        /// <param name="x">Input coordinate</param>
        /// <returns>Smoothed noise value</returns>
        private static float SmoothedNoise(float x)
        {
            int xInt = (int)x;
            float xFrac = x - xInt;

            return PerlinNoise.CosineInterpolate(Noise(xInt), Noise(xInt + 1), xFrac);
        }
        /// <summary>
        /// Computes smoothed 2D noise using cosine interpolation.
        /// </summary>
        /// <param name="x">X coordinate</param>
        /// <param name="y">Y coordinate</param>
        /// <returns>Smoothed noise value</returns>
        private static float SmoothedNoise(float x, float y)
        {
            // params
            int xInt = (int)x;
            int yInt = (int)y;
            float xFrac = x - xInt;
            float yFrac = y - yInt;

            // get four noise values
            float x0y0 = PerlinNoise.Noise(xInt, yInt);
            float x1y0 = PerlinNoise.Noise(xInt + 1, yInt);
            float x0y1 = PerlinNoise.Noise(xInt, yInt + 1);
            float x1y1 = PerlinNoise.Noise(xInt + 1, yInt + 1);

            // x interpolation
            float v1 = PerlinNoise.CosineInterpolate(x0y0, x1y0, xFrac);
            float v2 = PerlinNoise.CosineInterpolate(x0y1, x1y1, xFrac);

            // y interpolation
            return PerlinNoise.CosineInterpolate(v1, v2, yFrac);
        }
        /// <summary>
        /// Performs cosine interpolation between two values.
        /// </summary>
        /// <param name="x1">First value</param>
        /// <param name="x2">Second value</param>
        /// <param name="a">Interpolation factor</param>
        /// <returns>Interpolated value</returns>
        private static float CosineInterpolate(float x1, float x2, float a)
        {
            float f = (1 - Maths.Cos(a * Maths.Pi)) * 0.5f;

            return x1 * (1 - f) + x2 * f;
        }
        #endregion
    }
}
