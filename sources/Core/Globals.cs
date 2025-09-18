namespace UMapx.Core
{
    /// <summary>
    /// Used to configure the UMapx environment.
    /// </summary>
    public static class Globals
    {
        /// <summary>
        /// Use SIMD optimization or not.
        /// </summary>
        public static bool SIMD = false;
        /// <summary>
        /// Numbers format (works for Complex32 and Quaternion32).
        /// </summary>
        public static string DefaultFormat = "G6";
    }
}
