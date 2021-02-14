using System;

namespace UMapx.Core
{
    /// <summary>
    /// Uses for tic/toc operations.
    /// </summary>
    public static class TicToc
    {
        #region Class components
        /// <summary>
        /// Private member.
        /// </summary>
        private static int tic;
        /// <summary>
        /// Starts elapsing time.
        /// </summary>
        public static void Tic()
        {
            tic = Environment.TickCount;
        }
        /// <summary>
        /// Returns elapsed time in miliseconds.
        /// </summary>
        /// <returns>Int</returns>
        public static int Toc()
        {
            return Environment.TickCount - tic;
        }
        #endregion
    }
}
