using System;

namespace UMapx.Video
{
    /// <summary>
    /// Arguments for new depth event from video source.
    /// </summary>
    public class NewDepthEventArgs : EventArgs
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="NewDepthEventArgs"/> class.
        /// </summary>
        /// 
        /// <param name="depth">New depth</param>
        /// 
        public NewDepthEventArgs(ushort[,] depth)
        {
            Depth = depth;
        }

        /// <summary>
        /// Gets the depth.
        /// </summary>
        public ushort[,] Depth { get; private set; }
    }
}
