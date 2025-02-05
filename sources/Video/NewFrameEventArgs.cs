using System;
using SkiaDrawing;

namespace UMapx.Video
{
    /// <summary>
    /// Arguments for new frame event from video source.
    /// </summary>
    /// 
    public class NewFrameEventArgs : EventArgs
    {
        private readonly Bitmap frame;

        /// <summary>
        /// Initializes a new instance of the <see cref="NewFrameEventArgs"/> class.
        /// </summary>
        /// 
        /// <param name="frame">New frame.</param>
        /// 
        public NewFrameEventArgs(Bitmap frame)
        {
            this.frame = frame;
        }

        /// <summary>
        /// New frame from video source.
        /// </summary>
        /// 
        public Bitmap Frame
        {
            get { return frame; }
        }
    }
}
