namespace UMapx.Video
{
    /// <summary>
    /// Extension interface for handling Intel RealSense depth events.
    /// </summary>
    public interface IVideoDepthSource : IVideoSource
    {
        /// <summary>
        /// Handler of received frames
        /// </summary>
        public event NewDepthEventHandler NewDepth;
    }
}
