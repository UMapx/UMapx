namespace UMapx.Video
{
    /// <summary>
    /// Reason of finishing video playing.
    /// </summary>
    /// 
    /// <remarks><para>When video source class fire the <see cref="IVideoSource.PlayingFinished"/> event, they
    /// need to specify reason of finishing video playing. For example, it may be end of stream reached.</para></remarks>
    /// 
    public enum ReasonToFinishPlaying
    {
        /// <summary>
        /// Video playing has finished because it end was reached.
        /// </summary>
        EndOfStreamReached,
        /// <summary>
        /// Video playing has finished because it was stopped by user.
        /// </summary>
        StoppedByUser,
        /// <summary>
        /// Video playing has finished because the device was lost (unplugged).
        /// </summary>
        DeviceLost,
        /// <summary>
        /// Video playing has finished because of some error happened the video source (camera, stream, file, etc.).
        /// A error reporting event usually is fired to provide error information.
        /// </summary>
        VideoSourceError
    }
}
