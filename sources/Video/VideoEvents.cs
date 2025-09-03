namespace UMapx.Video
{
    /// <summary>
    /// Delegate for new frame event handler.
    /// </summary>
    /// 
    /// <param name="sender">Sender object</param>
    /// <param name="eventArgs">Event arguments</param>
    /// 
    public delegate void NewFrameEventHandler(object sender, NewFrameEventArgs eventArgs);

    /// <summary>
    /// Delegate for video source error event handler.
    /// </summary>
    /// 
    /// <param name="sender">Sender object</param>
    /// <param name="eventArgs">Event arguments</param>
    /// 
    public delegate void VideoSourceErrorEventHandler(object sender, VideoSourceErrorEventArgs eventArgs);

    /// <summary>
    /// Delegate for playing finished event handler.
    /// </summary>
    /// 
    /// <param name="sender">Sender object</param>
    /// <param name="reason">Reason of finishing video playing</param>
    /// 
    public delegate void PlayingFinishedEventHandler(object sender, ReasonToFinishPlaying reason);

}
