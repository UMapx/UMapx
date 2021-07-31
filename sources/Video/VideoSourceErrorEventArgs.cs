using System;

namespace UMapx.Video
{
    /// <summary>
    /// Arguments for video source error event from video source.
    /// </summary>
    /// 
    public class VideoSourceErrorEventArgs : EventArgs
    {
        private readonly string _description;
        private readonly Exception _exception;

        /// <summary>
        /// Initializes a new instance of the <see cref="VideoSourceErrorEventArgs"/> class.
        /// </summary>
        /// <param name="description">Error description.</param>
        public VideoSourceErrorEventArgs(string description)
            : this(description, null) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="VideoSourceErrorEventArgs"/> class.
        /// </summary>
        /// 
        /// <param name="description">Error description.</param>
        /// <param name="exception">Error exception.</param>
        /// 
        public VideoSourceErrorEventArgs(string description, Exception exception)
        {
            _description = description;
            _exception = exception;
        }

        /// <summary>
        /// Video source error description.
        /// </summary>
        /// 
        public string Description
        {
            get { return _description; }
        }

        /// <summary>
        /// Video source exception causing the error
        /// </summary>
        public Exception Exception
        {
            get { return _exception; }
        }
    }
}
