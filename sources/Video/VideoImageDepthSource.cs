using System;
using SkiaDrawing;
using System.Timers;
using UMapx.Imaging;

namespace UMapx.Video
{
    /// <summary>
    /// Defines video image depth source.
    /// </summary>
    public class VideoImageDepthSource : IVideoDepthSource
    {
        #region Fields

        private readonly object _locker = new object();
        private VideoCapabilities _videoResolution;
        private VideoCapabilities _depthResolution;
        private Timer _timer;
        private readonly Bitmap _image;
        private readonly ushort[,] _depth;
        private int _framesReceived;
        private long _bytesReceived;

        #endregion

        #region Constructor

        /// <summary>
        /// Initializes video image depth source.
        /// </summary>
        /// <param name="image">Image</param>
        /// <param name="depth">Depth</param>
        public VideoImageDepthSource(Bitmap image, ushort[,] depth)
        {
            _image = image ?? throw new ArgumentNullException();
            _depth = depth ?? throw new ArgumentNullException();
            _videoResolution = 
                _depthResolution = 
                new VideoCapabilities(
                new Size(640, 480),
                30,
                30,
                32
            );
        }

        /// <summary>
        /// Elasped frame timer.
        /// </summary>
        /// <param name="sender">sender</param>
        /// <param name="e">e</param>
        private void OnElapsed(object sender, ElapsedEventArgs e)
        {
            lock (_locker)
            {
                // frame
                var frame = BitmapTransform.Resize(
                    _image,
                    new Size(_videoResolution.FrameSize.Width,
                    _videoResolution.FrameSize.Height));

                // depth
                var depth = DepthTransform.Resize(
                    _depth,
                    new Size(_depthResolution.FrameSize.Width,
                    _depthResolution.FrameSize.Height));

                OnNewFrame(frame);
                OnNewDepth(depth);
                frame?.Dispose();
            }
        }

        /// <summary>
        /// Called when video source gets new frame.
        /// </summary>
        /// <param name="frame">Frame</param>
        private void OnNewFrame(Bitmap frame)
        {
            _framesReceived++;
            _bytesReceived += frame.Width * frame.Height * (Image.GetPixelFormatSize(frame.PixelFormat) >> 3);
            NewFrame?.Invoke(this, new NewFrameEventArgs(frame));
        }

        /// <summary>
        /// Called when video source gets new depth.
        /// </summary>
        /// <param name="depth">Depth</param>
        private void OnNewDepth(ushort[,] depth)
        {
            NewDepth?.Invoke(this, new NewDepthEventArgs(depth));
        }

        #endregion

        #region Properties

        /// <summary>
        /// Video resolution to set.
        /// </summary>
        /// 
        /// <remarks><para>The property allows to set one of the video resolutions supported by the camera.
        /// Use <see cref="VideoCapabilities"/> property to get the list of supported video resolutions.</para>
        /// 
        /// <para><note>The property must be set before camera is started to make any effect.</note></para>
        /// 
        /// <para>Default value of the property is set to <see langword="null"/>, which means default video
        /// resolution is used.</para>
        /// </remarks>
        /// 
        public VideoCapabilities VideoResolution
        {
            get
            {
                return _videoResolution;
            }
            set
            {
                _videoResolution = value;
            }
        }

        /// <summary>
        /// Depth resolution to set.
        /// </summary>
        /// 
        /// <remarks><para>The property allows to set one of the depth resolutions supported by the camera.
        /// Use <see cref="VideoCapabilities"/> property to get the list of supported depth resolutions.</para>
        /// 
        /// <para><note>The property must be set before camera is started to make any effect.</note></para>
        /// 
        /// <para>Default value of the property is set to <see langword="null"/>, which means default depth
        /// resolution is used.</para>
        /// </remarks>
        /// 
        public VideoCapabilities DepthResolution
        {
            get
            {
                return _depthResolution;
            }
            set
            {
                _depthResolution = value;
            }
        }

        /// <summary>
        /// Returns video source.
        /// </summary>
        public virtual string Source
        {
            get
            {
                return "Video image depth source";
            }
        }

        /// <summary>
        /// Received frames count.
        /// </summary>
        /// 
        /// <remarks>Number of frames the video source provided from the moment of the last
        /// access to the property.
        /// </remarks>
        /// 
        public int FramesReceived
        {
            get
            {
                int frames = _framesReceived;
                _framesReceived = 0;
                return frames;
            }
        }

        /// <summary>
        /// Received bytes count.
        /// </summary>
        /// 
        /// <remarks>Number of bytes the video source provided from the moment of the last
        /// access to the property.
        /// </remarks>
        /// 
        public long BytesReceived
        {
            get
            {
                long bytes = _bytesReceived;
                _bytesReceived = 0;
                return bytes;
            }
        }

        /// <summary>
        /// State of the video source.
        /// </summary>
        /// 
        /// <remarks>Current state of video source object - running or not.</remarks>
        /// 
        public bool IsRunning { get; private set; }

        #endregion

        #region Events 

        /// <summary>
        /// Video depth action event handler.
        /// </summary>
        public event NewDepthEventHandler NewDepth;

        /// <summary>
        /// Video frame action event handler.
        /// </summary>
        public event NewFrameEventHandler NewFrame;

        /// <summary>
        /// Video source error event.
        /// </summary>
        /// 
        /// <remarks>This event is used to notify clients about any type of errors occurred in
        /// video source object, for example internal exceptions.</remarks>
        /// 
        public event VideoSourceErrorEventHandler VideoSourceError;

        /// <summary>
        /// Video playing finished event.
        /// </summary>
        /// 
        /// <remarks><para>This event is used to notify clients that the video playing has finished.</para>
        /// </remarks>
        /// 
        public event PlayingFinishedEventHandler PlayingFinished;

        #endregion

        #region Methods

        /// <summary>
        /// Start video source.
        /// </summary>
        /// 
        /// <remarks>Starts video source and return execution to caller. Video source
        /// object creates background thread and notifies about new frames with the
        /// help of <see cref="NewFrame"/> event.</remarks>
        /// 
        public void Start()
        {
            try
            {
                if (!IsRunning)
                {
                    _timer = new Timer();
                    _timer.Elapsed += OnElapsed;
                    _timer.Interval = 1000.0 / _videoResolution.MaximumFrameRate;
                    _timer?.Start();
                    IsRunning = true;
                }
            }
            catch (Exception ex)
            {
                VideoSourceError?.Invoke(this, new VideoSourceErrorEventArgs(ex.Message));
                throw;
            }
        }

        /// <summary>
        /// Signal video source to stop its work.
        /// </summary>
        /// 
        /// <remarks>Signals video source to stop its background thread, stop to
        /// provide new frames and free resources.</remarks>
        /// 
        public void SignalToStop()
        {
            if (IsRunning)
            {
                _timer?.Stop();
                _timer?.Dispose();
                IsRunning = false;
                PlayingFinished?.Invoke(this, ReasonToFinishPlaying.StoppedByUser);
            }
        }

        /// <summary>
        /// Stop video source.
        /// </summary>
        /// 
        /// <remarks>Not implemented</remarks>
        /// 
        [Obsolete]
        public void Stop()
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Wait for video source has stopped.
        /// </summary>
        /// 
        /// <remarks>Not implemented</remarks>
        [Obsolete]
        public void WaitForStop()
        {
            throw new NotImplementedException();
        }

        #endregion

        #region IDisposable

        private bool _disposed;

        /// <inheritdoc/>
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        /// <inheritdoc/>
        protected virtual void Dispose(bool disposing)
        {
            if (!_disposed)
            {
                if (disposing)
                {
                    _timer?.Dispose();
                    _image?.Dispose();
                }
                _disposed = true;
            }
        }

        /// <inheritdoc/>
        ~VideoImageDepthSource()
        {
            Dispose(false);
        }

        #endregion
    }
}
