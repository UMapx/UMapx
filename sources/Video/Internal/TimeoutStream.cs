namespace UMapx.Video
{
    using System;
    using System.IO;
    using System.Threading;

    /// <summary>
    /// Adds read and write deadlines to streams through cooperative asynchronous cancellation.
    /// Streams with native timeout support use their own synchronous timeout implementation.
    /// </summary>
    public class TimeoutStream : Stream
    {
        private const int DEFAULT_TIMEOUT_READ = 30000;
        private const int DEFAULT_TIMEOUT_WRITE = 30000;

        private readonly Stream _baseStream;
        private readonly CancellationToken _readCancellation;

        private int _readTimeout = DEFAULT_TIMEOUT_READ;
        private int _writeTimeout = DEFAULT_TIMEOUT_WRITE;

        /// <summary>
        /// Creates an instance of a TimeoutStream wrapper.
        /// </summary>
        /// <param name="stream">Stream which may not support read or write timeouts.</param>
        public TimeoutStream(Stream stream)
            : this(stream, CancellationToken.None)
        {
        }

        /// <summary>Allows a video source to cancel response reads when it stops.</summary>
        internal TimeoutStream(Stream stream, CancellationToken readCancellation)
        {
            _baseStream = stream ?? throw new ArgumentNullException(nameof(stream));
            _readCancellation = readCancellation;
#if NET35 || NET40
            throw new NotSupportedException();
#endif
        }

        /// <summary>
        /// Stream wrapped by TimeoutStream wrapper.
        /// </summary>
        public Stream BaseStream {
            get
            {
                return _baseStream;
            }
        }

        /// <summary>
        /// Pass-through property.
        /// </summary>
        public override bool CanRead {
            get
            {
                return _baseStream.CanRead;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public override bool CanSeek {
            get
            {
                return _baseStream.CanSeek;
            }
        }

        /// <summary>
        /// Pass-through property.
        /// </summary>
        public override bool CanWrite {
            get
            {
                return _baseStream.CanWrite;
            }
        }

        /// <summary>
        /// Pass-through property.
        /// </summary>
        public override long Length {
            get
            {
                return _baseStream.Length;
            }
        }

        /// <summary>
        /// Pass-through property.
        /// </summary>
        public override bool CanTimeout {
            get
            {
                return true;
            }
        }

        /// <summary>
        /// Gets or sets the read deadline in milliseconds, or -1 for no deadline.
        /// </summary>
        public override int ReadTimeout
        {
            get
            {
                return _readTimeout;
            }

            set
            {
                if (value <= 0 && value != Timeout.Infinite) throw new ArgumentOutOfRangeException(nameof(value));
                _readTimeout = value;
            }
        }

        /// <summary>
        /// Gets or sets the write deadline in milliseconds, or -1 for no deadline.
        /// </summary>
        public override int WriteTimeout
        {
            get
            {
                return _writeTimeout;
            }

            set
            {
                if (value <= 0 && value != Timeout.Infinite) throw new ArgumentOutOfRangeException(nameof(value));
                _writeTimeout = value;
            }
        }

        /// <summary>
        /// Pass-through property.
        /// </summary>
        public override long Position
        {
            get { return _baseStream.Position; }
            set { _baseStream.Position = value; }
        }

        /// <summary>
        /// Pass-through method.
        /// </summary>
        public override void Flush()
        {
            _baseStream.Flush();
        }

        /// <summary>
        /// Reads from base stream using a timeout.
        /// </summary>
        /// <param name="buffer">Buffer byte array.</param>
        /// <param name="offset">Offset.</param>
        /// <param name="count">Number of bytes to read.</param>
        /// <returns>The number of bytes read, or zero at the end of the stream.</returns>
        public override int Read(byte[] buffer, int offset, int count)
        {
#if !NET35 && !NET40
            if (_baseStream.CanRead && !_baseStream.CanTimeout)
            {
                // A deadline belongs to one operation. Disposing its timer prevents an idle
                // interval from cancelling the next read or a simultaneous write.
                using var source = _readCancellation.CanBeCanceled
                    ? CancellationTokenSource.CreateLinkedTokenSource(_readCancellation)
                    : new CancellationTokenSource();
                source.CancelAfter(_readTimeout);
                try
                {
                    return _baseStream.ReadAsync(buffer, offset, count, source.Token).GetAwaiter().GetResult();
                }
                catch (OperationCanceledException exception) when (source.IsCancellationRequested && !_readCancellation.IsCancellationRequested)
                {
                    throw new TimeoutException("The operation timed out.", exception);
                }
            }
            if (_baseStream.CanTimeout) _baseStream.ReadTimeout = _readTimeout;
            return _baseStream.Read(buffer, offset, count);
#else
            throw new NotSupportedException();
#endif
        }

        /// <summary>
        /// Pass-through method.
        /// </summary>
        /// <param name="offset"></param>
        /// <param name="origin"></param>
        /// <returns></returns>
        public override long Seek(long offset, SeekOrigin origin)
        {
            return _baseStream.Seek(offset, origin);
        }

        /// <summary>
        /// Pass-through method.
        /// </summary>
        /// <param name="value"></param>
        public override void SetLength(long value)
        {
            _baseStream.SetLength(value);
        }

        /// <summary>
        /// Write to base stream using a timeout.
        /// </summary>
        /// <param name="buffer">Buffer byte array.</param>
        /// <param name="offset">Offset.</param>
        /// <param name="count">Number of bytes to write.</param>
        public override void Write(byte[] buffer, int offset, int count)
        {
#if !NET35 && !NET40
            if (_baseStream.CanWrite && !_baseStream.CanTimeout)
            {
                using var source = new CancellationTokenSource();
                source.CancelAfter(_writeTimeout);
                try
                {
                    _baseStream.WriteAsync(buffer, offset, count, source.Token).GetAwaiter().GetResult();
                }
                catch (OperationCanceledException exception) when (source.IsCancellationRequested)
                {
                    throw new TimeoutException("The operation timed out.", exception);
                }
            }
            else
            {
                if (_baseStream.CanTimeout) _baseStream.WriteTimeout = _writeTimeout;
                _baseStream.Write(buffer, offset, count);
            }
#else
            throw new NotSupportedException();
#endif
        }
    }
}
