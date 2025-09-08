﻿namespace UMapx.Video
{
    using System;
    using System.IO;
    using System.Threading;
    using System.Threading.Tasks;

    /// <summary>
    /// Wrapper that enables streams which do not support read and write timeouts to timeout
    /// Requires .NET 4.5 or above
    /// Used for .NET Standard 2.0
    /// </summary>
    public class TimeoutStream : Stream
    {
        private const int DEFAULT_TIMEOUT_READ = 30000;
        private const int DEFAULT_TIMEOUT_WRITE = 30000;

        private Stream _baseStream;
#if !NET35 && !NET40
        private CancellationTokenSource _source;
#endif

        private int _readTimeout = DEFAULT_TIMEOUT_READ;
        private int _writeTimeout = DEFAULT_TIMEOUT_WRITE;

        /// <summary>
        /// Creates an instance of a TimeoutStream wrapper
        /// </summary>
        /// <param name="stream">Stream which may not support read or write timeouts</param>
        public TimeoutStream(Stream stream)
        {
            _baseStream = stream;
#if !NET35 && !NET40
            _source = new CancellationTokenSource();
#else
            throw new NotSupportedException();
#endif
        }

        /// <summary>
        /// Stream wrapped by TimeoutStream wrapper
        /// </summary>
        public Stream BaseStream {
            get
            {
                return _baseStream;
            }
        }

        /// <summary>
        /// Pass-through property
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
        /// Pass-through property
        /// </summary>
        public override bool CanWrite {
            get
            {
                return _baseStream.CanWrite;
            }
        }

        /// <summary>
        /// Pass-through property
        /// </summary>
        public override long Length {
            get
            {
                return _baseStream.Length;
            }
        }

        /// <summary>
        /// Pass-through property
        /// </summary>
        public override bool CanTimeout {
            get
            {
                return true;
            }
        }

        /// <summary>
        /// Value of TimeoutStream's read timeout
        /// </summary>
        public override int ReadTimeout
        {
            get
            {
                return _readTimeout;
            }

            set
            {
                _readTimeout = value;
            }
        }

        /// <summary>
        /// Value of TimeoutStream's write timeout
        /// </summary>
        public override int WriteTimeout
        {
            get
            {
                return _writeTimeout;
            }

            set
            {
                _writeTimeout = value;
            }
        }

        /// <summary>
        /// Pass-through property
        /// </summary>
        public override long Position
        {
            get { return _baseStream.Position; }
            set { _baseStream.Position = value; }
        }

        /// <summary>
        /// Pass-through method
        /// </summary>
        public override void Flush()
        {
            _baseStream.Flush();
        }

        /// <summary>
        /// Reads from base stream using a timeout
        /// </summary>
        /// <param name="buffer">Buffer byte array</param>
        /// <param name="offset">Offset</param>
        /// <param name="count">Number of bytes to read</param>
        /// <returns></returns>
        public override int Read(byte[] buffer, int offset, int count)
        {
#if !NET35 && !NET40
            int result;

            if (_baseStream.CanRead && !_baseStream.CanTimeout)
            {
                try
                {
                    _source.CancelAfter(_readTimeout);
                    Task<int> task = _baseStream.ReadAsync(buffer, offset, count,
                        _source.Token);
                    result = task.Result;
                }
                catch (AggregateException)
                {
                    _source = new CancellationTokenSource();
                    throw new TimeoutException("The operation timed out");
                }
            }
            else
            {
                result = _baseStream.Read(buffer, offset, count);
            }

            return result;
#else
            throw new NotSupportedException();
#endif
        }

        /// <summary>
        /// Pass-through method
        /// </summary>
        /// <param name="offset"></param>
        /// <param name="origin"></param>
        /// <returns></returns>
        public override long Seek(long offset, SeekOrigin origin)
        {
            return _baseStream.Seek(offset, origin);
        }

        /// <summary>
        /// Pass-through method
        /// </summary>
        /// <param name="value"></param>
        public override void SetLength(long value)
        {
            _baseStream.SetLength(value);
        }

        /// <summary>
        /// Write to base stream using a timeout
        /// </summary>
        /// <param name="buffer">Buffer byte array</param>
        /// <param name="offset">Offset</param>
        /// <param name="count">Number of bytes to write</param>
        public override void Write(byte[] buffer, int offset, int count)
        {
#if !NET35 && !NET40
            if (_baseStream.CanWrite && !_baseStream.CanTimeout)
            {
                try
                {
                    _source.CancelAfter(_readTimeout);
                    Task task = _baseStream.WriteAsync(buffer, offset, count, _source.Token);
                    task.Wait();
                }
                catch (AggregateException)
                {
                    _source?.Dispose();
                    _source = new CancellationTokenSource();
                    throw new TimeoutException("The operation timed out");
                }
            }
            else
            {
                _baseStream.Write(buffer, offset, count);
            }
#else
            throw new NotSupportedException();
#endif
        }
    }
}
