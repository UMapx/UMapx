namespace UMapx.Video
{
    using System;
    using System.Drawing;
    using System.IO;

    /// <summary>
    /// Handles functionality related to parsing a MJPEG stream.
    /// </summary>
    public class MJPEGStreamParser
    {
        private const int READ_SIZE = 1024;
        private const int BUFFER_SIZE = READ_SIZE * READ_SIZE;

        private byte[] _buffer;

        private int _position = 0;
        private int _totalReadBytes = 0;

        private int _imageHeaderIndex = -1;
        private int _imageBoundaryIndex = -1;

        private readonly byte[] _header;
        private readonly Boundary _boundary;

        /// <summary>
        /// Creates instance of MJPEG stream parser using a boundary and a JPEG magic header.
        /// </summary>
        /// <param name="boundary"></param>
        /// <param name="header"></param>
        /// <param name="bufferSize">Initial buffer capacity. The buffer grows to fit incoming frames.</param>
        public MJPEGStreamParser(Boundary boundary, byte[] header, int bufferSize = BUFFER_SIZE)
        {
            _header = header;
            _boundary = boundary;

            _buffer = new byte[bufferSize];
        }

        /// <summary>
        /// Content of the current byte array buffer. Reading may replace the array when it grows.
        /// </summary>
        public byte[] Content
        {
            get { return _buffer; }
        }

        private int RemainingBytes
        {
            get { return _totalReadBytes - _position; }
        }

        private bool HasStart
        {
            get { return _imageHeaderIndex != -1; }
        }

        private bool HasEnd
        {
            get { return _imageBoundaryIndex != -1; }
        }

        /// <summary>
        /// True if frame is detected using DetectFrame and not removed using RemoveFrame.
        /// </summary>
        public bool HasFrame
        {
            get { return HasStart && HasEnd; }
        }

        /// <summary>
        /// Reads byte content to internal buffer from a stream.
        /// </summary>
        /// <param name="stream"></param>
        /// <returns></returns>
        public int Read(Stream stream)
        {
            EnsureBufferCapacity();

            int readBytes = stream.Read(_buffer, _totalReadBytes, READ_SIZE);

            if (readBytes == 0)
                throw new ApplicationException();

            _totalReadBytes += readBytes;

            return readBytes;
        }

        /// <summary>
        /// Reserves space for the next read without discarding frame data or search positions.
        /// </summary>
        private void EnsureBufferCapacity()
        {
            int required = checked(_totalReadBytes + READ_SIZE);
            if (required > _buffer.Length)
            {
                int capacity = Math.Max(required, (int)Math.Min((long)_buffer.Length * 2, int.MaxValue));
                Array.Resize(ref _buffer, capacity);
            }
        }

        /// <summary>
        /// Detects if a frame is present in the internal buffer.
        /// </summary>
        public void DetectFrame()
        {
            if (!HasStart && CanRead(_header))
            {
                _imageHeaderIndex = FindHeader();

                if (HasStart)
                {
                    PositionAfterHeader();
                }
                else
                {
                    PositionAtEnd();
                }
            }

            while (HasStart && !HasEnd && CanRead(_boundary.HasValue ? (byte[])_boundary : _header))
            {
                _imageBoundaryIndex = FindBoundary();

                if (!HasEnd)
                {
                    PositionAtEnd();
                }
            }
        }

        /// <summary>
        /// Retrieves the frame from the internal buffer.
        /// </summary>
        /// <returns></returns>
        public Bitmap GetFrame()
        {
            if (HasFrame)
            {
                PositionAtImageEnd();

                int length = _imageBoundaryIndex - _imageHeaderIndex;
                Stream imageStream = new MemoryStream(_buffer, _imageHeaderIndex, length);
                return (Bitmap)Image.FromStream(imageStream);
            }
            else
            {
                throw new InvalidOperationException("No frame detected in buffer");
            }
        }

        /// <summary>
        /// Removes current frame from buffer.
        /// </summary>
        public void RemoveFrame()
        {
            if (HasFrame)
            {
                _position = _imageBoundaryIndex + _boundary.Length;
                Array.Copy(_buffer, _position, _buffer, 0, RemainingBytes);

                _totalReadBytes = RemainingBytes;
                _position = 0;

                _imageHeaderIndex = -1;
                _imageBoundaryIndex = -1;
            }
            else
            {
                throw new InvalidOperationException("No frame detected in buffer");
            }
        }

        /// <summary>
        /// Retains a possible partial marker at the end of the buffer for the next read.
        /// </summary>
        private void PositionAtEnd()
        {
            int markerLength = HasStart && _boundary.HasValue ? _boundary.Length : _header.Length;
            // A partial match can occupy at most markerLength - 1 trailing bytes.
            _position = Math.Max(0, _totalReadBytes - markerLength + 1);
        }

        /// <summary>
        /// Searches for the JPEG frame header within the buffer.
        /// </summary>
        /// <returns>Index of the header or -1.</returns>
        private int FindHeader()
        {
            return ByteArrayUtils.Find(_buffer, _header, _position, RemainingBytes);
        }

        /// <summary>
        /// Searches for the boundary marker in the buffer.
        /// </summary>
        /// <returns>Index of boundary or -1.</returns>
        private int FindBoundary()
        {
            byte[] imageDelimiter;

            if (_boundary.Length != 0)
            {
                imageDelimiter = (byte[])_boundary;
            }
            else
            {
                imageDelimiter = _header;
            }

            return ByteArrayUtils.Find(_buffer, imageDelimiter, _position, RemainingBytes);
        }

        internal int FindImageBoundary()
        {
            return ByteArrayUtils.Find(_buffer, (byte[])_boundary, 0, RemainingBytes);
        }

        /// <summary>
        /// Moves the current position to the end of the detected image.
        /// </summary>
        private void PositionAtImageEnd()
        {
            _position = _imageBoundaryIndex;
        }

        /// <summary>
        /// Advances the current position just after the JPEG header.
        /// </summary>
        private void PositionAfterHeader()
        {
            _position = _imageHeaderIndex + _header.Length;
        }

        internal bool CanRead(Boundary boundary)
        {
            byte[] target = (byte[])boundary;
            return CanRead(target);
        }

        internal bool CanRead(byte[] target)
        {
            return RemainingBytes != 0 && RemainingBytes >= target.Length;
        }
    }
}
