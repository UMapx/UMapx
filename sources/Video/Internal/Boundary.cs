﻿namespace UMapx.Video
{
    using System;
    using System.Text;
    using System.Net;

    /// <summary>
    /// Container for MJPEG stream boundaries
    /// </summary>
    public class Boundary
    {
        private readonly static Encoding _encoding = Encoding.ASCII;

        private readonly StringBuilder _builder;

        private bool _isChecked = false;

        /// <summary>
        /// Creates an empty boundary for e.g. octet streams
        /// </summary>
        public Boundary()
        {
            _builder = new StringBuilder();
        }

        /// <summary>
        /// Creates instance using a string as boundary for e.g. multipart streams
        /// </summary>
        /// <param name="boundary">Boundary string</param>
        public Boundary(string boundary)
        {
            _builder = new StringBuilder(boundary);
        }

        /// <summary>
        /// Boundary string content
        /// </summary>
        public string Content
        {
            get { return _builder.ToString(); }
        }

        /// <summary>
        /// Length of boundary string
        /// </summary>
        public int Length
        {
            get { return _builder.Length; }
        }

        /// <summary>
        /// True if boundary string length is non-zero
        /// </summary>
        public bool HasValue
        {
            get { return Length != 0; }
        }

        /// <summary>
        /// True if FixMalformedBoundary has been run
        /// </summary>
        public bool IsChecked
        {
            get { return _isChecked; }
            set { _isChecked = value; }
        }

        /// <summary>
        /// True if IsChecked is true and HasValue is true, or if HasValue is false
        /// </summary>
        public bool IsValid
        {
            get { return (IsChecked && HasValue) || !HasValue; }
        }

        /// <summary>
        /// Adds character before boundary content
        /// </summary>
        /// <param name="c"></param>
        public void Prepend(char c)
        {
            _builder.Insert(0, c);
        }

        /// <summary>
        /// Some IP cameras, like AirLink, claim that boundary is "myboundary",
        /// when it is really "--myboundary". This corrects the issue.
        /// </summary>
        /// <param name="streamParser"></param>
        public void FixMalformedBoundary(MJPEGStreamParser streamParser)
        {
            byte[] content = streamParser.Content;

            int boundaryIndex = streamParser.FindImageBoundary();

            if (boundaryIndex != -1)
            {
                for (int i = boundaryIndex - 1; i >= 0; i--)
                {
                    char ch = (char)content[i];

                    if (ch == '\n' || ch == '\r')
                    {
                        break;
                    }

                    Prepend(ch);
                }

                IsChecked = true;
            }
        }

        /// <summary>
        /// Creates boundary from WebResponse
        /// </summary>
        /// <param name="response">Source of boundary string</param>
        /// <returns>Boundary with string content</returns>
        public static Boundary FromResponse(WebResponse response)
        {
            string contentType = response.ContentType;

            Boundary boundary = new Boundary();

            if (IsMultipartContent(contentType))
            {
                int boundaryIndex = GetBoundaryIndex(contentType);

                if (boundaryIndex != -1)
                {
                    boundary = TrimBoundary(contentType, boundaryIndex);
                    boundary.IsChecked = false;
                }
            }
            else if (!IsOctetStream(contentType))
            {
                throw new ArgumentException("Invalid content type.");
            }

            return boundary;
        }

        /// <summary>
        /// Locates the boundary parameter position in a content type string.
        /// </summary>
        /// <param name="contentType">Content type header</param>
        /// <returns>Index of boundary parameter or -1.</returns>
        private static int GetBoundaryIndex(string contentType)
        {
            int boundaryIndex = contentType.IndexOf("boundary", 0);
            if (boundaryIndex != -1)
            {
                boundaryIndex = contentType.IndexOf("=", boundaryIndex + 8);
            }

            return boundaryIndex;
        }

        /// <summary>
        /// Extracts and trims the boundary value from a content type string.
        /// </summary>
        /// <param name="contentType">Content type header</param>
        /// <param name="boundaryIndex">Index of the boundary token</param>
        /// <returns>Parsed boundary instance.</returns>
        private static Boundary TrimBoundary(string contentType, int boundaryIndex)
        {
            string boundary = contentType.Substring(boundaryIndex + 1);
            string trimmedBoundary = boundary.Trim(' ', '"');

            return new Boundary(trimmedBoundary);
        }

        /// <summary>
        /// Determines whether the content type is multipart/mixed.
        /// </summary>
        /// <param name="contentType">Content type header</param>
        /// <returns><c>true</c> if multipart/mixed; otherwise, false.</returns>
        private static bool IsMultipartContent(string contentType)
        {
            return contentType.StartsWith("multipart") && contentType.Contains("mixed");
        }

        /// <summary>
        /// Determines whether the content type is an octet-stream.
        /// </summary>
        /// <param name="contentType">Content type header</param>
        /// <returns><c>true</c> if octet-stream; otherwise, false.</returns>
        private static bool IsOctetStream(string contentType)
        {
            return contentType.StartsWith("application/octet-stream");
        }

        /// <summary>
        /// Converts boundary to string
        /// </summary>
        /// <param name="boundary">Boundary string content</param>
        public static explicit operator string(Boundary boundary)
        {
            string content = null;

            if (boundary != null)
            {
                content = boundary.Content;
            }

            return content;
        }

        /// <summary>
        /// Converts boundary to byte array
        /// </summary>
        /// <param name="boundary">Boundary byte content</param>
        public static explicit operator byte[] (Boundary boundary)
        {
            byte[] content = null;

            if (boundary != null)
            {
                content = _encoding.GetBytes(boundary.Content);
            }

            return content;
        }
    }
}
