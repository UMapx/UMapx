namespace UMapx.Video
{
    using System;
    using System.Text;
    using System.Net;
    using System.Net.Http.Headers;

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
            if (response == null) throw new ArgumentNullException(nameof(response));
            // Unfold legacy response headers; other line breaks remain invalid HTTP syntax.
            string header = response.ContentType?.Replace("\r\n ", " ").Replace("\r\n\t", " ");
            if (!MediaTypeHeaderValue.TryParse(header, out MediaTypeHeaderValue contentType))
                throw new ArgumentException("Invalid content type", nameof(response));

            if (contentType.MediaType.Equals("multipart/x-mixed-replace", StringComparison.OrdinalIgnoreCase) ||
                contentType.MediaType.Equals("multipart/mixed", StringComparison.OrdinalIgnoreCase))
            {
                string boundary = null;
                foreach (var parameter in contentType.Parameters)
                {
                    if (!parameter.Name.Equals("boundary", StringComparison.OrdinalIgnoreCase)) continue;
                    if (boundary != null || parameter.Value == null)
                        throw new ArgumentException("Invalid boundary parameter", nameof(response));
                    boundary = UnquoteBoundary(parameter.Value);
                }
                return new Boundary(boundary ?? string.Empty);
            }

            if (contentType.MediaType.Equals("application/octet-stream", StringComparison.OrdinalIgnoreCase))
                return new Boundary();

            throw new ArgumentException("Invalid content type", nameof(response));
        }

        /// <summary>
        /// Removes HTTP quoting and quoted-pair escapes from a parsed boundary value.
        /// </summary>
        /// <param name="value">Token or quoted string validated by the header parser</param>
        /// <returns>Boundary string with its original case preserved</returns>
        private static string UnquoteBoundary(string value)
        {
            if (value.Length == 0 || value[0] != '"') return value;
            var boundary = new StringBuilder(value.Length - 2);
            for (int i = 1; i < value.Length - 1; i++)
            {
                if (value[i] == '\\') i++;
                boundary.Append(value[i]);
            }
            return boundary.ToString();
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
