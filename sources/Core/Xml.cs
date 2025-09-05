using System;
using System.IO;
using System.Xml.Serialization;

namespace UMapx.Core
{
    /// <summary>
    /// Used to work with xml serialization of objects.
    /// </summary>
    public static class Xml
    {
        #region Binaries
        /// <summary>
        /// Save data to the file.
        /// </summary>
        /// <param name="stream">Stream</param>
        /// <param name="o">Object</param>
        public static void Save(Stream stream, object o)
        {
            XmlSerializer xml = new XmlSerializer(o.GetType());
            xml.Serialize(stream, o);
        }
        /// <summary>
        /// Save data to the file.
        /// </summary>
        /// <param name="fileName">File name</param>
        /// <param name="o">Object</param>
        public static void Save(string fileName, object o)
        {
            using var stream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            XmlSerializer xml = new XmlSerializer(o.GetType());
            xml.Serialize(stream, o);
        }
        /// <summary>
        /// Load data from the file.
        /// </summary>
        /// <param name="stream">Stream</param>
        /// <param name="type">Type</param>
        public static object Open(Stream stream, Type type)
        {
            XmlSerializer xml = new XmlSerializer(type);
            return xml.Deserialize(stream);
        }
        /// <summary>
        /// Load data from the file.
        /// </summary>
        /// <param name="fileName">File name</param>
        /// <param name="type">Type</param>
        public static object Open(string fileName, Type type)
        {
            using var stream = new FileStream(fileName, FileMode.Open, FileAccess.Read);
            XmlSerializer xml = new XmlSerializer(type);
            return xml.Deserialize(stream);
        }
        #endregion
    }
}
