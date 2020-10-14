using System.IO;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;

namespace UMapx.Core
{
    /// <summary>
    /// Uses for binary serialization of objects.
    /// </summary>
    public static class Binary
    {
        #region Binaries
        /// <summary>
        /// Save data from the file.
        /// </summary>
        /// <param name="stream">Stream</param>
        /// <param name="o">Object</param>
        public static void Save(Stream stream, object o)
        {
            IFormatter bin = new BinaryFormatter();
            bin.Serialize(stream, o);
        }
        /// <summary>
        /// Save data from the file.
        /// </summary>
        /// <param name="fileName">File name</param>
        /// <param name="o">Object</param>
        public static void Save(string fileName, object o)
        {
            FileStream stream = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            IFormatter bin = new BinaryFormatter();
            bin.Serialize(stream, o);
            stream.Close();
            stream.Dispose();
        }
        /// <summary>
        /// Load data from the file.
        /// </summary>
        /// <param name="stream">Stream</param>
        public static object Load(Stream stream)
        {
            IFormatter bin = new BinaryFormatter();
            return bin.Deserialize(stream);
        }
        /// <summary>
        /// Load data from the file.
        /// </summary>
        /// <param name="fileName">File name</param>
        public static object Load(string fileName)
        {
            FileStream stream = new FileStream(fileName, FileMode.Open, FileAccess.Read);
            IFormatter bin = new BinaryFormatter();
            object graph = bin.Deserialize(stream);
            stream.Close();
            stream.Dispose();
            return graph;
        }
        #endregion
    }
}
