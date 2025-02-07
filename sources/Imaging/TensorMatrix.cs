using SkiaDrawing;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Uses to work with tensor matrices.
    /// </summary>
    public static class TensorMatrix
    {
        #region Tensor convert components
        /// <summary>
        /// Converts a Bitmap to an BGR tensor arrays.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="rgb">RGB or BGR</param>
        /// <returns>RGB tensor arrays</returns>
        public static byte[][] ToByteTensor(this Bitmap Data, bool rgb = false)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            byte[][] _ix = TensorMatrix.ToByteTensor(bmData, rgb);
            BitmapFormat.Unlock(Data, bmData);
            return _ix;
        }
        /// <summary>
        /// Converts a Bitmap to an BGR tensor arrays.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="rgb">RGB or BGR</param>
        /// <returns>RGB tensor arrays</returns>
        public unsafe static byte[][] ToByteTensor(this BitmapData bmData, bool rgb = false)
        {
            // params
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int shift = height * width;
            byte[] _ix0 = new byte[shift];
            byte[] _ix1 = new byte[shift];
            byte[] _ix2 = new byte[shift];
            int z = 0;

            // do job
            if (rgb)
            {
                for (int j = 0; j < height; j++)
                {
                    for (int i = 0; i < width; i++, z++)
                    {
                        int k, jstride = j * stride;
                        k = jstride + i * 4;

                        // transform
                        _ix0[z] = p[k + 2];
                        _ix1[z] = p[k + 1];
                        _ix2[z] = p[k + 0];
                    }
                }
            }
            else
            {
                for (int j = 0; j < height; j++)
                {
                    for (int i = 0; i < width; i++, z++)
                    {
                        int k, jstride = j * stride;
                        k = jstride + i * 4;

                        // transform
                        _ix0[z] = p[k + 0];
                        _ix1[z] = p[k + 1];
                        _ix2[z] = p[k + 2];
                    }
                }
            }

            // arrays
            return new byte[][] { _ix0, _ix1, _ix2 };
        }
        /// <summary>
        /// Converts a Bitmap to an BGR tensor arrays.
        /// </summary>
        /// <param name="bmData">Bitmap data in BGR terms</param>
        /// <param name="rgb">RGB or BGR</param>
        /// <returns>RGB tensor arrays</returns>
        public unsafe static byte[][] ToByteTensor(this float[][,] bmData, bool rgb = false)
        {
            // params
            int width = bmData[0].GetLength(1), height = bmData[1].GetLength(0);
            int shift = height * width;
            byte[] _ix0 = new byte[shift];
            byte[] _ix1 = new byte[shift];
            byte[] _ix2 = new byte[shift];
            int z = 0;

            // do job
            if (rgb)
            {
                for (int j = 0; j < height; j++)
                {
                    for (int i = 0; i < width; i++, z++)
                    {
                        // transform
                        _ix0[z] = Maths.Byte(255 * bmData[2][j, i]);
                        _ix1[z] = Maths.Byte(255 * bmData[1][j, i]);
                        _ix2[z] = Maths.Byte(255 * bmData[0][j, i]);
                    }
                }
            }
            else
            {
                for (int j = 0; j < height; j++)
                {
                    for (int i = 0; i < width; i++, z++)
                    {
                        // transform
                        _ix0[z] = Maths.Byte(255 * bmData[0][j, i]);
                        _ix1[z] = Maths.Byte(255 * bmData[1][j, i]);
                        _ix2[z] = Maths.Byte(255 * bmData[2][j, i]);
                    }
                }
            }

            // arrays
            return new byte[][] { _ix0, _ix1, _ix2 };
        }
        /// <summary>
        /// Converts a Bitmap to an BGR tensor arrays.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="rgb">RGB or BGR</param>
        /// <returns>RGB tensor arrays</returns>
        public static float[][] ToFloatTensor(this Bitmap Data, bool rgb = false)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            float[][] _ix = TensorMatrix.ToFloatTensor(bmData, rgb);
            BitmapFormat.Unlock(Data, bmData);
            return _ix;
        }
        /// <summary>
        /// Converts a Bitmap to an BGR tensor arrays.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="rgb">RGB or BGR</param>
        /// <returns>RGB tensor arrays</returns>
        public unsafe static float[][] ToFloatTensor(this BitmapData bmData, bool rgb = false)
        {
            // params
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int shift = height * width;
            float[] _ix0 = new float[shift];
            float[] _ix1 = new float[shift];
            float[] _ix2 = new float[shift];
            int z = 0;

            // do job
            if (rgb)
            {
                for (int j = 0; j < height; j++)
                {
                    for (int i = 0; i < width; i++, z++)
                    {
                        int k, jstride = j * stride;
                        k = jstride + i * 4;

                        // transform
                        _ix0[z] = p[k + 2];
                        _ix1[z] = p[k + 1];
                        _ix2[z] = p[k + 0];
                    }
                }
            }
            else
            {
                for (int j = 0; j < height; j++)
                {
                    for (int i = 0; i < width; i++, z++)
                    {
                        int k, jstride = j * stride;
                        k = jstride + i * 4;

                        // transform
                        _ix0[z] = p[k + 0];
                        _ix1[z] = p[k + 1];
                        _ix2[z] = p[k + 2];
                    }
                }
            }

            // arrays
            return new float[][] { _ix0, _ix1, _ix2 };
        }
        /// <summary>
        /// Converts a Bitmap to an BGR tensor arrays.
        /// </summary>
        /// <param name="bmData">Bitmap data in BGR terms</param>
        /// <param name="rgb">RGB or BGR</param>
        /// <returns>RGB tensor arrays</returns>
        public unsafe static float[][] ToFloatTensor(this float[][,] bmData, bool rgb = false)
        {
            // params
            int width = bmData[0].GetLength(1), height = bmData[1].GetLength(0);
            int shift = height * width;
            float[] _ix0 = new float[shift];
            float[] _ix1 = new float[shift];
            float[] _ix2 = new float[shift];
            int z = 0;

            // do job
            if (rgb)
            {
                for (int j = 0; j < height; j++)
                {
                    for (int i = 0; i < width; i++, z++)
                    {
                        // transform
                        _ix0[z] = 255 * bmData[2][j, i];
                        _ix1[z] = 255 * bmData[1][j, i];
                        _ix2[z] = 255 * bmData[0][j, i];
                    }
                }
            }
            else
            {
                for (int j = 0; j < height; j++)
                {
                    for (int i = 0; i < width; i++, z++)
                    {
                        // transform
                        _ix0[z] = 255 * bmData[0][j, i];
                        _ix1[z] = 255 * bmData[1][j, i];
                        _ix2[z] = 255 * bmData[2][j, i];
                    }
                }
            }

            // arrays
            return new float[][] { _ix0, _ix1, _ix2 };
        }
        /// <summary>
        /// Converts a BGR tensor arrays to Bitmap.
        /// </summary>
        /// <param name="tensor">Tensor arrays</param>
        /// <param name="width">Width</param>
        /// <param name="height">Height</param>
        /// <param name="rgb">RGB or BGR</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromByteTensor(this byte[][] tensor, int width, int height, bool rgb = false)
        {
            // params
            Bitmap Data = new Bitmap(width, height);
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int shift = height * width;
            int z = 0;

            // do job
            if (rgb)
            {
                for (int j = 0; j < height; j++)
                {
                    for (int i = 0; i < width; i++, z++)
                    {
                        int k, jstride = j * stride;
                        k = jstride + i * 4;

                        // transform
                        p[k + 2] = tensor[0][z];
                        p[k + 1] = tensor[1][z];
                        p[k + 0] = tensor[2][z];
                    }
                }
            }
            else
            {
                for (int j = 0; j < height; j++)
                {
                    for (int i = 0; i < width; i++, z++)
                    {
                        int k, jstride = j * stride;
                        k = jstride + i * 4;

                        // transform
                        p[k + 0] = tensor[0][z];
                        p[k + 1] = tensor[1][z];
                        p[k + 2] = tensor[2][z];
                    }
                }
            }

            // arrays
            BitmapFormat.Unlock(Data, bmData);
            return Data;
        }
        /// <summary>
        /// Converts a BGR tensor arrays to Bitmap.
        /// </summary>
        /// <param name="tensor">Tensor arrays</param>
        /// <param name="width">Width</param>
        /// <param name="height">Height</param>
        /// <param name="rgb">RGB or BGR</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromFloatTensor(this float[][] tensor, int width, int height, bool rgb = false)
        {
            // params
            Bitmap Data = new Bitmap(width, height);
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int shift = height * width;
            int z = 0;

            // do job
            if (rgb)
            {
                for (int j = 0; j < height; j++)
                {
                    for (int i = 0; i < width; i++, z++)
                    {
                        int k, jstride = j * stride;
                        k = jstride + i * 4;

                        // transform
                        p[k + 2] = (byte)tensor[0][z];
                        p[k + 1] = (byte)tensor[1][z];
                        p[k + 0] = (byte)tensor[2][z];
                    }
                }
            }
            else
            {
                for (int j = 0; j < height; j++)
                {
                    for (int i = 0; i < width; i++, z++)
                    {
                        int k, jstride = j * stride;
                        k = jstride + i * 4;

                        // transform
                        p[k + 0] = (byte)tensor[0][z];
                        p[k + 1] = (byte)tensor[1][z];
                        p[k + 2] = (byte)tensor[2][z];
                    }
                }
            }

            // arrays
            BitmapFormat.Unlock(Data, bmData);
            return Data;
        }
        #endregion
    }
}
