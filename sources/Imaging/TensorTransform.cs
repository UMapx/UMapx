namespace UMapx.Imaging
{
    /// <summary>
    /// Uses for editing and transforming tensors.
    /// </summary>
    public static class TensorTransform
    {
        #region Merge
        /// <summary>
        /// Merges image tensors to single tensor.
        /// </summary>
        /// <param name="image">RGB tensor arrays</param>
        /// <param name="slice">Slice or not</param>
        /// <returns>Byte array</returns>
        public static byte[] Merge(this byte[][] image, bool slice = false)
        {
            int count = image.Length;
            int length = image[0].GetLength(0);
            byte[] _ix = new byte[count * length];
            int z = 0;

            if (slice)
            {
                for (int i = 0; i < count; i++)
                {
                    for (int j = 0; j < length; j++)
                    {
                        _ix[z++] = image[i][j];
                    }
                }
            }
            else
            {
                for (int j = 0; j < length; j++)
                {
                    for (int i = 0; i < count; i++)
                    {
                        _ix[z++] = image[i][j];
                    }
                }
            }

            return _ix;
        }
        /// <summary>
        /// Merges image tensors to single tensor.
        /// </summary>
        /// <param name="image">RGB tensor arrays</param>
        /// <param name="slice">Slice or not</param>
        /// <returns>Float array</returns>
        public static float[] Merge(this float[][] image, bool slice = false)
        {
            int count = image.Length;
            int length = image[0].GetLength(0);
            float[] _ix = new float[count * length];
            int z = 0;

            if (slice)
            {
                for (int i = 0; i < count; i++)
                {
                    for (int j = 0; j < length; j++)
                    {
                        _ix[z++] = image[i][j];
                    }
                }
            }
            else
            {
                for (int j = 0; j < length; j++)
                {
                    for (int i = 0; i < count; i++)
                    {
                        _ix[z++] = image[i][j];
                    }
                }
            }

            return _ix;
        }
        #endregion

        #region Average
        /// <summary>
        /// Averages image tensors to single tensor.
        /// </summary>
        /// <param name="image">RGB tensor arrays</param>
        /// <returns>Byte array</returns>
        public static byte[] Average(this byte[][] image)
        {
            int count = image.Length;
            int length = image[0].GetLength(0);
            byte[] _ix = new byte[length];
            int z = 0;

            for (int j = 0; j < length; j++)
            {
                var value = 0;

                for (int i = 0; i < count; i++)
                {
                    value += image[i][j];
                }

                _ix[z++] = (byte)(value / count);
            }

            return _ix;
        }
        /// <summary>
        /// Averages image tensors to single tensor.
        /// </summary>
        /// <param name="image">RGB tensor arrays</param>
        /// <returns>Byte array</returns>
        public static float[] Average(this float[][] image)
        {
            int count = image.Length;
            int length = image[0].GetLength(0);
            float[] _ix = new float[length];
            int z = 0;

            for (int j = 0; j < length; j++)
            {
                var value = 0.0f;

                for (int i = 0; i < count; i++)
                {
                    value += image[i][j];
                }

                _ix[z++] = (float)(value / count);
            }

            return _ix;
        }
        #endregion

        #region Operator
        /// <summary>
        /// Implements operator function.
        /// </summary>
        /// <param name="image">RGB tensor arrays</param>
        /// <param name="b">Vector</param>
        /// <param name="tensorOperator">Operator</param>
        public static void Compute(this float[][] image, float[] b, ITensorOperator tensorOperator)
        {
            int count = image.Length;

            for (int i = 0; i < count; i++)
            {
                image[i] = tensorOperator(image[i], b[i]);
            }

            return;
        }
        /// <summary>
        /// Implements operator function.
        /// </summary>
        /// <param name="image">RGB tensor arrays</param>
        /// <param name="b">Value</param>
        /// <param name="tensorOperator">Operator</param>
        public static void Compute(this float[][] image, float b, ITensorOperator tensorOperator)
        {
            int count = image.Length;

            for (int i = 0; i < count; i++)
            {
                image[i] = tensorOperator(image[i], b);
            }

            return;
        }
        #endregion

        #region Delegate
        /// <summary>
        /// Tensor operator.
        /// </summary>
        /// <param name="a">Vector</param>
        /// <param name="b">Value</param>
        /// <returns></returns>
        public delegate float[] ITensorOperator(float[] a, float b);
        #endregion
    }
}
