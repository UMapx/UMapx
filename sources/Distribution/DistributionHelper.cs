using System;

namespace UMapx.Distribution
{
    internal static class DistributionHelper
    {
        public static double[] ToDouble(float[] values)
        {
            if (values == null)
            {
                return null;
            }

            double[] result = new double[values.Length];
            for (int i = 0; i < values.Length; i++)
            {
                result[i] = values[i];
            }
            return result;
        }

        public static float[] ToFloat(double[] values)
        {
            if (values == null)
            {
                return null;
            }

            float[] result = new float[values.Length];
            for (int i = 0; i < values.Length; i++)
            {
                result[i] = (float)values[i];
            }
            return result;
        }
    }
}
