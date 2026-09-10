using System;

namespace UMapx.Core
{
    public static partial class Matrice
    {
        private static double HermitianSquaredDifference(Complex32[] x, Complex32[] y)
        {
            if (x.Length != y.Length) throw new ArgumentException("Vector lengths must match.", nameof(y));
            if (x.Length < 2) return double.NaN;
            double sum = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double re = (double)x[i].Real - y[i].Real, im = (double)x[i].Imag - y[i].Imag;
                sum += re * re + im * im;
            }
            return sum / (x.Length - 1);
        }

        private static double HermitianVariance(Complex32[] values)
        {
            if (values.Length < 2) return double.NaN;
            double meanReal = 0, meanImag = 0;
            foreach (var value in values) { meanReal += value.Real; meanImag += value.Imag; }
            meanReal /= values.Length; meanImag /= values.Length;
            double sum = 0;
            foreach (var value in values)
            {
                double re = value.Real - meanReal, im = value.Imag - meanImag;
                sum += re * re + im * im;
            }
            return sum / (values.Length - 1);
        }
    }
}
