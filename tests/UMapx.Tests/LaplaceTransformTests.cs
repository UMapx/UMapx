using System;
using UMapx.Core;
using UMapx.Transform;
using Xunit;

namespace UMapx.Tests
{
    public class LaplaceTransformTests
    {

        [Theory]
        [InlineData(0.0)]
        [InlineData(0.35)]
        public void LaplaceTransformMatchesExpected(double sigma)
        {
            float sigmaF = (float)sigma;
            Complex32[] signal = CreateUnitStep(8);
            Complex32[] original = (Complex32[])signal.Clone();
            Complex32[] expected = ComputeExpected(original, sigmaF, normalized: true);
            var transform = new LaplaceTransform(sigmaF, normalized: true);

            Complex32[] actual = transform.Forward(signal);

            AssertComplexArrayEqual(expected, actual, 1e-4f);

            Complex32[] reconstructed = transform.Backward(actual);
            AssertComplexArrayEqual(original, reconstructed, 1e-4f);
        }

        [Theory]
        [InlineData(0.0)]
        [InlineData(0.35)]
        public void FastLaplaceTransformMatchesExpected(double sigma)
        {
            float sigmaF = (float)sigma;
            Complex32[] signal = CreateUnitStep(8);
            Complex32[] original = (Complex32[])signal.Clone();
            Complex32[] expected = ComputeExpected(original, sigmaF, normalized: true);
            var transform = new FastLaplaceTransform(sigmaF, normalized: true);

            Complex32[] actual = transform.Forward(signal);

            AssertComplexArrayEqual(expected, actual, 1e-4f);

            if (sigmaF == 0f)
            {
                var fft = new FastFourierTransform(normalized: true);
                Complex32[] fftResult = fft.Forward(original);
                AssertComplexArrayEqual(fftResult, actual, 1e-4f);
            }

            Complex32[] reconstructed = transform.Backward(actual);
            AssertComplexArrayEqual(original, reconstructed, 1e-4f);
        }

        private static Complex32[] CreateUnitStep(int length)
        {
            var result = new Complex32[length];
            for (int i = 0; i < length; i++)
            {
                result[i] = new Complex32(1f, 0f);
            }
            return result;
        }

        private static Complex32[] ComputeExpected(Complex32[] signal, float sigma, bool normalized)
        {
            int length = signal.Length;
            var result = new Complex32[length];
            double scale = normalized ? Math.Sqrt(length) : 1.0;

            for (int k = 0; k < length; k++)
            {
                double sumRe = 0.0;
                double sumIm = 0.0;

                for (int n = 0; n < length; n++)
                {
                    double weight = Math.Exp(-sigma * n);
                    double re = signal[n].Real * weight;
                    double im = signal[n].Imag * weight;
                    double angle = -2.0 * Math.PI * k * n / length;
                    double cos = Math.Cos(angle);
                    double sin = Math.Sin(angle);

                    sumRe += re * cos - im * sin;
                    sumIm += re * sin + im * cos;
                }

                result[k] = new Complex32((float)(sumRe / scale), (float)(sumIm / scale));
            }

            return result;
        }

        private static void AssertComplexArrayEqual(Complex32[] expected, Complex32[] actual, float tolerance)
        {
            Assert.Equal(expected.Length, actual.Length);
            for (int i = 0; i < expected.Length; i++)
            {
                Assert.True(Math.Abs(expected[i].Real - actual[i].Real) <= tolerance,
                    $"Real parts differ at index {i}: expected {expected[i].Real}, actual {actual[i].Real}");
                Assert.True(Math.Abs(expected[i].Imag - actual[i].Imag) <= tolerance,
                    $"Imaginary parts differ at index {i}: expected {expected[i].Imag}, actual {actual[i].Imag}");
            }
        }
    }
}
