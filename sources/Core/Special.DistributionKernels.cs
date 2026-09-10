using System;
using System.Numerics;

namespace UMapx.Core
{
    public static partial class Special
    {
        // Reuse the audited double kernels without rounding intermediate distribution values to float.
        internal static double DistributionLogGamma(double x) => GammaLog(x);
        internal static double DistributionLogBeta(double a, double b) => BetaLog(a, b);
        internal static double DistributionDigamma(double x) => Polygamma(new Complex(x, 0), false).Real;
        internal static double DistributionGamma(double a, double x, bool upper) => IncompleteGamma(a, x, upper, true);
        internal static double DistributionBeta(double a, double b, double x) => IncompleteBeta(a, b, x, true);
        internal static double DistributionErfc(double x) => ErfcValue(x);
        internal static double DistributionExpm1(double x) => Expm1(x);
    }
}
