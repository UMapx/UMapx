using UMapx.Core;
using UMapx.Decomposition;

namespace UMapx.AuditProbe;

// Run potentially nonterminating operations in a process the tests can kill.
public static class Program
{
    public static void Main(string[] args)
    {
        Console.WriteLine(args[0] switch
        {
            "IsPrimeInt" => Maths.IsPrime(int.Parse(args[1])),
            "IsPrimeLong" => Maths.IsPrime(long.Parse(args[1])),
            "FactorLong" => string.Join(",", Maths.Itf(long.Parse(args[1]))),
            "SchurZero" => SchurZero(int.Parse(args[1])),
            "EigenScale" => EigenScale(args[1]),
            _ => throw new ArgumentException("Unknown audit operation.")
        });
    }

    private static bool EigenScale(string argument)
    {
        var parts = argument.Split(',');
        float scale = float.Parse(parts[1], System.Globalization.CultureInfo.InvariantCulture);
        float eps = float.Parse(parts[2], System.Globalization.CultureInfo.InvariantCulture);
        var a = parts[0] switch
        {
            "Rotation" => new float[,] { { 0, -scale }, { scale, 0 } },
            "Symmetric" => new float[,] { { 2 * scale, scale }, { scale, 2 * scale } },
            _ => new float[,] { { scale, scale, 0 }, { 0, 2 * scale, scale }, { 0, 0, 3 * scale } }
        };
        var d = EVD.Decompose(a, eps);
        int n = a.GetLength(0);
        if (d.V.Cast<float>().Any(x => !float.IsFinite(x))) return false;
        if (d.D.Any(x => !float.IsFinite(x.Real) || !float.IsFinite(x.Imag))) return false;
        double[] expected = parts[0] == "Rotation" ? new[] { 1.0, 1.0 }
            : parts[0] == "Symmetric" ? new[] { 1.0, 3.0 } : new[] { 1.0, 2.0, 3.0 };
        var magnitudes = d.D.Select(z => System.Numerics.Complex.Abs((System.Numerics.Complex)z) / scale).OrderBy(x => x).ToArray();
        for (int j = 0; j < n; j++)
        {
            if (Math.Abs(magnitudes[j] - expected[j]) > 2e-5) return false;
            // Verify real storage of conjugate eigenvectors directly, independently of RealEigenvalueMatrix.
            double norm = 0, error = 0;
            for (int i = 0; i < n; i++)
            {
                double av = 0;
                for (int k = 0; k < n; k++) av += (a[i, k] / (double)scale) * d.V[k, j];
                double vd = d.V[i, j] * (d.D[j].Real / (double)scale);
                if (d.D[j].Imag > 0) vd -= d.V[i, j + 1] * (d.D[j].Imag / (double)scale);
                if (d.D[j].Imag < 0) vd -= d.V[i, j - 1] * (d.D[j].Imag / (double)scale);
                error += (av - vd) * (av - vd);
                norm += (double)d.V[i, j] * d.V[i, j];
            }
            if (!(norm > 0) || Math.Sqrt(error / norm) > 2e-5) return false;
        }
        return true;
    }

    private static bool SchurZero(int size)
    {
        var d=Schur.Decompose(new float[size,size],1e-7f);
        if(d.T.Cast<float>().Any(x=>x!=0))return false;
        for(int i=0;i<size;i++)for(int j=0;j<size;j++)
        {
            double dot=Enumerable.Range(0,size).Sum(k=>(double)d.Q[k,i]*d.Q[k,j]);
            if(!double.IsFinite(dot)||Math.Abs(dot-(i==j?1:0))>1e-5)return false;
        }
        return true;
    }
}
