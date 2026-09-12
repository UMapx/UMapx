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
            _ => throw new ArgumentException("Unknown audit operation.")
        });
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
