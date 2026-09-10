using UMapx.Core;
using UMapx.Decomposition;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Decomposition")]
public class DecompositionAuditTests
{
    private static float[,] Identity(int n)=>Diagonal(Enumerable.Repeat(1f,n).ToArray());
    private static float[,] Input(int n,string kind)
    {
        if(kind=="Zero")return new float[n,n];
        if(kind=="Identity")return Identity(n);
        if(kind=="Diagonal")return Diagonal(Enumerable.Range(1,n).Select(i=>(float)i).ToArray());
        var a=Matrix(n,n);
        if(kind=="Positive")
        {
            a=Product(Transpose(a),a);for(int i=0;i<n;i++)a[i,i]+=1;
        }
        else if(kind=="RankOne")for(int i=0;i<n;i++)for(int j=0;j<n;j++)a[i,j]=(i+1)*(j+1);
        return a;
    }

    private static void Orthonormal(float[,] q,bool columns=true)
    {
        var qt=Transpose(q);
        Close(Identity(q.GetLength(columns?1:0)),columns?Product(qt,q):Product(q,qt),.003f);
    }
    private static void Band(float[,] a,int lower,int upper)
    {
        for(int i=0;i<a.GetLength(0);i++)for(int j=0;j<a.GetLength(1);j++)
            if(i-j>lower||j-i>upper)Close(0,a[i,j],.001);
    }

    public static IEnumerable<object[]> Factorizations()
    {
        foreach(string name in new[]{"QR","QL","LQ","RQ","LU","LDU","Diagonal","Polar","Hessenberg","Schur","EVD","Arnoldi","GramSchmidt","Cholesky","LDL","UDL","Householder","Lanczos","LanczosFull"})
            foreach(int n in new[]{2,3,5})
                foreach(string kind in new[]{"Positive","Diagonal"})yield return new object[]{name,n,kind};
        foreach(string name in new[]{"QR","QL","LQ","RQ","LU","LDU","Diagonal","Polar","Hessenberg","Schur","EVD","Arnoldi","GramSchmidt"})
            foreach(int n in new[]{2,3,5})yield return new object[]{name,n,"General"};
        foreach(string name in new[]{"QR","Polar","Hessenberg","Schur","EVD"})
            foreach(string kind in new[]{"Zero","Identity","RankOne"})
                if(name!="Schur"||kind!="Zero")yield return new object[]{name,3,kind};
    }

    [Theory] [InlineData(2)] [InlineData(3)] [InlineData(5)]
    public async Task SchurDecompositionOfTheZeroMatrixTerminates(int n)
    {
        Assert.Equal("True",await AuditProcess.RunAsync("SchurZero",n.ToString()));
    }

    [Theory, MemberData(nameof(Factorizations))]
    public void FactorsReconstructTheInputAndSatisfyTheirStructuralContracts(string name,int n,string kind)
    {
        var a=Input(n,kind); var original=(float[,])a.Clone(); float[,] actual;
        switch(name)
        {
            case "QR": {var d=new QR(a);actual=Product(d.Q,d.R);Orthonormal(d.Q);Band(d.R,0,n);Assert.Equal(a.Length,d.H.Length);break;}
            case "QL": {var d=new QL(a);actual=Product(d.Q,d.L);Orthonormal(d.Q);Band(d.L,n,0);break;}
            case "LQ": {var d=new LQ(a);actual=Product(d.L,d.Q);Orthonormal(d.Q);Band(d.L,n,0);break;}
            case "RQ": {var d=new RQ(a);actual=Product(d.R,d.Q);Orthonormal(d.Q);Band(d.R,0,n);break;}
            case "LU": {var d=new LU(a);actual=Product(d.L,d.U);Band(d.L,n,0);Band(d.U,0,n);break;}
            case "LDU": {var d=new LDU(a);actual=Product(Product(d.L,Diagonal(d.D)),d.U);Band(d.L,n,0);Band(d.U,0,n);break;}
            case "Diagonal": {var d=new UMapx.Decomposition.Diagonal(a);actual=Product(d.B,Diagonal(d.D));break;}
            case "Polar": {var d=new Polar(a,100);actual=Product(d.U,d.P);Close(d.P,Transpose(d.P));break;}
            case "Hessenberg": {var d=new Hessenberg(a);actual=Product(Product(d.P,d.H),Transpose(d.P));Orthonormal(d.P);Band(d.H,1,n);break;}
            case "Schur": {var d=new Schur(a,1e-7f);actual=Product(Product(d.Q,d.T),Transpose(d.Q));Orthonormal(d.Q);Band(d.T,1,n);break;}
            case "EVD": {var d=new EVD(a,1e-7f);Close(Product(a,d.V),Product(d.V,d.R),.003f);Close(Enumerable.Range(0,n).Sum(i=>a[i,i]),d.D.Sum(z=>z.Real),.003);return;}
            case "Arnoldi": {var d=new Arnoldi(a);actual=Product(Product(d.Q,d.H),Transpose(d.Q));Orthonormal(d.Q);Band(d.H,1,n);break;}
            case "GramSchmidt": {var d=new GramSchmidt(a);Orthonormal(d.Q);actual=Product(d.Q,Product(Transpose(d.Q),a));Band(Product(Transpose(d.Q),a),0,n);break;}
            case "Cholesky": {var d=new Cholesky(a);actual=Product(d.L,d.U);Close(d.U,Transpose(d.L));Band(d.L,n,0);break;}
            case "LDL": {var d=new LDL(a);actual=Product(Product(d.L,Diagonal(d.D)),d.U);Close(d.U,Transpose(d.L));Band(d.L,n,0);break;}
            case "UDL": {var d=new UDL(a);actual=Product(Product(d.U,Diagonal(d.D)),d.L);Close(d.U,Transpose(d.L));Band(d.L,n,0);break;}
            case "Householder": {var d=new Householder(a);actual=Product(Product(d.H,d.T),Transpose(d.H));Orthonormal(d.H);Band(d.T,1,1);break;}
            default: {var d=new Lanczos(a,name=="LanczosFull");actual=Product(Product(d.Q,d.T),Transpose(d.Q));Orthonormal(d.Q);Band(d.T,1,1);break;}
        }
        Close(original,actual,.004f); Close(original,a);
    }

    public static IEnumerable<object[]> RectangularCases()
    {
        foreach(string name in new[]{"SVD","Bidiagonal"})foreach(var shape in new[]{(1,1),(1,4),(4,1),(3,5),(5,3),(4,4)})
            foreach(string kind in new[]{"General","Zero","RankOne"})yield return new object[]{name,shape.Item1,shape.Item2,kind};
    }

    [Theory, MemberData(nameof(RectangularCases))]
    public void RectangularAndRankDeficientMatricesRetainTheirInformation(string name,int m,int n,string kind)
    {
        var a=Matrix(m,n);
        if(kind=="Zero")a=new float[m,n];
        if(kind=="RankOne")for(int i=0;i<m;i++)for(int j=0;j<n;j++)a[i,j]=(i+1)*(j+1);
        if(name=="Bidiagonal")
        {
            var d=new Bidiagonal(a);Close(a,Product(Product(d.U,d.B),Transpose(d.V)),.002f);
            Orthonormal(d.U);Orthonormal(d.V);Band(d.B,0,1);
        }
        else
        {
            var d=new SVD(a,100);Close(a,Product(Product(d.U,Diagonal(d.S)),Transpose(d.V)),.003f);
            Assert.All(d.S,s=>Assert.True(float.IsFinite(s)&&s>=0));Orthonormal(d.U);Orthonormal(d.V);
            // The pseudoinverse is independently characterized by the four Penrose equations.
            var p=d.P;var ap=Product(a,p);var pa=Product(p,a);
            Close(a,Product(ap,a),.005f);Close(p,Product(pa,p),.005f);
            Close(ap,Transpose(ap),.005f);Close(pa,Transpose(pa),.005f);
        }
    }

    [Theory] [InlineData(2)] [InlineData(3)] [InlineData(5)]
    public void GeneralizedFactorizationsSatisfyBothMatrixEquations(int n)
    {
        var a=Input(n,"General");var b=Input(n,"Positive");
        var qz=new QZ(a,b,1e-7f);Orthonormal(qz.Q);Orthonormal(qz.Z);
        Close(a,Product(Product(qz.Q,qz.S),Transpose(qz.Z)),.002f);
        Close(b,Product(Product(qz.Q,qz.T),Transpose(qz.Z)),.002f);Band(qz.T,0,n);Band(qz.S,1,n);
        var e=new GEVD(a,b,1e-7f);Assert.False(e.IsSingular);
        Assert.True(e.V.Cast<float>().All(float.IsFinite),"GEVD returned nonfinite eigenvectors for a finite nonsingular matrix pair.");
        Close(Product(a,e.V),Product(Product(b,e.V),e.D),.003f);
        for(int i=0;i<n;i++)Close((System.Numerics.Complex)e.Alpha[i]/e.Beta[i],e.Eigenvalues[i],.001);
    }

    [Theory] [InlineData(2,2,2)] [InlineData(5,4,3)] [InlineData(7,6,4)]
    public void GeneralizedSvdReconstructsBothRectangularInputs(int m,int p,int n)
    {
        var a=Matrix(m,n);var b=Matrix(p,n);for(int i=0;i<n;i++)b[i,i]+=1;
        var d=new GSVD(a,b,100);
        Close(a,Product(Product(d.U1,Diagonal(d.S1)),d.X),.003f);
        Close(b,Product(Product(d.U2,Diagonal(d.S2)),d.X),.003f);
        Orthonormal(d.U1);Orthonormal(d.U2);
        for(int i=0;i<n;i++){Close(1,d.S1[i]*d.S1[i]+d.S2[i]*d.S2[i],.001);Close(1,d.Identity[i],.001);Close(d.S1[i]/d.S2[i],d.Gamma[i],.001);}
    }

    [Theory] [InlineData(2)] [InlineData(3)] [InlineData(5)]
    public void HouseholderReflectionZerosTheTailAndIsOrthogonal(int n)
    {
        var v=Enumerable.Range(0,n).Select(i=>i+.5f).ToArray();var d=new Householder(v);
        Orthonormal(d.H);Close(d.H,Transpose(d.H));
        var transformed=new float[n];for(int i=0;i<n;i++)for(int j=0;j<n;j++)transformed[i]+=d.H[i,j]*v[j];
        for(int i=1;i<n;i++)Close(0,transformed[i],1e-5);
        Close(Math.Sqrt(v.Sum(t=>(double)t*t)),Math.Abs(transformed[0]),1e-5);
    }

    [Fact]
    public void PowerIterationFindsTheDominantEigenvector()
    {
        var a=new float[,]{{4,1,0},{0,2,0},{0,0,1}};var d=new Power(a,100);var v=d.V;
        Close(1,Math.Sqrt(v.Sum(x=>(double)x*x)));Close(0,v[1],1e-5);Close(0,v[2],1e-5);Close(Diagonal(v),d.J);
    }

    [Theory] [InlineData(3,5)] [InlineData(5,3)]
    public void NonnegativeRankOneFactorizationRecoversANonnegativeProduct(int m,int n)
    {
        var a=new float[m,n];for(int i=0;i<m;i++)for(int j=0;j<n;j++)a[i,j]=(i+1)*(j+1);
        var d=new NMF(a,1,100);Assert.All(d.W.Cast<float>(),x=>Assert.True(float.IsFinite(x)&&x>=0));Assert.All(d.H.Cast<float>(),x=>Assert.True(float.IsFinite(x)&&x>=0));
        Close(a,Product(d.W,d.H),.001f);
    }
}
