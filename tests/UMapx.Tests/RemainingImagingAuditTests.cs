using System.Drawing;
using System.Runtime.Versioning;
using UMapx.Core;
using UMapx.Imaging;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[SupportedOSPlatform("windows")]
[Trait("Category","Imaging")]
public class RemainingImagingAuditTests
{
    [Theory] [InlineData(false)] [InlineData(true)]
    public void TensorPackingAndAveragingPreserveChannelOrder(bool sliced)
    {
        byte[][] bytes={new byte[]{1,10,40,100},new byte[]{2,20,50,200},new byte[]{3,30,60,250}};
        var floats=bytes.Select(c=>c.Select(v=>(float)v).ToArray()).ToArray();
        var expected=sliced?bytes.SelectMany(c=>c).ToArray():Enumerable.Range(0,4).SelectMany(i=>bytes.Select(c=>c[i])).ToArray();
        Assert.Equal(expected,bytes.Merge(sliced));Close(expected.Select(v=>(float)v).ToArray(),floats.Merge(sliced));
        var bmean=bytes.Average();var fmean=floats.Average();for(int i=0;i<4;i++){double mean=bytes.Sum(c=>(int)c[i])/3.0;Assert.Equal((byte)mean,bmean[i]);Close(mean,fmean[i]);}
        TensorTransform.Compute(floats,new[]{1f,2f,3f},(a,b)=>a.Select(v=>v+b).ToArray());
        TensorTransform.Compute(floats,2f,(a,b)=>a.Select(v=>v*b).ToArray());
        for(int c=0;c<3;c++)for(int i=0;i<4;i++)Close(2*(bytes[c][i]+c+1),floats[c][i]);
    }
    [Theory] [InlineData(false)] [InlineData(true)]
    public void ConstantColorTransferPreservesIdenticalImages(bool inverted)
    {
        using var source=new CanvasColor(16,12,Color.FromArgb(77,99,123)).Create();using var actual=(Bitmap)source.Clone();
        new ColorTransfer(0,inverted).Apply(actual,source);ImagingAuditTests.Same(source,actual,1,false);
    }
    [Fact]
    public void ChannelRotationHasTheExpectedPermutationAndOrderThree()
    {
        using var original=ImagingAuditTests.Pattern();using var actual=(Bitmap)original.Clone();var filter=new RotateChannel();filter.Apply(actual);
        for(int y=0;y<actual.Height;y++)for(int x=0;x<actual.Width;x++){var c=original.GetPixel(x,y);ImagingAuditTests.Pixel(Color.FromArgb(c.A,c.G,c.B,c.R),actual.GetPixel(x,y));}
        filter.Apply(actual);filter.Apply(actual);ImagingAuditTests.Same(original,actual);
    }
    public static IEnumerable<object[]> AnaglyphCases()=>Enum.GetValues<AnaglyphMode>().Select(m=>new object[]{m});
    [Theory] [MemberData(nameof(AnaglyphCases))]
    public void AnaglyphColorChannelsMatchTheDocumentedMatrices(AnaglyphMode mode)
    {
        using var left=ImagingAuditTests.Pattern();using var original=(Bitmap)left.Clone();using var right=new CanvasColor(left.Width,left.Height,Color.FromArgb(71,83,197)).Create();
        new StereoAnaglyph(mode).Apply(left,right);
        for(int y=0;y<left.Height;y++)for(int x=0;x<left.Width;x++)
        {
            var a=original.GetPixel(x,y);var b=right.GetPixel(x,y);double la=.299*a.R+.587*a.G+.114*a.B,lb=.299*b.R+.587*b.G+.114*b.B;
            (double r,double g,double bl)=mode.ToString() switch{"True"=>(la,0,lb),"Gray"=>(la,lb,lb),"Color"=>((double)a.R,b.G,b.B),"HalfColor"=>(la,b.G,b.B),"Optimized"=>(.7*a.G+.3*a.B,b.G,b.B),_=>throw new InvalidOperationException(mode.ToString())};
            ImagingAuditTests.Pixel(Color.FromArgb(a.A,(int)r,(int)g,(int)bl),left.GetPixel(x,y),1);
        }
    }
    [Theory] [InlineData(false,false)] [InlineData(false,true)] [InlineData(true,false)] [InlineData(true,true)]
    public void SelectiveGrayscaleKeepsOnlyHuesInsideTheRequestedArc(bool hsl,bool wrap)
    {
        using var actual=new Bitmap(3,1);actual.SetPixel(0,0,Color.Red);actual.SetPixel(1,0,Color.Lime);actual.SetPixel(2,0,Color.Blue);
        int min=wrap?300:60,max=wrap?60:180;IBitmapFilter filter=hsl?new HSLGrayscale(min,max):new HSBGrayscale(min,max);filter.Apply(actual);
        for(int x=0;x<3;x++){var c=actual.GetPixel(x,0);bool keep=wrap?x==0:x==1;if(keep)ImagingAuditTests.Pixel(x==0?Color.Red:Color.Lime,c,1);else{Assert.Equal(c.R,c.G);Assert.Equal(c.R,c.B);}}
    }
    public static IEnumerable<object[]> DiffusionCases()
    {
        foreach(var p in typeof(ErrorDiffusionDithering).GetProperties(System.Reflection.BindingFlags.Public|System.Reflection.BindingFlags.Static))
        foreach(int level in new[]{0,255})yield return new object[]{p.Name,level};
    }
    [Theory] [MemberData(nameof(DiffusionCases))]
    public void DiffusionKernelsPreserveBlackAndWhiteFixedPoints(string name,int intensity)
    {
        using var bitmap=new CanvasColor(17,13,Color.FromArgb(intensity,intensity,intensity)).Create();
        var filter=(ErrorDiffusionDithering)typeof(ErrorDiffusionDithering).GetProperty(name)!.GetValue(null)!;filter.Apply(bitmap);
        for(int y=0;y<bitmap.Height;y++)for(int x=0;x<bitmap.Width;x++)ImagingAuditTests.Pixel(Color.FromArgb(intensity,intensity,intensity),bitmap.GetPixel(x,y));
    }
    [Fact]
    public void FloydSteinbergDiffusionMatchesAnIndependentScalarRaster()
    {
        const int width=17,height=13;var expected=new int[height,width];using var actual=new Bitmap(width,height);
        for(int y=0;y<height;y++)for(int x=0;x<width;x++){int v=(x*31+y*17)%256;expected[y,x]=v;actual.SetPixel(x,y,Color.FromArgb(v,v,v));}
        for(int y=0;y<height;y++)for(int x=0;x<width;x++)
        {
            // The API's two-level palette changes at 129; isolate diffusion from palette design.
            int old=expected[y,x],value=old<129?0:255,error=old-value;expected[y,x]=value;
            foreach(var d in new[]{(1,0,7),(-1,1,3),(0,1,5),(1,1,1)}){int xx=x+d.Item1,yy=y+d.Item2;if(xx>=0&&xx<width&&yy<height)expected[yy,xx]=Math.Clamp((int)(expected[yy,xx]+error*d.Item3/16.0),0,255);}
        }
        ErrorDiffusionDithering.FloydSteinberg.Apply(actual);
        for(int y=0;y<height;y++)for(int x=0;x<width;x++)ImagingAuditTests.Pixel(Color.FromArgb(expected[y,x],expected[y,x],expected[y,x]),actual.GetPixel(x,y));
    }
    [Theory] [InlineData(RGBA.Red)] [InlineData(RGBA.Green)] [InlineData(RGBA.Blue)]
    public void ChromaKeysSeparateTheKeyColorFromAnotherPrimary(RGBA channel)
    {
        var key=channel==RGBA.Red?Color.Red:channel==RGBA.Green?Color.Lime:Color.Blue;var foreground=channel==RGBA.Red?Color.Blue:Color.Red;
        using var actual=new CanvasColor(32,24,key).Create();using(var g=Graphics.FromImage(actual))using(var brush=new SolidBrush(foreground))g.FillRectangle(brush,16,0,16,24);
        new Chromakey(channel).Apply(actual);for(int y=0;y<24;y++)for(int x=0;x<32;x++)ImagingAuditTests.Pixel(x<16?Color.Black:Color.White,actual.GetPixel(x,y));
    }
    [Fact]
    public void GridDisplacementPreservesUniformImages()
    {
        using var expected=new CanvasColor(32,24,Color.FromArgb(71,83,197)).Create();using var actual=(Bitmap)expected.Clone();new Grid(4,1).Apply(actual);ImagingAuditTests.Same(expected,actual);
    }
}
