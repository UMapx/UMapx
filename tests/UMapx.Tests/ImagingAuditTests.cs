using System.Runtime.Versioning;
using System.Drawing;
using System.Drawing.Imaging;
using System.Reflection;
using System.Runtime.InteropServices;
using UMapx.Core;
using UMapx.Imaging;
using Xunit;
using static UMapx.Tests.NumericAssert;
using InterpolationMode=UMapx.Core.InterpolationMode;

namespace UMapx.Tests;

[Trait("Category", "Imaging")]
[SupportedOSPlatform("windows")]
public class ImagingAuditTests
{
    internal static Bitmap Pattern(int width=9,int height=7,bool alpha=true)
    {
        var b=new Bitmap(width,height,PixelFormat.Format32bppArgb);
        for(int y=0;y<height;y++)for(int x=0;x<width;x++)b.SetPixel(x,y,Color.FromArgb(alpha?40+(17*x+23*y)%216:255,(47*x+13*y)%256,(11*x+53*y)%256,(89*x+7*y)%256));return b;
    }
    internal static void Pixel(Color expected,Color actual,int tolerance=0,bool alpha=true)
    {
        if(alpha)Assert.InRange(Math.Abs(expected.A-actual.A),0,tolerance);
        Assert.True(Math.Abs(expected.R-actual.R)<=tolerance&&Math.Abs(expected.G-actual.G)<=tolerance&&Math.Abs(expected.B-actual.B)<=tolerance,$"Expected RGBA {expected.R},{expected.G},{expected.B},{expected.A}; actual {actual.R},{actual.G},{actual.B},{actual.A}; tolerance {tolerance}.");
    }
    internal static void Same(Bitmap a,Bitmap b,int tolerance=0,bool alpha=true)
    {
        Assert.Equal(a.Size,b.Size);for(int y=0;y<a.Height;y++)for(int x=0;x<a.Width;x++)Pixel(a.GetPixel(x,y),b.GetPixel(x,y),tolerance,alpha);
    }
    static IBitmapFilter Neutral(string name,Space space=Space.RGB)=>name switch
    {
        "Brightness"=>new BrightnessCorrection(0,space),"Contrast"=>new ContrastCorrection(0,space),"Gamma"=>new GammaCorrection(1,space),
        "Linear"=>new LinearCorrection(0,space),"Levels"=>new LevelsCorrection(new RangeFloat(0,1),new RangeFloat(0,1),space),"ShiftCorrection"=>new ShiftCorrection(0,space),
        "RGB"=>new RGBFilter(0,0,0),"HSB"=>new HSBFilter(0,0,0),"HSL"=>new HSLFilter(0,0,0),"CMYK"=>new CMYKFilter(0,0,0,0),"YCbCr"=>new YCbCrFilter(0,0,0),
        "Saturation"=>new SaturationCorrection(0),"Vibrance"=>new VibranceCorrection(0),"Transparency"=>new TransparencyCorrection(0),
        "Convolution"=>new Convolution(new float[,]{{1}}),"Flip"=>new Flip(false,false),"Shift"=>new Shift(0,0),"BoxBlur"=>new BoxBlur(1),
        "Erosion"=>new Erosion(1),"Dilatation"=>new Dilatation(1),"Median"=>new Median(1),"Opening"=>new Opening(1),"Closing"=>new Closing(1),
        _=>throw new ArgumentException(name)
    };
    public static IEnumerable<object[]> NeutralCases()
    {
        foreach(var name in new[]{"Brightness","Contrast","Gamma","Linear","Levels","ShiftCorrection"})foreach(var space in Enum.GetValues<Space>().Where(s=>s!=Space.Grayscale))yield return new object[]{name,space};
        foreach(var name in new[]{"RGB","HSB","HSL","CMYK","YCbCr","Saturation","Vibrance","Transparency","Convolution","Flip","Shift","BoxBlur","Erosion","Dilatation","Median","Opening","Closing"})yield return new object[]{name,Space.RGB};
    }
    [Theory] [MemberData(nameof(NeutralCases))]
    public void NeutralFilterSettingsPreservePixels(string name,Space space)
    {
        using var source=Pattern();using var actual=(Bitmap)source.Clone();Neutral(name,space).Apply(actual);
        // HSB/HSL have whole-degree hue storage, allowing at most five byte levels of quantization.
        int tolerance=name is "HSB" or "HSL"||space is Space.HSB or Space.HSL?5:name is "CMYK" or "YCbCr"||space==Space.YCbCr?2:0;
        Same(source,actual,tolerance);
    }

    public static IEnumerable<object[]> ConversionCases()
    {
        foreach(string space in new[]{"RGB","HSB","HSL","YCbCr"})foreach(bool alpha in new[]{false,true})foreach(bool locked in new[]{false,true})yield return new object[]{space,alpha,locked};
    }
    [Theory] [MemberData(nameof(ConversionCases))]
    public void BitmapChannelConversionsPreserveShapeAlphaAndColor(string space,bool alpha,bool locked)
    {
        using var source=Pattern();float[][,] channels;
        if(locked){var data=source.Lock32bpp();try{channels=(float[][,])typeof(BitmapMatrix).GetMethod("To"+space,new[]{typeof(BitmapData),typeof(bool)})!.Invoke(null,new object[]{data,alpha})!;}finally{source.Unlock(data);}}
        else channels=(float[][,])typeof(BitmapMatrix).GetMethod("To"+space,new[]{typeof(Bitmap),typeof(bool)})!.Invoke(null,new object[]{source,alpha})!;
        Assert.Equal(alpha?4:3,channels.Length);foreach(var c in channels){Assert.Equal(source.Height,c.GetLength(0));Assert.Equal(source.Width,c.GetLength(1));Assert.All(c.Cast<float>(),v=>Assert.True(float.IsFinite(v)));}
        if(space=="RGB")for(int y=0;y<source.Height;y++)for(int x=0;x<source.Width;x++){var c=source.GetPixel(x,y);Close(c.B/255.0,channels[0][y,x]);Close(c.G/255.0,channels[1][y,x]);Close(c.R/255.0,channels[2][y,x]);}
        // BitmapMatrix deliberately stores B,G,R channel planes; TensorMatrix has an explicit RGB switch.
        using var restored=(Bitmap)typeof(BitmapMatrix).GetMethod("From"+space,new[]{typeof(float[][,])})!.Invoke(null,new object[]{channels})!;
        Same(source,restored,space is "HSL" or "HSB"?5:1,alpha);
        using var destination=new Bitmap(source.Width,source.Height,PixelFormat.Format32bppArgb);
        typeof(BitmapMatrix).GetMethod("From"+space,new[]{typeof(float[][,]),typeof(Bitmap)})!.Invoke(null,new object[]{channels,destination});Same(restored,destination,0);
    }

    [Theory] [InlineData(false,false)] [InlineData(true,false)] [InlineData(false,true)] [InlineData(true,true)]
    public void TensorConversionsMatchPixelChannelOrder(bool rgb,bool floating)
    {
        using var source=Pattern(alpha:false);var planes=source.ToRGB();
        if(floating)
        {
            var tensor=source.ToFloatTensor(rgb);var fromPlanes=planes.ToFloatTensor(rgb);
            Assert.Equal(3,tensor.Length);for(int c=0;c<3;c++)for(int y=0;y<source.Height;y++)for(int x=0;x<source.Width;x++)
            {var p=source.GetPixel(x,y);double v=(rgb?new[]{p.R,p.G,p.B}:new[]{p.B,p.G,p.R})[c];Close(v,tensor[c][y*source.Width+x]);Close(v,fromPlanes[c][y*source.Width+x]);}
            using var restored=tensor.FromFloatTensor(source.Width,source.Height,rgb);Same(source,restored,1);
        }
        else
        {
            var tensor=source.ToByteTensor(rgb);var fromPlanes=planes.ToByteTensor(rgb);Assert.Equal(3,tensor.Length);
            for(int c=0;c<3;c++)for(int y=0;y<source.Height;y++)for(int x=0;x<source.Width;x++){var p=source.GetPixel(x,y);var v=(rgb?new[]{p.R,p.G,p.B}:new[]{p.B,p.G,p.R})[c];Assert.Equal(v,tensor[c][y*source.Width+x]);Assert.Equal(v,fromPlanes[c][y*source.Width+x]);}
            using var restored=tensor.FromByteTensor(source.Width,source.Height,rgb);Same(source,restored);
        }
    }

    [Theory] [InlineData("Invert")] [InlineData("RGB")] [InlineData("Gamma")] [InlineData("Brightness")] [InlineData("Grayscale")]
    public void PixelwiseFiltersMatchIndependentChannelEquations(string name)
    {
        using var source=Pattern();using var actual=(Bitmap)source.Clone();IBitmapFilter filter=name switch{"Invert"=>new InvertChannels(Space.RGB),"RGB"=>new RGBFilter(17,-23,31),"Gamma"=>new GammaCorrection(2,Space.RGB),"Brightness"=>new BrightnessCorrection(.2f,Space.RGB),_=>new Grayscale(.25f,.5f,.25f)};
        filter.Apply(actual);
        for(int y=0;y<source.Height;y++)for(int x=0;x<source.Width;x++)
        {
            var p=source.GetPixel(x,y);int Byte(double v)=>Math.Clamp((int)v,0,255);int Channel(int v,int offset)=>name switch{"Invert"=>255-v,"RGB"=>Byte(v+offset),"Gamma"=>Byte(255*Math.Pow(v/255.0,2)),"Brightness"=>Byte(v+25.5),_=>Byte(.25*p.R+.5*p.G+.25*p.B)};
            Pixel(Color.FromArgb(p.A,Channel(p.R,17),Channel(p.G,-23),Channel(p.B,31)),actual.GetPixel(x,y),name is "Gamma" or "Brightness"?1:0);
        }
    }

    [Theory] [InlineData("Erosion")] [InlineData("Dilatation")] [InlineData("Median")]
    public void BitmapMorphologyMatchesIndependentNeighborhoodSorting(string name)
    {
        using var source=Pattern();using var actual=(Bitmap)source.Clone();IBitmapFilter filter=name switch{"Erosion"=>new Erosion(3,5),"Dilatation"=>new Dilatation(3,5),_=>new Median(3,5)};filter.Apply(actual);
        for(int y=0;y<source.Height;y++)for(int x=0;x<source.Width;x++)
        {
            int Channel(Func<Color,int> component){var v=new List<int>();for(int dy=-2;dy<=2;dy++)for(int dx=-1;dx<=1;dx++)v.Add(component(source.GetPixel(Math.Clamp(x+dx,0,source.Width-1),Math.Clamp(y+dy,0,source.Height-1))));v.Sort();return name=="Erosion"?v[0]:name=="Dilatation"?v[^1]:v[v.Count/2];}
            Pixel(Color.FromArgb(source.GetPixel(x,y).A,Channel(c=>c.R),Channel(c=>c.G),Channel(c=>c.B)),actual.GetPixel(x,y));
        }
    }

    [Theory] [InlineData(InterpolationMode.NearestNeighbor)] [InlineData(InterpolationMode.Bilinear)] [InlineData(InterpolationMode.Bicubic)]
    public void BitmapResizingToTheSameDimensionsIsAnIdentity(InterpolationMode mode)
    {
        using var source=Pattern();using var actual=new Bitmap(source.Width,source.Height,PixelFormat.Format32bppArgb);new Resize(source.Width,source.Height,mode).Apply(actual,source);Same(source,actual);
    }

    [Theory] [InlineData(false,false)] [InlineData(true,false)] [InlineData(false,true)] [InlineData(true,true)]
    public void BitmapFlipsMatchPixelCoordinates(bool horizontal,bool vertical)
    {
        using var source=Pattern();using var actual=(Bitmap)source.Clone();new Flip(horizontal,vertical).Apply(actual);
        for(int y=0;y<source.Height;y++)for(int x=0;x<source.Width;x++)Pixel(source.GetPixel(horizontal?source.Width-1-x:x,vertical?source.Height-1-y:y),actual.GetPixel(x,y));
    }

    [Fact]
    public void BitmapCroppingCopiesTheRequestedRectangle()
    {
        using var source=Pattern();using var actual=new Bitmap(4,3,PixelFormat.Format32bppArgb);new Crop(2,1,4,3).Apply(actual,source);
        for(int y=0;y<3;y++)for(int x=0;x<4;x++)Pixel(source.GetPixel(x+2,y+1),actual.GetPixel(x,y));
        // Graphics.DrawImage uses premultiplied alpha and can round color channels.
        using var transformed=BitmapTransform.Crop(source,new Rectangle(2,1,4,3));Same(actual,transformed,2);
    }

    [Theory] [InlineData(RGBA.Red)] [InlineData(RGBA.Green)] [InlineData(RGBA.Blue)] [InlineData(RGBA.Alpha)]
    public void HistogramsAndHistogramStatisticsMatchPixelCounts(RGBA channel)
    {
        using var source=Pattern();var samples=new List<int>();var histogram=new int[256];
        for(int y=0;y<source.Height;y++)for(int x=0;x<source.Width;x++){var p=source.GetPixel(x,y);int v=channel switch{RGBA.Red=>p.R,RGBA.Green=>p.G,RGBA.Blue=>p.B,_=>p.A};histogram[v]++;samples.Add(v);}
        Assert.Equal(histogram,source.Histogram(channel));Assert.Equal(histogram.Sum(),Statistics.Sum(histogram));double mean=samples.Average();Close(mean,Statistics.Mean(histogram));Close(Math.Sqrt(samples.Average(v=>Math.Pow(v-mean,2))),Statistics.StdDev(histogram));
        Close(-histogram.Where(v=>v>0).Sum(v=>v/(double)samples.Count*Math.Log2(v/(double)samples.Count)),Statistics.Entropy(histogram));
        var cdf=Statistics.CDF(histogram);int sum=0;for(int i=0;i<256;i++){sum+=histogram[i];Assert.Equal(sum,cdf[i]);}
    }

    [Theory] [InlineData(1)] [InlineData(3)] [InlineData(5)]
    public void HistogramMedianSelectsTheMiddleObservationForOddCounts(int count)
    {
        var histogram=new int[256];for(int i=0;i<count;i++)histogram[20+i*30]++;Assert.Equal(20+(count/2)*30,Statistics.Median(histogram));
    }

    public static IEnumerable<object[]> StrideCases()
    {
        foreach(string name in new[]{"Grayscale","Transparency","RGB","Saturation","Vibrance","Gamma","Invert","HSB","HSL","CMYK","YCbCr","Diffusion"})foreach(bool negative in new[]{false,true})yield return new object[]{name,negative};
    }
    [Theory] [MemberData(nameof(StrideCases))]
    public void PixelFiltersHonorBitmapDataStrideAndLeavePaddingUntouched(string name,bool negative)
    {
        // Managed guard areas keep every malformed linear traversal inside allocated memory.
        const int width=3,height=4,stride=20;var bytes=Enumerable.Repeat((byte)173,3*height*stride).ToArray();int origin=height*stride+(negative?(height-1)*stride:0),step=negative?-stride:stride;
        var active=new HashSet<int>();using var reference=new Bitmap(width,height,PixelFormat.Format32bppArgb);
        for(int y=0;y<height;y++)for(int x=0;x<width;x++){int k=origin+y*step+4*x;var c=Color.FromArgb(60+x*20,30+x*30,90+y*20,180-x*20);reference.SetPixel(x,y,c);var pixel=new[]{c.B,c.G,c.R,c.A};for(int ch=0;ch<4;ch++){bytes[k+ch]=pixel[ch];active.Add(k+ch);}}
        IBitmapFilter Filter()=>name switch{"Grayscale"=>new Grayscale(.25f,.5f,.25f),"Gamma"=>new GammaCorrection(2,Space.RGB),"Invert"=>new InvertChannels(Space.RGB),"Transparency"=>new TransparencyCorrection(.2f),"RGB"=>new RGBFilter(12,23,-15),"Saturation"=>new SaturationCorrection(40),"Vibrance"=>new VibranceCorrection(40),"HSB"=>new HSBFilter(20,.1f,.1f),"HSL"=>new HSLFilter(20,.1f,.1f),"CMYK"=>new CMYKFilter(.1f,.1f,.1f,.1f),_=>new YCbCrFilter(.1f,.02f,.03f)};
        IBitmapFilter SelectedFilter()=>name=="Diffusion"?ErrorDiffusionDithering.FloydSteinberg:Filter();
        SelectedFilter().Apply(reference);var handle=GCHandle.Alloc(bytes,GCHandleType.Pinned);
        try{var data=new BitmapData{Width=width,Height=height,Stride=step,PixelFormat=PixelFormat.Format32bppArgb,Scan0=IntPtr.Add(handle.AddrOfPinnedObject(),origin)};SelectedFilter().Apply(data);}finally{handle.Free();}
        for(int y=0;y<height;y++)for(int x=0;x<width;x++){int k=origin+y*step+4*x;Pixel(reference.GetPixel(x,y),Color.FromArgb(bytes[k+3],bytes[k+2],bytes[k+1],bytes[k]));}
        for(int i=0;i<bytes.Length;i++)if(!active.Contains(i))Assert.True(bytes[i]==173,$"Filter wrote into padding or guard byte {i}.");
    }

    [Fact]
    public void DepthFloatConversionPreservesSixteenBitSamples()
    {
        ushort[,] depth={{0,1,255,256},{1023,4096,32768,65535}};var normalized=depth.ToFloat();var restored=DepthMatrix.FromFloat(normalized);
        for(int y=0;y<2;y++)for(int x=0;x<4;x++){Close(depth[y,x]/65535.0,normalized[y,x]);Assert.Equal(depth[y,x],restored[y,x]);}
    }
}
