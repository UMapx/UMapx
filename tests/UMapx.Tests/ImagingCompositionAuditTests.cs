using System.Drawing;
using System.Runtime.Versioning;
using UMapx.Core;
using UMapx.Imaging;
using Xunit;
using static UMapx.Tests.NumericAssert;
using InterpolationMode=UMapx.Core.InterpolationMode;

namespace UMapx.Tests;

[SupportedOSPlatform("windows")]
[Trait("Category","Imaging")]
public class ImagingCompositionAuditTests
{
    public static IEnumerable<object[]> ChannelCases(){foreach(string op in new[]{"Show","Hide","Equalize"})foreach(var c in Enum.GetValues<RGBA>())yield return new object[]{op,c};}
    [Theory] [MemberData(nameof(ChannelCases))]
    public void ChannelFiltersUseTheRequestedChannel(string operation,RGBA channel)
    {
        using var source=ImagingAuditTests.Pattern();using var actual=(Bitmap)source.Clone();IBitmapFilter f=operation switch{"Show"=>new ShowChannel(channel),"Hide"=>new HideChannel(channel),_=>new EqualizeChannel(channel)};f.Apply(actual);
        for(int y=0;y<source.Height;y++)for(int x=0;x<source.Width;x++)
        {
            var p=source.GetPixel(x,y);var components=new[]{p.B,p.G,p.R,p.A};var expected=(byte[])components.Clone();int index=(int)channel;
            if(operation=="Hide")expected[index]=0;
            else if(operation=="Equalize"||channel==RGBA.Alpha)for(int i=0;i<3;i++)expected[i]=components[index];
            else for(int i=0;i<3;i++)if(i!=index)expected[i]=0;
            ImagingAuditTests.Pixel(Color.FromArgb(expected[3],expected[2],expected[1],expected[0]),actual.GetPixel(x,y));
        }
    }
    [Fact]
    public void AlphaMasksUseTheMaskIntensityAndPreserveColorChannels()
    {
        using var target=ImagingAuditTests.Pattern();using var mask=ImagingAuditTests.Pattern(alpha:false);using var original=(Bitmap)target.Clone();new AlphaChannel().Apply(target,mask);
        for(int y=0;y<target.Height;y++)for(int x=0;x<target.Width;x++){var c=original.GetPixel(x,y);var m=mask.GetPixel(x,y);ImagingAuditTests.Pixel(Color.FromArgb((m.R+m.G+m.B)/3,c.R,c.G,c.B),target.GetPixel(x,y));}
    }
    [Theory] [InlineData(0)] [InlineData(128)] [InlineData(255)]
    public void LayerMergingUsesSourceAlphaAndTheRequestedOffset(int opacity)
    {
        using var background=ImagingAuditTests.Pattern(9,7,false);using var original=(Bitmap)background.Clone();using var foreground=ImagingAuditTests.Pattern(4,3);new Merge(2,1,opacity).Apply(background,foreground);
        for(int y=0;y<7;y++)for(int x=0;x<9;x++)
        {
            var a=original.GetPixel(x,y);var expected=a;if(x>=2&&x<6&&y>=1&&y<4){var b=foreground.GetPixel(x-2,y-1);int alpha=(int)(b.A*opacity/255.0);int Mix(int av,int bv)=>(av*(255-alpha)+bv*alpha)/255;expected=Color.FromArgb(255,Mix(a.R,b.R),Mix(a.G,b.G),Mix(a.B,b.B));}ImagingAuditTests.Pixel(expected,background.GetPixel(x,y));
        }
    }
    [Theory] [InlineData(1f,-1f)] [InlineData(.25f,.75f)] [InlineData(1f,1f)]
    public void LinearImageOperationsMatchSaturatedArithmetic(float a,float b)
    {
        using var left=ImagingAuditTests.Pattern();using var original=(Bitmap)left.Clone();using var right=ImagingAuditTests.Pattern(alpha:false);new Operation(a,b).Apply(left,right);
        for(int y=0;y<left.Height;y++)for(int x=0;x<left.Width;x++){var p=original.GetPixel(x,y);var q=right.GetPixel(x,y);int Mix(int v,int w)=>Math.Clamp((int)(a*v+b*w),0,255);ImagingAuditTests.Pixel(Color.FromArgb(p.A,Mix(p.R,q.R),Mix(p.G,q.G),Mix(p.B,q.B)),left.GetPixel(x,y));}
    }
    [Theory] [InlineData(InterpolationMode.NearestNeighbor)] [InlineData(InterpolationMode.Bilinear)] [InlineData(InterpolationMode.Bicubic)]
    public void ZeroAngleRotationPreservesTheOriginalImage(InterpolationMode mode)
    {
        using var source=ImagingAuditTests.Pattern();using var actual=(Bitmap)source.Clone();new Rotate(0,mode).Apply(actual);ImagingAuditTests.Same(source,actual,1);
    }
    [Fact]
    public void IdentityPerspectiveWarpPreservesPixels()
    {
        using var source=ImagingAuditTests.Pattern(16,12);using var actual=(Bitmap)source.Clone();new PerspectiveWarp(new PointFloat(0,0),new PointFloat(16,0),new PointFloat(0,12),new PointFloat(16,12)).Apply(actual);ImagingAuditTests.Same(source,actual);
    }
    [Theory] [InlineData(false)] [InlineData(true)]
    public void ColorTransferFromTheSameImageIsAnIdentityEvenWithEqualColumnVariances(bool inverted)
    {
        using var source=new Bitmap(16,12);for(int y=0;y<12;y++)for(int x=0;x<16;x++)source.SetPixel(x,y,Color.FromArgb(20+10*y,30+8*y,40+6*y));
        using var actual=(Bitmap)source.Clone();new ColorTransfer(0,inverted).Apply(actual,source);ImagingAuditTests.Same(source,actual,1,false);
    }
    [Theory] [InlineData(false)] [InlineData(true)]
    public void MotionDetectionReportsTheFractionOfPixelsThatChanged(bool filtered)
    {
        using var first=new Bitmap(16,12);using(var g=Graphics.FromImage(first))g.Clear(Color.Black);using var second=(Bitmap)first.Clone();using(var g=Graphics.FromImage(second))g.FillRectangle(Brushes.White,0,0,8,12);
        using var detector=new MotionDetector(15,filtered);Close(0,detector.Apply(first));Close(0,detector.Apply(first));Close(.5,detector.Apply(second));detector.Reset();Close(0,detector.Apply(first));
    }
    [Fact]
    public void CanvasesAndColorReplacementProduceRequestedColors()
    {
        using var solid=new CanvasColor(32,24,Color.Red).Create();Assert.Equal(new Size(32,24),solid.Size);ImagingAuditTests.Pixel(Color.Red,solid.GetPixel(13,7));new ColorReplace(Color.Red,Color.Blue).Apply(solid);ImagingAuditTests.Pixel(Color.Blue,solid.GetPixel(13,7));
        using var gradient=new CanvasGradient(32,24,0,Color.Black,Color.White).Create();Assert.Equal(new Size(32,24),gradient.Size);int previous=-1;for(int x=0;x<32;x++){var c=gradient.GetPixel(x,12);Assert.Equal(c.R,c.G);Assert.Equal(c.R,c.B);Assert.True(c.R>=previous);previous=c.R;}
        using var mask=new CanvasColor(32,24,Color.White).Create();new MaskColorFilter(Color.Lime).Apply(solid,mask);ImagingAuditTests.Pixel(Color.Lime,solid.GetPixel(13,7));
    }
    [Fact]
    public void FusingIdenticalExposuresPreservesTheImage()
    {
        using var source=ImagingAuditTests.Pattern(32,24,false);using var copy=(Bitmap)source.Clone();using var actual=new ExposureFusion(2).Apply(source,copy);ImagingAuditTests.Same(source,actual,3,false);
    }
    [Fact]
    public void IdenticalStereoImagesHaveZeroDisparity()
    {
        using var source=ImagingAuditTests.Pattern(32,24,false);using var copy=(Bitmap)source.Clone();var actual=new StereoDisparity(5,3,1).Apply(source,copy);Assert.Equal(24,actual.GetLength(0));Assert.Equal(32,actual.GetLength(1));
        for(int y=4;y<20;y++)for(int x=8;x<28;x++)Close(0,actual[y,x]);
    }
    [Theory] [InlineData(1)] [InlineData(3)] [InlineData(5)]
    public void FractalNoiseEqualsTheSumOfItsOctaves(int octaves)
    {
        var noise=new PerlinNoise(octaves,.6f,.7f,.8f);
        foreach(float x in new[]{-3.25f,-.5f,0f,.25f,2.7f})
        {
            double expected1=0,expected2=0;float frequency=.7f,amplitude=.8f;
            for(int k=0;k<octaves;k++){var single=new PerlinNoise(1,.6f,frequency,amplitude);expected1+=single.Function(x);expected2+=single.Function2D(x,.3f);frequency*=2;amplitude*=.6f;}
            Close(expected1,noise.Function(x));Close(expected2,noise.Function2D(x,.3f));Close(0,new PerlinNoise(octaves,.6f,.7f,0).Function(x));
        }
    }
}
