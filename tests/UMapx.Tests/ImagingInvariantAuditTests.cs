using System.Runtime.Versioning;
using System.Drawing;
using UMapx.Core;
using UMapx.Imaging;
using UMapx.Transform;
using Xunit;

namespace UMapx.Tests;

[Trait("Category","Imaging")]
[SupportedOSPlatform("windows")]
public class ImagingInvariantAuditTests
{
    static IBitmapFilter Smoother(string name)=>name switch
    {
        "BoxBlur"=>new BoxBlur(3),"GaussianBlur"=>new GaussianBlur(3,5,1,1.5f),"Erosion"=>new Erosion(3),"Dilatation"=>new Dilatation(3),"Median"=>new Median(3),"Opening"=>new Opening(3),"Closing"=>new Closing(3),
        "Kuwahara"=>new Kuwahara(3),"OilPainting"=>new OilPainting(3),"SmartOilPainting"=>new SmartOilPainting(3,64,0),"Wiener"=>new Wiener(3),"Pixelate"=>new Pixelate(3),"Jitter"=>new Jitter(3),"Noise"=>new Noise(3),"Water"=>new Water(3),
        "PhotoFilter"=>new PhotoFilter(Color.Red,0),"YUVPhotoFilter"=>new YUVPhotoFilter(Color.Red,0),"Temperature"=>new TemperatureCorrection(6500,0),"Texturer"=>new Texturer(new float[32,48],0),
        "AdditiveNoise"=>new AdditiveNoise(0),"SaltAndPepper"=>new SaltAndPepper(0),"FlatField"=>new FlatFieldCorrection(3),"BitmapFilter"=>new BitmapFilter(new GuidedFilter(3),Space.RGB),
        _=>throw new ArgumentException(name)
    };
    public static IEnumerable<object[]> ConstantCases(){foreach(string n in new[]{"BoxBlur","GaussianBlur","Erosion","Dilatation","Median","Opening","Closing","Kuwahara","OilPainting","SmartOilPainting","Wiener","Pixelate","Jitter","Noise","Water","PhotoFilter","YUVPhotoFilter","Temperature","Texturer","AdditiveNoise","SaltAndPepper","FlatField","BitmapFilter"})foreach(int v in new[]{0,77,255})yield return new object[]{n,v};}
    [Theory] [MemberData(nameof(ConstantCases))]
    public void SmoothingAndNeutralEffectsPreserveUniformFields(string name,int value)
    {
        using var image=new Bitmap(48,32);using(var g=Graphics.FromImage(image))g.Clear(Color.FromArgb(255,value,value,value));Smoother(name).Apply(image);
        int tolerance=name is "SmartOilPainting"?3:1;
        for(int y=0;y<image.Height;y++)for(int x=0;x<image.Width;x++)ImagingAuditTests.Pixel(Color.FromArgb(255,value,value,value),image.GetPixel(x,y),tolerance,false);
    }
    static IBitmapFilter Enhancer(string name,Space space)=>name switch
    {
        "Brightness"=>new BrightnessCorrection(.3f,space),"Contrast"=>new ContrastCorrection(.3f,space),"Gamma"=>new GammaCorrection(2,space),"Exposure"=>new ExposureCorrection(128,space),
        "Linear"=>new LinearCorrection(.2f,space),"Levels"=>new LevelsCorrection(new RangeFloat(.1f,.9f),new RangeFloat(0,1),space),"Log"=>new LogCorrection(3,.2f,space),"Sin"=>new SinCorrection(.2f,space),"Cos"=>new CosCorrection(.2f,space),"Shift"=>new ShiftCorrection(.2f,space),
        "Quantization"=>new Quantization(4,space),"Threshold"=>new Threshold(.5f,space),"Invert"=>new InvertChannels(space),"HistogramStretch"=>new HistogramStretch(0,1,space),"ContrastEnhancement"=>new ContrastEnhancement(.3f,space),
        "LocalContrast"=>new LocalContrastEnhancement(3,space),"LocalInversion"=>new LocalContrastInversion(3,space),"LocalHistogramStretch"=>new LocalHistogramStretch(3,space),"LocalThreshold"=>new LocalThreshold(3,space),
        "Homomorphic"=>new HomomorphicEnhancement(3,space),"KsiContrast"=>new KsiContrastEnhancement(3,space),"Retinex"=>new SingleScaleRetinex(3,space),"ShadowsHighlights"=>new ShadowsHighlightsCorrection(3,space),
        _=>throw new ArgumentException(name)
    };
    public static IEnumerable<object[]> GrayCases(){foreach(string n in new[]{"Brightness","Contrast","Gamma","Exposure","Linear","Levels","Log","Sin","Cos","Shift","Quantization","Threshold","Invert","HistogramStretch","ContrastEnhancement","LocalContrast","LocalInversion","LocalHistogramStretch","LocalThreshold","Homomorphic","KsiContrast","Retinex","ShadowsHighlights"})foreach(var s in Enum.GetValues<Space>())yield return new object[]{n,s};}
    [Theory] [MemberData(nameof(GrayCases))]
    public void AchromaticEnhancementDoesNotIntroduceColorIntoGrayImages(string name,Space space)
    {
        using var image=new Bitmap(48,32);for(int y=0;y<32;y++)for(int x=0;x<48;x++){int v=(7*x+13*y)%256;image.SetPixel(x,y,Color.FromArgb(255,v,v,v));}
        Enhancer(name,space).Apply(image);for(int y=0;y<32;y++)for(int x=0;x<48;x++){var p=image.GetPixel(x,y);Assert.InRange(Math.Abs(p.R-p.G),0,2);Assert.InRange(Math.Abs(p.R-p.B),0,2);}
    }
    [Theory] [InlineData("TopHat")] [InlineData("BottomHat")] [InlineData("EdgeGlow")] [InlineData("Canny")] [InlineData("FreiChen")]
    public void EdgeDetectorsAndMorphologicalResidualsVanishOnAUniformField(string name)
    {
        using var image=new Bitmap(48,32);using(var g=Graphics.FromImage(image))g.Clear(Color.FromArgb(77,77,77));IBitmapFilter f=name switch{"TopHat"=>new TopHat(3),"BottomHat"=>new BottomHat(3),"EdgeGlow"=>new EdgeGlow(3),"Canny"=>new CannyEdgeDetector(),_=>new FreiChen()};f.Apply(image);
        // The interior does not depend on each operator's choice of outside-image boundary condition.
        for(int y=8;y<24;y++)for(int x=8;x<40;x++)ImagingAuditTests.Pixel(Color.Black,image.GetPixel(x,y),1,false);
    }
    [Theory] [InlineData("HistogramEqualization")] [InlineData("LocalHistogramEqualization")] [InlineData("CLAHE")] [InlineData("ToneDither")]
    public void HistogramAndDitheringFiltersRetainAchromaticChannels(string name)
    {
        using var image=new Bitmap(48,32);for(int y=0;y<32;y++)for(int x=0;x<48;x++){int v=(x*7+y*13)%256;image.SetPixel(x,y,Color.FromArgb(v,v,v));}
        IBitmapFilter f=name switch{"HistogramEqualization"=>new HistogramEqualization(),"LocalHistogramEqualization"=>new LocalHistogramEqualization(3),"CLAHE"=>new CLAHE(4,4),_=>ToneDiffusionDithering.Bayer()};f.Apply(image);
        for(int y=0;y<32;y++)for(int x=0;x<48;x++){var p=image.GetPixel(x,y);Assert.Equal(p.R,p.G);Assert.Equal(p.R,p.B);if(name=="ToneDither")Assert.Contains(p.R,new byte[]{0,255});}
    }
    [Theory] [InlineData(false)] [InlineData(true)]
    public void FlatHeightMapsProduceNormalsPerpendicularToTheImage(bool sobel)
    {
        using var image=new Bitmap(16,12);using(var g=Graphics.FromImage(image))g.Clear(Color.FromArgb(77,77,77));new NormalBump(2){UseSobel=sobel}.Apply(image);
        for(int y=2;y<10;y++)for(int x=2;x<14;x++)ImagingAuditTests.Pixel(Color.FromArgb(127,127,255),image.GetPixel(x,y),1,false);
    }
}
