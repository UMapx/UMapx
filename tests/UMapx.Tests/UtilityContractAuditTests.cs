using System.Drawing;
using System.Runtime.Versioning;
using UMapx.Core;
using UMapx.Imaging;
using UMapx.Video;
using Xunit;

namespace UMapx.Tests;

[Trait("Category","Contract")]
public class UtilityContractAuditTests
{
    [Fact]
    public void XmlRoundTripPreservesNumericArrayValues()
    {
        float[] expected={-1.25f,0,1e-30f,1e30f};using var stream=new MemoryStream();Xml.Save(stream,expected);stream.Position=0;Assert.Equal(expected,(float[])Xml.Open(stream,typeof(float[])));
    }
    [SupportedOSPlatform("windows")]
    [Fact]
    public void UnstartedVideoSourcesExposeConfigurationWithoutAccessingDevicesOrNetwork()
    {
        var jpeg=new JPEGStream("http://localhost.invalid/frame.jpg"){Login="audit",Password="test",RequestTimeout=1234,FrameInterval=75,SeparateConnectionGroup=true,PreventCaching=false};
        Assert.Equal(1234,jpeg.RequestTimeout);Assert.Equal(75,jpeg.FrameInterval);Assert.Equal("audit",jpeg.Login);Assert.Equal("test",jpeg.Password);Assert.True(jpeg.SeparateConnectionGroup);Assert.False(jpeg.PreventCaching);
        var mjpeg=new MJPEGStream("http://localhost.invalid/frames"){Login="audit",Password="test",RequestTimeout=1234,SeparateConnectionGroup=true};
        Assert.Equal(1234,mjpeg.RequestTimeout);Assert.Equal("audit",mjpeg.Login);Assert.Equal("test",mjpeg.Password);Assert.True(mjpeg.SeparateConnectionGroup);
        var screen=new ScreenCaptureStream(new Rectangle(0,0,16,12),75);Assert.Equal(new Rectangle(0,0,16,12),screen.Region);Assert.Equal(75,screen.FrameInterval);
        foreach(IVideoSource source in new IVideoSource[]{jpeg,mjpeg,screen}){Assert.False(source.IsRunning);Assert.Equal(0,source.FramesReceived);Assert.Equal(0,source.BytesReceived);source.SignalToStop();}
        var error=new VideoException("test failure");var info=new VideoSourceErrorEventArgs("frame unavailable",error);Assert.Equal("frame unavailable",info.Description);Assert.Same(error,info.Exception);Assert.Null(new VideoSourceErrorEventArgs("no frame").Exception);
    }
    [Theory] [InlineData(255)] [InlineData(256)] [InlineData(300)]
    public void DepthHistogramEqualizationCountsMoreThan65535PixelsWithoutOverflow(int side)
    {
        var depth=new ushort[side,side];var actual=depth.Equalize();foreach(ushort value in actual)Assert.Equal(ushort.MaxValue,value);
    }
}
