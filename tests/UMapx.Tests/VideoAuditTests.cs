using System.Runtime.Versioning;
using System.Drawing;
using System.Drawing.Imaging;
using System.Net;
using System.Reflection;
using System.Text;
using UMapx.Video;
using Xunit;

namespace UMapx.Tests;

[Trait("Category","Video")]
[SupportedOSPlatform("windows")]
public class VideoAuditTests
{
    sealed class Response(string contentType):WebResponse{public override string ContentType{get;set;}=contentType;}
    [Theory]
    [InlineData("multipart/x-mixed-replace; boundary=frame","frame")]
    [InlineData("multipart/x-mixed-replace; boundary=\"frame\"","frame")]
    [InlineData("multipart/x-mixed-replace; boundary=frame; charset=utf-8","frame")]
    [InlineData("Multipart/X-Mixed-Replace; Boundary=frame","frame")]
    [InlineData("application/octet-stream","")]
    public void MultipartBoundaryParsingHonorsMimeParameters(string contentType,string expected)
    {
        using var response=new Response(contentType);var boundary=Boundary.FromResponse(response);Assert.Equal(expected,boundary.Content);Assert.Equal(Encoding.ASCII.GetBytes(expected),(byte[])boundary);
    }
    [Fact]
    public void UnsupportedMimeTypesAreRejected()
    {
        using var response=new Response("text/html");Assert.Throws<ArgumentException>(()=>Boundary.FromResponse(response));
        var b=new Boundary("frame");Assert.False(b.IsChecked);b.Prepend('-');b.Prepend('-');b.IsChecked=true;Assert.True(b.IsValid);Assert.Equal("--frame",(string)b);
    }
    static byte[] Jpeg(Color color)
    {
        using var b=new Bitmap(8,6);using(var g=Graphics.FromImage(b))g.Clear(color);using var s=new MemoryStream();b.Save(s,ImageFormat.Jpeg);return s.ToArray();
    }
    sealed class ChunkStream(byte[] data,int chunk):MemoryStream(data)
    {
        public override int Read(byte[] buffer,int offset,int count)=>base.Read(buffer,offset,Math.Min(count,chunk));
    }
    [Theory] [InlineData(1)] [InlineData(2)] [InlineData(3)] [InlineData(7)] [InlineData(63)] [InlineData(1024)]
    public void MjpegFramesSurviveArbitraryTransportChunkBoundaries(int chunk)
    {
        var header=Encoding.ASCII.GetBytes("--frame\r\nContent-Type: image/jpeg\r\n\r\n");
        var data=header.Concat(Jpeg(Color.Red)).Concat(Encoding.ASCII.GetBytes("\r\n")).Concat(header).Concat(Jpeg(Color.Blue)).Concat(Encoding.ASCII.GetBytes("\r\n--frame--\r\n")).ToArray();
        var boundary=new Boundary("--frame"){IsChecked=true};var parser=new MJPEGStreamParser(boundary,new byte[]{255,216,255});using var input=new ChunkStream(data,chunk);var colors=new List<Color>();
        while(input.Position<input.Length)
        {
            parser.Read(input);parser.DetectFrame();
            while(parser.HasFrame){using var frame=parser.GetFrame();Assert.Equal(new Size(8,6),frame.Size);colors.Add(frame.GetPixel(4,3));parser.RemoveFrame();parser.DetectFrame();}
        }
        Assert.Equal(2,colors.Count);ImagingAuditTests.Pixel(Color.Red,colors[0],4,false);ImagingAuditTests.Pixel(Color.Blue,colors[1],4,false);
    }
    [Theory] [InlineData(1)] [InlineData(17)] [InlineData(71)]
    public void BytePatternSearchMatchesASimpleBoundedSearch(int seed)
    {
        var type=typeof(Boundary).Assembly.GetType("UMapx.Video.ByteArrayUtils")!;var find=type.GetMethod("Find")!;var compare=type.GetMethod("Compare")!;var rng=new Random(seed);var source=Enumerable.Range(0,73).Select(_=>(byte)rng.Next(4)).ToArray();
        foreach(int n in new[]{1,2,3,7})foreach(int start in new[]{0,3,30})foreach(int length in new[]{0,1,15,40})
        {
            var needle=Enumerable.Range(0,n).Select(_=>(byte)rng.Next(4)).ToArray();int expected=-1;for(int i=start;i<=start+length-n;i++)if(source.Skip(i).Take(n).SequenceEqual(needle)){expected=i;break;}
            Assert.Equal(expected,(int)find.Invoke(null,new object[]{source,needle,start,length})!);
            if(expected>=0)Assert.True((bool)compare.Invoke(null,new object[]{source,needle,expected})!);
        }
    }
    [Fact]
    public void TimeoutStreamCanWrapAnOrdinaryReadableStream()
    {
        using var input=new MemoryStream(new byte[]{3,5,7});using var stream=new TimeoutStream(input);var buffer=new byte[3];Assert.Equal(3,stream.Read(buffer,0,3));Assert.Equal(new byte[]{3,5,7},buffer);
    }
    [Theory] [InlineData(false)] [InlineData(true)]
    public async Task SyntheticVideoSourcesDeliverFramesAndStopWithoutExternalDevices(bool depth)
    {
        using var bitmap=ImagingAuditTests.Pattern(12,10,false);var received=new TaskCompletionSource<Size>(TaskCreationOptions.RunContinuationsAsynchronously);var depthReceived=new TaskCompletionSource<int>(TaskCreationOptions.RunContinuationsAsynchronously);
        var capabilities=new VideoCapabilities(new Size(12,10),10,10,32);IVideoSource source;
        if(depth){var d=new VideoImageDepthSource((Bitmap)bitmap.Clone(),new ushort[10,12]);d.VideoResolution=capabilities;d.DepthResolution=capabilities;d.NewDepth+=(_,e)=>depthReceived.TrySetResult(e.Depth.Length);source=d;}
        else{var v=new VideoImageSource((Bitmap)bitmap.Clone());v.VideoResolution=capabilities;source=v;}
        source.NewFrame+=(_,e)=>received.TrySetResult(e.Frame.Size);
        try{source.Start();Assert.True(source.IsRunning);Assert.Equal(new Size(12,10),await received.Task.WaitAsync(TimeSpan.FromSeconds(3)));if(depth)Assert.Equal(120,await depthReceived.Task.WaitAsync(TimeSpan.FromSeconds(3)));source.SignalToStop();Assert.False(source.IsRunning);Assert.True(source.FramesReceived>0);Assert.True(source.BytesReceived>0);}
        finally{source.SignalToStop();((IDisposable)source).Dispose();}
    }
}
