using UMapx.Core;
using UMapx.Imaging;
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
    [Theory] [InlineData(255)] [InlineData(256)] [InlineData(300)]
    public void DepthHistogramEqualizationCountsMoreThan65535PixelsWithoutOverflow(int side)
    {
        var depth=new ushort[side,side];var actual=depth.Equalize();foreach(ushort value in actual)Assert.Equal(ushort.MaxValue,value);
    }
}
