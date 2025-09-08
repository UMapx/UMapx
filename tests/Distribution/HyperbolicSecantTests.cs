using UMapx.Distribution;
using UMapx.Core;

namespace UMapx.Tests;

public class HyperbolicSecantTests
{
    [Fact]
    public void Entropy_Equals_Log4()
    {
        var dist = new HyperbolicSecant();
        Assert.Equal(Maths.Log(4f), dist.Entropy);
    }
}