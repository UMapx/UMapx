using UMapx.Video;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Video")]
public class TimeoutStreamRepairTests
{
    [Fact]
    public void MissingBaseStreamIsRejected()
        => Assert.Throws<ArgumentNullException>(() => new TimeoutStream(null!));

    [Theory]
    [InlineData(0)] [InlineData(-2)] [InlineData(int.MinValue)]
    public void InvalidDeadlinesAreRejectedWhenAssigned(int deadline)
    {
        using var input = new MemoryStream();
        using var stream = new TimeoutStream(input);
        Assert.Throws<ArgumentOutOfRangeException>(() => stream.ReadTimeout = deadline);
        Assert.Throws<ArgumentOutOfRangeException>(() => stream.WriteTimeout = deadline);
    }

    [Theory]
    [InlineData(false)] [InlineData(true)]
    public void CancellationUnrelatedToTheDeadlineIsPreserved(bool write)
    {
        using var cancelled = new CancellationTokenSource();
        cancelled.Cancel();
        using var input = new AsyncStream
        {
            ReadOperation = _ => Task.FromCanceled<int>(cancelled.Token),
            WriteOperation = _ => Task.FromCanceled(cancelled.Token)
        };
        using var stream = new TimeoutStream(input) { ReadTimeout = Timeout.Infinite, WriteTimeout = Timeout.Infinite };
        Assert.ThrowsAny<OperationCanceledException>(() =>
        {
            if (write) stream.Write(new byte[1], 0, 1); else stream.Read(new byte[1], 0, 1);
        });
    }

    [Theory]
    [InlineData(false)] [InlineData(true)]
    public async Task IdleTimeAfterSuccessDoesNotCancelTheNextOperation(bool write)
    {
        using var input = new AsyncStream();
        using var stream = new TimeoutStream(input) { ReadTimeout = 30, WriteTimeout = 30 };
        var buffer = new byte[1];
        if (write) stream.Write(buffer, 0, 1); else Assert.Equal(1, stream.Read(buffer, 0, 1));
        await Task.Delay(200);
        stream.ReadTimeout = stream.WriteTimeout = Timeout.Infinite;
        if (write) stream.Write(buffer, 0, 1); else Assert.Equal(1, stream.Read(buffer, 0, 1));
    }

    [Fact]
    public void WritesUseTheWriteDeadlineWhenReadsHaveNoDeadline()
    {
        using var input = new AsyncStream { WriteOperation = token => Task.Delay(500, token) };
        using var stream = new TimeoutStream(input) { ReadTimeout = Timeout.Infinite, WriteTimeout = 30 };
        Assert.Throws<TimeoutException>(() => stream.Write(new byte[1], 0, 1));
        input.WriteOperation = _ => Task.CompletedTask;
        stream.Write(new byte[1], 0, 1);
    }

    [Fact]
    public void WritesWithNoDeadlineIgnoreTheReadDeadline()
    {
        using var input = new AsyncStream { WriteOperation = token => Task.Delay(100, token) };
        using var stream = new TimeoutStream(input) { ReadTimeout = 30, WriteTimeout = Timeout.Infinite };
        stream.Write(new byte[1], 0, 1);
    }

    [Fact]
    public void ReadTimeoutDoesNotPoisonSubsequentReads()
    {
        using var input = new AsyncStream { ReadOperation = async token => { await Task.Delay(500, token); return 1; } };
        using var stream = new TimeoutStream(input) { ReadTimeout = 30 };
        Assert.Throws<TimeoutException>(() => stream.Read(new byte[1], 0, 1));
        input.ReadOperation = _ => Task.FromResult(1);
        Assert.Equal(1, stream.Read(new byte[1], 0, 1));
    }

    [Theory]
    [InlineData(false)] [InlineData(true)]
    public void TransportFailuresKeepTheirOriginalException(bool write)
    {
        var failure = new IOException("Synthetic transport failure.");
        using var input = new AsyncStream
        {
            ReadOperation = _ => Task.FromException<int>(failure),
            WriteOperation = _ => Task.FromException(failure)
        };
        using var stream = new TimeoutStream(input);
        var actual = Assert.Throws<IOException>(() =>
        {
            if (write) stream.Write(new byte[1], 0, 1); else stream.Read(new byte[1], 0, 1);
        });
        Assert.Same(failure, actual);
    }

    [Fact]
    public void NativeStreamTimeoutsReceiveTheConfiguredDeadlines()
    {
        using var input = new NativeTimeoutStream();
        using var stream = new TimeoutStream(input) { ReadTimeout = 123, WriteTimeout = 456 };
        stream.Read(new byte[1], 0, 1);
        stream.Write(new byte[1], 0, 1);
        Assert.Equal(123, input.ReadTimeoutAtRead);
        Assert.Equal(456, input.WriteTimeoutAtWrite);
    }

    private sealed class AsyncStream : MemoryStream
    {
        public Func<CancellationToken, Task<int>> ReadOperation { get; set; } = _ => Task.FromResult(1);
        public Func<CancellationToken, Task> WriteOperation { get; set; } = _ => Task.CompletedTask;

        public override Task<int> ReadAsync(byte[] buffer, int offset, int count, CancellationToken cancellationToken)
            => cancellationToken.IsCancellationRequested ? Task.FromCanceled<int>(cancellationToken) : ReadOperation(cancellationToken);

        public override Task WriteAsync(byte[] buffer, int offset, int count, CancellationToken cancellationToken)
            => cancellationToken.IsCancellationRequested ? Task.FromCanceled(cancellationToken) : WriteOperation(cancellationToken);
    }

    private sealed class NativeTimeoutStream : MemoryStream
    {
        public override bool CanTimeout => true;
        public override int ReadTimeout { get; set; } = Timeout.Infinite;
        public override int WriteTimeout { get; set; } = Timeout.Infinite;
        public int ReadTimeoutAtRead { get; private set; }
        public int WriteTimeoutAtWrite { get; private set; }
        public override int Read(byte[] buffer, int offset, int count) { ReadTimeoutAtRead = ReadTimeout; return 0; }
        public override void Write(byte[] buffer, int offset, int count) { WriteTimeoutAtWrite = WriteTimeout; }
    }
}
