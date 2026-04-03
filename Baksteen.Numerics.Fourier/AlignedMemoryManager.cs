namespace Baksteen.Numerics.Fourier;

using System;
using System.Buffers;
using System.Runtime.InteropServices;

public unsafe class AlignedMemoryManager<T> : MemoryManager<T> where T : struct
{
    private readonly void* _ptr;
    private readonly int _length;
    private bool _disposed;

    public AlignedMemoryManager(int length, int alignment)
    {
        var numBytes = (nuint)(length * Marshal.SizeOf<T>());
        _ptr = NativeMemory.AlignedAlloc(numBytes, (nuint)alignment);
        NativeMemory.Clear(_ptr, numBytes);
        _length = length;
    }

    public override Span<T> GetSpan() => new(_ptr, _length);

    public override MemoryHandle Pin(int elementIndex = 0)
        => new(((T*)_ptr) + elementIndex);

    public override void Unpin()
    {
        // Nothing to do — the whole array is pinned
    }

    protected override void Dispose(bool disposing)
    {
        if (!_disposed)
        {
            NativeMemory.AlignedFree(_ptr);
            _disposed = true;
        }
    }
}


/*
var manager = new PinnedArrayMemoryManager<byte>(1024);

Memory<byte> mem = manager.Memory;
Span<byte> span = mem.Span;

nuint addr = manager.Address;

Console.WriteLine($"Address: 0x{addr:X}");
*/