namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics.X86;

public class FFTAvxVectorizedK
{
    private AlignedMemoryManager<Complex> _alignedMemoryManager;
    private readonly Memory<Complex> _wtable;

    public FFTAvxVectorizedK(int length)
    {
        if (!Avx.IsSupported || !Fma.IsSupported || !Sse2.IsSupported)
        {
            throw new NotSupportedException("need avx, fma, sse2 support for this fft implementation");
        }

        if (!BitOperations.IsPow2(length))
        {
            throw new ArgumentException("fft not a power of two", nameof(length));
        }

        _alignedMemoryManager = new AlignedMemoryManager<Complex>(length / 2, Marshal.SizeOf<Complex>());
        _wtable = _alignedMemoryManager.Memory;

        for (int t = 0; t < _wtable.Length; t++)
        {
            _wtable.Span[t] = Complex.FromPolarCoordinates(1, -(Math.Tau * t) / length);
        }
    }

    public unsafe void FastFourierTransform(Span<Complex> data, bool isInverse)
    {
        if (!BitOperations.IsPow2(data.Length))
        {
            throw new ArgumentException("fft not a power of two", nameof(data));
        }

        FFTUtils.AssertAlignment(data, 256 / 8);

        Reorder.Shuffle(data);

        var butterfliesPerPart = 1;             // a single butterfly does 2 angles, +w and -w (=w+pi radians)
        var nrOfParts = data.Length >> 1;       // so the first layer is len/2 parts of single butterflies
        var partStride = 2;

        fixed (Complex* vptr = data)
        fixed (Complex* wptr = _wtable.Span)
        {
            if (nrOfParts > 0)
            {
                //Vectorized.AvxButterfliesOnes(vptr, nrOfParts);
                Vectorized.Sse2ButterfliesOnes(vptr, nrOfParts);
                butterfliesPerPart <<= 1;
                nrOfParts >>= 1;
                partStride <<= 1;
            }

            while (nrOfParts > 0)
            {
                var ptmp = vptr;

                for (var p = 0; p < nrOfParts; p++)
                {
                    Vectorized.AvxButterflies(ptmp, wptr, butterfliesPerPart, nrOfParts);
                    ptmp += partStride;
                }

                butterfliesPerPart <<= 1;
                nrOfParts >>= 1;
                partStride <<= 1;
            }
        }

        if (isInverse) FFTUtils.Scale(data);
    }
}
