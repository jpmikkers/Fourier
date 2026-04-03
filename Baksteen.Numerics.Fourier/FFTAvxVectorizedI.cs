namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

public class FFTAvxVectorizedI
{
    private AlignedMemoryManager<Complex> _alignedMemoryManager;
    private readonly Memory<Complex> _wtable;

    public FFTAvxVectorizedI(int length)
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
                var ptmp = vptr;

                // no mul needed in the first combination layer, every w is 1+0i and -1+0i
                for (var p = 0; p < nrOfParts; p++)
                {
                    ButterflyOne(ptmp, ptmp + 1);
                    ptmp += partStride;
                }

                butterfliesPerPart <<= 1;
                nrOfParts >>= 1;
                partStride <<= 1;
            }

            while (nrOfParts > 0)
            {
                var ptmp = vptr;

                for (var p = 0; p < nrOfParts; p++)
                {
                    Butterflies(ptmp, wptr, butterfliesPerPart, nrOfParts);
                    ptmp += partStride;
                }

                butterfliesPerPart <<= 1;
                nrOfParts >>= 1;
                partStride <<= 1;
            }
        }

        if (isInverse) FFTUtils.Scale(data);
    }

    private static unsafe void Butterflies(Complex* peven, Complex* pw, int numbtf, int wstride)
    {
        for (var t = 0; t < numbtf; t += 2)
        {
            var evens = Avx.LoadAlignedVector256((double*)peven);   // Avx.LoadVector256((double*)peven);
            var odds = Avx.LoadAlignedVector256((double*)(peven + numbtf)); // Avx.LoadVector256((double*)(peven + numbtf));
            var w = Vector256.Create(Sse2.LoadAlignedVector128((double*)pw), Sse2.LoadAlignedVector128((double*)(pw + wstride)));
            var oddw = Vectorized.ComplexMulAvx(odds, w);
            pw += (wstride << 1);

            var re = Avx.Add(evens, oddw);
            Avx.StoreAligned((double*)peven, re);

            var ro = Avx.Subtract(evens, oddw);
            Avx.StoreAligned((double*)(peven + numbtf), ro);

            peven += 2;
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private unsafe void ButterflyOne(Complex* peven, Complex* podd)
    {
        var even0 = Sse2.LoadAlignedVector128((double*)peven);
        var odd0 = Sse2.LoadAlignedVector128((double*)podd);
        var re0 = Sse2.Add(even0, odd0);
        var ro0 = Sse2.Subtract(even0, odd0);
        Sse2.StoreAligned((double*)peven, re0);
        Sse2.StoreAligned((double*)podd, ro0);
    }
}
