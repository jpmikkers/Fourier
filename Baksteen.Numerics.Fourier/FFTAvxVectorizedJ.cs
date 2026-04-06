namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

public class FFTAvxVectorizedJ
{
    private AlignedMemoryManager<Complex> _alignedMemoryManager;
    private readonly Memory<Complex> _wtable;

    public FFTAvxVectorizedJ(int length)
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
                    Vectorized.Sse2ButterflyOne(ptmp, ptmp + 1);
                    ptmp += partStride;
                }

                butterfliesPerPart <<= 1;
                nrOfParts >>= 1;
                partStride <<= 1;
            }

            if (isInverse)
            {
                while (nrOfParts > 1)
                {
                    var ptmp = vptr;

                    for (var p = 0; p < nrOfParts; p++)
                    {
                        Vectorized.AvxButterfliesConjugate(ptmp, wptr, butterfliesPerPart, nrOfParts);
                        ptmp += partStride;
                    }

                    butterfliesPerPart <<= 1;
                    nrOfParts >>= 1;
                    partStride <<= 1;
                }

                if (nrOfParts > 0)
                {
                    ButterfliesConjugateScaled(vptr, wptr, butterfliesPerPart, nrOfParts, Math.ScaleB(1.0, -BitOperations.Log2((uint)data.Length)));
                }
            }
            else
            {
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
        }
    }

    private static unsafe void ButterfliesConjugateScaled(Complex* peven, Complex* pw, int numbtf, int wstride, double scalefactor)
    {
        var scaler = Vector256.Create(scalefactor);

        for (var t = 0; t < numbtf; t += 2)
        {
            var evens = Avx.LoadAlignedVector256((double*)peven);   // Avx.LoadVector256((double*)peven);
            var odds = Avx.LoadAlignedVector256((double*)(peven + numbtf)); // Avx.LoadVector256((double*)(peven + numbtf));
            var w = Vector256.Create(Sse2.LoadAlignedVector128((double*)pw), Sse2.LoadAlignedVector128((double*)(pw + wstride)));

            var oddw = Vectorized.ComplexMulAvxConjugate(odds, w);
            pw += (wstride << 1);

            var re = Avx.Add(evens, oddw);
            Avx.StoreAligned((double*)peven, Avx.Multiply(re, scaler));

            var ro = Avx.Subtract(evens, oddw);
            Avx.StoreAligned((double*)(peven + numbtf), Avx.Multiply(ro, scaler));

            peven += 2;
        }
    }
}
