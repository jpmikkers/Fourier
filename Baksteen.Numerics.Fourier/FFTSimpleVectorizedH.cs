namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

public class FFTSimpleVectorizedH
{
    private AlignedMemoryManager<Complex> _alignedMemoryManager;
    private double[] _wtable;

    public FFTSimpleVectorizedH(int length)
    {
        if (!Sse2.IsSupported || !Sse3.IsSupported)
        {
            throw new NotSupportedException("need sse2 and sse3 support for this fft implementation");
        }

        if (!BitOperations.IsPow2(length))
        {
            throw new ArgumentException("fft not a power of two", nameof(length));
        }

        _wtable = new double[length];

        for (int t = 0; t < _wtable.Length; t++)
        {
            _wtable[t] = Math.Sin((Math.Tau * t) / length);
        }
    }

    public unsafe void FastFourierTransform(Span<Complex> data, bool isInverse)
    {
        if (!BitOperations.IsPow2(data.Length))
        {
            throw new ArgumentException("fft not a power of two", nameof(data));
        }

        ref Complex r = ref MemoryMarshal.GetReference(data);
        nint address = (nint)Unsafe.AsPointer(ref r);

        if ((address & (Marshal.SizeOf<Complex>() - 1)) != 0)
        {
            throw new ArgumentException("data not properly aligned for vectorized fft, length: " + data.Length, nameof(data));
        }

        Reorder.Shuffle(data);

        var butterfliesPerPart = 1;             // a single butterfly does 2 angles, +w and -w (=w+pi radians)
        var nrOfParts = data.Length >> 1;       // so the first layer is len/2 parts of single butterflies
        var partStride = 2;

        fixed (Complex* vptr = data)
        fixed (double* wptr = _wtable)
        {
            var wr = &wptr[_wtable.Length >> 2];                        // always cosine
            var dw = isInverse ? -(_wtable.Length >> 2) : (_wtable.Length >> 2);     // +sin for inverse, -sin for forward

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
                    Butterflies(ptmp, wr, dw, butterfliesPerPart, nrOfParts);
                    ptmp += partStride;
                }

                butterfliesPerPart <<= 1;
                nrOfParts >>= 1;
                partStride <<= 1;
            }
        }

        if (isInverse) FFTUtils.Scale(data);
    }

    private static unsafe void Butterflies(Complex* peven, double* wr, int deltaw, int numbtf, int wstride)
    {
        for (var t = 0; t < numbtf; t++)
        {
            var w = Vector128.Create(*wr, *(wr + deltaw));   // this is 5% slower
            var podd = peven + numbtf;
            var even0 = Sse2.LoadAlignedVector128((double*)peven);
            var odd0 = ComplexMulSse2b(
                Sse2.LoadAlignedVector128((double*)podd),
                w
            );
            var re0 = Sse2.Add(even0, odd0);
            var ro0 = Sse2.Subtract(even0, odd0);
            wr += wstride;
            Sse2.StoreAligned((double*)peven, re0);
            peven++;
            Sse2.StoreAligned((double*)podd, ro0);
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private unsafe void Butterfly(Complex* peven, Complex* podd, Vector128<double> w)
    {
        var even0 = Sse2.LoadAlignedVector128((double*)peven);
        var odd0 = ComplexMulSse2b(Sse2.LoadAlignedVector128((double*)podd), w);
        var re0 = Sse2.Add(even0, odd0);
        var ro0 = Sse2.Subtract(even0, odd0);
        Sse2.StoreAligned((double*)peven, re0);
        Sse2.StoreAligned((double*)podd, ro0);
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

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Vector128<double> ComplexMulSse2b(Vector128<double> ab, Vector128<double> cd)
    {
        // Multiplication:  (a + bi)(c + di) = (ac -bd) + (bc + ad)i
        var ac_bd = Sse2.Multiply(ab, cd);
        var bc_ad = Sse2.Multiply(Sse2.Shuffle(ab, ab, 0x01), cd);
        return Sse3.AddSubtract(Sse2.UnpackLow(ac_bd, bc_ad), Sse2.UnpackHigh(ac_bd, bc_ad));  // (ac - bd) (bc + ad)
    }

    // using AVX.Permute, not faster than ComplexMulSse2b
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Vector128<double> ComplexMulSse2c(Vector128<double> ab, Vector128<double> cd)
    {
        // Multiplication:  (a + bi)(c + di) = (ac -bd) + (bc + ad)i
        var ac_bd = Avx.Multiply(ab, cd);
        var bc_ad = Sse2.Multiply(Avx.Permute(ab, 0x01), cd);
        return Sse3.AddSubtract(Sse2.UnpackLow(ac_bd, bc_ad), Sse2.UnpackHigh(ac_bd, bc_ad));  // (ac - bd) (bc + ad)
    }
}
