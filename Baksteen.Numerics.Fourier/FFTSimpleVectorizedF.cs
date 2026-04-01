namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

public class FFTSimpleVectorizedF
{
    private AlignedMemoryManager<Complex> _alignedMemoryManager;
    private Memory<Complex> _wtable;

    public FFTSimpleVectorizedF(int length)
    {
        if (!Sse2.IsSupported || !Sse3.IsSupported)
        {
            throw new NotSupportedException("need sse2 and sse3 support for this fft implementation");
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

        ref Complex r = ref MemoryMarshal.GetReference(data);
        nint address = (nint)Unsafe.AsPointer(ref r);

        if ((address & (Marshal.SizeOf<Complex>() - 1)) != 0)
        {
            throw new ArgumentException("data not properly aligned for vectorized fft, length: " + data.Length, nameof(data));
        }

        Reorder.Shuffle(data);

        var butterfliesPerPart = 1;             // a single butterfly does 2 angles, +w and -w (=w+pi radians)
        var nrOfParts = data.Length >> 1;       // so the first layer is len/2 parts of single butterflies
        var rotationLookupIndex = 1;
        var rotationIndexStep = data.Length >> 1;

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
                    ptmp += 2;
                }

                butterfliesPerPart <<= 1;
                nrOfParts >>= 1;
                rotationLookupIndex++;
                rotationIndexStep >>= 1;
            }

            while (nrOfParts > 0)
            {
                if (butterfliesPerPart > nrOfParts)
                {
                    for (var p = 0; p < nrOfParts; p++)
                    {
                        var evenindex = p << rotationLookupIndex;
                        var oddindex = evenindex + butterfliesPerPart;
                        var peven = (vptr + evenindex);
                        var podd = (vptr + oddindex);
                        var pw = wptr;

                        ButterflyOne(peven, podd);

                        peven++;
                        podd++;
                        pw += rotationIndexStep;

                        for (var a = 1; a < butterfliesPerPart; a++)
                        {
                            Butterfly(peven, podd, Vector128.LoadAligned((double*)pw));
                            peven++;
                            podd++;
                            pw += rotationIndexStep;
                        }
                    }
                }
                else
                {
                    var pw = wptr;

                    for (var a = 0; a < butterfliesPerPart; a++)
                    {
                        var w = Vector128.LoadAligned((double*)pw);
                        for (var p = 0; p < nrOfParts; p++)
                        {
                            var peven = (vptr + (p << rotationLookupIndex) + a);
                            Butterfly(peven, peven + butterfliesPerPart, w);
                        }
                        pw += rotationIndexStep;
                    }
                }

                butterfliesPerPart <<= 1;
                nrOfParts >>= 1;
                rotationLookupIndex++;
                rotationIndexStep >>= 1;
            }
        }

        if (isInverse) FFTUtils.Scale(data);
    }

#if USEGENERICLOAD
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private unsafe void Butterfly(double* peven, double* podd, Vector128<double> w)
    {
        var even0 = Vector128.LoadAligned(peven);
        var odd0 = Vectorized.ComplexMulSse2(Vector128.LoadAligned(podd), w);
        var re0 = Sse2.Add(even0, odd0);
        var ro0 = Sse2.Subtract(even0, odd0);
        Vector128.StoreAligned(re0, peven);
        Vector128.StoreAligned(ro0, podd);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private unsafe void ButterflyOne(double* peven, double* podd)
    {
        var even0 = Vector128.LoadAligned(peven);
        var odd0 = Vector128.LoadAligned(podd);
        var re0 = Sse2.Add(even0, odd0);
        var ro0 = Sse2.Subtract(even0, odd0);
        Vector128.StoreAligned(re0, peven);
        Vector128.StoreAligned(ro0, podd);
    }
#else
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private unsafe void Butterfly(Complex* peven, Complex* podd, Vector128<double> w)
    {
        var even0 = Sse2.LoadAlignedVector128((double*)peven);
        var odd0 = Vectorized.ComplexMulSse2(w, Sse2.LoadAlignedVector128((double*)podd));
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

    //[MethodImpl(MethodImplOptions.AggressiveInlining)]
    //private unsafe void Butterfly(Complex* peven, Complex* podd, Vector128<double> w)
    //{
    //    var even0 = Sse2.LoadVector128((double*)peven);
    //    var odd0 = Vectorized.ComplexMulSse2(w, Sse2.LoadVector128((double*)podd));
    //    var re0 = Sse2.Add(even0, odd0);
    //    var ro0 = Sse2.Subtract(even0, odd0);
    //    Sse2.Store((double*)peven, re0);
    //    Sse2.Store((double*)podd, ro0);
    //}

    //[MethodImpl(MethodImplOptions.AggressiveInlining)]
    //private unsafe void ButterflyOne(Complex* peven, Complex* podd)
    //{
    //    var even0 = Sse2.LoadVector128((double*)peven);
    //    var odd0 = Sse2.LoadVector128((double*)podd);
    //    var re0 = Sse2.Add(even0, odd0);
    //    var ro0 = Sse2.Subtract(even0, odd0);
    //    Sse2.Store((double*)peven, re0);
    //    Sse2.Store((double*)podd, ro0);
    //}
#endif
}
