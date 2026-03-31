namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

public class FFTSimpleVectorizedE
{
    private AlignedMemoryManager<Complex> _alignedMemoryManager;
    private Memory<Complex> _wtable;

    public FFTSimpleVectorizedE(int length)
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

        Reorder.Shuffle(data);

        var butterfliesPerPart = 1;             // a single butterfly does 2 angles, +w and -w (=w+pi radians)
        var nrOfParts = data.Length >> 1;       // so the first layer is len/2 parts of single butterflies
        var rotationLookupIndex = 1;
        var rotationIndexStep = data.Length >> 1;

        var vspan = MemoryMarshal.Cast<Complex, Vector128<double>>(data);
        var wspan = MemoryMarshal.Cast<Complex, Vector128<double>>(_wtable.Span);

        fixed (Vector128<double>* vptr = vspan)
        fixed (Vector128<double>* wptr = wspan)
        {
            if (nrOfParts > 0)
            {
                var ptmp = (double*)vptr;

                // no mul needed in the first combination layer, every w is 1+0i and -1+0i
                for (var p = 0; p < nrOfParts; p++)
                {
                    var even0 = Sse2.LoadAlignedVector128(ptmp);
                    var odd0 = Sse2.LoadAlignedVector128(ptmp + 2);
                    var re0 = Sse2.Add(even0, odd0);
                    var ro0 = Sse2.Subtract(even0, odd0);
                    Sse2.StoreAligned(ptmp, re0);
                    Sse2.StoreAligned(ptmp + 2, ro0);
                    ptmp += 4;
                }

                butterfliesPerPart <<= 1;
                nrOfParts >>= 1;
                rotationLookupIndex++;
                rotationIndexStep >>= 1;
            }

            while (nrOfParts > 0)
            {
                for (var p = 0; p < nrOfParts; p++)
                {
                    var evenindex = p << rotationLookupIndex;
                    var oddindex = evenindex + butterfliesPerPart;
                    var peven = (double*)(vptr + evenindex);
                    var podd = (double*)(vptr + oddindex);
                    var pw = (double*)wptr;

                    var even0 = Sse2.LoadAlignedVector128(peven);
                    var odd0 = Sse2.LoadAlignedVector128(podd);
                    var re0 = Sse2.Add(even0, odd0);
                    var ro0 = Sse2.Subtract(even0, odd0);
                    Sse2.StoreAligned(peven, re0);
                    Sse2.StoreAligned(podd, ro0);

                    peven += 2;
                    podd += 2;
                    pw += (rotationIndexStep << 1);

                    for (var a = 1; a < butterfliesPerPart; a++)
                    {
                        // using LoadAligned is about 3% faster than Load
                        even0 = Sse2.LoadAlignedVector128(peven);
                        odd0 = Vectorized.ComplexMulSse2(Sse2.LoadAlignedVector128(podd), Sse2.LoadAlignedVector128(pw));

                        re0 = Sse2.Add(even0, odd0);
                        ro0 = Sse2.Subtract(even0, odd0);
                        Sse2.StoreAligned(peven, re0);
                        Sse2.StoreAligned(podd, ro0);

                        peven += 2;
                        podd += 2;
                        pw += (rotationIndexStep << 1);
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
}
