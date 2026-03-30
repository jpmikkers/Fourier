namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

public static class FFTSimpleVectorizedB
{
    private const int maxPrecomputedW = 16;

    public static unsafe void FastFourierTransform(Span<Complex> data, bool isInverse)
    {
        if (!Sse2.IsSupported || !Sse3.IsSupported)
        {
            throw new NotSupportedException("need sse2 and sse3 support for this fft implementation");
        }

        if (!BitOperations.IsPow2(data.Length))
        {
            throw new ArgumentException("fft not a power of two", nameof(data));
        }

        Reorder.Shuffle(data);

        var butterfliesPerPart = 1;             // a single butterfly does 2 angles, +w and -w (=w+pi radians)
        var nrOfParts = data.Length >> 1;       // so the first layer is len/2 parts of single butterflies
        var rotationLookupIndex = 1;
        var complexOne = Vector128.Create(1.0, 0.0);

        //#if EVENFASTER
        if (nrOfParts > 0)
        {
            // no mul needed in the first combination layer, every w is 1+0i and -1+0i
            for (var p = 0; p < nrOfParts; p++)
            {
                Butterflies.Butterfly(ref data[p << 1], ref data[(p << 1) + 1]);
            }
            butterfliesPerPart <<= 1;
            nrOfParts >>= 1;
            rotationLookupIndex++;
        }
        //#endif

        var wspan = stackalloc Vector128<double>[maxPrecomputedW];
        var vspan = MemoryMarshal.Cast<Complex, Vector128<double>>(data);

        fixed (Vector128<double>* vptr = vspan)
        {
            while (nrOfParts > 0)
            {
                var wr = FFTUtils.GetRotation(rotationLookupIndex, isInverse);
                var vwr = Vector128.Create(wr.Real, wr.Imaginary);

                if (butterfliesPerPart < maxPrecomputedW)
                {
                    for (var p = 0; p < nrOfParts; p++)
                    {
                        Sse2.Store((double*)wspan, complexOne);
                        var w = vwr;

                        for (var a = 1; a < butterfliesPerPart; a++)
                        {
                            Sse2.Store((double*)(wspan + a), w);
                            w = Vectorized.ComplexMulSse2(w, vwr);
                        }

                        var evenindex = p << rotationLookupIndex;
                        var oddindex = evenindex + butterfliesPerPart;
                        var peven = (double*)(vptr + evenindex);
                        var podd = (double*)(vptr + oddindex);

                        var even0 = Sse2.LoadVector128(peven);
                        var odd0 = Sse2.LoadVector128(podd);
                        var re0 = Sse2.Add(even0, odd0);
                        var ro0 = Sse2.Subtract(even0, odd0);

                        Sse2.Store(peven, re0);
                        Sse2.Store(podd, ro0);

                        peven += 2;
                        podd += 2;

                        for (var a = 1; a < butterfliesPerPart; a++)
                        {
                            even0 = Sse2.LoadVector128(peven);
                            odd0 = Vectorized.ComplexMulSse2(Sse2.LoadVector128(podd), Sse2.LoadVector128((double*)(wspan + a)));
                            re0 = Sse2.Add(even0, odd0);
                            ro0 = Sse2.Subtract(even0, odd0);

                            Sse2.Store(peven, re0);
                            Sse2.Store(podd, ro0);

                            peven += 2;
                            podd += 2;
                        }
                    }
                }
                else
                {
                    for (var p = 0; p < nrOfParts; p++)
                    {

                        var evenindex = p << rotationLookupIndex;
                        var oddindex = evenindex + butterfliesPerPart;
                        var peven = (double*)(vptr + evenindex);
                        var podd = (double*)(vptr + oddindex);

                        var even0 = Sse2.LoadVector128(peven);
                        var odd0 = Sse2.LoadVector128(podd);
                        var re0 = Sse2.Add(even0, odd0);
                        var ro0 = Sse2.Subtract(even0, odd0);

                        Sse2.Store(peven, re0);
                        Sse2.Store(podd, ro0);

                        peven += 2;
                        podd += 2;
                        var w = vwr;

                        for (var a = 1; a < butterfliesPerPart; a++)
                        {
                            even0 = Sse2.LoadVector128(peven);
                            odd0 = Vectorized.ComplexMulSse2(Sse2.LoadVector128(podd), w);
                            w = Vectorized.ComplexMulSse2(w, vwr);
                            re0 = Sse2.Add(even0, odd0);
                            ro0 = Sse2.Subtract(even0, odd0);

                            Sse2.Store(peven, re0);
                            Sse2.Store(podd, ro0);

                            peven += 2;
                            podd += 2;
                        }
                    }
                }

                butterfliesPerPart <<= 1;
                nrOfParts >>= 1;
                rotationLookupIndex++;
            }
        }

        if (isInverse) FFTUtils.Scale(data);
    }
}
