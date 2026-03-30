namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

public static class FFTSimpleVectorized
{
    public static void FastFourierTransform(Span<Complex> data, bool isInverse)
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

        var vspan = MemoryMarshal.Cast<Complex, Vector128<double>>(data);

        while (nrOfParts > 0)
        {
            var wr = FFTUtils.GetRotation(rotationLookupIndex, isInverse);
            var vwr = Vector128.Create(wr.Real, wr.Imaginary);

            for (var p = 0; p < nrOfParts; p++)
            {
                var w = complexOne;
                var evenindex = p << rotationLookupIndex;
                var oddindex = evenindex + butterfliesPerPart;

                for (var a = 0; a < butterfliesPerPart; a++)
                {
                    var even = vspan[evenindex];
                    var odd_w = Vectorized.ComplexMulSse2(vspan[oddindex], w);
                    vspan[evenindex] = Sse2.Add(even, odd_w);
                    vspan[oddindex] = Sse2.Subtract(even, odd_w);
                    w = Vectorized.ComplexMulSse2(w, vwr);
                    evenindex++;
                    oddindex++;
                }
            }

            butterfliesPerPart <<= 1;
            nrOfParts >>= 1;
            rotationLookupIndex++;
        }

        if (isInverse) FFTUtils.Scale(data);
    }
}
