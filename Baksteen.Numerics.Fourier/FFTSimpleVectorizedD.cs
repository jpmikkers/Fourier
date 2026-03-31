namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

public class FFTSimpleVectorizedD
{
    private Complex[] _wtable;

    public FFTSimpleVectorizedD(int length)
    {
        if (!Sse2.IsSupported || !Sse3.IsSupported)
        {
            throw new NotSupportedException("need sse2 and sse3 support for this fft implementation");
        }

        if (!BitOperations.IsPow2(length))
        {
            throw new ArgumentException("fft not a power of two", nameof(length));
        }

        _wtable = [.. Enumerable.Range(0, length/2)
            .Select(t => Complex.FromPolarCoordinates(1, -(Math.Tau * t) / length))];
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
            rotationIndexStep >>= 1;
        }
        //#endif

        var vspan = MemoryMarshal.Cast<Complex, Vector128<double>>(data);
        var wspan = MemoryMarshal.Cast<Complex, Vector128<double>>(_wtable);

        while (nrOfParts > 0)
        {
            for (var p = 0; p < nrOfParts; p++)
            {
                var evenindex = p << rotationLookupIndex;
                var oddindex = evenindex + butterfliesPerPart;
                var widx = 0;

                for (var a = 0; a < butterfliesPerPart; a++)
                {
                    // this is slower than sse2.load
                    var even0 = vspan[evenindex];
                    var odd0 = Vectorized.ComplexMulSse2(vspan[oddindex], wspan[widx]);
                    vspan[evenindex] = Sse2.Add(even0, odd0);
                    vspan[oddindex] = Sse2.Subtract(even0, odd0);
                    evenindex++;
                    oddindex++;
                    widx += rotationIndexStep;
                }
            }

            butterfliesPerPart <<= 1;
            nrOfParts >>= 1;
            rotationLookupIndex++;
            rotationIndexStep >>= 1;
        }

        if (isInverse) FFTUtils.Scale(data);
    }
}
