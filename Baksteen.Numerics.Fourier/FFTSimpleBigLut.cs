namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;

public class FFTSimpleBigLut
{
    private readonly Complex[] _wtable;

    public FFTSimpleBigLut(int length)
    {
        if (!BitOperations.IsPow2(length))
        {
            throw new ArgumentException("fft not a power of two", nameof(length));
        }

        _wtable = [.. Enumerable.Range(0, length/2)
            .Select(t => Complex.FromPolarCoordinates(1, -(Math.Tau * t) / length))];
    }

    public void FastFourierTransform(Span<Complex> data, bool isInverse)
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

        while (nrOfParts > 0)
        {
            for (var p = 0; p < nrOfParts; p++)
            {
                var evenindex = p << rotationLookupIndex;
                var oddindex = evenindex + butterfliesPerPart;

                Butterflies.Butterfly(ref data[evenindex], ref data[oddindex]);
                evenindex++;
                oddindex++;
                var wi = rotationIndexStep;

                for (var a = 1; a < butterfliesPerPart; a++)
                {
                    Butterflies.Butterfly(ref data[evenindex], ref data[oddindex], _wtable[wi]);
                    evenindex++;
                    oddindex++;
                    wi += rotationIndexStep;
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
