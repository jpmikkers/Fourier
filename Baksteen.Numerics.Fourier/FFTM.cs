namespace Baksteen.Numerics.Fourier;

using System;
using System.Numerics;

public class FFTM
{
    private const int _cachedStages = 5;

    private static readonly Complex[] _rotations;

    private Complex[][] _wcache_cw;
    private Complex[][] _wcache_ccw;

    public FFTM(int length)
    {
        int butterfliesperpart = 1;
        var rotationLookupIndex = 1;

        _wcache_cw = new Complex[_cachedStages][];
        _wcache_ccw = new Complex[_cachedStages][];

        for (int t = 0; t < _cachedStages; t++)
        {
            _wcache_cw[t] = new Complex[butterfliesperpart];
            _wcache_ccw[t] = new Complex[butterfliesperpart];

            var wr = _rotations[rotationLookupIndex];
            var wri = Complex.Conjugate(wr);

            var w = Complex.One;
            var wi = Complex.One;

            for (var a = 0; a < butterfliesperpart; a++)
            {
                _wcache_cw[t][a] = w;
                _wcache_ccw[t][a] = wi;

                w *= wr;
                wi *= wri;
            }

            rotationLookupIndex++;
            butterfliesperpart <<= 1;
        }
    }

    static FFTM()
    {
        // contains 2π/2^n for n=0-31, so π/1, π/2, π/4, π/8 etc..
        _rotations = [.. Enumerable.Range(0, 32)
            .Select(lg2 => Complex.FromPolarCoordinates(1, -Math.Tau / Math.Pow(2.0, lg2)))];
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
        var stage = 0;

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
            stage++;
        }

        while (nrOfParts > 0)
        {
            if (stage < _cachedStages)
            {
                Complex[] lut = isInverse ? _wcache_ccw[stage] : _wcache_cw[stage];
                var evenindex = 0;
                var oddindex = evenindex + butterfliesPerPart;

                for (var a = 0; a < butterfliesPerPart; a++)
                {
                    var w = lut[a];
                    var rli = 0;
                    for (var p = 0; p < nrOfParts; p++)
                    {
                        Butterflies.Butterfly(ref data[evenindex + rli], ref data[oddindex + rli], w);
                        rli += (1 << rotationLookupIndex);
                    }
                    evenindex++;
                    oddindex++;
                }
            }
            else
            {
                var wr = isInverse ? Complex.Conjugate(_rotations[rotationLookupIndex]) : _rotations[rotationLookupIndex];

                for (var p = 0; p < nrOfParts; p++)
                {
                    var evenindex = p << rotationLookupIndex;
                    var oddindex = evenindex + butterfliesPerPart;

                    Butterflies.Butterfly(ref data[evenindex], ref data[oddindex]);
                    var w = wr;
                    evenindex++;
                    oddindex++;

                    for (var a = 1; a < butterfliesPerPart; a++)
                    {
                        Butterflies.Butterfly(ref data[evenindex], ref data[oddindex], w);
                        w *= wr;
                        evenindex++;
                        oddindex++;
                    }
                }
            }

            butterfliesPerPart <<= 1;
            nrOfParts >>= 1;
            rotationLookupIndex++;
            stage++;
        }

        if (isInverse) FFTUtils.Scale(data);
    }
}
