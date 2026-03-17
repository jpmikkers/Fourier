namespace Fourier;

using System.Numerics;

/*
    NWaves license:

    Copyright(c) 2017 Tim

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files(the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/
/// <summary>
/// FFT based on the NWaves FFT (https://github.com/ar1st0crat/NWaves), but merges the forward and reverse transform into one and uses
/// the built in Complex type.
/// </summary>
public class Fft64
{
    /// <summary>
    /// Gets FFT size.
    /// </summary>
    public int Size => _fftSize;
    private readonly int _fftSize;

    /// <summary>
    /// Precomputed rotations
    /// </summary>
    private readonly Complex[] _lookupTable;

    /// <summary>
    /// Constructs FFT transformer with given <paramref name="fftSize"/>. FFT size must be a power of two.
    /// </summary>
    /// <param name="fftSize">FFT size</param>
    public Fft64(int fftSize = 512)
    {
        if (!BitOperations.IsPow2(fftSize))
        {
            throw new ArgumentException("fft not a power of two", nameof(fftSize));
        }

        _fftSize = fftSize;

        var tblSize = (int)Math.Log(fftSize, 2);

        _lookupTable = new Complex[tblSize];

        for (int i = 1, pos = 0; i < _fftSize; i *= 2, pos++)
        {
            _lookupTable[pos] = Complex.FromPolarCoordinates(1.0, -i * Math.Tau / _fftSize);
        }
    }

    /// <summary>
    /// Does Fast Fourier Transform in-place.
    /// </summary>
    /// <param name="re">Array of real parts</param>
    /// <param name="im">Array of imaginary parts</param>
    public void Direct(Complex[] data, bool isInverse)
    {
        var L = _fftSize;
        var M = _fftSize >> 1;
        var lookupIndex = 0;

        while (L >= 2)
        {
            var l = L >> 1;
            var phasor = Complex.One;
            var w = isInverse ? Complex.Conjugate(_lookupTable[lookupIndex]) : _lookupTable[lookupIndex];
            lookupIndex++;
            for (var j = 0; j < l; j++)
            {
                for (var i = j; i < _fftSize; i += L)
                {
                    var p = i + l;

                    var di = data[i];
                    var dp = data[p];

                    data[i] = (di + dp);
                    data[p] = (di - dp) * phasor;
                }
                phasor *= w;
            }
            L >>= 1;
        }

        for (int i = 0, j = 0; i < (_fftSize - 1); i++)
        {
            if (i > j)
            {
                (data[j], data[i]) = (data[i], data[j]);
            }
            var k = M;
            while (j >= k)
            {
                j -= k;
                k >>= 1;
            }
            j += k;
        }

        if (isInverse)
        {
            var scaleFactor = 1.0 / _fftSize;
            for (var i = 0; i < _fftSize; i++)
            {
                data[i] *= scaleFactor;
            }
        }
    }
}
