namespace Baksteen.Numerics.Fourier;

using System.Numerics;
using System.Runtime.CompilerServices;

internal static class Butterflies
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void Butterfly(ref Complex even, ref Complex odd, Complex w)
    {
        var odd_w = w * odd;
        (even, odd) = (even + odd_w, even - odd_w);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void Butterfly2(ref Complex even, ref Complex odd, Complex w)
    {
        var tmpe = even;
        var odd_w = w * odd;
        even = tmpe + odd_w;
        odd = tmpe - odd_w;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void Butterfly(ref Complex even, ref Complex odd)
    {
        (even, odd) = (even + odd, even - odd);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void Butterfly90CW(ref Complex even, ref Complex odd)
    {
        var odd_w = new Complex(odd.Imaginary, -odd.Real);  // 90 degree clockwise
        (even, odd) = (even + odd_w, even - odd_w);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void Butterfly90CCW(ref Complex even, ref Complex odd)
    {
        var odd_w = new Complex(-odd.Imaginary, odd.Real);  // 90 degree counterclockwise
        (even, odd) = (even + odd_w, even - odd_w);
    }

    // Multiply using the 3‑multiply algorithm
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static Complex KaratsubaMultiply(Complex x, Complex y)
    {
        var p = x.Real * y.Real;
        var q = x.Imaginary * y.Imaginary;
        var r = (x.Real + x.Imaginary) * (y.Real + y.Imaginary);
        return new Complex(p - q, r - p - q);
    }
}
