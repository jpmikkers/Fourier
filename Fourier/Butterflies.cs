namespace Fourier;

using System.Numerics;
using System.Runtime.CompilerServices;

internal static class Butterflies
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void Butterfly(ref Complex even, ref Complex odd, Complex w)
    {
        var odd_w = odd * w;
        (even, odd) = (even + odd_w, even - odd_w);
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
}
