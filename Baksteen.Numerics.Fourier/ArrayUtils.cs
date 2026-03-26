namespace Baksteen.Numerics.Fourier;

public static class ArrayUtils
{
    /// <summary>
    /// Reorders the array in-place using only swaps so that:
    /// - All originally even-indexed elements (0-based) appear in the first ⌈n/2⌉ positions,
    ///   preserving their original relative order.
    /// - All originally odd-indexed elements appear in the remaining positions,
    ///   preserving their original relative order.
    /// 
    /// The algorithm runs in O(n) swaps (each O(1)) and O(n) total time in practice
    /// (the inner while loop amortizes to a very small constant factor for this specific
    /// permutation; worst-case chain length is O(log n)).
    /// No extra space beyond a few integer variables is used.
    /// </summary>
    public static void ReorderEvenOddRelative<T>(Span<T> array)
    {
        if (array.Length < 2) return;

        var n = array.Length;
        var m = (n + 1) / 2;   // ⌈n/2⌉ – size of the first half

        // Helper to compute the original source index for target position i
        // (this defines the exact permutation that preserves relative order)
        static int Source(int i, int n, int m)
        {
            return i < m ? 2 * i : 2 * (i - m) + 1;
        }

        for (int i = 1; i < n; i++)
        {
            int j = Source(i, n, m);

            // Follow the permutation chain backwards until we reach a position
            // that has not yet been finalized (j >= i). This guarantees we
            // always pull an untouched original value.
            while (j < i)
            {
                j = Source(j, n, m);
            }

            // Perform the swap (self-swap is harmless and never occurs here)
            if (j != i)
            {
                (array[i], array[j]) = (array[j], array[i]);
            }
        }
    }

    public static void ReorderEvenOddRelative2<T>(Span<T> array)
    {
        if (array.Length < 2) return;

        // Helper to compute the original source index for target position i
        // (this defines the exact permutation that preserves relative order)
        static int Source(int i, int n)
        {
            var m = (n + 1) >> 1;   // ⌈n/2⌉ – size of the first half
            var j = i;
            do
            {
                j = j < m ? (j << 1) : ((j - m) << 1) + 1;
            }
            while (j < i);
            return j;
        }

        for (var i = 1; i < array.Length; i++)
        {
            var j = Source(i, array.Length);

            // Perform the swap (self-swap is harmless and never occurs here)
            if (j != i)
            {
                (array[i], array[j]) = (array[j], array[i]);
            }
        }
    }
}

#if NEVER
// -----------------------------
// Example usage & verification
// -----------------------------
class Program
{
    static void Main()
    {
        // Test case 1: even length (n=6)
        int[] even = { 0, 1, 2, 3, 4, 5 };           // values = original indices
        ArrayUtils.ReorderEvenOddRelative(even);
        Console.WriteLine("n=6 result : " + string.Join(" ", even));
        // Expected: 0 2 4 1 3 5   (evens first in order, odds second in order)

        // Test case 2: odd length (n=5)
        int[] odd = { 0, 1, 2, 3, 4 };
        ArrayUtils.ReorderEvenOddRelative(odd);
        Console.WriteLine("n=5 result : " + string.Join(" ", odd));
        // Expected: 0 2 4 1 3

        // Test case 3: small edge case
        int[] small = { 10, 20, 30 };
        ArrayUtils.ReorderEvenOddRelative(small);
        Console.WriteLine("n=3 result : " + string.Join(" ", small));
        // Expected: 10 30 20
    }
}
#endif