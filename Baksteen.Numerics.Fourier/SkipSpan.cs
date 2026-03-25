namespace Baksteen.Numerics.Fourier;

using System;
using System.Collections;
using System.Runtime.CompilerServices;

/// <summary>
/// Span wrapper that allows slicing a span into odd and even indexed elements.
/// </summary>
/// <typeparam name="T"></typeparam>
// see https://github.com/dotnet/runtime/blob/main/src/libraries/System.Private.CoreLib/src/System/Span.cs
public readonly ref struct SkipSpan<T>
{
    public readonly Span<T> _original;
    public readonly int _bitshift;

    public SkipSpan(Span<T> original)
    {
        _original = original;
        _bitshift = 0;
    }

    private SkipSpan(Span<T> original, int bitshift)
    {
        _original = original;
        _bitshift = bitshift;
    }

    //
    // Summary:
    //     Gets the element at the specified zero-based index.
    //
    // Parameters:
    //   index:
    //     The zero-based index of the element.
    //
    // Returns:
    //     The element at the specified index.
    //
    // Exceptions:
    //   T:System.IndexOutOfRangeException:
    //     index is less than zero or greater than or equal to System.Span`1.Length.
    public ref T this[int index]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            return ref _original[index << _bitshift];
        }
    }

    /// <summary>
    /// The number of items in the span.
    /// </summary>
    public int Length
    {
        get => (_original.Length + ((1 << _bitshift) - 1)) >> _bitshift;
    }

    public bool IsEmpty
    {
        get => Length == 0;
    }

    public T[] ToArray()
    {
        if (IsEmpty)
        {
            return [];
        }
        var destination = new T[Length];
        for (var t = 0; t < Length; t++)
        {
            destination[t] = _original[t << _bitshift];
        }
        return destination;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public SkipSpan<T> SliceOdds()
    {
        return new SkipSpan<T>(_original[(1 << _bitshift)..], _bitshift + 1);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public SkipSpan<T> SliceEvens()
    {
        return new SkipSpan<T>(_original, _bitshift + 1);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public SkipSpan<T> Slice(int start)
    {
        return new SkipSpan<T>(_original.Slice(start << _bitshift), _bitshift);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public SkipSpan<T> Slice(int start, int length)
    {
        //var s = start << _bitshift;
        //var l = ((start + length) << _bitshift) - s;
        return new SkipSpan<T>(_original.Slice(start << _bitshift, length << _bitshift), _bitshift);
    }

    /// <summary>Gets an enumerator for this span.</summary>
    public Enumerator GetEnumerator() => new Enumerator(this);

    /// <summary>Enumerates the elements of a <see cref="Span{T}"/>.</summary>
    public ref struct Enumerator : IEnumerator<T>
    {
        /// <summary>The span being enumerated.</summary>
        private readonly SkipSpan<T> _span;
        /// <summary>The next index to yield.</summary>
        private int _index;

        /// <summary>Initialize the enumerator.</summary>
        /// <param name="span">The span to enumerate.</param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal Enumerator(SkipSpan<T> span)
        {
            _span = span;
            _index = -1;
        }

        /// <summary>Advances the enumerator to the next element of the span.</summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool MoveNext()
        {
            int index = _index + 1;
            if (index < _span.Length)
            {
                _index = index;
                return true;
            }

            return false;
        }

        /// <summary>Gets the element at the current position of the enumerator.</summary>
        public ref T Current
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => ref _span[_index];
        }

        /// <inheritdoc />
        T IEnumerator<T>.Current => Current;

        /// <inheritdoc />
        object IEnumerator.Current => Current!;

        /// <inheritdoc />
        void IEnumerator.Reset() => _index = -1;

        /// <inheritdoc />
        void IDisposable.Dispose() { }
    }
}