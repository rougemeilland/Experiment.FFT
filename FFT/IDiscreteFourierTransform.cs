using System;
using System.Numerics;

namespace FFT
{
    internal interface IDiscreteFourierTransform
    {
        void Transform(ReadOnlySpan<Complex> source, Span<Complex> destination);
        void InverseTransform(ReadOnlySpan<Complex> source, Span<Complex> destination);
    }
}
