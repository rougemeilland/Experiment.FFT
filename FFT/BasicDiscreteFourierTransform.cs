using System;
using System.Diagnostics.CodeAnalysis;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace FFT
{
    public sealed class BasicDiscreteFourierTransform
        : DiscreteFourierTransform
    {
        public BasicDiscreteFourierTransform(Options options = Options.None)
            : base(options)
        {
        }

        protected override void TransformCore(ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<Complex> source, Span<Complex> destination)
        {
            System.Diagnostics.Debug.Assert(source.Length > 0);
            System.Diagnostics.Debug.Assert(source.Length == destination.Length);

            TransformCore(false, GetTwiddleFactor(false, source.Length), source, destination);
        }

        protected override void InverseTransformCore(ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<Complex> source, Span<Complex> destination)
        {
            System.Diagnostics.Debug.Assert(source.Length > 4);
            System.Diagnostics.Debug.Assert(source.Length == destination.Length);

            TransformCore(true, GetTwiddleFactor(true, source.Length), source, destination);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
        private static void TransformCore([ConstantExpected] bool inverse, ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<Complex> source, Span<Complex> destination)
        {
            System.Diagnostics.Debug.Assert(source.Length == destination.Length);

            for (var destinationIndex = 0; destinationIndex < destination.Length; ++destinationIndex)
            {
                var sum = Complex.Zero;
                for (var sourceIndex = 0; sourceIndex < source.Length; ++sourceIndex)
                    sum += source[sourceIndex] * twiddleFactor[(int)((uint)sourceIndex * (uint)destinationIndex % (uint)source.Length)];
                destination[destinationIndex] = inverse ? sum / source.Length : sum;
            }

            destination[source.Length..].Clear();
        }
    }
}
