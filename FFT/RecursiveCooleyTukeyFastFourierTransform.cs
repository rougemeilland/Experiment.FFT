using System;
using System.Buffers;
using System.Diagnostics.CodeAnalysis;
using System.Numerics;

namespace FFT
{
    public class RecursiveCooleyTukeyFastFourierTransform
        : DiscreteFourierTransform
    {
        public RecursiveCooleyTukeyFastFourierTransform(Options options = Options.None)
            : base(options)
        {
        }

        protected override void TransformCore(ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<Complex> source, Span<Complex> destination)
        {
            if (!int.IsPow2(source.Length))
                throw new InvalidOperationException();

            System.Diagnostics.Debug.Assert(source.Length > 0);
            System.Diagnostics.Debug.Assert(source.Length == destination.Length);

            TransformCore(false, twiddleFactor, source, destination);
        }

        protected override void InverseTransformCore(ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<Complex> source, Span<Complex> destination)
        {
            if (!int.IsPow2(source.Length))
                throw new InvalidOperationException();

            System.Diagnostics.Debug.Assert(source.Length > 0);
            System.Diagnostics.Debug.Assert(source.Length == destination.Length);

            TransformCore(true, twiddleFactor, source, destination);
        }

        private static void TransformCore([ConstantExpected] bool inverse, ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<Complex> source, Span<Complex> destination)
        {
            System.Diagnostics.Debug.Assert(twiddleFactor.Length > 0);
            System.Diagnostics.Debug.Assert(twiddleFactor.Length == source.Length);
            System.Diagnostics.Debug.Assert(twiddleFactor.Length == destination.Length);
            System.Diagnostics.Debug.Assert(int.IsPow2(twiddleFactor.Length) == true);

            TransformCore2(twiddleFactor, source, destination);
            if (inverse)
            {
                for (var index = 0; index < destination.Length; ++index)
                    destination[index] /= destination.Length;
            }
        }

        private static void TransformCore2(ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<Complex> source, Span<Complex> destination)
        {
            System.Diagnostics.Debug.Assert(twiddleFactor.Length > 0);
            System.Diagnostics.Debug.Assert(twiddleFactor.Length == source.Length);
            System.Diagnostics.Debug.Assert(twiddleFactor.Length == destination.Length);
            System.Diagnostics.Debug.Assert(int.IsPow2(twiddleFactor.Length) == true);

            if (source.Length == 1)
            {
                destination[0] = source[0];
                return;
            }

            var tempFactor = (TwiddleFactor[]?)null;
            var tempSource1 = (Complex[]?)null;
            var tempSource2 = (Complex[]?)null;
            var tempDestination1 = (Complex[]?)null;
            var tempDestination2 = (Complex[]?)null;
            try
            {
                var half = twiddleFactor.Length >> 1;
                tempFactor = ArrayPool<TwiddleFactor>.Shared.Rent(half);
                tempSource1 = ArrayPool<Complex>.Shared.Rent(half);
                tempSource2 = ArrayPool<Complex>.Shared.Rent(half);
                tempDestination1 = ArrayPool<Complex>.Shared.Rent(half);
                tempDestination2 = ArrayPool<Complex>.Shared.Rent(half);

                for (var index = 0; index < half; ++index)
                    tempFactor[index] = twiddleFactor[index << 1];
                for (var index = 0; index < half; ++index)
                    tempSource1[index] = source[index << 1];
                for (var index = 0; index < half; ++index)
                    tempSource2[index] = source[(index << 1) + 1];
                TransformCore2(tempFactor.AsSpan(0, half), tempSource1.AsSpan(0, half), tempDestination1.AsSpan(0, half));
                TransformCore2(tempFactor.AsSpan(0, half), tempSource2.AsSpan(0, half), tempDestination2.AsSpan(0, half));
                for (var index = 0; index < half; ++index)
                {
                    var v1 = tempDestination1[index];
                    var v2 = twiddleFactor[index] * tempDestination2[index];
                    destination[index] = v1 + v2;
                    destination[index + (destination.Length >> 1)] = v1 - v2;
                }
            }
            finally
            {
                if (tempFactor is not null)
                    ArrayPool<TwiddleFactor>.Shared.Return(tempFactor);
                if (tempSource1 is not null)
                    ArrayPool<Complex>.Shared.Return(tempSource1);
                if (tempSource2 is not null)
                    ArrayPool<Complex>.Shared.Return(tempSource2);
                if (tempDestination1 is not null)
                    ArrayPool<Complex>.Shared.Return(tempDestination1);
                if (tempDestination2 is not null)
                    ArrayPool<Complex>.Shared.Return(tempDestination2);
            }
        }
    }
}
