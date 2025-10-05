using System;
using System.Collections.Generic;
using System.Numerics;

namespace FFT
{
    public abstract class CooleyTukeyFastFourierTransform
        : DiscreteFourierTransform
    {
        private static readonly Dictionary<int, int[]> _reversedIndexTables = [];

        protected CooleyTukeyFastFourierTransform(Options options)
            : base(options)
        {
        }

        protected override void TransformCore(ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<Complex> source, Span<Complex> destination)
        {
            if (!int.IsPow2(source.Length))
                throw new InvalidOperationException();

            System.Diagnostics.Debug.Assert(source.Length > 0);
            System.Diagnostics.Debug.Assert(source.Length == destination.Length);

            var reversedIndexTable = GetReversedIndexTable(source.Length);
            for (var index = 0; index < destination.Length; ++index)
                destination[index] = source[reversedIndexTable[index]];
            TransformCore(twiddleFactor, reversedIndexTable, destination);
        }

        protected override void InverseTransformCore(ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<Complex> source, Span<Complex> destination)
        {
            if (!int.IsPow2(source.Length))
                throw new InvalidOperationException();

            System.Diagnostics.Debug.Assert(source.Length > 0);
            System.Diagnostics.Debug.Assert(source.Length == destination.Length);

            var reversedIndexTable = GetReversedIndexTable(source.Length);
            for (var index = 0; index < destination.Length; ++index)
                destination[index] = source[reversedIndexTable[index]] / source.Length;
            TransformCore(twiddleFactor, reversedIndexTable, destination);
        }

        protected abstract void TransformCore(ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<int> reversedIndexTable, Span<Complex> data);

        private static int[] GetReversedIndexTable(int length)
        {
            System.Diagnostics.Debug.Assert(int.IsPow2(length) == true);

            lock (_reversedIndexTables)
            {
                if (_reversedIndexTables.TryGetValue(length, out var reversedIndexTable))
                    return reversedIndexTable;
                reversedIndexTable = new int[length];
                if (length < (1 << 8))
                {
                    for (var index = 0; index < length; ++index)
                    {
                        var reversedIndex = ((0x55555555u & (uint)index) << 1) | ((0xaaaaaaaau & (uint)index) >> 1);
                        reversedIndex = ((0x33333333u & reversedIndex) << 2) | ((0xccccccccu & reversedIndex) >> 2);
                        reversedIndex = ((0x0f0f0f0fu & reversedIndex) << 4) | ((0xf0f0f0f0u & reversedIndex) >> 4);
                        reversedIndexTable[index] = (int)(reversedIndex >> (8 - (int)uint.TrailingZeroCount((uint)length)));
                    }
                }
                else if (length < (1 << 16))
                {
                    for (var index = 0; index < length; ++index)
                    {
                        var reversedIndex = ((0x55555555u & (uint)index) << 1) | ((0xaaaaaaaau & (uint)index) >> 1);
                        reversedIndex = ((0x33333333u & reversedIndex) << 2) | ((0xccccccccu & reversedIndex) >> 2);
                        reversedIndex = ((0x0f0f0f0fu & reversedIndex) << 4) | ((0xf0f0f0f0u & reversedIndex) >> 4);
                        reversedIndex = ((0x00ff00ffu & reversedIndex) << 8) | ((0xff00ff00u & reversedIndex) >> 8);
                        reversedIndexTable[index] = (int)(reversedIndex >> (16 - (int)uint.TrailingZeroCount((uint)length)));
                    }
                }
                else
                {
                    for (var index = 0; index < length; ++index)
                    {
                        var reversedIndex = ((0x55555555u & (uint)index) << 1) | ((0xaaaaaaaau & (uint)index) >> 1);
                        reversedIndex = ((0x33333333u & reversedIndex) << 2) | ((0xccccccccu & reversedIndex) >> 2);
                        reversedIndex = ((0x0f0f0f0fu & reversedIndex) << 4) | ((0xf0f0f0f0u & reversedIndex) >> 4);
                        reversedIndex = ((0x00ff00ffu & reversedIndex) << 8) | ((0xff00ff00u & reversedIndex) >> 8);
                        reversedIndex = ((0x0000ffffu & reversedIndex) << 16) | ((0xffff0000u & reversedIndex) >> 16);
                        reversedIndexTable[index] = (int)(reversedIndex >> (32 - (int)uint.TrailingZeroCount((uint)length)));
                    }
                }

                _reversedIndexTables.Add(length, reversedIndexTable);
                return reversedIndexTable;
            }
        }
    }
}
