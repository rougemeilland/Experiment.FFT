using System;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace FFT
{
    public sealed class Base4FastFourierTransform
        : CooleyTukeyFastFourierTransform
    {
        public Base4FastFourierTransform(Options options = Options.None)
            : base(options)
        {
        }

        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        protected override void TransformCore(ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<int> reversedIndexTable, Span<Complex> data)
        {
#if false
https://www.hostmath.com/Show.aspx?Code=D_N(k)%3D%5Csum_0%5E%7BN-1%7D%7BS%27(n)%5Comega_N%5E%7Bk%7BR_N(n)%7D%7D%7D%5C%5C%0A%0AD_N(k)%0A%3DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k)%0A%2B%5Comega%5E%7B2k%7DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k%2B%5Cfrac%7B1%7D%7B4%7DN)%20%0A%2B%5Comega%5E%7Bk%7DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k%2B%5Cfrac%7B2%7D%7B4%7DN)%20%0A%2B%5Comega%5E%7B3k%7DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k%2B%5Cfrac%7B3%7D%7B4%7DN)%20%5C%5C%0A%0AD_N(k%20%2B%20%5Cfrac%7B1%7D%7B4%7DN)%0A%3DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k)%0A-%5Comega%5E%7B2k%7DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k%2B%5Cfrac%7B1%7D%7B4%7DN)%20%0A-i%5Comega%5E%7Bk%7DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k%2B%5Cfrac%7B2%7D%7B4%7DN)%20%0A%2Bi%5Comega%5E%7B3k%7DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k%2B%5Cfrac%7B3%7D%7B4%7DN)%20%5C%5C%0A%0AD_N(k%20%2B%20%5Cfrac%7B2%7D%7B4%7DN)%0A%3DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k)%0A%2B%5Comega%5E%7B2k%7DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k%2B%5Cfrac%7B1%7D%7B4%7DN)%20%0A-%5Comega%5E%7Bk%7DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k%2B%5Cfrac%7B2%7D%7B4%7DN)%20%0A-%5Comega%5E%7B3k%7DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k%2B%5Cfrac%7B3%7D%7B4%7DN)%20%5C%5C%0A%0AD_N(k%20%2B%20%5Cfrac%7B3%7D%7B4%7DN)%0A%3DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k)%0A-%5Comega%5E%7B2k%7DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k%2B%5Cfrac%7B1%7D%7B4%7DN)%20%0A%2Bi%5Comega%5E%7Bk%7DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k%2B%5Cfrac%7B2%7D%7B4%7DN)%20%0A-i%5Comega%5E%7B3k%7DD_%7B%5Cfrac%7BN%7D%7B4%7D%7D(k%2B%5Cfrac%7B3%7D%7B4%7DN)%20%5C%5C%0A%0A
#endif
            System.Diagnostics.Debug.Assert(twiddleFactor.Length >= 4);
            System.Diagnostics.Debug.Assert(twiddleFactor.Length == reversedIndexTable.Length);
            System.Diagnostics.Debug.Assert(twiddleFactor.Length == data.Length);
            System.Diagnostics.Debug.Assert(int.IsPow2(twiddleFactor.Length) == true);

            var twiddleFactorIndexMask = twiddleFactor.Length - 1;
            var twiddleFactorIndexStep = twiddleFactor.Length >> 2;
            var quarterOfBlockSize = 1;
            var blockSize = quarterOfBlockSize << 2;
            while (blockSize <= data.Length)
            {
                for (var block = 0; block < data.Length; block += blockSize)
                {
                    var twiddleFactorIndex = 0;
                    for (var index0 = block; index0 < block + quarterOfBlockSize; ++index0)
                    {
                        var twiddleFactorIndex1 = twiddleFactorIndex;
                        var twiddleFactorIndex2 = twiddleFactorIndex1 + twiddleFactorIndex;
                        var twiddleFactorIndex3 = twiddleFactorIndex2 + twiddleFactorIndex;
                        var index1 = index0 + quarterOfBlockSize;
                        var index2 = index1 + quarterOfBlockSize;
                        var index3 = index2 + quarterOfBlockSize;
                        var t00 = data[index0];
                        var t01 = twiddleFactor[twiddleFactorIndex2] * data[index1];
                        var t02 = twiddleFactor[twiddleFactorIndex1] * data[index2];
                        var t03 = twiddleFactor[twiddleFactorIndex3] * data[index3];
                        var t10 = t00 + t01;
                        var t11 = t00 - t01;
                        var t12 = t02 + t03;
                        var t13 = t02 - t03;
                        data[index0] = t10 + t12;
                        data[index1] = new Complex(t11.Real + t13.Imaginary, t11.Imaginary - t13.Real); // t11 - i * t13
                        data[index2] = t10 - t12;
                        data[index3] = new Complex(t11.Real - t13.Imaginary, t11.Imaginary + t13.Real); // t11 + i * t13
                        twiddleFactorIndex += twiddleFactorIndexStep;
                    }
                }

                twiddleFactorIndexStep >>= 2;
                quarterOfBlockSize = blockSize;
                blockSize <<= 2;
            }

            if (quarterOfBlockSize < data.Length)
            {
                System.Diagnostics.Debug.Assert(quarterOfBlockSize == data.Length >> 1);

                var twiddleFactorIndex = 0;
                for (var index0 = 0; index0 < quarterOfBlockSize; ++index0)
                {
                    var w = twiddleFactor[twiddleFactorIndex & twiddleFactorIndexMask];
                    var index1 = index0 + quarterOfBlockSize;
                    var t0 = data[index0];
                    var t1 = w * data[index1];
                    data[index0] = t0 + t1;
                    data[index1] = t0 - t1;
                    ++twiddleFactorIndex;
                }
            }
        }
    }
}
