using System;
using System.Numerics;

namespace FFT
{
#if false
https://www.HostMath.com/Show.aspx?Code=D(N%2C%20n_0%2C%20k)%0A%3DD(%5Cfrac%7BN%7D%7B2%7D%2C%20n_0%2C%20k)%0A%2B%5Comega_N%5Ek%20D(%5Cfrac%7BN%7D%7B2%7D%2C%20n_0%2B%5Cfrac%7BN%7D%7B2%7D%2C%20k)%0A%5C%5C%0AD(N%2C%20n_0%2C%20k%20%2B%20%5Cfrac%7BN%7D%7B2%7D)%0A%3DD(%5Cfrac%7BN%7D%7B2%7D%2C%20n_0%2C%20k)%0A-%5Comega_N%5Ek%20D(%5Cfrac%7BN%7D%7B2%7D%2C%20n_0%2B%20%5Cfrac%7BN%7D%7B2%7D%2C%20k)%0A%5C%5C%0A--------------------%0A%5C%5C%0AD(N%2C%20n_0%2C%20k)%0A%3DD(%5Cfrac%7BN%7D%7B2%7D%2C%20n_0%2C%20k)%0A%2B(%5Comega_N%5Ek%20D(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%2B%5Cfrac%7BN%7D%7B2%7D%2C%20k)%0A%2B%5Comega_N%5E%7B3k%7DD(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%2B%5Cfrac%7B3%7D%7B4%7DN%2C%20k))%5C%5C%0A%0AD(N%2C%20n_0%2C%20k%2B%5Cfrac%7BN%7D%7B2%7D)%0A%3DD(%5Cfrac%7BN%7D%7B2%7D%2C%20n_0%2C%20k)%0A-(%5Comega_N%5Ek%20D(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%20%2B%20%5Cfrac%7BN%7D%7B2%7D%2C%20k)%0A%2B%5Comega_N%5E%7B3k%7D%20D(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%20%2B%20%5Cfrac%7B3%7D%7B4%7DN%2C%20k))%5C%5C%0A%0AD(N%2C%20n_0%2C%20k%2B%5Cfrac%7BN%7D%7B4%7D)%0A%3DD(%5Cfrac%7BN%7D%7B2%7D%2C%20n_0%2C%20k%2B%5Cfrac%7BN%7D%7B4%7D)%0A-i(%5Comega_N%5Ek%20D(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%2B%5Cfrac%7BN%7D%7B2%7D%2C%20k)%0A-%5Comega_N%5E%7B3k%7DD(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%2B%5Cfrac%7B3%7D%7B4%7DN%2C%20k))%0A%5C%5C%0A%0AD(N%2C%20n_0%2C%20k%2B%5Cfrac%7B3%7D%7B4%7DN)%0A%3DD(%5Cfrac%7BN%7D%7B2%7D%2C%20n_0%2C%20k%2B%5Cfrac%7BN%7D%7B4%7D)%0A%2Bi(%5Comega_N%5Ek%20D(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%2B%5Cfrac%7BN%7D%7B2%7D%2C%20k)%0A-%5Comega_N%5E%7B3k%7DD(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%2B%5Cfrac%7B3%7D%7B4%7DN%2C%20k))%0A%5C%5C%0A
#endif

    // TODO: Split Radix の実装がうまくいかない。漸化式の変形に問題あり。

    public sealed class SplitRadixFastFourierTransform
        : CooleyTukeyFastFourierTransform
    {
        public SplitRadixFastFourierTransform(Options options = Options.None)
            : base(options)
        {
        }

        protected override void TransformCore(ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<int> reversedIndexTable, Span<Complex> data)
        {
            System.Diagnostics.Debug.Assert(twiddleFactor.Length > 0);
            System.Diagnostics.Debug.Assert(twiddleFactor.Length == reversedIndexTable.Length);
            System.Diagnostics.Debug.Assert(twiddleFactor.Length == data.Length);
            System.Diagnostics.Debug.Assert(int.IsPow2(twiddleFactor.Length) == true);

            var quarterOfBlockSize = 1;
            var halfOfBlockSize = 2;
            var blockSize = 4;
            while (blockSize <= data.Length)
            {
                // N = blockSize

                for (var block = 0; block < data.Length; block += blockSize)
                {
                    for (var index0 = block; index0 < blockSize; ++index0)
                    {
                        var index1 = index0 + quarterOfBlockSize;
                        var index2 = index1 + quarterOfBlockSize;
                        var index3 = index2 + quarterOfBlockSize;

                        // N=blockSize/2, n0=block, k=index0 で基数2のFFT
                        // N=blockSize/4, n0=block+halfOfBlockSize, k=index2 で基数4のFFT
                        // N=blockSize/4, n0=block+halfOfBlockSize+quarterOfBlockSize, k=index3 で基数4のFFT
                    }
                }

                halfOfBlockSize <<= 2;
                quarterOfBlockSize <<= 2;
                blockSize <<= 2;
            }

            if (halfOfBlockSize == data.Length)
            {
                // N = halfOfBlockSize
                for (var index0 = 0; index0 < halfOfBlockSize; ++index0)
                {
                    var w = twiddleFactor[index0];
                    var index1 = index0 + quarterOfBlockSize;
                    var t1 = data[index0];
                    var t2 = w * data[index1];
                    data[index0] = t1 + t2;
                    data[index1] = t1 - t2;
                }
            }
        }
    }
}
