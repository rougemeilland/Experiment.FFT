using System;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace FFT
{
    public sealed class Base2FastFourierTransform
        : CooleyTukeyFastFourierTransform
    {
        public Base2FastFourierTransform(Options options = Options.None)
            : base(options)
        {
        }

        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        protected override void TransformCore(ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<int> reversedIndexTable, Span<Complex> data)
        {
            System.Diagnostics.Debug.Assert(twiddleFactor.Length >= 2);
            System.Diagnostics.Debug.Assert(twiddleFactor.Length == reversedIndexTable.Length);
            System.Diagnostics.Debug.Assert(twiddleFactor.Length == data.Length);
            System.Diagnostics.Debug.Assert(int.IsPow2(twiddleFactor.Length) == true);

#if false

https://www.hostmath.com/Show.aspx?Code=D_N(k)%3D%5Csum_0%5E%7BN-1%7D%7BS%27(n)%5Comega_N%5E%7Bk%7BR_N(n)%7D%7D%7D%5C%5C%0A%0AD_N(k)%0A%3D%5Csum_%7Bn%3D0%7D%5E%7B%5Cfrac%7BN%7D%7B2%7D-1%7D%7BS%27(n)%5Comega_N%5E%7Bk%7BR_N(n)%7D%7D%7D%0A%2B%5Csum_%7Bn%3D%5Cfrac%7BN%7D%7B2%7D%7D%5E%7BN-1%7D%7BS%27(n)%5Comega_N%5E%7Bk%7BR_N(n)%7D%7D%7D%0A%3D%5Csum_%7Bn%3D0%7D%5E%7B%5Cfrac%7BN%7D%7B2%7D-1%7D%7BS%27(n)%5Comega_%7B%5Cfrac%7BN%7D%7B2%7D%7D%5E%7Bk%7BR_%7B%5Cfrac%7BN%7D%7B2%7D%7D(n)%7D%7D%7D%0A%2B%20%5Comega_N%5Ek%20%5Csum_%7Bn%3D0%7D%5E%7B%5Cfrac%7BN%7D%7B2%7D-1%7D%7BS%27(n%20%2B%20%5Cfrac%7BN%7D%7B2%7D)%5Comega_%7B%5Cfrac%7BN%7D%7B2%7D%7D%5E%7Bk%7BR_%7B%5Cfrac%7BN%7D%7B2%7D%7D(n)%7D%7D%7D%0A%3DD_%7B%5Cfrac%7BN%7D%7B2%7D%7D(k)%0A%2B%5Comega_N%5Ek%20D_%7B%5Cfrac%7BN%7D%7B2%7D%7D(k%20%2B%20%5Cfrac%7BN%7D%7B2%7D)%2C%20k%E3%81%AF%E5%81%B6%E6%95%B0%5C%5C%0A%0AD_N(k%20%2B%20%5Cfrac%7BN%7D%7B2%7D)%0A%3D%5Csum_%7Bn%3D0%7D%5E%7B%5Cfrac%7BN%7D%7B2%7D-1%7D%7BS%27(n)%5Comega_N%5E%7B(k%20%2B%20%5Cfrac%7BN%7D%7B2%7D)%7BR_N(n)%7D%7D%7D%0A%2B%5Comega_N%5Ek%20%5Csum_%7Bn%3D%5Cfrac%7BN%7D%7B2%7D%7D%5E%7BN-1%7D%7BS%27(n)%5Comega_N%5E%7B(k%20%2B%20%5Cfrac%7BN%7D%7B2%7D)%7BR_N(n)%7D%7D%7D%0A%0A%3D%5Csum_%7Bn%3D0%7D%5E%7B%5Cfrac%7BN%7D%7B2%7D-1%7D%7BS%27(n)%5Comega_N%5E%7BkR_%7B%5Cfrac%7BN%7D%7B2%7D%7D(n)%7D%7D%0A-%20%5Comega_N%5Ek%20%5Csum_%7Bn%3D%5Cfrac%7BN%7D%7B2%7D%7D%5E%7BN-1%7D%7BS%27(n%20%2B%20%5Cfrac%7BN%7D%7B2%7D)%5Comega_N%5E%7Bk%7BR_%7B%5Cfrac%7BN%7D%7B2%7D%7D(n)%7D%7D%7D%0A%3DD_%7B%5Cfrac%7BN%7D%7B2%7D%7D(k)%0A-%5Comega_N%5Ek%20D_%7B%5Cfrac%7BN%7D%7B2%7D%7D(k%20%2B%20%5Cfrac%7BN%7D%7B2%7D)%2C%20k%E3%81%AF%E5%81%B6%E6%95%B0%5C%5C%0A%0AD_2(k)%3DS%27(k)%2B%20S%27(k%20%2B%201)%2C%20k%E3%81%AF%E5%81%B6%E6%95%B0%5C%5C%0AD_2(k%20%2B%201)%3DS%27(k)-%20S%27(k%20%2B%201)%2C%20k%E3%81%AF%E5%81%B6%E6%95%B0%5C%5C%0A

#endif
            var twiddleFactorIndexMask = twiddleFactor.Length - 1;
            var twiddleFactorIndexStep = twiddleFactor.Length >> 1;
            var halfOfBlockSize = 1;
            var blockSize = halfOfBlockSize << 1;
            while (blockSize <= data.Length)
            {
                for (var block = 0; block < data.Length; block += blockSize)
                {
                    var twiddleFactorIndex = 0;
                    for (var index0 = block; index0 < block + halfOfBlockSize; ++index0)
                    {
                        var w = twiddleFactor[twiddleFactorIndex & twiddleFactorIndexMask];
                        var index1 = index0 + halfOfBlockSize;
                        var t0 = data[index0];
                        var t1 = w * data[index1];
                        data[index0] = t0 + t1;
                        data[index1] = t0 - t1;
                        twiddleFactorIndex += twiddleFactorIndexStep;
                    }
                }

                twiddleFactorIndexStep >>= 1;
                halfOfBlockSize = blockSize;
                blockSize <<= 1;
            }
        }
    }
}
