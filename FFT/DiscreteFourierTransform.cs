using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace FFT
{
    public abstract class DiscreteFourierTransform
        : IDiscreteFourierTransform
    {
        [Flags]
        public enum Options
        {
            None = 0b0,
            OptimizeTwiddleFactorMultiplication = 0b1,
        }

        private const double _2Pi = double.Pi * 2;
        private const double _negative2Pi = -double.Pi * 2;

        private static readonly double _sinOneThirdOfPi = double.Sqrt(3) / 2; // sin(π/3)
        private static readonly Dictionary<int, TwiddleFactor[]> _twiddleFactors1 = [];
        private static readonly Dictionary<int, TwiddleFactor[]> _twiddleFactors2 = [];
        
        private readonly Options _options;

        protected DiscreteFourierTransform(Options options)
        {
            _options = options;
        }

        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
        public void Transform(ReadOnlySpan<Complex> source, Span<Complex> destination)
        {
            ArgumentOutOfRangeException.ThrowIfZero(source.Length);
            ArgumentOutOfRangeException.ThrowIfNotEqual(destination.Length, source.Length);

            switch (source.Length)
            {
                case 1:
                    destination[0] = source[0];
                    break;
                case 2:
                    destination[0] = source[0] + source[1];
                    destination[1] = source[0] - source[1];
                    break;
                case 3:
                {
                    var t1 = source[1] + source[2];
                    var t2 = source[0] - t1 / 2;
                    var t3 = _sinOneThirdOfPi * (source[1] - source[2]);
                    destination[0] = source[0] + t1;
                    destination[1] = new Complex(t2.Real + t3.Imaginary, t2.Imaginary - t3.Real); // t2 - i * t3
                    destination[2] = new Complex(t2.Real - t3.Imaginary, t2.Imaginary + t3.Real); // t2 + i * t3
                    break;
                }
#if false
                case 4:
                {
                    var t1 = source[0] + source[2];
                    var t2 = source[1] + source[3];
                    var t3 = source[0] - source[2];
                    var t4 = source[1] - source[3];
                    destination[0] = t1 + t2;
                    destination[1] = new Complex(t3.Real + t4.Imaginary, t3.Imaginary - t4.Real); // t3 - i * t4
                    destination[2] = t1 - t2;
                    destination[3] = new Complex(t3.Real - t4.Imaginary, t3.Imaginary + t4.Real); // t3 + i * t4
                    break;
                }
#endif
                default:
                    TransformCore(GetTwiddleFactor(false, source.Length), source, destination);
                    break;
            }
        }

        public void InverseTransform(ReadOnlySpan<Complex> source, Span<Complex> destination)
        {
            ArgumentOutOfRangeException.ThrowIfZero(source.Length);
            ArgumentOutOfRangeException.ThrowIfNotEqual(destination.Length, source.Length);

            switch (source.Length)
            {
                case 1:
                    break;
                case 2:
                    break;
                case 3:
                    break;
                case 4:
                    break;
                default:
                    InverseTransformCore(GetTwiddleFactor(true, source.Length), source, destination);
                    break;
            }
        }

        protected abstract void TransformCore(ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<Complex> source, Span<Complex> destination);
        protected abstract void InverseTransformCore(ReadOnlySpan<TwiddleFactor> twiddleFactor, ReadOnlySpan<Complex> source, Span<Complex> destination);

        protected ReadOnlySpan<TwiddleFactor> GetTwiddleFactor([ConstantExpected] bool inverse, int length)
        {
            System.Diagnostics.Debug.Assert(length > 0);

            var twiddleFactors = inverse ? _twiddleFactors1 : _twiddleFactors2;
            lock (twiddleFactors)
            {
                if (twiddleFactors.TryGetValue(length, out var newFactors))
                    return newFactors;

                var optimize = (_options & Options.OptimizeTwiddleFactorMultiplication) != Options.None;
                newFactors = new TwiddleFactor[length];
                for (var index = 0; index < newFactors.Length; ++index)
                    newFactors[index] = new TwiddleFactor(newFactors.Length, index, inverse, optimize);

#if DEBUG
                SelfTest(newFactors, (inverse ? _2Pi : _negative2Pi) / length);
#endif
                twiddleFactors.Add(length, newFactors);
                return newFactors;
            }

            static void SelfTest(TwiddleFactor[] factors, double theta)
            {
                for (var index = 0; index < factors.Length; ++index)
                {
                    var radian = theta * index;
                    var expected = new Complex(double.Cos(radian), double.Sin(radian));
                    var factor = factors[index];
                    var delta = (factor * Complex.One - expected).Magnitude;
                    if (delta > 1e-14)
                        throw new Exception();
                }
            }
        }
    }
}
