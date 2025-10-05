using System;
using System.Globalization;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace FFT
{
    public readonly struct TwiddleFactor
    {
        private readonly Func<Complex, Complex> _multiplicator;
        private readonly Complex _complexValue;

        public TwiddleFactor(int N, int n, bool inverse, bool optimize = false)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(N);
            ArgumentOutOfRangeException.ThrowIfNegative(n);
            ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(n, N);

            _complexValue = CreateComplexValue(N, n, inverse);
            _multiplicator = CreateMultiplicator(N, n, inverse, optimize);
            IsRealNumber = n == 0 || n * 2L == N;
#if DEBUG
            {
                var theta = (inverse ? 2d : -2d) * n / N;
                var expected = new Complex(double.CosPi(theta), double.SinPi(theta));
                var actual = this * Complex.One;
                if ((actual - expected).Magnitude > 1e-14)
                    throw new Exception();
            }

            {
                var theta = (inverse ? 2d : -2d) * ((double)n / N + 0.25);
                var expected = new Complex(double.CosPi(theta), double.SinPi(theta));
                var actual = this * (inverse ? Complex.ImaginaryOne : -Complex.ImaginaryOne);
                if ((actual - expected).Magnitude > 1e-14)
                    throw new Exception();
            }
#endif
        }

        public bool IsRealNumber { get; }

        [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
        public Complex Multiply(Complex value) => _multiplicator(value);

        [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
        public static Complex operator *(TwiddleFactor left, Complex right) => left.Multiply(right);

        [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
        public static Complex operator *(Complex left, TwiddleFactor right) => right.Multiply(left);

        public override string ToString() => _complexValue.ToString(CultureInfo.InvariantCulture);

        private static Func<Complex, Complex> CreateMultiplicator(int N, int n, bool inverse, bool optimize)
        {
            System.Diagnostics.Debug.Assert(N > 0);
            System.Diagnostics.Debug.Assert(n >= 0);
            System.Diagnostics.Debug.Assert(n < N);

            if (!inverse && n != 0)
            {
                n = N - n;
                inverse = true;
            }

            if (n == 0)
            {
                return v => v;
            }
            else if (n * 2L == N)
            {
                return v => -v;
            }
            else if (n * 4L == N)
            {
                return v => new Complex(-v.Imaginary, v.Real);
            }
            else if (n * 4L == N * 3L)
            {
                return v => new Complex(v.Imaginary, -v.Real);
            }
            else if (n * 8L == N)
            {
                var sqrt2 = double.Sqrt(2);
                return
                [MethodImpl(MethodImplOptions.AggressiveOptimization)]
                (v) => new Complex((v.Real - v.Imaginary) / sqrt2, (v.Real + v.Imaginary) / sqrt2);
            }
            else if (n * 8L == N * 3L)
            {
                var negativeSqrt2 = -double.Sqrt(2);
                return
                [MethodImpl(MethodImplOptions.AggressiveOptimization)]
                (v) => new Complex((v.Real + v.Imaginary) / negativeSqrt2, (-v.Real + v.Imaginary) / negativeSqrt2);
            }
            else if (n * 8L == N * 5L)
            {
                var negativeSqrt2 = -double.Sqrt(2);
                return
                [MethodImpl(MethodImplOptions.AggressiveOptimization)]
                (v) => new Complex((v.Real - v.Imaginary) / negativeSqrt2, (v.Real + v.Imaginary) / negativeSqrt2);
            }
            else if (n * 8L == N * 7L)
            {
                var sqrt2 = double.Sqrt(2);
                return
                [MethodImpl(MethodImplOptions.AggressiveOptimization)]
                (v) => new Complex((v.Real + v.Imaginary) / sqrt2, (-v.Real + v.Imaginary) / sqrt2);
            }
            else if (optimize)
            {
                var c = CreateComplexValue(N, n, inverse);
                var realPlusImaginary = c.Real + c.Imaginary;
                var realMinusImaginary = c.Real - c.Imaginary;
                return
                [MethodImpl(MethodImplOptions.AggressiveOptimization)]
                (v) =>
                {
                    var t = (v.Real + v.Imaginary) * c.Real;
                    return new Complex(t - v.Imaginary * realPlusImaginary, t - v.Real * realMinusImaginary);
                };
            }
            else
            {
                var c = CreateComplexValue(N, n, inverse);
                return
                [MethodImpl(MethodImplOptions.AggressiveOptimization)]
                (v) => c * v;
            }
        }

        private static Complex CreateComplexValue(int N, int n, bool inverse)
        {
            System.Diagnostics.Debug.Assert(N > 0);
            System.Diagnostics.Debug.Assert(n >= 0);
            System.Diagnostics.Debug.Assert(n < N);

            if (!inverse && n != 0)
                n = N - n;

            if (n == 0)
            {
                return new Complex(1, 0);
            }
            else if (n * 4L == N)
            {
                return new Complex(0, 1);
            }
            else if (n * 2L == N)
            {
                return new Complex(-1, 0);
            }
            else if (n * 4L == N * 3L)
            {
                return new Complex(0, -1);
            }
            else if (n * 8L < N)
            {
                var (sin, cos) = double.SinCosPi(2d * n / N);
                return new Complex(cos, sin);
            }
            else if (n * 8L < N * 3L)
            {
                var (sin, cos) = double.SinCosPi(2d * n / N - 0.5d);
                return new Complex(-sin, cos);
            }
            else if (n * 8L < N * 5L)
            {
                var (sin, cos) = double.SinCosPi(2d * n / N - 1d);
                return new Complex(-cos, -sin);
            }
            else if (n * 8L < N * 7L)
            {
                var (sin, cos) = double.SinCosPi(2d * n / N - 1.5d);
                return new Complex(sin, -cos);
            }
            else
            {
                var (sin, cos) = double.SinCosPi(2d * n / N - 2d);
                return new Complex(cos, sin);
            }
        }
    }
}
