//#define COMPUTE_FFT

using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.ComponentModel.Design;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;

#if COMPUTE_FFT
using FFT;
#endif

namespace Experiment.CUI
{
    internal sealed class Program
    {
        private class FftTreeElement
        {
            public FftTreeElement(uint length, uint index, string expression)
            {
                Length = length;
                Index = index;
                Expression = expression;
            }

            public uint Length{get;}
            public uint Index { get; }
            public string Expression { get; }
        }

        private static void Main()
        {
#if COMPUTE_FFT
            ComputeFft();
#else
            var m = new SquareMatrix<int>(32);
            for (var row = 0; row < m.Size; ++row)
            {
                for (var column = 0; column < m.Size; ++column)
                {
                    var v = row * column % 32;
                    m[row, column] = v;
                }
            }

            for (var column = 0; column < m.Size; ++column)
            {
                var column2 = (int)ReverseBits((uint)column, m.Size);
                if (column > column2)
                    m.SwapColumn(column, column2);
            }

            for (var row = 0; row < m.Size; ++row)
            {
                Console.WriteLine($"D({row}),,{
                    string.Join(",", 
                        Enumerable.Range(0, m.Size)
                        .Select(column => 
                            m[row, column] == 0
                            ? "1"
                            : m[row, column] == 16
                            ? "-1"
                            : m[row, column] > 16
                            ? $"-w^{m[row, column] - 16}"
                            : $"w^{m[row, column]}"
                            ))},,S'({row})");
            }

#endif

            Console.Beep();
            _ = Console.ReadLine();
        }

#if COMPUTE_FFT
        private static void ComputeFft()
        {
            foreach (var count in Enumerable.Range(1, 4).Append(8).Append(16).Where(count => count >= 4 && int.IsPow2(count)))
            {
                var source = Enumerable.Range(0, count).Select(n => new Complex(n, 0)).ToArray();
                var destination1 = new Complex[count];
                var destination2 = new Complex[count];
                var destination3 = new Complex[count];
                var destination4 = new Complex[count];
                var destination5 = new Complex[count];
                new BasicDiscreteFourierTransform().Transform(source, destination1);
                new RecursiveCooleyTukeyFastFourierTransform( DiscreteFourierTransform.Options.OptimizeTwiddleFactorMultiplication).Transform(source, destination2);
                new Base2FastFourierTransform(DiscreteFourierTransform.Options.OptimizeTwiddleFactorMultiplication).Transform(source, destination3);
                new SplitRadixFastFourierTransform(DiscreteFourierTransform.Options.OptimizeTwiddleFactorMultiplication).Transform(source, destination4);
                new Base4FastFourierTransform(DiscreteFourierTransform.Options.OptimizeTwiddleFactorMultiplication).Transform(source, destination5);
                Console.WriteLine($"{count}: RecursiveCooleyTukeyFastFourierTransform: {EvaluateError(destination1, destination2)}");
                Console.WriteLine($"{count}: Base2FastFourierTransform: {EvaluateError(destination1, destination3)}");
                Console.WriteLine($"{count}: SplitRadixFastFourierTransform: {EvaluateError(destination1, destination4)}");
                Console.WriteLine($"{count}: Base4FastFourierTransform: {EvaluateError(destination1, destination5)}");
            }

            static double EvaluateError(Complex[] values1, Complex[] values2)
            {
                var error =
                    double.Sqrt(
                        values1
                        .Zip(values2, (v1, v2) =>
                        {
                            var diff = v1 - v2;
                            return diff.Real * diff.Real + diff.Imaginary * diff.Imaginary;
                        })
                        .Sum(n => n) / values1.Length);
                return double.Abs(error) < 1e-14 ? 0d : error;
            }
        }
#endif

        private static void BuildFftTree(IDictionary<(uint lebgth, uint index), FftTreeElement> tree, uint length, uint index)
        {
#if false
https://www.HostMath.com/Show.aspx?Code=D(N%2C%20n_0%2C%20k)%0A%3DD(%5Cfrac%7BN%7D%7B2%7D%2C%20n_0%2C%20k)%0A%2B(%5Comega_N%5Ek%20D(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%2B%5Cfrac%7BN%7D%7B2%7D%2C%20k)%0A%2B%5Comega_N%5E%7B3k%7DD(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%2B%5Cfrac%7B3%7D%7B4%7DN%2C%20k))%5C%5C%0A%0AD(N%2C%20n_0%2C%20k%2B%5Cfrac%7BN%7D%7B2%7D)%0A%3DD(%5Cfrac%7BN%7D%7B2%7D%2C%20n_0%2C%20k)%0A-(%5Comega_N%5Ek%20D(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%20%2B%20%5Cfrac%7BN%7D%7B2%7D%2C%20k)%0A%2B%5Comega_N%5E%7B3k%7D%20D(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%20%2B%20%5Cfrac%7B3%7D%7B4%7DN%2C%20k))%5C%5C%0A%0AD(N%2C%20n_0%2C%20k%2B%5Cfrac%7BN%7D%7B4%7D)%0A%3DD(%5Cfrac%7BN%7D%7B2%7D%2C%20n_0%2C%20k%2B%5Cfrac%7BN%7D%7B4%7D)%0A-i(%5Comega_N%5Ek%20D(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%2B%5Cfrac%7BN%7D%7B2%7D%2C%20k)%0A-%5Comega_N%5E%7B3k%7DD(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%2B%5Cfrac%7B3%7D%7B4%7DN%2C%20k))%0A%5C%5C%0A%0AD(N%2C%20n_0%2C%20k%2B%5Cfrac%7B3%7D%7B4%7DN)%0A%3DD(%5Cfrac%7BN%7D%7B2%7D%2C%20n_0%2C%20k%2B%5Cfrac%7BN%7D%7B4%7D)%0A%2Bi(%5Comega_N%5Ek%20D(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%2B%5Cfrac%7BN%7D%7B2%7D%2C%20k)%0A-%5Comega_N%5E%7B3k%7DD(%5Cfrac%7BN%7D%7B4%7D%2C%20n_0%2B%5Cfrac%7B3%7D%7B4%7DN%2C%20k))%0A%5C%5C%0A
#endif
            if (!uint.IsPow2(length))
                throw new Exception();
            var exp = uint.TrailingZeroCount(length);
            if (exp == 0)
            {
                AddToTree(tree, new FftTreeElement(length, index, $"D{length}({index}) = S'({index})"));
            }
            else if (uint.IsEvenInteger(exp))
            {
                // length が 4 のべき乗である場合
                if (index < length / 4)
                {
                }
                else if (index < length / 2)
                {
                }
                else if (index < length / 4 * 3)
                {
                }
                else
                {
                }
            }
            else if (exp > 0)
            {
                // length が 2 のべき乗である場合
                if (index < length / 2)
                {
                    AddToTree(tree, new FftTreeElement(length, index, $"D{length}({index}) = D{length / 2}({index}){GetOmega(length, index)}D{length / 2}({index + length / 2})"));
                    BuildFftTree(tree, length / 2, index);
                    BuildFftTree(tree, length / 2, index + length / 2);
                }
                else
                {
                    AddToTree(tree, new FftTreeElement(length, index, $"D{length}({index}) = D{length / 2}({index - length / 2}){GetOmega(length, index)}D{length / 2}({index})"));
                    BuildFftTree(tree, length / 2, index - length / 2);
                    BuildFftTree(tree, length / 2, index);
                }
            }
            else
            {
                throw new Exception();
            }

            static void AddToTree(IDictionary<(uint lebgth, uint index), FftTreeElement> tree, FftTreeElement treeElement)
            {
                ???
            }

            static string GetOmega(uint length, uint index)
            {
                var numerator = index % length;
                var denominator = length;
                var gcd = BigInteger.GreatestCommonDivisor(numerator, denominator);
                numerator /= checked((uint)gcd);
                denominator /= checked((uint)gcd);
                if (index == 0)
                    return" + ";
                if (denominator == 2)
                    return " - ";
                if (denominator == 4)
                {
                    if (numerator == 1)
                        return " - i * ";
                    else if (numerator == 3)
                        return " + i * ";
                    else
                        throw new Exception();
                }

                if (numerator * 2 < denominator)
                    return $" + w^({numerator}/{denominator}) * ";
                else
                    return $" - w^({numerator - denominator - 2}/{denominator}) * ";
            }
        }
    }
}
