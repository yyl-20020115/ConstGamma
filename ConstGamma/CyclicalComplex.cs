using System.Numerics;
using System.Globalization;

namespace ConstGamma;

public readonly struct CyclicalComplex(double real, double imaginary, double semicycle = CyclicalComplex.DefaultSemicycle) : IEquatable<CyclicalComplex>, IFormattable
{
    private readonly double real = real;
    private readonly double imaginary = imaginary;
    private readonly double semicycle = semicycle > 1 ?
        semicycle : throw new ArgumentOutOfRangeException(nameof(semicycle), "Semicycle must be greater than 1.");

    private const double LOG_10_INV = 0.43429448190325;

    public const double DefaultSemicycle = 1.34078079299426e+154;//Math.Sqrt(double.MaxValue)
    public static readonly CyclicalComplex Zero = new(0.0, 0.0, DefaultSemicycle);
    public static readonly CyclicalComplex One = new(1.0, 0.0, DefaultSemicycle);
    public static readonly CyclicalComplex ImaginaryOne = new(0.0, 1.0, DefaultSemicycle);

    public double Real => this.real;
    public double Imaginary => this.imaginary;
    public double Semicycle => this.semicycle;
    public double SemicycleSquared => this.semicycle * this.semicycle;
    public double SemicycleInverse => 1.0 / this.semicycle;
    public double Cycle => this.SemicycleSquared + 1;
    public double Magnitude => Abs(this);
    public double Phase => Math.Atan2(this.imaginary, this.real);
    public double PhaseDegrees => Phase * (180.0 / Math.PI);
    public static CyclicalComplex FromPolarCoordinates(double magnitude, double phase, double semicycle = DefaultSemicycle)
    {
        NormalizePolars(ref magnitude, ref phase, semicycle);
        return new(magnitude * Math.Cos(phase), magnitude * Math.Sin(phase), semicycle);
    }
    public static CyclicalComplex Negate(CyclicalComplex value)
        => -value;

    public static CyclicalComplex Add(CyclicalComplex left, CyclicalComplex right)
        => left + right;

    public static CyclicalComplex Subtract(CyclicalComplex left, CyclicalComplex right)
        => left - right;

    public static CyclicalComplex Multiply(CyclicalComplex left, CyclicalComplex right)
        => left * right;

    public static CyclicalComplex Divide(CyclicalComplex dividend, CyclicalComplex divisor)
        => dividend / divisor;

    public static CyclicalComplex operator -(CyclicalComplex value)
        => new(-value.real, -value.imaginary, value.semicycle);

    public static double NormalizeCycle(ref CyclicalComplex left, ref CyclicalComplex right)
    {
        if (left.semicycle <= 1) throw new ArgumentOutOfRangeException(nameof(left), "Semicycle must be greater than 1.");
        if (right.semicycle <= 1) throw new ArgumentOutOfRangeException(nameof(right), "Semicycle must be greater than 1.");

        var newCycle = left.semicycle;
        // 计算新的周期:以大周期为基准，小周期向大周期靠拢
        if (left.semicycle == right.semicycle)
        {
            newCycle = left.semicycle;
        }
        else if (left.semicycle > right.semicycle)
        {
            // 如果左边的周期大于右边的周期，则需要将右边的周期调整为左边的周期
            right = Scale(right, (newCycle = left.semicycle) / right.semicycle);
        }
        else if (left.semicycle < right.semicycle)
        {
            // 如果右边的周期大于左边的周期，则需要将左边的周期调整为右边的周期
            left = Scale(left, (newCycle = right.semicycle) / left.semicycle);
        }
        return newCycle;
    }

    public static void NormalizeParts(ref double real, ref double imaginary, double semicycle)
    {
        if (semicycle <= 1) throw new ArgumentOutOfRangeException(nameof(semicycle), "Semicycle must be greater than 1.");
        // 如果实部超过周期，则需要将其调整到周期范围内
        if (Math.Abs(real) > semicycle)
        {
            var reminder = Math.IEEERemainder(real, semicycle);
            imaginary += (real - reminder) / semicycle;
            real = reminder;
        }
        var inverseSemicycle = 1.0 / semicycle;
        if (Math.Abs(real) < inverseSemicycle)
        {
            real = 0.0;
        }
    }
    public static void NormalizePolars(ref double magnitude, ref double phase, double semicycle)
    {
        if (semicycle <= 1) throw new ArgumentOutOfRangeException(nameof(semicycle), "Semicycle must be greater than 1.");
        // 如果实部超过周期，则需要将其调整到周期范围内
        if (Math.Abs(magnitude) > semicycle)
        {
            var reminder = Math.IEEERemainder(magnitude, semicycle);
            phase += ((magnitude - reminder) / semicycle) / (2.0 * Math.PI);
            magnitude = reminder;
        }
    }
    /// <summary>
    /// 当两者周期不同的时候，如何相加
    /// 因为Sqrt(-1)=Sqrt(-1*1) 也就是1和无限的几何平方根，
    /// 而无限在这里就是周期。当周期不同的时候，新的周期应当
    /// 为两个周期的几何平方根。在新周期的基础上再重新计算
    /// 两个复数，然后再相加，减法也是一样的。
    /// 周期必须为正数。
    /// </summary>
    /// <param name="left"></param>
    /// <param name="right"></param>
    /// <returns></returns>
    /// 
    public static CyclicalComplex operator +(CyclicalComplex left, CyclicalComplex right)
    {
        var newCycle = NormalizeCycle(ref left, ref right);
        var real = left.real + right.real;
        var imaginary = left.imaginary + right.imaginary;
        NormalizeParts(ref real, ref imaginary, newCycle);
        return new CyclicalComplex(real, imaginary, newCycle);
    }


    public static CyclicalComplex operator -(CyclicalComplex left, CyclicalComplex right)
    {
        var newCycle = NormalizeCycle(ref left, ref right);
        var real = left.real - right.real;
        var imaginary = left.imaginary - right.imaginary;
        NormalizeParts(ref real, ref imaginary, newCycle);
        return new CyclicalComplex(real, imaginary, newCycle);
    }


    public static CyclicalComplex operator *(CyclicalComplex left, CyclicalComplex right)
    {
        var newCycle = NormalizeCycle(ref left, ref right);
        var real = left.real * right.real - left.imaginary * right.imaginary;
        var imaginary = left.imaginary * right.real + left.real * right.imaginary;
        NormalizeParts(ref real, ref imaginary, newCycle);
        return new(real, imaginary, newCycle);
    }


    public static CyclicalComplex operator /(CyclicalComplex left, CyclicalComplex right)
    {
        var newCycle = NormalizeCycle(ref left, ref right);
        var real = left.real;
        var imaginary = left.imaginary;
        var real2 = right.real;
        var imaginary2 = right.imaginary;
        if (Math.Abs(imaginary2) < Math.Abs(real2))
        {
            var n = imaginary2 / real2;
            real = (real + imaginary * n) / (real2 + imaginary2 * n);
            imaginary = (imaginary - real * n) / (real2 + imaginary2 * n);
            NormalizeParts(ref real, ref imaginary, newCycle);
            return new(real, imaginary, newCycle);
        }
        else
        {
            var n = real2 / imaginary2;
            real = (imaginary + real * n) / (imaginary2 + real2 * n);
            imaginary = (-real + imaginary * n) / (imaginary2 + real2 * n);
            NormalizeParts(ref real, ref imaginary, newCycle);
            return new(real, imaginary, newCycle);
        }
    }

    public static double ToIntegral(CyclicalComplex value)
        => value.real + value.imaginary * value.semicycle
        ;

    public static double Abs(CyclicalComplex value)
    {
        if (double.IsInfinity(value.real) || double.IsInfinity(value.imaginary)) return 0.0;

        var nr = Abs(value.real);
        var ni = Abs(value.imaginary);
        if (nr > ni)
        {
            double m = ni / nr;
            return nr * Math.Sqrt(1.0 + m * m);
        }
        if (ni == 0.0)
        {
            return nr;
        }
        var n = nr / ni;
        return ni * Math.Sqrt(1.0 + n * n);
    }
    //共轭
    public static CyclicalComplex Conjugate(CyclicalComplex value)
        => new(value.real, -value.imaginary, value.semicycle);

    //倒数
    public static CyclicalComplex Reciprocal(CyclicalComplex value)
        => value.real == 0.0 && value.imaginary == 0.0 ? Zero : One / value;

    //相等
    public static bool operator ==(CyclicalComplex left, CyclicalComplex right)
        => left.real == right.real && left.imaginary == right.imaginary && left.semicycle == right.semicycle;

    //不等
    public static bool operator !=(CyclicalComplex left, CyclicalComplex right)
        => left.real != right.real || left.imaginary != right.imaginary || left.semicycle != right.semicycle;


    public override bool Equals(object? o)
        => o is CyclicalComplex complex && this == complex
        ;


    public bool Equals(CyclicalComplex value)
        => this.real.Equals(value.real)
        && this.imaginary.Equals(value.imaginary)
        && this.semicycle.Equals(this.semicycle)
        ;


    public static implicit operator CyclicalComplex(short value)
        => new(value, 0.0, DefaultSemicycle);


    public static implicit operator CyclicalComplex(int value)
        => new(value, 0.0, DefaultSemicycle);


    public static implicit operator CyclicalComplex(long value)
        => new(value, 0.0, DefaultSemicycle);

    public static implicit operator CyclicalComplex(ushort value)
        => new(value, 0.0, DefaultSemicycle);

    public static implicit operator CyclicalComplex(uint value)
        => new(value, 0.0, DefaultSemicycle);

    public static implicit operator CyclicalComplex(ulong value)
        => new(value, 0.0, DefaultSemicycle);

    public static implicit operator CyclicalComplex(sbyte value)
        => new(value, 0.0, DefaultSemicycle);


    public static implicit operator CyclicalComplex(byte value)
        => new(value, 0.0, DefaultSemicycle);


    public static implicit operator CyclicalComplex(float value)
        => new((double)value, 0.0, DefaultSemicycle);


    public static implicit operator CyclicalComplex(double value)
        => new((double)value, 0.0, DefaultSemicycle);


    public static explicit operator CyclicalComplex(BigInteger value)
        => new((double)value, 0.0, DefaultSemicycle);


    public override string ToString()
        => $"({this.real}, {this.imaginary}, {this.semicycle})";

    public string ToString(string format)
        => $"({this.real.ToString(format, CultureInfo.CurrentCulture)}, {this.imaginary.ToString(format, CultureInfo.CurrentCulture)}, {this.semicycle.ToString(format, CultureInfo.CurrentCulture)})";


    public string ToString(IFormatProvider provider) => string.Format(provider, "({0}, {1}, {2})",

            this.real,
            this.imaginary,
            this.semicycle
        );

    public string ToString(string format, IFormatProvider provider) => string.Format(provider, "({0}, {1}, {2})",

            this.real.ToString(format, provider),
            this.imaginary.ToString(format, provider),
            this.semicycle.ToString(format, provider)
        );

    public override int GetHashCode()
    {
        const int Factor = 99999997;
        int n = this.real.GetHashCode() % Factor;
        int hashCode = this.imaginary.GetHashCode();
        return n ^ hashCode;
    }


    public static CyclicalComplex Sin(CyclicalComplex value)
    {
        var real = value.real;
        var imaginary = value.imaginary;
        real = Math.Sin(real) * Math.Cosh(imaginary);
        imaginary = Math.Cos(real) * Math.Sinh(imaginary);
        NormalizeParts(ref real, ref imaginary, value.semicycle);
        return new(real, imaginary, value.semicycle);
    }

    public static CyclicalComplex Asin(CyclicalComplex value)
    {
        var io = new CyclicalComplex(0.0, 1.0, value.semicycle);
        var ro = new CyclicalComplex(1.0, 0.0, value.semicycle);
        return -io * Log(io * value + Sqrt(ro - value * value));
    }

    public static CyclicalComplex Sinh(CyclicalComplex value)
    {
        var real = value.real;
        var imaginary = value.imaginary;
        real = Math.Sinh(real) * Math.Cos(imaginary);
        imaginary = Math.Cosh(real) * Math.Sin(imaginary);
        NormalizeParts(ref real, ref imaginary, value.semicycle);
        return new(real, imaginary, value.semicycle);
    }

    public static CyclicalComplex Cos(CyclicalComplex value)
    {
        var real = value.real;
        var imaginary = value.imaginary;
        real = Math.Cos(real) * Math.Cosh(imaginary);
        imaginary = -(Math.Sin(real) * Math.Sinh(imaginary));
        NormalizeParts(ref real, ref imaginary, value.semicycle);
        return new(real, imaginary, value.semicycle);
    }

    public static CyclicalComplex Cosh(CyclicalComplex value)
    {
        var real = value.real;
        var imaginary = value.imaginary;
        real = Math.Cosh(real) * Math.Cos(imaginary);
        imaginary = Math.Sinh(real) * Math.Sin(imaginary);
        NormalizeParts(ref real, ref imaginary, value.semicycle);
        return new(real, imaginary, value.semicycle);
    }

    public static CyclicalComplex Acos(CyclicalComplex value)
    {
        var io = new CyclicalComplex(0.0, 1.0, value.semicycle);
        var ro = new CyclicalComplex(1.0, 0.0, value.semicycle);
        return -io * Log(value + io * Sqrt(ro - value * value));
    }

    public static CyclicalComplex Tan(CyclicalComplex value)
        => Sin(value) / Cos(value);

    public static CyclicalComplex Tanh(CyclicalComplex value)
        => Sinh(value) / Cosh(value);

    public static CyclicalComplex Atan(CyclicalComplex value)
    {
        var right = new CyclicalComplex(2.0, 0.0, value.semicycle);
        var io = new CyclicalComplex(0.0, 1.0, value.semicycle);
        var ro = new CyclicalComplex(1.0, 0.0, value.semicycle);
        return io / right * (Log(ro - io * value) - Log(ro + io * value));
    }

    public static CyclicalComplex Log(CyclicalComplex value)
    {
        var real = Math.Log(Abs(value));
        var imaginary = Math.Atan2(value.imaginary, value.real);
        NormalizeParts(ref real, ref imaginary, value.semicycle);
        return new(real, imaginary, value.semicycle);
    }

    public static CyclicalComplex Log(CyclicalComplex value, double baseValue)
        => Log(value) / Log(baseValue)
        ;

    public static CyclicalComplex Log10(CyclicalComplex value)
        => Scale(Log(value), LOG_10_INV)
        ;

    public static CyclicalComplex Exp(CyclicalComplex value)
    {
        var num = Math.Exp(value.real);
        var real = num * Math.Cos(value.imaginary);
        var imaginary = num * Math.Sin(value.imaginary);
        NormalizeParts(ref real, ref imaginary, value.semicycle);
        return new(real, imaginary, value.semicycle);
    }

    public static CyclicalComplex Sqrt(CyclicalComplex value)
        => FromPolarCoordinates(Math.Sqrt(value.Magnitude), value.Phase / 2.0, value.semicycle)
        ;

    public static CyclicalComplex Pow(CyclicalComplex value, CyclicalComplex power)
    {
        if (power == Zero)
        {
            return One;
        }
        else if (value == Zero)
        {
            return Zero;
        }
        var real = value.real;
        var imaginary = value.imaginary;
        var real2 = power.real;
        var imaginary2 = power.imaginary;
        var num1 = Abs(value);
        var num2 = Math.Atan2(imaginary, real);
        var num3 = real2 * num2 + imaginary2 * Math.Log(num1);
        var num4 = Math.Pow(num1, real2) * Math.Pow(Math.E, -imaginary2 * num2);
        real = num4 * Math.Cos(num3);
        imaginary = num4 * Math.Sin(num3);
        NormalizeParts(ref real, ref imaginary, value.semicycle);
        return new(real, imaginary, value.semicycle);
    }

    public static CyclicalComplex Pow(CyclicalComplex value, double power)
        => Pow(value, new CyclicalComplex(power, 0.0, value.semicycle))
        ;

    private static CyclicalComplex Scale(CyclicalComplex value, double factor)
    {
        var real = factor * value.real;
        var imaginary = factor * value.imaginary;
        NormalizeParts(ref real, ref imaginary, value.semicycle);
        return new(real, imaginary, value.semicycle);
    }
}
