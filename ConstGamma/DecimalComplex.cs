using System.Globalization;

namespace System.Numerics;

public readonly struct DecimalComplex(decimal real, decimal imaginary) : IEquatable<DecimalComplex>, IFormattable
{
    private readonly decimal real = real;
    private readonly decimal imaginary = imaginary;
    private const decimal LOG_10_INV = 0.43429448190325m;

    public decimal Real => this.real;

    public decimal Imaginary => this.imaginary;

    public decimal Magnitude => Abs(this);

    public decimal Phase => Atan2(this.imaginary, this.real);

    public static DecimalComplex FromPolarCoordinates(decimal magnitude, decimal phase)
        => new(magnitude * Cos(phase), magnitude * Sin(phase));

    public static DecimalComplex Negate(DecimalComplex value)
        => -value;


    public static DecimalComplex Add(DecimalComplex left, DecimalComplex right)
        => left + right;


    public static DecimalComplex Subtract(DecimalComplex left, DecimalComplex right)
        => left - right;


    public static DecimalComplex Multiply(DecimalComplex left, DecimalComplex right)
        => left * right;


    public static DecimalComplex Divide(DecimalComplex dividend, DecimalComplex divisor)
        => dividend / divisor;


    public static DecimalComplex operator -(DecimalComplex value)
        => new(-value.real, -value.imaginary);


    public static DecimalComplex operator +(DecimalComplex left, DecimalComplex right)
        => new(left.real + right.real, left.imaginary + right.imaginary);


    public static DecimalComplex operator -(DecimalComplex left, DecimalComplex right)
        => new(left.real - right.real, left.imaginary - right.imaginary);


    public static DecimalComplex operator *(DecimalComplex left, DecimalComplex right)
    {
        var real = left.real * right.real - left.imaginary * right.imaginary;
        var imaginary = left.imaginary * right.real + left.real * right.imaginary;
        return new(real, imaginary);
    }


    public static DecimalComplex operator /(DecimalComplex left, DecimalComplex right)
    {
        var real = left.real;
        var imaginary = left.imaginary;
        var real2 = right.real;
        var imaginary2 = right.imaginary;
        if (Math.Abs(imaginary2) < Math.Abs(real2))
        {
            var num = imaginary2 / real2;
            return new DecimalComplex((real + imaginary * num) / (real2 + imaginary2 * num), (imaginary - real * num) / (real2 + imaginary2 * num));
        }
        var num2 = real2 / imaginary2;
        return new DecimalComplex((imaginary + real * num2) / (imaginary2 + real2 * num2), (-real + imaginary * num2) / (imaginary2 + real2 * num2));
    }

    public static decimal Abs(DecimalComplex value)
    {
        if (IsInfinity(value.real) || IsInfinity(value.imaginary)) return Decimal.Zero; ;

        var nr = Abs(value.real);
        var ni = Abs(value.imaginary);
        if (nr > ni)
        {
            decimal m = ni / nr;
            return nr * Sqrt(1.0m + m * m);
        }
        if (ni == 0.0m)
        {
            return nr;
        }
        var n = nr / ni;
        return ni * Sqrt(1.0m + n * n);
    }
    public static DecimalComplex Conjugate(DecimalComplex value)
        => new(value.real, -value.imaginary);


    public static DecimalComplex Reciprocal(DecimalComplex value)
        => value.real == 0.0m && value.imaginary == 0.0m ? DecimalComplex.Zero : DecimalComplex.One / value;


    public static bool operator ==(DecimalComplex left, DecimalComplex right)
        => left.real == right.real && left.imaginary == right.imaginary;


    public static bool operator !=(DecimalComplex left, DecimalComplex right)
        => left.real != right.real || left.imaginary != right.imaginary;


    public override bool Equals(object? o)
        => o is DecimalComplex complex && this == complex;


    public bool Equals(DecimalComplex value)
        => this.real.Equals(value.real) && this.imaginary.Equals(value.imaginary);


    public static implicit operator DecimalComplex(short value)
        => new(value, 0.0m);


    public static implicit operator DecimalComplex(int value)
        => new(value, 0.0m);


    public static implicit operator DecimalComplex(long value)
        => new(value, 0.0m);

    [CLSCompliant(false)]

    public static implicit operator DecimalComplex(ushort value)
        => new((decimal)value, 0.0m);

    [CLSCompliant(false)]

    public static implicit operator DecimalComplex(uint value)
        => new(value, 0.0m);

    [CLSCompliant(false)]

    public static implicit operator DecimalComplex(ulong value)
        => new(value, 0.0m);

    [CLSCompliant(false)]

    public static implicit operator DecimalComplex(sbyte value)
        => new(value, 0.0m);


    public static implicit operator DecimalComplex(byte value)
        => new((decimal)value, 0.0m);


    public static implicit operator DecimalComplex(float value)
        => new((decimal)value, 0.0m);


    public static implicit operator DecimalComplex(double value)
        => new((decimal)value, 0.0m);


    public static explicit operator DecimalComplex(BigInteger value)
        => new((decimal)value, 0.0m);


    public static explicit operator DecimalComplex(decimal value)
        => new(value, 0.0m);


    public override string ToString() => string.Format(CultureInfo.CurrentCulture, "({0}, {1})",

            this.real,
            this.imaginary
        );


    public string ToString(string format) => string.Format(CultureInfo.CurrentCulture, "({0}, {1})",

            this.real.ToString(format, CultureInfo.CurrentCulture),
            this.imaginary.ToString(format, CultureInfo.CurrentCulture)
        );


    public string ToString(IFormatProvider provider) => string.Format(provider, "({0}, {1})",

            this.real,
            this.imaginary
        );


    public string ToString(string format, IFormatProvider provider) => string.Format(provider, "({0}, {1})",

            this.real.ToString(format, provider),
            this.imaginary.ToString(format, provider)
        );


    public override int GetHashCode()
    {
        const int Factor = 99999997;
        int n = this.real.GetHashCode() % Factor;
        int hashCode = this.imaginary.GetHashCode();
        return n ^ hashCode;
    }


    public static DecimalComplex Sin(DecimalComplex value)
    {
        var real = value.real;
        var imaginary = value.imaginary;
        return new(Sin(real) * Cosh(imaginary), Cos(real) * Sinh(imaginary));
    }


    public static DecimalComplex Sinh(DecimalComplex value)
    {
        var real = value.real;
        var imaginary = value.imaginary;
        return new(Sinh(real) * Cos(imaginary), Cosh(real) * Sin(imaginary));
    }


    public static DecimalComplex Asin(DecimalComplex value)
        => -ImaginaryOne * Log(ImaginaryOne * value + Sqrt(One - value * value));


    public static DecimalComplex Cos(DecimalComplex value)
    {
        var real = value.real;
        var imaginary = value.imaginary;
        return new (Cos(real) * Cosh(imaginary), -(Sin(real) * Sinh(imaginary)));
    }

    public static DecimalComplex Cosh(DecimalComplex value)
    {
        var real = value.real;
        var imaginary = value.imaginary;
        return new (Cosh(real) * Cos(imaginary), Sinh(real) * Sin(imaginary));
    }


    public static DecimalComplex Acos(DecimalComplex value)
        => -ImaginaryOne * Log(value + ImaginaryOne * Sqrt(One - value * value));


    public static DecimalComplex Tan(DecimalComplex value)
        => Sin(value) / Cos(value);


    public static DecimalComplex Tanh(DecimalComplex value)
        => Sinh(value) / Cosh(value);


    public static DecimalComplex Atan(DecimalComplex value)
    {
        var right = new DecimalComplex(2.0m, 0.0m);
        return ImaginaryOne / right * (Log(One - ImaginaryOne * value) - Log(One + ImaginaryOne * value));
    }

    public static bool IsInfinity(decimal value) 
        => Math.Abs(value) == decimal.MaxValue;

    public static decimal Log(decimal value, decimal baseValue) 
        => Log(value) / Log(baseValue);

    public static DecimalComplex Log(DecimalComplex value) 
        => new (Log(Abs(value)), Atan2(value.imaginary, value.real));

    public static DecimalComplex Log(DecimalComplex value, decimal baseValue) 
        => new (Log(Abs(value)), Atan2(value.imaginary, value.real));


    public static DecimalComplex Log10(DecimalComplex value)
        => Scale(Log(value), LOG_10_INV);

    public static DecimalComplex Exp(DecimalComplex value)
    {
        var num = Exp(value.real);
        var real = num * Cos(value.imaginary);
        var imaginary = num * Sin(value.imaginary);
        return new DecimalComplex(real, imaginary);
    }

    public static DecimalComplex Sqrt(DecimalComplex value) 
        => FromPolarCoordinates(Sqrt(value.Magnitude), value.Phase / 2.0m);


    public static DecimalComplex Pow(DecimalComplex value, DecimalComplex power)
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
        var num2 = Atan2(imaginary, real);
        var num3 = real2 * num2 + imaginary2 * Log(num1);
        var num4 = Pow(num1, real2) * Pow(
            2.718281828459045m, -imaginary2 * num2);
        return new DecimalComplex(num4 * Cos(num3), num4 * Sin(num3));
    }

    public static DecimalComplex Pow(DecimalComplex value, decimal power) 
        => Pow(value, new DecimalComplex(power, 0.0m));

    private static DecimalComplex Scale(DecimalComplex value, decimal factor)
    {
        var real = factor * value.real;
        var imaginary = factor * value.imaginary;
        return new DecimalComplex(real, imaginary);
    }

    public static readonly DecimalComplex Zero = new (0.0m, 0.0m);
    public static readonly DecimalComplex One = new (1.0m, 0.0m);
    public static readonly DecimalComplex ImaginaryOne = new (0.0m, 1.0m);

    public const decimal E = 2.71828182845904523536028747135m;
    public const decimal Pi = 3.14159265358979323846264338327m;
    public const decimal HalfPi = 3.14159265358979323846264338327m / 2.0m;
    public const decimal PI_OVER_4 = Pi / 4.0m;
    public const decimal THREE_PI_OVER_4 = 3.0m * Pi / 4.0m;
    public const decimal EPSILON = 1e-4m; // 精度控制
    public const int MAX_ITERATIONS = 8; // 最大迭代次数

    public static decimal Abs(decimal n) => Math.Abs(n);
    public static decimal Sinh(decimal value) => (Exp(value) - Exp(value)) / 2.0m;
    public static decimal Cosh(decimal value) => (Exp(value) + Exp(value)) / 2.0m;
    public static decimal Exp(decimal value) => Pow(E, value);
    public static decimal Sqrt(decimal value, decimal e = 1e-50m)
    {
        if (value < 0) return decimal.MinValue;
        var x = value;
        var y = (x + value / x) / 2;
        while (Math.Abs(x - y) > e)
        {
            x = y;
            y = (x + value / x) / 2;
        }
        return x;
    }
    public static decimal Sin(decimal value, decimal eps = EPSILON)
    {
        var sum = 0.0m;
        var n = 0L;
        var factorial = 1.0m;
        decimal term;
        while (true)
        {
            term = (n % 2 == 0 ? 1 : -1) * Pow(value, n) / factorial;
            sum += term;
            if (Math.Abs(term) < eps)
                break;
            n += 2;
            factorial *= (n + 1) * (n + 2);
        }
        return sum;
    }

    public static decimal Cos(decimal value, decimal eps = EPSILON) 
        => Sin(value + HalfPi, eps);

    public static decimal Atan2(decimal y, decimal x, int iteration = MAX_ITERATIONS, decimal eps = EPSILON)
    {
        if (x == 0.0m && y == 0.0m) return 0.0m; // 处理(0,0)的情况

        var angle = 0.0m;
        var tx = x;
        var ty = y;
        var invert = false;

        // 处理x为负的情况
        if (x < 0m)
        {
            tx = -x;
            ty = -y;
            invert = true;
        }

        // 处理y为负的情况
        if (y < 0)
        {
            angle = Pi;
            tx = -y;
            ty = x;
        }

        // 处理x为0的情况，即y的正负决定了角度在哪个象限
        if (x == 0) return y > 0 ? PI_OVER_4 : THREE_PI_OVER_4;

        // 泰勒迭代过程
        var factor = PI_OVER_4 / (1 << iteration); // 每次迭代的步长
        for (int i = 0; i < iteration; i++)
        {
            var txTy = tx + ty;
            var tyTx = ty - tx;
            if (Math.Abs(tyTx) > eps) // 检查是否足够小以避免除以零错误
            {
                var ratio = tyTx / txTy; // 新的比率逼近tan(θ)
                angle += factor * Math.Sign(tyTx); // 根据符号调整角度方向
                ty = ratio * tx - ty; // 更新y坐标以逼近正确的角度值
                tx = (ratio * ty) + tx; // 更新x坐标以逼近正确的角度值
            }
            else break; // 如果比率太小，则停止迭代以避免数值问题
        }

        return invert ? Pi - angle : angle; // 根据原始输入调整角度

    }
    public static decimal Pow(decimal x, decimal y, int iterations = MAX_ITERATIONS, bool skip_check = false)
    {
        switch (x)
        {
            case 0.0m when y < 0:
                throw new ArgumentException("Cannot compute zero to a negative power.");
            case 0.0m:
                return 0.0m;
        }
        switch (y)
        {
            case 0.0m:
                return 1.0m;
            case 1.0m:
                return x;
            case -1.0m:
                return 1.0m / x;
            case 0.5m:
                return Sqrt(x);
            case -0.5m:
                return 1.0m / Sqrt(x);
        }
        if (!skip_check)
        {
            switch (y)
            {
                case > 0 and < 1:
                    return 1.0m / Pow(1.0m / x, y, iterations, true);
                case < 0 and > -1:
                    return Pow(1.0m / x, -y, iterations, true);
            }
        }
        {
            var result = x;
            var n = (long)y;
            var fraction = y - n;
            for (int i = 1; i < n; i++) result *= x;
            //TODO
            var approx = result;
            for (int i = 0; i < iterations; i++)
            {
                approx = (approx + x / approx) / 2.0m;
            }
            return approx;
        }
    }
    public static decimal Log(decimal x, int iterations = MAX_ITERATIONS)
    {
        if (x <= 0) throw new ArgumentException("x must be greater than 0");
        if (x == 1) return 0;

        var sum = 0.0m;
        for (int n = 1; n <= iterations; n++)
        {
            var term = Pow(x - 1.0m, n) * ((n + 1) % 2 == 1 ? -1.0m : 1.0m) / n;
            sum += term;
        }
        return 2.0m * sum;
    }
}
