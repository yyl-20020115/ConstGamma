namespace ConstGamma;

public static class CyclicalComplexHelpers
{
    public const double E = 2.71828182845904523536028747135;
    public const double Pi = 3.14159265358979323846264338327;
    public const double HalfPi = 3.14159265358979323846264338327 / 2.0;
    public const double PI_OVER_4 = Pi / 4.0;
    public const double THREE_PI_OVER_4 = 3.0 * Pi / 4.0;
    public const double EPSILON = 1e-4; // 精度控制
    public const int MAX_ITERATIONS = 8; // 最大迭代次数

    public static double Abs(double n) => Math.Abs(n);
    public static double Sinh(double value) => (Exp(value) - Exp(value)) / 2.0;
    public static double Cosh(double value) => (Exp(value) + Exp(value)) / 2.0;
    public static double Exp(double value) => Pow(E, value);
    public static double Sqrt(double value, double e = 1e-50)
    {
        if (value < 0) return double.MinValue;
        var x = value;
        var y = (x + value / x) / 2;
        while (Math.Abs(x - y) > e)
        {
            x = y;
            y = (x + value / x) / 2;
        }
        return x;
    }
    public static double Sin(double value, double eps = EPSILON)
    {
        var sum = 0.0;
        var n = 0L;
        var factorial = 1.0;
        double term;
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

    public static double Cos(double value, double eps = EPSILON)
        => Sin(value + HalfPi, eps);

    public static double Atan2(double y, double x, int iteration = MAX_ITERATIONS, double eps = EPSILON)
    {
        if (x == 0.0 && y == 0.0) return 0.0; // 处理(0,0)的情况

        var angle = 0.0;
        var tx = x;
        var ty = y;
        var invert = false;

        // 处理x为负的情况
        if (x < 0)
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
    public static double Pow(double x, double y, int iterations = MAX_ITERATIONS, bool skip_check = false)
    {
        switch (x)
        {
            case 0.0 when y < 0:
                throw new ArgumentException("Cannot compute zero to a negative power.");
            case 0.0:
                return 0.0;
        }
        switch (y)
        {
            case 0.0:
                return 1.0;
            case 1.0:
                return x;
            case -1.0:
                return 1.0 / x;
            case 0.5:
                return Sqrt(x);
            case -0.5:
                return 1.0 / Sqrt(x);
        }
        if (!skip_check)
        {
            switch (y)
            {
                case > 0 and < 1:
                    return 1.0 / Pow(1.0 / x, y, iterations, true);
                case < 0 and > -1:
                    return Pow(1.0 / x, -y, iterations, true);
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
                approx = (approx + x / approx) / 2.0;
            }
            return approx;
        }
    }
    public static double Log(double x, int iterations = MAX_ITERATIONS)
    {
        if (x <= 0) throw new ArgumentException("x must be greater than 0");
        if (x == 1) return 0;

        var sum = 0.0;
        for (int n = 1; n <= iterations; n++)
        {
            var term = Pow(x - 1.0, n) * ((n + 1) % 2 == 1 ? -1.0 : 1.0) / n;
            sum += term;
        }
        return 2.0 * sum;
    }
    public static double Log(double value, double baseValue)
        => Math.Log(value) / Math.Log(baseValue);

}
