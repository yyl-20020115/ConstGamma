using System.Numerics;

namespace ConstGamma;

public class Program
{
    public static readonly Complex[] NontrivalZeros = [
        new (0.5, -14.134725141734693730968711),
        new (0.5, -21.02203964204755334282958),
        new (0.5, -25.01082562396885707143764),
        new (0.5, -28.57182697383300413620860),
        new (0.5, -31.46260704941565719677549),
        new (0.5, -34.13428610227788947322655),
        new (0.5, -36.62904395886886437039029),
        new (0.5, -38.9997672119303999877777),
        new (0.5, -41.2425926398458654312819),
        new (0.5, -43.4167801777518566141529)
        ];


    public const double Gamma = 0.577215664901532860606512090082402431042159335;
    public const double Pi = Math.PI;

    static double CalculateS(int n)
        => Enumerable.Range(1, n).Aggregate(0.0, (a, b) => a + 1.0 / b);
    static double CalculateT(int n)
        => (Pi / 8.0 * n) + Gamma;

    static string FormatEquation(int n)
        => string.Join(" + ", Enumerable.Range(1, n).Select(i => $"1/{i}"));


    /// <summary>
    /// ZetaBasic
    /// Z=sum(1/n^s)
    /// </summary>
    /// <param name="n"></param>
    /// <param name="s"></param>
    /// <returns></returns>
    static double ZetaBasic(long n, double s)
    {
        var sum = 0.0;
        for (long i = 1; i <= n; i++)
        {
            sum += 1.0 / Math.Pow(i, s);
        }
        return sum;
    }
    /// <summary>
    /// ZetaExponential
    /// Z=sum(e^(-s*ln(n)))
    /// verified
    /// </summary>
    /// <param name="n"></param>
    /// <param name="s"></param>
    /// <returns></returns>
    static double ZetaExponential(long n, double s)
    {
        var sum = 0.0;
        for (long i = 1; i <= n; i++)
        {
            sum += Math.Pow(Math.E, -s * Math.Log(i));
        }
        return sum;
    }
    /// <summary>
    /// SumKs
    /// S=sum(1.0/(n*i+1)..1.0/(n*i+i))
    /// verified
    /// </summary>
    /// <param name="n"></param>
    /// <param name="imaginary"></param>
    /// <returns></returns>
    static double SumKs(long n, long imaginary)
    {
        var sum = 0.0;
        for (long j = 1; j <= imaginary; j++)
        {
            sum += 1.0 / (n * imaginary + j);
        }
        return sum;
    }
    /// <summary>
    /// ZetaRecursive
    /// Verified
    /// </summary>
    /// <param name="n"></param>
    /// <param name="s"></param>
    /// <returns></returns>
    static double ZetaRecursive(long n, double s)
    {
        //init z: z_(i+1)= (1+1/n)
        var z = Math.Pow(1.0 + 1.0 / n, s);
        for (long j = n; j >= 1; j--)
        {
            z = 1.0 + Math.Pow(Math.E, -s * SumKs(j, n)) * z;
        }
        return z;
    }

    /// <summary>
    /// ZetaRecursiveExtended
    /// failed
    /// </summary>
    /// <param name="n"></param>
    /// <param name="s"></param>
    /// <returns></returns>
    static double ZetaRecursiveExtended(long n, double s)
    {
        double zn = 1.0;
        for (long i = n; i >= 1; i--)
        {
            zn = 1.0 + Math.Pow((1.0 * i) / (i + 1), s) * zn;
        }
        return zn;
    }

    static double Factorial(long n)
    {
        double f = 1.0;
        for (long i = n; i >= 1; i--)
        {
            f *= i;
        }
        return f;
    }
    static double ZetaRecursiveFactor(long n, double s)
    {
        var sum = 0.0;
        for (long i = 1; i <= n; i++)
        {
            sum += Math.Pow(Factorial(i - 1) / Factorial(i), s);
        }
        return sum;
    }

    static Complex ZetaComplex(double b, double error)
    {
        var sum = Complex.One;
        var bp = b * Math.PI;
        var s = new Complex(0.5, bp);
        for (var i = 1.0; ; i++)
        {
            var v = Complex.Pow(i, -s);
            if (v.Magnitude < error)
            {
                break;
            }
            sum += v;
        }
        return sum;
    }
    static string GetSign(double d) => d >= 0 ? "+" : "";

    class ComplexComparer : IComparer<(Complex, double)>
    {
        public int Compare((Complex, double) x, (Complex, double) y)
        {
            return Math.Abs(x.Item1.Real)>=Math.Abs(y.Item1.Real)?1:0;
        }
    }
    static bool GetZetaResults(double min = -16, long steps = 1600, double error = 5e-3)
    {
        var ms = new List<(Complex c,double b)>();
        for (double b = 0.0; b >= min; b += min / steps)
        {
            var n = ZetaComplex(b, error);
            var z = n / n.Magnitude;

            var m = Math.Abs(n.Real) < Math.Abs(z.Real) ? n : z;
            ms.Add((m,b));
        }

        ms.Sort(new ComplexComparer());
        for (int i = 0; i < ms.Count; i++)
        {
            var m = ms[i].c;
            var b = ms[i].b;
            if (Math.Abs(m.Real) < 1e-2)
            {
                Console.WriteLine($"T:ZetaZero = {GetSign(m.Real)}{m.Real:N8}{GetSign(m.Imaginary)}{m.Imaginary:N8}i, b={GetSign(b)}{b:N8}");
            }
            else
            {
                Console.WriteLine($"F:ZetaZero = {GetSign(m.Real)}{m.Real:N8}{GetSign(m.Imaginary)}{m.Imaginary:N8}i, b={GetSign(b)}{b:N8}");
            }
        }

        return true;
    }
    static void Main(string[] args)
    {
        ArgumentNullException.ThrowIfNull(args);
        for (int i = 0; i < NontrivalZeros.Length; i++)
        {
            Console.WriteLine($"NontrivalZeros[{i}] = {NontrivalZeros[i]}, Imaginary/Pi={NontrivalZeros[i].Imaginary / Math.PI}");
        }

        GetZetaResults();

        var n = 4000;
        var nr = 2200;
        var fr = 171;

        var s = 2;
        var zo = ZetaBasic(n, s);
        var ze = ZetaExponential(n, s);
        var zr = ZetaRecursive(nr, s);

        var zre = ZetaRecursiveExtended(n, s);
        var zrf = ZetaRecursiveFactor(fr, s);

        Console.WriteLine($"ZetaOrdinary    ({n}, {s}) = {zo}");
        Console.WriteLine($"ZetaExponential ({n}, {s}) = {ze}");
        Console.WriteLine($"ZetaRecursive   ({nr}, {s}) = {zr}");
        Console.WriteLine($"ZetaRecursiveExtended   ({n}, {s}) = {zre}");
        Console.WriteLine($"ZetaRecursiveFactor     ({fr}, {s}) = {zrf}");


        if (false)
        {
            const int max = 10;
            for (int i = 1; i < max; i++)
            {
                var _s = CalculateS(i);
                var _t = CalculateT(i);
                Console.WriteLine(
                    $"i={i},\tdiff={_s - _t},\tt = {Pi}/8 +{Gamma} = {_t},\ts = {FormatEquation(i)}={_s}");
            }
        }
    }
}
