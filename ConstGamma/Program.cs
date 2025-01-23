namespace ConstGamma;

public class Program
{
    public const double Gamma = 0.577215664901532860606512090082402431042159335;
    public const double Pi = Math.PI;

    static double CalculateS(int n)
        => Enumerable.Range(1, n).Aggregate(0.0, (a, b) => a + 1.0 / b);
    static double CalculateT(int n)
        => (Pi / 8.0 * n) + Gamma;

    static string FormatEquation(int n)
        => string.Join(" + ", Enumerable.Range(1, n).Select(i => $"1/{i}"));


    static double ZetaOrdinary(long n, long s)
    {
        var sum = 0.0;
        for (long i = 1; i <= n; i++)
        {
            sum += 1.0 / Math.Pow(i, s);
        }
        return sum;
    }
    static double ZetaExponetial(long n, long s)
    {
        var sum = 0.0;
        for (long i = 1; i <= n; i++)
        {
            sum += Math.Pow(Math.E, -s * Math.Log(i));
        }
        return sum;
    }
    static double SumKs(long n, long i)
    {
        var sum = 0.0;
        for (long j = 1; j <= i; j++)
        {
            sum += 1.0 / (n * i + j);
        }
        return sum;
    }
    static double ZetaRecursive(long n, long s)
    {
        var zn = 0.0;
        for (long j = n; j >= 1; j--)
        {
            zn = 1.0 + Math.Pow(Math.E, -s * SumKs(j, n)) * zn;
        }
        return zn;
    }

    static double ZetaRecursiveExtended(long n, long s)
    {
        double zn = 0;
        for (long i = n; i >= 1; i--)
        {
            zn = 1.0 + Math.Pow(Math.E, s * Math.Log(i)) * zn;
        }
        return zn;
    }

    static double ZetaRecursiveLogExtended(long n, long s)
    {
        double zn = n;
        for (long i = n; i >= 1; i--)
        {
            zn = 1.0 + Math.Pow(Math.E, s * Math.Log(i)) * zn;
        }
        return zn;
    }

    static void Main(string[] args)
    {
        ArgumentNullException.ThrowIfNull(args);
        var n = 1000;
        var s = 2;
        var zo = ZetaOrdinary(n, s);
        var ze = ZetaExponetial(n, s);
        var zr = ZetaRecursive(n, s);

        var zre = ZetaRecursiveExtended(n, s);
        var zle = ZetaRecursiveLogExtended(n, s);

        Console.WriteLine($"ZetaOrdinary\t({n}, {s}) = {zo}");
        Console.WriteLine($"ZetaExponetial\t({n}, {s}) = {ze}");
        Console.WriteLine($"ZetaRecursive\t({n}, {s}) = {zr}");

        Console.WriteLine($"ZetaRecursiveExtended\t({n}, {s}) = {zre}");
        Console.WriteLine($"ZetaRecursiveLogExtended\t({n}, {s}) = {zle}");

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
