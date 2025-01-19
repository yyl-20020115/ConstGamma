using System.Text;

namespace ConstGamma;

public class Program
{
    public const double Gamma = 0.57721_56649_01532_86060_65120_90082_40243_10421_59335;
    public const double Pi = Math.PI;

    static double CalculateS(int n) 
        => Enumerable.Range(1, n).Aggregate(0.0, (a, b) => a + 1.0 / b);
    static double CalculateT(int n) 
        => (Pi / 8.0 * n) + Gamma;

    static string FormatSEquation(int n)
    {
        var builder = new StringBuilder();
        for(int i = 1; i <= n; i++)
        {
            builder.Append($"1/{i}");
            if (i < n) builder.Append(" + ");
        }
        return builder.ToString();
    }
    static void Main(string[] args)
    {
        const int n = 10;
        for(int i = 1; i < n; i++)
        {
            double s = CalculateS(i);
            double t = CalculateT(i);
            Console.WriteLine($"i={i},\tdiff={s - t},\tt = {Pi}/8 +{Gamma} = {t},\ts = {FormatSEquation(i)}={s}");
        }

    }
}
