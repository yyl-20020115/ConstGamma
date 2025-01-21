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

    static void Main(string[] args)
    {
        ArgumentNullException.ThrowIfNull(args);

        const int n = 10;
        for(int i = 1; i < n; i++)
        {
            var s = CalculateS(i);
            var t = CalculateT(i);
            Console.WriteLine(
                $"i={i},\tdiff={s - t},\tt = {Pi}/8 +{Gamma} = {t},\ts = {FormatEquation(i)}={s}");
        }
    }
}
