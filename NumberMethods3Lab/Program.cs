using System;
using System.Text;
using static System.Math;

namespace NumberMethods
{
    public class Program
    {
        const double DefaultEpsilon = 1e-2;
        const double DefaultStep = 1e-8;
        const int DefaultMaxIterations = 100;

        static double CentralDifference(Func<double, double, double> f, double x, double y, double h, string variable) =>
            variable switch
            {
                "x" => (f(x + h, y) - f(x - h, y)) / (2 * h),
                "y" => (f(x, y + h) - f(x, y - h)) / (2 * h),
                _ => throw new ArgumentException("Invalid variable")
            };

        public static (double x, double y) Newton2D(List<Func<double, double, double>> system, double x0, double y0,
            double eps = DefaultEpsilon, int maxIter = DefaultMaxIterations, double h = DefaultStep, bool verbose = false)
        {
            double x = x0, y = y0;

            for (int iter = 0; iter < maxIter; iter++)
            {
                var fValues = system.Select(f => f(x, y)).ToList();
                var jacobian = system
                    .Select(f => (
                        dx: CentralDifference(f, x, y, h, "x"),
                        dy: CentralDifference(f, x, y, h, "y")
                    )).ToList();

                var (f1, f2) = (fValues[0], fValues[1]);
                var (df1x, df1y) = jacobian[0];
                var (df2x, df2y) = jacobian[1];

                double det = df1x * df2y - df1y * df2x;
                if (Abs(det) < 1e-14)
                    throw new InvalidOperationException("Детермінант Якобіана дорівнює нулю");

                double dx = (-f1 * df2y + f2 * df1y) / det;
                double dy = (-df1x * f2 + df2x * f1) / det;

                (x, y) = (x + dx, y + dy);

                if (verbose)
                    Console.WriteLine($"{iter + 1}-ітерація: x = {x}, y = {y}");

                if (Abs(dx) < eps && Abs(dy) < eps)
                    return (x, y);
            }

            throw new Exception("Метод Ньютона не збігся за задану кількість ітерацій.");
        }

        static List<Func<double, double, double>> GetPredefinedEquations()
        {
            return new List<Func<double, double, double>>
            {
                (x, y) => Sin(x) + Sqrt(2 * Pow(y, 3)) - 4,
                (x, y) => Tan(x) - Pow(y, 2) + 4
            };
        }

        public static void Main()
        {
            Console.OutputEncoding = Encoding.UTF8;
            var equations = GetPredefinedEquations();
            var (x, y) = Newton2D(equations, 3.17, 2, verbose: true);

            Console.WriteLine($"\nРозв'язок: x = {x}, y = {y}");
            for (int i = 0; i < equations.Count; i++)
                Console.WriteLine($"Перевірка {i + 1}: f = {equations[i](x, y)}");
        }
    }
}
