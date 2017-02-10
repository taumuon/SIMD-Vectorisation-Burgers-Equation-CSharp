using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Burgers
{
    class Program
    {
        static Vector<double> dummy;
        static Program() { dummy = Vector<double>.One; }

        private static double BurgersAnalytical(double t, double x, double nu)
        {
            return -2 * nu * (-(-8 * t + 2 * x) * Math.Exp(-Math.Pow((-4 * t + x), 2) / (4 * nu * (t + 1))) / (4 * nu * (t + 1)) - (-8 * t + 2 * x - 12.5663706143592) * Math.Exp(-Math.Pow(-4 * t + x - 6.28318530717959, 2) / (4 * nu * (t + 1))) / (4 * nu * (t + 1))) / (Math.Exp(-Math.Pow(-4 * t + x - 6.28318530717959, 2) / (4 * nu * (t + 1))) + Math.Exp(-Math.Pow(-4 * t + x, 2) / (4 * nu * (t + 1)))) + 4;
        }

        private const string path = @"C:\FILES\Learning\mooc\numerical methods\BurgersCode\output";

        private static double[] linspace(double first, double last, int num)
        {
            var step = (last - first) / (double)num;
            return Enumerable.Range(0, num).Select(v => (v * step) + first).ToArray();
        }

        private static double[] GetAnalytical(double[] x, double t, double nu)
        {
            double[] u = new double[x.Length];

            for (int i = 0; i < x.Length; ++i)
            {
                u[i] = BurgersAnalytical(t, x[i], nu);
            }

            return u;
        }

        private static double[] GetCalculated(int nt, int nx, double dx, double dt, double nu, double[] u)
        {
            for (int tStep = 0; tStep < nt; tStep++)
            {
                double[] un = new double[nx];
                Array.Copy(u, un, u.Length);

                for (int i = 1; i < nx - 1; i++)
                {
                    u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) + Math.Pow(nu * dt / dx, 2.0) *
                            (un[i + 1] - 2 * un[i] + un[i - 1]);
                }

                u[0] = un[0] - un[0] * dt / dx * (un[0] - un[nx - 1]) + Math.Pow(nu * dt / dx, 2.0) *
                            (un[1] - 2 * un[0] + un[nx - 1]);

                u[nx - 1] = un[nx - 1] - un[nx - 1] * dt / dx * (un[nx - 1] - un[nx - 2]) + Math.Pow(nu * dt / dx, 2.0) *
                            (un[0] - 2 * un[nx - 1] + un[nx - 2]);
            }

            return u;
        }

        // Reduce new array allocation and copying, ping-pong between them
        private static double[] GetCalculated2(int nt, int nx, double dx, double dt, double nu, double[] initial)
        {
            double[] u = new double[nx];
            double[] un = new double[nx];
            Array.Copy(initial, un, un.Length);

            for (int tStep = 0; tStep < nt; tStep++)
            {
                for (int i = 1; i < nx - 1; i++)
                {
                    u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) + Math.Pow(nu * dt / dx, 2.0) *
                            (un[i + 1] - 2 * un[i] + un[i - 1]);
                }

                u[0] = un[0] - un[0] * dt / dx * (un[0] - un[nx - 1]) + Math.Pow(nu * dt / dx, 2.0) *
                            (un[1] - 2 * un[0] + un[nx - 1]);

                u[nx - 1] = un[nx - 1] - un[nx - 1] * dt / dx * (un[nx - 1] - un[nx - 2]) + Math.Pow(nu * dt / dx, 2.0) *
                            (un[0] - 2 * un[nx - 1] + un[nx - 2]);

                double[] swap = u;
                u = un;
                un = swap;
            }

            return un;
        }

        // Pull calculation of (nu * dt / dx)^2 out into a variable
        private static double[] GetCalculated3(int nt, int nx, double dx, double dt, double nu, double[] initial)
        {
            double[] u = new double[nx];
            double[] un = new double[nx];
            Array.Copy(initial, un, un.Length);

            double factor = Math.Pow(nu * dt / dx, 2.0);

            for (int tStep = 0; tStep < nt; tStep++)
            {
                for (int i = 1; i < nx - 1; i++)
                {
                    u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) + factor *
                            (un[i + 1] - 2 * un[i] + un[i - 1]);
                }

                u[0] = un[0] - un[0] * dt / dx * (un[0] - un[nx - 1]) + factor *
                            (un[1] - 2 * un[0] + un[nx - 1]);

                u[nx - 1] = un[nx - 1] - un[nx - 1] * dt / dx * (un[nx - 1] - un[nx - 2]) + factor *
                            (un[0] - 2 * un[nx - 1] + un[nx - 2]);

                double[] swap = u;
                u = un;
                un = swap;
            }

            return un;
        }

        private static double[] GetCalculated4(int nt, int nx, double dx, double dt, double nu, double[] initial)
        {
            var nx2 = nx + (Vector<double>.Length - (nx % Vector<double>.Length));

            double[] u = new double[nx2];
            double[] un = new double[nx2];
            Array.Copy(initial, un, initial.Length);

            double factor = Math.Pow(nu * dt / dx, 2.0);

            for (int tStep = 0; tStep < nt; tStep++)
            {
                for (int i = 1; i < nx2 - 1; i += Vector<double>.Length)
                {
                    var vectorIn0 = new Vector<double>(un, i);
                    var vectorInPrev = new Vector<double>(un, i - 1);
                    var vectorInNext = new Vector<double>(un, i + 1);

                    var vectorOut = vectorIn0 - vectorIn0 * (dt / dx) * (vectorIn0 - vectorInPrev) + factor *
                        (vectorInNext - 2.0 * vectorIn0 + vectorInPrev);

                    vectorOut.CopyTo(u, i);
                }

                //for (int i = 1; i < nx - 1; i++)
                //{
                //    u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) + factor *
                //            (un[i + 1] - 2 * un[i] + un[i - 1]);
                //}

                u[0] = un[0] - un[0] * dt / dx * (un[0] - un[nx - 1]) + factor *
                            (un[1] - 2 * un[0] + un[nx - 1]);

                u[nx - 1] = un[nx - 1] - un[nx - 1] * dt / dx * (un[nx - 1] - un[nx - 2]) + factor *
                            (un[0] - 2 * un[nx - 1] + un[nx - 2]);

                double[] swap = u;
                u = un;
                un = swap;
            }

            return un;
        }

        static void Main(string[] args)
        {
            if (!VectorMath.IsHardwareAccelerated)
            {
                Console.WriteLine("Not hardware accelerated!");
            }
            else
            {
                Console.WriteLine("Vector<double>.Length: " + Vector<double>.Length);
            }

            int nx = 10001;
            int nt = 1000;
            double dx = 2.0 * Math.PI / (nx - 1.0);
            double nu = 0.07;
            double dt = dx * nu;

            double[] x = linspace(0.0, 2.0 * Math.PI, nx);

            double[] initial = GetAnalytical(x, 0.0, nu);

            // double[] justToWarmupJit = GetCalculated(1, nx, dx, dt, nu, initial);
            double[] justToWarmupJit = GetCalculated4(1, nx, dx, dt, nu, initial);

            var stopwatch = new System.Diagnostics.Stopwatch();
            stopwatch.Start();
            // double[] calculated = GetCalculated(nt, nx, dx, dt, nu, initial);
            double[] calculated = GetCalculated4(nt, nx, dx, dt, nu, initial);
            stopwatch.Stop();

            double[] analytical = GetAnalytical(x, nt * dt, nu);

            var fileName = System.IO.Path.Combine(path, "out.dat");
            using (var writer = new StreamWriter(fileName))
            {
                writer.WriteLine(string.Format("#x\tanalytical\tcalculated"));
                for (int i = 0; i < nx; ++i)
                {
                    writer.WriteLine(string.Format("{0}\t{1}\t{2}", x[i], analytical[i], calculated[i]));
                }
            }

            Console.WriteLine("Elapsed: " + stopwatch.ElapsedMilliseconds);

        }
    }
}
