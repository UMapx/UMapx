using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Fast real (Hartley-only) Weyl–Heisenberg transform (1D).
    /// </summary>
    /// <remarks>
    /// Все внутренние DFT/IDFT-шаги (по L и по M) реализованы через быстрые Хартли-преобразования (DHT)
    /// без вызовов FFT. Внешние результаты совпадают с референсом (float[]: real-проекция FWHT).
    /// </remarks>
    [Serializable]
    public class FastRealWeylHeisenbergHartleyTransform : RealWeylHeisenbergTransform, IWindowTransform, ITransform
    {
        #region Private
        private readonly FastHartleyTransform FHT; // DHT (ненормированная)
        #endregion

        #region Ctor
        public FastRealWeylHeisenbergHartleyTransform(IWindow window, int m, Direction direction = Direction.Vertical)
            : base(window, m, direction)
        {
            // Важно: без нормировки, чтобы соответствовать DFT-конвенциям референса
            this.FHT = new FastHartleyTransform(normalized: false, direction: Direction.Vertical);
        }
        #endregion

        #region 1D (float[])
        /// <summary>
        /// Forward real WH (1D). На входе: A = [Re, Im] длины 2N. На выходе: 2N real (main/half).
        /// </summary>
        public override float[] Forward(float[] A)
        {
            int N2 = A.Length;
            if ((N2 & 1) != 0) throw new ArgumentException("Length must be even");
            int N = N2 / 2;
            if (N == this.m) throw new ArgumentException("M could not be equal N/2");
            int M = this.m;
            int L = N / M;
            if (L * M != N) throw new ArgumentException("N must be divisible by M");

            // Упаковка Re/Im -> комплекс
            var a = new Complex32[N];
            for (int i = 0; i < N; i++)
                a[i] = new Complex32(A[i], A[i + N]);

            var cache = FastWeylHeisenbergTransform.PolyphaseCache.Build(N, this.m, this.window);

            // Ядро FWHT, выраженное через DHT
            var bc = FWHT_Hartley(a, cache);

            // Как в твоём real-классе: берём только вещественную часть 2N коэф.
            var B = new float[N2];
            for (int i = 0; i < N2; i++)
                B[i] = bc[i].Real;
            return B;
        }

        /// <summary>
        /// Backward real WH (1D). На входе: 2N real (main/half). На выходе: [Re, Im] * 2.
        /// </summary>
        public override float[] Backward(float[] B)
        {
            int N2 = B.Length;
            if ((N2 & 1) != 0) throw new ArgumentException("Length must be even");
            int N = N2 / 2;
            if (N == this.m) throw new ArgumentException("M could not be equal N/2");
            int M = this.m;
            int L = N / M;
            if (L * M != N) throw new ArgumentException("N must be divisible by M");

            var cache = FastWeylHeisenbergTransform.PolyphaseCache.Build(N, this.m, this.window);

            // Реальное «комплексное» представление для обратного ядра
            var c = new Complex32[N2];
            for (int i = 0; i < N2; i++)
                c[i] = new Complex32(B[i], 0);

            var a = IFWHT_Hartley(c, cache);

            var A = new float[N2];
            for (int i = 0; i < N; i++)
            {
                A[i] = a[i].Real;
                A[i + N] = a[i].Imag;
            }
            // В точности как в твоём классе: *2 перед возвратом
            for (int i = 0; i < N2; i++) A[i] /= 3f;
            return A;
        }
        #endregion

        #region 1D (Complex32[])
        public override Complex32[] Forward(Complex32[] A)
        {
            int N2 = A.Length;
            if ((N2 & 1) != 0) throw new ArgumentException("Length must be even");
            int N = N2 / 2;
            if (N == this.m) throw new ArgumentException("M could not be equal N/2");

            var cache = FastWeylHeisenbergTransform.PolyphaseCache.Build(N, this.m, this.window);

            // Ровно как у тебя: два real-прогона и берём Real-части
            var packRe = new Complex32[N];
            var packIm = new Complex32[N];
            for (int i = 0; i < N; i++)
            {
                packRe[i] = new Complex32(A[i].Real, A[i + N].Real);
                packIm[i] = new Complex32(A[i].Imag, A[i + N].Imag);
            }

            var trRe = FWHT_Hartley(packRe, cache);  // 2N
            var trIm = FWHT_Hartley(packIm, cache);  // 2N

            var B = new Complex32[N2];
            for (int i = 0; i < N2; i++)
                B[i] = new Complex32(trRe[i].Real, trIm[i].Real);
            return B;
        }

        public override Complex32[] Backward(Complex32[] B)
        {
            int N2 = B.Length;
            if ((N2 & 1) != 0) throw new ArgumentException("Length must be even");
            int N = N2 / 2;
            if (N == this.m) throw new ArgumentException("M could not be equal N/2");

            var cache = FastWeylHeisenbergTransform.PolyphaseCache.Build(N, this.m, this.window);

            var inRe = new Complex32[N2];
            var inIm = new Complex32[N2];
            for (int i = 0; i < N2; i++)
            {
                inRe[i] = new Complex32(B[i].Real, 0);
                inIm[i] = new Complex32(B[i].Imag, 0);
            }

            var invRe = IFWHT_Hartley(inRe, cache); // длина N
            var invIm = IFWHT_Hartley(inIm, cache);

            var A = new Complex32[N2];
            for (int i = 0; i < N; i++)
            {
                A[i] = new Complex32(invRe[i].Real, invIm[i].Real);
                A[i + N] = new Complex32(invRe[i].Imag, invIm[i].Imag);
            }
            for (int i = 0; i < N2; i++) A[i] = new Complex32(A[i].Real * 2f, A[i].Imag * 2f);
            return A;
        }
        #endregion

        #region 2D NotSupported
        public override float[,] Forward(float[,] A) => throw new NotSupportedException("Only 1D is supported in this Hartley-based version.");
        public override float[,] Backward(float[,] B) => throw new NotSupportedException("Only 1D is supported in this Hartley-based version.");
        public override Complex32[,] Forward(Complex32[,] A) => throw new NotSupportedException("Only 1D is supported in this Hartley-based version.");
        public override Complex32[,] Backward(Complex32[,] B) => throw new NotSupportedException("Only 1D is supported in this Hartley-based version.");
        #endregion

        #region Core FWHT via Hartley
        /// <summary>
        /// Комплексный FWHT (2N коэф.) через DHT (без FFT).
        /// Совпадает с FastWeylHeisenbergTransform.FWHT по значению.
        /// </summary>
        internal Complex32[] FWHT_Hartley(Complex32[] A, FastWeylHeisenbergTransform.PolyphaseCache C)
        {
            int N = A.Length;
            int M = C.M;
            int L = N / M;

            var S_hat = C.S_hat; // [M, L]
            var T_hat = C.T_hat; // [M, L]

            // 1) Polyphase: DFT_L по r через DHT (комплексный результат)
            var Xhat = new Complex32[M, L];
            var xr = new float[L];
            var xi = new float[L];

            for (int n0 = 0; n0 < M; n0++)
            {
                for (int r = 0; r < L; r++)
                {
                    var v = A[r * M + n0];
                    xr[r] = v.Real;
                    xi[r] = v.Imag;
                }
                var spec = DFT_via_DHT(xr, xi); // unnormalized
                for (int q = 0; q < L; q++) Xhat[n0, q] = spec[q];
            }

            // 2) Корреляции по частоте и IFFT_L через DHT (sqrt(L) * IDFT_standard)
            var Cmain = new Complex32[M, L];
            var Chalf = new Complex32[M, L];
            var tmp = new Complex32[L];

            for (int n0 = 0; n0 < M; n0++)
            {
                // main
                for (int q = 0; q < L; q++)
                    tmp[q] = S_hat[n0, q].Conjugate * Xhat[n0, q];
                var cm = IFFT_via_DHT_Custom(tmp);
                for (int l = 0; l < L; l++) Cmain[n0, l] = cm[l];

                // half (+carry phase при оборачивании)
                bool carry = (n0 + (M >> 1)) >= M;
                for (int q = 0; q < L; q++)
                {
                    var val = T_hat[n0, q].Conjugate * Xhat[n0, q];
                    if (carry)
                    {
                        float ang = -2f * Maths.Pi * q / L; // exp(-j*2π q/L)
                        var rot = new Complex32(Maths.Cos(ang), Maths.Sin(ang));
                        val *= rot;
                    }
                    tmp[q] = val;
                }
                var ch = IFFT_via_DHT_Custom(tmp);
                for (int l = 0; l < L; l++) Chalf[n0, l] = ch[l];
            }

            // 3) DFT_M по n0 через DHT и четверть-фаза
            var B = new Complex32[2 * N];
            float gain = Maths.Sqrt(M) / Maths.Sqrt(N);

            var colr = new float[M];
            var coli = new float[M];

            for (int l = 0; l < L; l++)
            {
                // main
                for (int n0 = 0; n0 < M; n0++) { colr[n0] = Cmain[n0, l].Real; coli[n0] = Cmain[n0, l].Imag; }
                var Ymain = DFT_via_DHT(colr, coli);

                // half
                for (int n0 = 0; n0 < M; n0++) { colr[n0] = Chalf[n0, l].Real; coli[n0] = Chalf[n0, l].Imag; }
                var Yhalf = DFT_via_DHT(colr, coli);

                for (int k = 0; k < M; k++)
                {
                    var phase = PhasePlusPiOver2(k);   // e^{+jπk/2}
                    var P = phase * Ymain[k];
                    var Q = phase * Yhalf[k];

                    int u = l * M + k;
                    B[u + 0] = P * gain;                 // main
                    B[u + N] = -Complex32.I * Q * gain;  // half (−j)
                }
            }

            return B;
        }

        /// <summary>
        /// Обратный FWHT через DHT (вход: 2N коэф. main+half).
        /// </summary>
        internal Complex32[] IFWHT_Hartley(Complex32[] B, FastWeylHeisenbergTransform.PolyphaseCache C)
        {
            int N = C.N;
            int M = C.M;
            int L = C.L;
            if (B.Length != 2 * N) throw new ArgumentException("Expect 2N coefficients (main + half)");

            var S_hat = C.S_hat;
            var T_hat = C.T_hat;

            // 1) Распаковка по k для каждого l: снимаем четверть-фазу и делаем IFFT_M (через DHT)
            var Cmain = new Complex32[M, L];
            var Chalf = new Complex32[M, L];

            var Y_main = new Complex32[M];
            var Y_half = new Complex32[M];

            float invGain = Maths.Sqrt(N) / Maths.Sqrt(M);

            for (int l = 0; l < L; l++)
            {
                for (int k = 0; k < M; k++)
                {
                    int u = l * M + k;
                    var bMain = B[u + 0];
                    var bHalf = B[u + N];

                    var P = bMain * invGain;
                    var Q = Complex32.I * bHalf * invGain; // т.к. вперёд был −j*Q*gain

                    var phase = PhasePlusPiOver2(k);
                    var phaseConj = new Complex32(phase.Real, -phase.Imag);

                    Y_main[k] = phaseConj * P;
                    Y_half[k] = phaseConj * Q;
                }

                // IFFT_M со стандартной 1/M нормой, затем *M (как в референсе)
                var cmain_col = IFFT_via_DHT_Standard(Y_main);
                var chalf_col = IFFT_via_DHT_Standard(Y_half);
                for (int n0 = 0; n0 < M; n0++)
                {
                    Cmain[n0, l] = cmain_col[n0] * M;
                    Chalf[n0, l] = chalf_col[n0] * M;
                }
            }

            // 2) Аджойнт частотных корреляций (по l)
            var Xhat = new Complex32[M, L];

            var xr = new float[L];
            var xi = new float[L];

            for (int n0 = 0; n0 < M; n0++)
            {
                // FFT_L{Cmain}
                for (int l = 0; l < L; l++) { xr[l] = Cmain[n0, l].Real; xi[l] = Cmain[n0, l].Imag; }
                var Y1 = DFT_via_DHT(xr, xi);

                // FFT_L{Chalf}
                for (int l = 0; l < L; l++) { xr[l] = Chalf[n0, l].Real; xi[l] = Chalf[n0, l].Imag; }
                var Y2 = DFT_via_DHT(xr, xi);

                bool carry = (n0 + (M >> 1)) >= M;
                for (int q = 0; q < L; q++)
                {
                    var sum = S_hat[n0, q] * Y1[q];
                    if (carry)
                    {
                        float ang = 2f * Maths.Pi * q / L; // exp(+j*2π q/L)
                        var rot = new Complex32(Maths.Cos(ang), Maths.Sin(ang));
                        sum += T_hat[n0, q] * (rot * Y2[q]);
                    }
                    else
                    {
                        sum += T_hat[n0, q] * Y2[q];
                    }
                    Xhat[n0, q] = sum;
                }
            }

            // 3) Обратная полифазная композиция: IFFT_L через DHT (sqrt(L) * IDFT_standard)
            var Arec = new Complex32[N];
            for (int n0 = 0; n0 < M; n0++)
            {
                var tmp = new Complex32[L];
                for (int q = 0; q < L; q++) tmp[q] = Xhat[n0, q];
                var time = IFFT_via_DHT_Custom(tmp);
                for (int r = 0; r < L; r++) Arec[r * M + n0] = time[r];
            }

            return Arec;
        }
        #endregion

        #region DFT via DHT helpers
        /// <summary>Комплексный unnormalized DFT через две DHT над Re/Im.</summary>
        private Complex32[] DFT_via_DHT(float[] xr, float[] xi)
        {
            int n = xr.Length;
            var Hx = FHT.Forward((float[])xr.Clone());
            var Hy = FHT.Forward((float[])xi.Clone());

            var Z = new Complex32[n];
            for (int k = 0; k < n; k++)
            {
                int k2 = (n - k) % n;
                float Xr = 0.5f * (Hx[k] + Hx[k2]);
                float Xi = 0.5f * (Hx[k2] - Hx[k]);
                float Yr = 0.5f * (Hy[k] + Hy[k2]);
                float Yi = 0.5f * (Hy[k2] - Hy[k]);
                Z[k] = new Complex32(Xr - Yi, Xi + Yr);
            }
            return Z;
        }

        /// <summary>Стандартная (1/n) комплексная IDFT через две DHT.</summary>
        private Complex32[] IFFT_via_DHT_Standard(Complex32[] Z)
        {
            int n = Z.Length;
            var Zr = new float[n];
            var Zi = new float[n];
            for (int k = 0; k < n; k++) { Zr[k] = Z[k].Real; Zi[k] = Z[k].Imag; }

            var Hx = new float[n];
            var Hy = new float[n];

            for (int k = 0; k <= n / 2; k++)
            {
                int k2 = (n - k) % n;
                float A = Zr[k] + Zi[k];
                float B = Zr[k] - Zi[k];
                float C = Zr[k2] + Zi[k2];
                float D = Zr[k2] - Zi[k2];

                Hx[k] = 0.5f * (B + C);
                Hx[k2] = 0.5f * (A + D);
                Hy[k] = 0.5f * (A - D);
                Hy[k2] = 0.5f * (C - B);
            }

            var x = FHT.Forward(Hx);
            var y = FHT.Forward(Hy);

            float invN = 1f / n;
            var time = new Complex32[n];
            for (int i = 0; i < n; i++)
                time[i] = new Complex32(x[i] * invN, y[i] * invN);

            return time;
        }

        /// <summary>
        /// Кастомная «IFFT» как в референсном коде: sqrt(n) * IDFT_standard.
        /// Совпадает с их FFT(true), где применяется 1/sqrt(n) после «unnorm» обратного ДПФ.
        /// </summary>
        private Complex32[] IFFT_via_DHT_Custom(Complex32[] Z)
        {
            int n = Z.Length;
            var std = IFFT_via_DHT_Standard(Z);
            float s = Maths.Sqrt(n);
            for (int i = 0; i < n; i++) std[i] *= s;
            return std;
        }

        /// <summary>Четверть-периодный множитель e^{+jπk/2}.</summary>
        private static Complex32 PhasePlusPiOver2(int k)
        {
            switch (k & 3)
            {
                case 0: return new Complex32(+1, 0);
                case 1: return new Complex32(0, +1);
                case 2: return new Complex32(-1, 0);
                default: return new Complex32(0, -1);
            }
        }
        #endregion
    }
}
