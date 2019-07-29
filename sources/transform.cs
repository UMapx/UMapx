// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;
using System.Threading.Tasks;
using UMapx.Core;
using UMapx.Window;

namespace UMapx.Transform
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                              UMAPX.TRANSFORM
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
    // **************************************************************************

    #region Orthogonal transforms
    /// <summary>
    /// Определяет преобразование Уолша-Адамара.
    /// <remarks>
    /// Дискретное преобразование Уолша-Адамара - это одно из первых дискретных ортогональных преобразований. Оно активно используется для сжатия сигналов, размерность которых
    /// равна степени 2 (например 64, 128, 256 и так далее).
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// http://kibia.ru/teachers/kreindelin/pdf/2.pdf
    /// </remarks>
    /// </summary>
    public class WalshHadamardTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Нормализованное преобразование или нет.
        /// </summary>
        private bool normalized;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="normalized">Нормализированное преобразование или нет</param>
        /// <param name="direction">Направление обработки</param>
        public WalshHadamardTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.normalized = normalized; this.direction = direction;
        }
        /// <summary>
        /// Нормализированное преобразование или нет.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.normalized;
            }
            set
            {
                this.normalized = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        #endregion

        #region Walsh-Hadamard static components
        /// <summary>
        /// Реализует построение матрицы Уолша-Адамара.
        /// </summary>
        /// <param name="powOf2">Степень двойки</param>
        /// <returns>Матрица</returns>
        public static double[,] Hadamard(int powOf2)
        {
            // Количество необходимых итераций:
            int iterations = powOf2 - 1;
            // Пораждающая матрица Адамара:
            double[,] hadamard = WalshHadamardTransform.Hadamard();

            if (iterations > 0)
            {
                double[,] H = WalshHadamardTransform.Hadamard();

                for (int i = 0; i < iterations; i++)
                {
                    H = Matrice.Kronecker(H, hadamard);
                }
                return H;
            }

            return hadamard;
        }
        /// <summary>
        /// Реализует построение матрицы Уолша-Адамара [2 x 2].
        /// </summary>
        /// <returns>Матрица</returns>
        public static double[,] Hadamard()
        {
            return (new double[2, 2] { { 1, 1 }, { 1, -1 } });
        }
        #endregion

        #region Walsh-Hadamard Transform
        /// <summary>
        /// Прямое дискретное преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            int N = A.Length;

            if (!Maths.IsPower(N, 2))
                throw new Exception("Размерность сигнала должна быть степенью 2");

            int n = (int)Maths.Log2(N);
            double[,] U = WalshHadamardTransform.Hadamard(n);
            double[] B = Matrice.Dot(A, U);

            if (normalized)
            {
                B = Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            int N = B.Length;

            if (!Maths.IsPower(N, 2))
                throw new Exception("Размерность сигнала должна быть степенью 2");

            int n = (int)Maths.Log2(N);
            double[,] U = WalshHadamardTransform.Hadamard(n);
            double[] A = Matrice.Dot(B, U.Transponate());

            if (normalized)
            {
                A = Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);

            if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                throw new Exception("Размерность сигнала должна быть степенью 2");

            double[,] U = WalshHadamardTransform.Hadamard((int)Maths.Log2(N));
            double[,] V = WalshHadamardTransform.Hadamard((int)Maths.Log2(M));
            double[,] B;

            if (direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Transponate());
                B = normalized ? B.Div(Math.Sqrt(N * M)) : B;
            }
            else if (direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = normalized ? B.Div(Math.Sqrt(N)) : B;
            }
            else
            {
                B = A.Dot(V.Transponate());
                B = normalized ? B.Div(Math.Sqrt(M)) : B;
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);

            if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                throw new Exception("Размерность сигнала должна быть степенью 2");

            double[,] U = WalshHadamardTransform.Hadamard((int)Maths.Log2(N));
            double[,] V = WalshHadamardTransform.Hadamard((int)Maths.Log2(M));
            double[,] A;

            if (direction == Direction.Both)
            {
                A = U.Transponate().Dot(B).Dot(V);
                A = normalized ? A.Div(Math.Sqrt(N * M)) : A;
            }
            else if (direction == Direction.Vertical)
            {
                A = U.Transponate().Dot(B);
                A = normalized ? A.Div(Math.Sqrt(N)) : A;
            }
            else
            {
                A = B.Dot(V);
                A = normalized ? A.Div(Math.Sqrt(M)) : A;
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;

            if (!Maths.IsPower(N, 2))
                throw new Exception("Размерность сигнала должна быть степенью 2");

            int n = (int)Maths.Log2(N);
            double[,] U = WalshHadamardTransform.Hadamard(n);
            Complex[] B = Matrice.Dot(A, U);

            if (normalized)
            {
                B = Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;

            if (!Maths.IsPower(N, 2))
                throw new Exception("Размерность сигнала должна быть степенью 2");

            int n = (int)Maths.Log2(N);
            double[,] U = WalshHadamardTransform.Hadamard(n);
            Complex[] A = Matrice.Dot(B, U.Transponate());

            if (normalized)
            {
                A = Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);

            if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                throw new Exception("Размерность сигнала должна быть степенью 2");

            double[,] U = WalshHadamardTransform.Hadamard((int)Maths.Log2(N));
            double[,] V = WalshHadamardTransform.Hadamard((int)Maths.Log2(M));
            Complex[,] B;

            if (direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Transponate());
                B = normalized ? B.Div(Math.Sqrt(N * M)) : B;
            }
            else if (direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = normalized ? B.Div(Math.Sqrt(N)) : B;
            }
            else
            {
                B = A.Dot(V.Transponate());
                B = normalized ? B.Div(Math.Sqrt(M)) : B;
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);

            if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                throw new Exception("Размерность сигнала должна быть степенью 2");

            double[,] U = WalshHadamardTransform.Hadamard((int)Maths.Log2(N));
            double[,] V = WalshHadamardTransform.Hadamard((int)Maths.Log2(M));
            Complex[,] A;

            if (direction == Direction.Both)
            {
                A = U.Transponate().Dot(B).Dot(V);
                A = normalized ? A.Div(Math.Sqrt(N * M)) : A;
            }
            else if (direction == Direction.Vertical)
            {
                A = U.Transponate().Dot(B);
                A = normalized ? A.Div(Math.Sqrt(N)) : A;
            }
            else
            {
                A = B.Dot(V);
                A = normalized ? A.Div(Math.Sqrt(M)) : A;
            }

            return A;
        }
        #endregion
    }
    /// <summary>
    /// Определяет дискретное косинусоидальное преобразование.
    /// <remarks>
    /// Дискретное консинусное преобразование (ДКП) тесно связано с дискретным преобразованием Фурье и является гомоморфизмом его векторного пространства.
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Discrete_cosine_transform
    /// </remarks>
    /// </summary>
    public class CosineTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует дискретное косинусоидальное преобразование.
        /// </summary>
        /// <param name="direction">Направление обработки</param>
        public CosineTransform(Direction direction = Direction.Vertical) 
        {
            this.direction = direction;
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        #endregion

        #region Cosine static components
        /// <summary>
        /// Реализует построение матрицы косинусного преобразования.
        /// </summary>
        /// <param name="n">Размер матрицы</param>
        /// <returns>Матрица</returns>
        public static double[,] Cosine(int n)
        {
            int j, i;
            double[,] H = new double[n, n];
            double c = Maths.Pi / (2.0 * n);
            double g1 = Maths.Sqrt(1.0 / n);
            double g2 = Math.Sqrt(2.0 / n);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = (i == 0) ? g1 : Math.Cos((2 * j + 1) * i * c) * g2;
                }
            }

            return (H);
        }
        #endregion

        #region Cosine Discrete Transform
        /// <summary>
        /// Прямое дискретное косинусоидальное преобразование.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            int N = A.Length;
            double[,] U = CosineTransform.Cosine(N);
            return Matrice.Dot(A, U);
        }
        /// <summary>
        /// Обратное дискретное косинусоидальное преобразование.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            int N = B.Length;
            double[,] U = CosineTransform.Cosine(N);
            return Matrice.Dot(B, U.Transponate());
        }
        /// <summary>
        /// Прямое дискретное косинусоидальное преобразование.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            double[,] U = CosineTransform.Cosine(N);
            double[,] V = CosineTransform.Cosine(M);

            if (direction == Direction.Both)
            {
                return U.Dot(A).Dot(V.Transponate());
            }
            else if (direction == Direction.Vertical)
            {
                return U.Dot(A);
            }
            return A.Dot(V.Transponate());
        }
        /// <summary>
        /// Обратное дискретное косинусоидальное преобразование.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            double[,] U = CosineTransform.Cosine(N);
            double[,] V = CosineTransform.Cosine(M);

            if (direction == Direction.Both)
            {
                return U.Transponate().Dot(B).Dot(V);
            }
            else if (direction == Direction.Vertical)
            {
                return U.Transponate().Dot(B);
            }
            return B.Dot(V);
        }
        /// <summary>
        /// Прямое дискретное косинусоидальное преобразование.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            double[,] U = CosineTransform.Cosine(N);
            return Matrice.Dot(A, U);
        }
        /// <summary>
        /// Обратное дискретное косинусоидальное преобразование.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;
            double[,] U = CosineTransform.Cosine(N);
            return Matrice.Dot(B, U.Transponate());
        }
        /// <summary>
        /// Прямое дискретное косинусоидальное преобразование.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            double[,] U = CosineTransform.Cosine(N);
            double[,] V = CosineTransform.Cosine(M);

            if (direction == Direction.Both)
            {
                return U.Dot(A).Dot(V.Transponate());
            }
            else if (direction == Direction.Vertical)
            {
                return U.Dot(A);
            }
            return A.Dot(V.Transponate());
        }
        /// <summary>
        /// Обратное дискретное косинусоидальное преобразование.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            double[,] U = CosineTransform.Cosine(N);
            double[,] V = CosineTransform.Cosine(M);

            if (direction == Direction.Both)
            {
                return U.Transponate().Dot(B).Dot(V);
            }
            else if (direction == Direction.Vertical)
            {
                return U.Transponate().Dot(B);
            }
            return B.Dot(V);
        }
        #endregion
    }
    /// <summary>
    /// Определяет дискретное синусоидальное преобразование.
    /// <remarks>
    /// Дискретное ортогональное синусное преобразование было предложено Джейном в качестве аппроксимации преобразования Карунена-Лоэва для марковского процесса.
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// http://sernam.ru/book_prett1.php?id=91
    /// </remarks>
    /// </summary>
    public class SineTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует дискретное синусоидальное преобразование.
        /// </summary>
        /// <param name="direction">Направление обработки</param>
        public SineTransform(Direction direction = Direction.Vertical) 
        {
            this.direction = direction;
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        #endregion

        #region Sine static components
        /// <summary>
        /// Реализует построение матрицы синусного преобразования.
        /// </summary>
        /// <param name="n">Размер матрицы</param>
        /// <returns>Матрица</returns>
        public static double[,] Sine(int n)
        {
            int j, i;
            double[,] H = new double[n, n];
            double n1 = n + 1;
            double scale = Maths.Sqrt(2.0 / n1);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = Math.Sin(Maths.Pi * (j + 1) * (i + 1) / n1) * scale;
                }
            }

            return (H);
        }
        #endregion

        #region Sine Discrete Transform
        /// <summary>
        /// Прямое дискретное синусоидальное преобразование.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            int N = A.Length;
            double[,] U = SineTransform.Sine(N);
            return Matrice.Dot(A, U);
        }
        /// <summary>
        /// Обратное дискретное синусоидальное преобразование.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            int N = B.Length;
            double[,] U = SineTransform.Sine(N);
            return Matrice.Dot(B, U.Transponate());
        }
        /// <summary>
        /// Прямое дискретное синусоидальное преобразование.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            double[,] U = SineTransform.Sine(N);
            double[,] V = SineTransform.Sine(M);

            if (direction == Direction.Both)
            {
                return U.Dot(A).Dot(V.Transponate());
            }
            else if (direction == Direction.Vertical)
            {
                return U.Dot(A);
            }
            return A.Dot(V.Transponate());
        }
        /// <summary>
        /// Обратное дискретное синусоидальное преобразование.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            double[,] U = SineTransform.Sine(N);
            double[,] V = SineTransform.Sine(M);

            if (direction == Direction.Both)
            {
                return U.Transponate().Dot(B).Dot(V);
            }
            else if (direction == Direction.Vertical)
            {
                return U.Transponate().Dot(B);
            }
            return B.Dot(V);
        }
        /// <summary>
        /// Прямое дискретное синусоидальное преобразование.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            double[,] U = SineTransform.Sine(N);
            return Matrice.Dot(A, U);
        }
        /// <summary>
        /// Обратное дискретное синусоидальное преобразование.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;
            double[,] U = SineTransform.Sine(N);
            return Matrice.Dot(B, U.Transponate());
        }
        /// <summary>
        /// Прямое дискретное синусоидальное преобразование.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            double[,] U = SineTransform.Sine(N);
            double[,] V = SineTransform.Sine(M);

            if (direction == Direction.Both)
            {
                return U.Dot(A).Dot(V.Transponate());
            }
            else if (direction == Direction.Vertical)
            {
                return U.Dot(A);
            }
            return A.Dot(V.Transponate());
        }
        /// <summary>
        /// Обратное дискретное синусоидальное преобразование.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            double[,] U = SineTransform.Sine(N);
            double[,] V = SineTransform.Sine(M);

            if (direction == Direction.Both)
            {
                return U.Transponate().Dot(B).Dot(V);
            }
            else if (direction == Direction.Vertical)
            {
                return U.Transponate().Dot(B);
            }
            return B.Dot(V);
        }
        #endregion
    }
    /// <summary>
    /// Определяет дискретное преобразование Хартли.
    /// <remarks>
    /// Данное дискретное ортогональное преобразование служит своего рода заменой дискретного преобразования Фурье для вещественных сигналов, однако, может применяться
    /// и для анализа комплексных.
    /// Кроме того, в данном программном модуле реализовано быстрое преобразование Хартли (БПХ), вычисляемое через быстрое преобразование Фурье (алгоритм Кули-Тьюки).
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Discrete_Hartley_transform
    /// </remarks>
    /// </summary>
    public class HartleyTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Нормализованное преобразование или нет.
        /// </summary>
        private bool normalized;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует дискретное преобразование Хартли.
        /// </summary>
        /// <param name="normalized">Нормализованное преобразование или нет</param>
        /// <param name="direction">Направление обработки</param>
        public HartleyTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.normalized = normalized; this.direction = direction;
        }
        /// <summary>
        /// Нормализированное преобразование или нет.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.normalized;
            }
            set
            {
                this.normalized = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        #endregion

        #region Hartley static components
        /// <summary>
        /// Реализует построение матрицы преобразования Хартли.
        /// </summary>
        /// <param name="n">Размер матрицы</param>
        /// <returns>Матрица</returns>
        public static double[,] Hartley(int n)
        {
            int j, i;
            double[,] H = new double[n, n];
            double s = (2.0f * Maths.Pi) / n;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = Special.Cas(s * i * j);
                }
            }

            return (H);
        }
        #endregion

        #region Hartley Discrete Transform
        /// <summary>
        /// Прямое дискретное преобразование Хартли.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            int N = A.Length;
            double[,] U = HartleyTransform.Hartley(N);
            double[] B = Matrice.Dot(A, U);

            if (normalized)
            {
                B = Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Хартли.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            int N = B.Length;
            double[,] U = HartleyTransform.Hartley(N);
            double[] A = Matrice.Dot(B, U.Transponate());

            if (normalized)
            {
                A = Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Хартли.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            double[,] U = HartleyTransform.Hartley(N);
            double[,] V = HartleyTransform.Hartley(M);
            double[,] B;

            if (direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Transponate());
                B = normalized ? B.Div(Math.Sqrt(N * M)) : B;
            }
            else if (direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = normalized ? B.Div(Math.Sqrt(N)) : B;
            }
            else
            {
                B = A.Dot(V.Transponate());
                B = normalized ? B.Div(Math.Sqrt(M)) : B;
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Хартли.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            double[,] U = HartleyTransform.Hartley(N);
            double[,] V = HartleyTransform.Hartley(M);
            double[,] A;

            if (direction == Direction.Both)
            {
                A = U.Transponate().Dot(B).Dot(V);
                A = normalized ? A.Div(Math.Sqrt(N * M)) : A;
            }
            else if (direction == Direction.Vertical)
            {
                A = U.Transponate().Dot(B);
                A = normalized ? A.Div(Math.Sqrt(N)) : A;
            }
            else
            {
                A = B.Dot(V);
                A = normalized ? A.Div(Math.Sqrt(M)) : A;
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Хартли.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            double[,] U = HartleyTransform.Hartley(N);
            Complex[] B = Matrice.Dot(A, U);

            if (normalized)
            {
                B = Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Хартли.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;
            double[,] U = HartleyTransform.Hartley(N);
            Complex[] A = Matrice.Dot(B, U.Transponate());

            if (normalized)
            {
                A = Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Хартли.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            double[,] U = HartleyTransform.Hartley(N);
            double[,] V = HartleyTransform.Hartley(M);
            Complex[,] B;

            if (direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Transponate());
                B = normalized ? B.Div(Math.Sqrt(N * M)) : B;
            }
            else if (direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = normalized ? B.Div(Math.Sqrt(N)) : B;
            }
            else
            {
                B = A.Dot(V.Transponate());
                B = normalized ? B.Div(Math.Sqrt(M)) : B;
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Хартли.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            double[,] U = HartleyTransform.Hartley(N);
            double[,] V = HartleyTransform.Hartley(M);
            Complex[,] A;

            if (direction == Direction.Both)
            {
                A = U.Transponate().Dot(B).Dot(V);
                A = normalized ? A.Div(Math.Sqrt(N * M)) : A;
            }
            else if (direction == Direction.Vertical)
            {
                A = U.Transponate().Dot(B);
                A = normalized ? A.Div(Math.Sqrt(N)) : A;
            }
            else
            {
                A = B.Dot(V);
                A = normalized ? A.Div(Math.Sqrt(M)) : A;
            }

            return A;
        }
        #endregion
    }
    /// <summary>
    /// Определяет дискретное преобразование Фурье.
    /// <remarks>
    /// Это одно из преобразований Фурье, широко применяемых в алгоритмах цифровой обработки сигналов. Дискретные преобразования Фурье помогают решать дифференциальные уравнения 
    /// в частных производных и выполнять такие операции, как свёртки. Дискретные преобразования Фурье также активно используются в статистике, при анализе временных рядов.
    /// Кроме того, в данном программном модуле реализовано быстрое преобразование Фурье (алгоритм Кули-Тьюки).
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Discrete_Fourier_transform
    /// </remarks>
    /// </summary>
    public class FourierTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Нормализованное преобразование или нет.
        /// </summary>
        private bool normalized;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует дискретное преобразование Фурье.
        /// </summary>
        /// <param name="normalized">Нормализированное преобразование или нет</param>
        /// <param name="direction">Направление обработки</param>
        public FourierTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.normalized = normalized; this.direction = direction;
        }
        /// <summary>
        /// Нормализированное преобразование или нет.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.normalized;
            }
            set
            {
                this.normalized = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        #endregion

        #region Fourier static components
        /// <summary>
        /// Реализует построение матрицы Фурье.
        /// </summary>
        /// <param name="n">Размер матрицы</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Fourier(int n)
        {
            Complex[,] H = new Complex[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = Maths.Exp(-2 * Maths.Pi * Maths.I * i * j / n);
                }
            }
            return H;
        }
        #endregion

        #region Fourier Discrete Transform
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            Complex[,] U = FourierTransform.Fourier(N);
            Complex[] B = Matrice.Dot(A, U);

            if (normalized)
            {
                B = Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;
            Complex[,] U = FourierTransform.Fourier(N);
            Complex[] A = Matrice.Dot(B, U.Hermitian());

            if (normalized)
            {
                A = Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            Complex[,] U = FourierTransform.Fourier(N);
            Complex[,] V = FourierTransform.Fourier(M);
            Complex[,] B;

            if (direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Hermitian());
                B = normalized ? B.Div(Math.Sqrt(N * M)) : B;
            }
            else if (direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = normalized ? B.Div(Math.Sqrt(N)) : B;
            }
            else
            {
                B = A.Dot(V.Hermitian());
                B = normalized ? B.Div(Math.Sqrt(M)) : B;
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            Complex[,] U = FourierTransform.Fourier(N);
            Complex[,] V = FourierTransform.Fourier(M);
            Complex[,] A;

            if (direction == Direction.Both)
            {
                A = U.Hermitian().Dot(B).Dot(V);
                A = normalized ? A.Div(Math.Sqrt(N * M)) : A;
            }
            else if (direction == Direction.Vertical)
            {
                A = U.Hermitian().Dot(B);
                A = normalized ? A.Div(Math.Sqrt(N)) : A;
            }
            else
            {
                A = B.Dot(V);
                A = normalized ? A.Div(Math.Sqrt(M)) : A;
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
    /// <summary>
    /// Определяет дискретное преобразование Лапласа.
    /// <remarks>
    /// Данный класс реализует дискретный эквивалент непрерывного преобразования Лапласа.
    /// Непрерывное преобразование Лапласа тесно связано с преобразованием Фурье. Иными словами, преобразование Фурье эквивалентно двустороннему преобразованию Лапласа с комплексным 
    /// аргументом: s = iω. Эта связь между преобразованиями часто используется для того, чтобы определить частотный спектр сигнала.
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Laplace_transform
    /// </remarks>
    /// </summary>
    public class LaplaceTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Среднеквадратическое отклонение.
        /// </summary>
        private double sigma;
        /// <summary>
        /// Нормализованное преобразование или нет.
        /// </summary>
        private bool normalized;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует дискретное преобразование Лапласа.
        /// </summary>
        /// <param name="sigma">Среднеквадратическое отклонение (0, 1)</param>
        /// <param name="normalized">Нормализированное преобразование или нет</param>
        /// <param name="direction">Направление обработки</param>
        public LaplaceTransform(double sigma = 0.0005, bool normalized = true, Direction direction = Direction.Vertical)
        {
            Sigma = sigma; this.normalized = normalized; this.direction = direction;
        }
        /// <summary>
        /// Получает или задает значение среднеквадратического отклонения (0, 1).
        /// <remarks>
        /// В случае, если σ = 0, то преобразование Лапласа принимает вид преобразования Фурье.
        /// </remarks>
        /// </summary>
        public double Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Нормализированное преобразование или нет.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.normalized;
            }
            set
            {
                this.normalized = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        #endregion

        #region Laplace static components
        /// <summary>
        /// Реализует построение матрицы Лапласа.
        /// </summary>
        /// <param name="n">Размер матрицы</param>
        /// <param name="sigma">Среднеквадратическое отклонение (0, 1)</param>
        /// <param name="backward">Возвращать матрицу обратного преобразования или нет</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Laplace(int n, double sigma, bool backward = false)
        {
            Complex[,] H = new Complex[n, n];
            double factor;
            int i, j;

            // inverse matrix or not?
            if (backward)
            {
                for (i = 0; i < n; i++)
                {
                    factor = Maths.Exp(-sigma * i);

                    for (j = 0; j < n; j++)
                    {
                        H[i, j] = Maths.Exp(-2 * Maths.Pi * Maths.I * j / n * i) / factor;
                    }
                }
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    factor = Maths.Exp(-sigma * i);

                    for (j = 0; j < n; j++)
                    {
                        H[i, j] = factor * Maths.Exp(-2 * Maths.Pi * Maths.I * j / n * i);
                    }
                }
            }
            return H;
        }
        #endregion

        #region Laplace Discrete Transform
        /// <summary>
        /// Прямое дискретное преобразование Лапласа.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            Complex[,] U = LaplaceTransform.Laplace(N, sigma);
            Complex[] B = Matrice.Dot(A, U);

            if (normalized)
            {
                B = Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Лапласа.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;
            Complex[,] U = LaplaceTransform.Laplace(N, sigma, true);
            Complex[] A = Matrice.Dot(B, U.Hermitian());

            if (normalized)
            {
                A = Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Лапласа.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            Complex[,] U = LaplaceTransform.Laplace(N, sigma);
            Complex[,] V = LaplaceTransform.Laplace(M, sigma);
            Complex[,] B;

            if (direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Hermitian());
                B = normalized ? B.Div(Math.Sqrt(N * M)) : B;
            }
            else if (direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = normalized ? B.Div(Math.Sqrt(N)) : B;
            }
            else
            {
                B = A.Dot(V.Hermitian());
                B = normalized ? B.Div(Math.Sqrt(M)) : B;
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Лапласа.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            Complex[,] U = LaplaceTransform.Laplace(N, sigma, true);
            Complex[,] V = LaplaceTransform.Laplace(M, sigma, true);
            Complex[,] A;

            if (direction == Direction.Both)
            {
                A = U.Hermitian().Dot(B).Dot(V);
                A = normalized ? A.Div(Math.Sqrt(N * M)) : A;
            }
            else if (direction == Direction.Vertical)
            {
                A = U.Hermitian().Dot(B);
                A = normalized ? A.Div(Math.Sqrt(N)) : A;
            }
            else
            {
                A = B.Dot( V);
                A = normalized ? A.Div(Math.Sqrt(M)) : A;
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
    /// <summary>
    /// Определяет преобразование Гильберта.
    /// <remarks>
    /// Суть преобразования заключается в том, что при прямом преобразовании Гильберта происходит усиление положительных частот и обнуление отрицательных.
    /// В свою очередь, обратное преобразование Гильберта вычисляется путем отображения прямого преобразования: H^–1{h(t)} = –H{h(t)}.
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Hilbert_transform
    /// </remarks>
    /// </summary>
    public class HilbertTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Преобразование Фурье.
        /// </summary>
        private FourierTransform FFT;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует преобразование Гильберта.
        /// </summary>
        /// <param name="normalized">Нормализированное преобразование или нет</param>
        /// <param name="direction">Направление обработки</param>
        public HilbertTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.FFT = new FourierTransform(normalized, direction);
        }
        /// <summary>
        /// Нормализированное преобразование или нет.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.FFT.Normalized;
            }
            set
            {
                this.FFT.Normalized = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.FFT.Direction;
            }
            set
            {
                this.FFT.Direction = value;
            }
        }
        #endregion

        #region Hilbert Discrete Transform
        /// <summary>
        /// Прямое дискретное преобразование Гильберта.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;

            // Перевод массива вещественных чисел в массив комплексных
            // и быстрое преобразование Фурье
            Complex[] F = FFT.Forward(A);

            // Перегруппировка:
            HilbertTransform.hilbertf(F, N);

            // Обратное преобразование Фурье:
            F = FFT.Backward(F);

            // Группировка:
            return HilbertTransform.hilbertb(A, F, N);
        }
        /// <summary>
        /// Обратное дискретное преобразование Гильберта.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            // Обратное преобразование Гильберта вычисляется путем
            // отображения прямого преобразования:
            // H^–1{h(t)} = –H{h(t)}

            int N = B.Length, i;
            Complex[] A = new Complex[N];

            for (i = 0; i < N; i++)
            {
                A[i] = new Complex(B[i].Re, 0);
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Гильберта.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            Complex[,] B = (Complex[,])A.Clone();
            int N = A.GetLength(0);
            int M = A.GetLength(1);

            if (Direction == Direction.Both)
            {
                // 2-dimension horizontal Hilbert transform:
                FFT.Direction = Direction.Horizontal;
                B = FFT.Forward(B);

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    HilbertTransform.hilbertf(row, M);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                B = FFT.Backward(B);

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    Complex[] num = new Complex[M];

                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                        num[j] = A[i, j];
                    }

                    num = HilbertTransform.hilbertb(num, row, M);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = num[j];
                    }
                });


                // 2-dimension vertical Hilbert transform:
                FFT.Direction = Direction.Vertical;
                B = FFT.Forward(B);

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    HilbertTransform.hilbertf(col, N);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });

                B = FFT.Backward(B);

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    Complex[] num = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                        num[i] = A[i, j];
                    }

                    num = HilbertTransform.hilbertb(num, col, N);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = num[i];
                    }
                });

                FFT.Direction = Direction.Both;
            }
            else if (Direction == Direction.Vertical)
            {
                // 2-dimension vertical Hilbert transform:
                B = FFT.Forward(B);

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    HilbertTransform.hilbertf(col, N);

                    for (i = 0; i < M; i++)
                    {
                        B[i, j] = col[i];
                    }
                });

                B = FFT.Backward(B);

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    Complex[] num = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                        num[i] = A[i, j];
                    }

                    num = HilbertTransform.hilbertb(num, col, N);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = num[i];
                    }
                });
            }
            else
            {
                // 2-dimension horizontal Hilbert transform:
                B = FFT.Forward(B);

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    HilbertTransform.hilbertf(row, M);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                B = FFT.Backward(B);

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    Complex[] num = new Complex[M];

                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                        num[j] = A[i, j];
                    }

                    num = HilbertTransform.hilbertb(num, row, M);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = num[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Гильберта.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = (Complex[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (Direction == Direction.Both)
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                }
                );
            }
            else if (Direction == Direction.Vertical)
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                }
                );
            }
            else
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Реализует перегруппировку спектра по Гильберту.
        /// </summary>
        /// <param name="f">Спектр</param>
        /// <param name="n">Размерность</param>
        internal static void hilbertf(Complex[] f, int n)
        {
            // Усиление положительных частот и обнуление отрицательных:
            int n2 = n / 2;

            for (int i = 0; i < n2; i++)
            {
                f[i     ] *= 2.0;
                f[i + n2]  = Complex.Zero;
            }
            return;
        }
        /// <summary>
        /// Реализует группировку по Гильберту.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="f">Спектр, упорядоченный по Гильберту</param>
        /// <param name="n">Размерность</param>
        /// <returns>Одномерный массив</returns>
        internal static Complex[] hilbertb(Complex[] a, Complex[] f, int n)
        {
            Complex[] B = new Complex[n];

            // Обратное преобразование в массив вещественных чисел:
            for (int i = 0; i < n; i++)
            {
                B[i] = new Complex(a[i].Re, f[i].Im);
            }

            return B;
        }
        #endregion
    }
    #endregion

    #region Fast orthogonal transforms
    /// <summary>
    /// Определяет быстрое преобразование Уолша-Адамара.
    /// <remarks>
    /// Данная оптимизация преобразования Уолша-Адамара сложностью O(N^2) позволяет произвести вычисления за O(Nlog(N)). Алгоритм быстрого преобразования 
    /// Уолша-Адамара имеет схожую структуру с быстрым преобразованием Фурье (алгоритм Кули-Тьюки).
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// http://www.mathworks.com/matlabcentral/fileexchange/6879-fast-walsh-hadamard-transform
    /// </remarks>
    /// </summary>
    public class FastWalshHadamardTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        /// <summary>
        /// Нормализованное преобразование или нет.
        /// </summary>
        private bool normalized;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует быстрое преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="normalized">Нормализированное преобразование или нет</param>
        /// <param name="direction">Направление обработки</param>
        public FastWalshHadamardTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.normalized = normalized; this.direction = direction;
        }
        /// <summary>
        /// Нормализированное преобразование или нет.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.normalized;
            }
            set
            {
                this.normalized = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        #endregion

        #region Fast Walsh-Hadamard Discrete Transform
        /// <summary>
        /// Прямое быстрое преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            int N = A.Length;
            if (!Maths.IsPower(N, 2))
                throw new Exception("Размер сигнала должен быть равен степени 2.");

            double[] B = (double[])A.Clone();
            fwht(B);

            if (normalized)
            {
                B = Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Обратное быстрое преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            int N = B.Length;
            if (!Maths.IsPower(N, 2))
                throw new Exception("Размер сигнала должен быть равен степени 2.");

            double[] A = (double[])B.Clone();
            fwht(B);

            if (normalized)
            {
                A = Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Прямое быстрое преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            double[,] B = (double[,])A.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    fwht(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                }
                );

                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    fwht(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    B = Matrice.Div(B, Math.Sqrt(N * M));
                }
            }
            else if (direction == Direction.Vertical)
            {

                if (!Maths.IsPower(N, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    fwht(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    B = Matrice.Div(B, Math.Sqrt(N));
                }
            }
            else
            {
                if (!Maths.IsPower(M, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    fwht(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                if (normalized == true)
                {
                    B = Matrice.Div(B, Math.Sqrt(M));
                }
            }

            return B;
        }
        /// <summary>
        /// Обратное быстрое преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            double[,] A = (double[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    fwht(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    fwht(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                }
                );

                if (normalized == true)
                {
                    A = Matrice.Div(A, Math.Sqrt(N * M));
                }
            }
            else if (direction == Direction.Vertical)
            {
                if (!Maths.IsPower(N, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    fwht(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    A = Matrice.Div(A, Math.Sqrt(N));
                }
            }
            else
            {
                if (!Maths.IsPower(M, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    fwht(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });

                if (normalized == true)
                {
                    A = Matrice.Div(A, Math.Sqrt(M));
                }
            }

            return A;
        }
        /// <summary>
        /// Прямое быстрое преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            if (!Maths.IsPower(N, 2))
                throw new Exception("Размер сигнала должен быть равен степени 2.");

            Complex[] B = (Complex[])A.Clone();
            fwht(B);

            if (normalized)
            {
                B = Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Обратное быстрое преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;
            if (!Maths.IsPower(N, 2))
                throw new Exception("Размер сигнала должен быть равен степени 2.");

            Complex[] A = (Complex[])B.Clone();
            fwht(B);

            if (normalized)
            {
                A = Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Прямое быстрое преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            Complex[,] B = (Complex[,])A.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    fwht(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                }
                );

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    fwht(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    B = Matrice.Div(B, Math.Sqrt(N * M));
                }
            }
            else if (direction == Direction.Vertical)
            {
                if (!Maths.IsPower(N, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    fwht(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    B = Matrice.Div(B, Math.Sqrt(N));
                }
            }
            else
            {
                if (!Maths.IsPower(M, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    fwht(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                if (normalized == true)
                {
                    B = Matrice.Div(B, Math.Sqrt(M));
                }
            }

            return B;
        }
        /// <summary>
        /// Обратное быстрое преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = (Complex[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    fwht(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    fwht(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                }
                );

                if (normalized == true)
                {
                    A = Matrice.Div(A, Math.Sqrt(N * M));
                }
            }
            else if (direction == Direction.Vertical)
            {
                if (!Maths.IsPower(N, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    fwht(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    A = Matrice.Div(A, Math.Sqrt(N));
                }
            }
            else
            {
                if (!Maths.IsPower(M, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    fwht(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });

                if (normalized == true)
                {
                    A = Matrice.Div(A, Math.Sqrt(M));
                }
            }

            return A;
        }
        #endregion

        #region Private data
        /// <summary>
        /// Реализует быстрое преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        private void fwht(double[] data)
        {
            int N = data.Length;
            int log2N = (int)Maths.Log2(N);
            double x_even, x_odd;

            int k0 = N, k1 = 1, k2 = N / 2;
            int x, y, z, i, j, l;

            for (x = 0; x < log2N; x++)
            {
                l = 0;

                for (y = 0; y < k1; y++, l += k0)
                {
                    for (z = 0; z < k2; z++)
                    {
                        i = z + l; j = i + k2;

                        x_even = data[i];
                        x_odd  = data[j];

                        data[i] = x_even + x_odd;
                        data[j] = x_even - x_odd;
                    }
                }

                k0 /= 2; k1 *= 2; k2 /= 2;
            }
            return;
        }
        /// <summary>
        /// Реализует быстрое преобразование Уолша-Адамара.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        private void fwht(Complex[] data)
        {
            int N = data.Length;
            int log2N = (int)Maths.Log2(N);
            Complex x_even, x_odd;

            int k0 = N, k1 = 1, k2 = N / 2;
            int x, y, z, i, j, l;

            for (x = 0; x < log2N; x++)
            {
                l = 0;

                for (y = 0; y < k1; y++, l += k0)
                {
                    for (z = 0; z < k2; z++)
                    {
                        i = z + l; j = i + k2;

                        x_even = data[i];
                        x_odd  = data[j];

                        data[i] = x_even + x_odd;
                        data[j] = x_even - x_odd;
                    }
                }

                k0 /= 2; k1 *= 2; k2 /= 2;
            }
            return;
        }
        #endregion
    }
    /// <summary>
    /// Определяет быстрое преобразование Фурье по алгоритму Кули-Тюки.
    /// <remarks>
    /// Данная оптимизация алгоритма Кули-Тьюки для вычисления преобразования Фурье предназначена для анализа крупных сигналов.
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
    /// </remarks>
    /// </summary>
    public class FastFourierTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        /// <summary>
        /// Нормализованное преобразование или нет.
        /// </summary>
        private bool normalized;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует быстрое преобразование Фурье по алгоритму Кули-Тюки.
        /// </summary>
        /// <param name="normalized">Нормализированное преобразование или нет</param>
        /// <param name="direction">Направление обработки</param>
        public FastFourierTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.normalized = normalized; this.direction = direction;
        }
        /// <summary>
        /// Нормализированное преобразование или нет.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.normalized;
            }
            set
            {
                this.normalized = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        #endregion

        #region Fast Fourier Transform
        /// <summary>
        /// Прямое дискретное быстрое преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            if (!Maths.IsPower(N, 2))
                throw new Exception("Размерность сигнала должна быть степенью 2");

            Complex[] B = (Complex[])A.Clone();
            fft(B);

            if (normalized == true)
            {
                B = Matrice.Div(B, Math.Sqrt(N));
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное быстрое преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;
            if (!Maths.IsPower(N, 2))
                throw new Exception("Размерность сигнала должна быть степенью 2");

            Complex[] A = (Complex[])B.Clone();
            ifft(A);

            if (normalized == true)
            {
                A = Matrice.Div(A, Math.Sqrt(N));
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное быстрое преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            Complex[,] B = (Complex[,])A.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    ifft(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                }
                );

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    fft(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    B = Matrice.Div(B, Math.Sqrt(N * M));
                }
            }
            else if (direction == Direction.Vertical)
            {
                if (!Maths.IsPower(N, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    fft(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    B = Matrice.Div(B, Math.Sqrt(N));
                }
            }
            else
            {
                if (!Maths.IsPower(M, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    ifft(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                if (normalized == true)
                {
                    B = Matrice.Div(B, Math.Sqrt(M));
                }
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное быстрое преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = (Complex[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                if (!Maths.IsPower(N, 2) || !Maths.IsPower(M, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    ifft(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    fft(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                }
                );

                if (normalized == true)
                {
                    A = Matrice.Div(A, Math.Sqrt(N * M));
                }
            }
            else if (direction == Direction.Vertical)
            {
                if (!Maths.IsPower(N, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    ifft(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                if (normalized == true)
                {
                    A = Matrice.Div(A, Math.Sqrt(N));
                }
            }
            else
            {
                if (!Maths.IsPower(M, 2))
                    throw new Exception("Размерность сигнала должна быть степенью 2");

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    fft(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });

                if (normalized == true)
                {
                    A = Matrice.Div(A, Math.Sqrt(M));
                }
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Private data
        private const int minLength = 2;
        private const int maxLength = 16384;
        private const int minBits = 1;
        private const int maxBits = 14;
        private static int[][] reversedBits = new int[maxBits][];
        private static Complex[,][] complexRotation = new Complex[maxBits, 2][];
        #endregion

        #region Private voids
        /// <summary>
        /// Прямое дискретное быстрое преобразование Фурье.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        private static void fft(Complex[] data)
        {
            int n = data.Length;
            int m = Log2(n);

            // reorder data first
            FastFourierTransform.ReorderData(data);

            // compute FFT
            int tn = 1, tm, k, i, even, odd;
            Complex[] rotation;
            Complex t, ce, co;
            double tr, ti;

            for (k = 1; k <= m; k++)
            {
                rotation = FastFourierTransform.ForwardComplexRotation(k);
                tm = tn; tn <<= 1;

                for (i = 0; i < tm; i++)
                {
                    t = rotation[i];

                    for (even = i; even < n; even += tn)
                    {
                        odd = even + tm;
                        ce = data[even];
                        co = data[odd];

                        tr = co.Re * t.Re - co.Im * t.Im;
                        ti = co.Re * t.Im + co.Im * t.Re;

                        data[even].Re += tr;
                        data[even].Im += ti;

                        data[odd].Re = ce.Re - tr;
                        data[odd].Im = ce.Im - ti;
                    }
                }
            }
        }
        /// <summary>
        /// Обратное дискретное быстрое преобразование Фурье.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        private static void ifft(Complex[] data)
        {
            int n = data.Length;
            int m = Log2(n);

            // reorder data first
            FastFourierTransform.ReorderData(data);

            // compute FFT
            int tn = 1, tm, k, i, even, odd;
            Complex[] rotation;
            Complex t, ce, co;
            double tr, ti;

            for (k = 1; k <= m; k++)
            {
                rotation = FastFourierTransform.BackwardComplexRotation(k);
                tm = tn; tn <<= 1;

                for (i = 0; i < tm; i++)
                {
                    t = rotation[i];

                    for (even = i; even < n; even += tn)
                    {
                        odd = even + tm;
                        ce = data[even];
                        co = data[odd];

                        tr = co.Re * t.Re - co.Im * t.Im;
                        ti = co.Re * t.Im + co.Im * t.Re;

                        data[even].Re += tr;
                        data[even].Im += ti;

                        data[odd].Re = ce.Re - tr;
                        data[odd].Im = ce.Im - ti;
                    }
                }
            }
        }
        /// <summary>
        /// Получает массив с указателями на члены данных, которые должны быть заменены перед БПФ.
        /// </summary>
        /// <param name="numberOfBits">Количество битов</param>
        /// <returns>Массив</returns>
        private static int[] GetReversedBits(int numberOfBits)
        {
            if ((numberOfBits < minBits) || (numberOfBits > maxBits))
                throw new ArgumentOutOfRangeException();

            // check if the array is already calculated
            if (reversedBits[numberOfBits - 1] == null)
            {
                int n = Pow2(numberOfBits);
                int[] rBits = new int[n];
                int i, j, oldBits, newBits;

                // calculate the array
                for (i = 0; i < n; i++)
                {
                    oldBits = i;
                    newBits = 0;

                    for (j = 0; j < numberOfBits; j++)
                    {
                        newBits = (newBits << 1) | (oldBits & 1);
                        oldBits = (oldBits >> 1);
                    }
                    rBits[i] = newBits;
                }
                reversedBits[numberOfBits - 1] = rBits;
            }
            return reversedBits[numberOfBits - 1];
        }
        /// <summary>
        /// Получает прямое вращение комплексного числа.
        /// </summary>
        /// <param name="numberOfBits">Количество битов</param>
        /// <returns>Массив</returns>
        private static Complex[] ForwardComplexRotation(int numberOfBits)
        {
            int directionIndex = 0;

            // check if the array is already calculated
            if (complexRotation[numberOfBits - 1, directionIndex] == null)
            {
                int n = 1 << (numberOfBits - 1), i;
                double uR = 1.0;
                double uI = 0.0;
                double angle = -System.Math.PI / n;
                double wR = System.Math.Cos(angle);
                double wI = System.Math.Sin(angle);
                double t;
                Complex[] rotation = new Complex[n];

                for (i = 0; i < n; i++)
                {
                    rotation[i] = new Complex(uR, uI);
                    t = uR * wI + uI * wR;
                    uR = uR * wR - uI * wI;
                    uI = t;
                }

                complexRotation[numberOfBits - 1, directionIndex] = rotation;
            }
            return complexRotation[numberOfBits - 1, directionIndex];
        }
        /// <summary>
        /// Получает обратное вращение комплексного числа.
        /// </summary>
        /// <param name="numberOfBits">Количество битов</param>
        /// <returns>Массив</returns>
        private static Complex[] BackwardComplexRotation(int numberOfBits)
        {
            int directionIndex = 1;

            // check if the array is already calculated
            if (complexRotation[numberOfBits - 1, directionIndex] == null)
            {
                int n = 1 << (numberOfBits - 1), i;
                double uR = 1.0;
                double uI = 0.0;
                double angle = System.Math.PI / n;
                double wR = System.Math.Cos(angle);
                double wI = System.Math.Sin(angle);
                double t;
                Complex[] rotation = new Complex[n];

                for (i = 0; i < n; i++)
                {
                    rotation[i] = new Complex(uR, uI);
                    t = uR * wI + uI * wR;
                    uR = uR * wR - uI * wI;
                    uI = t;
                }

                complexRotation[numberOfBits - 1, directionIndex] = rotation;
            }
            return complexRotation[numberOfBits - 1, directionIndex];
        }
        /// <summary>
        /// Переупорядочивает данные для использования БПФ.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        private static void ReorderData(Complex[] data)
        {
            int length = data.Length;
            int[] rBits = GetReversedBits(Log2(length));
            Complex t;
            int i, s;

            for (i = 0; i < length; i++)
            {
                s = rBits[i];

                if (s > i)
                {
                    t = data[i];
                    data[i] = data[s];
                    data[s] = t;
                }
            }
        }
        /// <summary>
        /// Вычисляет степень двойки.
        /// </summary>
        /// <param name="power">Степень</param>
        /// <returns>Целое число</returns>
        private static int Pow2(int power)
        {
            return ((power >= 0) && (power <= 30)) ? (1 << power) : 0;
        }
        /// <summary>
        /// Выисляет логарифм по основанию 2.
        /// </summary>
        /// <param name="x">Целое число</param>
        /// <returns>Целое число</returns>
        private static int Log2(int x)
        {
            return (int)(Math.Log10(x) / 0.30102999566398);
        }
        #endregion
    }
    /// <summary>
    /// Определяет быстрое преобразование Хартли.
    /// <remarks>
    /// Алгоритм быстрого преобразования Хартли (БПХ) построен на основе быстрого преобразования Фурье (алгоритм Кули-Тьюки).
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Discrete_Hartley_transform
    /// </remarks>
    /// </summary>
    public class FastHartleyTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Преобразование Фурье.
        /// </summary>
        private FastFourierTransform FFT;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует быстрое преобразование Хартли.
        /// </summary>
        /// <param name="normalized">Нормализированное преобразование или нет</param>
        /// <param name="direction">Направление обработки</param>
        public FastHartleyTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(normalized, direction);
        }
        /// <summary>
        /// Нормализированное преобразование или нет.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.FFT.Normalized;
            }
            set
            {
                this.FFT.Normalized = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.FFT.Direction;
            }
            set
            {
                this.FFT.Direction = value;
            }
        }
        #endregion

        #region Fast Hartley Discrete Transform
        /// <summary>
        /// Прямое дискретное быстрое преобразование Хартли.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            Complex[] B = Matrice.ToComplex(A);
            B = FFT.Forward(B);

            int length = A.Length, i;
            double[] Hk = new double[length];

            for (i = 0; i < length; i++)
            {
                Hk[i] = B[i].Re + B[i].Im;
            }

            return Hk;
        }
        /// <summary>
        /// Обратное дискретное быстрое преобразование Хартли.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            Complex[] A = Matrice.ToComplex(B);
            A = FFT.Backward(A);

            int length = B.Length, i;
            double[] Hk = new double[length];

            for (i = 0; i < length; i++)
            {
                Hk[i] = A[i].Re - A[i].Im;
            }

            return Hk;
        }
        /// <summary>
        /// Прямое дискретное быстрое преобразование Хартли.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            Complex[,] B = Matrice.ToComplex(A);
            B = FFT.Forward(B);

            int width = A.GetLength(1), height = A.GetLength(0);
            double[,] Hk = new double[height, width];
            int i, j;

            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    Hk[i, j] = B[i, j].Re + B[i, j].Im;
                }
            }

            return Hk;
        }
        /// <summary>
        /// Обратное дискретное быстрое преобразование Хартли.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            Complex[,] A = Matrice.ToComplex(B);
            A = FFT.Backward(A);

            int width = B.GetLength(1), height = B.GetLength(0);
            double[,] Hk = new double[height, width];
            int i, j;

            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    Hk[i, j] = A[i, j].Re - A[i, j].Im;
                }
            }

            return Hk;
        }
        /// <summary>
        /// Прямое дискретное быстрое преобразование Хартли.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int length = A.Length, i;
            Complex[] B = FFT.Forward(A);
            Complex[] Hk = new Complex[length];

            for (i = 0; i < length; i++)
            {
                Hk[i] = B[i].Re + B[i].Im;
            }

            return Hk;
        }
        /// <summary>
        /// Обратное дискретное быстрое преобразование Хартли.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int length = B.Length, i;
            Complex[] A = FFT.Backward(B);
            Complex[] Hk = new Complex[length];

            for (i = 0; i < length; i++)
            {
                Hk[i] = A[i].Re - A[i].Im;
            }

            return Hk;
        }
        /// <summary>
        /// Прямое дискретное быстрое преобразование Хартли.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            Complex[,] B = FFT.Forward(A);
            int width = A.GetLength(1), height = A.GetLength(0);
            Complex[,] Hk = new Complex[height, width];
            int i, j;

            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    Hk[i, j] = B[i, j].Re + B[i, j].Im;
                }
            }

            return Hk;
        }
        /// <summary>
        /// Обратное дискретное быстрое преобразование Хартли.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = FFT.Backward(B);
            int width = B.GetLength(1), height = B.GetLength(0);
            Complex[,] Hk = new Complex[height, width];
            int i, j;

            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    Hk[i, j] = A[i, j].Re - A[i, j].Im;
                }
            }

            return Hk;
        }
        #endregion
    }
    /// <summary>
    /// Определяет быстрое преобразование Гильберта.
    /// <remarks>
    /// Данная вычислительная оптимизация реализована за счет быстрого преобразования Фурье (алгоритм Кули-Тьюки), сложность которого составляет O(Nlog(N)).
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Hilbert_transform
    /// </remarks>
    /// </summary>
    public class FastHilbertTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Преобразование Фурье.
        /// </summary>
        private FastFourierTransform FFT;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует быстрое преобразование Гильберта.
        /// </summary>
        /// <param name="normalized">Нормализированное преобразование или нет</param>
        /// <param name="direction">Направление обработки</param>
        public FastHilbertTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(normalized, Direction.Both);
            this.direction = direction;
        }
        /// <summary>
        /// Нормализированное преобразование или нет.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.FFT.Normalized;
            }
            set
            {
                this.FFT.Normalized = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        #endregion

        #region Fast Hilbert Discrete Transform
        /// <summary>
        /// Прямое дискретное быстрое преобразование Гильберта.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;

            // Перевод массива вещественных чисел в массив комплексных
            // и быстрое преобразование Фурье
            Complex[] F = FFT.Forward(A);

            // Перегруппировка:
            HilbertTransform.hilbertf(F, N);

            // Обратное преобразование Фурье:
            F = FFT.Backward(F);

            // Группировка:
            return HilbertTransform.hilbertb(A, F, N);
        }
        /// <summary>
        /// Обратное дискретное быстрое преобразование Гильберта.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            // Обратное преобразование Гильберта вычисляется путем
            // отображения прямого преобразования:
            // H^–1{h(t)} = –H{h(t)}

            int N = B.Length, i;
            Complex[] A = new Complex[N];

            for (i = 0; i < N; i++)
            {
                A[i] = new Complex(B[i].Re, 0);
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное быстрое преобразование Гильберта.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            Complex[,] B = (Complex[,])A.Clone();
            int N = A.GetLength(0);
            int M = A.GetLength(1);

            if (Direction == Direction.Both)
            {
                // 2-dimension horizontal Hilbert transform:
                FFT.Direction = Direction.Horizontal;
                B = FFT.Forward(B);

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    HilbertTransform.hilbertf(row, M);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                B = FFT.Backward(B);

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    Complex[] num = new Complex[M];

                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                        num[j] = A[i, j];
                    }

                    num = HilbertTransform.hilbertb(num, row, M);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = num[j];
                    }
                });


                // 2-dimension vertical Hilbert transform:
                FFT.Direction = Direction.Vertical;
                B = FFT.Forward(B);

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    HilbertTransform.hilbertf(col, N);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });

                B = FFT.Backward(B);

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    Complex[] num = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                        num[i] = A[i, j];
                    }

                    num = HilbertTransform.hilbertb(num, col, N);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = num[i];
                    }
                });

                FFT.Direction = Direction.Both;
            }
            else if (Direction == Direction.Vertical)
            {
                // 2-dimension vertical Hilbert transform:
                B = FFT.Forward(B);

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    HilbertTransform.hilbertf(col, N);

                    for (i = 0; i < M; i++)
                    {
                        B[i, j] = col[i];
                    }
                });

                B = FFT.Backward(B);

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    Complex[] num = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                        num[i] = A[i, j];
                    }

                    num = HilbertTransform.hilbertb(num, col, N);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = num[i];
                    }
                });
            }
            else
            {
                // 2-dimension horizontal Hilbert transform:
                B = FFT.Forward(B);

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    HilbertTransform.hilbertf(row, M);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                B = FFT.Backward(B);

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    Complex[] num = new Complex[M];

                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                        num[j] = A[i, j];
                    }

                    num = HilbertTransform.hilbertb(num, row, M);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = num[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное быстрое преобразование Гильберта.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = (Complex[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (Direction == Direction.Both)
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                }
                );
            }
            else if (Direction == Direction.Vertical)
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                }
                );
            }
            else
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
    /// <summary>
    /// Определяет быстрое дискретное косинусное преобразование.
    /// </summary>
    public class FastCosineTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Преобразование Фурье.
        /// </summary>
        private FastFourierTransform FFT;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует быстрое дискретное косинусное преобразование.
        /// </summary>
        /// <param name="direction">Направление обработки</param>
        public FastCosineTransform(Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(true, Direction.Both);
            this.direction = direction;
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        #endregion

        #region Fast Cosine Discrete Transform
        /// <summary>
        /// Прямое быстрое дискретное косинусное преобразование.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            int N = A.Length, N2 = N * 2, i, k;
            Complex[] B = new Complex[N2];

            for (i = 0; i < N; i++)
            {
                B[i] = A[i];
            }

            // Преобразование Фурье:
            B = FFT.Forward(B);

            double[] C = new double[N];
            Complex c = -Maths.I * Maths.Pi / N2;

            for (k = 0; k < N; k++)
            {
                C[k] = 2.0 * (B[k] * Maths.Exp(c * k)).Re;
            }
            C[0] = C[0] / Math.Sqrt(2); // DCT-I форма

            return C;
        }
        /// <summary>
        /// Обратное быстрое дискретное косинусное преобразование.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            int N = B.Length, N2 = N * 2, i, k;
            Complex[] A = new Complex[N2];
            double Bk, temp, c = Maths.Pi / N2;

            B[0] /= Math.Sqrt(2); // DCT-I форма

            for (k = 0; k < N; k++)
            {
                Bk = B[k];
                temp = k * c;
                A[k] = new Complex(Bk * Math.Cos(temp), Bk * Math.Sin(temp));
            }

            A = FFT.Backward(A);
            double[] C = new double[N];

            for (i = 0; i < N; i++)
            {
                C[i] = A[i].Re * 2;
            }

            return C;
        }
        /// <summary>
        /// Прямое быстрое дискретное косинусное преобразование.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            double[,] B = (double[,])A.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = Forward(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                }
                );

                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = Forward(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = Forward(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = Forward(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Обратное быстрое дискретное косинусное преобразование.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            double[,] A = (double[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                }
                );
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
    /// <summary>
    /// Определяет быстрое дискретное синусное преобразование.
    /// </summary>
    public class FastSineTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Преобразование Фурье.
        /// </summary>
        private FastFourierTransform FFT;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует быстрое дискретное синусное преобразование.
        /// </summary>
        /// <param name="direction">Направление обработки</param>
        public FastSineTransform(Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(true, Direction.Both);
            this.direction = direction;
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        #endregion

        #region Fast Sine Discrete Transform
        /// <summary>
        /// Прямое быстрое дискретное синусное преобразование.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            int N = A.Length, N2 = N * 2, i, k;
            Complex[] B = new Complex[N2];

            for (i = 0; i < N; i++)
            {
                B[i] = A[i];
            }

            // Преобразование Фурье:
            B = FFT.Forward(B);

            double[] C = new double[N];
            Complex c = -Maths.I * Maths.Pi / N;

            for (k = 0; k < N; k++)
            {
                C[k] = 2.0 * (B[k] * Maths.Exp(c * k)).Im;
            }

            // Перестановка:
            C[0] = A[N - 1] / Math.Sqrt(2);

            return C;
        }
        /// <summary>
        /// Обратное быстрое дискретное синусное преобразование.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            int N = B.Length, N2 = N * 2, i, k;
            Complex[] A = new Complex[N2];
            double Bk, temp, c = Maths.Pi / N;

            for (k = 0; k < N; k++)
            {
                Bk = B[k];
                temp = k * c;
                A[k] = new Complex(Bk * Math.Cos(temp), Bk * Math.Sin(temp));
            }

            // Преобразование Фурье:
            A = FFT.Backward(A);
            double[] C = new double[N];

            for (i = 0; i < N; i++)
            {
                C[i] = -A[i].Im * 2.0;
            }

            // Перестановка:
            C[N - 1] = B[0] * Math.Sqrt(2);

            return C;
        }
        /// <summary>
        /// Прямое быстрое дискретное синусное преобразование.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            double[,] B = (double[,])A.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = Forward(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                }
                );

                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = Forward(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = Forward(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = Forward(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Обратное быстрое дискретное синусное преобразование.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            double[,] A = (double[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                }
                );
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
    /// <summary>
    /// Определяет быстрое дискретное преобразование Лапласа.
    /// <remarks>
    /// Данный класс реализует дискретный эквивалент непрерывного преобразования Лапласа.
    /// Непрерывное преобразование Лапласа тесно связано с преобразованием Фурье. Иными словами, преобразование Фурье эквивалентно двустороннему преобразованию Лапласа с комплексным 
    /// аргументом: s = iω. Эта связь между преобразованиями часто используется для того, чтобы определить частотный спектр сигнала.
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Laplace_transform
    /// </remarks>
    /// </summary>
    public class FastLaplaceTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Преобразование Фурье.
        /// </summary>
        private FastFourierTransform FFT;
        /// <summary>
        /// Среднеквадратическое отклонение.
        /// </summary>
        private double sigma;
        /// <summary>
        /// Нормализованное преобразование или нет.
        /// </summary>
        private bool normalized;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует быстрое дискретное преобразование Лапласа.
        /// </summary>
        /// <param name="sigma">Среднеквадратическое отклонение (0, 1)</param>
        /// <param name="normalized">Нормализированное преобразование или нет</param>
        /// <param name="direction">Направление обработки</param>
        public FastLaplaceTransform(double sigma = 0.0005, bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(true, Direction.Vertical);
            Sigma = sigma; this.normalized = normalized; this.direction = direction;
        }
        /// <summary>
        /// Получает или задает значение среднеквадратического отклонения (0, 1).
        /// <remarks>
        /// В случае, если σ = 0, то преобразование Лапласа принимает вид преобразования Фурье.
        /// </remarks>
        /// </summary>
        public double Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Нормализированное преобразование или нет.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.normalized;
            }
            set
            {
                this.normalized = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        #endregion

        #region Laplace Discrete Transform
        /// <summary>
        /// Прямое дискретное преобразование Лапласа.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            // Fourier transform:
            Complex[] B = FFT.Forward(A);

            // Fourier to Laplace transform:
            laplace(B, sigma);
            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Лапласа.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            // Laplace to Fourier transform:
            Complex[] A = (Complex[])B.Clone();
            invlaplace(A, sigma);

            // Fourier transform:
            return FFT.Backward(A);
        }
        /// <summary>
        /// Прямое быстрое дискретное преобразование Лапласа.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            Complex[,] B = (Complex[,])A.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = Forward(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                }
                );

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = Forward(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = Forward(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = Forward(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Обратное быстрое дискретное преобразование Лапласа.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = (Complex[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                }
                );

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                }
                );
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;
                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Прямое преобразование Лапласа.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="sigma">Сигма</param>
        internal static void laplace(Complex[] v, double sigma = 0.003)
        {
            // forward laplacian transform
            // for discrete signal:
            int N = v.Length;

            for (int i = 0; i < N; i++)
            {
                v[i] = Math.Exp(-sigma * i) * v[i];
            }

            return;
        }
        /// <summary>
        /// Обратное преобразование Лапласа.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="sigma">Сигма</param>
        internal static void invlaplace(Complex[] v, double sigma = 0.003)
        {
            // inverse laplacian transform
            // for discrete signal:
            int N = v.Length;

            for (int i = 0; i < N; i++)
            {
                v[i] = v[i] / Math.Exp(-sigma * i);
            }

            return;
        }
        #endregion
    }
    #endregion

    #region Short-time Fourier analysis
    /// <summary>
    /// Определяет быстрое оконное преобразование Фурье.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Short-time_Fourier_transform
    /// </remarks>
    /// </summary>
    public class FastShortTimeFourierTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Преобразование Фурье.
        /// </summary>
        private FastFourierTransform FFT;
        /// <summary>
        /// Оконная функция.
        /// </summary>
        private IWindow window;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        /// <summary>
        /// Коэффициенты оконной функции.
        /// </summary>
        private double[] coefs;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует быстрое оконное преобразование Фурье.
        /// </summary>
        /// <param name="function">Оконная функция</param>
        /// <param name="normalized">Нормализированное преобразование или нет</param>
        /// <param name="direction">Направление обработки</param>
        public FastShortTimeFourierTransform(IWindow function, bool normalized = true, Direction direction = Direction.Vertical)
        {
            // fourier transform initialization:
            this.FFT = new FastFourierTransform(normalized, direction);
            Direction = direction;
            Function = function;

            // sampling window function:
            this.coefs = function.GetWindow().Add(1e-64);
        }
        /// <summary>
        /// Нормализированное преобразование или нет.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.FFT.Normalized;
            }
            set
            {
                this.FFT.Normalized = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        /// <summary>
        /// Получает или задает оконную функцию.
        /// </summary>
        public IWindow Function
        {
            get { return window; }
            set { window = value; }
        }
        #endregion

        #region Short-Time Fourier Discrete Transform
        /// <summary>
        /// Прямое дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;

            // Быстрое преобразование Фурье:
            Complex[] B = FFT.Forward(A);

            // Прямое оконное преобразование:
            ShortTimeFourierTransform.shortfourierf(B, N, coefs);

            return B;
        }
        /// <summary>
        /// Обратное дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;

            // Быстрое обратное преобразование Фурье:
            Complex[] A = (Complex[])B.Clone();

            // Обратное оконное преобразование:
            ShortTimeFourierTransform.shortfourierb(A, N, coefs);

            A = FFT.Backward(A);
            return A;
        }
        /// <summary>
        /// Прямое дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            // Fourier transform:
            Complex[,] B = FFT.Forward(A);
            int N = A.GetLength(0);
            int M = A.GetLength(1);

            if (direction == Direction.Both)
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    ShortTimeFourierTransform.shortfourierf(row, M, coefs);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    ShortTimeFourierTransform.shortfourierf(col, N, coefs);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });

            }
            else if (direction == Direction.Vertical)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    ShortTimeFourierTransform.shortfourierf(col, N, coefs);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    ShortTimeFourierTransform.shortfourierf(row, M, coefs);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = (Complex[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    ShortTimeFourierTransform.shortfourierb(col, N, coefs);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });

                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    ShortTimeFourierTransform.shortfourierb(row, M, coefs);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    ShortTimeFourierTransform.shortfourierb(col, N, coefs);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });
            }
            else
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    ShortTimeFourierTransform.shortfourierb(row, M, coefs);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }

            A = FFT.Backward(A);

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконное преобразование Фурье.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Short-time_Fourier_transform
    /// </remarks>
    /// </summary>
    public class ShortTimeFourierTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Преобразование Фурье.
        /// </summary>
        private FourierTransform FFT;
        /// <summary>
        /// Оконная функция.
        /// </summary>
        private IWindow window;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        /// <summary>
        /// Коэффициенты оконной функции.
        /// </summary>
        private double[] coefs;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует оконное преобразование Фурье.
        /// </summary>
        /// <param name="function">Оконная функция</param>
        /// <param name="normalized">Нормализированное преобразование или нет</param>
        /// <param name="direction">Направление обработки</param>
        public ShortTimeFourierTransform(IWindow function, bool normalized = true, Direction direction = Direction.Vertical)
        {
            // fourier transform initialization:
            this.FFT = new FourierTransform(normalized, direction);
            Direction = direction;
            Function = function;

            // sampling window function:
            this.coefs = function.GetWindow().Add(1e-64);
        }
        /// <summary>
        /// Нормализированное преобразование или нет.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.FFT.Normalized;
            }
            set
            {
                this.FFT.Normalized = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        /// <summary>
        /// Получает или задает оконную функцию.
        /// </summary>
        public IWindow Function
        {
            get { return window; }
            set { window = value; }
        }
        #endregion

        #region Short-Time Fourier Discrete Transform
        /// <summary>
        /// Прямое дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;

            // Быстрое преобразование Фурье:
            Complex[] B = FFT.Forward(A);

            // Прямое оконное преобразование:
            ShortTimeFourierTransform.shortfourierf(B, N, coefs);

            return B;
        }
        /// <summary>
        /// Обратное дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;

            // Быстрое обратное преобразование Фурье:
            Complex[] A = (Complex[])B.Clone();

            // Обратное оконное преобразование:
            ShortTimeFourierTransform.shortfourierb(A, N, coefs);

            A = FFT.Backward(A);
            return A;
        }
        /// <summary>
        /// Прямое дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            // Fourier transform:
            Complex[,] B = FFT.Forward(A);
            int N = A.GetLength(0);
            int M = A.GetLength(1);

            if (direction == Direction.Both)
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    ShortTimeFourierTransform.shortfourierf(row, M, coefs);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    ShortTimeFourierTransform.shortfourierf(col, N, coefs);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });

            }
            else if (direction == Direction.Vertical)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    ShortTimeFourierTransform.shortfourierf(col, N, coefs);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    ShortTimeFourierTransform.shortfourierf(row, M, coefs);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = (Complex[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    ShortTimeFourierTransform.shortfourierb(col, N, coefs);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });

                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    ShortTimeFourierTransform.shortfourierb(row, M, coefs);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    ShortTimeFourierTransform.shortfourierb(col, N, coefs);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });
            }
            else
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    ShortTimeFourierTransform.shortfourierb(row, M, coefs);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }

            A = FFT.Backward(A);

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Прямое оконное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <param name="N">Размерность</param>
        /// <param name="coefs">Оконная функция</param>
        internal static void shortfourierf(Complex[] B, int N, double[] coefs)
        {
            int time = coefs.Length; for (int i = 0; i < N; i++) { B[i] *= coefs[Maths.Mod(i, time)]; } return;
        }
        /// <summary>
        /// Обратное оконное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <param name="N">Размерность</param>
        /// <param name="coefs">Оконная функция</param>
        internal static void shortfourierb(Complex[] B, int N, double[] coefs)
        {
            int time = coefs.Length; for (int i = 0; i < N; i++) { B[i] /= coefs[Maths.Mod(i, time)]; } return;
        }
        #endregion
    }
    #endregion

    #region Gabor analysis
    /// <summary>
    /// Определяет группу ортогональных базисов и дискретных преобразований Вейля-Гейзенберга.
    /// <remarks>
    /// Базисы Вейля-Гейзенберга используются для получения частотно-временных характеристик сигнала.
    /// Более подробную информацию можно найти на сайте:
    /// https://elibrary.ru/item.asp?id=29767333
    /// </remarks>
    /// </summary>
    public class WeylHeisenbergTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Количество сдвигов по частоте.
        /// </summary>
        private int m;
        /// <summary>
        /// Среднеквадратическое отклонение гауссиана.
        /// </summary>
        private double sigma;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует группу ортогональных базисов и преобразований Вейля-Гейзенберга.
        /// </summary>
        /// <param name="m">Количество сдвигов по частоте [4, N].<remarks>
        /// </remarks>Четное число</param>
        /// <param name="sigma">Среднеквадратическое отклонение (0, +inf)</param>
        /// <param name="direction">Направление обработки</param>
        public WeylHeisenbergTransform(int m = 8, double sigma = 0.0025, Direction direction = Direction.Vertical)
        {
            Sigma = sigma; M = m; Direction = direction;
        }
        /// <summary>
        /// Получает или задает количество сдвигов по частоте [4, N].
        /// <remarks>
        /// Четное число.
        /// </remarks>
        /// </summary>
        public int M
        {
            get
            {
                return this.m;
            }
            set
            {
                if (value <= 2 || !Maths.IsEven(value))
                    throw new Exception("Неверное значение аргумента");

                this.m = value;
            }
        }
        /// <summary>
        /// Получает или задает значение среднеквадратического отклонения (0, +inf).
        /// </summary>
        public double Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        #endregion

        #region Weyl-Heisenberg static components
        /// <summary>
        /// Возвращает комплексную матрицу базиса Габора.
        /// <remarks>
        /// Размерность матрицы [N, N], где N = M * L.
        /// Данный базис не является ортогональным.
        /// </remarks>
        /// </summary>
        /// <param name="M">Количество сдвигов по частоте</param>
        /// <param name="L">Количество сдвигов по времени</param>
        /// <param name="sigma">Среднеквадратическое отклонение</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Gabor(int M, int L, double sigma)
        {
            int N = M * L, N2 = N * 2;                                  // Определение параметров,
            double[] g0 = Gabor(N, sigma);                           // вычисление формирующей WH-функции,
            Complex[,] G = new Complex[N, N2];                          // комплексный прямоугольный сигнальный базис,
            Complex c = 2 * Maths.Pi * Maths.I;                         // комплексный коэффициент преобразования,
            double a = M / 2.0 - 1.0;                                   // вычисление оптимального фазового параметра.

            Parallel.For(0, N, n =>                                     // Использование параллельных циклов.
            {
                double phase = n - a / 2.0;                             // Фазовый сдвиг.
                int k, l, u, i, j;
                Complex exp, psi;

                for (k = 0; k < M; k++)
                {
                    exp = Maths.Exp(c * k / M * phase);                 // Экспоненциальный коэффициент.

                    for (l = 0; l < L; l++)
                    {
                        u = l * M + k;                                  // Элемент матрицы,
                        i = Maths.Mod(n - l * M, N);                    // сдвиг по времени,
                        j = Maths.Mod(n + M / 2 - l * M, N);            // сдвиг по частоте,

                        psi = new Complex(
                            (g0[i] * exp).Re,                           // Комлпексный сигнальный базис,
                            (Maths.I * g0[j] * exp).Re);                // < Ψ Re, Ψ Im >

                        G[n, u] = psi;
                    }
                }
            }
            );

            return G;
        }
        /// <summary>
        /// Возвращает комплексную матрицу оптимального базиса Вейля-Гейзенберга.
        /// <remarks>
        /// Размерность матрицы [N, N], где N = M * L.
        /// </remarks>
        /// </summary>
        /// <param name="M">Количество сдвигов по частоте</param>
        /// <param name="L">Количество сдвигов по времени</param>
        /// <param name="sigma">Среднеквадратическое отклонение</param>
        /// <returns>Матрица</returns>
        public static Complex[,] WeylHeisenberg(int M, int L, double sigma)
        {
            int N = M * L;                                              // Определение параметров,
            double[] g0 = Matrice.Zak(Gabor(N, sigma), M);           // вычисление формирующей WH-функции,
            Complex[,] G = new Complex[N, N];                           // комплексный прямоугольный сигнальный базис,
            Complex c = 2 * Maths.Pi * Maths.I;                         // комплексный коэффициент преобразования,
            double a = M / 2.0;                                         // вычисление оптимального фазового параметра.

            Parallel.For(0, N, n =>                                     // Использование параллельных циклов.
            {
                double phase = n - a / 2.0;                             // Фазовый сдвиг.
                int k, l, u, i, j;
                Complex exp, psi;

                for (k = 0; k < M; k++)
                {
                    exp = Maths.Exp(c * k / M * phase);                 // Экспоненциальный коэффициент.

                    for (l = 0; l < L; l++)
                    {
                        u = l * M + k;                                  // Элемент матрицы,
                        i = Maths.Mod(n - l * M, N);                    // сдвиг по времени,
                        j = Maths.Mod(n + M / 2 - l * M, N);            // сдвиг по частоте,

                        psi = new Complex(
                            (g0[i] * exp).Re,                           // Комлпексный сигнальный базис,
                            (Maths.I * g0[j] * exp).Re);                // < Ψ Re, Ψ Im >

                        G[n, u] = psi;
                    }
                }
            }
            );

            return G;
        }
        /// <summary>
        /// Возвращает комплексную матрицу базиса Вейля-Гейзенберга.
        /// <remarks>
        /// Размерность матрицы [N, N], где N = M * L.
        /// </remarks>
        /// </summary>
        /// <param name="g0">Формирующая функция</param>
        /// <param name="M">Количество сдвигов по частоте</param>
        /// <returns>Матрица</returns>
        public static Complex[,] WeylHeisenberg(double[] g0, int M)
        {
            int N = g0.Length, L = N / M;                               // Определение параметров,
            Complex[,] G = new Complex[N, N];                           // комплексный прямоугольный сигнальный базис,
            Complex c = 2 * Maths.Pi * Maths.I;                         // комплексный коэффициент преобразования,
            double a = M / 2.0;                                         // вычисление оптимального фазового параметра.

            Parallel.For(0, N, n =>                                     // Использование параллельных циклов.
            {
                double phase = n - a / 2.0;                             // Фазовый сдвиг.
                int k, l, u, i, j;
                Complex exp, psi;

                for (k = 0; k < M; k++)
                {
                    exp = Maths.Exp(c * k / M * phase);                 // Экспоненциальный коэффициент.

                    for (l = 0; l < L; l++)
                    {
                        u = l * M + k;                                  // Элемент матрицы,
                        i = Maths.Mod(n - l * M, N);                    // сдвиг по времени,
                        j = Maths.Mod(n + M / 2 - l * M, N);            // сдвиг по частоте,

                        psi = new Complex(
                            (g0[i] * exp).Re,                           // Комлпексный сигнальный базис,
                            (Maths.I * g0[j] * exp).Re);                // < Ψ Re, Ψ Im >

                        G[n, u] = psi;
                    }
                }
            }
            );

            return G;
        }
        /// <summary>
        /// Возвращает вектор значений усеченного гауссиана.
        /// </summary>
        /// <param name="N">Количество отсчетов функции</param>
        /// <param name="sigma">Среднеквадратическое отклонение</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Gabor(int N, double sigma)
        {
            // Рассматривается только случай сопряженной N-симметрии.
            int N2 = N / 2, i; double xmin = -N2, xmax = N2 + 1;
            double[] g0 = new double[N];
            double t = xmin;

            for (i = 0; i < N; i++, t++)
            {
                g0[i] = Special.Gabor(t, sigma);
            }
            return g0;
        }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Прямое дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            int L = N / m;

            if (L <= 0)
                throw new Exception("Количество сдвигов по частоте определено неверно");

            Complex[,] U = WeylHeisenbergTransform.WeylHeisenberg(m, L, sigma);
            Complex[] B = Matrice.Dot(A, U.Hermitian());

            // Выделение действительной части от преобразования:
            int length = B.Length;
            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;
            int L = N / m;

            if (L <= 0)
                throw new Exception("Количество сдвигов по частоте определено неверно");

            Complex[,] U = WeylHeisenbergTransform.WeylHeisenberg(m, L, sigma);
            Complex[] A = Matrice.Dot(B, U);

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            int l0 = N / m, l1 = M / m;

            if (l0 <= 0 || l1 <= 0)
                throw new Exception("Количество сдвигов по частоте определено неверно");

            Complex[,] U = WeylHeisenbergTransform.WeylHeisenberg(m, l0, sigma);
            Complex[,] V = WeylHeisenbergTransform.WeylHeisenberg(m, l1, sigma);
            Complex[,] B;

            if (direction == Direction.Both)
            {
                B = U.Hermitian().Dot(A).Dot(V);
            }
            else if (direction == Direction.Vertical)
            {
                B = U.Hermitian().Dot(A);
            }
            else
            {
                B = A.Dot(V);
            }
            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            int l0 = N / m, l1 = M / m;

            if (l0 <= 0 || l1 <= 0)
                throw new Exception("Количество сдвигов по частоте определено неверно");

            Complex[,] U = WeylHeisenbergTransform.WeylHeisenberg(m, l0, sigma);
            Complex[,] V = WeylHeisenbergTransform.WeylHeisenberg(m, l1, sigma);
            Complex[,] A;

            if (direction == Direction.Both)
            {
                A = U.Dot(B).Dot(V.Hermitian());
            }
            else if (direction == Direction.Vertical)
            {
                A = U.Dot(B);
            }
            else
            {
                A = B.Dot(V.Hermitian());
            }
            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс быстрого преобразования Вейля-Гейзенберга.
    /// <remarks>
    /// Класс представляет вычислительно эффективную реализацию одномерного и двумерного дикретных ортогональных
    /// преобразований Вейля-Гейзенберга.
    /// Более подробную информацию можно найти на сайте:
    /// https://elibrary.ru/title_about.asp?id=58245
    /// </remarks>
    /// </summary>
    public class FastWeylHeisenbergTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Преобразование Фурье.
        /// <remarks>
        /// Используется как вспомогательный компонент.
        /// </remarks>
        /// </summary>
        private static FastFourierTransform FFT = new FastFourierTransform(false, Direction.Horizontal);
        /// <summary>
        /// Количество сдвигов по частоте.
        /// </summary>
        private int m;
        /// <summary>
        /// Среднеквадратическое отклонение гауссиана.
        /// </summary>
        private double sigma;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует быстрое преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="m">Количество сдвигов по частоте [4, N].<remarks>
        /// </remarks>Четное число</param>
        /// <param name="sigma">Среднеквадратическое отклонение (0, +inf)</param>
        /// <param name="direction">Направление обработки</param>
        public FastWeylHeisenbergTransform(int m = 8, double sigma = 0.0025, Direction direction = Direction.Vertical)
        {
            Sigma = sigma; M = m; Direction = direction;
        }
        /// <summary>
        /// Получает или задает количество сдвигов по частоте [4, N].
        /// <remarks>
        /// Четное число.
        /// </remarks>
        /// </summary>
        public int M
        {
            get
            {
                return this.m;
            }
            set
            {
                if (value <= 2 || !Maths.IsEven(value))
                    throw new Exception("Неверное значение аргумента");

                this.m = value;
            }
        }
        /// <summary>
        /// Получает или задает значение среднеквадратического отклонения (0, +inf).
        /// </summary>
        public double Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Прямое дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            double[] g0 = Matrice.Zak(WeylHeisenbergTransform.Gabor(A.Length, sigma), m);
            return FastWeylHeisenbergTransform.WHT(A, g0, m);
        }
        /// <summary>
        /// Обратное дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            double[] g0 = Matrice.Zak(WeylHeisenbergTransform.Gabor(B.Length, sigma), m);
            return FastWeylHeisenbergTransform.IWHT(B, g0, m);
        }
        /// <summary>
        /// Прямое дискретное быстрое преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            Complex[,] B = (Complex[,])A.Clone();
            int N = B.GetLength(0), M = B.GetLength(1);

            double[] g0 = Matrice.Zak(WeylHeisenbergTransform.Gabor(N, sigma), m);
            double[] g1 = Matrice.Zak(WeylHeisenbergTransform.Gabor(M, sigma), m);

            if (direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = WHT(row, g1, m);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = WHT(col, g0, m);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = WHT(col, g0, m);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = WHT(row, g1, m);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное быстрое преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = (Complex[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            double[] g0 = Matrice.Zak(WeylHeisenbergTransform.Gabor(N, sigma), m);
            double[] g1 = Matrice.Zak(WeylHeisenbergTransform.Gabor(M, sigma), m);

            if (direction == Direction.Both)
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    col = IWHT(col, g0, m);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    row = IWHT(row, g1, m);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    col = IWHT(col, g0, m);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    row = IWHT(row, g1, m);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Public static components
        /// <summary>
        /// Прямое быстрое дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="input">Одномерный массив</param>
        /// <param name="g0">Формирующая WH-функция</param>
        /// <param name="M">Количество сдвигов по частоте</param>
        /// <returns>Одноменый массив</returns>
        public static Complex[] WHT(Complex[] input, double[] g0, int M)
        {
            // Функция реализует быстрый алгоритм прямого преобразования Вейля-Гейзенберга, 
            // изложенный в следующих статьях:
            // A. Vahlin, "EFFICIENT ALGORITHMS FOR MODULATION AND DEMODULATION IN OFDM-SYSTEMS" [1].
            // В.М. Асирян, В.П. Волчков, "ЭФФЕКТИВНАЯ РЕАЛИЗАЦИЯ ПРЯМОГО ПРЕОБРАЗОВАНИЯ ВЕЙЛЯ-ГЕЙЗЕНБЕРГА" [2].
            // Алгоритм является вычислительно эффективным при больших значениях M.

            // Параметры преобразования и выходной сигнал:
            int N = input.Length, L = N / M, M2 = M / 2, M4 = M2 / 2;
            Complex[] output = new Complex[N];
            Complex[] exp = FastWeylHeisenbergTransform.GetRotation(M);

            // Вспомогательные матрицы и переменные:
            Complex[,] s0 = new Complex[M, L];
            Complex[,] a0 = new Complex[M, L];
            Complex[,] b0 = new Complex[M, L];
            Complex[,] A0 = new Complex[L, M];
            Complex[,] B0 = new Complex[L, M];
            Complex[,] A1 = new Complex[L, M2];
            Complex[,] B1 = new Complex[L, M2];
            Complex c1re, c2re;
            Complex c1im, c2im;
            int k, i, j, u, n, m, l;

            // 1. Формирование матриц перестановок:
            for (m = 0; m < M; m++)
            {
                for (n = 0; n < L; n++)
                {
                    u = n * M;
                    i = Maths.Mod(m + M4 + u, N);
                    j = Maths.Mod(m - M4 + u, N);
                    k = Maths.Mod(-m - M4 - u, N);

                    s0[m, n] = input[k];
                    a0[m, n] = g0[i];
                    b0[m, n] = g0[j];
                }
            }

            // 2. Вычисление циклической свертки для матриц
            // полученных в пункте 1:
            for (l = 0; l < L; l++)
            {
                for (n = 0; n < L; n++)
                {
                    k = Maths.Mod(n - l, L);

                    for (m = 0; m < M; m++)
                    {
                        A0[l, m] += a0[m, n] * s0[m, k];
                        B0[l, m] += b0[m, n] * s0[m, k];
                    }
                }
            }

            // Переменные замены:
            Complex x, y, z, w;

            // 3. Вычисление новых последовательностей:
            for (l = 0; l < L; l++)
            {
                for (m = 0; m < M2; m++)
                {
                    // Замена переменных для действительной части:
                    x = A0[l, m];
                    y = A0[l, m + M2];
                    z = A0[l, M2 - m].Conjugate;
                    w = A0[l, Maths.Mod(M - m, M)].Conjugate;

                    // Вычисление выражений:
                    c1re = x + y + z + w;
                    c2re = x - y - z + w;

                    // Аналогичная замена переменных для мнимой части:
                    x = B0[l, m];
                    y = B0[l, m + M2];
                    z = B0[l, M2 - m].Conjugate;
                    w = B0[l, Maths.Mod(M - m, M)].Conjugate;

                    // Вычисление выражений:
                    c1im = x + y - z - w;
                    c2im = x - y + z - w;

                    // Вычисление элементов матриц:
                    A1[l, m] = 1.0 / (2) * (c1re + Maths.I * c2re * exp[m]);
                    B1[l, m] = 1.0 / (2 * Maths.I) * (c1im + Maths.I * c2im * exp[m]);
                }
            }

            // 4. Быстрое обратное M/2- точечное преобразование Фурье по строкам:
            A1 = FFT.Backward(Matrice.Conjugate(A1));
            B1 = FFT.Backward(Matrice.Conjugate(B1));

            // 5. Формирование комплексного вектора сигнала:
            for (k = 0; k < M2; k++)
            {
                for (l = 0; l < L; l++)
                {
                    // Четные и нечетные индексы фильтра:
                    i = l * M + 2 * k;
                    j = l * M + 2 * k + 1;

                    // Замена переменных:
                    x = A1[l, k];
                    y = B1[l, k];

                    // Вычисление элементов вектора:
                    output[i] = x.Re + Maths.I * y.Re;
                    output[j] = x.Im + Maths.I * y.Im;
                }
            }

            // Результат прямого WH-преобразования.
            return output;
        }
        /// <summary>
        /// Обратное быстрое дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="input">Одномерный массив</param>
        /// <param name="g0">Формирующая WH-функция</param>
        /// <param name="M">Количество сдвигов по частоте</param>
        /// <returns>Одноменый массив</returns>
        public static Complex[] IWHT(Complex[] input, double[] g0, int M)
        {
            // Функция реализует быстрый алгоритм обратного преобразования Вейля-Гейзенберга, 
            // изложенный в следующих статьях:
            // A. Vahlin, "EFFICIENT ALGORITHMS FOR MODULATION AND DEMODULATION IN OFDM-SYSTEMS" [1].
            // В.М. Асирян, В.П. Волчков, "ЭФФЕКТИВНАЯ РЕАЛИЗАЦИЯ ПРЯМОГО ПРЕОБРАЗОВАНИЯ ВЕЙЛЯ-ГЕЙЗЕНБЕРГА" [2].
            // Алгоритм является вычислительно эффективным при больших значениях M.

            // Вспомогательные матрицы и переменные:
            int N = input.Length, L = N / M, M2 = M / 2, M4 = M2 / 2;
            Complex[] output = new Complex[N];
            Complex[,] A1 = new Complex[L, M];
            Complex[,] B1 = new Complex[L, M];
            Complex[] exp = FastWeylHeisenbergTransform.GetRotation(M);
            Complex s;
            int n, k, l;

            // 1. Формирование вещественных матриц сдвигов:
            for (k = 0; k < M; k++)
            {
                for (l = 0; l < L; l++)
                {
                    // Замена переменных:
                    s = input[k + l * M];

                    // Вычисление элементов матриц:
                    A1[l, k] = s.Re;
                    B1[l, k] = s.Im;
                }
            }

            // 2. Вычисление матриц:
            Complex[,] Za = new Complex[L, M2];
            Complex[,] Zb = new Complex[L, M2];

            for (k = 0; k < M2; k++)
            {
                for (l = 0; l < L; l++)
                {
                    Za[l, k] = A1[l, k * 2] + Maths.I * A1[l, k * 2 + 1];
                    Zb[l, k] = B1[l, k * 2] + Maths.I * B1[l, k * 2 + 1];
                }
            }

            // 3. Быстрое обратное M/2- точечное преобразование Фурье.
            Za = Matrice.Conjugate(FFT.Backward(Za));
            Zb = Matrice.Conjugate(FFT.Backward(Zb));

            // Переменные замены:
            Complex a0, a1, b0, b1;
            Complex x, y, u, v;

            // 4. Формирование новых матриц:
            for (k = 0; k < M2; k++)
            {
                for (l = 0; l < L; l++)
                {
                    // Для матрицы A(l, k)
                    // Замена переменных:
                    a0 = Za[l, k]; a1 = Za[l, Maths.Mod(M - k, M2)].Conjugate;

                    // Вычисление выражений:
                    x = 1.0 / (2) * (a0 + a1);
                    y = 1.0 / (2 * Maths.I) * (a0 - a1);
                    y *= exp[k];

                    // Формирование матрицы:
                    A1[l, k] = x + y;
                    A1[l, k + M2] = x - y;


                    // Для матрицы B(l, k)
                    // Замена переменных:
                    b0 = Zb[l, k]; b1 = Zb[l, Maths.Mod(M - k, M2)].Conjugate;

                    // Вычисление выражений:
                    u = 1.0 / (2) * (b0 + b1);
                    v = 1.0 / (2 * Maths.I) * (b0 - b1);
                    v *= exp[k];

                    // Формирование матрицы:
                    B1[l, k] = u + v;
                    B1[l, k + M2] = u - v;
                }
            }

            // 5. Формирование выходного сигнала:
            for (l = 0; l < L; l++)
            {
                for (n = 0; n < N; n++)
                {
                    // Вычисление выражения:
                    output[n] += A1[l, Maths.Mod(n - M4, M)] * g0[Maths.Mod(n - l * M, N)] - Maths.I
                               * B1[l, Maths.Mod(n - M4, M)] * g0[Maths.Mod(n - l * M + M2, N)];
                }
            }

            // Результат обратного WH-преобразования.
            return output;
        }
        /// <summary>
        /// Возвращает одномерный массив фазовых поворотов.
        /// </summary>
        /// <param name="M">Количество сдвигов по частоте</param>
        /// <returns>Одномерный массив</returns>
        private static Complex[] GetRotation(int M)
        {
            int M2 = M / 2;
            Complex[] phase = new Complex[M2];
            for (int k = 0; k < M2; k++)
            {
                phase[k] = Maths.Exp(Maths.I * 2 * Math.PI / M * k);
            }
            return phase;
        }
        #endregion
    }
    #endregion

    #region Pyramid transforms
    /// <summary>
    /// Определяет класс представления сигнала в виде пирамиды Лапласа.
    /// Более подробную информацию можно найти на сайте:
    /// http://www.cs.toronto.edu/~jepson/csc320/notes/pyramids.pdf
    /// </summary>
    public class LaplacianPyramidTransform : IPyramidTransform
    {
        #region Private data
        int levels;
        #endregion

        #region Pyramid components
        /// <summary>
        /// Инициализирует класс представления сигнала в виде пирамиды Лапласа.
        /// </summary>
        /// <param name="levels">Количество уровней представления</param>
        public LaplacianPyramidTransform(int levels)
        {
            this.Levels = levels;
        }
        /// <summary>
        /// Инициализирует класс представления сигнала в виде пирамиды Лапласа.
        /// </summary>
        public LaplacianPyramidTransform()
        {
            this.Levels = int.MaxValue;
        }
        /// <summary>
        /// Получает или задает количество уровней представления.
        /// </summary>
        public int Levels
        {
            get
            {
                return this.levels;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.levels = value;
            }

        }
        #endregion

        #region Apply voids
        // **************************************************
        //            Laplacian Pyramid Transform
        // **************************************************
        // ORIGINALS: Burt, P., and Adelson, E. H.
        // IEEE Transactions on Communication, COM-31:532-540 
        // (1983).
        // Designed by Asiryan Valeriy (c), 2015-2019
        // Moscow, Russia.
        // **************************************************

        /// <summary>
        /// Прямое пирамидоидальное преобразование Лапласа.
        /// </summary>
        /// <param name="data">Двумерный массив</param>
        /// <returns>Пирамида</returns>
        public double[][,] Forward(double[,] data)
        {
            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = (int)Math.Min((Math.Log(Math.Min(r, c)) / Math.Log(2)), levels);
            double[][,] lapl = new double[nlev][,];
            double[,] I, J = data;

            for (int i = 0; i < nlev - 1; i++)
            {
                I = GaussianPyramidTransform.downsample(J);
                lapl[i] = GaussianPyramidTransform.sub(J, GaussianPyramidTransform.upsample(I));
                J = I;
            }

            lapl[nlev - 1] = J;
            return lapl;
        }
        /// <summary>
        /// Прямое пирамидоидальное преобразование Лапласа.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        /// <returns>Пирамида</returns>
        public double[][] Forward(double[] data)
        {
            int r = data.Length;
            int nlev = (int)Math.Min((Math.Log(r) / Math.Log(2)), levels);

            double[][] lapl = new double[nlev][];
            double[] I, J = data;

            for (int i = 0; i < nlev - 1; i++)
            {
                I = GaussianPyramidTransform.downsample(J);
                lapl[i] = GaussianPyramidTransform.sub(J, GaussianPyramidTransform.upsample(I));
                J = I;
            }

            lapl[nlev - 1] = J;
            return lapl;
        }
        /// <summary>
        /// Обратное пирамидоидальное преобразование Лапласа.
        /// </summary>
        /// <param name="pyramid">Пирамида</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[][,] pyramid)
        {
            int nlev = pyramid.Length - 1;
            double[,] I = pyramid[nlev];

            for (int i = nlev - 1; i >= 0; i--)
            {
                I = GaussianPyramidTransform.add(pyramid[i], GaussianPyramidTransform.upsample(I));
            }

            return I;
        }
        /// <summary>
        /// Обратное пирамидоидальное преобразование Лапласа.
        /// </summary>
        /// <param name="pyramid">Пирамида</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[][] pyramid)
        {
            int nlev = pyramid.Length;
            double[] I = pyramid[nlev];

            for (int i = nlev - 1; i >= 0; i--)
            {
                I = GaussianPyramidTransform.add(pyramid[i], GaussianPyramidTransform.upsample(I));
            }

            return I;
        }
        /// <summary>
        /// Прямое пирамидоидальное преобразование Лапласа.
        /// </summary>
        /// <param name="data">Двумерный массив</param>
        /// <returns>Пирамида</returns>
        public Complex[][,] Forward(Complex[,] data)
        {
            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = (int)Math.Min((Math.Log(Math.Min(r, c)) / Math.Log(2)), levels);
            Complex[][,] lapl = new Complex[nlev][,];
            Complex[,] I, J = data;

            for (int i = 0; i < nlev - 1; i++)
            {
                I = GaussianPyramidTransform.downsample(J);
                lapl[i] = GaussianPyramidTransform.sub(J, GaussianPyramidTransform.upsample(I));
                J = I;
            }

            lapl[nlev - 1] = J;
            return lapl;
        }
        /// <summary>
        /// Прямое пирамидоидальное преобразование Лапласа.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        /// <returns>Пирамида</returns>
        public Complex[][] Forward(Complex[] data)
        {
            int r = data.Length;
            int nlev = (int)Math.Min((Math.Log(r) / Math.Log(2)), levels);

            Complex[][] lapl = new Complex[nlev][];
            Complex[] I, J = data;

            for (int i = 0; i < nlev - 1; i++)
            {
                I = GaussianPyramidTransform.downsample(J);
                lapl[i] = GaussianPyramidTransform.sub(J, GaussianPyramidTransform.upsample(I));
                J = I;
            }

            lapl[nlev - 1] = J;
            return lapl;
        }
        /// <summary>
        /// Обратное пирамидоидальное преобразование Лапласа.
        /// </summary>
        /// <param name="pyramid">Пирамида</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[][,] pyramid)
        {
            int nlev = pyramid.Length - 1;
            Complex[,] I = pyramid[nlev];

            for (int i = nlev - 1; i >= 0; i--)
            {
                I = GaussianPyramidTransform.add(pyramid[i], GaussianPyramidTransform.upsample(I));
            }

            return I;
        }
        /// <summary>
        /// Обратное пирамидоидальное преобразование Лапласа.
        /// </summary>
        /// <param name="pyramid">Пирамида</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[][] pyramid)
        {
            int nlev = pyramid.Length;
            Complex[] I = pyramid[nlev];

            for (int i = nlev - 1; i >= 0; i--)
            {
                I = GaussianPyramidTransform.add(pyramid[i], GaussianPyramidTransform.upsample(I));
            }

            return I;
        }
        #endregion

        #region Gaussian pyramid to Laplacian pyramid
        /// <summary>
        /// Прямое пирамидоидальное преобразование Лапласа.
        /// </summary>
        /// <param name="data">Пирамида Гусса</param>
        /// <returns>Пирамида</returns>
        public double[][,] Forward(double[][,] data)
        {
            int nlev = data.Length;
            double[][,] lapl = new double[nlev][,];

            for (int i = 1; i < nlev; i++)
            {
                lapl[i - 1] = GaussianPyramidTransform.sub(data[i - 1], GaussianPyramidTransform.upsample(data[i]));
            }

            lapl[nlev - 1] = data[nlev - 1];
            return lapl;
        }
        /// <summary>
        /// Прямое пирамидоидальное преобразование Лапласа.
        /// </summary>
        /// <param name="data">Пирамида Гусса</param>
        /// <returns>Пирамида</returns>
        public double[][] Forward(double[][] data)
        {
            int nlev = data.Length;
            double[][] lapl = new double[nlev][];

            for (int i = 1; i < nlev; i++)
            {
                lapl[i - 1] = GaussianPyramidTransform.sub(data[i - 1], GaussianPyramidTransform.upsample(data[i]));
            }

            lapl[nlev - 1] = data[nlev - 1];
            return lapl;
        }
        /// <summary>
        /// Прямое пирамидоидальное преобразование Лапласа.
        /// </summary>
        /// <param name="data">Пирамида Гусса</param>
        /// <returns>Пирамида</returns>
        public Complex[][,] Forward(Complex[][,] data)
        {
            int nlev = data.Length;
            Complex[][,] lapl = new Complex[nlev][,];

            for (int i = 1; i < nlev; i++)
            {
                lapl[i - 1] = GaussianPyramidTransform.sub(data[i - 1], GaussianPyramidTransform.upsample(data[i]));
            }

            lapl[nlev - 1] = data[nlev - 1];
            return lapl;
        }
        /// <summary>
        /// Прямое пирамидоидальное преобразование Лапласа.
        /// </summary>
        /// <param name="data">Пирамида Гусса</param>
        /// <returns>Пирамида</returns>
        public Complex[][] Forward(Complex[][] data)
        {
            int nlev = data.Length;
            Complex[][] lapl = new Complex[nlev][];

            for (int i = 1; i < nlev; i++)
            {
                lapl[i - 1] = GaussianPyramidTransform.sub(data[i - 1], GaussianPyramidTransform.upsample(data[i]));
            }

            lapl[nlev - 1] = data[nlev - 1];
            return lapl;
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс представления сигнала в виде пирамиды Гаусса.
    /// Более подробную информацию можно найти на сайте:
    /// http://www.cs.toronto.edu/~jepson/csc320/notes/pyramids.pdf
    /// </summary>
    public class GaussianPyramidTransform : IPyramidTransform
    {
        #region Private data
        int levels;
        #endregion

        #region Pyramid components
        /// <summary>
        /// Инициализирует класс представления сигнала в виде пирамиды Гаусса.
        /// </summary>
        public GaussianPyramidTransform()
        {
            this.Levels = int.MaxValue;
        }
        /// <summary>
        /// Инициализирует класс представления сигнала в виде пирамиды Гаусса.
        /// </summary>
        /// <param name="levels">Количество уровней представления</param>
        public GaussianPyramidTransform(int levels)
        {
            this.Levels = levels;
        }
        /// <summary>
        /// Получает или задает количество уровней представления.
        /// </summary>
        public int Levels
        {
            get
            {
                return this.levels;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.levels = value;
            }

        }
        #endregion

        #region Apply voids
        // **************************************************
        //            Gaussian Pyramid Transform
        // **************************************************
        // ORIGINALS: Burt, P., and Adelson, E. H.
        // IEEE Transactions on Communication, COM-31:532-540 
        // (1983).
        // Designed by Asiryan Valeriy (c), 2015-2019
        // Moscow, Russia.
        // **************************************************

        /// <summary>
        /// Прямое пирамидоидальное преобразование Гаусса.
        /// </summary>
        /// <param name="data">Двумерный массив</param>
        /// <returns>Пирамида</returns>
        public double[][,] Forward(double[,] data)
        {
            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = (int)Math.Min((Math.Log(Math.Min(r, c))
                / Math.Log(2)), levels);

            double[][,] pyr = new double[nlev][,];
            double[,] dummy = (double[,])data.Clone();

            for (int i = 0; i < nlev; i++)
            {
                pyr[i] = dummy;
                dummy = downsample(dummy);
            }

            return pyr;
        }
        /// <summary>
        /// Прямое пирамидоидальное преобразование Гаусса.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        /// <returns>Пирамида</returns>
        public double[][] Forward(double[] data)
        {
            int r = data.Length;
            int nlev = (int)Math.Min((Math.Log(r) / Math.Log(2)), levels);

            double[][] pyr = new double[nlev][];
            double[] dummy = (double[])data.Clone();

            for (int i = 0; i < nlev; i++)
            {
                pyr[i] = dummy;
                dummy = downsample(dummy);
            }

            return pyr;
        }
        /// <summary>
        /// Обратное пирамидоидальное преобразование Гаусса (не поддерживается).
        /// </summary>
        /// <param name="pyramid">Пирамида</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[][,] pyramid)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное пирамидоидальное преобразование Гаусса (не поддерживается).
        /// </summary>
        /// <param name="pyramid">Пирамида</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[][] pyramid)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое пирамидоидальное преобразование Гаусса.
        /// </summary>
        /// <param name="data">Двумерный массив</param>
        /// <returns>Пирамида</returns>
        public Complex[][,] Forward(Complex[,] data)
        {
            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = (int)Math.Min((Math.Log(Math.Min(r, c))
                / Math.Log(2)), levels);

            Complex[][,] pyr = new Complex[nlev][,];
            Complex[,] dummy = (Complex[,])data.Clone();

            for (int i = 0; i < nlev; i++)
            {
                pyr[i] = dummy;
                dummy = downsample(dummy);
            }

            return pyr;
        }
        /// <summary>
        /// Прямое пирамидоидальное преобразование Гаусса.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        /// <returns>Пирамида</returns>
        public Complex[][] Forward(Complex[] data)
        {
            int r = data.Length;
            int nlev = (int)Math.Min((Math.Log(r) / Math.Log(2)), levels);

            Complex[][] pyr = new Complex[nlev][];
            Complex[] dummy = (Complex[])data.Clone();

            for (int i = 0; i < nlev; i++)
            {
                pyr[i] = dummy;
                dummy = downsample(dummy);
            }

            return pyr;
        }
        /// <summary>
        /// Обратное пирамидоидальное преобразование Гаусса (не поддерживается).
        /// </summary>
        /// <param name="pyramid">Пирамида</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[][,] pyramid)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное пирамидоидальное преобразование Гаусса (не поддерживается).
        /// </summary>
        /// <param name="pyramid">Пирамида</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[][] pyramid)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Static methods
        /// <summary>
        /// Возвращает одномерный фильтр Гаусса.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public static double[] Filter
        {
            get
            {
                return new double[] { .0625, .25, .375, .25, .0625 };
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Увеличивает размерность в два раза. 
        /// </summary>
        /// <param name="u">Матрица</param>
        /// <returns>Матрица</returns>
        internal static double[,] upsample(double[,] u)
        {
            int r = u.GetLength(0), c = u.GetLength(1);
            int n = r * 2, m = c * 2;
            int i, j, k, l;
            double[,] v = new double[n, m];

            for (k = 0, i = 0; i < r; i++, k += 2)
            {
                for (l = 0, j = 0; j < c; j++, l += 2)
                {
                    v[k + 1, l] = u[i, j];
                    v[k, l + 1] = u[i, j];
                    v[k, l] = u[i, j];
                    v[k + 1, l + 1] = u[i, j];
                }
            }

            MeanFilter.boxf(v, Direction.Both, 4);

            return v;
        }
        /// <summary>
        /// Увеличивает размерность в два раза. 
        /// </summary>
        /// <param name="u">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        internal static double[] upsample(double[] u)
        {
            int r = u.GetLength(0);
            int n = r * 2;
            int i, k;
            double[] v = new double[n];

            for (k = 0, i = 0; i < r; i++, k += 2)
            {
                v[k] = u[i];
                v[k + 1] = u[i];
            }

            MeanFilter.boxf(v, r, 4);

            return v;
        }
        /// <summary>
        /// Уменьшает размерность в два раза. 
        /// </summary>
        /// <param name="u">Матрица</param>
        /// <returns>Матрица</returns>
        internal static double[,] downsample(double[,] u)
        {
            int r = u.GetLength(0);
            int c = u.GetLength(1);
            int n = (r + 1) / 2, m = (c + 1) / 2;
            int i, j, k, l;
            double[,] v = new double[n, m];

            for (k = 0, i = 0; i < r; i += 2, k++)
            {
                for (l = 0, j = 0; j < c; j += 2, l++)
                {
                    v[k, l] = u[i, j];
                }
            }

            MeanFilter.boxf(v, Direction.Both, 4);

            return v;
        }
        /// <summary>
        /// Уменьшает размерность в два раза. 
        /// </summary>
        /// <param name="u">Матрица</param>
        /// <returns>Матрица</returns>
        internal static double[] downsample(double[] u)
        {
            int r = u.Length;
            int n = (r + 1) / 2;
            int i, k;
            double[] v = new double[n];

            for (k = 0, i = 0; i < r; i += 2, k++)
            {
                v[k] = u[i];
            }

            MeanFilter.boxf(v, r, 4);

            return v;
        }
        /// <summary>
        /// Складывает две матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        internal static double[,] add(double[,] m, double[,] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            int mr = (int)Math.Min(m.GetLength(1), n.GetLength(1));
            double[,] H = new double[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] + n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Складывает два вектора.
        /// </summary>
        /// <param name="m">Одномерный массив</param>
        /// <param name="n">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        internal static double[] add(double[] m, double[] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            double[] v = new double[ml];
            int i;

            for (i = 0; i < ml; i++)
            {
                v[i] = m[i] + n[i];
            }
            return v;
        }
        /// <summary>
        /// Вычитает две матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        internal static double[,] sub(double[,] m, double[,] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            int mr = (int)Math.Min(m.GetLength(1), n.GetLength(1));
            double[,] H = new double[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] - n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Вычитает два вектора.
        /// </summary>
        /// <param name="m">Одномерный массив</param>
        /// <param name="n">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        internal static double[] sub(double[] m, double[] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            double[] v = new double[ml];
            int i;

            for (i = 0; i < ml; i++)
            {
                v[i] = m[i] - n[i];
            }
            return v;
        }
        /// <summary>
        /// Увеличивает размерность в два раза. 
        /// </summary>
        /// <param name="u">Матрица</param>
        /// <returns>Матрица</returns>
        internal static Complex[,] upsample(Complex[,] u)
        {
            int r = u.GetLength(0), c = u.GetLength(1);
            int n = r * 2, m = c * 2;
            int i, j, k, l;
            Complex[,] v = new Complex[n, m];

            for (k = 0, i = 0; i < r; i++, k += 2)
            {
                for (l = 0, j = 0; j < c; j++, l += 2)
                {
                    v[k + 1, l] = u[i, j];
                    v[k, l + 1] = u[i, j];
                    v[k, l] = u[i, j];
                    v[k + 1, l + 1] = u[i, j];
                }
            }

            MeanFilter.boxf(v, Direction.Both, 4);

            return v;
        }
        /// <summary>
        /// Увеличивает размерность в два раза. 
        /// </summary>
        /// <param name="u">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        internal static Complex[] upsample(Complex[] u)
        {
            int r = u.GetLength(0);
            int n = r * 2;
            int i, k;
            Complex[] v = new Complex[n];

            for (k = 0, i = 0; i < r; i++, k += 2)
            {
                v[k] = u[i];
                v[k + 1] = u[i];
            }

            MeanFilter.boxf(v, r, 4);

            return v;
        }
        /// <summary>
        /// Уменьшает размерность в два раза. 
        /// </summary>
        /// <param name="u">Матрица</param>
        /// <returns>Матрица</returns>
        internal static Complex[,] downsample(Complex[,] u)
        {
            int r = u.GetLength(0);
            int c = u.GetLength(1);
            int n = (r + 1) / 2, m = (c + 1) / 2;
            int i, j, k, l;
            Complex[,] v = new Complex[n, m];

            for (k = 0, i = 0; i < r; i += 2, k++)
            {
                for (l = 0, j = 0; j < c; j += 2, l++)
                {
                    v[k, l] = u[i, j];
                }
            }

            MeanFilter.boxf(v, Direction.Both, 4);

            return v;
        }
        /// <summary>
        /// Уменьшает размерность в два раза. 
        /// </summary>
        /// <param name="u">Матрица</param>
        /// <returns>Матрица</returns>
        internal static Complex[] downsample(Complex[] u)
        {
            int r = u.Length;
            int n = (r + 1) / 2;
            int i, k;
            Complex[] v = new Complex[n];

            for (k = 0, i = 0; i < r; i += 2, k++)
            {
                v[k] = u[i];
            }

            MeanFilter.boxf(v, r, 4);

            return v;
        }
        /// <summary>
        /// Складывает две матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        internal static Complex[,] add(Complex[,] m, Complex[,] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            int mr = (int)Math.Min(m.GetLength(1), n.GetLength(1));
            Complex[,] H = new Complex[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] + n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Складывает два вектора.
        /// </summary>
        /// <param name="m">Одномерный массив</param>
        /// <param name="n">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        internal static Complex[] add(Complex[] m, Complex[] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            Complex[] v = new Complex[ml];
            int i;

            for (i = 0; i < ml; i++)
            {
                v[i] = m[i] + n[i];
            }
            return v;
        }
        /// <summary>
        /// Вычитает две матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        internal static Complex[,] sub(Complex[,] m, Complex[,] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            int mr = (int)Math.Min(m.GetLength(1), n.GetLength(1));
            Complex[,] H = new Complex[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] - n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Вычитает два вектора.
        /// </summary>
        /// <param name="m">Одномерный массив</param>
        /// <param name="n">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        internal static Complex[] sub(Complex[] m, Complex[] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            Complex[] v = new Complex[ml];
            int i;

            for (i = 0; i < ml; i++)
            {
                v[i] = m[i] - n[i];
            }
            return v;
        }
        #endregion
    }
    #endregion
    
    #region Processing filters
    /// <summary>
    /// Определяет фильтр частот.
    /// </summary>
    public class FrequencyFilter : IFilter
    {
        #region Private data
        /// <summary>
        /// Минимальные и максимальные частоты.
        /// </summary>
        protected RangeInt frequencyRange = new RangeInt(0, 64);
        #endregion

        #region Filter components
        /// <summary>
        /// Получает или задает диапазон частот.
        /// </summary>
        public RangeInt FrequencyRange
        {
            get
            {
                return frequencyRange;
            }
            set
            {
                frequencyRange = value;
            }
        }
        /// <summary>
        /// Инициализирует фильтр частот.
        /// </summary>
        public FrequencyFilter() { }
        /// <summary>
        /// Инициализирует фильтр частот.
        /// </summary>
        /// <param name="frequencyRange">Диапазон частот</param>
        public FrequencyFilter(RangeInt frequencyRange)
        {
            this.frequencyRange = frequencyRange;
        }
        /// <summary>
        /// Инициализирует фильтр частот.
        /// </summary>
        /// <param name="min">Минимальная частота</param>
        /// <param name="max">Максимальная частота</param>
        public FrequencyFilter(int min, int max)
        {
            this.frequencyRange = new RangeInt(min, max);
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Матрица</param>
        public void Apply(double[,] data)
        {
            int height = data.GetLength(0);
            int width = data.GetLength(1);

            // Половины измерений:
            int hh = height >> 1;
            int hw = width >> 1;

            // Минимальные и максимальные частоты:
            int min = frequencyRange.Min;
            int max = frequencyRange.Max;

            int i, j, x, y, d;

            for (i = 0; i < height; i++)
            {
                y = i - hh;

                for (j = 0; j < width; j++)
                {
                    x = j - hw;
                    d = (int)Math.Sqrt(x * x + y * y);

                    // Значения за пределами фильтра.
                    if ((d > max) || (d < min))
                    {
                        data[i, j] = 0;
                    }
                }
            }
            return;
        }
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        public void Apply(double[] data)
        {
            int length = data.Length;
            int hh = length >> 1;

            // Минимальные и максимальные частоты:
            int min = frequencyRange.Min;
            int max = frequencyRange.Max;

            int i, d;

            for (i = 0; i < length; i++)
            {
                d = i - hh;

                // Значения за пределами фильтра.
                if ((d > max) || (d < min))
                {
                    data[i] = 0;
                }
            }

            return;
        }
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Матрица</param>
        public void Apply(Complex[,] data)
        {
            int height = data.GetLength(0);
            int width = data.GetLength(1);

            // Половины измерений:
            int hh = height >> 1;
            int hw = width >> 1;

            // Минимальные и максимальные частоты:
            int min = frequencyRange.Min;
            int max = frequencyRange.Max;

            int i, j, x, y, d;

            for (i = 0; i < height; i++)
            {
                y = i - hh;

                for (j = 0; j < width; j++)
                {
                    x = j - hw;
                    d = (int)Math.Sqrt(x * x + y * y);

                    // Значения за пределами фильтра.
                    if ((d > max) || (d < min))
                    {
                        data[i, j].Re = 0;
                        data[i, j].Im = 0;
                    }
                }
            }
            return;
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        public void Apply(Complex[] data)
        {
            int length = data.Length;
            int hh = length >> 1;

            // Минимальные и максимальные частоты:
            int min = frequencyRange.Min;
            int max = frequencyRange.Max;

            int i, d;

            for (i = 0; i < length; i++)
            {
                d = i - hh;

                // Значения за пределами фильтра.
                if ((d > max) || (d < min))
                {
                    data[i].Re = 0;
                    data[i].Im = 0;
                }
            }

            return;
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр сжатия по пороговому значению.
    /// <remarks>
    /// Суть фильтра заключается в обнулении всех тех значений сигнала, модуль которых лежит ниже определенного значения - порога.
    /// </remarks>
    /// </summary>
    public class CompressFilter : IFilter
    {
        #region Private data
        /// <summary>
        /// Пороговое значение.
        /// </summary>
        private double threshold = 0;
        /// <summary>
        /// Тип сжатия.
        /// </summary>
        private Compress compresstype = Compress.Abs;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр сжатия по пороговому значению.
        /// </summary>
        public CompressFilter() { }
        /// <summary>
        /// Инициализирует фильтр сжатия по пороговому значению.
        /// </summary>
        /// <param name="threshold">Пороговое значение</param>
        /// <param name="compresstype">Тип сжатия</param>
        public CompressFilter(double threshold, Compress compresstype = Compress.Abs)
        {
            this.threshold = threshold;
            this.compresstype = compresstype;
        }
        /// <summary>
        /// Получает или задает тип сжатия.
        /// </summary>
        public Compress CompressType
        {
            get
            {
                return this.compresstype;
            }
            set
            {
                this.compresstype = value;
            }
        }
        /// <summary>
        /// Получает или задает пороговое значение.
        /// </summary>
        public double Threshold
        {
            get
            {
                return this.threshold;
            }
            set
            {
                this.threshold = value;
            }
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Одомерный массив</param>
        public void Apply(double[] data)
        {
            int length = data.Length;
            int i;

            if (this.compresstype == Compress.Abs)
            {
                for (i = 0; i < length; i++)
                {
                    if (Math.Abs(data[i]) < threshold)
                    {
                        data[i] = 0;
                    }
                }
            }
            else if (this.compresstype == Compress.Over)
            {
                for (i = 0; i < length; i++)
                {
                    if (data[i] > threshold)
                    {
                        data[i] = 0;
                    }
                }
            }
            else
            {
                for (i = 0; i < length; i++)
                {
                    if (data[i] < threshold)
                    {
                        data[i] = 0;
                    }
                }
            }
            return;
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Одомерный массив</param>
        public void Apply(Complex[] data)
        {
            int length = data.Length;
            int i;

            if (this.compresstype == Compress.Abs)
            {
                for (i = 0; i < length; i++)
                {
                    if (Maths.Abs(data[i]) < threshold)
                    {
                        data[i] = 0;
                    }
                }
            }
            else if (this.compresstype == Compress.Over)
            {
                for (i = 0; i < length; i++)
                {
                    if (data[i].Re > threshold)
                    {
                        data[i].Re = 0;
                    }
                    if (data[i].Im > threshold)
                    {
                        data[i].Im = 0;
                    }
                }
            }
            else
            {
                for (i = 0; i < length; i++)
                {
                    if (data[i].Re < threshold)
                    {
                        data[i].Re = 0;
                    }
                    if (data[i].Im < threshold)
                    {
                        data[i].Im = 0;
                    }
                }
            }

            return;
        }
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Матрица</param>
        public void Apply(double[,] data)
        {
            int width = data.GetLength(1);
            int height = data.GetLength(0);
            int i, j;

            if (this.compresstype == Compress.Abs)
            {
                for (i = 0; i < height; i++)
                {
                    for (j = 0; j < width; j++)
                    {
                        if (Math.Abs(data[i, j]) < threshold)
                        {
                            data[i, j] = 0;
                        }
                    }
                }
            }
            else if (this.compresstype == Compress.Over)
            {
                for (i = 0; i < height; i++)
                {
                    for (j = 0; j < width; j++)
                    {
                        if (data[i, j] > threshold)
                        {
                            data[i, j] = 0;
                        }
                    }
                }
            }
            else
            {
                for (i = 0; i < height; i++)
                {
                    for (j = 0; j < width; j++)
                    {
                        if (data[i, j] < threshold)
                        {
                            data[i, j] = 0;
                        }
                    }
                }
            }
            return;
        }
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Матрица</param>
        public void Apply(Complex[,] data)
        {
            int width = data.GetLength(1);
            int height = data.GetLength(0);
            int i, j;

            if (this.compresstype == Compress.Abs)
            {
                for (i = 0; i < height; i++)
                {
                    for (j = 0; j < width; j++)
                    {
                        if (Maths.Abs(data[i, j]) < threshold)
                        {
                            data[i, j] = 0;
                        }
                    }
                }
            }
            else if (this.compresstype == Compress.Over)
            {
                for (i = 0; i < height; i++)
                {
                    for (j = 0; j < width; j++)
                    {
                        if (data[i, j].Re > threshold)
                        {
                            data[i, j] = 0;
                        }
                        if (data[i, j].Im > threshold)
                        {
                            data[i, j] = 0;
                        }
                    }
                }
            }
            else
            {
                for (i = 0; i < height; i++)
                {
                    for (j = 0; j < width; j++)
                    {
                        if (data[i, j].Re < threshold)
                        {
                            data[i, j] = 0;
                        }
                        if (data[i, j].Im < threshold)
                        {
                            data[i, j] = 0;
                        }
                    }
                }
            }
            return;
        }
        #endregion

        #region Compress type
        /// <summary>
        /// Определяет тип сжатия.
        /// </summary>
        public enum Compress
        {
            #region Types
            /// <summary>
            /// Сжатие по модулю.
            /// </summary>
            Abs,
            /// <summary>
            /// Сжатие значений меньше порогового.
            /// </summary>
            Under,
            /// <summary>
            /// Сжатие значений больше порогового.
            /// </summary>
            Over,
            #endregion
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр локального усреднения.
    /// <remarks>
    /// Фильтр вычисляет средние значения в локальных областях сигнала.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Box_blur
    /// </remarks>
    /// </summary>
    public class MeanFilter : IFilter, IBlendFilter
    {
        #region Private data
        /// <summary>
        /// Радиус фильтра.
        /// </summary>
        private int r;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        /// <summary>
        /// Множитель.
        /// </summary>
        private double factor;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр локального усреднения.
        /// </summary>
        /// <param name="radius">Радиус фильтра (>1)</param>
        /// <param name="direction">Направление обработки</param>
        /// <param name="factor">Множитель [-1, 1]</param>
        public MeanFilter(int radius, Direction direction = Direction.Both, double factor = -1.0)
        {
            this.Radius = radius;
            this.direction = direction;
            this.factor = factor;
        }
        /// <summary>
        /// Получает или задает значение радиуса фильтра.
        /// </summary>
        public int Radius
        {
            get
            {
                return this.r;
            }
            set
            {
                if (value < 1)
                    throw new Exception("Неверное значение аргумента");

                this.r = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        /// <summary>
        /// Получает или задает значение множителя [-1, 1].
        /// </summary>
        public double Factor
        {
            get
            {
                return this.factor;
            }
            set
            {
                this.factor = value;
            }
        }
        #endregion

        #region Public apply voids
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Матрица</param>
        public void Apply(double[,] data)
        {
            if (this.factor != 0)
            {
                int N = data.GetLength(0);
                int M = data.GetLength(1);
                int i, j;

                double[,] copy = (double[,])data.Clone();
                MeanFilter.boxf(copy, this.direction, this.r);

                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < M; j++)
                    {
                        data[i, j] = (1.0 + this.factor) * (data[i, j] - copy[i, j]) + copy[i, j];
                    }
                }

                return;
            }

            return;
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        public void Apply(double[] data)
        {
            int l = data.Length;

            if (this.factor != 0)
            {
                double[] copy = (double[])data.Clone();
                MeanFilter.boxf(copy, l, r);

                for (int i = 0; i < l; i++)
                {
                    data[i] = (1.0 + this.factor) * (data[i] - copy[i]) + copy[i];
                }

                return;
            }

            return;
        }
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Матрица</param>
        public void Apply(Complex[,] data)
        {
            if (this.factor != 0)
            {
                int N = data.GetLength(0);
                int M = data.GetLength(1);
                int i, j;

                Complex[,] copy = (Complex[,])data.Clone();
                MeanFilter.boxf(copy, this.direction, this.r);

                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < M; j++)
                    {
                        data[i, j] = (1.0 + this.factor) * (data[i, j] - copy[i, j]) + copy[i, j];
                    }
                }

                return;
            }

            return;
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        public void Apply(Complex[] data)
        {
            int l = data.Length;

            if (this.factor != 0)
            {
                double[] copy = (double[])data.Clone();
                MeanFilter.boxf(copy, l, r);

                for (int i = 0; i < l; i++)
                {
                    data[i] = (1.0 + this.factor) * (data[i] - copy[i]) + copy[i];
                }

                return;
            }

            return;
        }
        #endregion

        #region Blender apply voids
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Набор матриц</param>
        /// <returns>Матрица</returns>
        public double[,] Apply(double[][,] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0), c = data[0].GetLength(1);
            double[,] sum = new double[r, c], cur, low;
            int i, j, k;

            // process
            for (i = 0; i < length; i++)
            {
                // lowpass filter:
                cur = data[i];
                low = (double[,])cur.Clone();
                MeanFilter.boxf(low, this.direction, this.r);

                for (j = 0; j < r; j++)
                {
                    for (k = 0; k < c; k++)
                    {
                        // summarize high and low frequencies:
                        sum[j, k] += (1.0 + this.factor) / length * (cur[j, k] - low[j, k]) + (low[j, k] / length);
                    }
                }
            }

            return sum;
        }
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Набор матриц</param>
        /// <returns>Матрица</returns>
        public Complex[,] Apply(Complex[][,] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0), c = data[0].GetLength(1);
            Complex[,] sum = new Complex[r, c], cur, low;
            int i, j, k;

            // process
            for (i = 0; i < length; i++)
            {
                // lowpass filter:
                cur = data[i];
                low = (Complex[,])cur.Clone();
                MeanFilter.boxf(low, this.direction, this.r);

                for (j = 0; j < r; j++)
                {
                    for (k = 0; k < c; k++)
                    {
                        // summarize high and low frequencies:
                        sum[j, k] += (1.0 + this.factor) / length * (cur[j, k] - low[j, k]) + (low[j, k] / length);

                    }
                }
            }

            return sum;
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Набор векторов</param>
        /// <returns>Одномерный массив</returns>
        public double[] Apply(double[][] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0);
            double[] sum = new double[r], cur, low;
            int i, j;

            // process
            for (i = 0; i < length; i++)
            {
                // lowpass filter:
                cur = data[i];
                low = (double[])cur.Clone();
                MeanFilter.boxf(low, length, this.r);

                for (j = 0; j < r; j++)
                {
                    // summarize high and low frequencies:
                    sum[j] += (1.0 + this.factor) / length * (cur[j] - low[j]) + (low[j] / length);
                }
            }

            return sum;
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Набор векторов</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Apply(Complex[][] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0);
            Complex[] sum = new Complex[r], cur, low;
            int i, j;

            // process
            for (i = 0; i < length; i++)
            {
                // lowpass filter:
                cur = data[i];
                low = (Complex[])cur.Clone();
                MeanFilter.boxf(low, length, this.r);

                for (j = 0; j < r; j++)
                {
                    // summarize high and low frequencies:
                    sum[j] += (1.0 + this.factor) / length * (cur[j] - low[j]) + (low[j] / length);
                }
            }

            return sum;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Матрица</param>
        /// <param name="direction">Направление обработки</param>
        /// <param name="r">Радиус обработки</param>
        internal static void boxf(double[,] data, Direction direction, int r)
        {
            int N = data.GetLength(0);
            int M = data.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = data[i, j];
                    }

                    MeanFilter.boxf(row, M, r);

                    for (j = 0; j < M; j++)
                    {
                        data[i, j] = row[j];
                    }
                }
                );

                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = data[i, j];
                    }

                    MeanFilter.boxf(col, N, r);

                    for (i = 0; i < N; i++)
                    {
                        data[i, j] = col[i];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = data[i, j];
                    }

                    MeanFilter.boxf(col, N, r);

                    for (i = 0; i < N; i++)
                    {
                        data[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = data[i, j];
                    }

                    MeanFilter.boxf(row, M, r);

                    for (j = 0; j < M; j++)
                    {
                        data[i, j] = row[j];
                    }
                });
            }

            return;
        }
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Матрица</param>
        /// <param name="direction">Направление обработки</param>
        /// <param name="r">Радиус обработки</param>
        internal static void boxf(Complex[,] data, Direction direction, int r)
        {
            int N = data.GetLength(0);
            int M = data.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = data[i, j];
                    }

                    MeanFilter.boxf(row, M, r);

                    for (j = 0; j < M; j++)
                    {
                        data[i, j] = row[j];
                    }
                }
                );

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = data[i, j];
                    }

                    MeanFilter.boxf(col, N, r);

                    for (i = 0; i < N; i++)
                    {
                        data[i, j] = col[i];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = data[i, j];
                    }

                    MeanFilter.boxf(col, N, r);

                    for (i = 0; i < N; i++)
                    {
                        data[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = data[i, j];
                    }

                    MeanFilter.boxf(row, M, r);

                    for (j = 0; j < M; j++)
                    {
                        data[i, j] = row[j];
                    }
                });
            }

            return;
        }
        /// <summary>
        /// Реализует одномерный фильтр арифметического локального усреднения.
        /// </summary>
        /// <param name="input">Одномерный массив</param>
        /// <param name="l">Длина сигнала</param>
        /// <param name="r">Ядро фильтра</param>
        internal static void boxf(double[] input, int l, int r)
        {
            // Исключение по размерности:
            if (l == 1) return;
            int h = r >= l ? l - 1 : r;

            // Определение параметров фильтра:
            int v = h >> 1;
            int dl = l - v;
            double s = 0;
            int i;

            // Вычисление глобальной суммы [0, h):
            for (i = 0; i < h; i++)
            {
                s += input[i];
            }
            // Вычисление фильтра на отрезке [0, v):
            for (i = 0; i < v; i++)
            {
                input[i] = s / h;
            }
            // Вычисление фильтра на отрезке [v, l-v):
            for (i = v; i < dl; i++)
            {
                s = s - input[i - v] + input[i + v];
                input[i] = s / h;
            }
            // Вычисление фильтра на отрезке [l-v, l):
            for (i = dl; i < l; i++)
            {
                s = s - input[i - v] + input[i];
                input[i] = s / h;
            }

            return;
        }
        /// <summary>
        /// Реализует одномерный фильтр арифметического локального усреднения.
        /// </summary>
        /// <param name="input">Одномерный массив</param>
        /// <param name="l">Длина сигнала</param>
        /// <param name="r">Ядро фильтра</param>
        internal static void boxf(Complex[] input, int l, int r)
        {
            // Исключение по размерности:
            if (l == 1) return;
            int h = r >= l ? l - 1 : r;

            // Определение параметров фильтра:
            int v = h >> 1;
            int dl = l - v;
            Complex s = 0;
            int i;

            // Вычисление глобальной суммы [0, h):
            for (i = 0; i < h; i++)
            {
                s += input[i];
            }
            // Вычисление фильтра на отрезке [0, v):
            for (i = 0; i < v; i++)
            {
                input[i] = s / h;
            }
            // Вычисление фильтра на отрезке [v, l-v):
            for (i = v; i < dl; i++)
            {
                s = s - input[i - v] + input[i + v];
                input[i] = s / h;
            }
            // Вычисление фильтра на отрезке [l-v, l):
            for (i = dl; i < l; i++)
            {
                s = s - input[i - v] + input[i];
                input[i] = s / h;
            }

            return;
        }
        /// <summary>
        /// Реализует двумерный фильтр математического усреднения.
        /// </summary>
        /// <param name="data">Набор матриц</param>
        /// <returns>Матрица</returns>
        internal static double[,] meanf(double[][,] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0), c = data[0].GetLength(1);
            double[,] sum = new double[r, c];
            double[,] wei = new double[r, c];
            int i, j, k;

            // process
            for (i = 0; i < length; i++)
            {
                for (j = 0; j < r; j++)
                {
                    for (k = 0; k < c; k++)
                    {
                        // summarize all signals:
                        sum[j, k] += data[i][j, k] / length;
                    }
                }
            }

            return sum;
        }
        /// <summary>
        /// Реализует двумерный фильтр математического усреднения.
        /// </summary>
        /// <param name="data">Набор матриц</param>
        /// <returns>Матрица</returns>
        internal static Complex[,] meanf(Complex[][,] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0), c = data[0].GetLength(1);
            Complex[,] sum = new Complex[r, c];
            int i, j, k;

            // process
            for (i = 0; i < length; i++)
            {
                for (j = 0; j < r; j++)
                {
                    for (k = 0; k < c; k++)
                    {
                        // summarize all signals:
                        sum[j, k] += data[i][j, k] / length;
                    }
                }
            }

            return sum;
        }
        /// <summary>
        /// Реализует одномерный фильтр математического усреднения.
        /// </summary>
        /// <param name="data">Набор векторов</param>
        /// <returns>Одномерный массив</returns>
        internal static double[] meanf(double[][] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0);
            double[] sum = new double[r];
            int i, j;

            // process
            for (i = 0; i < length; i++)
            {
                for (j = 0; j < r; j++)
                {
                    // summarize all signals:
                    sum[j] += data[i][j] / length;
                }
            }

            return sum;
        }
        /// <summary>
        /// Реализует одномерный фильтр математического усреднения.
        /// </summary>
        /// <param name="data">Набор векторов</param>
        /// <returns>Одномерный массив</returns>
        internal static Complex[] meanf(Complex[][] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0);
            Complex[] sum = new Complex[r];
            int i, j;

            // process
            for (i = 0; i < length; i++)
            {
                for (j = 0; j < r; j++)
                {
                    // summarize all signals:
                    sum[j] += data[i][j] / length;
                }
            }

            return sum;
        }
        #endregion
    }
    /// <summary>
    /// Определяет направленный фильтр.
    /// <remarks>
    /// Данный фильтр представляет собой вычислительно эффективный аналог билатерального фильтра (bilateral filter).
    /// Более подробную информацию можно найти на сайте:
    /// http://kaiminghe.com/eccv10/index.html
    /// </remarks>
    /// </summary>
    public class GuidedFilter : IFilter, IBlendFilter
    {
        #region Private data
        /// <summary>
        /// Погрешность.
        /// </summary>
        private double eps;
        /// <summary>
        /// Множитель.
        /// </summary>
        private double factor;
        /// <summary>
        /// Радиус фильтра.
        /// </summary>
        private int radius;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует направленный фильтр.
        /// </summary>
        /// <param name="radius">Радиус фильтра (>1)</param>
        /// <param name="eps">Погрешность (0, 1)</param>
        /// <param name="factor">Множитель [-1, 1]</param>
        public GuidedFilter(int radius, double eps = 0.025, double factor = -1.0)
        {
            this.radius = radius;
            this.Eps = eps;
            this.factor = factor;
        }
        /// <summary>
        /// Получает или задает значение радиуса фильтра.
        /// </summary>
        public int Radius
        {
            get
            {
                return this.radius;
            }
            set
            {
                this.radius = value;
            }
        }
        /// <summary>
        /// Получает или задает значение погрешности (0, 1).
        /// <remarks>
        /// Оптимальное значение ε = 0.025.
        /// </remarks>
        /// </summary>
        public double Eps
        {
            get
            {
                return this.eps;
            }
            set
            {
                this.eps = Maths.Double(value);
            }
        }
        /// <summary>
        /// Получает или задает значение множителя [-1, 1].
        /// </summary>
        public double Factor
        {
            get
            {
                return this.factor;
            }
            set
            {
                this.factor = value;
            }
        }
        #endregion

        #region Static voids
        /// <summary>
        /// Создает направленный фильтр с заданными параметрами для билатерального фильтра.
        /// </summary>
        /// <param name="s">Значение первой сигмы σs</param>
        /// <param name="r">Значение второй сигмы σr</param>
        /// <returns>Направленный фильтр</returns>
        public static GuidedFilter FromBilateral(int s, double r)
        {
            return new GuidedFilter(s, r * r);
        }
        #endregion

        #region Public apply voids
        /// <summary>
        /// Реализует одномерный направленный фильтр.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        public void Apply(double[] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int i;

                // guided filter
                double[] copy = (double[])data.Clone();
                GuidedFilter.guidedfilter(copy, this.radius, this.eps);

                // process
                for (i = 0; i < l0; i++)
                    data[i] = (1.0 + this.factor) * (data[i] - copy[i]) + copy[i];
            }

            return;
        }
        /// <summary>
        /// Реализует двумерный направленный фильтр.
        /// </summary>
        /// <param name="data">Матрица</param>
        public void Apply(double[,] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int l1 = data.GetLength(1);
                int i, j;

                // guided filter
                double[,] copy = (double[,])data.Clone();
                GuidedFilter.guidedfilter(copy, this.radius, this.eps);

                // process
                for (i = 0; i < l0; i++)
                    for (j = 0; j < l1; j++)
                        data[i, j] = (1.0 + this.factor) * (data[i, j] - copy[i, j]) + copy[i, j];
            }

            return;
        }
        /// <summary>
        /// Реализует одномерный направленный фильтр.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        public void Apply(Complex[] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int i;

                // guided filter
                Complex[] copy = (Complex[])data.Clone();
                GuidedFilter.guidedfilter(copy, this.radius, this.eps);

                // process
                for (i = 0; i < l0; i++)
                    data[i] = (1.0 + this.factor) * (data[i] - copy[i]) + copy[i];
            }

            return;
        }
        /// <summary>
        /// Реализует двумерный направленный фильтр.
        /// </summary>
        /// <param name="data">Матрица</param>
        public void Apply(Complex[,] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int l1 = data.GetLength(1);
                int i, j;

                // guided filter
                Complex[,] copy = (Complex[,])data.Clone();
                GuidedFilter.guidedfilter(copy, this.radius, this.eps);

                // process
                for (i = 0; i < l0; i++)
                    for (j = 0; j < l1; j++)
                        data[i, j] = (1.0 + this.factor) * (data[i, j] - copy[i, j]) + copy[i, j];
            }

            return;
        }
        #endregion

        #region Blender apply voids
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Набор матриц</param>
        /// <returns>Матрица</returns>
        public double[,] Apply(double[][,] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0), c = data[0].GetLength(1);
            double[,] sum = new double[r, c], cur, low;
            int i, j, k;

            // process
            for (i = 0; i < length; i++)
            {
                // lowpass filter:
                cur = data[i];
                low = (double[,])cur.Clone();
                GuidedFilter.guidedfilter(low, this.radius, this.eps);

                for (j = 0; j < r; j++)
                {
                    for (k = 0; k < c; k++)
                    {
                        // summarize high and low frequencies:
                        sum[j, k] += (1.0 + this.factor) / length * (cur[j, k] - low[j, k]) + low[j, k] / length;

                    }
                }
            }

            return sum;
        }
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Набор матриц</param>
        /// <returns>Матрица</returns>
        public Complex[,] Apply(Complex[][,] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0), c = data[0].GetLength(1);
            Complex[,] sum = new Complex[r, c], cur, low;
            int i, j, k;

            // process
            for (i = 0; i < length; i++)
            {
                // lowpass filter:
                cur = data[i];
                low = (Complex[,])cur.Clone();
                GuidedFilter.guidedfilter(low, this.radius, this.eps);

                for (j = 0; j < r; j++)
                {
                    for (k = 0; k < c; k++)
                    {
                        // summarize high and low frequencies:
                        sum[j, k] += (1.0 + this.factor) / length * (cur[j, k] - low[j, k]) + (low[j, k] / length);

                    }
                }
            }

            return sum;
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Набор векторов</param>
        /// <returns>Одномерный массив</returns>
        public double[] Apply(double[][] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0);
            double[] sum = new double[r], cur, low;
            int i, j;

            // process
            for (i = 0; i < length; i++)
            {
                // lowpass filter:
                cur = data[i];
                low = (double[])cur.Clone();
                GuidedFilter.guidedfilter(low, this.radius, this.eps);

                for (j = 0; j < r; j++)
                {
                    // summarize high and low frequencies:
                    sum[j] += (1.0 + this.factor) / length * (cur[j] - low[j]) + (low[j] / length);
                }
            }

            return sum;
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Набор векторов</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Apply(Complex[][] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0);
            Complex[] sum = new Complex[r], cur, low;
            int i, j;

            // process
            for (i = 0; i < length; i++)
            {
                // lowpass filter:
                cur = data[i];
                low = (Complex[])cur.Clone();
                GuidedFilter.guidedfilter(low, this.radius, this.eps);

                for (j = 0; j < r; j++)
                {
                    // summarize high and low frequencies:
                    sum[j] += (1.0 + this.factor) / length * (cur[j] - low[j]) + (low[j] / length);
                }
            }

            return sum;
        }
        #endregion

        #region Private voids
        // **************************************************
        //                   GUIDED FILTER
        // **************************************************
        // ORIGINALS: Kaiming He, Jian Sun, and Xiaoou Tang.
        // Designed by Asiryan Valeriy (c), 2015-2018
        // Moscow, Russia.
        // **************************************************

        /// <summary>
        /// Guided filer function.
        /// </summary>
        /// <param name="input">Input signal</param>
        /// <param name="r">Filter size</param>
        /// <param name="eps">Epsilon (0, 1)</param>
        internal static void guidedfilter(double[,] input, int r, double eps)
        {
            // Input signal properties:
            int l0 = input.GetLength(0), l1 = input.GetLength(1), i, j;

            // Calculating μ(I) and μ(I^2):
            double[,] x = (double[,])input.Clone();
            double[,] y = Matrice.Pow(input, 2.0);

            // Applying fast box filter:
            MeanFilter.boxf(x, Direction.Both, r);
            MeanFilter.boxf(y, Direction.Both, r);

            // Calculating cov(I):
            // This is the covariance of input in each local patch:
            double[,] c = new double[l0, l1];
            for (i = 0; i < l0; i++)
                for (j = 0; j < l1; j++)
                    c[i, j] = y[i, j] - x[i, j] * x[i, j];

            // Calculating μ(a) and μ(b):
            double[,] a = new double[l0, l1];
            double[,] b = new double[l0, l1];
            for (i = 0; i < l0; i++)
                for (j = 0; j < l1; j++)
                {
                    a[i, j] = c[i, j] / (c[i, j] + eps);
                    b[i, j] = x[i, j] - a[i, j] * x[i, j];
                }

            // Applying fast box filter:
            MeanFilter.boxf(a, Direction.Both, r);
            MeanFilter.boxf(b, Direction.Both, r);

            // Calculating μ(a) * I + μ(b):
            double[,] q = new double[l0, l1];
            for (i = 0; i < l0; i++)
                for (j = 0; j < l1; j++)
                    input[i, j] = a[i, j] * input[i, j] + b[i, j];

            return;
        }
        /// <summary>
        /// Guided filer function.
        /// </summary>
        /// <param name="input">Input signal</param>
        /// <param name="r">Filter size</param>
        /// <param name="eps">Epsilon (0, 1)</param>
        internal static void guidedfilter(Complex[,] input, int r, double eps)
        {
            // Input signal properties:
            int l0 = input.GetLength(0), l1 = input.GetLength(1), i, j;

            // Calculating μ(I) and μ(I^2):
            Complex[,] x = (Complex[,])input.Clone();
            Complex[,] y = Matrice.Pow(input, 2.0);

            // Applying fast box filter:
            MeanFilter.boxf(x, Direction.Both, r);
            MeanFilter.boxf(y, Direction.Both, r);

            // Calculating cov(I):
            // This is the covariance of input in each local patch:
            Complex[,] c = new Complex[l0, l1];
            for (i = 0; i < l0; i++)
                for (j = 0; j < l1; j++)
                    c[i, j] = y[i, j] - x[i, j] * x[i, j];

            // Calculating μ(a) and μ(b):
            Complex[,] a = new Complex[l0, l1];
            Complex[,] b = new Complex[l0, l1];
            for (i = 0; i < l0; i++)
                for (j = 0; j < l1; j++)
                {
                    a[i, j] = c[i, j] / (c[i, j] + eps);
                    b[i, j] = x[i, j] - a[i, j] * x[i, j];
                }

            // Applying fast box filter:
            MeanFilter.boxf(a, Direction.Both, r);
            MeanFilter.boxf(b, Direction.Both, r);

            // Calculating μ(a) * I + μ(b):
            Complex[,] q = new Complex[l0, l1];
            for (i = 0; i < l0; i++)
                for (j = 0; j < l1; j++)
                    input[i, j] = a[i, j] * input[i, j] + b[i, j];

            return;
        }
        /// <summary>
        /// Guided filer function.
        /// </summary>
        /// <param name="input">Input signal</param>
        /// <param name="r">Filter size</param>
        /// <param name="eps">Epsilon (0, 1)</param>
        internal static void guidedfilter(double[] input, int r, double eps)
        {
            // Input signal properties:
            int length = input.Length, i;

            // Calculating μ(I) and μ(I^2):
            double[] x = (double[])input.Clone();
            double[] y = Matrice.Pow(input, 2.0);

            // Applying fast box filter:
            MeanFilter.boxf(x, length, r);
            MeanFilter.boxf(y, length, r);

            // Calculating cov(I):
            // This is the covariance of input in each local patch:
            double[] c = new double[length];
            for (i = 0; i < length; i++)
                c[i] = y[i] - x[i] * x[i];

            // Calculating μ(a) and μ(b):
            double[] a = new double[length];
            double[] b = new double[length];
            for (i = 0; i < length; i++)
            {
                a[i] = c[i] / (c[i] + eps);
                b[i] = x[i] - a[i] * x[i];
            }

            // Applying fast box filter:
            MeanFilter.boxf(a, length, r);
            MeanFilter.boxf(b, length, r);

            // Calculating μ(a) * I + μ(b):
            double[] q = new double[length];
            for (i = 0; i < length; i++)
                input[i] = a[i] * input[i] + b[i];

            return;
        }
        /// <summary>
        /// Guided filer function.
        /// </summary>
        /// <param name="input">Input signal</param>
        /// <param name="r">Filter size</param>
        /// <param name="eps">Epsilon (0, 1)</param>
        internal static void guidedfilter(Complex[] input, int r, double eps)
        {
            // Input signal properties:
            int length = input.Length, i;

            // Calculating μ(I) and μ(I^2):
            Complex[] x = (Complex[])input.Clone();
            Complex[] y = Matrice.Pow(input, 2.0);

            // Applying fast box filter:
            MeanFilter.boxf(x, length, r);
            MeanFilter.boxf(y, length, r);

            // Calculating cov(I):
            // This is the covariance of input in each local patch:
            Complex[] c = new Complex[length];
            for (i = 0; i < length; i++)
                c[i] = y[i] - x[i] * x[i];

            // Calculating μ(a) and μ(b):
            Complex[] a = new Complex[length];
            Complex[] b = new Complex[length];
            for (i = 0; i < length; i++)
            {
                a[i] = c[i] / (c[i] + eps);
                b[i] = x[i] - a[i] * x[i];
            }

            // Applying fast box filter:
            MeanFilter.boxf(a, length, r);
            MeanFilter.boxf(b, length, r);

            // Calculating μ(a) * I + μ(b):
            Complex[] q = new Complex[length];
            for (i = 0; i < length; i++)
                input[i] = a[i] * input[i] + b[i];

            return;
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр на основе пирамиды Лапласа.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// http://www.cs.toronto.edu/~jepson/csc320/notes/pyramids.pdf
    /// </remarks>
    /// </summary>
    public class LaplacianPyramidFilter : IFilter, IBlendFilter
    {
        #region Private data
        private LaplacianPyramidTransform lap;
        private double factor;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр на основе пирамиды Лапласа.
        /// </summary>
        /// <param name="levels">Количество уровней представления</param>
        /// <param name="factor">Множитель [-1, 1]</param>
        public LaplacianPyramidFilter(int levels = 10, double factor = -1.0)
        {
            this.lap = new LaplacianPyramidTransform(levels);
            this.factor = factor;
        }
        /// <summary>
        /// Получает или задает количество уровней представления.
        /// </summary>
        public int Levels
        {
            get
            {
                return lap.Levels;
            }
            set
            {
                lap.Levels = value;
            }
        }
        /// <summary>
        /// Получает или задает значение множителя [-1, 1].
        /// </summary>
        public double Factor
        {
            get
            {
                return this.factor;
            }
            set
            {
                this.factor = value;
            }
        }
        #endregion

        #region Blender apply voids
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Набор матриц</param>
        /// <returns>Матрица</returns>
        public double[,] Apply(double[][,] data)
        {
            int length = data.Length;

            if (length == 0) return null;

            int r = data[0].GetLength(0), c = data[0].GetLength(1);
            double[][,] zero = lap.Forward(new double[r, c]), curr;
            int nlev = zero.Length;

            // process
            for (int j = 0; j < length; j++)
            {
                curr = lap.Forward(data[j]);

                for (int i = 0; i < nlev - 1; i++)
                {
                    zero[i] = Matrice.Add(zero[i], curr[i].Mul((1.0 + this.factor) / length));
                }

                zero[nlev - 1] = Matrice.Add(zero[nlev - 1], curr[nlev - 1].Div(length));
            }

            return lap.Backward(zero);
        }
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Набор матриц</param>
        /// <returns>Матрица</returns>
        public Complex[,] Apply(Complex[][,] data)
        {
            int length = data.Length;

            if (length == 0) return null;

            int r = data[0].GetLength(0), c = data[0].GetLength(1);
            Complex[][,] zero = lap.Forward(new Complex[r, c]), curr;
            int nlev = zero.Length;

            // process
            for (int j = 0; j < length; j++)
            {
                curr = lap.Forward(data[j]);

                for (int i = 0; i < nlev - 1; i++)
                {
                    zero[i] = Matrice.Add(zero[i], curr[i].Mul((1.0 + this.factor) / length));
                }

                zero[nlev - 1] = Matrice.Add(zero[nlev - 1], curr[nlev - 1].Div(length));
            }

            return lap.Backward(zero);
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Набор векторов</param>
        /// <returns>Вектор</returns>
        public double[] Apply(double[][] data)
        {
            int length = data.Length;

            if (length == 0) return null;

            int r = data[0].GetLength(0);
            double[][] zero = lap.Forward(new double[r]), curr;
            int nlev = zero.Length;

            // process
            for (int j = 0; j < length; j++)
            {
                curr = lap.Forward(data[j]);

                for (int i = 0; i < nlev - 1; i++)
                {
                    zero[i] = Matrice.Add(zero[i], curr[i].Mul((1.0 + this.factor) / length));
                }

                zero[nlev - 1] = Matrice.Add(zero[nlev - 1], curr[nlev - 1].Div(length));
            }

            return lap.Backward(zero);
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Набор векторов</param>
        /// <returns>Вектор</returns>
        public Complex[] Apply(Complex[][] data)
        {
            int length = data.Length;

            if (length == 0) return null;

            int r = data[0].GetLength(0);
            Complex[][] zero = lap.Forward(new Complex[r]), curr;
            int nlev = zero.Length;

            // process
            for (int j = 0; j < length; j++)
            {
                curr = lap.Forward(data[j]);

                for (int i = 0; i < nlev - 1; i++)
                {
                    zero[i] = Matrice.Add(zero[i], curr[i].Mul((1.0 + this.factor) / length));
                }

                zero[nlev - 1] = Matrice.Add(zero[nlev - 1], curr[nlev - 1].Div(length));
            }

            return lap.Backward(zero);
        }
        #endregion

        #region Apply voids
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Матрица</param>
        public void Apply(double[,] data)
        {
            // forward pyramid transform
            double[][,] pA = lap.Forward(data);

            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = pA.Length - 1, i, j;

            for (i = 0; i < nlev; i++)
            {
                pA[i] = Matrice.Mul(pA[i], 1.0 + this.factor);
            }

            // backward pyramid transform
            double[,] dummy = lap.Backward(pA);

            for (i = 0; i < r; i++)
            {
                for (j = 0; j < c; j++)
                {
                    data[i, j] = dummy[i, j];
                }
            }

            return;
        }
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Матрица</param>
        public void Apply(Complex[,] data)
        {
            // forward pyramid transform
            Complex[][,] pA = lap.Forward(data);

            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = pA.Length - 1, i, j;

            for (i = 0; i < nlev; i++)
            {
                pA[i] = Matrice.Mul(pA[i], 1.0 + this.factor);
            }

            // backward pyramid transform
            Complex[,] dummy = lap.Backward(pA);

            for (i = 0; i < r; i++)
            {
                for (j = 0; j < c; j++)
                {
                    data[i, j] = dummy[i, j];
                }
            }

            return;
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        public void Apply(double[] data)
        {
            // forward pyramid transform
            double[][] pA = lap.Forward(data);

            int r = data.GetLength(0);
            int nlev = pA.Length - 1, i;

            for (i = 0; i < nlev; i++)
            {
                pA[i] = Matrice.Mul(pA[i], 1.0 + this.factor);
            }

            // backward pyramid transform
            double[] dummy = lap.Backward(pA);

            for (i = 0; i < r; i++)
            {
                data[i] = dummy[i];
            }

            return;
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        public void Apply(Complex[] data)
        {
            // forward pyramid transform
            Complex[][] pA = lap.Forward(data);

            int r = data.GetLength(0);
            int nlev = pA.Length - 1, i;

            for (i = 0; i < nlev; i++)
            {
                pA[i] = Matrice.Mul(pA[i], 1.0 + this.factor);
            }

            // backward pyramid transform
            Complex[] dummy = lap.Backward(pA);

            for (i = 0; i < r; i++)
            {
                data[i] = dummy[i];
            }

            return;
        }
        #endregion
    }
    /// <summary>
    /// Определяет локальный фильтр Лапласа.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://people.csail.mit.edu/sparis/publi/2011/siggraph/
    /// </remarks>
    /// </summary>
    public class LocalLaplacianFilter : IFilter, IBlendFilter
    {
        #region Private data
        /// <summary>
        /// Параметр сигма.
        /// </summary>
        protected double sigma;
        /// <summary>
        /// Множитель.
        /// </summary>
        protected double factor;
        /// <summary>
        /// Количество отсчетов.
        /// </summary>
        protected int n;
        /// <summary>
        /// Количество уровней
        /// </summary>
        protected int levels;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует локальный фильтр Лапласа.
        /// </summary>
        /// <param name="sigma">σ-параметр</param>
        /// <param name="n">Количество отсчетов</param>
        /// <param name="levels">Количество уровней</param>
        /// <param name="factor">Множитель [-1, 1]</param>
        public LocalLaplacianFilter(double sigma = 0.05, int n = 10, int levels = 10, double factor = -1.0)
        {
            this.Sigma = sigma;
            this.N = n;
            this.Levels = levels;
            this.Factor = factor;
        }
        /// <summary>
        /// Получает или задает значение σ-параметра.
        /// </summary>
        public double Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                this.sigma = Maths.Double(value);
            }
        }
        /// <summary>
        /// Получает или задает значение множителя.
        /// </summary>
        public double Factor
        {
            get
            {
                return this.factor;
            }
            set
            {
                this.factor = value;
            }
        }
        /// <summary>
        /// Получает или задает количество отсчетов.
        /// </summary>
        public int N
        {
            get
            {
                return this.n;
            }
            set
            {
                this.n = Math.Max(value, 0);
            }
        }
        /// <summary>
        /// Получает или задает количество уровней.
        /// </summary>
        public int Levels
        {
            get
            {
                return this.levels;
            }
            set
            {
                this.levels = value;
            }
        }
        #endregion

        #region Blender apply voids
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Набор матриц</param>
        /// <returns>Матрица</returns>
        public double[,] Apply(double[][,] data)
        {
            // local laplcian filter blending
            double[,] output = MeanFilter.meanf(data);
            LocalLaplacianFilter.llfilter(output, this.sigma, this.factor, this.n, this.levels);
            return output;
        }
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Набор матриц</param>
        /// <returns>Матрица</returns>
        public Complex[,] Apply(Complex[][,] data)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Набор векторов</param>
        /// <returns>Одномерный массив</returns>
        public double[] Apply(double[][] data)
        {
            // local laplcian filter blending
            double[] output = MeanFilter.meanf(data);
            LocalLaplacianFilter.llfilter(output, this.sigma, this.factor, this.n, this.levels);
            return output;
        }
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Набор векторов</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Apply(Complex[][] data)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Apply voids
        /// <summary>
        /// Реализует двумерный локальный фильтр Лапласа.
        /// </summary>
        /// <param name="data">Матрица</param>
        public void Apply(double[,] data)
        {
            llfilter(data, this.sigma, this.factor, this.n, this.levels);
            return;
        }
        /// <summary>
        /// Реализует одномерный локальный фильтр Лапласа.
        /// </summary>
        /// <param name="data">Матрица</param>
        public void Apply(double[] data)
        {
            llfilter(data, this.sigma, this.factor, this.n, this.levels);
            return;
        }
        /// <summary>
        /// Реализует двумерный локальный фильтр Лапласа.
        /// </summary>
        /// <param name="data">Матрица</param>
        public void Apply(Complex[,] data)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Реализует одномерный локальный фильтр Лапласа.
        /// </summary>
        /// <param name="data">Матрица</param>
        public void Apply(Complex[] data)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Private voids
        // **************************************************
        //            Local Laplacian Filter
        // **************************************************
        // This function implements edge-aware detail and 
        // tone manipulation as described in:
        // "Fast and Robust Pyramid-based Image Processing"
        // Mathieu Aubry, Sylvain Paris, Samuel W. Hasinoff, 
        // Jan Kautz, and Fredo Durand.
        // MIT technical report, November 2011.
        // 
        // Designed by Asiryan Valeriy (c), 2015-2019
        // Moscow, Russia.
        // **************************************************

        /// <summary>
        /// Local laplacian filter.
        /// </summary>
        /// <param name="input">Input data</param>
        /// <param name="sigma">Sigma</param>
        /// <param name="factor">Factor</param>
        /// <param name="n">Number of steps</param>
        /// <param name="levels">Levels</param>
        /// <returns>Output data</returns>
        internal static void llfilter(double[,] input, double sigma, double factor, int n, int levels)
        {
            // exception
            if (factor == 0)
                return;

            // data
            int height = input.GetLength(0);
            int width = input.GetLength(1);
            int y, x, level, length = 256;
            double step = 1.0 / n;
            double min = 0.0, max = 1.0;

            // pyramids
            int n_levels = (int)Math.Min((Math.Log(Math.Min(height, width)) / Math.Log(2)), levels);
            LaplacianPyramidTransform lpt = new LaplacianPyramidTransform(n_levels);
            GaussianPyramidTransform gpt = new GaussianPyramidTransform(n_levels);

            double[][,] input_gaussian_pyr = gpt.Forward(input);
            double[][,] output_laplace_pyr = lpt.Forward(input_gaussian_pyr);
            double[][,] temp_laplace_pyr;
            double[,] I_temp, I_gaus, I_outp;
            double[] T;

            // do job
            for (double i = min; i <= max; i += step)
            {
                height = input.GetLength(0); width = input.GetLength(1);
                I_temp = new double[height, width];
                T = Rem(sigma, factor, i, length);

                // remapping function
                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++)
                    {
                        I_temp[y, x] = T[Maths.Byte(input[y, x] * (length - 1))];
                    }
                }

                temp_laplace_pyr = lpt.Forward(I_temp);
                T = Rec(i, step, length);

                // pyramid reconstruction
                for (level = 0; level < n_levels; level++)
                {
                    I_gaus = input_gaussian_pyr[level];
                    I_temp = temp_laplace_pyr[level];
                    I_outp = output_laplace_pyr[level];
                    height = I_outp.GetLength(0);
                    width = I_outp.GetLength(1);

                    for (y = 0; y < height; y++)
                    {
                        for (x = 0; x < width; x++)
                        {
                            I_outp[y, x] += T[Maths.Byte(I_gaus[y, x] * (length - 1))] * I_temp[y, x];
                        }
                    }

                    output_laplace_pyr[level] = I_outp;
                }
            }

            // backward transform
            I_outp = lpt.Backward(output_laplace_pyr);
            height = input.GetLength(0);
            width = input.GetLength(1);

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++)
                {
                    input[y, x] = I_outp[y, x];
                }
            }

            return;
        }
        /// <summary>
        /// Local laplacian filter.
        /// </summary>
        /// <param name="input">Input data</param>
        /// <param name="sigma">Sigma</param>
        /// <param name="factor">Factor</param>
        /// <param name="n">Number of steps</param>
        /// <param name="levels">Levels</param>
        /// <returns>Output data</returns>
        internal static void llfilter(double[] input, double sigma, double factor, int n, int levels)
        {
            // exception
            if (factor == 0)
                return;

            // data
            int height = input.GetLength(0);
            int y, level, length = 256;
            double step = 1.0 / n;
            double min = 0.0, max = 1.0;

            // pyramids
            int n_levels = (int)Math.Min((Math.Log(height) / Math.Log(2)), levels);
            LaplacianPyramidTransform lpt = new LaplacianPyramidTransform(n_levels);
            GaussianPyramidTransform gpt = new GaussianPyramidTransform(n_levels);

            double[][] input_gaussian_pyr = gpt.Forward(input);
            double[][] output_laplace_pyr = lpt.Forward(input_gaussian_pyr);
            double[][] temp_laplace_pyr;
            double[] I_temp, I_gaus, I_outp;
            double[] T;

            // do job
            for (double i = min; i <= max; i += step)
            {
                height = input.GetLength(0);
                I_temp = new double[height];
                T = Rem(sigma, factor, i, length);

                // remapping function
                for (y = 0; y < height; y++)
                {
                    I_temp[y] = T[Maths.Byte(input[y] * (length - 1))];
                }

                temp_laplace_pyr = lpt.Forward(I_temp);
                T = Rec(i, step, length);

                // pyramid reconstruction
                for (level = 0; level < n_levels; level++)
                {
                    I_gaus = input_gaussian_pyr[level];
                    I_temp = temp_laplace_pyr[level];
                    I_outp = output_laplace_pyr[level];
                    height = I_outp.GetLength(0);

                    for (y = 0; y < height; y++)
                    {
                        I_outp[y] += T[Maths.Byte(I_gaus[y] * (length - 1))] * I_temp[y];
                    }

                    output_laplace_pyr[level] = I_outp;
                }
            }

            // backward transform
            I_outp = lpt.Backward(output_laplace_pyr);
            height = input.GetLength(0);

            for (y = 0; y < height; y++)
            {
                input[y] = I_outp[y];
            }

            return;
        }

        /// <summary>
        /// Reconstruct function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="i">Increment</param>
        /// <param name="step">Step</param>
        /// <returns>Function</returns>
        internal static double Rec(double x, double i, double step)
        {
            double y = Math.Abs(x - i);
            return y < step ? (1.0 - y / step) : 0;
        }
        /// <summary>
        /// Reconstruct function.
        /// </summary>
        /// <param name="i">Increment</param>
        /// <param name="step">Step</param>
        /// <param name="length">Length of table</param>
        /// <returns>Table</returns>
        internal static double[] Rec(double i, double step, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = LocalLaplacianFilter.Rec(x / (double)length, i, step);
            }
            return table;
        }
        /// <summary>
        /// Remapping function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="sigma">Sigma</param>
        /// <param name="factor">Factor</param>
        /// <param name="i">Increment</param>
        /// <returns>Function</returns>
        internal static double Rem(double x, double sigma, double factor, double i)
        {
            double z = 2 * sigma * sigma;
            double y = x - i;
            return factor * y * Math.Exp(-y * y / z);
        }
        /// <summary>
        /// Remapping function.
        /// </summary>
        /// <param name="sigma">Sigma</param>
        /// <param name="factor">Factor</param>
        /// <param name="i">Increment</param>
        /// <param name="length">Length of table</param>
        /// <returns>Table</returns>
        internal static double[] Rem(double sigma, double factor, double i, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = LocalLaplacianFilter.Rem(x / (double)length, sigma, factor, i);
            }
            return table;
        }
        #endregion
    }
    #endregion

    #region Interfaces
    /// <summary>
    /// Определяет общий интерфейс дискретных преобразований.
    /// </summary>
    public interface ITransform
    {
        #region Interface
        /// <summary>
        /// Прямое преобразование.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        double[] Forward(double[] data);
        /// <summary>
        /// Прямое преобразование.
        /// </summary>
        /// <param name="data">Двумерный массив</param>
        /// <returns>Одномерный массив</returns>
        double[,] Forward(double[,] data);
        /// <summary>
        /// Обратное преобразование.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        double[] Backward(double[] data);
        /// <summary>
        /// Обратное преобразование.
        /// </summary>
        /// <param name="data">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        double[,] Backward(double[,] data);
        /// <summary>
        /// Прямое преобразование.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        Complex[] Forward(Complex[] data);
        /// <summary>
        /// Прямое преобразование.
        /// </summary>
        /// <param name="data">Двумерный массив</param>
        /// <returns>Одномерный массив</returns>
        Complex[,] Forward(Complex[,] data);
        /// <summary>
        /// Обратное преобразование.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        Complex[] Backward(Complex[] data);
        /// <summary>
        /// Обратное преобразование.
        /// </summary>
        /// <param name="data">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        Complex[,] Backward(Complex[,] data);
        #endregion
    }
    /// <summary>
    /// Определяет общий вид пирамидоидальных преобразований.
    /// </summary>
    public interface IPyramidTransform
    {
        #region Interface
        /// <summary>
        /// Прямое пирамидоидальное преобразование.
        /// </summary>
        /// <param name="data">Двумерный массив</param>
        /// <returns>Пирамида</returns>
        double[][,] Forward(double[,] data);
        /// <summary>
        /// Обратное пирамидоидальное преобразование.
        /// </summary>
        /// <param name="pyramid">Пирамида</param>
        /// <returns>Двумерный массив</returns>
        double[,] Backward(double[][,] pyramid);
        /// <summary>
        /// Прямое пирамидоидальное преобразование.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        /// <returns>Пирамида</returns>
        double[][] Forward(double[] data);
        /// <summary>
        /// Обратное пирамидоидальное преобразование.
        /// </summary>
        /// <param name="pyramid">Пирамида</param>
        /// <returns>Одномерный массив</returns>
        double[] Backward(double[][] pyramid);
        /// <summary>
        /// Прямое пирамидоидальное преобразование.
        /// </summary>
        /// <param name="data">Двумерный массив</param>
        /// <returns>Пирамида</returns>
        Complex[][,] Forward(Complex[,] data);
        /// <summary>
        /// Обратное пирамидоидальное преобразование.
        /// </summary>
        /// <param name="pyramid">Пирамида</param>
        /// <returns>Двумерный массив</returns>
        Complex[,] Backward(Complex[][,] pyramid);
        /// <summary>
        /// Прямое пирамидоидальное преобразование.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        /// <returns>Пирамида</returns>
        Complex[][] Forward(Complex[] data);
        /// <summary>
        /// Обратное пирамидоидальное преобразование.
        /// </summary>
        /// <param name="pyramid">Пирамида</param>
        /// <returns>Одномерный массив</returns>
        Complex[] Backward(Complex[][] pyramid);
        #endregion
    }
    /// <summary>
    /// Определяет общий вид фильтров.
    /// </summary>
    public interface IFilter
    {
        #region Interface
        /// <summary>
        /// Применяет фильтр к вектору вещественных чисел.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        void Apply(double[] data);
        /// <summary>
        /// Применяет фильтр к матрице вещественных чисел.
        /// </summary>
        /// <param name="data">Матрица</param>
        void Apply(double[,] data);
        /// <summary>
        /// Применяет фильтр к вектору комплексных чисел.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        void Apply(Complex[] data);
        /// <summary>
        /// Применяет фильтр к матрице комплексных чисел.
        /// </summary>
        /// <param name="data">Матрица</param>
        void Apply(Complex[,] data);
        #endregion
    }
    /// <summary>
    /// Определяет общий вид фильтров смешивания.
    /// </summary>
    public interface IBlendFilter
    {
        #region Interface
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Набор матриц</param>
        /// <returns>Матрица</returns>
        double[,] Apply(double[][,] data);
        /// <summary>
        /// Реализует двумерный фильтр.
        /// </summary>
        /// <param name="data">Набор матриц</param>
        /// <returns>Матрица</returns>
        Complex[,] Apply(Complex[][,] data);
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Набор векторов</param>
        /// <returns>Одномерный массив</returns>
        double[] Apply(double[][] data);
        /// <summary>
        /// Реализует одномерный фильтр.
        /// </summary>
        /// <param name="data">Набор векторов</param>
        /// <returns>Одномерный массив</returns>
        Complex[] Apply(Complex[][] data);
        #endregion
    }
    #endregion
}