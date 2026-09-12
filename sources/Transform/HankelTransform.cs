using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Hankel transform.
    /// </summary>
    /// <remarks>
    /// NOT RECOMMENDED.
    /// There is no fast O(N log N) algorithm for this transform.
    /// 
    /// More information can be found on the website:
    /// <see href="https://en.wikipedia.org/wiki/Hankel_transform"/>.
    /// </remarks>
    [Serializable]
    public class HankelTransform : TransformBaseMatrixFloat, ITransform
    {
        #region Private data
        /// <summary>
        /// Param.
        /// </summary>
        private int a;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Hankel transform.
        /// </summary>
        /// <param name="a">Param.</param>
        /// <param name="direction">Processing direction.</param>
        public HankelTransform(int a = 0, Direction direction = Direction.Vertical)
        {
            this.Direction = direction;
            this.a = a;
        }
        /// <summary>
        /// Gets or sets the param of transform.
        /// </summary>
        public int A
        {
            get
            { 
                return this.a; 
            }
            set
            { 
                this.a = value; 
            }
        }
        #endregion

        #region Hankel static components
        /// <summary>
        /// Implements the construction of the Hankel transform matrix.
        /// </summary>
        /// <param name="N">Size.</param>
        /// <param name="a">Param.</param>
        /// <returns>Matrix.</returns>
        /// <exception cref="ArgumentException">Exception.</exception>
        public static float[,] Matrix(int N, int a)
        {
            if (N <= 0 || a < 0) throw new ArgumentException("Arguments could not be negative");

            float[] j = BesselZerosJ(a, N + 1);

            float[] Jp = new float[N + 1 + 1];

            for (int k = 1; k <= N + 1; k++)
                Jp[k] = Special.J(j[k], a + 1);

            float denom = j[N + 1];

            var T = new float[N, N];

            for (int m = 1; m <= N; m++)
            {
                for (int n = 1; n <= N; n++)
                {
                    float arg = (float)((double)j[m] * j[n] / denom);
                    float num = Special.J((float)arg, a);
                    float val = 2.0f / denom * num / (Jp[m] * Jp[n]);
                    T[m - 1, n - 1] = val;
                }
            }
            return T;
        }
        /// <summary>
        /// Finds consecutive positive zeros of the integer-order Bessel function J.
        /// </summary>
        /// <param name="a">Nonnegative integer order.</param>
        /// <param name="count">Number of positive zeros required.</param>
        /// <returns>Increasing roots in entries 1 through count; entry 0 is unused.</returns>
        /// <exception cref="ArithmeticException">The roots cannot be separated at single precision.</exception>
        private static float[] BesselZerosJ(int a, int count)
        {
            float[] roots = new float[count + 1];
            // The first positive zero lies above the order. Scanning with pi/4
            // cannot skip consecutive zeros at nonnegative integer orders and
            // excludes the zero at the origin when a > 0.
            float left = a, fLeft = Special.J(left, a);
            int found = 0;
            int scanLimit = checked(16 * count + 128 + (int)(16 * Math.Pow(a, 1.0 / 3)));
            for (int scan = 0; scan < scanLimit && found < count; scan++)
            {
                float right = (float)(left + Math.PI / 4);
                float fRight = Special.J(right, a);
                if (right <= left || float.IsNaN(fRight) || float.IsNaN(fLeft)) break;
                if (fRight == 0 || (fLeft != 0 && (fLeft < 0) != (fRight < 0)))
                {
                    float lo = left, hi = right, flo = fLeft;
                    if (fRight != 0)
                    {
                        for (int iteration = 0; iteration < 32; iteration++)
                        {
                            float mid = (float)(((double)lo + hi) / 2);
                            if (mid == lo || mid == hi) break;
                            float fm = Special.J(mid, a);
                            if (fm == 0) { lo = hi = mid; break; }
                            if ((fm < 0) == (flo < 0)) { lo = mid; flo = fm; }
                            else hi = mid;
                        }
                    }
                    else lo = hi;
                    roots[++found] = (float)(((double)lo + hi) / 2);
                }
                left = right;
                fLeft = fRight;
            }
            if (found != count) throw new ArithmeticException("Unable to bracket distinct positive Bessel zeros.");
            return roots;
        }
        #endregion

        #region Hankel Transform
        /// <inheritdoc/>
        protected override float[,] TransformationMatrix(int n)
        {
            return Matrix(n, this.a);
        }
        #endregion
    }
}
