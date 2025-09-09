using System;
using UMapx.Core;

namespace UMapx.Distance
{
    /// <summary>
    /// Defines Minkowski distance.
    /// </summary>
    public class Minkowski : DistanceBase, IDistance
    {
        #region Private data
        private float order = 1.0f;
        #endregion

        #region Constuctor and properties
        /// <summary>
        /// Initializes Minkowski distance.
        /// </summary>
        /// <param name="order">Order (0, +inf)</param>
        public Minkowski(float order)
        {
            this.order = order;
        }
        /// <summary>
        /// Gets or sets order (0, +inf).
        /// </summary>
        public float Order
        {
            get
            {
                return order;
            }
            set
            {
                if (value < 0)
                    throw new ArgumentException("Order cannot be less or equal 0");

                order = value;
            }
        }
        #endregion

        #region Minkowski distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override float Compute(float[] p, float[] q)
        {
            int n = p.Length;
            float sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += Maths.Pow(Math.Abs(p[i] - q[i]), order);
            }
            return Maths.Pow(sum, 1.0f / order);
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public override Complex32 Compute(Complex32[] p, Complex32[] q)
        {
            int n = p.Length;
            float sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += Maths.Pow(Maths.Abs(p[i] - q[i]), order);
            }
            return Maths.Pow(sum, 1.0f / order);
        }
        #endregion
    }
}
