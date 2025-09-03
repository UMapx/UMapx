using System;
using System.Drawing;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Ported from gimp. See gimpmatrix.c and gimp-transform-utils.c.
    /// </summary>
    internal class Float3x3
    {
        #region Private data

        private readonly float[][] coefficients;

        #endregion

        #region Static

        /// <summary>
        /// Identity matrix (cached).
        /// </summary>
        private static readonly Float3x3 Identity;

        /// <summary>
        /// Static initializer builds the identity matrix.
        /// </summary>
        static Float3x3()
        {
            Identity = new Float3x3();
            Identity.coefficients[0][0] = 1;
            Identity.coefficients[1][1] = 1;
            Identity.coefficients[2][2] = 1;
        }

        #endregion

        #region Constructor

        /// <summary>
        /// Creates a zero-initialized 3×3 matrix.
        /// </summary>
        public Float3x3()
        {
            coefficients = new float[3][];
            for (int i = 0; i < 3; i++)
            {
                coefficients[i] = new float[3];
            }
        }

        #endregion

        #region Public methods

        /// <summary>
        /// Builds a perspective transform that maps a source axis-aligned
        /// rectangle to an arbitrary quadrilateral (t1..t4), then normalizes
        /// the source by translating to the rectangle origin and scaling
        /// into [0,1]×[0,1] prior to the warp.
        /// </summary>
        /// <param name="rectangle">Source rectangle in input coordinates</param>
        /// <param name="t1">Destination corner corresponding to source (0,0)</param>
        /// <param name="t2">Destination corner corresponding to source (1,0)</param>
        /// <param name="t3">Destination corner corresponding to source (1,1)</param>
        /// <param name="t4">Destination corner corresponding to source (0,1)</param>
        /// <returns>
        /// Homography matrix H such that dest ≈ H · norm(src), with a perspective divide
        /// </returns>
        /// <remarks>
        /// If the destination quad is (nearly) affine, a simplified affine mapping is used.
        /// Otherwise the full projective warp is constructed. Based on GIMP’s homography derivation.
        /// </remarks>
        public static Float3x3 Perspective(Rectangle rectangle, PointFloat t1, PointFloat t2, PointFloat t3, PointFloat t4)
        {

            var scalex = 1.0f;
            var scaley = 1.0f;

            if (rectangle.Width > 0)
                scalex = 1.0f / rectangle.Width;

            if (rectangle.Height > 0)
                scaley = 1.0f / rectangle.Height;

            var matrix = Float3x3.Identity;
            matrix = matrix.Translate(-rectangle.X, -rectangle.Y);
            matrix = matrix.Scale(scalex, scaley);

            var trafo = new Float3x3();
            {

                var t_x1 = t1.X;
                var t_y1 = t1.Y;
                var t_x2 = t2.X;
                var t_y2 = t2.Y;
                var t_x3 = t3.X;
                var t_y3 = t3.Y;
                var t_x4 = t4.X;
                var t_y4 = t4.Y;
                var dx1 = t_x2 - t_x4;
                var dx2 = t_x3 - t_x4;
                var dx3 = t_x1 - t_x2 + t_x4 - t_x3;
                var dy1 = t_y2 - t_y4;
                var dy2 = t_y3 - t_y4;
                var dy3 = t_y1 - t_y2 + t_y4 - t_y3;

                /*  Is the mapping affine?  */
                var epsilon = 1e-4f;
                if (EpsilonEquals(dx3, 0, epsilon) && EpsilonEquals(dy3, 0, epsilon))
                {
                    trafo.coefficients[0][0] = t_x2 - t_x1;
                    trafo.coefficients[0][1] = t_x4 - t_x2;
                    trafo.coefficients[0][2] = t_x1;
                    trafo.coefficients[1][0] = t_y2 - t_y1;
                    trafo.coefficients[1][1] = t_y4 - t_y2;
                    trafo.coefficients[1][2] = t_y1;
                    trafo.coefficients[2][0] = 0.0f;
                    trafo.coefficients[2][1] = 0.0f;
                }
                else
                {
                    var det1 = dx3 * dy2 - dy3 * dx2;
                    var det2 = dx1 * dy2 - dy1 * dx2;

                    trafo.coefficients[2][0] = EpsilonEquals(det2, 0, epsilon) ? 1.0f : det1 / det2;

                    det1 = dx1 * dy3 - dy1 * dx3;

                    trafo.coefficients[2][1] = EpsilonEquals(det2, 0, epsilon) ? 1.0f : det1 / det2;

                    trafo.coefficients[0][0] = t_x2 - t_x1 + trafo.coefficients[2][0] * t_x2;
                    trafo.coefficients[0][1] = t_x3 - t_x1 + trafo.coefficients[2][1] * t_x3;
                    trafo.coefficients[0][2] = t_x1;

                    trafo.coefficients[1][0] = t_y2 - t_y1 + trafo.coefficients[2][0] * t_y2;
                    trafo.coefficients[1][1] = t_y3 - t_y1 + trafo.coefficients[2][1] * t_y3;
                    trafo.coefficients[1][2] = t_y1;
                }
                trafo.coefficients[2][2] = 1.0f;
            }

            return trafo.Multiply(matrix);
        }

        /// <summary>
        /// Computes the determinant of the 3×3 matrix.
        /// </summary>
        /// <returns>Determinant value.</returns>
        /// <remarks>
        /// Implementation expands along the first column.
        /// </remarks>
        public float Determinant()
        {
            var matrix = this;
            var determinant = (matrix.coefficients[0][0]
                               * (matrix.coefficients[1][1] * matrix.coefficients[2][2] - matrix.coefficients[1][2] * matrix.coefficients[2][1]));
            determinant -= (matrix.coefficients[1][0]
                            * (matrix.coefficients[0][1] * matrix.coefficients[2][2] - matrix.coefficients[0][2] * matrix.coefficients[2][1]));
            determinant += (matrix.coefficients[2][0]
                            * (matrix.coefficients[0][1] * matrix.coefficients[1][2] - matrix.coefficients[0][2] * matrix.coefficients[1][1]));

            return determinant;
        }

        /// <summary>
        /// Returns the inverse matrix using the classical adjugate / determinant formula.
        /// </summary>
        /// <returns>Inverse matrix.</returns>
        /// <remarks>
        /// <para>
        /// Numerical note: if the matrix is singular or nearly singular (|det| ≈ 0),
        /// this operation is unstable. Consider checking the determinant magnitude and
        /// handling the case appropriately.
        /// </para>
        /// <para>
        /// ⚠ NOTE in original code: the conditional
        /// <c>if (EpsilonEquals(det, 0, 1e-5f)) det = 1.0f / det;</c>
        /// is suspicious: typically you would <b>throw</b> when det ≈ 0, and otherwise
        /// set <c>invDet = 1.0f / det</c>. Kept as-is to preserve behavior.
        /// </para>
        /// </remarks>
        public Float3x3 Invert()
        {
            var inv = new Float3x3();
            var matrix = this;

            var det = Determinant();

            if (EpsilonEquals(det, 0, 1e-5f))
                det = 1.0f / det;

            inv.coefficients[0][0] = (matrix.coefficients[1][1] * matrix.coefficients[2][2] -
                                 matrix.coefficients[1][2] * matrix.coefficients[2][1]) * det;

            inv.coefficients[1][0] = -(matrix.coefficients[1][0] * matrix.coefficients[2][2] -
                                 matrix.coefficients[1][2] * matrix.coefficients[2][0]) * det;

            inv.coefficients[2][0] = (matrix.coefficients[1][0] * matrix.coefficients[2][1] -
                                 matrix.coefficients[1][1] * matrix.coefficients[2][0]) * det;

            inv.coefficients[0][1] = -(matrix.coefficients[0][1] * matrix.coefficients[2][2] -
                                 matrix.coefficients[0][2] * matrix.coefficients[2][1]) * det;

            inv.coefficients[1][1] = (matrix.coefficients[0][0] * matrix.coefficients[2][2] -
                                 matrix.coefficients[0][2] * matrix.coefficients[2][0]) * det;

            inv.coefficients[2][1] = -(matrix.coefficients[0][0] * matrix.coefficients[2][1] -
                                 matrix.coefficients[0][1] * matrix.coefficients[2][0]) * det;

            inv.coefficients[0][2] = (matrix.coefficients[0][1] * matrix.coefficients[1][2] -
                                 matrix.coefficients[0][2] * matrix.coefficients[1][1]) * det;

            inv.coefficients[1][2] = -(matrix.coefficients[0][0] * matrix.coefficients[1][2] -
                                 matrix.coefficients[0][2] * matrix.coefficients[1][0]) * det;

            inv.coefficients[2][2] = (matrix.coefficients[0][0] * matrix.coefficients[1][1] -
                                 matrix.coefficients[0][1] * matrix.coefficients[1][0]) * det;

            return inv;
        }

        /// <summary>
        /// Applies the homography to a 2D point (with a perspective divide).
        /// </summary>
        /// <param name="point">Input point</param>
        /// <returns>Transformed point after perspective division</returns>
        /// <remarks>
        /// If the homogeneous w becomes (near) zero, the result is undefined.
        /// This implementation guards only against exact zero; consider using an epsilon.
        /// </remarks>
        public PointFloat TransformPoint(PointFloat point)
        {
            var x = point.X;
            var y = point.Y;
            var matrix = this;

            var w = matrix.coefficients[2][0] * x + matrix.coefficients[2][1] * y + matrix.coefficients[2][2];

            if (w == 0.0)
                w = 1.0f;
            else
                w = 1.0f / w;

            var newx = (matrix.coefficients[0][0] * x +
                     matrix.coefficients[0][1] * y +
                     matrix.coefficients[0][2]) * w;
            var newy = (matrix.coefficients[1][0] * x +
                     matrix.coefficients[1][1] * y +
                     matrix.coefficients[1][2]) * w;
            return new PointFloat(newx, newy);
        }

        /// <summary>
        /// Post-multiplies the current matrix by a translation T(x, y).
        /// </summary>
        /// <param name="x">Translation along X</param>
        /// <param name="y">Translation along Y</param>
        /// <returns>New matrix equal to <c>this * T(x,y)</c></returns>
        public Float3x3 Translate(float x, float y)
        {
            var matrix = Clone();
            var g = matrix.coefficients[2][0];
            var h = matrix.coefficients[2][1];
            var i = matrix.coefficients[2][2];
            matrix.coefficients[0][0] += x * g;
            matrix.coefficients[0][1] += x * h;
            matrix.coefficients[0][2] += x * i;
            matrix.coefficients[1][0] += y * g;
            matrix.coefficients[1][1] += y * h;
            matrix.coefficients[1][2] += y * i;
            return matrix;
        }

        /// <summary>
        /// Post-multiplies the current matrix by a non-uniform scale S(x, y).
        /// </summary>
        /// <param name="x">Scale on X</param>
        /// <param name="y">Scale on Y</param>
        /// <returns>New matrix equal to <c>this * S(x,y)</c></returns>
        public Float3x3 Scale(float x, float y)
        {
            var matrix = Clone();
            matrix.coefficients[0][0] *= x;
            matrix.coefficients[0][1] *= x;
            matrix.coefficients[0][2] *= x;

            matrix.coefficients[1][0] *= y;
            matrix.coefficients[1][1] *= y;
            matrix.coefficients[1][2] *= y;

            return matrix;
        }

        /// <summary>
        /// Deep copy of the matrix.
        /// </summary>
        /// <returns>Cloned matrix</returns>
        public Float3x3 Clone()
        {
            var r = new Float3x3();
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    r.coefficients[i][j] = coefficients[i][j];
                }
            }
            return r;
        }

        /// <summary>
        /// Matrix product: <c>this * other</c>.
        /// </summary>
        /// <param name="other">Right-hand matrix</param>
        /// <returns>Product matrix</returns>
        public Float3x3 Multiply(Float3x3 other)
        {
            var tmp = new Float3x3();

            var matrix1 = this;
            var matrix2 = other;

            for (var i = 0; i < 3; i++)
            {
                var t1 = matrix1.coefficients[i][0];
                var t2 = matrix1.coefficients[i][1];
                var t3 = matrix1.coefficients[i][2];

                for (var j = 0; j < 3; j++)
                {
                    tmp.coefficients[i][j] = t1 * matrix2.coefficients[0][j];
                    tmp.coefficients[i][j] += t2 * matrix2.coefficients[1][j];
                    tmp.coefficients[i][j] += t3 * matrix2.coefficients[2][j];
                }
            }
            return tmp;
        }

        #endregion

        #region Private methods

        /// <summary>
        /// Absolute-ε equality check for floats: |f - other| ≤ ε.
        /// </summary>
        /// <param name="f">First value</param>
        /// <param name="other">Second value</param>
        /// <param name="epsilon">Tolerance</param>
        /// <returns>True if within epsilon</returns>
        internal static bool EpsilonEquals(float f, float other, float epsilon)
        {
            return Math.Abs(f - other) <= epsilon;
        }

        #endregion
    }
}