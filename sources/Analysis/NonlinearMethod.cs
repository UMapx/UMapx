using System;

namespace UMapx.Analysis
{
    /// <summary>
    /// Method for solving a nonlinear equation.
    /// </summary>
    [Serializable]
    public enum NonlinearMethod
    {
        /// <summary>
        /// Bisection method.
        /// </summary>
        Bisection,
        /// <summary>
        /// Chord method.
        /// </summary>
        Chord,
        /// <summary>
        /// Secant method.
        /// </summary>
        Secant,
        /// <summary>
        /// False position method.
        /// </summary>
        FalsePosition,
    }
}
