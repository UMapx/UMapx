using System;

namespace UMapx.Analysis
{
    /// <summary>
    /// Differentiation method
    /// </summary>
    [Serializable]
    public enum DifferentialMethod
    {
        /// <summary>
        /// Euler method.
        /// </summary>
        Euler,
        /// <summary>
        /// The second-order Runge-Kutta method.
        /// </summary>
        RungeKutta2,
        /// <summary>
        /// Fourth-order Runge-Kutta method.
        /// </summary>
        RungeKutta4,
        /// <summary>
        /// Felberg's method.
        /// </summary>
        Fehlberg,
    }
}
