using System;

namespace UMapx.Analysis
{
    /// <summary>
    /// Integration method.
    /// </summary>
    [Serializable]
    public enum IntegrationMethod
    {
        /// <summary>
        /// Rectangle method.
        /// </summary>
        Rectangle,
        /// <summary>
        /// Midpoint method.
        /// </summary>
        Midpoint,
        /// <summary>
        /// Trapezoidal method.
        /// </summary>
        Trapezoidal,
        /// <summary>
        /// Simpson method.
        /// </summary>
        Simpson,
        /// <summary>
        /// Romberg method.
        /// </summary>
        Romberg,
    }
}
