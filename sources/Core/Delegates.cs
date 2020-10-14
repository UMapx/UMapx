namespace UMapx.Core
{
    /// <summary>
    /// Defines the delegate of a continuous function that depends on a single argument.
    /// </summary>
    /// <param name="x">Argument</param>
    /// <returns>Double precision floating point number</returns>
    public delegate double IDouble(double x);
    /// <summary>
    /// Defines the delegate of a continuous function that depends on a single argument.
    /// </summary>
    /// <param name="x">Argument</param>
    /// <returns>Complex number</returns>
    public delegate Complex IComplex(Complex x);
    /// <summary>
    /// Defines the delegate of a continuous function that depends on two arguments.
    /// </summary>
    /// <param name="x">First argument</param>
    /// <param name="y">Second argument</param>
    /// <returns>Double precision floating point number</returns>
    public delegate double IDoubleMesh(double x, double y);
    /// <summary>
    /// Defines the delegate of a continuous function that depends on two arguments.
    /// </summary>
    /// <param name="x">First argument</param>
    /// <param name="y">Second argument</param>
    /// <returns>Complex number</returns>
    public delegate Complex IComplexMesh(Complex x, Complex y);
}
