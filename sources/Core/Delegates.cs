namespace UMapx.Core
{
    /// <summary>
    /// Defines the delegate of a continuous function that depends on a single argument.
    /// </summary>
    /// <param name="x">Value</param>
    /// <returns>Value</returns>
    public delegate float IFloat(float x);
    /// <summary>
    /// Defines the delegate of a continuous function that depends on a single argument.
    /// </summary>
    /// <param name="x">Value</param>
    /// <returns>Complex number</returns>
    public delegate ComplexF IComplex32(ComplexF x);
    /// <summary>
    /// Defines the delegate of a continuous function that depends on two arguments.
    /// </summary>
    /// <param name="x">First argument</param>
    /// <param name="y">Second argument</param>
    /// <returns>Value</returns>
    public delegate float IFloatMesh(float x, float y);
    /// <summary>
    /// Defines the delegate of a continuous function that depends on two arguments.
    /// </summary>
    /// <param name="x">First argument</param>
    /// <param name="y">Second argument</param>
    /// <returns>Complex number</returns>
    public delegate ComplexF IComplex32Mesh(ComplexF x, ComplexF y);
}
