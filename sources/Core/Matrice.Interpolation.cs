using System;

namespace UMapx.Core
{
    public static partial class Matrice
    {
        private static void RotationCoefficients(float angle, out double cosine, out double sine)
        {
            double reduced = angle % 360.0;
            if (reduced < 0) reduced += 360;
            if (reduced == 0) { cosine = 1; sine = 0; }
            else if (reduced == 90) { cosine = 0; sine = -1; }
            else if (reduced == 180) { cosine = -1; sine = 0; }
            else if (reduced == 270) { cosine = 0; sine = 1; }
            else
            {
                double radians = -reduced * Math.PI / 180;
                cosine = Math.Cos(radians); sine = Math.Sin(radians);
            }
        }
    }
}
