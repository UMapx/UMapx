using System;
using System.Text.RegularExpressions;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a class of string operations.
    /// </summary>
    internal class StringOptions
    {
        #region String voids
        /// <summary>
        /// Complex number format.
        /// </summary>
        public static string[] C
        {
            get
            {
                return new string[] { "", "i" };
            }
        }
        /// <summary>
        /// Quaternion format.
        /// </summary>
        public static string[] Q
        {
            get
            {
                return new string[] { "i", "j", "k", "" };
            }
        }
        /// <summary>
        /// The function of converting an array of numbers to a string.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="format">Format string</param>
        /// <param name="symbol">String array</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string Disp(float[] v, string format, string[] symbol)
        {
            int length = v.Length, i;
            int start = -1;

            for (i = 0; i < length; i++)
            {
                if (v[i] != 0)
                {
                    start = i;
                    break;
                }
            }

            if (start != -1)
            {
                string result = Disp(v[start], format, true, symbol[start]);

                for (i = start + 1; i < length; i++)
                {
                    result += Disp(v[i], format, false, symbol[i]);
                }
                return result;
            }
            return "0";
        }
        /// <summary>
        /// The function of converting number to a string
        /// </summary>
        /// <param name="v">Number</param>
        /// <param name="format">Format string</param>
        /// <param name="s">First in a row or not</param>
        /// <param name="symbol">Symbol</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static string Disp(double v, string format, bool s, string symbol)
        {
            if (v == 0)
            {
                return "";
            }
            else if (v < 0)
            {
                return (s) ? "-" + (-v).ToString(format) + symbol : " - " + (-v).ToString(format) + symbol;
            }
            return (s) ? v.ToString(format) + symbol : " + " + v.ToString(format) + symbol;
        }
        /// <summary>
        /// Defines a general method for casting the original row to the matrix form.
        /// </summary>
        /// <param name="s">Input string</param>
        /// <returns>String array</returns>
        public static string[] Matpar(string s)
        {
            // example: s = "[1,2,3,4]".
            // Regex options:
            Regex regex = new Regex(@"\[(?<matrice>.*)]", RegexOptions.None);
            Match match = regex.Match(s);

            // success?
            if (match.Success)
            {
                // get new string:
                return Regex.Split(match.Result("${matrice}").Replace(",", "|"), ";");
            }
            throw new Exception("The input string was in the wrong format");
        }
        /// <summary>
        /// Translates the original string to complex number.
        /// <remarks>
        /// Example: "1 + 2i", "0.321 + 11i", ".1i".
        /// </remarks>
        /// </summary>
        /// <param name="s">Input string</param>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public static Complex Compar(string s)
        {
            string u = s.Replace(" ", "");
            int i, k = 0;
            int length = u.Length;
            string re = string.Empty + u[0];
            string im = string.Empty;

            if (!u.Contains("i"))
            {
                for (i = 1; i < length; i++)
                {
                    if (u[i] != '+' && u[i] != '-')
                    {
                        re += u[i]; k++;
                    }
                    else break;
                }

                return new Complex(float.Parse(re), 0);
            }
            else
            {
                if (u == "i") return new Complex(0, 1);

                for (i = 1; i < length; i++)
                {
                    if (u[i] != '+' && u[i] != '-')
                    {
                        re += u[i]; k++;
                    }
                    else break;
                }

                if (k != length - 1)
                {
                    int k1 = k + 1, k2 = k + 2; im += u[k1];

                    if (u[k2] == 'i') return new Complex(float.Parse(re), float.Parse(im + '1'));

                    for (i = k2; i < length; i++)
                    {
                        if (u[i] != 'i')
                        {
                            im += u[i];
                        }
                        else break;
                    }
                    return new Complex(float.Parse(re), float.Parse(im));
                }
                else
                {
                    return new Complex(0, float.Parse(re.Replace("i", "")));
                }
            }
        }
        #endregion
    }
}
