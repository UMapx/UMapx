using System;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the levels correction filter.
    /// <remarks>
    /// Filter usage example:
    /// https://digital-photography-school.com/using-levels-photoshop-image-correct-color-contrast/
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LevelsCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private RangeDouble input;
        private RangeDouble output;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the levels correction filter.
        /// </summary>
        /// <param name="input">Input channel values</param>
        /// <param name="output">Output channel values</param>
        /// <param name="space">Color space</param>
        public LevelsCorrection(RangeDouble input, RangeDouble output, Space space)
        {
            Input = input; Output = output; this.Space = space;
        }
        /// <summary>
        /// Initializes the levels correction filter.
        /// </summary>
        public LevelsCorrection()
        {
            Input = new RangeDouble(0, 1);
            Output = new RangeDouble(0, 1);
        }
        /// <summary>
        /// Gets or sets input channel values.
        /// </summary>
        public RangeDouble Input
        {
            get
            {
                return this.input;
            }
            set
            {
                this.input = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Gets or sets output channel values.
        /// </summary>
        public RangeDouble Output
        {
            get
            {
                return this.output;
            }
            set
            {
                this.output = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Levels(input, output, 256);
        }
        #endregion
    }
}
