namespace EngineLayer.CrosslinkSearch
{
    public class MatchedIonInfo
    {
        #region Public Constructors

        public MatchedIonInfo(int length)
        {
            MatchedIonMz = new double[length];
            MatchedIonIntensity = new double[length];
            MatchedIonName = new string[length];
            MatchedIonIntensityRank = new int[length];

            PredictedIonMZ = new double[length];
            PredictedIonIntensity = new double[length];
            PredictedIonName = new string[length];
        }

        #endregion Public Constructors

        #region Public Properties

        public double[] MatchedIonMz { get; set; }
        public string[] MatchedIonName { get; set; }
        public double[] MatchedIonIntensity { get; set; }
        public int[] MatchedIonIntensityRank { get; set; }

        public double[] PredictedIonMZ { get; set; }
        public string[] PredictedIonName{ get; set; }
        public double[] PredictedIonIntensity { get; set; }

        #endregion Public Properties
    }
}