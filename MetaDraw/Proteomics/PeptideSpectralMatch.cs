using Chemistry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class PeptideSpectralMatch
    {
        #region Private Fields

        private const double tolForDoubleResolution = 1e-6;
        
        private Dictionary<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>> compactPeptides = new Dictionary<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>>();

        #endregion Private Fields

        #region Public Fields

        public const double tolForScoreDifferentiation = 1e-9;

        #endregion Public Fields

        #region Public Constructors

        public PeptideSpectralMatch() { }

        public PeptideSpectralMatch(CompactPeptideBase peptide, int notch, double score, int scanIndex, IScan scan)
        {
            this.ScanIndex = scanIndex;
            this.FullFilePath = scan.FullFilePath;
            this.ScanNumber = scan.OneBasedScanNumber;
            this.PrecursorScanNumber = scan.OneBasedPrecursorScanNumber;
            this.ScanRetentionTime = scan.RetentionTime;
            this.ScanExperimentalPeaks = scan.NumPeaks;
            this.TotalIonCurrent = scan.TotalIonCurrent;
            this.ScanPrecursorCharge = scan.PrecursorCharge;
            this.ScanPrecursorMonoisotopicPeakMz = scan.PrecursorMonoisotopicPeakMz;
            this.ScanPrecursorMass = scan.PrecursorMass;
            AddOrReplace(peptide, score, notch, true);
            this.AllScores = new List<int>(new int[(int)Math.Floor(score) + 1]);
            this.AllScores[AllScores.Count - 1]++;
            MatchedIonDictOnlyMatches = new Dictionary<ProductType, double[]>();
            ProductMassErrorDa = new Dictionary<ProductType, double[]>();
            ProductMassErrorPpm = new Dictionary<ProductType, double[]>();
        }

        #endregion Public Constructors

        #region Public Properties

        public ChemicalFormula ModsChemicalFormula { get; private set; }
        public int ScanNumber { get; set; }
        public int? PrecursorScanNumber { get; set; }
        public double ScanRetentionTime { get; }
        public int ScanExperimentalPeaks { get; }
        public double TotalIonCurrent { get; }
        public int ScanPrecursorCharge { get; set; }
        public double ScanPrecursorMonoisotopicPeakMz { get; set; }
        public double ScanPrecursorMass { get; }
        public string FullFilePath { get; }
        public int ScanIndex { get; }
        public IEnumerable<KeyValuePair<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>>> CompactPeptides { get { return compactPeptides.AsEnumerable(); } }
        public int NumDifferentCompactPeptides { get { return compactPeptides.Count; } }
        public FdrInfo FdrInfo { get; private set; }
        public double Score { get; private set; }
        public double DeltaScore { get; private set; }
        public double RunnerUpScore { get; set; }
        public bool IsDecoy { get; set; }
        public string FullSequence { get; set; }
        public int? Notch { get; private set; }
        public string BaseSequence { get; set; }
        public int? PeptideLength { get; private set; }
        public int? OneBasedStartResidueInProtein { get; private set; }
        public int? OneBasedEndResidueInProtein { get; private set; }
        public double? PeptideMonisotopicMass { get; set; }
        public int? ProteinLength { get; private set; }
        public List<double> LocalizedScores { get; internal set; }
        public Dictionary<ProductType, double[]> MatchedIonDictOnlyMatches { get; internal set; }
        public string ProteinAccesion { get; private set; }
        public string Organism { get; private set; }
        public Dictionary<string, int> ModsIdentified { get; private set; }
        public Dictionary<ProductType, double[]> ProductMassErrorDa { get; internal set; }
        public Dictionary<ProductType, double[]> ProductMassErrorPpm { get; internal set; }

        public List<int> AllScores { get; set; }

        public double[] Features
        {
            get
            {
                return new[] { Math.Round(Score), Score - Math.Round(Score) };
            }
        }

        #endregion Public Properties

        #region Public Methods

        public void AddOrReplace(CompactPeptideBase compactPeptide, double score, int notch, bool reportAllAmbiguity)
        {
            if (score - Score > tolForScoreDifferentiation) //if new score beat the old score, overwrite it
            {
                compactPeptides = new Dictionary<CompactPeptideBase, Tuple<int, HashSet<PeptideWithSetModifications>>>
                {
                    { compactPeptide, new  Tuple<int, HashSet<PeptideWithSetModifications>>(notch,null)}
                };
                if (Score - RunnerUpScore > tolForScoreDifferentiation)
                {
                    RunnerUpScore = Score;
                }
                Score = score;
            }
            else if (score - Score > -tolForScoreDifferentiation && reportAllAmbiguity) //else if the same score and ambiguity is allowed
            {
                compactPeptides[compactPeptide] = new Tuple<int, HashSet<PeptideWithSetModifications>>(notch, null);
            }
            else if (Score - RunnerUpScore > tolForScoreDifferentiation)
            {
                RunnerUpScore = score;
            }
        }

        #endregion Public Methods
        
    }
}