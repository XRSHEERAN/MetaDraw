﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using System.Text.RegularExpressions;
using Proteomics;
using MzLibUtil;
using System.Windows;

namespace MetaDraw
{
    public class TsvResultReader
    {
        public static List<PsmCross> ReadTsv(string filePath)
        {
            List<PsmCross> PSMs = new List<PsmCross>();

            Dictionary<string, int> ids = new Dictionary<string, int>();

            //Empty array. Is this correct?
            string[][] resultArray = new string[][] { };
            try
            {
               resultArray = File.ReadLines(filePath).Select(p => p.Split('\t')).ToArray();
            }
            catch (Exception)
            {
                MessageBox.Show("Please check the file.");
                return null;
            }
            

            for (int i = 0; i < resultArray[0].Length; i++)
            {
                    ids.Add(resultArray[0][i], i);         
            }

            for (int i = 1; i < resultArray.Length; i++)
            {
                PsmCross PSM = new PsmCross();
                PSM.ScanNumber = Convert.ToInt32(resultArray[i][ids["Scan Number"]]);
                var baseSeq = resultArray[i][ids["Base Sequence"]];
                var chargeState = resultArray[i][ids["Precursor Charge"]];
                var fullSeq = resultArray[i][ids["Full Sequence"]];
                if (baseSeq.Contains("|"))
                {
                    baseSeq = baseSeq.Split('|').First();
                    fullSeq = fullSeq.Split('|').First();
                }
                PSM.BaseSequence = baseSeq;
                PSM.ScanPrecursorCharge = Convert.ToInt32(chargeState);
                PSM.FullSequence = fullSeq;
                PSM.IsDecoy = (resultArray[i][ids["Decoy"]]=="N" ? false :true);
                PSM.PeptideMonisotopicMass = Convert.ToDouble(resultArray[i][ids["Precursor Mass"]]);

                var mods = GetMods(PSM);

                Dictionary<int, ModificationWithMass> allModsOneIsNterminus = new Dictionary<int, ModificationWithMass>();

                foreach (var item in mods)
                {
                    //I don't really know why use as here.
                    var theMod = GlobalVariables.AllModsKnown.Where(p => p.id == item.Value).First() as ModificationWithMass;                

                    allModsOneIsNterminus.Add(item.Key, theMod);

                }

                PepWithSetModForCompactPep pepWithSetModForCompactPep = new PepWithSetModForCompactPep();
                pepWithSetModForCompactPep.allModsOneIsNterminus = allModsOneIsNterminus;
                pepWithSetModForCompactPep.BaseSequence = PSM.BaseSequence;
                pepWithSetModForCompactPep.Length = PSM.BaseSequence.Length;
                pepWithSetModForCompactPep.MonoisotopicMass = (double)PSM.PeptideMonisotopicMass;

                var compactPeptide = new CompactPeptide(pepWithSetModForCompactPep, TerminusType.None);

                PSM.CompactPeptide = compactPeptide;
                PSMs.Add(PSM);
               
            }

            return PSMs;
        }

        public static Dictionary<int, string> GetMods(PeptideSpectralMatch PSM)
        {
            Dictionary<int, string> modIds = new Dictionary<int, string>();

            var fullseq = PSM.FullSequence;

            string[] parts = fullseq.Split('[', ']');

            int pos = 0;

            Regex regex = new Regex(@":([A-Za-z\s]+)");

            for (int i = 0; i < parts.Length; i++)
            {
                if (parts[i].Contains(":"))
                {
                    Match match = regex.Match(parts[i]);
                    modIds.Add(pos, match.Groups[1].Value);
                }
                else
                {
                    pos += parts[i].Length;
                }
            }

            return modIds;
        }
    }
}
