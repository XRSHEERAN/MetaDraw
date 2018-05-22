using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using System.Reflection;

namespace MetaDraw
{
    public class Pdeep
    {        
        public static void pDeepSearchInput(List<PsmCross> psmCrosses)
        {
            //Get psms of interest, take base sequences and charge states, and generate example.py
            string modelName = "./h5/2-layer-BiLSTM-QE-M-M.h5"; //trypsin HCD
            //string modelName = "./h5/ZR_FirstPass.h5"; //e001318-20
            string outputPath = Path.Combine(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location), "input.py");       
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(outputPath))
            {
                file.WriteLine("from model.featurize import Ion2Vector");
                file.WriteLine("from keras.models import load_model");
                file.WriteLine("from model import fragmentation_config as fconfig");
                file.WriteLine("import numpy as np");
                file.WriteLine("config = fconfig.HCD_Config()");
                file.WriteLine("pep2vec = Ion2Vector(modini = config.modini, ion_types = config.ion_types, cleave_site_num = config.time_step, max_charge = config.max_ion_charge, prev = 1, next = 1, from_CTerm = config.from_CTerm)");
                file.WriteLine("def output_predict_with_iontype(pred, peptide, charge):");
                file.WriteLine("    pred_charge = charge-1 if charge <= config.max_ion_charge else config.max_ion_charge");
                file.WriteLine("    output = {}");
                file.WriteLine("    for i in range(len(config.ion_types)):");
                file.WriteLine("        it = config.ion_types[i]");
                file.WriteLine("        for ch in range(1, pred_charge+1):");
                file.WriteLine("            output['{}+{}'.format(ch, it)] = pred[:len(peptide)-1, i*config.max_ion_charge + ch-1]");
                file.WriteLine("    return output");
                file.WriteLine("pdeep = load_model('" + modelName + "')");
                file.WriteLine("mod_info = ''");
                //file.WriteLine("path = '" + Path.Combine(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location), "output.txt'"));
                file.WriteLine("path = '/mnt/e/GitHub/MetaDraw/MetaDraw/bin/x64/Debug/output.txt'");
                file.WriteLine("file = open(path,'w')");
                file.WriteLine("newline = " + '"' + "\\n" + '"');
                foreach (var psm in psmCrosses)
                {
                    if (psm.BaseSequence.Length <= 50)
                    {
                        file.WriteLine("peptide = " + "'" + psm.BaseSequence + "'");
                        if (psm.ScanPrecursorCharge != 1)
                            file.WriteLine("charge = " + psm.ScanPrecursorCharge);
                        else
                            file.WriteLine("charge = " + 2);
                        file.WriteLine("x, __ = pep2vec.FeaturizeOnePeptide(peptide, mod_info, charge)");
                        file.WriteLine("pred = pdeep.predict(np.array([x]))[0,:,:]");
                        file.WriteLine("for key, value in output_predict_with_iontype(pred, peptide, charge).items():");
                        file.WriteLine("    print('{}: {}'.format(key, value))");
                        file.WriteLine("for key, value in output_predict_with_iontype(pred, peptide, charge).items():");
                        file.WriteLine("    file.write('{}: {}'.format(key, value))");
                        file.WriteLine("    file.write(newline)");
                        file.WriteLine("file.write(newline)");
                    }
                }
                file.WriteLine("file.close()");
            }
        }

        public static Dictionary<int, double[]> pDeepParser(List<PsmCross> psmCrosses)
        {
            //pdeep has a length limitation of 20? ash Zach
            //import pDeep output
            string outputPath = Path.Combine(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location), "output.txt");
            string[] pDeepLines = File.ReadAllLines(outputPath);
            int pDeepIndex = 0;
            Dictionary<int, double[]> predictedIntensities = new Dictionary<int, double[]>();
            foreach (var psmCross in psmCrosses)
            {
                //get pDeep for this psm
                //assumes two lists for charge 2-, 4 lists for 3+
                //y+1,y+2,b+1,b+2
                List<List<double>> pDeepIntensities = new List<List<double>>();
                string line = "";
                if (pDeepLines[pDeepIndex].Equals(""))
                {
                    pDeepIndex++;
                }
                for (; pDeepIndex < pDeepLines.Length; pDeepIndex++)
                {
                    if (pDeepLines[pDeepIndex].Contains(']'))
                    {
                        line += pDeepLines[pDeepIndex];
                        string originalLine = line;
                        line = line.Replace("[ ", "[");
                        while (line.Contains(" ]"))
                            line = line.Replace(" ]", "]");
                        line = line.Replace("]", "[");
                        line = line.Split('[')[1];
                        while (line.Contains("  "))
                            line = line.Replace("  ", " ");
                        //   try
                        //  {
                        List<string> lineArray = line.Split(' ').ToList();
                        List<double> values = new List<double>();
                        foreach (string s in lineArray)
                        {
                            if (s.Length == 0)
                                continue;
                            double d = Convert.ToDouble(s.Split('e')[0]);
                            if (s.Contains("e"))
                            {
                                string iPre = s.Split('e')[1].Replace("0", "");
                                if (iPre.Length == 1)
                                    iPre = "0";
                                int i = Convert.ToInt16(iPre);
                                while (i > 0)
                                {
                                    i--;
                                    d = d * 10;
                                }
                                while (i < 0)
                                {
                                    i++;
                                    d = d * 0.1;
                                }

                            }
                            values.Add(d);
                        }
                        pDeepIntensities.Add(values);
                        //   }
                        //   catch { }
                        line = "";
                    }
                    else if (pDeepLines[pDeepIndex].Length > 0)
                    {
                        line += pDeepLines[pDeepIndex];
                    }
                    else
                    {
                        //pDeepIndex++;
                        break;
                    }
                }

                //pdeep is n to c, so reverse y such that they match mz order
                for (int i = 0; i < pDeepIntensities.Count / 2; i++)
                {
                    pDeepIntensities[i].Reverse();
                }
                for (int i = 0; i < pDeepIntensities.Count; i++)
                {
                    List<double> currentlist = pDeepIntensities[i];
                    for (int j = 0; j < currentlist.Count; j++)
                        if (currentlist[j] == 0)
                            currentlist[j] = 0.001;
                }
                pDeepIntensities[0].RemoveAt(0);
                predictedIntensities.Add(psmCross.ScanNumber, pDeepIntensities[0].Concat(pDeepIntensities[1]).ToArray());                      
            }
            return predictedIntensities;
        }
    }
}
