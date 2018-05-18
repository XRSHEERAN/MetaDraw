﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Proteomics;

namespace MetaDraw
{
    public class PepWithSetModForCompactPep
    {
        public int Length { get; set; }

        public double MonoisotopicMass { get; set; }

        public string BaseSequence { get; set; }

        public Dictionary<int, ModificationWithMass> allModsOneIsNterminus;
    }
}
