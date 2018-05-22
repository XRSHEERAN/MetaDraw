using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using Proteomics;
using MzLibUtil;
using UsefulProteomicsDatabases;
using System.IO;
using EngineLayer.CrosslinkSearch;

namespace MetaDraw
{
    [TestFixture]
    class Test
    {
        [SetUpFixture]
        public class MySetUpClass
        {
            private const string elementsLocation = @"Data\elements.dat";

            [OneTimeSetUp]
            public static void Setup()
            {
                Loaders.LoadElements(Path.Combine(TestContext.CurrentContext.TestDirectory, elementsLocation));
            }

        }


        [TestCase]
        public void TestMods()
        {
            PsmCross PSM = new PsmCross();
            PSM.BaseSequence = "DMHGDSEYNIMFGPDICGPGTK";
            PSM.FullSequence = "DM[Common Variable:Oxidation of M]HGDSEYNIM[Common Variable:Oxidation of M]FGPDIC[Common Fixed:Carbamidomethyl of C]GPGTK";
            PSM.PeptideMonisotopicMass = 2472.0032;

            var mods = TsvResultReader.GetMods(PSM);

            Dictionary<int, ModificationWithMass> allModsOneIsNterminus = new Dictionary<int, ModificationWithMass>();

            foreach (var item in mods)
            {
                ModificationMotif.TryGetMotif(item.Value, out ModificationMotif modificationMotif);

                ModificationWithMass mod = new ModificationWithMass("mod", null, modificationMotif, TerminusLocalization.Any, 10);

                allModsOneIsNterminus.Add(item.Key, mod);

            }

            PepWithSetModForCompactPep pepWithSetModForCompactPep = new PepWithSetModForCompactPep();
            pepWithSetModForCompactPep.allModsOneIsNterminus = allModsOneIsNterminus;
            pepWithSetModForCompactPep.BaseSequence = PSM.BaseSequence;
            pepWithSetModForCompactPep.Length = PSM.BaseSequence.Length;
            pepWithSetModForCompactPep.MonoisotopicMass = (double)PSM.PeptideMonisotopicMass;

            var compactPeptide = new CompactPeptide(pepWithSetModForCompactPep, TerminusType.None);

            PSM.CompactPeptide = compactPeptide;

            var lp = new List<ProductType> { ProductType.BnoB1ions, ProductType.Y};
            Tolerance productMassTolerance = new PpmTolerance(10);

            var pmm = PsmCross.XlCalculateTotalProductMassesForSingle(PSM, lp, false);

            var matchedIonMassesListPositiveIsMatch = new MatchedIonInfo(pmm.ProductMz.Length);
            //double pmmScore = PsmCross.XlMatchIons(theScan.TheScan, productMassTolerance, pmm.ProductMz, pmm.ProductName, matchedIonMassesListPositiveIsMatch);
        }

        [Test]
        public void TestpDeepParser()
        {
            Pdeep.pDeepParser(null);
        }
    }
}
