using System;
using System.Linq;
using MetaDraw;

namespace EngineLayer
{
    [Serializable]
    public class CompactPeptide : CompactPeptideBase
    {
        public CompactPeptide(PeptideWithSetModifications peptideWithSetModifications, TerminusType terminusType)
        {
            NTerminalMasses = null;
            CTerminalMasses = null;
            if (terminusType == TerminusType.None || terminusType == TerminusType.N)
            {
                NTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, 0, 0, 1).ToArray();
            }
            if (terminusType == TerminusType.None || terminusType == TerminusType.C)
            {
                CTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, 0, peptideWithSetModifications.Length + 1, -1).ToArray();
            }
            MonoisotopicMassIncludingFixedMods = peptideWithSetModifications.MonoisotopicMass;
        }

        public CompactPeptide(PepWithSetModForCompactPep peptideWithSetModifications, TerminusType terminusType)
        {
            NTerminalMasses = null;
            CTerminalMasses = null;
            if (terminusType == TerminusType.None || terminusType == TerminusType.N)
            {
                NTerminalMasses = ComputeFollowingFragmentMassesByPepCom(peptideWithSetModifications, 0, 0, 1).ToArray();
            }
            if (terminusType == TerminusType.None || terminusType == TerminusType.C)
            {
                CTerminalMasses = ComputeFollowingFragmentMassesByPepCom(peptideWithSetModifications, 0, peptideWithSetModifications.Length + 1, -1).ToArray();
            }
            MonoisotopicMassIncludingFixedMods = peptideWithSetModifications.MonoisotopicMass;
        }
    }
}