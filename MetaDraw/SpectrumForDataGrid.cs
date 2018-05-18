using System.IO;

namespace MetaDraw
{
    class SpectrumForDataGrid
    {
        public SpectrumForDataGrid(int scanNum)
        {
            ScanNum = scanNum;
        }

        public int ScanNum { get; set; }
    }
}
