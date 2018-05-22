using System;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Collections.ObjectModel;
using System.IO;
using ViewModels;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using System.Collections.Generic;
using MzLibUtil;
using System.Text.RegularExpressions;
using CMD;

namespace MetaDraw
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private readonly ObservableCollection<RawDataForDataGrid> spectraFilesObservableCollection = new ObservableCollection<RawDataForDataGrid>();
        private readonly ObservableCollection<RawDataForDataGrid> resultFilesObservableCollection = new ObservableCollection<RawDataForDataGrid>();
        private MainViewModel mainViewModel = null;
        private PdeepModelView pdeepModelView = null;
        private Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = null;
        private List<PsmCross> PSMs = null;
        private Dictionary<int, double[]> predictedIntensities = null;
        private readonly ObservableCollection<SpectrumForDataGrid> spectrumNumsObservableCollection = new ObservableCollection<SpectrumForDataGrid>();

        public MainWindow()
        {

            InitializeComponent();

            mainViewModel = new MainViewModel();

            plotView.DataContext = mainViewModel;

            pdeepModelView = new PdeepModelView();

            plotViewPdeep.DataContext = pdeepModelView;

            dataGridMassSpectraFiles.DataContext = spectraFilesObservableCollection;

            dataGridResultFiles.DataContext = resultFilesObservableCollection;

            dataGridScanNums.DataContext = spectrumNumsObservableCollection;

            Title = "MetaDraw: version " + GlobalVariables.MetaMorpheusVersion;

            UpdateOutputFolderTextbox();

        }

        private void Window_Drop(object sender, DragEventArgs e)
        {
            if (true)
            {
                string[] files = ((string[])e.Data.GetData(DataFormats.FileDrop)).OrderBy(p => p).ToArray();

                if (files != null)
                {
                    foreach (var draggedFilePath in files)
                    {
                        if (Directory.Exists(draggedFilePath))
                        {
                            foreach (string file in Directory.EnumerateFiles(draggedFilePath, "*.*", SearchOption.AllDirectories))
                            {
                                AddAFile(file);
                            }
                        }
                        else
                        {
                            AddAFile(draggedFilePath);
                        }
                       dataGridMassSpectraFiles.CommitEdit(DataGridEditingUnit.Row, true);
                       dataGridMassSpectraFiles.Items.Refresh();
                    }
                }

            }
        }

        private void btnClearFiles_Click(object sender, RoutedEventArgs e)
        {
            spectraFilesObservableCollection.Clear();
        }

        private void btnAddFiles_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Spectra Files(*.raw;*.mzML)|*.raw;*.mzML",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openFileDialog1.ShowDialog() == true)
                foreach (var rawDataFromSelected in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    AddAFile(rawDataFromSelected);
                }
            dataGridMassSpectraFiles.Items.Refresh();
        }

        private void AddAFile(string draggedFilePath)
        {
            // this line is NOT used because .xml.gz (extensions with two dots) mess up with Path.GetExtension
            //var theExtension = Path.GetExtension(draggedFilePath).ToLowerInvariant();

            // we need to get the filename before parsing out the extension because if we assume that everything after the dot
            // is the extension and there are dots in the file path (i.e. in a folder name), this will mess up
            var filename = Path.GetFileName(draggedFilePath);
            var theExtension = filename.Substring(filename.IndexOf(".")).ToLowerInvariant();

            switch (theExtension)
            {
                case ".raw":
                case ".mzml":
                    RawDataForDataGrid zz = new RawDataForDataGrid(draggedFilePath);
                    if (!SpectraFileExists(spectraFilesObservableCollection, zz)) { spectraFilesObservableCollection.Add(zz); }
                    break;
                case ".pep.XML":
                case ".pep.xml":
                    break;
                case ".psmtsv":
                case ".tsv":
                    RawDataForDataGrid resultFileDataGrid = new RawDataForDataGrid(draggedFilePath);
                    if (!SpectraFileExists(resultFilesObservableCollection, resultFileDataGrid)) { resultFilesObservableCollection.Add(resultFileDataGrid); }
                    break;
                default:
                    break;
            }
        }

        private bool SpectraFileExists(ObservableCollection<RawDataForDataGrid> rDOC, RawDataForDataGrid zzz)
        {
            foreach (RawDataForDataGrid rdoc in rDOC)
                if (rdoc.FileName == zzz.FileName) { return true; }
            return false;
        }

        private void btnRun_Click(object sender, RoutedEventArgs e)
        {
            if (!spectraFilesObservableCollection.Any())
            {
                return;
            }

            DrawEngine drawEngine = new DrawEngine(spectraFilesObservableCollection.Where(b => b.Use).First().FilePath, txtBoxOutputFolder.Text);

            drawEngine.Run();

            arrayOfMs2ScansSortedByMass = drawEngine.arrayOfMs2ScansSortedByMass;

        }

        private void UpdateOutputFolderTextbox()
        {
            if (spectraFilesObservableCollection.Any())
            {
                // if current output folder is blank and there is a spectra file, use the spectra file's path as the output path
                if (string.IsNullOrWhiteSpace(txtBoxOutputFolder.Text))
                {
                    var pathOfFirstSpectraFile = Path.GetDirectoryName(spectraFilesObservableCollection.Where(p => p.Use).First().FilePath);
                    txtBoxOutputFolder.Text = Path.Combine(pathOfFirstSpectraFile, @"$DATETIME");
                }
                // else do nothing (do not override if there is a path already there; might clear user-defined path)
            }
            else
            {
                // no spectra files; clear the output folder from the GUI
                txtBoxOutputFolder.Clear();
            }
        }

        private void btnDraw_Click(object sender, RoutedEventArgs e)
        {
            btnRun.IsEnabled = false;

            mainViewModel.Model.InvalidatePlot(true);
            pdeepModelView.PdeepModel.InvalidatePlot(true);

            int x = Convert.ToInt32(txtScanNum.Text);

            UpdateModel(x);

        }

        private void btnReadResultFile_Click(object sender, RoutedEventArgs e)
        {
            var resultFilePath = resultFilesObservableCollection.Where(b => b.Use).First().FilePath;
            PSMs = TsvResultReader.ReadTsv(resultFilePath);
            foreach (var item in PSMs)
            {
                spectrumNumsObservableCollection.Add(new SpectrumForDataGrid(item.ScanNumber));
            }
            dataGridScanNums.Items.Refresh();
        }

        private void Row_DoubleClick(object sender, System.Windows.Input.MouseButtonEventArgs e)
        {

            if (sender != null)
            {
                try
                {
                    Regex regex = new Regex(@"\d+");
                    int x = Convert.ToInt32(regex.Match(sender.ToString()).Value);
                    UpdateModel(x);

                }
                catch (Exception)
                {
                    MessageBox.Show("Please check the data.");
                }
            }
        }

        private void UpdateModel(int x)
        {
            Ms2ScanWithSpecificMass msScanForDraw = arrayOfMs2ScansSortedByMass.Where(p => p.OneBasedScanNumber == x).First();

            //mainViewModel.UpdateScanModel(msScanForDraw);

            PsmCross psmCross = PSMs.Where(p => p.ScanNumber == x).First();

            var lp = new List<ProductType> { ProductType.BnoB1ions, ProductType.Y };
            Tolerance productMassTolerance = new PpmTolerance(20);

            var pmm = PsmCross.XlCalculateTotalProductMassesForSingle(psmCross, lp, false);

            var matchedIonMassesListPositiveIsMatch = new MatchedIonInfo(pmm.ProductMz.Length);

            double pmmScore = PsmCross.XlMatchIons(msScanForDraw.TheScan, productMassTolerance, pmm.ProductMz, pmm.ProductName, matchedIonMassesListPositiveIsMatch);

            matchedIonMassesListPositiveIsMatch.PredictedIonName = psmCross.ProductMassesMightHaveDuplicatesAndNaNs(lp).ProductName;
            matchedIonMassesListPositiveIsMatch.PredictedIonMZ = psmCross.ProductMassesMightHaveDuplicatesAndNaNs(lp).ProductMz;
            matchedIonMassesListPositiveIsMatch.PredictedIonIntensity = predictedIntensities[psmCross.ScanNumber];

            psmCross.MatchedIonInfo = matchedIonMassesListPositiveIsMatch;


            mainViewModel.UpdateCrosslinkModelForSingle(msScanForDraw, psmCross);

            pdeepModelView.UpdataModelForPdeep(psmCross);

        }

        private void btnPdeep_Click(object sender, RoutedEventArgs e)
        {
            Pdeep.pDeepSearchInput(PSMs);
            CMD.Program.Main(new string[] { });
            predictedIntensities = Pdeep.pDeepParser(PSMs);
        }
    }
}
