﻿using System.Linq;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using OxyPlot.Annotations;
using EngineLayer.CrosslinkSearch;
using EngineLayer;
using System.ComponentModel;

namespace ViewModels
{
    public class MainViewModel : INotifyPropertyChanged
    {
        private PlotModel privateModel;

        public PlotModel Model
        {
            get
            {
                return this.privateModel;
            }
            set
            {
                this.privateModel = value;
                NotifyPropertyChanged("Model");
            }
        }

        public event PropertyChangedEventHandler PropertyChanged;

        protected void NotifyPropertyChanged(string propertyName)
        {
            PropertyChangedEventHandler handler = PropertyChanged;
            if (handler != null)
            {
                handler(this, new PropertyChangedEventArgs(propertyName));
            }
        }
    

        public MainViewModel()
        {
            // Create the plot model
            var tmp = new PlotModel { Title = "Spectrum Annotation", Subtitle = "using OxyPlot" };

            // Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            this.privateModel = tmp;
        }

        public void UpdateScanModel(Ms2ScanWithSpecificMass MsScanForDraw)
        {
            var x = MsScanForDraw.TheScan.MassSpectrum.XArray;
            var y = MsScanForDraw.TheScan.MassSpectrum.YArray;

            string scanNum = MsScanForDraw.OneBasedScanNumber.ToString();


            PlotModel model = new PlotModel { Title = "Spectrum anotation of Scan " + scanNum, DefaultFontSize = 15 };
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "m/z", Minimum = 0, Maximum = x.Max() * 1.02 });
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Intensity(counts)", Minimum = 0, Maximum = y.Max() * 1.2 });

            LineSeries[] s0 = new LineSeries[x.Length];
            LineSeries[] s1 = new LineSeries[x.Length];
            LineSeries[] s2 = new LineSeries[x.Length];

            //Draw the ms/ms scan peaks
            for (int i = 0; i < x.Length; i++)
            {
                s0[i] = new LineSeries();
                s0[i].Color = OxyColors.DimGray;
                s0[i].StrokeThickness = 0.15;
                s0[i].Points.Add(new DataPoint(x[i], 0));
                s0[i].Points.Add(new DataPoint(x[i], y[i]));
                model.Series.Add(s0[i]);
            }

            // Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            this.Model = model;

        }

        public void UpdateCrosslinkModelForSingle(Ms2ScanWithSpecificMass MsScanForDraw, PsmCross psmParentsForDraw)
        {
            var x = MsScanForDraw.TheScan.MassSpectrum.XArray;
            var y = MsScanForDraw.TheScan.MassSpectrum.YArray;

            string scanNum = psmParentsForDraw.ScanNumber.ToString();
            string sequence1 = psmParentsForDraw.FullSequence + "-" + psmParentsForDraw.XlPos.ToString();

            var matchedIonDic1 = psmParentsForDraw.MatchedIonInfo;

            PlotModel model = new PlotModel { Title = "Spectrum anotation of Scan " + scanNum + " for Crosslinked Peptide", DefaultFontSize = 15 };
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "m/z", Minimum = 0, Maximum = x.Max() * 1.02 });
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Intensity(counts)", Minimum = 0, Maximum = y.Max() * 1.2 });
            var textAnnoSeq1 = new TextAnnotation() { };
            textAnnoSeq1.FontSize = 9; textAnnoSeq1.TextColor = OxyColors.Red; textAnnoSeq1.StrokeThickness = 0; textAnnoSeq1.TextPosition = new DataPoint(x.Max() / 2, y.Max() * 1.15); textAnnoSeq1.Text = sequence1;
            model.Annotations.Add(textAnnoSeq1);

            LineSeries[] s0 = new LineSeries[x.Length];
            LineSeries[] s1 = new LineSeries[x.Length];


            //Draw the ms/ms scan peaks
            for (int i = 0; i < x.Length; i++)
            {
                s0[i] = new LineSeries();
                s0[i].Color = OxyColors.DimGray;
                s0[i].StrokeThickness = 0.15;
                s0[i].Points.Add(new DataPoint(x[i], 0));
                s0[i].Points.Add(new DataPoint(x[i], y[i]));
                model.Series.Add(s0[i]);
            }
            //Draw the ms/ms scan matched peaks

            for (int i = 0; i < matchedIonDic1.MatchedIonMz.Length; i++)
            {
                OxyColor ionColor = OxyColors.Red;
                if (matchedIonDic1.MatchedIonMz[i] > 0)
                {
                    s1[i] = new LineSeries();
                    s1[i].Color = ionColor;
                    s1[i].StrokeThickness = 0.2;
                    s1[i].Points.Add(new DataPoint(matchedIonDic1.MatchedIonMz[i] + 1.007277, 0));
                    s1[i].Points.Add(new DataPoint(matchedIonDic1.MatchedIonMz[i] + 1.007277, matchedIonDic1.MatchedIonIntensity[i]));

                    var textAnno1 = new TextAnnotation();
                    textAnno1.FontSize = 6;
                    textAnno1.TextColor = ionColor;
                    textAnno1.StrokeThickness = 0;
                    textAnno1.TextPosition = s1[i].Points[1];
                    textAnno1.Text = matchedIonDic1.MatchedIonMz[i].ToString("f3");

                    var textAnno2 = new TextAnnotation();
                    textAnno2.FontSize = 6;
                    textAnno2.TextColor = ionColor;
                    textAnno2.StrokeThickness = 0;
                    textAnno2.TextPosition = new DataPoint(s1[i].Points[1].X, s1[i].Points[1].Y + y.Max() * 0.02);
                    textAnno2.Text = matchedIonDic1.MatchedIonName[i];

                    model.Annotations.Add(textAnno1);
                    model.Annotations.Add(textAnno2);
                    model.Series.Add(s1[i]);
                }

            }

            // Axes are created automatically if they are not defined

            // Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            this.Model = model;
        }

        public void UpdateCrosslinkModel(Ms2ScanWithSpecificMass MsScanForDraw, PsmCross psmParentsForDraw)
        {
            var x = MsScanForDraw.TheScan.MassSpectrum.XArray;
            var y = MsScanForDraw.TheScan.MassSpectrum.YArray;

            string scanNum = psmParentsForDraw.ScanNumber.ToString();
            string sequence1 = psmParentsForDraw.FullSequence + "-" + psmParentsForDraw.XlPos.ToString();
            string sequence2 = psmParentsForDraw.BetaPsmCross.FullSequence + "-" + psmParentsForDraw.BetaPsmCross.XlPos.ToString();

            var matchedIonDic1 = psmParentsForDraw.MatchedIonInfo;
            var matchedIonDic2 = psmParentsForDraw.BetaPsmCross.MatchedIonInfo;

            PlotModel model = new PlotModel { Title = "Spectrum anotation of Scan " + scanNum + " for Crosslinked Peptide", DefaultFontSize = 15 };
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "m/z", Minimum = 0, Maximum = x.Max() * 1.02 });
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Intensity(counts)", Minimum = 0, Maximum = y.Max() * 1.2 });
            var textAnnoSeq1 = new TextAnnotation() { };
            textAnnoSeq1.FontSize = 9; textAnnoSeq1.TextColor = OxyColors.Red; textAnnoSeq1.StrokeThickness = 0; textAnnoSeq1.TextPosition = new DataPoint(x.Max() / 2, y.Max() * 1.15); textAnnoSeq1.Text = sequence1;
            var textAnnoSeq2 = new TextAnnotation() { };
            textAnnoSeq2.FontSize = 9; textAnnoSeq2.TextColor = OxyColors.Blue; textAnnoSeq2.StrokeThickness = 0; textAnnoSeq2.TextPosition = new DataPoint(x.Max() / 2, y.Max() * 1.1); textAnnoSeq2.Text = sequence2;
            model.Annotations.Add(textAnnoSeq1);
            model.Annotations.Add(textAnnoSeq2);

            LineSeries[] s0 = new LineSeries[x.Length];
            LineSeries[] s1 = new LineSeries[x.Length];
            LineSeries[] s2 = new LineSeries[x.Length];

            //Draw the ms/ms scan peaks
            for (int i = 0; i < x.Length; i++)
            {
                s0[i] = new LineSeries();
                s0[i].Color = OxyColors.DimGray;
                s0[i].StrokeThickness = 0.15;
                s0[i].Points.Add(new DataPoint(x[i], 0));
                s0[i].Points.Add(new DataPoint(x[i], y[i]));
                model.Series.Add(s0[i]);
            }
            //Draw the ms/ms scan matched peaks

            for (int i = 0; i < matchedIonDic1.MatchedIonMz.Length; i++)
            {
                OxyColor ionColor = OxyColors.Red;
                if (matchedIonDic1.MatchedIonMz[i] > 0)
                {
                    s1[i] = new LineSeries();
                    s1[i].Color = ionColor;
                    s1[i].StrokeThickness = 0.2;
                    s1[i].Points.Add(new DataPoint(matchedIonDic1.MatchedIonMz[i] + 1.007277, 0));
                    s1[i].Points.Add(new DataPoint(matchedIonDic1.MatchedIonMz[i] + 1.007277, matchedIonDic1.MatchedIonIntensity[i]));

                    var textAnno1 = new TextAnnotation();
                    textAnno1.FontSize = 6;
                    textAnno1.TextColor = ionColor;
                    textAnno1.StrokeThickness = 0;
                    textAnno1.TextPosition = s1[i].Points[1];
                    textAnno1.Text = matchedIonDic1.MatchedIonMz[i].ToString("f3");

                    var textAnno2 = new TextAnnotation();
                    textAnno2.FontSize = 6;
                    textAnno2.TextColor = ionColor;
                    textAnno2.StrokeThickness = 0;
                    textAnno2.TextPosition = new DataPoint(s1[i].Points[1].X, s1[i].Points[1].Y + y.Max() * 0.02);
                    textAnno2.Text = matchedIonDic1.MatchedIonName[i];

                    model.Annotations.Add(textAnno1);
                    model.Annotations.Add(textAnno2);
                    model.Series.Add(s1[i]);
                }

            }

            for (int i = 0; i < matchedIonDic2.MatchedIonMz.Length; i++)
            {
                OxyColor ionColor = OxyColors.Blue;
                if (matchedIonDic2.MatchedIonMz[i] > 0)
                {
                    s2[i] = new LineSeries();
                    s2[i].Color = ionColor;
                    s2[i].StrokeThickness = 0.2;
                    s2[i].Points.Add(new DataPoint(matchedIonDic2.MatchedIonMz[i] + 1.007277, 0));
                    s2[i].Points.Add(new DataPoint(matchedIonDic2.MatchedIonMz[i] + 1.007277, matchedIonDic2.MatchedIonIntensity[i]));

                    var textAnno1 = new TextAnnotation();
                    textAnno1.FontSize = 6;
                    textAnno1.TextColor = ionColor;
                    textAnno1.StrokeThickness = 0;
                    textAnno1.TextPosition = s2[i].Points[1];
                    textAnno1.Text = matchedIonDic2.MatchedIonMz[i].ToString("f3");

                    var textAnno2 = new TextAnnotation();
                    textAnno2.FontSize = 6;
                    textAnno2.TextColor = ionColor;
                    textAnno2.StrokeThickness = 0;
                    textAnno2.TextPosition = new DataPoint(s2[i].Points[1].X, s2[i].Points[1].Y + y.Max() * 0.02);
                    textAnno2.Text = matchedIonDic2.MatchedIonName[i];

                    model.Annotations.Add(textAnno1);
                    model.Annotations.Add(textAnno2);
                    model.Series.Add(s2[i]);
                }
            }

            // Axes are created automatically if they are not defined

            // Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            this.Model = model;
        }

    }
}
