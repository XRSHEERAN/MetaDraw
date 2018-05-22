using System.Linq;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using OxyPlot.Annotations;
using EngineLayer.CrosslinkSearch;
using EngineLayer;
using System.ComponentModel;

namespace ViewModels
{
    public class PdeepModelView : INotifyPropertyChanged
    {
        private PlotModel pDeepModel;

        public PlotModel PdeepModel
        {
            get
            {
                return this.pDeepModel;
            }
            set
            {
                this.pDeepModel = value;
                NotifyPropertyChanged("PdeepModel");
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

        public PdeepModelView()
        {
            // Create the plot model
            var tmp = new PlotModel { Title = "pDeep Spectrum Annotation", Subtitle = "using OxyPlot" };

            // Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            this.PdeepModel = tmp;
        }

        public void UpdataModelForPdeep(PsmCross psmParentsForDraw)
        {
            var x = psmParentsForDraw.MatchedIonInfo.PredictedIonMZ;
            var y = psmParentsForDraw.MatchedIonInfo.PredictedIonIntensity;

            var xM = psmParentsForDraw.MatchedIonInfo.MatchedIonMz;
            var yM = psmParentsForDraw.MatchedIonInfo.MatchedIonIntensity;

            string scanNum = psmParentsForDraw.ScanNumber.ToString();
            string sequence1 = psmParentsForDraw.FullSequence + "-" + psmParentsForDraw.XlPos.ToString();

            PlotModel model = new PlotModel { Title = "pDeep Spectrum anotation of Scan " + scanNum + " for Peptide", DefaultFontSize = 15 };
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "m/z", Minimum = 0, Maximum = x.Max() * 1.02 });
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Intensity(counts)", Minimum = -y.Max()*1.2, Maximum = y.Max() * 1.2 });
            var textAnnoSeq1 = new TextAnnotation() { };
            textAnnoSeq1.FontSize = 12; textAnnoSeq1.TextColor = OxyColors.Red; textAnnoSeq1.StrokeThickness = 0; textAnnoSeq1.TextPosition = new DataPoint(x.Max() / 2, y.Max() * 1.15); textAnnoSeq1.Text = sequence1;
            model.Annotations.Add(textAnnoSeq1);

            LineSeries s0 = new LineSeries();
            s0.Color = OxyColors.Black;
            s0.StrokeThickness = 1;
            s0.Points.Add(new DataPoint(0, 0));
            s0.Points.Add(new DataPoint(x.Max(), 0));
            model.Series.Add(s0);

            LineSeries[] s1 = new LineSeries[x.Length];

            for (int i = 0; i < x.Length; i++)
            {
                OxyColor ionColor = OxyColors.Red;
                if (x[i] > 0)
                {
                    s1[i] = new LineSeries();
                    s1[i].Color = ionColor;
                    s1[i].StrokeThickness = 0.75;
                    s1[i].Points.Add(new DataPoint(x[i] + 1.007277, 0));
                    s1[i].Points.Add(new DataPoint(x[i] + 1.007277, y[i]));

                    var textAnno1 = new TextAnnotation();
                    textAnno1.FontSize = 9;
                    textAnno1.TextColor = ionColor;
                    textAnno1.StrokeThickness = 0;
                    textAnno1.TextPosition = s1[i].Points[1];
                    textAnno1.Text = x[i].ToString("f3");

                    var textAnno2 = new TextAnnotation();
                    textAnno2.FontSize = 9;
                    textAnno2.TextColor = ionColor;
                    textAnno2.StrokeThickness = 0;
                    textAnno2.TextPosition = new DataPoint(s1[i].Points[1].X, s1[i].Points[1].Y + y.Max() * 0.02);
                    textAnno2.Text = psmParentsForDraw.MatchedIonInfo.PredictedIonName[i];

                    model.Annotations.Add(textAnno1);
                    model.Annotations.Add(textAnno2);
                    model.Series.Add(s1[i]);
                }
            }

            //Matched fragment ions
            LineSeries[] s2 = new LineSeries[xM.Length];

            for (int i = 0; i < xM.Length; i++)
            {
                OxyColor ionColor = OxyColors.Blue;
                if (xM[i] > 0)
                {
                    s2[i] = new LineSeries();
                    s2[i].Color = ionColor;
                    s2[i].StrokeThickness = 0.75;
                    s2[i].Points.Add(new DataPoint(xM[i] + 1.007277, 0));
                    s2[i].Points.Add(new DataPoint(xM[i] + 1.007277, - yM[i]/yM.Max()));

                    var textAnno1 = new TextAnnotation();
                    textAnno1.FontSize = 9;
                    textAnno1.TextColor = ionColor;
                    textAnno1.StrokeThickness = 0;
                    textAnno1.TextPosition = s2[i].Points[1];
                    textAnno1.Text = xM[i].ToString("f3");

                    var textAnno2 = new TextAnnotation();
                    textAnno2.FontSize = 9;
                    textAnno2.TextColor = ionColor;
                    textAnno2.StrokeThickness = 0;
                    textAnno2.TextPosition = new DataPoint(s2[i].Points[1].X, s2[i].Points[1].Y - yM.Max()/ yM.Max() * 0.02);
                    textAnno2.Text = psmParentsForDraw.MatchedIonInfo.PredictedIonName[i];

                    model.Annotations.Add(textAnno1);
                    model.Annotations.Add(textAnno2);
                    model.Series.Add(s2[i]);
                }
            }

            // Axes are created automatically if they are not defined

            // Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            this.PdeepModel = model;
        }

    }
}
