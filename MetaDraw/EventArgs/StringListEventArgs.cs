using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MetaDraw
{
    class StringListEventArgs : EventArgs
    {
        #region Public Constructors

        public StringListEventArgs(List<string> stringList)
        {
            this.StringList = stringList;
        }

        #endregion Public Constructors

        #region Public Properties

        public IEnumerable<string> StringList { get; }

        #endregion Public Properties
    }
}
