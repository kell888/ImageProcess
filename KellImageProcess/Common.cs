using System;
using System.Collections.Generic;
using System.Text;
using System.Configuration;

namespace KellImageProcess
{
    /// <summary>
    /// 公共配置静态类
    /// </summary>
    public static class Common
    {
        /// <summary>
        /// 误差
        /// </summary>
        public static double Eps
        {
            get
            {
                double eps = 0.000001;
                string epsStr = ConfigurationManager.AppSettings["eps"];
                if (!string.IsNullOrEmpty(epsStr))
                {
                    double R;
                    if (double.TryParse(epsStr, out R))
                        eps = R;
                }
                return eps;
            }
        }
        /// <summary>
        /// 扫描位图时分四区域还是九区域
        /// </summary>
        public static int FourOrNine
        {
            get
            {
                int fourOrNine = 4;
                string fourOrNineStr = ConfigurationManager.AppSettings["fourOrNine"];
                if (!string.IsNullOrEmpty(fourOrNineStr))
                {
                    int R;
                    if (int.TryParse(fourOrNineStr, out R))
                        fourOrNine = R;
                }
                return fourOrNine;
            }
        }

    }
}
