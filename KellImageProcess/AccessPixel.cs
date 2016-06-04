using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.Windows.Forms;

namespace KellImageProcess
{
    /// <summary>
    /// 像素访问处理类
    /// </summary>
    public class AccessPixel
    {
        int BPP = 4;
        Bitmap bmp;
        int width;
        int height;
        static uint[] histogram = new uint[256];
        static byte currentThreshold;
        static List<Point> tmpPoints;
        static bool isLinkCZ = false;
        static readonly object lockObj = new object();

        /// <summary>
        /// 默认构造函数
        /// </summary>
        public AccessPixel()
        {
            tmpPoints = new List<Point>();
        }
        /// <summary>
        /// 初始化位图构造函数
        /// </summary>
        /// <param name="bitmap"></param>
        public AccessPixel(Bitmap bitmap)
        {
            tmpPoints = new List<Point>();
            if (bitmap == null)
            {
                throw new KellImageProcessException("The Initial bitmap is null !");
            }
            bmp = (Bitmap)bitmap.Clone();
            BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            width = bitmap.Width;
            height = bitmap.Height;
            currentThreshold = GetThresholdAuto(bmp);
        }

        /// <summary>
        /// 将指定位图的背景以指定的背景颜色及颜色容差进行透明化
        /// </summary>
        /// <param name="bmp">指定位图</param>
        /// <param name="backgroundColor">指定的背景颜色</param>
        /// <param name="tolerance">颜色容差，默认为0</param>
        /// <returns></returns>
        public static Bitmap TransparentBackground(Bitmap bmp, Color backgroundColor, byte tolerance)
        {
            if (bmp == null)
                return null;
            int width = bmp.Width;
            int height = bmp.Height;
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            if (BPP < 3)
                return null;
            Bitmap b1 = new Bitmap(width, height, PixelFormat.Format32bppArgb);
            unsafe
            {
                BitmapData data = bmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, bmp.PixelFormat);
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * width;

                BitmapData data1 = b1.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, PixelFormat.Format32bppArgb);
                byte* p1 = (byte*)data1.Scan0;
                int stride1 = data1.Stride;
                int offset1 = stride1 - 4 * width;

                for (int i = 0; i < width; i++)
                {
                    for (int j = 0; j < height; j++)
                    {
                        p1[0] = p[0];
                        p1[1] = p[1];
                        p1[2] = p[2];

                        double distance = Math.Sqrt((backgroundColor.R - p[2]) * (backgroundColor.R - p[2]) + (backgroundColor.G - p[1]) * (backgroundColor.G - p[1]) + (backgroundColor.B - p[0]) * (backgroundColor.B - p[0]));
                        if (distance <= tolerance)
                        {
                            p1[3] = 0;
                        }
                        else
                        {
                            if (BPP > 3)
                                p1[3] = p[3];
                            else
                                p1[3] = 255;
                        }
                        p += BPP;
                        p1 += 4;
                    }
                    p += offset;
                    p1 += offset1;
                }
                bmp.UnlockBits(data);
                b1.UnlockBits(data1);
            }
            return b1;
        }

        /// <summary>
        /// 获取指定带有位置信息像素数组的反色信息
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static List<PointColor> GetReverseColor(List<PointColor> ps)
        {
            List<PointColor> ret = new List<PointColor>();
            foreach (PointColor p in ps)
            {
                int R = p.Color.R ^ 0xFF;
                int G = p.Color.G ^ 0xFF;
                int B = p.Color.B ^ 0xFF;
                PointColor pc = new PointColor();
                pc.Location = p.Location;
                pc.Color = Color.FromArgb(R, G, B);
                ret.Add(pc);
            }
            return ret;
        }
        /// <summary>
        /// 去除点集中的异常点，且插值趋势点(注：errPs.Count不一定等于interPs.Count，即errPs.Count>=interPs.Count)
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="badDegree">点集粒度，常用范围[3,5,7,9,11]，粒度越小插值越细，最佳的经验值为[3,5]</param>
        /// <param name="factor">波动系数，常用范围[0,1]，系数越小插值趋势点越逼近原始点，最佳的经验值范围[0.01,0.1]</param>
        /// <param name="containsInterPs"></param>
        /// <param name="errPs"></param>
        /// <param name="interPs"></param>
        /// <returns></returns>
        public static List<Point> RemoveExceptionPointAndInterpolation(List<Point> ps, int badDegree, double factor, bool containsInterPs, out List<Point> errPs, out List<Point> interPs)
        {
            List<Point> ps1 = new List<Point>();
            interPs = new List<Point>();
            errPs = new List<Point>();
            if (ps.Count > 2)
            {
                errPs = GetBadPoints(ps, badDegree, factor);
                double[] xs = new double[ps.Count - errPs.Count];
                double[] ys = new double[ps.Count - errPs.Count];
                int valid = 0;
                for (int i = 0; i < ps.Count; i++)
                {
                    if (!errPs.Contains(ps[i]))
                    {
                        xs[valid] = ps[i].X;
                        ys[valid] = ps[i].Y;
                        ps1.Add(ps[i]);
                        valid++;
                    }
                }
                for (int i = 0; i < errPs.Count; i++)
                {
                    Point p1, p2;
                    Point cz = Point.Empty;
                    try
                    {
                        cz = GetShortDisP(errPs[i], errPs, ps1, out p1, out p2);
                    }
                    catch
                    {
                        break;
                    }
                    double t = cz.X;
                    double yt = 0;
                    try
                    {
                        yt = CSharpAlgorithm.Algorithm.Interpolation.GetValuePqs(ps.Count - errPs.Count, xs, ys, t);
                    }
                    catch
                    {
                        yt = cz.Y;
                    }
                    if (!double.IsNaN(yt) && Math.Abs(yt) < short.MaxValue)//所以errPs.Count不一定等于interPs.Count！
                    {
                        Point p = new Point(Convert.ToInt32(Math.Round(t)), Convert.ToInt32(Math.Round(yt)));
                        interPs.Add(p);
                        if (containsInterPs)
                            ps1.Add(p);
                    }
                }
            }
            else if (ps.Count == 1 || ps.Count == 2)
            {
                ps1.Add(ps[0]);
                ps1.Add(ps[1]);
            }
            return ps1;
        }

        /// <summary>
        /// 求指定异常点到两个临近点的垂足，且返回该两个临近点
        /// </summary>
        /// <param name="p"></param>
        /// <param name="errPs"></param>
        /// <param name="goodPs"></param>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static Point GetShortDisP(Point p, List<Point> errPs, List<Point> goodPs, out Point p1, out Point p2)
        {
            p1 = Point.Empty;
            p2 = Point.Empty;
            if (goodPs.Count > 0)
            {
                //是否需要移除正常点集中与异常点p位置相同的点？
                //while (goodPs.Contains(p))
                //{
                //    goodPs.Remove(p);
                //}
                Point sp = goodPs[0];
                double sd = Distance2D(p, sp);
                Point lastSp = sp;
                for (int i = 1; i < goodPs.Count; i++)
                {
                    double dis = Distance2D(p, goodPs[i]);
                    if (dis < sd)
                    {
                        sd = dis;
                        lastSp = sp;
                        sp = goodPs[i];
                    }
                    else
                    {
                        if (dis < Distance2D(p, lastSp))
                        {
                            lastSp = goodPs[i];
                        }
                    }
                }
                p1 = lastSp;
                p2 = sp;
                if (p1.X == p2.X && p1.Y == p2.Y)//异常点在正常点的两侧（即：只有一个临近点，所以要插值出第二个临近点）
                {
                    Rectangle rect = RectangleOnMorePoint(goodPs);
                    Point goodCenter = new Point(rect.X + rect.Width / 2, rect.Y + rect.Height / 2);
                    double maxDis = Distance2D(errPs[0], goodCenter);
                    Point longDisP = errPs[0];
                    for (int i = 1; i < errPs.Count; i++)
                    {
                        double dis = Distance2D(errPs[i], goodCenter);
                        if (dis > maxDis)
                        {
                            maxDis = dis;
                            longDisP = errPs[i];
                        }
                    }
                    p2 = longDisP;
                }
                return Pvint(p, p1, p2);
            }
            else
            {
                throw new Exception("正常点集为空！");
            }
        }

        /// <summary>
        /// 过直线外一点的在直线上的垂足
        /// </summary>
        /// <param name="pt">直线外的一点</param>
        /// <param name="lt1">直线上的一点</param>
        /// <param name="lt2">直线上的另一点</param>
        /// <returns></returns>
        public static Point Pvint(Point pt, Point lt1, Point lt2)
        {
            Point t = new Point();
            double x = 0.0;
            double y = 0.0;
            if (lt1.X != lt2.X)
            {
                if (lt1.Y != lt2.Y)
                {
                    double k1 = 0.0;
                    double k2 = 0.0;
                    k1 = (lt2.Y - lt1.Y) / ((double)lt2.X - lt1.X);
                    k2 = -1 / k1;//垂直线的斜率

                    x = ((lt1.Y * lt2.X - lt2.Y * lt1.X) / (lt2.X - lt1.X) + k2 * pt.X - pt.Y) / (k2 - k1);
                    y = k2 * x - k2 * pt.X + pt.Y;
                }
                else
                {
                    x = pt.X;//与X平行
                    y = lt1.Y;
                }
            }
            else
            {
                x = lt1.X;//与Y平行
                y = pt.Y;
            }
            t.X = (int)x;
            t.Y = (int)y;
            return t;
        }

        /// <summary>
        /// 获取指定点集中的异常点(集)
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="badDegree">点集粒度，常用范围[3,5,7,9,11]，粒度越小插值越细，最佳的经验值为[3,5]</param>
        /// <param name="factor">波动系数，常用范围[0,1]，系数越小插值趋势点越逼近原始点，最佳的经验值范围[0.01,0.1]</param>
        /// <returns></returns>
        public static List<Point> GetBadPoints(List<Point> ps, int badDegree, double factor)
        {
            List<Point> err = new List<Point>();
            //以下是最为关键的代码：
            int points = ps.Count / badDegree;
            Point[] tmp = new Point[badDegree];
            Hashtable ht = new Hashtable();
            int max = 0;
            double current = 0;
            List<int> currentIndexs = new List<int>();
            List<int> goodIndexs = new List<int>();
            List<DictionaryEntry> des = new List<DictionaryEntry>();
            for (int i = 0; i < points; i++)
            {
                ps.CopyTo(badDegree * i, tmp, 0, badDegree);
                double ds = GetDSN(tmp, tmp[tmp.Length / 2], 2);
                ht.Add(i, ds);
                if (Math.Abs(ds - current) > factor)
                {
                    for (int j = 0; j < des.Count; j++)
                    {
                        if (Math.Abs((double)des[j].Key - current) <= factor)
                            des[j] = new DictionaryEntry(current, (int)des[j].Value + currentIndexs.Count);
                        else
                            des.Add(new DictionaryEntry(current, currentIndexs.Count));
                    }
                    current = ds;
                    currentIndexs.Clear();
                    currentIndexs.Add(i);
                }
                else
                {
                    currentIndexs.Add(i);
                }
            }
            double maxCountDs = 0;
            foreach (DictionaryEntry de in des)
            {
                if (max < (int)de.Value)
                {
                    max = (int)de.Value;
                    maxCountDs = (double)de.Key;
                }
            }
            foreach (object ind in ht.Keys)
            {
                if (Math.Abs((double)ht[ind] - maxCountDs) <= factor)
                    goodIndexs.Add((int)ind);
            }
            foreach (object key in ht.Keys)
            {
                int i = (int)key;
                if (!goodIndexs.Contains(i))
                {
                    for (int j = 0; j < badDegree; j++)
                    {
                        err.Add(ps[badDegree * i + j]);
                    }
                }
            }
            return err;
        }

        /// <summary>
        /// 求ceof系数决定函数所对应x自变量的函数值
        /// </summary>
        /// <param name="ceof"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double F(double[] ceof, double x)
        {
            double sum = 0;
            for (int i = ceof.Length - 1; i >= 0; i--)
            {
                sum += ceof[i] * Math.Pow(x, i);
            }
            return sum;
        }

        /// <summary>
        /// First-Derivative
        /// </summary>
        /// <param name="ceof"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double GetDS(double[] ceof, double x)
        {
            double e = 0.000000001;
            double sum = 0;
            sum = F(ceof, x + e) - F(ceof, x);
            double ds = sum / e;
            return ds;
        }

        /// <summary>
        /// N-Derivative
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="p"></param>
        /// <param name="js"></param>
        /// <returns></returns>
        public static double GetDSN(List<Point> ps, Point p, int js)
        {
            double[] dCeof = new double[4];
            while (js > 1)
            {
                dCeof = GetDerivativeFunction(ps);
                js--;
            }
            return GetDS(dCeof, p.X);
        }

        /// <summary>
        /// 求导函数，返回导函数的系数数组(三次函数的系数数组，一维且长度为4)
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static double[] GetDerivativeFunction(List<Point> ps)
        {
            if (ps.Count > 1)
            {
                int points = ps.Count;
                double[] xs = new double[points];
                double[] ys = new double[points];
                double[] s = new double[5];
                for (int i = 0; i < points; i++)
                {
                    xs[i] = ps[i].X;
                    ys[i] = ps[i].Y;
                }
                CSharpAlgorithm.Algorithm.Interpolation.GetValueAkima(points, xs, ys, ps[0].X, s, -1);
                for (int i = 0; i < points; i++)
                {
                    double x = ps[i].X;
                    double[] ceof = new double[4];
                    for (int j = 0; j < 4; j++)
                    {
                        ceof[j] = s[j];
                    }
                    ys[i] = GetDS(ceof, x);
                }
                double y = CSharpAlgorithm.Algorithm.Interpolation.GetValueAkima(points, xs, ys, ps[0].X, s, -1);
                double[] dCeof = new double[4];
                for (int i = 0; i < 4; i++)
                {
                    dCeof[i] = s[i];
                }
                return dCeof;
            }
            else
            {
                throw new KellImageProcessException("点数少于2点！");
            }
        }

        /// <summary>
        /// N-Derivative
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="p"></param>
        /// <param name="js"></param>
        /// <returns></returns>
        public static double GetDSN(Point[] ps, Point p, int js)
        {
            double[] dCeof = new double[4];
            while (js > 1)
            {
                dCeof = GetDerivativeFunction(ps);
                js--;
            }
            return GetDS(dCeof, p.X);
        }

        /// <summary>
        /// 求导函数，返回导函数的系数数组(三次函数的系数数组，一维且长度为4)
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static double[] GetDerivativeFunction(Point[] ps)
        {
            if (ps.Length > 1)
            {
                int points = ps.Length;
                double[] xs = new double[points];
                double[] ys = new double[points];
                double[] s = new double[5];
                for (int i = 0; i < points; i++)
                {
                    xs[i] = ps[i].X;
                    ys[i] = ps[i].Y;
                }
                CSharpAlgorithm.Algorithm.Interpolation.GetValueAkima(points, xs, ys, ps[0].X, s, -1);
                for (int i = 0; i < points; i++)
                {
                    double x = ps[i].X;
                    double[] ceof = new double[4];
                    for (int j = 0; j < 4; j++)
                    {
                        ceof[j] = s[j];
                    }
                    ys[i] = GetDS(ceof, x);
                }
                double y = CSharpAlgorithm.Algorithm.Interpolation.GetValueAkima(points, xs, ys, ps[0].X, s, -1);
                double[] dCeof = new double[4];
                for (int i = 0; i < 4; i++)
                {
                    dCeof[i] = s[i];
                }
                return dCeof;
            }
            else
            {
                throw new KellImageProcessException("点数少于2点！");
            }
        }
        /// <summary>
        /// N-Derivative
        /// </summary>
        /// <param name="ceof"></param>
        /// <param name="x"></param>
        /// <param name="js"></param>
        /// <param name="range"></param>
        /// <returns></returns>
        public static double GetDSN(double[] ceof, double x, int js, double range)
        {
            double[] dCeof = ceof;
            while (js > 1)
            {
                dCeof = GetDerivativeFunction(dCeof, x - range, x + range);
                js--;
            }
            return GetDS(dCeof, x);
        }

        /// <summary>
        /// 求导函数，返回导函数的系数数组(三次函数的系数数组，一维且长度为4)
        /// </summary>
        /// <param name="ceof"></param>
        /// <param name="start"></param>
        /// <param name="end"></param>
        /// <returns></returns>
        public static double[] GetDerivativeFunction(double[] ceof, double start, double end)
        {
            int points = 11;
            double[] ys = new double[points];
            double x = 0;
            double step = (end - start) / 10;
            for (int i = 0; i < points; i++)
            {
                x = start + step * i;
                ys[i] = GetDS(ceof, x);
            }
            double[] s = new double[5];
            double y = CSharpAlgorithm.Algorithm.Interpolation.GetValueAkima(points, start, step, ys, start, s, -1);
            double[] dCeof = new double[4];
            for (int i = 0; i < 4; i++)
            {
                dCeof[i] = s[i];
            }
            return dCeof;
        }

        /// <summary>
        /// 对点集进行排序
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="rect"></param>
        /// <param name="near">临近系数，取2或3时效果最佳</param>
        /// <param name="firP">指定的路径开始点</param>
        /// <param name="moreWay">如果为true，则返回所有连通路径中的最长者，否则返回一条合并所有临近点的总路径</param>
        /// <param name="linkWays">所有连通路径的点集数组的列表</param>
        /// <returns></returns>
        public static List<Point> SortLinkPs(List<Point> ps, Rectangle rect, double near, Point firP, bool moreWay, out List<List<Point>> linkWays)
        {
            List<Point> retPs = new List<Point>();
            linkWays = new List<List<Point>>();
            if (ps != null && ps.Count > 0)
            {
                Point firstP = firP;
                if (firP == Point.Empty)
                    firstP = FindFirstP(ps, rect);
                retPs.Add(firstP);
                Point current = firstP;
                List<Point> pss = RemoveListElement(ps, firstP);
                while (pss.Count > 0)
                {
                    List<Point> moreW = new List<Point>();
                    foreach (Point p in pss)
                    {
                        if (ThreadCommon.LinkedPoints(current, p, near))
                        {
                            moreW.Add(p);
                            retPs.Add(p);
                            current = p;
                            pss = RemoveListElement(pss, p);
                        }
                    }
                    linkWays.Add(moreW);
                    if (pss.Count == 0)
                        break;
                    Point newP;
                    if (moreWay)
                        newP = FindFirstP(pss, rect);
                    else
                        newP = FindNearP(current, pss);
                    moreW.Add(newP);
                    retPs.Add(newP);
                    current = newP;
                    pss = RemoveListElement(pss, newP);
                }
            }
            if (moreWay)
            {
                int max = 0;
                List<Point> lenMax = new List<Point>();
                foreach (List<Point> pps in linkWays)
                {
                    if (pps.Count > max)
                    {
                        lenMax = pps;
                        max = pps.Count;
                    }
                }
                return lenMax;
            }
            else
            {
                return retPs;
            }
        }
        /// <summary>
        /// 寻找点集路径的开始点
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="rect"></param>
        /// <returns></returns>
        public static Point FindFirstP(List<Point> ps, Rectangle rect)
        {
            if (ps != null && ps.Count > 0)
            {
                Point firstP = ps[0];
                foreach (Point p in ps)
                {
                    if (p.X == rect.Left)
                    {
                        firstP = p;
                        //MessageBox.Show("Left");
                        break;
                    }
                    if (p.Y == rect.Top)
                    {
                        firstP = p;
                        //MessageBox.Show("Top");
                        break;
                    }
                    if (p.X == rect.Right)
                    {
                        firstP = p;
                        //MessageBox.Show("Right");
                        break;
                    }
                    if (p.Y == rect.Bottom)
                    {
                        firstP = p;
                        //MessageBox.Show("Bottom");
                        break;
                    }
                }
                return firstP;
            }
            throw new KellImageProcessException("点集为空，无法求开始点！");
        }
        /// <summary>
        /// 寻找点集路径的最近点
        /// </summary>
        /// <param name="p"></param>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static Point FindNearP(Point p, List<Point> ps)
        {
            if (ps != null && ps.Count > 0)
            {
                Point nearP = ps[0];
                double minDis = Math.Sqrt((p.X - nearP.X) * (p.X - nearP.X) + (p.Y - nearP.Y) * (p.Y - nearP.Y));
                for (int i = 1; i < ps.Count; i++)
                {
                    double tmp = Math.Sqrt((p.X - ps[i].X) * (p.X - ps[i].X) + (p.Y - ps[i].Y) * (p.Y - ps[i].Y));
                    if (minDis > tmp)
                    {
                        minDis = tmp;
                        nearP = ps[i];
                    }
                }
                return nearP;
            }
            throw new KellImageProcessException("点集为空，无法求最近点！");
        }

        private static List<Point> RemoveListElement(List<Point> ps, Point p)
        {
            List<Point> pss = new List<Point>();
            foreach (Point pe in ps)
                pss.Add(pe);
            pss.Remove(p);
            return pss;
        }

        /// <summary>
        /// 根据指定区域获取转弯点
        /// </summary>
        /// <param name="b"></param>
        /// <param name="rect"></param>
        /// <param name="threshold"></param>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static Point GetCornerPoint(Bitmap b, Rectangle rect, byte threshold, out List<Point> ps)
        {
            //threshold = OTSUThresholdValue(FastClipBitmap(b, rect));
            ps = GetEdgeByRect(b, rect, threshold);
            List<List<Point>> lps;
            ps = SortLinkPs(ps, new Rectangle(rect.X + 1, rect.Y + 1, rect.Width - 1, rect.Height - 1), 3, Point.Empty, false, out lps);
            Algorithm sf = new Algorithm();
            Point p;
            if (sf.GetCornerPoint(ps, out p))
                return p;
            else
                return Point.Empty;
        }
        /// <summary>
        /// 根据指定路径获取转弯点
        /// </summary>
        /// <param name="b"></param>
        /// <param name="path"></param>
        /// <param name="threshold"></param>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static Point GetCornerPoint(Bitmap b, GraphicsPath path, byte threshold, out List<Point> ps)
        {
            //threshold = OTSUThresholdValue(FastClipBitmap(b, path));
            ps = GetEdgeByPath(b, path, threshold);
            Rectangle rect = Rectangle.Intersect(new Rectangle(0, 0, b.Width, b.Height), Rectangle.Truncate(path.GetBounds()));
            List<List<Point>> lps;
            ps = SortLinkPs(ps, new Rectangle(rect.X + 1, rect.Y + 1, rect.Width - 1, rect.Height - 1), 3, Point.Empty, false, out lps);
            Algorithm sf = new Algorithm();
            Point p;
            if (sf.GetCornerPoint(ps, out p))
                return p;
            else
                return Point.Empty;
        }
        /// <summary>
        /// 获取三维点集所确定的立方体的对顶线的长度
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static double GetMaxDistanceOf3DpointsCube(List<Point3> ps)
        {
            double dis = 0;
            if (ps != null && ps.Count > 0)
            {
                float minX = float.MaxValue, minY = float.MaxValue, minZ = float.MaxValue;
                float maxX = float.MinValue, maxY = float.MinValue, maxZ = float.MinValue;
                foreach (Point3 p in ps)
                {
                    if (p.X < minX)
                    {
                        minX = p.X;
                    }
                    if (p.Y < minY)
                    {
                        minY = p.Y;
                    }
                    if (p.Z < minZ)
                    {
                        minZ = p.Z;
                    }
                    if (p.X > maxX)
                    {
                        maxX = p.X;
                    }
                    if (p.Y > maxY)
                    {
                        maxY = p.Y;
                    }
                    if (p.Z > maxZ)
                    {
                        maxZ = p.Z;
                    }
                }
                float dx = maxX - minX;
                float dy = maxY - minY;
                float dz = maxZ - minZ;
                dis = Math.Sqrt(dx * dx + dy * dy + dz * dz);
            }
            return dis;
        }
        /// <summary>
        /// 获取位图中指定路径下的平均颜色
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="path"></param>
        public static Color GetAvgColorByPath(Bitmap bmp, GraphicsPath path)
        {
            Color c = Color.Empty;
            PointF Loc = path.GetBounds().Location;
            Bitmap b = FastClipBitmap(bmp, path);
            BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadOnly, b.PixelFormat);
            unsafe
            {
                int sumR = 0;
                int sumG = 0;
                int sumB = 0;
                int count = 0;
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * b.Width;
                for (int y = 0; y < b.Height; y++)
                {
                    for (int x = 0; x < b.Width; x++)
                    {
                        if (path.IsVisible(x + Loc.X, y + Loc.Y))
                        {
                            sumR += p[2];
                            sumG += p[1];
                            sumB += p[0];
                            count++;
                        }
                        p += BPP;
                    }
                    p += offset;
                }
                if (count > 0)
                    c = Color.FromArgb(sumR / count, sumG / count, sumB / count);
            }
            b.UnlockBits(data);
            return c;
        }
        /// <summary>
        /// 获取位图中指定路径下的平均颜色，并返回该路径下颜色的色系范围
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="path"></param>
        /// <param name="needSmooth"></param>
        /// <param name="range"></param>
        public static Color GetAvgColorByPath(Bitmap bmp, GraphicsPath path, bool needSmooth, out int range)
        {
            Color c = Color.Empty;
            PointF Loc = path.GetBounds().Location;
            Bitmap b = FastClipBitmap(bmp, path);

            if (needSmooth)
                b = Smooth(b);

            BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadOnly, b.PixelFormat);
            List<Point3> ps = new List<Point3>();
            unsafe
            {
                int sumR = 0;
                int sumG = 0;
                int sumB = 0;
                int count = 0;
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * b.Width;
                for (int y = 0; y < b.Height; y++)
                {
                    for (int x = 0; x < b.Width; x++)
                    {
                        if (path.IsVisible(x + Loc.X, y + Loc.Y))
                        {
                            sumR += p[2];
                            sumG += p[1];
                            sumB += p[0];
                            ps.Add(new Point3(p[2], p[1], p[0]));
                            count++;
                        }
                        p += BPP;
                    }
                    p += offset;
                }
                if (count > 0)
                    c = Color.FromArgb(sumR / count, sumG / count, sumB / count);
            }
            b.UnlockBits(data);
            range = (int)GetMaxDistanceOf3DpointsCube(ps);
            return c;
        }
        /// <summary>
        /// 获取位图中指定区域下的平均颜色
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="rect"></param>
        public static Color GetAvgColorByRect(Bitmap bmp, Rectangle rect)
        {
            Color c = Color.Empty;
            Bitmap b = FastClipBitmap(bmp, rect);
            BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadOnly, b.PixelFormat);
            unsafe
            {
                int sumR = 0;
                int sumG = 0;
                int sumB = 0;
                int count = b.Width * b.Height;
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * b.Width;
                for (int y = 0; y < b.Height; y++)
                {
                    for (int x = 0; x < b.Width; x++)
                    {
                        sumR += p[2];
                        sumG += p[1];
                        sumB += p[0];
                        p += BPP;
                    }
                    p += offset;
                }
                if (count > 0)
                    c = Color.FromArgb(sumR / count, sumG / count, sumB / count);
            }
            b.UnlockBits(data);
            return c;
        }
        /// <summary>
        /// 获取位图中指定区域下的平均颜色，并返回该区域下颜色的色系范围
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="rect"></param>
        /// <param name="needSmooth"></param>
        /// <param name="range"></param>
        public static Color GetAvgColorByRect(Bitmap bmp, Rectangle rect, bool needSmooth, out int range)
        {
            Color c = Color.Empty;
            Bitmap b = FastClipBitmap(bmp, rect);

            if (needSmooth)
                b = Smooth(b);

            BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadOnly, b.PixelFormat);
            List<Point3> ps = new List<Point3>();
            unsafe
            {
                int sumR = 0;
                int sumG = 0;
                int sumB = 0;
                int count = b.Width * b.Height;
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * b.Width;
                for (int y = 0; y < b.Height; y++)
                {
                    for (int x = 0; x < b.Width; x++)
                    {
                        sumR += p[2];
                        sumG += p[1];
                        sumB += p[0];
                        ps.Add(new Point3(p[2], p[1], p[0]));
                        p += BPP;
                    }
                    p += offset;
                }
                if (count > 0)
                    c = Color.FromArgb(sumR / count, sumG / count, sumB / count);
            }
            b.UnlockBits(data);
            range = (int)GetMaxDistanceOf3DpointsCube(ps);
            return c;
        }
        /// <summary>
        /// 判断位图中指定路径下是否存在超出指定色系范围内的NG像素点，采用Distance方法，较快
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="path"></param>
        /// <param name="needSmooth"></param>
        /// <param name="avgColor"></param>
        /// <param name="tolerance"></param>
        /// <param name="step"></param>
        /// <returns></returns>
        public static bool HaveNGcolorPointByPath(Bitmap bmp, GraphicsPath path, bool needSmooth, Color avgColor, int tolerance, int step)
        {
            if (step < 1)
            {
                throw new KellImageProcessException("步长必须是一个像素以上！");
            }
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            if (BPP < 3)
            {
                throw new KellImageProcessException("必须是24位以上的RGB位图！");
            }
            bool ng = false;
            PointF Loc = path.GetBounds().Location;
            Bitmap b = FastClipBitmap(bmp, path);

            if (needSmooth)
                b = Smooth(b);

            BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadOnly, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * b.Width;
                for (int y = 0; y < b.Height; y += step)
                {
                    for (int x = 0; x < b.Width; x += step)
                    {
                        if (path.IsVisible(x + Loc.X, y + Loc.Y))
                        {
                            double distance = Math.Sqrt((avgColor.R - p[2]) * (avgColor.R - p[2]) + (avgColor.G - p[1]) * (avgColor.G - p[1]) + (avgColor.B - p[0]) * (avgColor.B - p[0]));
                            if (distance > tolerance)
                            {
                                ng = true;
                                break;
                            }
                        }
                        p += BPP;
                    }
                    if (ng)
                        break;
                    p += offset;
                }
            }
            b.UnlockBits(data);
            return ng;
        }
        /// <summary>
        /// 判断位图中指定区域下是否存在超出指定色系范围内的NG像素点，采用Distance方法，较快
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="rect"></param>
        /// <param name="needSmooth"></param>
        /// <param name="avgColor"></param>
        /// <param name="tolerance"></param>
        /// <param name="step"></param>
        /// <returns></returns>
        public static bool HaveNGcolorPointByRect(Bitmap bmp, Rectangle rect, bool needSmooth, Color avgColor, int tolerance, int step)
        {
            if (step < 1)
            {
                throw new KellImageProcessException("步长必须是一个像素以上！");
            }
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            if (BPP < 3)
            {
                throw new KellImageProcessException("必须是24位以上的RGB位图！");
            }
            bool ng = false;
            Bitmap b = FastClipBitmap(bmp, rect);

            if (needSmooth)
                b = Smooth(b);

            BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadOnly, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * b.Width;
                for (int y = 0; y < b.Height; y += step)
                {
                    for (int x = 0; x < b.Width; x += step)
                    {
                        double distance = Math.Sqrt((avgColor.R - p[2]) * (avgColor.R - p[2]) + (avgColor.G - p[1]) * (avgColor.G - p[1]) + (avgColor.B - p[0]) * (avgColor.B - p[0]));
                        if (distance > tolerance)
                        {
                            ng = true;
                            break;
                        }
                        p += BPP;
                    }
                    if (ng)
                        break;
                    p += offset;
                }
            }
            b.UnlockBits(data);
            return ng;
        }
        /// <summary>
        /// 判断位图中指定路径下是否存在超出指定色系范围内的NG像素点，采用Contains方法，较慢
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="path"></param>
        /// <param name="needSmooth"></param>
        /// <param name="cs"></param>
        /// <param name="step"></param>
        /// <returns></returns>
        public static bool HaveNGcolorPointByPath(Bitmap bmp, GraphicsPath path, bool needSmooth, List<int> cs, int step)
        {
            if (step < 1)
            {
                throw new KellImageProcessException("步长必须是一个像素以上！");
            }
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            if (BPP < 3)
            {
                throw new KellImageProcessException("必须是24位以上的RGB位图！");
            }
            bool ng = false;
            PointF Loc = path.GetBounds().Location;
            Bitmap b = FastClipBitmap(bmp, path);

            if (needSmooth)
                b = Smooth(b);

            BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadOnly, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * b.Width;
                for (int y = 0; y < b.Height; y += step)
                {
                    for (int x = 0; x < b.Width; x += step)
                    {
                        if (path.IsVisible(x + Loc.X, y + Loc.Y))
                        {
                            if (!cs.Contains(Color.FromArgb(p[2], p[1], p[0]).ToArgb()))
                            {
                                ng = true;
                                break;
                            }
                        }
                        p += BPP;
                    }
                    if (ng)
                        break;
                    p += offset;
                }
            }
            b.UnlockBits(data);
            return ng;
        }
        /// <summary>
        /// 判断位图中指定区域下是否存在超出指定色系范围内的NG像素点，采用Contains方法，较慢
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="rect"></param>
        /// <param name="needSmooth"></param>
        /// <param name="cs"></param>
        /// <param name="step"></param>
        /// <returns></returns>
        public static bool HaveNGcolorPointByRect(Bitmap bmp, Rectangle rect, bool needSmooth, List<int> cs, int step)
        {
            if (step < 1)
            {
                throw new KellImageProcessException("步长必须是一个像素以上！");
            }
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            if (BPP < 3)
            {
                throw new KellImageProcessException("必须是24位以上的RGB位图！");
            }
            bool ng = false;
            Bitmap b = FastClipBitmap(bmp, rect);

            if (needSmooth)
                b = Smooth(b);

            BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadOnly, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * b.Width;
                for (int y = 0; y < b.Height; y += step)
                {
                    for (int x = 0; x < b.Width; x += step)
                    {
                        if (!cs.Contains(Color.FromArgb(p[2], p[1], p[0]).ToArgb()))
                        {
                            ng = true;
                            break;
                        }
                        p += BPP;
                    }
                    if (ng)
                        break;
                    p += offset;
                }
            }
            b.UnlockBits(data);
            return ng;
        }
        /// <summary>
        /// 获取指定颜色和色差的相似颜色数组(色系)
        /// </summary>
        /// <param name="c"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public static List<int> GetNearestColors(Color c, int tolerance)
        {
            if (tolerance > 443)
            {
                throw new KellImageProcessException("色差超出最大的范围！");
            }

            List<int> cs = new List<int>();

            int minR = c.R - tolerance;
            if (minR < 0) minR = 0;
            int maxR = c.R + tolerance;
            if (maxR > 255) maxR = 255;

            int minG = c.G - tolerance;
            if (minG < 0) minG = 0;
            int maxG = c.G + tolerance;
            if (maxG > 255) maxG = 255;

            int minB = c.B - tolerance;
            if (minB < 0) minB = 0;
            int maxB = c.B + tolerance;
            if (maxB > 255) maxB = 255;

            for (int r = minR; r < maxR + 1; r++)
            {
                for (int g = minG; g < maxG + 1; g++)
                {
                    for (int b = minB; b < maxB + 1; b++)
                    {
                        double distance = Math.Sqrt((c.R - r) * (c.R - r) + (c.G - g) * (c.G - g) + (c.B - b) * (c.B - b));
                        if (distance <= tolerance)
                        {
                            cs.Add(Color.FromArgb(r, g, b).ToArgb());
                        }
                    }
                }
            }
            return cs;
        }
        /// <summary>
        /// 获取在位图中按照指定颜色和色差所确定颜色系内的所有像素点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="c"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public static List<Point> GetPointsByNearestColor(Bitmap bmp, Color c, int tolerance)
        {
            if (tolerance > 443)
            {
                throw new KellImageProcessException("色差超出最大的范围！");
            }
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            if (BPP < 3)
            {
                throw new KellImageProcessException("要求24位以上的位图！");
            }
            List<Point> ps = new List<Point>();
            BitmapData data = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * bmp.Width;
                for (int y = 0; y < bmp.Height; y++)
                {
                    for (int x = 0; x < bmp.Width; x++)
                    {
                        double distance = Math.Sqrt((c.R - p[2]) * (c.R - p[2]) + (c.G - p[1]) * (c.G - p[1]) + (c.B - p[0]) * (c.B - p[0]));
                        if (distance <= tolerance)
                        {
                            ps.Add(new Point(x, y));
                        }
                        p += BPP;
                    }
                    p += offset;
                }
            }
            bmp.UnlockBits(data);

            return ps;
        }
        /// <summary>
        /// 获取在位图中按照指定颜色和色差所确定颜色系内的所有像素点，并返回只含该颜色区域的位图
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="c"></param>
        /// <param name="tolerance"></param>
        /// <param name="backgroundColor"></param>
        /// <param name="outBmp"></param>
        /// <returns></returns>
        public static List<Point> GetPointsByNearestColorByContains(Bitmap bmp, Color c, int tolerance, Color backgroundColor, out Bitmap outBmp)
        {
            outBmp = null;
            if (tolerance > 443)
            {
                throw new KellImageProcessException("色差超出最大的范围！");
            }
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            if (BPP < 3)
            {
                throw new KellImageProcessException("要求24位以上的位图！");
            }
            List<Point> ps = new List<Point>();
            outBmp = new Bitmap(bmp.Width, bmp.Height, bmp.PixelFormat);
            BitmapData data = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
            BitmapData data1 = outBmp.LockBits(new Rectangle(0, 0, outBmp.Width, outBmp.Height), ImageLockMode.WriteOnly, outBmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                byte* p1 = (byte*)data1.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * bmp.Width;
                for (int y = 0; y < bmp.Height; y++)
                {
                    for (int x = 0; x < bmp.Width; x++)
                    {
                        double distance = Math.Sqrt((c.R - p[2]) * (c.R - p[2]) + (c.G - p[1]) * (c.G - p[1]) + (c.B - p[0]) * (c.B - p[0]));
                        if (distance <= tolerance)
                        {
                            ps.Add(new Point(x, y));
                            p1[2] = p[2];
                            p1[1] = p[1];
                            p1[0] = p[0];
                        }
                        else
                        {
                            p1[2] = backgroundColor.R;
                            p1[1] = backgroundColor.G;
                            p1[0] = backgroundColor.B;
                        }
                        p += BPP;
                        p1 += BPP;
                    }
                    p += offset;
                    p1 += offset;
                }
            }
            bmp.UnlockBits(data);
            outBmp.UnlockBits(data1);

            return ps;
        }
        /// <summary>
        /// 获取在位图中按照指定颜色和色差所确定颜色系内的所有像素点，并返回去除该颜色区域的位图
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="c"></param>
        /// <param name="tolerance"></param>
        /// <param name="insteadColor">
        /// <param name="outBmp"></param></param>
        /// <returns></returns>
        public static List<Point> GetPointsByNearestColorByRemove(Bitmap bmp, Color c, int tolerance, Color insteadColor, out Bitmap outBmp)
        {
            outBmp = null;
            if (tolerance > 443)
            {
                throw new KellImageProcessException("色差超出最大的范围！");
            }
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            if (BPP < 3)
            {
                throw new KellImageProcessException("要求24位以上的位图！");
            }
            List<Point> ps = new List<Point>();
            outBmp = new Bitmap(bmp.Width, bmp.Height, bmp.PixelFormat);
            BitmapData data = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
            BitmapData data1 = outBmp.LockBits(new Rectangle(0, 0, outBmp.Width, outBmp.Height), ImageLockMode.WriteOnly, outBmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                byte* p1 = (byte*)data1.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * bmp.Width;
                for (int y = 0; y < bmp.Height; y++)
                {
                    for (int x = 0; x < bmp.Width; x++)
                    {
                        double distance = Math.Sqrt((c.R - p[2]) * (c.R - p[2]) + (c.G - p[1]) * (c.G - p[1]) + (c.B - p[0]) * (c.B - p[0]));
                        if (distance <= tolerance)
                        {
                            ps.Add(new Point(x, y));
                            p1[2] = insteadColor.R;
                            p1[1] = insteadColor.G;
                            p1[0] = insteadColor.B;
                        }
                        else
                        {
                            p1[2] = p[2];
                            p1[1] = p[1];
                            p1[0] = p[0];
                        }
                        p += BPP;
                        p1 += BPP;
                    }
                    p += offset;
                    p1 += offset;
                }
            }
            bmp.UnlockBits(data);
            outBmp.UnlockBits(data1);

            return ps;
        }
        /// <summary>
        /// 获取两个位图的差异点集，并返回差异位图(保留原图)
        /// </summary>
        /// <param name="b1"></param>
        /// <param name="b2"></param>
        /// <param name="tolerance"></param>
        /// <param name="fillColor"></param>
        /// <param name="diffBmp"></param>
        /// <returns></returns>
        public static List<Point> GetDifferentFrom2BmpWithOriginImage(Bitmap b1, Bitmap b2, int tolerance, Color fillColor, out Bitmap diffBmp)
        {
            diffBmp = null;
            if (tolerance > 443)
            {
                throw new KellImageProcessException("色差超出最大的范围！");
            }
            int BPP1 = Image.GetPixelFormatSize(b1.PixelFormat) / 8;
            int BPP2 = Image.GetPixelFormatSize(b2.PixelFormat) / 8;
            int BPP = Math.Min(BPP1, BPP2);
            PixelFormat pf;
            if (BPP1 < BPP2)
                pf = b1.PixelFormat;
            else
                pf = b2.PixelFormat;
            if (BPP < 3)
            {
                throw new KellImageProcessException("要求24位以上的位图！");
            }
            Bitmap bb1, bb2;
            Rectangle rec1 = new Rectangle(0, 0, b1.Width, b1.Height);
            Rectangle rec2 = new Rectangle(0, 0, b2.Width, b2.Height);
            if (BPP1 != BPP2)
            {
                if (BPP1 > BPP2)
                {
                    bb1 = b1.Clone(rec1, pf);
                    bb2 = (Bitmap)b2.Clone();
                }
                else
                {
                    bb1 = (Bitmap)b1.Clone();
                    bb2 = b2.Clone(rec2, pf);
                }
            }
            else
            {
                bb1 = (Bitmap)b1.Clone();
                bb2 = (Bitmap)b2.Clone();
            }
            Rectangle rect = Rectangle.Union(rec1, rec2);
            List<Point> ps = new List<Point>();
            diffBmp = new Bitmap(rect.Width, rect.Height, pf);
            BitmapData d1 = bb1.LockBits(rec1, ImageLockMode.ReadOnly, bb1.PixelFormat);
            BitmapData d2 = bb2.LockBits(rec2, ImageLockMode.ReadOnly, bb2.PixelFormat);
            BitmapData data = diffBmp.LockBits(new Rectangle(0, 0, diffBmp.Width, diffBmp.Height), ImageLockMode.WriteOnly, diffBmp.PixelFormat);
            unsafe
            {
                byte* p1 = (byte*)d1.Scan0;
                byte* p2 = (byte*)d2.Scan0;
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * diffBmp.Width;
                for (int y = 0; y < rect.Height; y++)
                {
                    for (int x = 0; x < rect.Width; x++)
                    {
                        Point pp = new Point(x, y);
                        if (rec1.Contains(pp) && rec2.Contains(pp))
                        {
                            p1 = (byte*)d1.Scan0 + d1.Stride * y + BPP * x;
                            p2 = (byte*)d2.Scan0 + d2.Stride * y + BPP * x;
                            double distance = Math.Sqrt((p1[2] - p2[2]) * (p1[2] - p2[2]) + (p1[1] - p2[1]) * (p1[1] - p2[1]) + (p1[0] - p2[0]) * (p1[0] - p2[0]));
                            if (distance <= tolerance)
                            {
                                p[2] = p1[2];
                                p[1] = p1[1];
                                p[0] = p1[0];
                            }
                            else
                            {
                                ps.Add(new Point(x, y));
                                p[2] = fillColor.R;
                                p[1] = fillColor.G;
                                p[0] = fillColor.B;
                            }
                        }
                        else
                        {
                            ps.Add(new Point(x, y));
                            p[2] = fillColor.R;
                            p[1] = fillColor.G;
                            p[0] = fillColor.B;
                        }
                        if (BPP > 3)
                            p[3] = 255;
                        p += BPP;
                    }
                    p += offset;
                }
            }
            bb1.UnlockBits(d1);
            bb2.UnlockBits(d2);
            diffBmp.UnlockBits(data);

            return ps;
        }
        /// <summary>
        /// 获取两个位图的差异点集，并返回差异二值位图(不保留原图)
        /// </summary>
        /// <param name="b1"></param>
        /// <param name="b2"></param>
        /// <param name="tolerance"></param>
        /// <param name="diffBmp"></param>
        /// <returns></returns>
        public static List<Point> GetDifferentFrom2Bmp(Bitmap b1, Bitmap b2, int tolerance, out Bitmap diffBmp)
        {
            diffBmp = null;
            if (tolerance > 443)
            {
                throw new KellImageProcessException("色差超出最大的范围！");
            }
            int BPP1 = Image.GetPixelFormatSize(b1.PixelFormat) / 8;
            int BPP2 = Image.GetPixelFormatSize(b2.PixelFormat) / 8;
            int BPP = Math.Min(BPP1, BPP2);
            PixelFormat pf;
            if (BPP1 < BPP2)
                pf = b1.PixelFormat;
            else
                pf = b2.PixelFormat;
            if (BPP < 3)
            {
                throw new KellImageProcessException("要求24位以上的位图！");
            }
            Bitmap bb1, bb2;
            Rectangle rec1 = new Rectangle(0, 0, b1.Width, b1.Height);
            Rectangle rec2 = new Rectangle(0, 0, b2.Width, b2.Height);
            if (BPP1 != BPP2)
            {
                if (BPP1 > BPP2)
                {
                    bb1 = b1.Clone(rec1, pf);
                    bb2 = (Bitmap)b2.Clone();
                }
                else
                {
                    bb1 = (Bitmap)b1.Clone();
                    bb2 = b2.Clone(rec2, pf);
                }
            }
            else
            {
                bb1 = (Bitmap)b1.Clone();
                bb2 = (Bitmap)b2.Clone();
            }
            Rectangle rect = Rectangle.Union(rec1, rec2);
            List<Point> ps = new List<Point>();
            diffBmp = new Bitmap(rect.Width, rect.Height, pf);
            BitmapData d1 = bb1.LockBits(rec1, ImageLockMode.ReadOnly, bb1.PixelFormat);
            BitmapData d2 = bb2.LockBits(rec2, ImageLockMode.ReadOnly, bb2.PixelFormat);
            BitmapData data = diffBmp.LockBits(new Rectangle(0, 0, diffBmp.Width, diffBmp.Height), ImageLockMode.WriteOnly, diffBmp.PixelFormat);
            unsafe
            {
                byte* p1 = (byte*)d1.Scan0;
                byte* p2 = (byte*)d2.Scan0;
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * diffBmp.Width;
                for (int y = 0; y < rect.Height; y++)
                {
                    for (int x = 0; x < rect.Width; x++)
                    {
                        Point pp = new Point(x, y);
                        if (rec1.Contains(pp) && rec2.Contains(pp))
                        {
                            p1 = (byte*)d1.Scan0 + d1.Stride * y + BPP * x;
                            p2 = (byte*)d2.Scan0 + d2.Stride * y + BPP * x;
                            double distance = Math.Sqrt((p1[2] - p2[2]) * (p1[2] - p2[2]) + (p1[1] - p2[1]) * (p1[1] - p2[1]) + (p1[0] - p2[0]) * (p1[0] - p2[0]));
                            if (distance <= tolerance)
                            {
                                p[2] = 0;
                                p[1] = 0;
                                p[0] = 0;
                            }
                            else
                            {
                                ps.Add(new Point(x, y));
                                p[2] = 255;
                                p[1] = 255;
                                p[0] = 255;
                            }
                        }
                        else
                        {
                            ps.Add(new Point(x, y));
                            p[2] = 255;
                            p[1] = 255;
                            p[0] = 255;
                        }
                        if (BPP > 3)
                            p[3] = 255;
                        p += BPP;
                    }
                    p += offset;
                }
            }
            bb1.UnlockBits(d1);
            bb2.UnlockBits(d2);
            diffBmp.UnlockBits(data);

            return ps;
        }
        /// <summary>
        /// 边缘增强(半阀值化)
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="threshold">阈值[0,255]</param>
        /// <returns></returns>
        public static Bitmap EdgeEnhance(Bitmap b, int threshold)
        {
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            Bitmap dstImage = (Bitmap)b.Clone();

            BitmapData srcData = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, b.PixelFormat);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, b.PixelFormat);

            // 图像实际处理区域
            // 不考虑最左 1 列和最右 1 列
            // 不考虑最上 1 行和最下 1 行
            int rectTop = 1;
            int rectBottom = height - 1;
            int rectLeft = 1;
            int rectRight = width - 1;

            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;

                int stride = srcData.Stride;
                int offset = stride - width * BPP;

                int pixel = 0;
                int maxPixel = 0;


                // 指向第 1 行
                src += stride;
                dst += stride;
                for (int y = rectTop; y < rectBottom; y++)
                {
                    // 指向每行第 1 列像素
                    src += BPP;
                    dst += BPP;

                    for (int x = rectLeft; x < rectRight; x++)
                    {
                        // Alpha
                        if (BPP > 3)
                            dst[3] = src[3];

                        // 处理 B, G, R 三分量
                        for (int i = 0; i < 3; i++)
                        {
                            // 右上-左下
                            maxPixel = src[i - stride + BPP] - src[i + stride - BPP];
                            if (maxPixel < 0) maxPixel = -maxPixel;

                            // 左上-右下
                            pixel = src[i - stride - BPP] - src[i + stride + BPP];
                            if (pixel < 0) pixel = -pixel;
                            if (pixel > maxPixel) maxPixel = pixel;

                            // 上-下
                            pixel = src[i - stride] - src[i + stride];
                            if (pixel < 0) pixel = -pixel;
                            if (pixel > maxPixel) maxPixel = pixel;

                            // 左-右
                            pixel = src[i - BPP] - src[i + BPP];
                            if (pixel < 0) pixel = -pixel;
                            if (pixel > maxPixel) maxPixel = pixel;

                            // 进行阈值判断
                            if (maxPixel < threshold) maxPixel = 0;

                            dst[i] = (byte)maxPixel;
                        }

                        // 向后移一像素
                        src += BPP;
                        dst += BPP;
                    } // x

                    // 移向下一行
                    // 这里得注意要多移 1 列，因最右边还有 1 列不必处理
                    src += offset + BPP;
                    dst += offset + BPP;
                } // y
            }

            b.UnlockBits(srcData);
            dstImage.UnlockBits(dstData);

            return dstImage;
        } // end of EdgeEnhance

        /// <summary>
        /// 边缘增强
        /// </summary>
        /// <param name="b">位图流</param>
        /// <returns></returns>
        public static Bitmap EdgeEnhance(Bitmap b)
        {
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            Bitmap dstImage = (Bitmap)b.Clone();

            BitmapData srcData = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, b.PixelFormat);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, b.PixelFormat);

            // 图像实际处理区域
            // 不考虑最左 1 列和最右 1 列
            // 不考虑最上 1 行和最下 1 行
            int rectTop = 1;
            int rectBottom = height - 1;
            int rectLeft = 1;
            int rectRight = width - 1;

            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;

                int stride = srcData.Stride;
                int offset = stride - width * BPP;

                int pixel = 0;
                int maxPixel = 0;


                // 指向第 1 行
                src += stride;
                dst += stride;
                for (int y = rectTop; y < rectBottom; y++)
                {
                    // 指向每行第 1 列像素
                    src += BPP;
                    dst += BPP;

                    for (int x = rectLeft; x < rectRight; x++)
                    {
                        // Alpha
                        if (BPP > 3)
                            dst[3] = src[3];

                        // 处理 B, G, R 三分量
                        for (int i = 0; i < 3; i++)
                        {
                            // 右上-左下
                            maxPixel = src[i - stride + BPP] - src[i + stride - BPP];
                            if (maxPixel < 0) maxPixel = -maxPixel;

                            // 左上-右下
                            pixel = src[i - stride - BPP] - src[i + stride + BPP];
                            if (pixel < 0) pixel = -pixel;
                            if (pixel > maxPixel) maxPixel = pixel;

                            // 上-下
                            pixel = src[i - stride] - src[i + stride];
                            if (pixel < 0) pixel = -pixel;
                            if (pixel > maxPixel) maxPixel = pixel;

                            // 左-右
                            pixel = src[i - BPP] - src[i + BPP];
                            if (pixel < 0) pixel = -pixel;
                            if (pixel > maxPixel) maxPixel = pixel;

                            dst[i] = (byte)maxPixel;
                        }

                        // 向后移一像素
                        src += BPP;
                        dst += BPP;
                    } // x

                    // 移向下一行
                    // 这里得注意要多移 1 列，因最右边还有 1 列不必处理
                    src += offset + BPP;
                    dst += offset + BPP;
                } // y
            }

            b.UnlockBits(srcData);
            dstImage.UnlockBits(dstData);

            return dstImage;
        } // end of EdgeEnhance
        /// <summary>
        /// 自动获取指定位图的前背景灰度阀值
        /// </summary>
        /// <param name="bmp"></param>
        /// <returns></returns>
        private static byte GetThresholdAuto(Bitmap bmp)
        {
            byte threshold = OTSUThresholdValue(bmp);

            return threshold;
        }
        /// <summary>
        /// OTSU大津法求阈值
        /// </summary>
        /// <param name="bitmap"></param>
        /// <returns></returns>
        public static byte OTSUThresholdValue(Bitmap bitmap)
        {
            int BPP = Image.GetPixelFormatSize(bitmap.PixelFormat) / 8;
            int wight = bitmap.Width;
            int height = bitmap.Height;
            Rectangle rectangle = new Rectangle(0, 0, wight, height);
            BitmapData bitmapdata = bitmap.LockBits(rectangle, ImageLockMode.ReadOnly, bitmap.PixelFormat);
            int stride = bitmapdata.Stride;
            int offset = stride - rectangle.Width * BPP;
            double[] graybyte = new double[256];//保存灰度值的大小

            #region 统计灰度值
            unsafe
            {
                byte* p = (byte*)bitmapdata.Scan0;
                for (int i = 0; i < rectangle.Height; i++)
                {
                    for (int j = 0; j < rectangle.Width; j++)
                    {
                        byte t = (byte)((19661 * p[2] + 38666 * p[1] + 7209 * p[0]) >> 16);
                        graybyte[t] += 1;//对应的灰度数加一
                        p += BPP;
                    }
                    p += offset;
                }
            }
            #endregion

            bitmap.UnlockBits(bitmapdata);

            int u = 0;
            double sumGray = 0.0;
            double sumNoGray = 0;
            for (int i = 0; i < 256; i++)
            {
                sumGray += graybyte[i] * i;//所有的灰度值之和
                sumNoGray += graybyte[i];//直方图的个元素个数
            }
            double sumBigGray = 0, sumSallGray = 0, sumBigNo = 0, sumSallNo = 0, maxvalue = double.MinValue;
            for (int i = 0; i < 256; i++)
            {
                sumSallGray += graybyte[i] * i;
                sumSallNo += graybyte[i];
                sumBigGray = sumGray - sumSallGray;
                sumBigNo = sumNoGray - sumSallNo;
                if (sumBigGray == 0)
                    break;
                double ProbabilitySmall = sumSallGray / sumSallNo;
                double ProbabilityBig = sumBigGray / sumBigNo;
                double temp = sumBigNo * sumSallNo * (ProbabilityBig - ProbabilitySmall) * (ProbabilityBig - ProbabilitySmall);
                if (temp > maxvalue)
                {
                    maxvalue = temp;
                    u = i;
                }
            }
            return (byte)u;
        }
        /// <summary>
        /// 最佳法求阈值
        /// </summary>
        /// <param name="bitmap"></param>
        /// <returns></returns>
        public static byte TheBestThresholdValue(Bitmap bitmap)
        {
            int BPP = Image.GetPixelFormatSize(bitmap.PixelFormat) / 8;
            int wight = bitmap.Width;
            int height = bitmap.Height;
            Rectangle rectangle = new Rectangle(0, 0, wight, height);
            BitmapData bitmapdata = bitmap.LockBits(rectangle, ImageLockMode.ReadOnly, bitmap.PixelFormat);
            int stride = bitmapdata.Stride;
            int offset = stride - rectangle.Width * BPP;
            double[] graybyte = new double[256];//保存灰度值的大小

            #region 统计灰度值
            unsafe
            {

                byte* p = (byte*)bitmapdata.Scan0;
                for (int i = 0; i < rectangle.Height; i++)
                {
                    for (int j = 0; j < rectangle.Width; j++)
                    {
                        byte t = (byte)((19661 * p[2] + 38666 * p[1] + 7209 * p[0]) >> 16);
                        graybyte[t] += 1;//对应的灰度数加一
                        p += BPP;
                    }
                    p += offset;
                }
            }
            #endregion

            bitmap.UnlockBits(bitmapdata);

            byte maxgray = byte.MinValue;
            byte mingray = byte.MaxValue;

            #region 找出最大的灰度值与最小的灰度值
            for (int i = 0; i < 256; i++)
            {
                if (graybyte[i] == 0)
                    continue;
                else
                {
                    if (i > maxgray)
                        maxgray = (byte)i;
                    if (i < mingray)
                        mingray = (byte)i;
                }
            }
            #endregion

            byte averagemaxandmin = (byte)((maxgray + mingray) / 2);

            // 迭代法求阈值 求以平均值为中心的两部的平均值
            byte newaer = 0;
            while (newaer != averagemaxandmin)
            {
                averagemaxandmin = newaer;
                //第一部份
                double sum1 = 0.0;
                double su1 = 0;
                for (int i = mingray; i < averagemaxandmin; i++)
                {
                    sum1 += i * graybyte[i];
                    su1 += graybyte[i];

                }
                byte aver1 = (byte)(sum1 / su1);

                //第二部份
                double sum2 = 0.0;
                double su2 = 0;
                for (int i = averagemaxandmin; i < maxgray; i++)
                {
                    sum2 += i * graybyte[i];
                    su2 += graybyte[i];
                }
                byte aver2 = (byte)(sum2 / su2);
                newaer = (byte)((aver1 + aver2) / 2);
            }

            return (byte)averagemaxandmin;
        }
        /// <summary>
        /// 获取当前位图前背景的灰度阀值
        /// </summary>
        public static byte CurrentThreshold
        {
            get
            {
                return currentThreshold;
            }
        }
        /// <summary>
        /// 访问或处理的位图
        /// </summary>
        public Bitmap Bitmap
        {
            get
            {
                return bmp;
            }
            set
            {
                bmp = (Bitmap)value.Clone();
                width = bmp.Width;
                height = bmp.Height;
            }
        }
        /// <summary>
        /// 直方图各灰度使用量数组，元素值的取值范围[0, 100](百分比%)
        /// </summary>
        public uint[] HistogramValues
        {
            get
            {
                return histogram;
            }
        }

        /// <summary>
        /// 在指定的路径下寻找最值
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="IsY"></param>
        /// <param name="MinMax"></param>
        /// <param name="path"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <returns></returns>
        public static int GetTopPoint(Bitmap bmp, bool IsY, int MinMax, GraphicsPath path, int blkORwht, byte threshold, bool needGrayExt)
        {
            List<Point> ps = GetKeyPs(bmp, path, blkORwht, threshold, needGrayExt, true, true, false);
            if (ps.Count == 0)
                return 0;
            int minMax = 0;
            if (IsY)
            {
                if (MinMax == 0)
                {
                    minMax = ps[0].Y;
                    for (int i = 1; i < ps.Count; i++)
                    {
                        if (minMax > ps[i].Y)
                            minMax = ps[i].Y;
                    }
                }
                else
                {
                    minMax = ps[0].Y;
                    for (int i = 1; i < ps.Count; i++)
                    {
                        if (minMax < ps[i].Y)
                            minMax = ps[i].Y;
                    }
                }
            }
            else
            {
                if (MinMax == 0)
                {
                    minMax = ps[0].X;
                    for (int i = 1; i < ps.Count; i++)
                    {
                        if (minMax > ps[i].X)
                            minMax = ps[i].X;
                    }
                }
                else
                {
                    minMax = ps[0].X;
                    for (int i = 1; i < ps.Count; i++)
                    {
                        if (minMax < ps[i].X)
                            minMax = ps[i].X;
                    }
                }
            }
            return minMax;
        }

        /// <summary>
        /// 将指定的位图转化为亚像素大位图
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="quasi">正整数</param>
        /// <returns></returns>
        public static Bitmap QuasiPixel(Bitmap bmp, uint quasi)
        {
            if (bmp == null)
                return null;
            Bitmap bm = (Bitmap)bmp.Clone();
            Rectangle rect = new Rectangle(0, 0, bm.Width * (int)quasi, bm.Height * (int)quasi);
            Bitmap b = new Bitmap(rect.Width, rect.Height, bm.PixelFormat);
            using (Graphics g = Graphics.FromImage(b))
            {
                Rectangle re = new Rectangle(0, 0, bm.Width, bm.Height);
                g.InterpolationMode = InterpolationMode.HighQualityBicubic;
                g.DrawImage(bm, rect, re, GraphicsUnit.Pixel);
            }
            bm.Dispose();
            return b;
        }

        /// <summary>
        /// 获取原位图中指定像素的亚像素位图
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="quasi">正整数</param>
        /// <returns></returns>
        public static Bitmap QuasiPixel(Bitmap bmp, int x, int y, uint quasi)
        {
            if (bmp == null)
                return null;
            GraphicsUnit gu = GraphicsUnit.Pixel;
            RectangleF r = bmp.GetBounds(ref gu);
            if (!r.Contains(new PointF(x, y)))
                return null;
            Bitmap retBmp = new Bitmap((int)quasi, (int)quasi, bmp.PixelFormat);
            Bitmap b = QuasiPixel(FastClipBitmap(bmp, new Rectangle(x - 4, y - 4, 9, 9)), quasi);
            PixelFormat pf = b.PixelFormat;
            int BPP = Image.GetPixelFormatSize(pf) / 8;
            BitmapData srcData = b.LockBits(new Rectangle(4 * retBmp.Width, 4 * retBmp.Width, retBmp.Width, retBmp.Width), ImageLockMode.ReadOnly, pf);
            BitmapData dstData = retBmp.LockBits(new Rectangle(0, 0, retBmp.Width, retBmp.Width), ImageLockMode.WriteOnly, pf);
            int stride1 = srcData.Stride;
            int offset1 = stride1 - BPP * retBmp.Width;
            int stride2 = dstData.Stride;
            int offset2 = stride2 - BPP * retBmp.Width;
            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;
                for (int yy = 0; yy < quasi; yy++)
                {
                    for (int xx = 0; xx < quasi; xx++)
                    {
                        for (int i = 0; i < BPP; i++)
                        {
                            dst[i] = src[i];
                        } // i
                        src += BPP;
                        dst += BPP;
                    } // xx
                    src += offset1;
                    dst += offset2;
                } // yy
            }
            b.UnlockBits(srcData);
            retBmp.UnlockBits(dstData);
            b.Dispose();
            return retBmp;
        }

        /// <summary>
        /// 以指定放大系数放大位图并且格式化
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="scale"></param>
        /// <returns></returns>
        public static Bitmap ScalePixelByFormat(Bitmap bmp, uint scale)
        {
            if (bmp == null)
                return null;
            Bitmap bm = (Bitmap)bmp.Clone();
            int width = bm.Width * (int)scale;
            int height = bm.Height * (int)scale;
            Bitmap b = new Bitmap(width, height, bm.PixelFormat);
            byte[,] aveR = new byte[bm.Width, bm.Height];
            byte[,] aveG = new byte[bm.Width, bm.Height];
            byte[,] aveB = new byte[bm.Width, bm.Height];

            PixelFormat pf = bm.PixelFormat;
            int BPP = Image.GetPixelFormatSize(pf) / 8;
            BitmapData srcData = bm.LockBits(new Rectangle(0, 0, bm.Width, bm.Height), ImageLockMode.ReadOnly, pf);
            int stride1 = srcData.Stride;
            int offset1 = stride1 - BPP * bm.Width;
            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                for (int yy = 0; yy < bm.Height; yy++)
                {
                    for (int xx = 0; xx < bm.Width; xx++)
                    {
                        aveR[xx, yy] = src[2];
                        aveG[xx, yy] = src[1];
                        aveB[xx, yy] = src[0];
                        src += BPP;
                    } // xx
                    src += offset1;
                } // yy
            }
            bm.UnlockBits(srcData);
            bm.Dispose();
            BitmapData dstData = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, pf);
            int stride2 = dstData.Stride;
            int offset2 = stride2 - BPP * b.Width;
            unsafe
            {
                byte* dst = (byte*)dstData.Scan0;
                for (int yy = 0; yy < height; yy++)
                {
                    for (int xx = 0; xx < width; xx++)
                    {
                        int divx = (int)(xx / scale);
                        int divy = (int)(yy / scale);
                        dst[2] = aveR[divx, divy];
                        dst[1] = aveG[divx, divy];
                        dst[0] = aveB[divx, divy];
                        dst += BPP;
                    } // xx
                    dst += offset2;
                } // yy
            }
            b.UnlockBits(dstData);
            return b;
        }

        /// <summary>
        /// 以指定放大系数放大位图不用格式化
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="scale"></param>
        /// <returns></returns>
        public static Bitmap ScalePixelNoFormat(Bitmap bmp, uint scale)
        {
            if (bmp == null)
                return null;
            Bitmap bm = (Bitmap)bmp.Clone();
            Rectangle rect = new Rectangle(0, 0, bm.Width * (int)scale, bm.Height * (int)scale);
            Bitmap b = new Bitmap(rect.Width, rect.Height, bm.PixelFormat);
            using (Graphics g = Graphics.FromImage(b))
            {
                Rectangle re = new Rectangle(0, 0, bm.Width, bm.Height);
                g.InterpolationMode = InterpolationMode.HighQualityBicubic;
                g.DrawImage(bm, rect, re, GraphicsUnit.Pixel);
            }
            bm.Dispose();
            return b;
        }

        /// <summary>
        /// 快速剪裁指定区域的位图
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="rect"></param>
        /// <returns></returns>
        public static Bitmap FastClipBitmap(Bitmap bmp, Rectangle rect)
        {
            //Bitmap srcBmp = (Bitmap)bmp.Clone();
            PixelFormat pf = bmp.PixelFormat;
            int BPP = Image.GetPixelFormatSize(pf) / 8;
            GraphicsUnit gu = GraphicsUnit.Pixel;
            Rectangle rect1 = Rectangle.Intersect(Rectangle.Truncate(bmp.GetBounds(ref gu)), rect);
            if (rect1 == Rectangle.Empty || rect1.Width == 0 || rect1.Height == 0)
                return null;
            Bitmap dstImage = new Bitmap(rect1.Width, rect1.Height, pf);
            BitmapData srcData = bmp.LockBits(new Rectangle(rect1.X, rect1.Y, rect1.Width, rect1.Height), ImageLockMode.ReadOnly, pf);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, rect1.Width, rect1.Height), ImageLockMode.WriteOnly, pf);
            int stride1 = srcData.Stride;
            int offset1 = stride1 - BPP * rect1.Width;
            int stride2 = dstData.Stride;
            int offset2 = stride2 - BPP * rect1.Width;
            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;
                for (int y = rect1.Y; y < rect1.Height + rect1.Y; y++)
                {
                    for (int x = rect1.X; x < rect1.Width + rect1.X; x++)
                    {
                        for (int i = 0; i < BPP; i++)
                        {
                            dst[i] = src[i];
                        } // i
                        src += BPP;
                        dst += BPP;
                    } // x
                    src += offset1;
                    dst += offset2;
                } // y
            }
            bmp.UnlockBits(srcData);
            dstImage.UnlockBits(dstData);
            return dstImage;
        }
        /// <summary>
        /// 合并景深位图(返回的位图尺寸为最大的那张大小)
        /// Bright比Gray速度快
        /// </summary>
        /// <param name="bmps">位图数组</param>
        /// <returns></returns>
        public static Bitmap MergeJSByBright(List<Bitmap> bmps)
        {
            if (bmps.Count < 2)
            {
                return bmps[0];
            }
            else
            {
                Color[,] maxTrenchants = GetTheMaxTrenchantPixelsByBright(bmps);
                Rectangle rect = new Rectangle(0, 0, maxTrenchants.GetLength(0), maxTrenchants.GetLength(1));
                Bitmap bmp = new Bitmap(rect.Width, rect.Height, bmps[0].PixelFormat);
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                BitmapData data = bmp.LockBits(rect, ImageLockMode.WriteOnly, bmp.PixelFormat);
                unsafe
                {
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride - BPP * rect.Width;
                    for (int y = 0; y < rect.Height; y++)
                    {
                        for (int x = 0; x < rect.Width; x++)
                        {
                            Color color = maxTrenchants[x, y];
                            if (BPP > 3)
                                p[3] = color.A;
                            p[2] = color.R;
                            p[1] = color.G;
                            p[0] = color.B;
                            p += BPP;
                        }  // x
                        p += offset;
                    } // y
                }
                bmp.UnlockBits(data);
                return bmp;
            }
        }

        /// <summary>
        /// 获取位图数组中最清晰的像素颜色数组
        /// </summary>
        /// <param name="bmps">位图数组</param>
        /// <returns></returns>
        private static Color[,] GetTheMaxTrenchantPixelsByBright(List<Bitmap> bmps)
        {
            if (bmps.Count < 1)
                return null;
            Size maxSize = new Size(bmps[0].Width, bmps[0].Height);
            for (int i = 1; i < bmps.Count; i++)
            {
                if (bmps[i].Width > maxSize.Width && bmps[i].Height > maxSize.Height)
                {
                    maxSize = bmps[i].Size;
                }
            }
            float[,] maxTrenchants = new float[maxSize.Width, maxSize.Height];
            Color[,] colors = new Color[maxSize.Width, maxSize.Height];
            byte R, G, B;
            float bright = 0;
            float max = 0, min = 1F;
            unsafe
            {
                Bitmap bmp = bmps[0];
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                Rectangle rect = new Rectangle(0, 0, bmp.Width, bmp.Height);
                BitmapData data = bmp.LockBits(rect, ImageLockMode.ReadOnly, bmp.PixelFormat);
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * rect.Width;
                p += stride;
                for (int y = 1; y < rect.Height - 1; y++)
                {
                    p += BPP;
                    for (int x = 1; x < rect.Width - 1; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];
                        Color color = Color.FromArgb(R, G, B);
                        bright = color.GetBrightness();
                        byte* l = p - BPP;
                        byte* u = p - stride;
                        byte* r = p + BPP;
                        byte* d = p + stride;
                        byte lR = l[2];
                        byte lG = l[1];
                        byte lB = l[0];
                        float lb = Color.FromArgb(lR, lG, lB).GetBrightness();
                        byte uR = u[2];
                        byte uG = u[1];
                        byte uB = u[0];
                        float ub = Color.FromArgb(uR, uG, uB).GetBrightness();
                        byte rR = r[2];
                        byte rG = r[1];
                        byte rB = r[0];
                        float rb = Color.FromArgb(rR, rG, rB).GetBrightness();
                        byte dR = d[2];
                        byte dG = d[1];
                        byte dB = d[0];
                        float db = Color.FromArgb(dR, dG, dB).GetBrightness();
                        float lbAbs = Math.Abs(lb - bright);
                        float ubAbs = Math.Abs(ub - bright);
                        float rbAbs = Math.Abs(rb - bright);
                        float dbAbs = Math.Abs(db - bright);
                        max = GetMax(lbAbs, ubAbs, rbAbs, dbAbs);
                        min = GetMin(lbAbs, ubAbs, rbAbs, dbAbs);
                        maxTrenchants[x, y] = max - min;
                        colors[x, y] = color;
                        p += BPP;
                    }  // x
                    p += BPP + offset;
                } // y
                bmp.UnlockBits(data);
                for (int i = 1; i < bmps.Count; i++)
                {
                    max = 0;
                    min = 1F;
                    Bitmap bm = bmps[i];
                    BPP = Image.GetPixelFormatSize(bm.PixelFormat) / 8;
                    Rectangle rec = new Rectangle(0, 0, bm.Width, bm.Height);
                    BitmapData dat = bm.LockBits(rec, ImageLockMode.ReadOnly, bm.PixelFormat);
                    byte* p1 = (byte*)dat.Scan0;
                    int stride1 = data.Stride;
                    int offset1 = stride1 - BPP * rec.Width;
                    p1 += stride1;
                    for (int y = 1; y < rec.Height - 1; y++)
                    {
                        p1 += BPP;
                        for (int x = 1; x < rec.Width - 1; x++)
                        {
                            R = p1[2];
                            G = p1[1];
                            B = p1[0];
                            Color color = Color.FromArgb(R, G, B);
                            bright = color.GetBrightness();
                            byte* l = p1 - BPP;
                            byte* u = p1 - stride1;
                            byte* r = p1 + BPP;
                            byte* d = p1 + stride1;
                            byte lR = l[2];
                            byte lG = l[1];
                            byte lB = l[0];
                            float lb = Color.FromArgb(lR, lG, lB).GetBrightness();
                            byte uR = u[2];
                            byte uG = u[1];
                            byte uB = u[0];
                            float ub = Color.FromArgb(uR, uG, uB).GetBrightness();
                            byte rR = r[2];
                            byte rG = r[1];
                            byte rB = r[0];
                            float rb = Color.FromArgb(rR, rG, rB).GetBrightness();
                            byte dR = d[2];
                            byte dG = d[1];
                            byte dB = d[0];
                            float db = Color.FromArgb(dR, dG, dB).GetBrightness();
                            float lbAbs = Math.Abs(lb - bright);
                            float ubAbs = Math.Abs(ub - bright);
                            float rbAbs = Math.Abs(rb - bright);
                            float dbAbs = Math.Abs(db - bright);
                            max = GetMax(lbAbs, ubAbs, rbAbs, dbAbs);
                            min = GetMin(lbAbs, ubAbs, rbAbs, dbAbs);
                            if (maxTrenchants[x, y] < max - min)
                            {
                                maxTrenchants[x, y] = max - min;
                                colors[x, y] = color;
                            }
                            p1 += BPP;
                        }  // x
                        p1 += BPP + offset1;
                    } // y
                    bm.UnlockBits(dat);
                }
            }
            return colors;
        }

        /// <summary>
        /// 合并景深位图(返回的位图尺寸为最大的那张大小)
        /// Gray速度没有Bright快
        /// </summary>
        /// <param name="bmps">位图数组</param>
        /// <returns></returns>
        public static Bitmap MergeJSByGray(List<Bitmap> bmps)
        {
            if (bmps.Count == 1)
            {
                return bmps[0];
            }
            else
            {
                Color[,] maxTrenchants = GetTheMaxTrenchantPixelsByGray(bmps);
                Rectangle rect = new Rectangle(0, 0, maxTrenchants.GetLength(0), maxTrenchants.GetLength(1));
                Bitmap bmp = new Bitmap(rect.Width, rect.Height, bmps[0].PixelFormat);
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                BitmapData data = bmp.LockBits(rect, ImageLockMode.WriteOnly, bmp.PixelFormat);
                unsafe
                {
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride - BPP * rect.Width;
                    for (int y = 0; y < rect.Height; y++)
                    {
                        for (int x = 0; x < rect.Width; x++)
                        {
                            Color color = maxTrenchants[x, y];
                            if (BPP > 3)
                                p[3] = color.A;
                            p[2] = color.R;
                            p[1] = color.G;
                            p[0] = color.B;
                            p += BPP;
                        }  // x
                        p += offset;
                    } // y
                }
                bmp.UnlockBits(data);
                return bmp;
            }
        }
        /// <summary>
        /// 获取位图数组中最清晰的像素颜色数组
        /// </summary>
        /// <param name="bmps">位图数组</param>
        /// <returns></returns>
        private static Color[,] GetTheMaxTrenchantPixelsByGray(List<Bitmap> bmps)
        {
            if (bmps.Count < 1)
                return null;
            Size maxSize = new Size(bmps[0].Width, bmps[0].Height);
            for (int i = 1; i < bmps.Count; i++)
            {
                if (bmps[i].Width > maxSize.Width && bmps[i].Height > maxSize.Height)
                {
                    maxSize = bmps[i].Size;
                }
            }
            byte[,] maxTrenchants = new byte[maxSize.Width, maxSize.Height];
            Color[,] colors = new Color[maxSize.Width, maxSize.Height];
            byte R, G, B;
            byte gray = 128;
            byte max = 0, min = 255;
            unsafe
            {
                Bitmap bmp = bmps[0];
                Rectangle rect = new Rectangle(0, 0, bmp.Width, bmp.Height);
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                BitmapData data = bmp.LockBits(rect, ImageLockMode.ReadOnly, bmp.PixelFormat);
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * rect.Width;
                p += stride;
                for (int y = 1; y < rect.Height - 1; y++)
                {
                    p += BPP;
                    for (int x = 1; x < rect.Width - 1; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];
                        gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                        byte* l = p - BPP;
                        byte* u = p - stride;
                        byte* r = p + BPP;
                        byte* d = p + stride;
                        byte lR = l[2];
                        byte lG = l[1];
                        byte lB = l[0];
                        byte lb = (byte)((19661 * lR + 38666 * lG + 7209 * lB) >> 16);
                        byte uR = u[2];
                        byte uG = u[1];
                        byte uB = u[0];
                        byte ub = (byte)((19661 * uR + 38666 * uG + 7209 * uB) >> 16);
                        byte rR = r[2];
                        byte rG = r[1];
                        byte rB = r[0];
                        byte rb = (byte)((19661 * rR + 38666 * rG + 7209 * rB) >> 16);
                        byte dR = d[2];
                        byte dG = d[1];
                        byte dB = d[0];
                        byte db = (byte)((19661 * dR + 38666 * dG + 7209 * dB) >> 16);
                        byte lbAbs = (byte)Math.Abs(lb - gray);
                        byte ubAbs = (byte)Math.Abs(ub - gray);
                        byte rbAbs = (byte)Math.Abs(rb - gray);
                        byte dbAbs = (byte)Math.Abs(db - gray);
                        max = GetMax(lbAbs, ubAbs, rbAbs, dbAbs);
                        min = GetMin(lbAbs, ubAbs, rbAbs, dbAbs);
                        maxTrenchants[x, y] = (byte)(max - min);
                        colors[x, y] = Color.FromArgb(R, G, B);
                        p += BPP;
                    }  // x
                    p += BPP + offset;
                } // y
                bmp.UnlockBits(data);
                for (int i = 1; i < bmps.Count; i++)
                {
                    max = 0;
                    min = 255;
                    Bitmap bm = bmps[i];
                    Rectangle rec = new Rectangle(0, 0, bm.Width, bm.Height);
                    BPP = Image.GetPixelFormatSize(bm.PixelFormat) / 8;
                    BitmapData dat = bm.LockBits(rec, ImageLockMode.ReadOnly, bm.PixelFormat);
                    byte* p1 = (byte*)dat.Scan0;
                    int stride1 = data.Stride;
                    int offset1 = stride1 - BPP * rec.Width;
                    p1 += stride1;
                    for (int y = 1; y < rec.Height - 1; y++)
                    {
                        p1 += BPP;
                        for (int x = 1; x < rec.Width - 1; x++)
                        {
                            R = p1[2];
                            G = p1[1];
                            B = p1[0];
                            gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                            byte* l = p1 - BPP;
                            byte* u = p1 - stride1;
                            byte* r = p1 + BPP;
                            byte* d = p1 + stride1;
                            byte lR = l[2];
                            byte lG = l[1];
                            byte lB = l[0];
                            byte lb = (byte)((19661 * lR + 38666 * lG + 7209 * lB) >> 16);
                            byte uR = u[2];
                            byte uG = u[1];
                            byte uB = u[0];
                            byte ub = (byte)((19661 * uR + 38666 * uG + 7209 * uB) >> 16);
                            byte rR = r[2];
                            byte rG = r[1];
                            byte rB = r[0];
                            byte rb = (byte)((19661 * rR + 38666 * rG + 7209 * rB) >> 16);
                            byte dR = d[2];
                            byte dG = d[1];
                            byte dB = d[0];
                            byte db = (byte)((19661 * dR + 38666 * dG + 7209 * dB) >> 16);
                            byte lbAbs = (byte)Math.Abs(lb - gray);
                            byte ubAbs = (byte)Math.Abs(ub - gray);
                            byte rbAbs = (byte)Math.Abs(rb - gray);
                            byte dbAbs = (byte)Math.Abs(db - gray);
                            max = GetMax(lbAbs, ubAbs, rbAbs, dbAbs);
                            min = GetMin(lbAbs, ubAbs, rbAbs, dbAbs);
                            if (maxTrenchants[x, y] < max - min)
                            {
                                maxTrenchants[x, y] = (byte)(max - min);
                                colors[x, y] = Color.FromArgb(R, G, B);
                            }
                            p1 += BPP;
                        }  // x
                        p1 += BPP + offset1;
                    } // y
                    bm.UnlockBits(dat);
                }
            }
            return colors;
        }

        private static float GetMin(float lbAbs, float ubAbs, float rbAbs, float dbAbs)
        {
            float tmp1 = lbAbs < ubAbs ? lbAbs : ubAbs;
            float tmp2 = rbAbs < dbAbs ? rbAbs : dbAbs;
            return tmp1 < tmp2 ? tmp1 : tmp2;
        }

        private static float GetMax(float lbAbs, float ubAbs, float rbAbs, float dbAbs)
        {
            float tmp1 = lbAbs > ubAbs ? lbAbs : ubAbs;
            float tmp2 = rbAbs > dbAbs ? rbAbs : dbAbs;
            return tmp1 > tmp2 ? tmp1 : tmp2;
        }

        private static byte GetMin(byte lbAbs, byte ubAbs, byte rbAbs, byte dbAbs)
        {
            byte tmp1 = lbAbs < ubAbs ? lbAbs : ubAbs;
            byte tmp2 = rbAbs < dbAbs ? rbAbs : dbAbs;
            return tmp1 < tmp2 ? tmp1 : tmp2;
        }

        private static byte GetMax(byte lbAbs, byte ubAbs, byte rbAbs, byte dbAbs)
        {
            byte tmp1 = lbAbs > ubAbs ? lbAbs : ubAbs;
            byte tmp2 = rbAbs > dbAbs ? rbAbs : dbAbs;
            return tmp1 > tmp2 ? tmp1 : tmp2;
        }

        /// <summary>
        /// 按指定的裁剪路径对图像进行裁剪
        /// </summary>
        /// <param name="b">原始图像</param>
        /// <param name="region">裁剪区域</param>
        /// <returns></returns>
        public static Bitmap Crop(Bitmap b, Region region)
        {
            Graphics g = System.Drawing.Graphics.FromImage(b);

            // 获取区域边界
            RectangleF validRect = region.GetBounds(g);

            int x = (int)validRect.X;
            int y = (int)validRect.Y;
            int width = (int)validRect.Width;
            int height = (int)validRect.Height;

            // 对区域进行平移
            Region validRegion = region;
            validRegion.Translate(-x, -y);

            Bitmap dstImage = new Bitmap(width, height, b.PixelFormat);
            Graphics dstGraphics = System.Drawing.Graphics.FromImage(dstImage);

            // 设置剪辑区域
            dstGraphics.SetClip(validRegion, CombineMode.Replace);

            // 绘图
            dstGraphics.DrawImage(b, new Rectangle(0, 0, width, height),
              validRect, GraphicsUnit.Pixel);

            g.Dispose();
            dstGraphics.Dispose();

            return dstImage;
        } // end of Crop

        /// <summary>
        /// 缩放滤波　类似开运算的方法，利用先缩小再放大的不可逆性来达到钝化的效果
        /// </summary>
        /// <param name="b"></param>
        /// <param name="insteadBoundColor"></param>
        /// <returns></returns>
        public static Bitmap Filtering(Bitmap b, Color insteadBoundColor)
        {
            Bitmap bmp = new Bitmap(b.Width, b.Height, b.PixelFormat);
            Bitmap tmp = new Bitmap(b.Width / 2, b.Height / 2, b.PixelFormat);
            using (Graphics g = Graphics.FromImage(tmp))
            {
                g.InterpolationMode = InterpolationMode.NearestNeighbor;
                g.DrawImage(b, new Rectangle(0, 0, tmp.Width, tmp.Height), new Rectangle(0, 0, b.Width, b.Height), GraphicsUnit.Pixel);
            }
            using (Graphics g = Graphics.FromImage(bmp))
            {
                g.InterpolationMode = InterpolationMode.High;// QualityBicubic;
                g.DrawImage(tmp, new Rectangle(0, 0, bmp.Width, bmp.Height), new Rectangle(0, 0, tmp.Width, tmp.Height), GraphicsUnit.Pixel);
                //以下是去边框
                g.DrawRectangle(new Pen(insteadBoundColor), new Rectangle(0, 0, b.Width, b.Height));
                g.DrawLine(new Pen(insteadBoundColor), new Point(1, b.Height - 1), new Point(b.Width - 1, b.Height - 1));
                g.DrawLine(new Pen(insteadBoundColor), new Point(b.Width - 1, 1), new Point(b.Width - 1, b.Height - 1));
                g.DrawLine(new Pen(insteadBoundColor), new Point(1, b.Height - 2), new Point(b.Width - 1, b.Height - 2));
                g.DrawLine(new Pen(insteadBoundColor), new Point(b.Width - 2, 1), new Point(b.Width - 2, b.Height - 1));
            }
            bmp = GrayExtend(bmp);
            return bmp;
        }

        /// <summary>
        /// 对图像进行锐化处理
        /// </summary>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Bitmap Sharpen(Bitmap b)
        {
            //    0 -1  0
            //   -1  5 -1
            //    0 -1  0 / 1
            Matrix3x3 m = new Matrix3x3();
            m.Init(0);
            m.Center = 5;
            m.TopMid = m.MidLeft = m.MidRight = m.BottomMid = -1;
            m.Scale = 1;
            return m.Convolute(b);
        }

        /// <summary>
        /// 对图像进行加强锐化处理
        /// </summary>
        /// <param name="b">位图流</param>
        /// <returns></returns>
        public static Bitmap SharpenMore(Bitmap b)
        {
            //   -1 -1 -1
            //   -1  9 -1
            //   -1 -1 -1 / 1
            Matrix3x3 m = new Matrix3x3();
            m.Init(-1);
            m.Center = 9;
            return m.Convolute(b);
        } // end of SharpenMore

        /// <summary>
        /// 快速算法的中值滤波，返回滤波后的彩图(未经灰度化或二值化的位图)，缺点：效果不明显！
        /// </summary>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Bitmap Filtering(Bitmap b)
        {
            return FilterNxN(b, 3);
        }

        /// <summary>
        /// 利用迭代开运算的方法进行图像滤波，返回二值位图
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="blkORwht">黑检测还是白检测</param>
        /// <param name="threshold"></param>
        /// <param name="degree">迭代次数，范围[1, 255]，默认为1次</param>
        /// <returns></returns>
        public static Bitmap Filtering(Bitmap b, int blkORwht, byte threshold, byte degree)
        {
            Bitmap bmp = Opening(b, blkORwht, threshold, false);
            int cnt = 1;
            while (cnt < degree)
            {
                bmp = Opening(bmp, blkORwht, threshold, true);
                cnt++;
            }
            if (blkORwht == 0)
            {
                using (Graphics g = Graphics.FromImage(bmp))
                {
                    g.DrawRectangle(Pens.White, new Rectangle(0, 0, bmp.Width - 1, bmp.Height - 1));
                }
            }
            return bmp;
        }
        /// <summary>
        /// N×N 窗口中值滤波
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="N">滤波窗口大小，N 为奇数</param>
        /// <returns></returns>
        public static Bitmap FilterNxN(Bitmap b, int N)
        {
            // 如果 N 为偶数，则变 N 为奇数
            if (N % 2 == 0) N++;

            // N×N 窗口数字序列
            //byte[] sequence = new byte[N * N];

            // 窗口半径
            int radius = N / 2;

            int width = b.Width;
            int height = b.Height;

            //Bitmap srcImage = (Bitmap)b.Clone();
            Bitmap dstImage = (Bitmap)b.Clone();
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData srcData = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, b.PixelFormat);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, b.PixelFormat);

            // 图像实际处理区域
            int rectTop = radius;
            int rectBottom = height - radius;
            int rectLeft = radius;
            int rectRight = width - radius;

            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;

                int stride = srcData.Stride;
                int offset = stride - width * BPP;

                // 移向最顶行，即第 radius 行
                src += stride * rectTop;
                dst += stride * rectTop;
                for (int y = rectTop; y < rectBottom; y++)
                {
                    // 移向最左列，即每行第 radius 列
                    src += BPP * rectLeft;
                    dst += BPP * rectLeft;

                    for (int x = rectLeft; x < rectRight; x++)
                    {
                        // Alpha
                        if (BPP > 3)
                            dst[3] = src[3];

                        // 处理 B, G, R 三分量
                        for (int i = 0; i < 3; i++)
                        {
                            List<byte> sequence = new List<byte>();
                            // 收集 N×N 窗口内数字序列
                            for (int m = -radius; m <= radius; m++)
                            {
                                for (int n = -radius; n <= radius; n++)
                                {
                                    //sequence[(m + radius) * N + (n + radius)] = src[i + m * stride + n * BPP];
                                    sequence.Add(src[i + m * stride + n * BPP]);
                                } // n
                            } // m

                            // 根据用户指定的统计方法计算数字序列
                            dst[i] = MathLogic.SearchMid(sequence);
                        } // i


                        // 向后移一像素
                        src += BPP;
                        dst += BPP;
                    } // x

                    // 移向下一行
                    // 这里得注意要多移 radius 列，因最右边还有 radius 列不必处理
                    src += (offset + BPP * radius);
                    dst += (offset + BPP * radius);
                } // y
            }

            b.UnlockBits(srcData);//GDI+中发生一般性错误 ？！！！已解决，是因为b的访问模式为ReadOnly引起的，该为ReadWrite就OK
            dstImage.UnlockBits(dstData);

            return dstImage;
        } // end of FilterNxN
        /// <summary>
        /// 利用迭代闭运算的方法进行图像融合
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="blkORwht">黑检测还是白检测</param>
        /// <param name="threshold"></param>
        /// <param name="degree">迭代次数，范围[2, 255]，默认为2次</param>
        /// <returns></returns>
        public static Bitmap Thawing(Bitmap b, int blkORwht, byte threshold, byte degree)
        {
            Bitmap bmp = Closing(Closing(b, blkORwht, threshold, false), blkORwht, threshold, true);
            int cnt = 2;
            while (cnt < degree)
            {
                bmp = Closing(bmp, blkORwht, threshold, true);
                cnt++;
            }
            return bmp;
        }

        /// <summary>
        /// 开运算
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="blkORwht">黑检测还是白检测</param>
        /// <param name="threshold"></param>
        /// <param name="graied">已经灰度化</param>
        /// <returns></returns>
        private static Bitmap Opening(Bitmap b, int blkORwht, byte threshold, bool graied)
        {
            // 先腐蚀，后膨胀
            b = ErosionCross(b, blkORwht, threshold);
            b = DilationCross(b, blkORwht, threshold);

            return b;
        } // end of Opening

        /// <summary>
        /// 闭运算
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="blkORwht">黑检测还是白检测</param>
        /// <param name="threshold"></param>
        /// <param name="graied">已经灰度化</param>
        /// <returns></returns>
        private static Bitmap Closing(Bitmap b, int blkORwht, byte threshold, bool graied)
        {
            // 先膨胀，后腐蚀
            b = DilationCross(b, blkORwht, threshold);
            b = ErosionCross(b, blkORwht, threshold);

            return b;
        } // end of Closing


        /// <summary>
        /// 阈值化
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="threshold">阈值</param>
        /// <param name="needGrayExtend">是否需要灰度拉伸</param>
        /// <returns></returns>
        public static Bitmap Thresholding(Bitmap b, byte threshold, bool needGrayExtend)
        {
            byte[,] GrayArray = BinaryArray(b, threshold, needGrayExtend);
            Bitmap dstImage = BinaryImage(GrayArray, b.PixelFormat, Color.Black, Color.White);
            return dstImage;
        } // end of Thresholding 
        /// <summary>
        /// 指定区域阈值化
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="rect">指定区域</param>
        /// <param name="threshold">阈值</param>
        /// <param name="needGrayExtend">是否需要灰度拉伸</param>
        /// <returns></returns>
        public static Bitmap Thresholding(Bitmap b, Rectangle rect, byte threshold, bool needGrayExtend)
        {
            Bitmap bb = AccessPixel.FastClipBitmap((Bitmap)b.Clone(), rect);
            byte[,] GrayArray = BinaryArray(bb, threshold, needGrayExtend);
            Bitmap dstImage = BinaryImage(GrayArray, bb.PixelFormat, Color.Black, Color.White);
            bb.Dispose();
            return dstImage;
        } // end of Thresholding 
        /// <summary>
        /// 十字型膨胀
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public static Bitmap DilationCross(Bitmap b, int blkORwht, byte threshold)
        {
            // 先将原始二值图转化为二维数组
            byte[,] srcGray = Image2Array(b);

            // 进行十字型膨胀处理
            byte[,] dstGray = DilationCross(srcGray, blkORwht, threshold);

            // 转换为灰度图像
            return Array2Image(dstGray, b.PixelFormat);
        } // end of DilationCross

        /// <summary>
        /// 在没有灰度化位图的指定区域的十字型膨胀，注意在调用时记得返回的位图位置在path的左上角(Location = rect.Location)
        /// </summary>
        /// <param name="b">没有灰度化的位图</param>
        /// <param name="rect">指定区域</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public static Bitmap DilationCross(Bitmap b, Rectangle rect, int blkORwht, byte threshold)
        {
            Bitmap bb = AccessPixel.FastClipBitmap((Bitmap)b.Clone(), rect);
            // 先将原始二值图转化为二维数组
            byte[,] srcGray = Image2Array(bb);

            // 进行十字型膨胀处理
            byte[,] dstGray = DilationCross(srcGray, blkORwht, threshold);

            bb.Dispose();
            // 转换为灰度图像
            return Array2Image(dstGray, b.PixelFormat);
        } // end of DilationCross
        /// <summary>
        /// 十字型膨胀
        /// </summary>
        /// <param name="src">二值化数组</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        private static byte[,] DilationCross(byte[,] src, int blkORwht, byte threshold)
        {
            int width = src.GetLength(0);
            int height = src.GetLength(1);

            // 初始化目标数组
            byte[,] dst = null;
            if (blkORwht == 0)
                dst = InitArray(width, height, 255);
            else
                dst = InitArray(width, height, 0);

            // 3*3 结构元素
            // 1 0 1
            // 0 0 0
            // 1 0 1
            // 由于使用 3*3 的结构元素，为防止越界，所以不处理最上、下、左、右四边的像素
            int topRect = 1;
            int bottomRect = height - 1;
            int leftRect = 1;
            int rightRect = width - 1;

            for (int y = topRect; y < bottomRect; y++)
            {
                for (int x = leftRect; x < rightRect; x++)
                {
                    // 假设目标图像中当前点为白色
                    byte c = 255;
                    if (blkORwht != 0)
                        c = 0;
                    // 如果原图像中 (-1,0), (0,0), (1,0), (0,-1), (0, 1) 五个点之一有黑点，
                    // 则将目标图像中的当前点，即(0,0)点赋予黑色
                    for (int i = -1; i <= 1; i++)
                    {
                        for (int j = -1; j <= 1; j++)
                        {
                            // 根据 3*3 结构元素，
                            // 不必判断当前像素左上、右上、左下、右下四角上的四点
                            if ((i + 2) % 2 == 1 && (j + 2) % 2 == 1)
                                continue;
                            if (blkORwht == 0)
                            {
                                if (src[x + i, y + j] <= threshold)
                                {
                                    c = 0;
                                    break;
                                }
                            }
                            else
                            {
                                if (src[x + i, y + j] >= threshold)
                                {
                                    c = 255;
                                    break;
                                }
                            }
                        } // j
                    } // i

                    dst[x, y] = c;
                } // x
            } // y

            return dst;
        } // end of DilationCross
        /// <summary>
        /// 十字型腐蚀
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public static Bitmap ErosionCross(Bitmap b, int blkORwht, byte threshold)
        {
            // 先将原始二值图转化为二维数组
            byte[,] srcGray = Image2Array(b);

            // 进行十字型腐蚀处理
            byte[,] dstGray = ErosionCross(srcGray, blkORwht, threshold);

            // 转换为灰度图像
            return Array2Image(dstGray, b.PixelFormat);
        } // end of ErosionCross

        /// <summary>
        /// 在没有灰度化位图的指定区域的十字型腐蚀，注意在调用时记得返回的位图位置在path的左上角(Location = rect.Location)
        /// </summary>
        /// <param name="b">没有灰度化的位图</param>
        /// <param name="rect">指定区域</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public static Bitmap ErosionCross(Bitmap b, Rectangle rect, int blkORwht, byte threshold)
        {
            Bitmap bb = AccessPixel.FastClipBitmap((Bitmap)b.Clone(), rect);
            // 先将原始二值图转化为二维数组
            byte[,] srcGray = Image2Array(bb);

            // 进行十字型腐蚀处理
            byte[,] dstGray = ErosionCross(srcGray, blkORwht, threshold);

            bb.Dispose();
            // 转换为灰度图像
            return Array2Image(dstGray, b.PixelFormat);
        } // end of ErosionCross
        /// <summary>
        /// 十字型腐蚀
        /// </summary>
        /// <param name="src">二值化数组</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        private static byte[,] ErosionCross(byte[,] src, int blkORwht, byte threshold)
        {
            int width = src.GetLength(0);
            int height = src.GetLength(1);

            // 初始化目标数组
            byte[,] dst = null;
            if (blkORwht == 0)
                dst = InitArray(width, height, 255);
            else
                dst = InitArray(width, height, 0);

            // 3*3 结构元素
            // 1 0 1
            // 0 0 0
            // 1 0 1
            // 由于使用 3*3 的结构元素，为防止越界，所以不处理最上、下、左、右四边的像素
            int topRect = 1;
            int bottomRect = height - 1;
            int leftRect = 1;
            int rightRect = width - 1;

            for (int y = topRect; y < bottomRect; y++)
            {
                for (int x = leftRect; x < rightRect; x++)
                {
                    // 假设目标图像中当前点为黑色
                    byte c = 0;
                    if (blkORwht != 0)
                        c = 255;
                    // 如果原图像中的 (-1,0), (0,0), (1,0), (0,-1), (0, 1) 五个点之一有白点，
                    // 则将目标图像中的当前点，即(0,0)点赋予白色
                    for (int i = -1; i <= 1; i++)
                    {
                        for (int j = -1; j <= 1; j++)
                        {
                            // 根据 3*3 十字型结构元素，
                            // 不必判断当前像素左上、右上、左下、右下四角上的四点
                            if ((i + 2) % 2 == 1 && (j + 2) % 2 == 1)
                                continue;
                            if (blkORwht == 0)
                            {
                                if (src[x + i, y + j] >= threshold)
                                {
                                    c = 255;
                                    break;
                                }
                            }
                            else
                            {
                                if (src[x + i, y + j] <= threshold)
                                {
                                    c = 0;
                                    break;
                                }
                            }
                        } // j
                    } // i

                    dst[x, y] = c;
                } // x
            } // y

            return dst;
        } // end of ErosionCross
        /// <summary>
        /// 细化
        /// </summary>
        /// <param name="b">位图流</param>
        /// <returns></returns>
        public static Bitmap Thinning(Bitmap b)
        {
            // 先将原始二值图转化为二维数组
            byte[,] srcGray = Image2Array(b);

            // 进行细化处理
            byte[,] dstGray = Thinning(srcGray);

            // 转换为灰度图像
            return Array2Image(dstGray, b.PixelFormat);
        } // end of Thinning


        /// <summary>
        /// 粗化
        /// </summary>
        /// <param name="b">位图流</param>
        /// <returns></returns>
        public static Bitmap Thickening(Bitmap b)
        {
            // 先将原始二值图转化为二维数组
            byte[,] srcGray = Image2Array(b);

            // 进行粗化处理
            byte[,] dstGray = Thickening(srcGray);

            // 转换为灰度图像
            return Array2Image(dstGray, b.PixelFormat);
        } // end of Thickening

        /// <summary>
        /// 细化
        /// </summary>
        /// <param name="src">二值化数组</param>
        /// <returns></returns>
        private static byte[,] Thinning(byte[,] src)
        {
            int width = src.GetLength(0);
            int height = src.GetLength(1);

            // 需要细化吗？
            bool needThinning = true;

            byte[,] S = new byte[5, 5];
            int num = 0;

            // 5*5 结构元素
            // 由于使用 5*5 的结构元素，为防止越界，所以不处理外围的 2 行、 2 列像素
            int topRect = 2;
            int bottomRect = height - 2;
            int leftRect = 2;
            int rightRect = width - 2;

            while (needThinning)
            {
                needThinning = false;

                // 初始化目标数组为 255，即白色
                byte[,] dst = InitArray(width, height, 255);

                for (int y = topRect; y < bottomRect; y++)
                {
                    for (int x = leftRect; x < rightRect; x++)
                    {
                        // 如果原图像中当前点为白色，则跳过
                        if (src[x, y] > 127)
                            continue;

                        // 获得当前点相邻的 5*5 区域内像素值，
                        // 注意：白色用 0 代表，黑色用 1 代表
                        for (int i = -2; i <= 2; i++)
                        {
                            for (int j = -2; j <= 2; j++)
                            {
                                if (src[x + j, y + i] > 127)
                                    S[i + 2, j + 2] = 0; // 白色
                                else
                                    S[i + 2, j + 2] = 1; // 黑色
                            } // j
                        } // i


                        // 判断条件 1 是否成立
                        num = S[1, 1] + S[1, 2] + S[1, 3] +
                              S[2, 1] + S[2, 3] +
                              S[3, 1] + S[3, 2] + S[3, 3];
                        if (!(num >= 2 && num <= 6))
                        {
                            dst[x, y] = 0;
                            continue;
                        }


                        // 判断条件 2 是否成立
                        num = 0;

                        if (S[1, 2] == 0 && S[1, 1] == 1) num++;
                        if (S[1, 1] == 0 && S[2, 1] == 1) num++;
                        if (S[2, 1] == 0 && S[3, 1] == 1) num++;
                        if (S[3, 1] == 0 && S[3, 2] == 1) num++;
                        if (S[3, 2] == 0 && S[3, 3] == 1) num++;
                        if (S[3, 3] == 0 && S[2, 3] == 1) num++;
                        if (S[2, 3] == 0 && S[1, 3] == 1) num++;
                        if (S[1, 3] == 0 && S[1, 2] == 1) num++;

                        if (!(num == 1))
                        {
                            dst[x, y] = 0;
                            continue;
                        }


                        // 判断条件 3 是否成立
                        if (!(S[1, 2] * S[2, 1] * S[2, 3] == 0))
                        {
                            num = 0;

                            if (S[0, 2] == 0 && S[0, 1] == 1) num++;
                            if (S[0, 1] == 0 && S[1, 1] == 1) num++;
                            if (S[1, 1] == 0 && S[2, 1] == 1) num++;
                            if (S[2, 1] == 0 && S[2, 2] == 1) num++;
                            if (S[2, 2] == 0 && S[2, 3] == 1) num++;
                            if (S[2, 3] == 0 && S[1, 3] == 1) num++;
                            if (S[1, 3] == 0 && S[0, 3] == 1) num++;
                            if (S[0, 3] == 0 && S[0, 2] == 1) num++;

                            if (num == 1)
                            {
                                dst[x, y] = 0;
                                continue;
                            }
                        }


                        // 判断条件 4 是否成立
                        if (!(S[1, 2] * S[2, 1] * S[3, 2] == 0))
                        {
                            num = 0;

                            if (S[1, 1] == 0 && S[1, 0] == 1) num++;
                            if (S[1, 0] == 0 && S[2, 0] == 1) num++;
                            if (S[2, 0] == 0 && S[3, 0] == 1) num++;
                            if (S[3, 0] == 0 && S[3, 1] == 1) num++;
                            if (S[3, 1] == 0 && S[3, 2] == 1) num++;
                            if (S[3, 2] == 0 && S[2, 2] == 1) num++;
                            if (S[2, 2] == 0 && S[1, 2] == 1) num++;
                            if (S[1, 2] == 0 && S[1, 1] == 1) num++;

                            if (num == 1)
                            {
                                dst[x, y] = 0;
                                continue;
                            }
                        }


                        // 如果条件均满足，则删除该点
                        dst[x, y] = 255;

                        // 继续细化！
                        needThinning = true;
                    } // x
                } // y

                // 复制细化后的图像作为新一轮细化对象
                src = (byte[,])dst.Clone();
            } // while

            return src;
        } // end of Thinning

        /// <summary>
        /// 粗化
        /// </summary>
        /// <param name="src">二值化数组</param>
        /// <returns></returns>
        private static byte[,] Thickening(byte[,] src)
        {
            int width = src.GetLength(0);
            int height = src.GetLength(1);

            // 对已二值化像素各颜色分量进行求补
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    src[x, y] ^= 255;
                } // x
            } // y

            return Thinning(src);
        } // end of Thickening

        /// <summary>
        /// 获取经过灰度拉伸后的位图(灰度图)
        /// </summary>
        /// <param name="srcBmp"></param>
        /// <param name="graied"></param>
        /// <returns></returns>
        public static Bitmap GrayExtend(Bitmap srcBmp, bool graied)
        {
            Bitmap dstBmp = null;
            if (srcBmp != null)
            {
                //先灰度化
                if (!graied)
                    dstBmp = Gray((Bitmap)srcBmp.Clone());
                else
                    dstBmp = (Bitmap)srcBmp.Clone();
                int width = dstBmp.Width;
                int height = dstBmp.Height;
                int BPP = Image.GetPixelFormatSize(dstBmp.PixelFormat) / 8;
                BitmapData data = dstBmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, dstBmp.PixelFormat);
                unsafe
                {
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride - width * BPP;
                    for (int y = 0; y < height; y++)
                    {
                        for (int x = 0; x < width; x++)
                        {
                            //因为已经经过灰度化，所以p[0] = p[1] = p[2] = gray
                            byte gray = p[0];
                            p[0] = p[1] = p[2] = GetGrayExValue(gray);
                            p += BPP;
                        } //  x
                        p += offset;
                    } // y
                }
                dstBmp.UnlockBits(data);
            }
            return dstBmp;
        }
        /// <summary>
        /// 获取经过拉伸后的位图
        /// </summary>
        /// <param name="srcBmp"></param>
        /// <returns></returns>
        public static Bitmap GrayExtend(Bitmap srcBmp)
        {
            Bitmap dstBmp = null;
            if (srcBmp != null)
            {
                dstBmp = (Bitmap)srcBmp.Clone();
                int width = dstBmp.Width;
                int height = dstBmp.Height;
                int BPP = Image.GetPixelFormatSize(dstBmp.PixelFormat) / 8;
                BitmapData data = dstBmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, dstBmp.PixelFormat);
                unsafe
                {
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride - width * BPP;
                    for (int y = 0; y < height; y++)
                    {
                        for (int x = 0; x < width; x++)
                        {
                            p[0] = GetGrayExValue(p[0]);
                            p[1] = GetGrayExValue(p[1]);
                            p[2] = GetGrayExValue(p[2]);
                            p += BPP;
                        } //  x
                        p += offset;
                    } // y
                }
                dstBmp.UnlockBits(data);
            }
            return dstBmp;
        }
        /// <summary>
        /// 拉伸范围[0,60)[60,180)[180,255)[255] --> [0,30)[30,220)[220,255)[255]
        /// </summary>
        /// <param name="gray"></param>
        /// <returns></returns>
        private static byte GetGrayExValue(byte gray)
        {
            byte ret = 0;
            if (gray < 60)
            {
                ret = (byte)(gray / 2);//最大值:30
            }
            else if (gray >= 60 && gray < 180)
            {
                if (gray >= 60 && gray <= 85)
                    ret = 30;
                else
                    ret = (byte)(2 * (gray - 85) + 30);//最大值:218
            }
            else if (gray >= 180 && gray < 255)
            {
                if (gray >= 180 && gray < 186)
                    ret = 219;
                else
                    ret = (byte)((gray - 186) / 2 + 220);//最大值:254
            }
            else if (gray == 255)
            {
                ret = 255;
            }
            return ret;
        }

        /// <summary>
        /// 将二维数组转换为灰度位图流
        /// </summary>
        /// <param name="GrayArray">灰度数组</param>
        /// <param name="pf"></param>
        /// <returns></returns>
        public static Bitmap Array2Image(byte[,] GrayArray, PixelFormat pf)
        {
            int width = GrayArray.GetLength(0);
            int height = GrayArray.GetLength(1);
            Bitmap b = new Bitmap(width, height, pf);
            int BPP = Image.GetPixelFormatSize(pf) / 8;
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        if (BPP <= 3)
                        {
                            for (int i = 0; i < BPP; i++)
                                p[i] = GrayArray[x, y];
                        }
                        else
                        {
                            for (int i = 0; i < 3; i++)
                                p[i] = GrayArray[x, y];
                            p[3] = 255;
                        }

                        p += BPP;
                    } //  x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            return b;
        } // end of Array2Image
        /// <summary>
        /// 轮廓提取
        /// </summary>
        /// <param name="b">二值图数据数组</param>
        /// <returns></returns>
        private static byte[,] ContourPick(byte[,] b)
        {
            int width = b.GetLength(0);
            int height = b.GetLength(1);

            byte[,] dst = new byte[width, height];

            // 不考虑图像最外围一圈
            int topRect = 1;
            int bottomRect = height - 1;
            int leftRect = 1;
            int rightRect = width - 1;

            for (int y = topRect; y < bottomRect; y++)
            {
                for (int x = leftRect; x < rightRect; x++)
                {
                    if (b[x, y] == 0)
                    {
                        int sum =
                          b[x - 1, y - 1] + b[x, y - 1] + b[x + 1, y - 1] +
                          b[x - 1, y] + b[x + 1, y] +
                          b[x - 1, y + 1] + b[x, y + 1] + b[x + 1, y + 1];

                        // 如果周围八点全为黑点，则置当前点为白点
                        if (sum == 0)
                            dst[x, y] = 255;
                    }
                    else
                    {
                        dst[x, y] = 255;
                    }
                } // x
            } // y

            return dst;
        } // end of ContourPick
        /// <summary>
        /// 轮廓跟踪
        /// </summary>
        /// <param name="Sign">区域标记数组</param>
        /// <returns></returns>
        private static ushort[,] ContourTrace(ushort[,] Sign)
        {
            int width = Sign.GetLength(0);
            int height = Sign.GetLength(1);

            // 面积及区域数
            int[] Area = ImageArea(Sign);
            int areaNumber = Area.Length;

            // 图像边界
            ushort[,] Boundary = new ushort[width, height];

            // 找到起始点或回到起始点？
            bool findStartPoint = false;

            // 扫描到一个边界点？
            bool findPoint;

            // 起始边界点及当前边界点
            Point Start = new Point(0, 0);
            Point Current = new Point(0, 0);

            // 起始点在左上方，扫描方向为逆时针
            Point[] Direct = new Point[8];
            Direct[0].X = -1; Direct[0].Y = 1;    // SW
            Direct[1].X = 0; Direct[1].Y = 1;    // S
            Direct[2].X = 1; Direct[2].Y = 1;    // SE
            Direct[3].X = 1; Direct[3].Y = 0;    // E
            Direct[4].X = 1; Direct[4].Y = -1;   // NE
            Direct[5].X = 0; Direct[5].Y = -1;   // N
            Direct[6].X = -1; Direct[6].Y = -1;   // NW
            Direct[7].X = -1; Direct[7].Y = 0;    // W

            // 先找到最左上方的边界点
            for (int sign = 1; sign < areaNumber; sign++)
            {
                findStartPoint = false;

                for (int y = 0; y < height && !findStartPoint; y++)
                {
                    for (int x = 0; x < width && !findStartPoint; x++)
                    {
                        if (Sign[x, y] == sign)
                        {
                            findStartPoint = true;

                            // 记录下当前边界起始点
                            Start = new Point(x, y);
                            Boundary[x, y] = (ushort)sign;
                        }
                    } // x
                } // y

                // 当前扫描方向
                int direction = 0;

                // 从初始点开始扫描
                Current = Start;

                // 开始跟踪边界
                findStartPoint = false;
                while (!findStartPoint)
                {
                    // 扫描次数
                    int searchNumber = 0;

                    findPoint = false;
                    while (!findPoint)
                    {
                        // 沿扫描方向查看一个像素
                        int x = Current.X + Direct[direction].X;
                        int y = Current.Y + Direct[direction].Y;

                        if ((x >= 0 && x < width && y >= 0 && y < height) && (Sign[x, y] == sign))
                        {
                            // 找到边界点
                            findPoint = true;

                            Current.X += Direct[direction].X;
                            Current.Y += Direct[direction].Y;

                            // 回到起始点，边界扫描完毕
                            if (Current.X == Start.X && Current.Y == Start.Y)
                            {
                                findStartPoint = true;
                                break;
                            }

                            Boundary[Current.X, Current.Y] = (ushort)sign;

                            // 扫描方向顺时针旋转两格
                            direction -= 2;
                            direction += 8; // 避免出现负数
                            direction %= 8;
                        }
                        else
                        {
                            // 扫描方向逆时针旋转一格
                            direction++;
                            direction %= 8;

                            searchNumber++;
                            if (searchNumber > 8)
                            {
                                // 为孤立点
                                findStartPoint = true;
                                findPoint = true;
                            }
                        }

                    } // findPoint
                } // findStartPoint
            } // sign

            return Boundary;
        } // end of ContourTrace
        /// <summary>
        /// 按指定的区域号绘制出对应的区域
        /// </summary>
        /// <param name="b">二值位图流</param>
        /// <param name="graied">已经灰度化</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="Region">区域号</param>
        /// <param name="showContour">指定bool值，是显示轮廓线，否显示区域块</param>
        /// <param name="color">轮廓线或区域块的颜色</param>
        /// <param name="ps">返回轮廓线或区域块内的点集</param>
        /// <returns></returns>
        public static Bitmap ImageRegion(Bitmap b, bool graied, int blkORwht, byte threshold, bool needGrayExt, ushort[] Region, bool showContour, Color color, out Point[] ps)
        {
            Bitmap bb = (Bitmap)b.Clone();
            List<Point> al = new List<Point>();
            // 进行区域标记
            ushort[,] Sign = ImageSign(bb, blkORwht, threshold, needGrayExt);

            // 按轮廓线进行显示
            if (showContour)
                Sign = ContourTrace(Sign);

            int len = Region.Length;

            int width = bb.Width;
            int height = bb.Height;

            int BPP = Image.GetPixelFormatSize(bb.PixelFormat) / 8;
            BitmapData data = bb.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, bb.PixelFormat);

            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;

                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        ushort sign = Sign[x, y];
                        bool showRegion = false;

                        for (int i = 0; i < len; i++)
                        {
                            if (sign == Region[i])
                            {
                                showRegion = true;
                                break;
                            }
                        } // i

                        // 绘制区域
                        if (showRegion)
                        {
                            p[0] = color.B;
                            p[1] = color.G;
                            p[2] = color.R;
                            al.Add(new Point(x, y));
                        }

                        p += BPP;
                    } // x

                    p += offset;
                } // y
            }

            bb.UnlockBits(data);
            ps = new Point[al.Count];
            for (int i = 0; i < al.Count; i++)
            {
                ps[i] = (Point)al[i];
            }
            return bb;
        } // end of ImageRegion
        /// <summary>
        /// 轮廓跟踪，速度较慢，且得到的轮廓信息很单纯
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExtend"></param>
        /// <returns></returns>
        public static Bitmap ContourTrace(Bitmap b, int blkORwht, byte threshold, bool needGrayExtend)
        {
            Bitmap bm = (Bitmap)b.Clone();
            Bitmap bm1 = (Bitmap)b.Clone();

            // 进行区域标记
            ushort[,] Sign = ImageSign(bm, blkORwht, threshold, needGrayExtend);
            bm.Dispose();
            // 轮廓跟踪
            ushort[,] Boundary = ContourTrace(Sign);

            int width = bm1.Width;
            int height = bm1.Height;
            int BPP = Image.GetPixelFormatSize(bm1.PixelFormat) / 8;

            BitmapData data = bm1.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, bm1.PixelFormat);

            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;

                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        if (Boundary[x, y] != 0)
                        {
                            p[0] = p[1] = p[2] = 0;
                        }
                        else
                        {
                            p[0] = p[1] = p[2] = 255;
                        }

                        p += BPP;
                    } // x

                    p += offset;
                } // y
            }

            bm1.UnlockBits(data);

            return bm1;
        } // end of ContourTrace

        /// <summary>
        /// 指定区域中轮廓跟踪，速度较慢，且得到的轮廓信息很单纯
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExtend"></param>
        /// <param name="rect"></param>
        /// <returns></returns>
        public static Bitmap ContourTrace(Bitmap b, int blkORwht, byte threshold, bool needGrayExtend, Rectangle rect)
        {
            Bitmap bm = (Bitmap)b.Clone();
            Bitmap bm1 = (Bitmap)b.Clone();
            // 进行区域标记
            ushort[,] Sign = ImageSign(bm, rect, blkORwht, threshold, needGrayExtend);
            bm.Dispose();
            // 轮廓跟踪
            ushort[,] Boundary = ContourTrace(Sign);

            int width = rect.Width;
            int height = rect.Height;
            int BPP = Image.GetPixelFormatSize(bm1.PixelFormat) / 8;

            BitmapData data = bm1.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, bm1.PixelFormat);

            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;

                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        if (Boundary[x, y] != 0)
                        {
                            p[0] = p[1] = p[2] = 0;
                        }
                        else
                        {
                            p[0] = p[1] = p[2] = 255;
                        }

                        p += BPP;
                    } // x

                    p += offset;
                } // y
            }

            bm1.UnlockBits(data);

            return bm1;
        } // end of ContourTrace

        /// <summary>
        /// 轮廓提取，速度较快，且得到的轮廓信息很详细，且包括外框
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="threshold">阀值</param>
        /// <param name="needGrayExtend">是否需要灰度拉伸</param>
        /// <returns></returns>
        public static Bitmap ContourPick(Bitmap b, byte threshold, bool needGrayExtend)
        {
            // 将原始二值图转化为二维数组
            Bitmap bm = Bitize((Bitmap)b.Clone(), threshold, needGrayExtend);
            byte[,] srcGray = Image2Array(bm);
            bm.Dispose();
            // 轮廓提取
            byte[,] dstGray = ContourPick(srcGray);

            // 转换为灰度图像
            return Array2Image(dstGray, b.PixelFormat);
        } // end of ContourPick

        /// <summary>
        /// 指定区域中轮廓提取，速度较快，且得到的轮廓信息很详细，且包括外框
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="threshold">阀值</param>
        /// <param name="needGrayExtend">是否需要灰度拉伸</param>
        /// <param name="rect"></param>
        /// <returns></returns>
        public static Bitmap ContourPick(Bitmap b, byte threshold, bool needGrayExtend, Rectangle rect)
        {
            // 将原始二值图转化为二维数组
            Bitmap bm = FastClipBitmap((Bitmap)b.Clone(), rect);
            bm = Bitize(bm, threshold, needGrayExtend);
            byte[,] srcGray = Image2Array(bm);
            bm.Dispose();
            // 轮廓提取
            byte[,] dstGray = ContourPick(srcGray);

            // 转换为灰度图像
            return Array2Image(dstGray, b.PixelFormat);
        } // end of ContourPick

        /// <summary>
        /// 区域周长
        /// </summary>
        /// <param name="Sign">二值图像标记数组</param>
        /// <returns></returns>
        private static int[] ImagePerimeter(ushort[,] Sign)
        {
            // 进行轮廓跟踪，划出轮廓线
            Sign = ContourTrace(Sign);

            int width = Sign.GetLength(0);
            int height = Sign.GetLength(1);

            // 找出最大标记号，即找出区域数
            int max = 0;
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    if (Sign[x, y] > max)
                        max = Sign[x, y];
                } // x
            } // y

            // 周长统计数组
            int[] Perimeter = new int[max + 1];

            // 计算区域周长
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    Perimeter[Sign[x, y]]++;
                } // x
            } // y

            return Perimeter;
        } // end of ImagePerimeter
        /// <summary>
        /// 获取每个区域的周长信息
        /// </summary>
        /// <param name="b">二值位图流</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExtend"></param>
        /// <returns></returns>
        public static int[] ImagePerimeter(Bitmap b, int blkORwht, byte threshold, bool needGrayExtend)
        {
            // 进行区域标记
            ushort[,] Sign = ImageSign((Bitmap)b.Clone(), blkORwht, threshold, needGrayExtend);

            // 区域周长
            int[] Perimeter = ImagePerimeter(Sign);

            return Perimeter;
        } // end of ImagePerimeter
        /// <summary>
        /// 获取指定区域中每个区域的周长信息
        /// </summary>
        /// <param name="b">二值位图流</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExtend"></param>
        /// <param name="rect"></param>
        /// <returns></returns>
        public static int[] ImagePerimeter(Bitmap b, int blkORwht, byte threshold, bool needGrayExtend, Rectangle rect)
        {
            // 进行区域标记
            ushort[,] Sign = ImageSign((Bitmap)b.Clone(), rect, blkORwht, threshold, needGrayExtend);

            // 区域周长
            int[] Perimeter = ImagePerimeter(Sign);

            return Perimeter;
        } // end of ImagePerimeter
        /// <summary>
        /// 区域面积
        /// </summary>
        /// <param name="Sign">二值图像标记数组</param>
        /// <returns></returns>
        private static int[] ImageArea(ushort[,] Sign)
        {
            int width = Sign.GetLength(0);
            int height = Sign.GetLength(1);

            // 找出最大标记号，即找出区域数
            int max = 0;
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    if (Sign[x, y] > max)
                        max = Sign[x, y];
                } // x
            } // y

            // 面积统计数组
            int[] Area = new int[max + 1];

            // 计算区域面积
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    Area[Sign[x, y]]++;
                } // x
            } // y

            return Area;
        } // end of ImageArea
        /// <summary>
        /// 标记一幅没有灰度化的位图图像
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="blkORwht">黑边还是白边检测</param>
        /// <param name="threshold">阀值</param>
        /// <param name="needGrayExtend">是否需要灰度拉伸</param>
        /// <returns></returns>
        public static ushort[,] ImageSign(Bitmap bmp, int blkORwht, byte threshold, bool needGrayExtend)
        {
            byte bw = 255;
            if (blkORwht == 0)
                bw = 0;
            Bitmap bm = Bitize((Bitmap)bmp.Clone(), threshold, needGrayExtend);
            byte[,] b = Image2Array(bm);
            bm.Dispose();
            int width = b.GetLength(0);
            int height = b.GetLength(1);

            // 标记号，最多可以标记 65536 个不同的连通区域
            // 注意标记号从 1 开始依次递增，标记号 0 代表背景
            ushort signNo = 1;

            // 用堆栈记录所有空标记
            System.Collections.Stack Seat = new System.Collections.Stack();

            // 二值图像连通区域标识，存储的是区域标识，而非图像数据
            ushort[,] Sign = new ushort[width, height];


            // 初始化最顶行标记
            for (int x = 0; x < width; x++)
            {
                if (b[x, 0] != bw) continue;

                // 处理所有连续的黑点
                while (x < width && b[x, 0] == bw)
                {
                    Sign[x, 0] = signNo;
                    x++;
                } // while

                signNo++;
            } // x


            // 处理最左列及最右列标记
            for (int y = 1; y < height; y++)
            {
                // 第左列
                if (b[0, y] == bw)
                {
                    if (b[0, y - 1] == bw)
                        Sign[0, y] = Sign[0, y - 1];
                    else
                        Sign[0, y] = signNo++;
                }

                // 最右列
                if (b[width - 1, y] == bw)
                {
                    if (b[width - 1, y - 1] == bw)
                        Sign[width - 1, y] = Sign[width - 1, y - 1];
                    else
                        Sign[width - 1, y] = signNo++;
                }
            } // y


            // 上面已经处理了图像最顶行、最左列、最右列，
            // 故下面只处理排除这三行后的主图像区
            int topRect = 1;
            int bottomRect = height;
            int leftRect = 1;
            int rightRect = width - 1;

            // 从左到右开始标记
            for (int y = topRect; y < bottomRect; y++)
            {
                for (int x = leftRect; x < rightRect; x++)
                {
                    // 如果当前点不为黑点，则跳过不处理
                    if (b[x, y] != bw) continue;

                    // 右上
                    if (b[x + 1, y - 1] == bw)
                    {
                        // 将当前点置为与右上点相同的标记
                        ushort sign = Sign[x, y] = Sign[x + 1, y - 1];

                        // 当左前点为黑点，且左前点的标记与右上点的标记不同时
                        if (b[x - 1, y] == bw && Sign[x - 1, y] != sign)
                        {
                            // 进栈：记录左前点标记，因为该标记将被替换掉
                            Seat.Push(Sign[x - 1, y]);

                            // 用右上点的标记号替换掉所有与左前点标记号相同的点
                            ReplaceSign(ref Sign, Sign[x - 1, y], sign);
                        }

                          // 当左上点为黑点，且左上点的标记与右上点的标记不同时
                        else if (b[x - 1, y - 1] == bw && Sign[x - 1, y - 1] != sign)
                        {
                            // 进栈：记录左上点标记，因为该标记将被替换掉
                            Seat.Push(Sign[x - 1, y - 1]);

                            // 用右上点的标记号替换掉所有与左上点标记号相同的点
                            ReplaceSign(ref Sign, Sign[x - 1, y - 1], sign);
                        }
                    } // 右上完

                        // 正上
                    else if (b[x, y - 1] == bw)
                    {
                        // 将当前点置为与正上点相同的标记
                        Sign[x, y] = Sign[x, y - 1];
                    }

                      // 左上
                    else if (b[x - 1, y - 1] == bw)
                    {
                        // 将当前点置为与左上点相同的标记
                        Sign[x, y] = Sign[x - 1, y - 1];
                    }

                      // 左前
                    else if (b[x - 1, y] == bw)
                    {
                        // 将当前点置为与左前点相同的标记
                        Sign[x, y] = Sign[x - 1, y];
                    }

                      // 右上、正上、左上及左前四点均不为黑点，即表示新区域开始
                    else
                    {
                        // 避免区域数超过 0xFFFF
                        if (signNo >= 0xFFFF)
                            return Sign;

                        // 如果堆栈里无空标记，
                        // 则使用新标记，否则使用该空标记
                        if (Seat.Count == 0)
                            Sign[x, y] = signNo++;
                        else
                            Sign[x, y] = (ushort)Seat.Pop(); // 出栈
                    }

                } // x
            } // y


            // 堆栈里存在空标记，则使用完没有空标记
            while (Seat.Count > 0)
            {
                ReplaceSign(ref Sign, (ushort)(--signNo), (ushort)Seat.Pop());
            } // while

            return Sign;
        } // end of ImageSign

        /// <summary>
        /// 在没有灰度化位图中获取每个区域的面积信息
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExtend">是否需要灰度拉伸</param>
        /// <returns></returns>
        public static int[] ImageArea(Bitmap b, int blkORwht, byte threshold, bool needGrayExtend)
        {
            // 进行区域标记
            ushort[,] Sign = ImageSign((Bitmap)b.Clone(), blkORwht, threshold, needGrayExtend);

            // 区域面积
            int[] Area = ImageArea(Sign);

            return Area;
        } // end of ImageArea

        /// <summary>
        /// 在没有灰度化位图的指定区域下获取每个区域的面积信息
        /// </summary>
        /// <param name="b">二值位图流</param>
        /// <param name="rect">指定区域</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExtend"></param>
        /// <returns></returns>
        public static int[] ImageArea(Bitmap b, Rectangle rect, int blkORwht, byte threshold, bool needGrayExtend)
        {
            // 进行区域标记
            ushort[,] Sign = ImageSign((Bitmap)b.Clone(), rect, blkORwht, threshold, needGrayExtend);

            // 区域面积
            int[] Area = ImageArea(Sign);

            return Area;
        } // end of ImageArea

        /// <summary>
        /// 在没有灰度化位图的指定路径下获取每个区域的面积信息
        /// </summary>
        /// <param name="b">二值位图流</param>
        /// <param name="path">指定路径</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExtend">是否需要灰度拉伸</param>
        /// <returns></returns>
        public static int[] ImageArea(Bitmap b, GraphicsPath path, int blkORwht, byte threshold, bool needGrayExtend)
        {
            // 进行区域标记
            ushort[,] Sign = ImageSign((Bitmap)b.Clone(), path, blkORwht, threshold, needGrayExtend);

            // 区域面积
            int[] Area = ImageArea(Sign);

            return Area;
        } // end of ImageArea

        /// <summary>
        /// 用新的标记号替换掉标记数组中旧的标记号
        /// </summary>
        /// <param name="Sign">二值图像标记数组</param>
        /// <param name="srcSign">原始标记号</param>
        /// <param name="dstSign">目标标记号</param>
        private static void ReplaceSign(ref ushort[,] Sign, ushort srcSign, ushort dstSign)
        {
            int width = Sign.GetLength(0);
            int height = Sign.GetLength(1);

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    if (Sign[x, y] == srcSign)
                        Sign[x, y] = dstSign;
                } // x
            } // y
        } // end of ReplaceSign

        /// <summary>
        /// 在指定的区域中标记一幅位图图像
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="rect">指定的区域</param>
        /// <param name="blkORwht">黑边还是白边检测</param>
        /// <param name="threshold">阀值</param>
        /// <param name="needGrayExtend">是否需要灰度拉伸</param>
        /// <returns></returns>
        public static ushort[,] ImageSign(Bitmap bmp, Rectangle rect, int blkORwht, byte threshold, bool needGrayExtend)
        {
            /*
            byte bw = 255;
            if (blkORwht == 0)
                bw = 0;
            // 将原始二值图转化为二维数组*/
            Bitmap outBmp = FastClipBitmap((Bitmap)bmp.Clone(), rect);
            return ImageSign(outBmp, blkORwht, threshold, needGrayExtend);
            /*outBmp = Bitize(outBmp, threshold, needGray, needGrayExtend);
            byte[,] b = Image2Array(outBmp);
            

            int X = 0;
            int Y = 0;
            int width = outBmp.Width;
            int height = outBmp.Height;
            outBmp.Dispose();
            // 标记号，最多可以标记 65536 个不同的连通区域
            // 注意标记号从 1 开始依次递增，标记号 0 代表背景
            ushort signNo = 1;

            // 用堆栈记录所有空标记
            System.Collections.Stack Seat = new System.Collections.Stack();

            // 二值图像连通区域标识，存储的是区域标识，而非图像数据
            ushort[,] Sign = new ushort[width, height];


            // 初始化最顶行标记
            for (int x = X; x < width; x++)
            {
                if (b[x, Y] != bw) continue;

                // 处理所有连续的黑点
                while (x < width && b[x, Y] == bw)
                {
                    Sign[x, Y] = signNo;
                    x++;
                } // while
                signNo++;
            } // x


            // 处理最左列及最右列标记
            for (int y = Y + 1; y < height; y++)
            {
                // 最左列
                if (b[X, y] == bw)
                {
                    if (b[X, y - 1] == bw)
                        Sign[X, y] = Sign[X, y - 1];
                    else
                        Sign[X, y] = signNo++;
                }
                // 最右列
                if (b[width - 1, y] == bw)
                {
                    if (b[width - 1, y - 1] == bw)
                        Sign[width - 1, y] = Sign[width - 1, y - 1];
                    else
                        Sign[width - 1, y] = signNo++;
                }
            } // y


            // 上面已经处理了图像最顶行、最左列、最右列，
            // 故下面只处理排除这三行后的主图像区
            int topRect = Y + 1;
            int bottomRect = height;
            int leftRect = X + 1;
            int rightRect = width - 1;

            // 从左到右开始标记
            for (int y = topRect; y < bottomRect; y++)
            {
                for (int x = leftRect; x < rightRect; x++)
                {
                    // 如果当前点不为黑点，则跳过不处理
                    if (b[x, y] != bw) continue;

                    // 右上
                    if (b[x + 1, y - 1] == bw)
                    {
                        // 将当前点置为与右上点相同的标记
                        ushort sign = Sign[x, y] = Sign[x + 1, y - 1];

                        // 当左前点为黑点，且左前点的标记与右上点的标记不同时
                        if (b[x - 1, y] == bw && Sign[x - 1, y] != sign)
                        {
                            // 进栈：记录左前点标记，因为该标记将被替换掉
                            Seat.Push(Sign[x - 1, y]);

                            // 用右上点的标记号替换掉所有与左前点标记号相同的点
                            ReplaceSign(ref Sign, Sign[x - 1, y], sign);
                        }

                          // 当左上点为黑点，且左上点的标记与右上点的标记不同时
                        else if (b[x - 1, y - 1] == bw && Sign[x - 1, y - 1] != sign)
                        {
                            // 进栈：记录左上点标记，因为该标记将被替换掉
                            Seat.Push(Sign[x - 1, y - 1]);

                            // 用右上点的标记号替换掉所有与左上点标记号相同的点
                            ReplaceSign(ref Sign, Sign[x - 1, y - 1], sign);
                        }
                    } // 右上完

                    // 正上
                    else if (b[x, y - 1] == bw)
                    {
                        // 将当前点置为与正上点相同的标记
                        Sign[x, y] = Sign[x, y - 1];
                    }

                  // 左上
                    else if (b[x - 1, y - 1] == bw)
                    {
                        // 将当前点置为与左上点相同的标记
                        Sign[x, y] = Sign[x - 1, y - 1];
                    }

                  // 左前
                    else if (b[x - 1, y] == bw)
                    {
                        // 将当前点置为与左前点相同的标记
                        Sign[x, y] = Sign[x - 1, y];
                    }

                  // 右上、正上、左上及左前四点均不为黑点，即表示新区域开始
                    else
                    {
                        // 避免区域数超过 0xFFFF
                        if (signNo >= 0xFFFF)
                            return Sign;
                        // 如果堆栈里无空标记，
                        // 则使用新标记，否则使用该空标记
                        if (Seat.Count == 0)
                            Sign[x, y] = signNo++;
                        else
                            Sign[x, y] = (ushort)Seat.Pop(); // 出栈
                    }

                } // x
            } // y


            // 堆栈里存在空标记，则使用完没有空标记
            while (Seat.Count > 0)
            {
                ReplaceSign(ref Sign, (ushort)(--signNo), (ushort)Seat.Pop());
            } // while

            return Sign;*/
        } // end of ImageSign

        /// <summary>
        /// 在指定的路径中标记一幅位图图像
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="path">指定的路径</param>
        /// <param name="blkORwht">黑边还是白边检测</param>
        /// <param name="threshold">阀值</param>
        /// <param name="needGrayExtend">是否需要灰度拉伸</param>
        /// <returns></returns>
        private static ushort[,] ImageSign(Bitmap bmp, GraphicsPath path, int blkORwht, byte threshold, bool needGrayExtend)
        {
            byte bw = 255;
            if (blkORwht == 0)
                bw = 0;
            // 将原始二值图转化为二维数组
            Rectangle rect = Rectangle.Truncate(path.GetBounds());
            Bitmap bm = FastClipBitmap((Bitmap)bmp.Clone(), rect);
            bm = Bitize(bm, threshold, needGrayExtend);
            byte[,] b = Image2Array(bm);

            int X = 0;
            int Y = 0;
            int width = bm.Width;
            int height = bm.Height;
            bm.Dispose();
            // 标记号，最多可以标记 65536 个不同的连通区域
            // 注意标记号从 1 开始依次递增，标记号 0 代表背景
            ushort signNo = 1;

            // 用堆栈记录所有空标记
            System.Collections.Stack Seat = new System.Collections.Stack();

            // 二值图像连通区域标识，存储的是区域标识，而非图像数据
            ushort[,] Sign = new ushort[width, height];


            // 初始化最顶行标记
            for (int x = X; x < width; x++)
            {
                if (b[x, Y] != bw) continue;
                // 处理所有连续的黑点
                while (x < width && b[x, Y] == bw)
                {
                    if (path.IsVisible(x + rect.X, Y + rect.Y))
                    {
                        Sign[x, Y] = signNo;
                    }
                    x++;
                } // while
                signNo++;
            } // x


            // 处理最左列及最右列标记
            for (int y = Y + 1; y < height; y++)
            {
                if (path.IsVisible(X + rect.X, y + rect.Y))
                {
                    // 最左列
                    if (b[X, y] == bw)
                    {
                        if (b[X, y - 1] == bw)
                            Sign[X, y] = Sign[X, y - 1];
                        else
                            Sign[X, y] = signNo++;
                    }
                }
                if (path.IsVisible(width - 1 + rect.X, y + rect.Y))
                {
                    // 最右列
                    if (b[width - 1, y] == bw)
                    {
                        if (b[width - 1, y - 1] == bw)
                            Sign[width - 1, y] = Sign[width - 1, y - 1];
                        else
                            Sign[width - 1, y] = signNo++;
                    }
                }
            } // y


            // 上面已经处理了图像最顶行、最左列、最右列，
            // 故下面只处理排除这三行后的主图像区
            int topRect = Y + 1;
            int bottomRect = height + X;
            int leftRect = X + 1;
            int rightRect = width + Y - 1;

            // 从左到右开始标记
            for (int y = topRect; y < bottomRect; y++)
            {
                for (int x = leftRect; x < rightRect; x++)
                {
                    // 如果当前点不为黑点，则跳过不处理
                    if (b[x, y] != bw) continue;

                    // 右上
                    if (b[x + 1, y - 1] == bw)
                    {
                        if (path.IsVisible(x + 1 + rect.X, y - 1 + rect.Y))
                        {
                            // 将当前点置为与右上点相同的标记
                            ushort sign = Sign[x, y] = Sign[x + 1, y - 1];

                            // 当左前点为黑点，且左前点的标记与右上点的标记不同时
                            if (b[x - 1, y] == bw && Sign[x - 1, y] != sign)
                            {
                                // 进栈：记录左前点标记，因为该标记将被替换掉
                                Seat.Push(Sign[x - 1, y]);

                                // 用右上点的标记号替换掉所有与左前点标记号相同的点
                                ReplaceSign(ref Sign, Sign[x - 1, y], sign);
                            }

                              // 当左上点为黑点，且左上点的标记与右上点的标记不同时
                            else if (b[x - 1, y - 1] == bw && Sign[x - 1, y - 1] != sign)
                            {
                                // 进栈：记录左上点标记，因为该标记将被替换掉
                                Seat.Push(Sign[x - 1, y - 1]);

                                // 用右上点的标记号替换掉所有与左上点标记号相同的点
                                ReplaceSign(ref Sign, Sign[x - 1, y - 1], sign);
                            }
                        }
                    } // 右上完

                        // 正上
                    else if (b[x, y - 1] == bw)
                    {
                        if (path.IsVisible(x + rect.X, y - 1 + rect.Y))
                        {
                            // 将当前点置为与正上点相同的标记
                            Sign[x, y] = Sign[x, y - 1];
                        }
                    }

                      // 左上
                    else if (b[x - 1, y - 1] == bw)
                    {
                        if (path.IsVisible(x - 1 + rect.X, y - 1 + rect.Y))
                        {
                            // 将当前点置为与左上点相同的标记
                            Sign[x, y] = Sign[x - 1, y - 1];
                        }
                    }

                      // 左前
                    else if (b[x - 1, y] == bw)
                    {
                        if (path.IsVisible(x - 1 + rect.X, y + rect.Y))
                        {
                            // 将当前点置为与左前点相同的标记
                            Sign[x, y] = Sign[x - 1, y];
                        }
                    }

                      // 右上、正上、左上及左前四点均不为黑点，即表示新区域开始
                    else
                    {
                        // 避免区域数超过 0xFFFF
                        if (signNo >= 0xFFFF)
                            return Sign;
                        if (path.IsVisible(x + rect.X, y + rect.Y))
                        {
                            // 如果堆栈里无空标记，
                            // 则使用新标记，否则使用该空标记
                            if (Seat.Count == 0)
                                Sign[x, y] = signNo++;
                            else
                                Sign[x, y] = (ushort)Seat.Pop(); // 出栈
                        }
                    }

                } // x
            } // y


            // 堆栈里存在空标记，则使用完没有空标记
            while (Seat.Count > 0)
            {
                ReplaceSign(ref Sign, (ushort)(--signNo), (ushort)Seat.Pop());
            } // while

            return Sign;
        } // end of ImageSign
        /// <summary>
        /// 在指定的位图中以指定的颜色绘画封闭区域的边框（或填充整块区域），并返回结果图
        /// </summary>
        /// <param name="b"></param>
        /// <param name="rect"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExtend"></param>
        /// <param name="Region"></param>
        /// <param name="fill"></param>
        /// <param name="color"></param>
        /// <param name="containBounds">是否包括边框</param>
        public static Bitmap GetSquareContourOrFillInternal(Bitmap b, Rectangle rect, int blkORwht, byte threshold, bool needGrayExtend, ushort[] Region, bool fill, Color color, bool containBounds)
        {
            Bitmap bb = (Bitmap)b.Clone();
            Bitmap bm = (Bitmap)b.Clone();
            // 进行区域标记
            ushort[,] Sign = ImageSign(bb, rect, blkORwht, threshold, needGrayExtend);
            bb.Dispose();
            // 按轮廓线进行显示
            if (!fill)
                Sign = ContourTrace(Sign);

            int len = Region.Length;

            GraphicsUnit gu = GraphicsUnit.Pixel;
            Rectangle rec = Rectangle.Intersect(rect, Rectangle.Truncate(b.GetBounds(ref gu)));
            int X = rec.X;
            int Y = rec.Y;
            int width = rec.Width;
            int height = rec.Height;
            int BPP = Image.GetPixelFormatSize(bm.PixelFormat) / 8;

            BitmapData data = bm.LockBits(new Rectangle(X, Y, width, height), ImageLockMode.ReadWrite, bm.PixelFormat);

            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;

                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        ushort sign = Sign[x, y];
                        bool showRegion = false;
                        Point pp = new Point(x + X, y + Y);
                        for (int i = 0; i < len; i++)
                        {
                            if (containBounds)
                            {
                                if (sign == Region[i])
                                {
                                    showRegion = true;
                                    break;
                                }
                            }
                            else
                            {
                                if (sign == Region[i] && pp.X != rec.X && pp.X != rec.Right - 1 && pp.Y != rec.Y && pp.Y != rec.Bottom - 1)
                                {
                                    showRegion = true;
                                    break;
                                }
                            }
                        } // i

                        // 绘制区域
                        if (showRegion)
                        {
                            p[0] = color.B;
                            p[1] = color.G;
                            p[2] = color.R;
                        }
                        p += BPP;
                    } // x

                    p += offset;
                } // y
            }
            bm.UnlockBits(data);
            return bm;
        }

        /// <summary>
        /// 在指定的位图中获取指定区域的边框（或填充整块区域）的点集
        /// </summary>
        /// <param name="b"></param>
        /// <param name="rect"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExtend"></param>
        /// <param name="Region"></param>
        /// <param name="fill"></param>
        /// <param name="containBounds">是否包括边框</param>
        public static List<Point> GetSquareContourOrFillInternal(Bitmap b, Rectangle rect, int blkORwht, byte threshold, bool needGrayExtend, ushort[] Region, bool fill, bool containBounds)
        {
            List<Point> ps = new List<Point>();
            Bitmap bb = (Bitmap)b.Clone();
            // 进行区域标记
            ushort[,] Sign = ImageSign(bb, rect, blkORwht, threshold, needGrayExtend);

            // 按轮廓线进行显示
            if (!fill)
                Sign = ContourTrace(Sign);

            int len = Region.Length;

            GraphicsUnit gu = GraphicsUnit.Pixel;
            Rectangle rec = Rectangle.Intersect(rect, Rectangle.Truncate(b.GetBounds(ref gu)));
            int X = rec.X;
            int Y = rec.Y;
            int width = rec.Width;
            int height = rec.Height;

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    ushort sign = Sign[x, y];
                    bool showRegion = false;
                    Point p = new Point(x + X, y + Y);
                    for (int i = 0; i < len; i++)
                    {
                        if (containBounds)
                        {
                            if (sign == Region[i])
                            {
                                showRegion = true;
                                break;
                            }
                        }
                        else
                        {
                            if (sign == Region[i] && p.X != rec.X && p.X != rec.Right - 1 && p.Y != rec.Y && p.Y != rec.Bottom - 1)
                            {
                                showRegion = true;
                                break;
                            }
                        }
                    }

                    // 绘制区域
                    if (showRegion)
                    {
                        ps.Add(new Point(x + X, y + Y));
                    }
                }
            }
            return ps;
        }
        /// <summary>
        /// 在指定的位图中获取指定区域的边框（或填充整块区域）的点集，以及其他非指定区域的点集(即otherSquares)
        /// </summary>
        /// <param name="b"></param>
        /// <param name="rect"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExtend"></param>
        /// <param name="Region"></param>
        /// <param name="fill"></param>
        /// <param name="otherSquares"></param>
        /// <param name="containBounds">是否包含边框</param>
        public static List<Point> GetSquareContourOrFillInternal(Bitmap b, Rectangle rect, int blkORwht, byte threshold, bool needGrayExtend, ushort[] Region, bool fill, out List<Point> otherSquares, bool containBounds)
        {
            List<Point> al = new List<Point>();
            List<Point> al1 = new List<Point>();
            Bitmap bb = (Bitmap)b.Clone();
            // 进行区域标记
            ushort[,] Sign = ImageSign(bb, rect, blkORwht, threshold, needGrayExtend);
            int[] area = ImageArea(Sign);
            // 按轮廓线进行显示
            if (!fill)
                Sign = ContourTrace(Sign);
            List<ushort> others = new List<ushort>();
            //int accessed = 0;
            for (int s = 1; s < area.Length; s++)//去掉第一个背景区域s从1开始
            {
                //if (s == accessed)
                //continue;
                bool same = false;
                for (int r = 0; r < Region.Length; r++)
                {
                    //if (s != Region[r])
                    //{
                    //    others.Add((ushort)s);
                    //    accessed = Region[r];
                    //}
                    if (s == Region[r])
                        same = true;
                }
                if (!same)
                    others.Add((ushort)s);
            }
            ushort[] OtherRegion = others.ToArray();

            GraphicsUnit gu = GraphicsUnit.Pixel;
            Rectangle rec = Rectangle.Intersect(rect, Rectangle.Truncate(b.GetBounds(ref gu)));
            int X = rec.X;
            int Y = rec.Y;
            int width = rec.Width;
            int height = rec.Height;

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    ushort sign = Sign[x, y];
                    Point p = new Point(x + X, y + Y);
                    for (int i = 0; i < Region.Length; i++)
                    {
                        if (containBounds)
                        {
                            if (sign == Region[i])
                            {
                                al.Add(new Point(x + X, y + Y));
                                break;
                            }
                        }
                        else
                        {
                            if (sign == Region[i] && p.X != rec.X && p.X != rec.Right - 1 && p.Y != rec.Y && p.Y != rec.Bottom - 1)
                            {
                                al.Add(new Point(x + X, y + Y));
                                break;
                            }
                        }
                    } // i
                    for (int i = 0; i < OtherRegion.Length; i++)
                    {
                        if (sign == OtherRegion[i])
                        {
                            al1.Add(new Point(x + X, y + Y));
                            break;
                        }
                    }
                }
            }
            otherSquares = al1;
            return al;
        }
        /// <summary>
        /// 在指定的位图中以指定的颜色绘画封闭路径的边框（或填充整个路径），并返回结果图
        /// </summary>
        /// <param name="b"></param>
        /// <param name="path"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExtend"></param>
        /// <param name="Region"></param>
        /// <param name="fill"></param>
        /// <param name="color"></param>
        /// <param name="containBounds">是否包括边框</param>
        public static Bitmap GetSquareContourOrFillInternal(Bitmap b, GraphicsPath path, int blkORwht, byte threshold, bool needGrayExtend, ushort[] Region, bool fill, Color color, bool containBounds)
        {
            Bitmap bb = (Bitmap)b.Clone();
            Bitmap bm = (Bitmap)b.Clone();
            // 进行区域标记
            ushort[,] Sign = ImageSign(bb, path, blkORwht, threshold, needGrayExtend);
            bb.Dispose();
            // 按轮廓线进行显示
            if (!fill)
                Sign = ContourTrace(Sign);

            int len = Region.Length;

            GraphicsUnit gu = GraphicsUnit.Pixel;
            Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(path.GetBounds()), Rectangle.Truncate(b.GetBounds(ref gu)));
            int X = rec.X;
            int Y = rec.Y;
            int width = rec.Width;
            int height = rec.Height;
            int BPP = Image.GetPixelFormatSize(bm.PixelFormat) / 8;

            BitmapData data = bm.LockBits(new Rectangle(X, Y, width, height), ImageLockMode.ReadWrite, bm.PixelFormat);

            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;

                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        ushort sign = Sign[x, y];
                        bool showRegion = false;
                        Point pp = new Point(x + X, y + Y);
                        for (int i = 0; i < len; i++)
                        {
                            if (containBounds)
                            {
                                if (sign == Region[i])
                                {
                                    showRegion = true;
                                    break;
                                }
                            }
                            else
                            {
                                if (sign == Region[i] && pp.X != rec.X && pp.X != rec.Right - 1 && pp.Y != rec.Y && pp.Y != rec.Bottom - 1)
                                {
                                    showRegion = true;
                                    break;
                                }
                            }
                        } // i

                        // 绘制区域
                        if (showRegion)
                        {
                            p[0] = color.B;
                            p[1] = color.G;
                            p[2] = color.R;
                        }
                        p += BPP;
                    } // x

                    p += offset;
                } // y
            }
            bm.UnlockBits(data);
            return bm;
        }

        /// <summary>
        /// 在指定的位图中获取指定路径的边框（或填充整个路径）的点集
        /// </summary>
        /// <param name="b"></param>
        /// <param name="path"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExtend"></param>
        /// <param name="Region"></param>
        /// <param name="fill"></param>
        /// <param name="containBounds">是否包含边框</param>
        public static List<Point> GetSquareContourOrFillInternal(Bitmap b, GraphicsPath path, int blkORwht, byte threshold, bool needGrayExtend, ushort[] Region, bool fill, bool containBounds)
        {
            List<Point> ps = new List<Point>();
            Bitmap bb = (Bitmap)b.Clone();
            // 进行区域标记
            ushort[,] Sign = ImageSign(bb, path, blkORwht, threshold, needGrayExtend);

            // 按轮廓线进行显示
            if (!fill)
                Sign = ContourTrace(Sign);

            int len = Region.Length;

            GraphicsUnit gu = GraphicsUnit.Pixel;
            Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(path.GetBounds()), Rectangle.Truncate(b.GetBounds(ref gu)));
            int X = rec.X;
            int Y = rec.Y;
            int width = rec.Width;
            int height = rec.Height;

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    ushort sign = Sign[x, y];
                    bool showRegion = false;
                    Point p = new Point(x + X, y + Y);
                    for (int i = 0; i < len; i++)
                    {
                        if (containBounds)
                        {
                            if (sign == Region[i])
                            {
                                showRegion = true;
                                break;
                            }
                        }
                        else
                        {
                            if (sign == Region[i] && p.X != rec.X && p.X != rec.Right - 1 && p.Y != rec.Y && p.Y != rec.Bottom - 1)
                            {
                                showRegion = true;
                                break;
                            }
                        }
                    }

                    // 绘制区域
                    if (showRegion)
                    {
                        ps.Add(new Point(x + X, y + Y));
                    }
                }
            }
            return ps;
        }
        /// <summary>
        /// 在指定的位图中获取指定路径的边框（或填充整个路径）的点集，以及其他非指定路径的点集(即otherSquares)
        /// </summary>
        /// <param name="b"></param>
        /// <param name="path"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExtend"></param>
        /// <param name="Region"></param>
        /// <param name="fill"></param>
        /// <param name="otherSquares"></param>
        /// <param name="containBounds">是否包括边框</param>
        public static List<Point> GetSquareContourOrFillInternal(Bitmap b, GraphicsPath path, int blkORwht, byte threshold, bool needGrayExtend, ushort[] Region, bool fill, out List<Point> otherSquares, bool containBounds)
        {
            List<Point> al = new List<Point>();
            List<Point> al1 = new List<Point>();
            Bitmap bb = (Bitmap)b.Clone();
            // 进行区域标记
            ushort[,] Sign = ImageSign(bb, path, blkORwht, threshold, needGrayExtend);
            int[] area = ImageArea(Sign);
            // 按轮廓线进行显示
            if (!fill)
                Sign = ContourTrace(Sign);
            List<ushort> others = new List<ushort>();
            //int accessed = 0;
            for (int s = 1; s < area.Length; s++)//去掉第一个背景区域s从1开始
            {
                //if (s == accessed)
                //continue;
                bool same = false;
                for (int r = 0; r < Region.Length; r++)
                {
                    //if (s != Region[r])
                    //{
                    //    others.Add((ushort)s);
                    //    accessed = Region[r];
                    //}
                    if (s == Region[r])
                        same = true;
                }
                if (!same)
                    others.Add((ushort)s);
            }
            ushort[] OtherRegion = others.ToArray();
            Rectangle rect = Rectangle.Truncate(path.GetBounds());
            GraphicsUnit gu = GraphicsUnit.Pixel;
            Rectangle rec = Rectangle.Intersect(rect, Rectangle.Truncate(b.GetBounds(ref gu)));
            int X = rec.X;
            int Y = rec.Y;
            int width = rec.Width;
            int height = rec.Height;

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    ushort sign = Sign[x, y];
                    Point p = new Point(x + X, y + Y);
                    for (int i = 0; i < Region.Length; i++)
                    {
                        if (containBounds)
                        {
                            if (sign == Region[i])
                            {
                                al.Add(p);
                                break;
                            }
                        }
                        else
                        {
                            if (sign == Region[i] && p.X != rec.X && p.X != rec.Right - 1 && p.Y != rec.Y && p.Y != rec.Bottom - 1)
                            {
                                al.Add(p);
                                break;
                            }
                        }
                    } // i
                    for (int i = 0; i < OtherRegion.Length; i++)
                    {
                        if (sign == OtherRegion[i])
                        {
                            al1.Add(p);
                            break;
                        }
                    }
                }
            }
            otherSquares = al1;
            return al;
        }
        /// <summary>
        /// 获取窗体的外观位图
        /// </summary>
        /// <param name="form">窗体</param>
        /// <param name="showAll">是否包括标题栏和边框</param>
        /// <returns></returns>
        public static Bitmap GetFormImage(Form form, bool showAll)
        {
            int edge = 0;
            int x = 0;
            int y = 0;
            int wid = form.Width;
            int hei = form.Height;
            if (!showAll)
            {
                if (form.FormBorderStyle != FormBorderStyle.None)
                {
                    edge = (form.Width - form.ClientSize.Width) / 2;
                    x = edge + 1;
                    y = form.Height - form.ClientSize.Height - edge + 1;
                    wid = form.ClientSize.Width - 1;
                    hei = form.ClientSize.Height - 1;
                }
            }
            Bitmap bmp = new Bitmap(form.Width, form.Height);
            form.DrawToBitmap(bmp, new Rectangle(0, 0, bmp.Width, bmp.Height));
            Rectangle dstRect = new Rectangle(x, y, wid, hei);
            Bitmap b = new Bitmap(wid, hei);
            using (Graphics g = Graphics.FromImage(b))
            {
                g.DrawImage(bmp, new Rectangle(0, 0, b.Width, b.Height), dstRect, GraphicsUnit.Pixel);
            }
            bmp.Dispose();
            return b;
        }

        /// <summary>
        /// 获取控件的外观位图
        /// </summary>
        /// <param name="control">控件</param>
        /// <param name="showAll">是否包括边框</param>
        /// <returns></returns>
        public static Bitmap GetControlImage(Control control, bool showAll)
        {
            int edge = 0;
            int x = 0;
            int y = 0;
            int wid = control.Width;
            int hei = control.Height;
            if (!showAll)
            {
                if (control.ClientSize.Width != control.Width || control.ClientSize.Height != control.Height)
                {
                    edge = (control.Width - control.ClientSize.Width) / 2;
                    x = edge + 1;
                    y = control.Height - control.ClientSize.Height - edge + 1;
                    wid = control.ClientSize.Width - 1;
                    hei = control.ClientSize.Height - 1;
                }
            }
            Bitmap bmp = new Bitmap(control.Width, control.Height);
            control.DrawToBitmap(bmp, new Rectangle(0, 0, bmp.Width, bmp.Height));
            Rectangle dstRect = new Rectangle(x, y, wid, hei);
            Bitmap b = new Bitmap(wid, hei);
            using (Graphics g = Graphics.FromImage(b))
            {
                g.DrawImage(bmp, new Rectangle(0, 0, b.Width, b.Height), dstRect, GraphicsUnit.Pixel);
            }
            bmp.Dispose();
            return b;
        }

        /// <summary>
        /// 取得一个图片中非透明色部分的区域
        /// </summary>
        /// <param name="b">取其区域的图片</param>
        /// <param name="transparentColor">透明色</param>
        /// <returns>图片中非透明色部分的区域</returns>
        public static Region GetRegionFromBmp(Bitmap b, Color transparentColor)
        {
            int width = b.Width;
            int height = b.Height;
            Region rgn = new Region();
            rgn.MakeEmpty();
            bool isTransRgn;//前一个点是否在透明区
            Color curColor;//当前点的颜色
            Rectangle curRect = new Rectangle();
            curRect.Height = 1;
            unsafe
            {
                int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
                BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, b.PixelFormat);
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - BPP * width;
                //逐像素扫描这个图片，找出非透明色部分区域并合并起来。
                for (int y = 0; y < height; y++)
                {
                    isTransRgn = true;
                    for (int x = 0; x < width; x++)
                    {
                        curColor = Color.FromArgb(p[2], p[1], p[0]);
                        if (curColor == transparentColor || x == width - 1)//如果遇到透明色或行尾
                        {
                            if (isTransRgn == false)//退出有效区
                            {
                                curRect.Width = x - curRect.X;
                                rgn.Union(curRect);
                            }
                        }
                        else//非透明色
                        {
                            if (isTransRgn == true)//进入有效区
                            {
                                curRect.X = x;
                                curRect.Y = y;
                            }
                        }//if curColor
                        isTransRgn = curColor == transparentColor;
                        p += BPP;
                    }//for x
                    p += offset;
                }//for y
                b.UnlockBits(data);
            }
            return rgn;
        }

        /// <summary>
        /// 在指定的区域中无方向检测出线段，成功检测返回true，否则返回false，无问题！
        /// </summary>
        /// <param name="rgn"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        public bool GetLineByRng(Region rgn, int blkORwht, byte threshold, bool needGrayExt, out PointF p1, out PointF p2)
        {
            bool flag = false;
            List<PointF> al = new List<PointF>();
            p1 = PointF.Empty;
            p2 = PointF.Empty;
            if (bmp != null)
            {
                Graphics g = Graphics.FromImage((Image)bmp.Clone());
                RectangleF rect = rgn.GetBounds(g);
                Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(rect), new Rectangle(0, 0, width, height));
                int X = rec.X;
                int Y = rec.Y;
                int Width = rec.Width + rec.X;
                int Height = rec.Height + rec.Y;

                Bitmap b = Thresholding(bmp, rec, threshold, needGrayExt);
                Bitmap edgeImage = GetBmpByRoberts(b);
                b.Dispose();
                unsafe
                {
                    int BPP = Image.GetPixelFormatSize(edgeImage.PixelFormat) / 8;
                    BitmapData data = edgeImage.LockBits(new Rectangle(0, 0, rec.Width, rec.Height), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                    byte* p = (byte*)data.Scan0;
                    int stride = data.Stride;
                    int offset = stride - BPP * rec.Width;
                    for (int j = Y; j < Height; j++)
                    {
                        for (int i = X; i < Width; i++)
                        {
                            PointF pp = new PointF(i, j);
                            if (rgn.IsVisible(pp))
                            {
                                if (blkORwht == 0)
                                {
                                    if (p[0] == 0)
                                        al.Add(pp);
                                }
                                else
                                {
                                    if (p[0] == 255)
                                        al.Add(pp);
                                }
                            }
                            p += BPP;
                        }
                        p += offset;
                    }
                    edgeImage.UnlockBits(data);
                }
                edgeImage.Dispose();
                //if (ps != null && ps.Length > 0)
                if (al.Count > 0)
                {
                    Algorithm NiHeXian = new Algorithm();
                    NiHeXian.NiHeLine(al, out p1, out p2);
                    flag = true;
                }
            }
            return flag;
        }
        /// <summary>
        /// 从外到内或由内而外的方向枚举
        /// </summary>
        public enum InOutDirection
        {
            /// <summary>
            /// 从外到内
            /// </summary>
            In,
            /// <summary>
            /// 由内而外
            /// </summary>
            Out,
            /// <summary>
            /// 其它方向
            /// </summary>
            Other
        }
        /// <summary>
        /// 在指定的区域中有方向检测出圆，成功检测返回true，否则返回false
        /// </summary>
        /// <param name="rgn"></param>
        /// <param name="dir">如果可以尽量采用In方向（即从外到内），准确度会更高</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="center"></param>
        /// <param name="rand"></param>
        public bool GetCircleByRng(Region rgn, InOutDirection dir, int blkORwht, byte threshold, bool needGrayExt, out PointF center, out double rand)
        {
            bool flag = false;
            center = PointF.Empty;
            rand = 0;
            if (bmp != null)
            {
                Graphics g = Graphics.FromImage((Image)bmp.Clone());
                RectangleF rect = rgn.GetBounds(g);
                Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(rect), new Rectangle(0, 0, width, height));
                int X = rec.X;
                int Y = rec.Y;
                int Width = rec.Width + rec.X;
                int Height = rec.Height + rec.Y;

                Bitmap b = Thresholding(bmp, rec, threshold, needGrayExt);
                Bitmap edgeImage = GetBmpByRoberts(b);
                int BPP = Image.GetPixelFormatSize(edgeImage.PixelFormat) / 8;
                b.Dispose();
                List<PointF> al1 = new List<PointF>();
                List<PointF> al2 = new List<PointF>();
                unsafe
                {
                    if (dir == InOutDirection.In)
                    {
                        BitmapData data1 = edgeImage.LockBits(new Rectangle(0, 0, rec.Width / 2, rec.Height), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        byte* p1 = (byte*)data1.Scan0;
                        int stride1 = data1.Stride;
                        int offset = stride1 - BPP * rec.Width / 2;
                        for (int j = Y; j < Y + rec.Height; j++)
                        {
                            p1 = (byte*)data1.Scan0 + (j - Y) * stride1;
                            for (int i = X; i < X + rec.Width / 2; i++)
                            {
                                Point pp = new Point(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p1[0] == 0)
                                        {
                                            al1.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p1[0] == 255)
                                        {
                                            al1.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p1 += BPP;
                            }
                            p1 += offset;
                        }
                        edgeImage.UnlockBits(data1);

                        BitmapData data2 = edgeImage.LockBits(new Rectangle(rec.Width / 2, 0, rec.Width / 2, rec.Height), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        byte* p2 = (byte*)data2.Scan0 + BPP * rec.Width;
                        int stride2 = data2.Stride;
                        for (int j = Y; j < Y + rec.Height; j++)
                        {
                            p2 = (byte*)data2.Scan0 + BPP * rec.Width + (j - Y) * stride2;
                            for (int i = X + rec.Width - 1; i > X + rec.Width / 2; i--)
                            {
                                Point pp = new Point(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p2[0] == 0)
                                        {
                                            al2.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p2[0] == 255)
                                        {
                                            al2.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p2 -= BPP;
                            }
                            p2 += stride2 + BPP * rec.Width / 2;
                        }
                        edgeImage.UnlockBits(data2);
                    }
                    else if (dir == InOutDirection.Out)
                    {
                        BitmapData data3 = edgeImage.LockBits(new Rectangle(0, 0, rec.Width / 2, rec.Height), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        byte* p3 = (byte*)data3.Scan0 + BPP * rec.Width;
                        int stride3 = data3.Stride;
                        for (int j = Y; j < Y + rec.Height; j++)
                        {
                            p3 = (byte*)data3.Scan0 + BPP * rec.Width + (j - Y) * stride3;
                            for (int i = X + rec.Width - 1; i > X + rec.Width / 2; i--)
                            {
                                Point pp = new Point(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p3[0] == 0)
                                        {
                                            al1.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p3[0] == 255)
                                        {
                                            al1.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p3 -= BPP;
                            }
                            p3 += stride3 + BPP * rec.Width / 2;
                        }
                        edgeImage.UnlockBits(data3);

                        BitmapData data4 = edgeImage.LockBits(new Rectangle(rec.Width / 2, 0, rec.Width / 2, rec.Height), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        byte* p4 = (byte*)data4.Scan0;
                        int stride4 = data4.Stride;
                        int offset = stride4 - BPP * rec.Width / 2;
                        for (int j = Y; j < Y + rec.Height; j++)
                        {
                            p4 = (byte*)data4.Scan0 + (j - Y) * stride4;
                            for (int i = X; i < X + rec.Width / 2; i++)
                            {
                                Point pp = new Point(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p4[0] == 0)
                                        {
                                            al2.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p4[0] == 255)
                                        {
                                            al2.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p4 += BPP;
                            }
                            p4 += offset;
                        }
                        edgeImage.UnlockBits(data4);

                    }
                }
                edgeImage.Dispose();
                PointF[] ps1 = al1.ToArray();
                PointF[] ps2 = al2.ToArray();
                List<PointF> al = new List<PointF>();
                //PointF[] ps = new PointF[ps1.Length + ps2.Length];
                for (int i = 0; i < ps1.Length; i++)
                {
                    //ps[i] = ps1[i]; 
                    al.Add(ps1[i]);
                }
                for (int i = 0; i < ps2.Length; i++)
                {
                    //ps[ps1.Length + i] = ps2[i];
                    al.Add(ps2[i]);
                }
                Algorithm NiHeYuan = new Algorithm();

                NiHeYuan.NiHeCircle(al, out center, out rand);
                flag = true;
            }
            return flag;
        }
        /// <summary>
        /// 在指定的区域中有方向检测出弧，成功检测返回true，否则返回false
        /// </summary>
        /// <param name="rgn"></param>
        /// <param name="dir"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="center"></param>
        /// <param name="rand"></param>
        /// <param name="startA"></param>
        /// <param name="sweepA"></param>
        public bool GetArcByRng(Region rgn, InOutDirection dir, int blkORwht, byte threshold, bool needGrayExt, out PointF center, out double rand, out float startA, out float sweepA)//, out PointF pp1, out PointF pp2, out PointF pp3)
        {
            bool flag = false;
            center = PointF.Empty;
            rand = 0;
            startA = 0;
            sweepA = 0;// pp1 = PointF.Empty; pp2 = PointF.Empty; pp3 = PointF.Empty;
            if (bmp != null)
            {
                Graphics g = Graphics.FromImage((Image)bmp.Clone());
                RectangleF rect = rgn.GetBounds(g);
                Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(rect), new Rectangle(0, 0, width, height));
                int X = rec.X;
                int Y = rec.Y;
                //int Width = rec.Width + rec.X;
                //int Height = rec.Height + rec.Y;

                Bitmap b = Thresholding(bmp, rec, threshold, needGrayExt);
                Bitmap edgeImage = GetBmpByRoberts(b);
                int BPP = Image.GetPixelFormatSize(edgeImage.PixelFormat) / 8;
                b.Dispose();
                List<PointF> al1 = new List<PointF>();
                List<PointF> al2 = new List<PointF>();
                unsafe
                {
                    if (dir == InOutDirection.In)
                    {
                        BitmapData data1 = edgeImage.LockBits(new Rectangle(0, 0, rec.Width / 2, rec.Height), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        byte* p1 = (byte*)data1.Scan0;
                        int stride1 = data1.Stride;
                        int offset = stride1 - BPP * rec.Width / 2;
                        for (int j = Y; j < Y + rec.Height; j++)
                        {
                            p1 = (byte*)data1.Scan0 + (j - Y) * stride1;
                            for (int i = X; i < X + rec.Width / 2; i++)
                            {
                                Point pp = new Point(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p1[0] == 0)
                                        {
                                            al1.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p1[0] == 255)
                                        {
                                            al1.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p1 += BPP;
                            }
                            p1 += offset;
                        }
                        edgeImage.UnlockBits(data1);

                        BitmapData data2 = edgeImage.LockBits(new Rectangle(rec.Width / 2, 0, rec.Width / 2, rec.Height), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        byte* p2 = (byte*)data2.Scan0 + BPP * rec.Width / 2;
                        int stride2 = data2.Stride;
                        for (int j = Y; j < Y + rec.Height; j++)
                        {
                            p2 = (byte*)data2.Scan0 + BPP * rec.Width / 2 + (j - Y) * stride2;
                            for (int i = X + rec.Width - 1; i > X + rec.Width / 2; i--)
                            {
                                Point pp = new Point(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p2[0] == 0)
                                        {
                                            al2.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p2[0] == 255)
                                        {
                                            al2.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p2 -= BPP;
                            }
                            p2 += stride2 + BPP * rec.Width / 2;
                        }
                        edgeImage.UnlockBits(data2);
                    }
                    else if (dir == InOutDirection.Out)
                    {
                        BitmapData data3 = edgeImage.LockBits(new Rectangle(0, 0, rec.Width / 2, rec.Height), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        byte* p3 = (byte*)data3.Scan0 + BPP * rec.Width / 2;
                        int stride3 = data3.Stride;
                        for (int j = Y; j < Y + rec.Height; j++)
                        {
                            p3 = (byte*)data3.Scan0 + BPP * rec.Width / 2 + (j - Y) * stride3;
                            for (int i = X + rec.Width - 1; i > X + rec.Width / 2; i--)
                            {
                                Point pp = new Point(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p3[0] == 0)
                                        {
                                            al1.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p3[0] == 255)
                                        {
                                            al1.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p3 -= BPP;
                            }
                            p3 += stride3 + BPP * rec.Width / 2;
                        }
                        edgeImage.UnlockBits(data3);

                        BitmapData data4 = edgeImage.LockBits(new Rectangle(rec.Width / 2, 0, rec.Width / 2, rec.Height), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        byte* p4 = (byte*)data4.Scan0;
                        int stride4 = data4.Stride;
                        int offset = stride4 - BPP * rec.Width / 2;
                        for (int j = Y; j < Y + rec.Height; j++)
                        {
                            p4 = (byte*)data4.Scan0 + (j - Y) * stride4;
                            for (int i = X; i < X + rec.Width / 2; i++)
                            {
                                Point pp = new Point(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p4[0] == 0)
                                        {
                                            al2.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p4[0] == 255)
                                        {
                                            al2.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p4 += BPP;
                            }
                            p4 += offset;
                        }
                        edgeImage.UnlockBits(data4);
                    }
                }
                edgeImage.Dispose();
                PointF[] ps1 = al1.ToArray();
                PointF[] ps2 = al2.ToArray();
                List<PointF> al = new List<PointF>();
                //PointF[] ps = new PointF[ps1.Length + ps2.Length];
                for (int i = 0; i < ps1.Length; i++)
                {
                    //ps[i] = ps1[i];
                    al.Add(ps1[i]);
                }
                for (int i = 0; i < ps2.Length; i++)
                {
                    //ps[ps1.Length + i] = ps2[i];
                    al.Add(ps2[i]);
                }
                if (al.Count > 0)
                {
                    //float x, y;
                    PointF startP, middleP, endP;

                    Algorithm NiHeHu = new Algorithm();
                    NiHeHu.NiHeCircle(al, out center, out rand);
                    GetArc3P1(ps1, ps2, out startP, out middleP, out endP);
                    //pp1 = startP;
                    //pp2 = middleP;
                    //pp3 = endP;
                    GetArcPara(center, startP, middleP, endP, out startA, out sweepA);
                    //PointF[] points = new PointF[3];
                    //points[0] = startP;
                    //points[1] = middleP;
                    //points[2] = endP;
                    //ClassLibrary1.MyClass HuClass = new ClassLibrary1.MyClass();
                    //HuClass.HuaHu(x, y, points, rand, out sweepA, out startA);
                    //NiHeHu.Calculate(ps, out x, out y, out rand);
                    /*//屏蔽掉的方法有问题...
                    RectangleF rectF;
                    GetArc(ps, out rectF, out startA, out sweepA);有问题...
                    x = rectF.X + rectF.Width / 2;
                    y = rectF.Y + rectF.Height / 2;
                    rand = (rectF.Width + rectF.Height) / 4;
                    GetArc3P(rectF, startA, sweepA, out startP, out middleP, out endP);
                    */
                    flag = true;
                }
            }
            return flag;
        }
        /// <summary>
        /// 在指定的区域中有方向检测出矩形，成功检测返回true，否则返回false
        /// </summary>
        /// <param name="rgn"></param>
        /// <param name="dir"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="pp1"></param>
        /// <param name="pp2"></param>
        /// <param name="pp3"></param>
        /// <param name="pp4"></param>
        /// <param name="rectPs"></param>
        /// <returns></returns>
        public bool GetRectByRng(Region rgn, InOutDirection dir, int blkORwht, byte threshold, bool needGrayExt, out PointF pp1, out PointF pp2, out PointF pp3, out PointF pp4, out List<PointF> rectPs)
        {
            bool flag = false;
            pp1 = PointF.Empty;
            pp2 = PointF.Empty;
            pp3 = PointF.Empty;
            pp4 = PointF.Empty; rectPs = null;
            if (bmp != null)
            {
                Graphics g = Graphics.FromImage((Image)bmp.Clone());
                RectangleF rect = rgn.GetBounds(g);
                Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(rect), new Rectangle(0, 0, width, height));
                int X = rec.X;
                int Y = rec.Y;
                //int Width = rec.Width + rec.X;
                //int Height = rec.Height + rec.Y;

                Bitmap b = Thresholding(bmp, rec, threshold, needGrayExt);
                Bitmap edgeImage = GetBmpByRoberts(b);
                b.Dispose();
                List<PointF> al1 = new List<PointF>();
                List<PointF> al2 = new List<PointF>();
                List<PointF> al3 = new List<PointF>();
                List<PointF> al4 = new List<PointF>();
                unsafe
                {
                    #region//左右分割
                    if (dir == InOutDirection.In)
                    {
                        BitmapData data1 = edgeImage.LockBits(new Rectangle(0, 0, rec.Width / 2, rec.Height), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        byte* p1 = (byte*)data1.Scan0;
                        int stride1 = data1.Stride;
                        int offset = stride1 - BPP * rec.Width / 2;
                        for (int j = Y; j < Y + rec.Height; j++)
                        {
                            p1 = (byte*)data1.Scan0 + (j - Y) * stride1;
                            for (int i = X; i < X + rec.Width / 2; i++)
                            {
                                PointF pp = new PointF(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p1[0] == 0)
                                        {
                                            al1.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p1[0] == 255)
                                        {
                                            al1.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p1 += BPP;
                            }
                            p1 += offset;
                        }
                        edgeImage.UnlockBits(data1);

                        BitmapData data2 = edgeImage.LockBits(new Rectangle(rec.Width / 2, 0, rec.Width / 2, rec.Height), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        byte* p2 = (byte*)data2.Scan0;
                        int stride2 = data2.Stride;
                        for (int j = Y; j < Y + rec.Height; j++)
                        {
                            p2 = (byte*)data2.Scan0 + (j - Y) * stride2;
                            for (int i = X + rec.Width - 1; i > X + rec.Width / 2; i--)
                            {
                                PointF pp = new PointF(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p2[0] == 0)
                                        {
                                            al2.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p2[0] == 255)
                                        {
                                            al2.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p2 -= BPP;
                            }
                            p2 += stride2 + BPP * rec.Width / 2;
                        }
                        edgeImage.UnlockBits(data2);
                    }
                    else if (dir == InOutDirection.Out)
                    {
                        BitmapData data3 = edgeImage.LockBits(new Rectangle(0, 0, rec.Width / 2, rec.Height), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        byte* p3 = (byte*)data3.Scan0 + BPP * rec.Width / 2;
                        int stride3 = data3.Stride;
                        for (int j = Y; j < Y + rec.Height; j++)
                        {
                            p3 = (byte*)data3.Scan0 + BPP * rec.Width / 2 + (j - Y) * stride3;
                            for (int i = X + rec.Width - 1; i > X + rec.Width / 2; i--)
                            {
                                PointF pp = new PointF(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p3[0] == 0)
                                        {
                                            al1.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p3[0] == 255)
                                        {
                                            al1.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p3 -= BPP;
                            }
                            p3 += stride3 + BPP * rec.Width / 2;
                        }
                        edgeImage.UnlockBits(data3);

                        BitmapData data4 = edgeImage.LockBits(new Rectangle(rec.Width / 2, 0, rec.Width / 2, rec.Height), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        byte* p4 = (byte*)data4.Scan0;
                        int stride4 = data4.Stride;
                        int offset = stride4 - BPP * rec.Width / 2;
                        for (int j = Y; j < Y + rec.Height; j++)
                        {
                            p4 = (byte*)data4.Scan0 + (j - Y) * stride4;
                            for (int i = X; i < X + rec.Width / 2; i++)
                            {
                                PointF pp = new PointF(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p4[0] == 0)
                                        {
                                            al2.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p4[0] == 255)
                                        {
                                            al2.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p4 += BPP;
                            }
                            p4 += offset;
                        }
                        edgeImage.UnlockBits(data4);
                    }
                    #endregion
                    #region//上下分割
                    if (dir == InOutDirection.In)
                    {
                        BitmapData data5 = edgeImage.LockBits(new Rectangle(0, 0, rec.Width, rec.Height / 2), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        int stride5 = data5.Stride;
                        byte* p5 = (byte*)data5.Scan0;
                        for (int i = X; i < X + rec.Width; i++)
                        {
                            p5 = (byte*)data5.Scan0 + (i - X) * BPP;
                            for (int j = Y; j < Y + rec.Height / 2; j++)
                            {
                                PointF pp = new PointF(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p5[0] == 0)
                                        {
                                            al3.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p5[0] == 255)
                                        {
                                            al3.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p5 += stride5;
                            }
                            p5 += BPP - stride5 * rec.Height / 2;
                        }
                        edgeImage.UnlockBits(data5);

                        BitmapData data6 = edgeImage.LockBits(new Rectangle(0, rec.Height / 2, rec.Width, rec.Height / 2), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        int stride6 = data6.Stride;
                        byte* p6 = (byte*)data6.Scan0 + stride6 * rec.Height / 2;
                        for (int i = X; i < X + rec.Width; i++)
                        {
                            p6 = (byte*)data6.Scan0 + stride6 * rec.Height / 2 + (i - X) * BPP;
                            for (int j = Y + rec.Height / 2 - 1; j > Y; j--)
                            {
                                PointF pp = new PointF(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p6[0] == 0)
                                        {
                                            al4.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p6[0] == 255)
                                        {
                                            al4.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p6 -= stride6;
                            }
                            p6 += BPP;
                        }
                        edgeImage.UnlockBits(data6);
                    }
                    else if (dir == InOutDirection.Out)
                    {
                        BitmapData data7 = edgeImage.LockBits(new Rectangle(0, 0, rec.Width, rec.Height / 2), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        int stride7 = data7.Stride;
                        byte* p7 = (byte*)data7.Scan0 + stride7 * rec.Height / 2;
                        for (int i = X; i < X + rec.Width; i++)
                        {
                            p7 = (byte*)data7.Scan0 + stride7 * rec.Height / 2 + (i - X) * BPP;
                            for (int j = Y + rec.Height / 2 - 1; j > Y; j--)
                            {
                                PointF pp = new PointF(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p7[0] == 0)
                                        {
                                            al3.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p7[0] == 255)
                                        {
                                            al3.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p7 -= stride7;
                            }
                            p7 += BPP;
                        }
                        edgeImage.UnlockBits(data7);

                        BitmapData data8 = edgeImage.LockBits(new Rectangle(0, rec.Height / 2, rec.Width, rec.Height / 2), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                        int stride8 = data8.Stride;
                        byte* p8 = (byte*)data8.Scan0;
                        for (int i = X; i < X + rec.Width; i++)
                        {
                            p8 = (byte*)data8.Scan0 + (i - X) * BPP;
                            for (int j = Y; j < Y + rec.Height / 2; j++)
                            {
                                PointF pp = new PointF(i, j);
                                if (rgn.IsVisible(pp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p8[0] == 0)
                                        {
                                            al4.Add(pp);
                                            break;
                                        }
                                    }
                                    else
                                    {
                                        if (p8[0] == 255)
                                        {
                                            al4.Add(pp);
                                            break;
                                        }
                                    }
                                }
                                p8 += stride8;
                            }
                            p8 += BPP - stride8 * rec.Height / 2;
                        }
                        edgeImage.UnlockBits(data8);
                    }
                    #endregion
                }
                edgeImage.Dispose();
                PointF[] ps1 = al1.ToArray();
                PointF[] ps2 = al2.ToArray();
                PointF[] ps3 = al3.ToArray();
                PointF[] ps4 = al4.ToArray();
                List<PointF> al = new List<PointF>();
                for (int i = 0; i < ps1.Length; i++)
                {
                    al.Add(ps1[i]);
                }
                for (int i = 0; i < ps2.Length; i++)
                {
                    if (!al.Contains(ps2[i]))
                        al.Add(ps2[i]);
                }
                for (int i = 0; i < ps3.Length; i++)
                {
                    if (!al.Contains(ps3[i]))
                        al.Add(ps3[i]);
                }
                for (int i = 0; i < ps4.Length; i++)
                {
                    if (!al.Contains(ps4[i]))
                        al.Add(ps4[i]);
                }
                PointF[] ps = al.ToArray();
                rectPs = al;
                GetRectangle4P(ps, out pp1, out pp2, out pp3, out pp4);

                if (pp1 != PointF.Empty || pp2 != PointF.Empty || pp3 != PointF.Empty || pp4 != PointF.Empty)
                    flag = true;
            }
            return flag;
        }
        /// <summary>
        /// 从检测到的点集中选取矩形的四点
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <param name="p4"></param>
        public static void GetRectangle4P(PointF[] ps, out PointF p1, out PointF p2, out PointF p3, out PointF p4)
        {
            p1 = PointF.Empty;
            p2 = PointF.Empty;
            p3 = PointF.Empty;
            p4 = PointF.Empty;
            PointF[] pp = Get4PFromRectangle(ps);
            if (pp == null)
                return;
            p1 = pp[0];
            p2 = pp[1];
            p3 = pp[2];
            p4 = pp[3];
            /*
            int vcount = ps.Length;
            if (vcount < 4)
            {
                throw new KellImageProcessException("检测到的点数少于四点，无法取矩形！");
            }
            PointF p = ps[0];
            System.Collections.ArrayList al = new System.Collections.ArrayList();
            for (int i = 1; i < vcount; i++) // 寻找凸顶点
            {
                if (ps[i].Y < p.Y || (ps[i].Y == p.Y && ps[i].X < p.X))
                {
                    p = ps[i];
                    al.Add(i);
                }
            }
            //MessageBox.Show("共有" + al.Count.ToString() + "个拐点。");
            if (al.Count > 3)
            {
                p1 = new Point2(ps[(int)al[0]].X, ps[(int)al[0]].Y);
                p2 = new Point2(ps[(int)al[1]].X, ps[(int)al[1]].Y);
                p3 = new Point2(ps[(int)al[2]].X, ps[(int)al[2]].Y);
                p4 = new Point2(ps[(int)al[3]].X, ps[(int)al[3]].Y);
            }
            */
        }

        /// <summary>
        /// 求平面向量的点积，并可以根据返回值判断向量的夹角是锐是直还是钝
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p0"></param>
        /// <returns></returns>
        public static double Dotmultiply(PointF p1, PointF p2, PointF p0)
        {
            return ((p1.X - p0.X) * (p2.X - p0.X) + (p1.Y - p0.Y) * (p2.Y - p0.Y));
        }

        /// <summary>
        /// 根据矩形点集获取四个顶点（驻点）
        /// </summary>
        /// <param name="rectangle"></param>
        /// <returns></returns>
        public static PointF[] Get4PFromRectangle(PointF[] rectangle)
        {
            PointF[] fourPs = null;
            if (rectangle.Length < 4)
            {
                throw new KellImageProcessException("检测点少于四个，无法取矩形！");
            }
            System.Collections.ArrayList al = new System.Collections.ArrayList();
            for (int i = 0; i < rectangle.Length - 3; i += 3)
            {
                //PointF[] deris = new PointF[3];
                //deris[0] = rectangle[i];
                //deris[1] = rectangle[i + 1];
                //deris[2] = rectangle[i + 2];
                //if (Derivative2(deris) != 0)//找出拐点(二阶导数不为零即为拐点)
                PointF a = new PointF(rectangle[i].X, rectangle[i].Y);
                PointF b = new PointF(rectangle[i].X, rectangle[i].Y);
                PointF c = new PointF(rectangle[i].X, rectangle[i].Y);
                if (Math.Abs(Dotmultiply(a, c, b)) < Common.Eps) // 说明b点是直角拐角处
                {
                    al.Add(b);
                }
            }
            //PointF[] deris1 = new PointF[3];
            //deris1[0] = rectangle[rectangle.Length - 1];
            //deris1[1] = rectangle[0];
            //deris1[2] = rectangle[1];
            //if (Derivative2(deris1) != 0)//找出rectangle的首部，看是否为拐点。
            PointF a1 = new PointF(rectangle[rectangle.Length - 1].X, rectangle[rectangle.Length - 1].Y);
            PointF b1 = new PointF(rectangle[0].X, rectangle[0].Y);
            PointF c1 = new PointF(rectangle[1].X, rectangle[1].Y);
            if (Math.Abs(Dotmultiply(a1, c1, b1)) < Common.Eps) // 说明b1点是直角拐角处
            {
                al.Add(b1);
            }
            //PointF[] deris2 = new PointF[3];
            //deris2[0] = rectangle[rectangle.Length - 2];
            //deris2[1] = rectangle[rectangle.Length - 1];
            //deris2[2] = rectangle[0];
            //if (Derivative2(deris2) != 0)//找出rectangle的尾部，看是否为拐点。
            PointF a2 = new PointF(rectangle[rectangle.Length - 2].X, rectangle[rectangle.Length - 2].Y);
            PointF b2 = new PointF(rectangle[rectangle.Length - 1].X, rectangle[rectangle.Length - 1].Y);
            PointF c2 = new PointF(rectangle[0].X, rectangle[0].Y);
            if (Math.Abs(Dotmultiply(a2, c2, b2)) < Common.Eps) // 说明b2点是直角拐角处
            {
                al.Add(b2);
            }
            if (al.Count < 4)
            {
                if (al.Count == 3)//如果只有三点，就智能地添加第四点
                {
                    al.Add(Rect4thP((PointF)al[0], (PointF)al[1], (PointF)al[2]));
                }
            }
            else if (al.Count > 4)
            {
                PointF firstP = (PointF)al[0];
                int second = al.Count / 4;
                PointF secondP = (PointF)al[second];
                int third = al.Count * 2 / 4;
                PointF thirdP = (PointF)al[third];
                int fourth = al.Count * 3 / 4;
                PointF fourthP = (PointF)al[fourth];
                al.Clear();
                al.Add(firstP);
                al.Add(secondP);
                al.Add(thirdP);
                al.Add(fourthP);
            }
            if (al.Count < 4)
            {
                throw new KellImageProcessException("拐点少于四个，无法取矩形！");
            }
            fourPs = new PointF[4];
            for (int i = 0; i < 4; i++)
            {
                fourPs[i] = (PointF)al[i];
            }
            return fourPs;
        }

        /// <summary>
        ///  判断象限
        /// </summary>
        /// <param name="a1">圆心X坐标</param>
        /// <param name="b1">圆心Y坐标</param>
        /// <param name="a2">任意X坐标</param>
        /// <param name="b2">任意Y坐标</param>
        /// <param name="rand">半径</param>
        /// <returns>返回任意点在圆中的角度值</returns>
        public float PanDuanXiangXian(float a1, float b1, float a2, float b2, double rand)
        {

            float a, b, c, d, e;
            //a 为余玄值。
            a = Convert.ToSingle((a2 - a1) / rand);
            //b为正玄值。
            b = Convert.ToSingle((b1 - b2) / rand);
            c = 0;
            d = 0;
            e = 0;

            if (a > 0)
            {
                if (b > 0)
                {
                    //第一象限 
                    c = Convert.ToSingle(Math.Acos(a));
                    d = 360 - c * 180 / (float)Math.PI;

                }
                else
                {
                    //第四象限
                    e = Convert.ToSingle(Math.Asin(b));
                    d = -e * 180 / (float)Math.PI;
                }
            }
            if (a < 0)
            {
                if (b > 0)
                {
                    //第二象限  
                    c = Convert.ToSingle(Math.Acos(-a));
                    d = 180 + c * 180 / (float)Math.PI;
                }
                else
                {
                    //第三象限   
                    e = Convert.ToSingle(Math.Asin(-b));
                    d = 180 - e * 180 / (float)Math.PI;
                }
            }
            return d;
        }

        /// <summary>
        ///  获取任意点somePoint相对于指定点center的方向角
        /// </summary>
        /// <param name="center">指定点</param>
        /// <param name="somePoint">任意点</param>
        /// <returns>返回圆上的任意点在圆中的角度值，逆时针方向为正</returns>
        public float GetFangXiangJiao(PointF center, PointF somePoint)
        {
            float a1 = center.X, b1 = center.Y, a2 = somePoint.X, b2 = somePoint.Y;
            float a, b, c, d, e;
            //a 为余玄值。
            a = Convert.ToSingle(a2 - a1);
            //b为正玄值。
            b = Convert.ToSingle(b1 - b2);
            c = 0;
            d = 0;
            e = 0;

            if (a > 0)
            {
                if (b > 0)
                {
                    //第一象限 
                    c = Convert.ToSingle(Math.Acos(a));
                    d = 360 - c * 180 / (float)Math.PI;
                }
                else
                {
                    //第四象限
                    e = Convert.ToSingle(Math.Asin(b));
                    d = -e * 180 / (float)Math.PI;
                }
            }
            if (a < 0)
            {
                if (b > 0)
                {
                    //第二象限  
                    c = Convert.ToSingle(Math.Acos(-a));
                    d = 180 + c * 180 / (float)Math.PI;
                }
                else
                {
                    //第三象限   
                    e = Convert.ToSingle(Math.Asin(-b));
                    d = 180 - e * 180 / (float)Math.PI;
                }
            }
            return d;
        }

        private void GetArcPara(PointF center, PointF startP, PointF middleP, PointF endP, out float startA, out float sweepA)
        {
            double Rand;
            PointF cc;
            GetCircle(startP, middleP, endP, out cc, out Rand);
            float X1 = center.X;
            float Y1 = center.Y;
            float start, mid, end, sweep;
            float x1, x2, x3, y1, y2, y3;
            x1 = Convert.ToSingle(startP.X);
            x2 = Convert.ToSingle(middleP.X);
            x3 = Convert.ToSingle(endP.X);

            y1 = Convert.ToSingle(startP.Y);
            y2 = Convert.ToSingle(middleP.Y);
            y3 = Convert.ToSingle(endP.Y);

            start = PanDuanXiangXian(X1, Y1, x1, y1, Rand);
            mid = PanDuanXiangXian(X1, Y1, x2, y2, Rand);
            end = PanDuanXiangXian(X1, Y1, x3, y3, Rand);
            double leng = 0;
            double squar = 0;
            sweep = start - end;
            if (start < mid)
            {
                if (end < mid)
                {
                    if (start < end)
                    {
                        leng = 2 * Math.PI * Rand * (-sweep) / 360;
                        squar = Math.PI * Rand * Rand * (-sweep) / 360;
                        sweep = -(360 + sweep);// 是负角。 
                    }
                    else
                    {
                        leng = 2 * Math.PI * Rand * sweep / 360;
                        squar = Math.PI * Rand * Rand * sweep / 360;
                        sweep = 360 - sweep;//sweep 是正角。
                    }
                }
                else
                {
                    sweep = -sweep;//sweep是负角
                    leng = 2 * Math.PI * Rand * sweep / 360;
                    squar = Math.PI * Rand * Rand * sweep / 360;

                }
            }
            else
            {
                if (end > mid)
                {
                    if (start < end)
                    {
                        leng = 2 * Math.PI * Rand * (-sweep) / 360;
                        squar = Math.PI * Rand * Rand * (-sweep) / 360;
                        sweep = -(360 + sweep);//sweep是负角。
                    }
                    else
                    {
                        leng = 2 * Math.PI * Rand * sweep / 360;
                        squar = Math.PI * Rand * Rand * sweep / 360;
                        sweep = 360 - sweep;// sweep是正角。
                    }
                }
                else
                {
                    leng = 2 * Math.PI * Rand * sweep / 360;
                    squar = Math.PI * Rand * Rand * sweep / 360;
                    sweep = -sweep;//sweep 是正角
                }
            }
            sweepA = sweep;
            startA = start;
        }
        /// <summary>
        /// 返回点p以点o为圆心逆时针旋转alpha(单位：弧度)后所在的位置
        /// </summary>
        /// <param name="o"></param>
        /// <param name="alpha"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public static PointF Rotate(PointF o, double alpha, PointF p)
        {
            PointF tp = PointF.Empty;
            p.X -= o.X;
            p.Y -= o.Y;
            tp.X = p.X * (float)Math.Cos(alpha) - p.Y * (float)Math.Sin(alpha) + o.X;
            tp.Y = p.Y * (float)Math.Cos(alpha) + p.X * (float)Math.Sin(alpha) + o.Y;
            return tp;
        }
        /// <summary>
        /// 根据弧的参数，获取弧上的标准三点（即：始点p1、中点p2、终点p3）
        /// </summary>
        /// <param name="rect"></param>
        /// <param name="startA"></param>
        /// <param name="sweepA"></param>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        public static void GetArc3P(RectangleF rect, float startA, float sweepA, out PointF p1, out PointF p2, out PointF p3)
        {
            double start = 2 * Math.PI - startA * Math.PI / 180;
            double sweep = sweepA * Math.PI / 180;
            double end = start - sweep;
            PointF arcCenter = new PointF(rect.X + rect.Width / 2, (rect.Y + rect.Height / 2));
            //弧的半径默认为 rect.Width / 2
            p1 = new PointF((float)(arcCenter.X + Math.Cos(start) * rect.Width / 2), (float)(arcCenter.Y + Math.Sin(start) * rect.Width / 2));
            p3 = new PointF((float)(arcCenter.X + Math.Cos(end) * rect.Width / 2), (float)(arcCenter.Y + Math.Sin(end) * rect.Width / 2));
            PointF pp1 = new PointF(p1.X, p1.Y);
            PointF centerP = Rotate(arcCenter, -(sweep / 2), pp1);//取半角（sweep / 2），逆时针故取负(-)
            p2 = centerP;
        }
        /// <summary>
        /// 多点定弧
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="rect"></param>
        /// <param name="startA"></param>
        /// <param name="sweepA"></param>
        public static void GetArc(Point[] ps, out Rectangle rect, out float startA, out float sweepA)
        {
            rect = Rectangle.Empty;
            startA = 0;
            sweepA = 0;
            if (ps.Length < 3)
            {
                throw new KellImageProcessException("检测到的点数少于三点，无法取弧！");
            }
            Point p;
            double ar;
            Point ps1 = ps[0];
            Point ps2 = ps[ps.Length / 2];
            Point ps3 = ps[ps.Length - 1];
            GetCircle(ps1, ps2, ps3, out p, out ar);
            rect = new Rectangle((int)(p.X - ar), (int)(p.Y - ar), (int)(2 * ar), (int)(2 * ar));
            double a1, a2, mida, mina, maxa;
            a1 = Math.Atan2(p.Y - ps1.Y, ps1.X - p.X);
            a2 = Math.Atan2(p.Y - ps3.Y, ps3.X - p.X);
            mida = Math.Atan2(p.Y - ps2.Y, ps2.X - p.X);
            if (mida < 0)
                mida += 2 * Math.PI;
            if (a1 < 0)
                a1 += 2 * Math.PI;
            if (a2 < 0)
                a2 += 2 * Math.PI;
            mina = a1 < a2 ? a1 : a2;
            maxa = a1 > a2 ? a1 : a2;
            startA = (float)(maxa * 180 / Math.PI);
            if (mida < mina || mida > maxa)//关键所在
                startA = (float)(mina * 180 / Math.PI);
            startA = 360 - startA;
            sweepA = (float)GetAnyJJ(ps1, ps2, ps3);
        }
        /// <summary>
        /// 多点定弧
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="rect"></param>
        /// <param name="startA"></param>
        /// <param name="sweepA"></param>
        public static void GetArc(PointF[] ps, out RectangleF rect, out float startA, out float sweepA)
        {
            rect = RectangleF.Empty;
            startA = 0;
            sweepA = 0;
            if (ps.Length < 3)
            {
                throw new KellImageProcessException("检测到的点数少于三点，无法取弧！");
            }
            PointF p;
            double ar;
            PointF ps1 = ps[0];
            PointF ps2 = ps[ps.Length / 2];
            PointF ps3 = ps[ps.Length - 1];
            GetCircle(ps1, ps2, ps3, out p, out ar);
            rect = new RectangleF((float)(p.X - ar), (float)(p.Y - ar), (float)(2 * ar), (float)(2 * ar));
            double a1, a2, mida, mina, maxa;
            a1 = Math.Atan2(p.Y - ps1.Y, ps1.X - p.X);
            a2 = Math.Atan2(p.Y - ps3.Y, ps3.X - p.X);
            mida = Math.Atan2(p.Y - ps2.Y, ps2.X - p.X);
            if (mida < 0)
                mida += 2 * Math.PI;
            if (a1 < 0)
                a1 += 2 * Math.PI;
            if (a2 < 0)
                a2 += 2 * Math.PI;
            mina = a1 < a2 ? a1 : a2;
            maxa = a1 > a2 ? a1 : a2;
            startA = (float)(maxa * 180 / Math.PI);
            if (mida < mina || mida > maxa)//关键所在
                startA = (float)(mina * 180 / Math.PI);
            startA = 360 - startA;
            sweepA = (float)GetAnyJJ(ps1, ps2, ps3);
        }
        /// <summary>
        /// 三点定圆
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <param name="cc"></param>
        /// <param name="cr"></param>
        public static void GetCircle(Point p1, Point p2, Point p3, out Point cc, out double cr)
        {
            cc = Point.Empty;
            cr = 0;
            float a1 = p1.X * (p2.Y - p3.Y) + p2.X * (p3.Y - p1.Y) + p3.X * (p1.Y - p2.Y);
            float b1 = (p1.X * p1.X + p1.Y * p1.Y - p2.X * p2.X - p2.Y * p2.Y) / 2;
            float c1 = (p1.X * p1.X + p1.Y * p1.Y - p3.X * p3.X - p3.Y * p3.Y) / 2;
            float d = b1 * (p1.Y - p3.Y) - c1 * (p1.Y - p2.Y);
            float e = c1 * (p1.X - p2.X) - b1 * (p1.X - p3.X);
            if (a1 != 0)//当a1 = 0 时三点共线，无法成圆。
            {
                cc = new Point((int)(d / a1), (int)(e / a1));
                cr = Math.Sqrt((cc.X - p1.X) * (cc.X - p1.X) + (cc.Y - p1.Y) * (cc.Y - p1.Y));
            }
        }
        /// <summary>
        /// 三点定圆
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <param name="cc"></param>
        /// <param name="cr"></param>
        public static void GetCircle(PointF p1, PointF p2, PointF p3, out PointF cc, out double cr)
        {
            cc = PointF.Empty;
            cr = 0;
            float a1 = p1.X * (p2.Y - p3.Y) + p2.X * (p3.Y - p1.Y) + p3.X * (p1.Y - p2.Y);
            float b1 = (p1.X * p1.X + p1.Y * p1.Y - p2.X * p2.X - p2.Y * p2.Y) / 2;
            float c1 = (p1.X * p1.X + p1.Y * p1.Y - p3.X * p3.X - p3.Y * p3.Y) / 2;
            float d = b1 * (p1.Y - p3.Y) - c1 * (p1.Y - p2.Y);
            float e = c1 * (p1.X - p2.X) - b1 * (p1.X - p3.X);
            if (a1 != 0)//当a1 = 0 时三点共线，无法成圆。
            {
                cc = new PointF(d / a1, e / a1);
                cr = Math.Sqrt((cc.X - p1.X) * (cc.X - p1.X) + (cc.Y - p1.Y) * (cc.Y - p1.Y));
            }
        }
        /// <summary>
        /// 多点拟合圆
        /// </summary>
        /// <param name="pts">点集,点集数大于三</param>
        /// <param name="cen"></param>
        /// <param name="r"></param>
        /// <returns></returns>
        public static void GetCircle(List<Point> pts, out PointF cen, out double r)
        {
            cen = PointF.Empty;
            r = 0;
            if (pts.Count < 3)
            {
                throw new KellImageProcessException("点数少于3，无法拟合圆！");
            }
            try
            {
                double X1 = 0;
                double Y1 = 0;
                double X2 = 0;
                double Y2 = 0;
                double X3 = 0;
                double Y3 = 0;
                double X1Y1 = 0;
                double X1Y2 = 0;
                double X2Y1 = 0;
                int m_nNum = pts.Count;
                for (int i = 0; i < m_nNum; i++)
                {
                    X1 = X1 + pts[i].X;
                    Y1 = Y1 + pts[i].Y;
                    X2 = X2 + pts[i].X * pts[i].X;
                    Y2 = Y2 + pts[i].Y * pts[i].Y;
                    X3 = X3 + pts[i].X * pts[i].X * pts[i].X;
                    Y3 = Y3 + pts[i].Y * pts[i].Y * pts[i].Y;
                    X1Y1 = X1Y1 + pts[i].X * pts[i].Y;
                    X1Y2 = X1Y2 + pts[i].X * pts[i].Y * pts[i].Y;
                    X2Y1 = X2Y1 + pts[i].X * pts[i].X * pts[i].Y;
                }
                double C, D, E, G, H, N;
                double a, b, c;
                N = m_nNum;
                C = N * X2 - X1 * X1;
                D = N * X1Y1 - X1 * Y1;
                E = N * X3 + N * X1Y2 - (X2 + Y2) * X1;
                G = N * Y2 - Y1 * Y1;
                H = N * X2Y1 + N * Y3 - (X2 + Y2) * Y1;
                a = (H * D - E * G) / (C * G - D * D);
                b = (H * C - E * D) / (D * D - G * C);
                c = -(a * X1 + b * Y1 + X2 + Y2) / N;
                double A, B, R;
                A = a / (-2);
                B = b / (-2);
                R = Math.Sqrt(a * a + b * b - 4 * c) / 2;
                PointF po = new PointF();
                po.X = (float)A;
                po.Y = (float)B;
                cen = po;
                r = R;
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }
        /// <summary>
        /// 多点拟合圆
        /// </summary>
        /// <param name="pts">点集,点集数大于三</param>
        /// <param name="cen"></param>
        /// <param name="r"></param>
        /// <returns></returns>
        public static void GetCircle(List<PointF> pts, out PointF cen, out double r)
        {
            cen = PointF.Empty;
            r = 0;
            if (pts.Count < 3)
            {
                throw new KellImageProcessException("点数少于3，无法拟合圆！");
            }
            try
            {
                double X1 = 0;
                double Y1 = 0;
                double X2 = 0;
                double Y2 = 0;
                double X3 = 0;
                double Y3 = 0;
                double X1Y1 = 0;
                double X1Y2 = 0;
                double X2Y1 = 0;
                int m_nNum = pts.Count;
                for (int i = 0; i < m_nNum; i++)
                {
                    X1 = X1 + pts[i].X;
                    Y1 = Y1 + pts[i].Y;
                    X2 = X2 + pts[i].X * pts[i].X;
                    Y2 = Y2 + pts[i].Y * pts[i].Y;
                    X3 = X3 + pts[i].X * pts[i].X * pts[i].X;
                    Y3 = Y3 + pts[i].Y * pts[i].Y * pts[i].Y;
                    X1Y1 = X1Y1 + pts[i].X * pts[i].Y;
                    X1Y2 = X1Y2 + pts[i].X * pts[i].Y * pts[i].Y;
                    X2Y1 = X2Y1 + pts[i].X * pts[i].X * pts[i].Y;
                }
                double C, D, E, G, H, N;
                double a, b, c;
                N = m_nNum;
                C = N * X2 - X1 * X1;
                D = N * X1Y1 - X1 * Y1;
                E = N * X3 + N * X1Y2 - (X2 + Y2) * X1;
                G = N * Y2 - Y1 * Y1;
                H = N * X2Y1 + N * Y3 - (X2 + Y2) * Y1;
                a = (H * D - E * G) / (C * G - D * D);
                b = (H * C - E * D) / (D * D - G * C);
                c = -(a * X1 + b * Y1 + X2 + Y2) / N;
                double A, B, R;
                A = a / (-2);
                B = b / (-2);
                R = Math.Sqrt(a * a + b * b - 4 * c) / 2;
                PointF po = new PointF();
                po.X = (float)A;
                po.Y = (float)B;
                cen = po;
                r = R;
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }

        /// <summary>
        /// 点数组转换
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static List<PointF> PointList2PointFList(List<Point> ps)
        {
            List<PointF> retPs = new List<PointF>();
            for (int i = 0; i < ps.Count; i++)
            {
                retPs.Add(new PointF(ps[i].X, ps[i].Y));
            }
            return retPs;
        }

        /// <summary>
        /// 点数组转换
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static List<Point> PointFList2PointList(List<PointF> ps)
        {
            List<Point> retPs = new List<Point>();
            for (int i = 0; i < ps.Count; i++)
            {
                retPs.Add(Point.Round(ps[i]));
            }
            return retPs;
        }

        /// <summary>
        /// 三点决定的角，有顺序性
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <returns></returns>
        public static float GetAnyJJ(PointF p1, PointF p2, PointF p3)
        {
            float jj = 0;
            PointF c;
            //GetCircleCenter(p1, p2, p3, out c);
            double cr;
            GetCircle(p1, p2, p3, out c, out cr);
            //关键所在：
            if (c != null)
            {
                double j1, j2, midj, minj, maxj;
                j1 = Math.Atan2(c.Y - p1.Y, p1.X - c.X);
                j2 = Math.Atan2(c.Y - p3.Y, p3.X - c.X);
                midj = Math.Atan2(c.Y - p2.Y, p2.X - c.X);
                if (midj < 0)
                    midj += 2 * Math.PI;
                if (j1 < 0)
                    j1 += 2 * Math.PI;
                if (j2 < 0)
                    j2 += 2 * Math.PI;
                minj = j1 < j2 ? j1 : j2;
                maxj = j1 > j2 ? j1 : j2;
                jj = (float)Math.Abs((maxj - minj) * 180 / Math.PI);
                if (midj < minj || midj > maxj)
                    jj = 360 - jj;
            }
            return jj;
        }
        /// <summary>
        /// 由两点获取方向角
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static double GetAngle(Point p1, Point p2)
        {
            double angle = Math.Atan2(p1.Y - p2.Y, p2.X - p1.X);
            if (angle < 0)
                angle += 2 * Math.PI;
            return angle;
        }
        /// <summary>
        /// 由两点获取方向角
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static double GetAngle(PointF p1, PointF p2)
        {
            double angle = Math.Atan2(p1.Y - p2.Y, p2.X - p1.X);
            if (angle < 0)
                angle += 2 * Math.PI;
            return angle;
        }
        /// <summary>
        /// 获取弧的始点、中点、终点
        /// </summary>
        /// <param name="PointList1">左半区域检测到的点集</param>
        /// <param name="PointList2">右半区域检测到的点集</param>
        /// <param name="StartP"></param>
        /// <param name="MiddleP"></param>
        /// <param name="EndP"></param>
        public void GetArc3P1(PointF[] PointList1, PointF[] PointList2, out PointF StartP, out PointF MiddleP, out PointF EndP)
        {
            StartP = PointF.Empty;
            MiddleP = PointF.Empty;
            EndP = PointF.Empty;
            if (PointList1.Length == 0)
            {//ok
                if (PointList2.Length > 2)
                {
                    StartP = PointList2[0];
                    MiddleP = PointList2[PointList2.Length / 2];
                    EndP = PointList2[PointList2.Length - 1];
                }
                else
                {
                    throw new KellImageProcessException("达不到画弧要求");
                }
            }
            else if (PointList1.Length == 1)
            {
                //ok
                if (PointList2.Length == 2)
                {
                    StartP = PointList1[0];
                    if (PointList2[0].X > PointList2[1].X)
                    {
                        MiddleP = PointList2[1];
                        EndP = PointList2[0];
                    }
                    else if (PointList2[0].X == PointList2[1].X)
                    {
                        if (PointList2[1].Y > PointList1[0].Y)
                        {
                            MiddleP = PointList2[0];
                            EndP = PointList2[1];
                        }
                        else
                        {
                            MiddleP = PointList2[1];
                            EndP = PointList2[0];
                        }
                    }
                    else
                    {
                        MiddleP = PointList2[0];
                        EndP = PointList2[1];
                    }
                }
                //ok
                else if (PointList2.Length > 2)
                {
                    StartP = PointList1[0];
                    if (PointList2[0].X > PointList2[PointList2.Length - 1].X)
                    {
                        MiddleP = PointList2[PointList2.Length / 2];
                        EndP = PointList2[0];
                    }
                    else
                    {
                        MiddleP = PointList2[PointList2.Length / 2];
                        EndP = PointList2[PointList2.Length - 1];
                    }
                }
                else
                {
                    throw new KellImageProcessException("达不到画弧要求");
                }
            }
            else if (PointList1.Length == 2)
            {
                //ok
                if (PointList2.Length == 1)
                {
                    EndP = PointList2[0];
                    if (PointList1[0].X > PointList1[1].X)
                    {
                        StartP = PointList1[1];
                        MiddleP = PointList1[0];
                    }
                    else if (PointList1[0].X == PointList1[1].X)
                    {
                        if (PointList1[1].Y > PointList2[0].Y)
                        {
                            MiddleP = PointList1[0];
                            EndP = PointList1[1];
                        }
                        else
                        {
                            MiddleP = PointList1[1];
                            EndP = PointList1[0];
                        }
                    }
                    else
                    {
                        StartP = PointList1[0];
                        MiddleP = PointList1[1];
                    }
                }//
                else if (PointList2.Length >= 2)
                {
                    if (PointList1[0].X > PointList1[1].X)
                    {
                        StartP = PointList1[1];
                        MiddleP = PointList2[PointList2.Length / 2];
                        EndP = PointList2[PointList2.Length - 1];
                    }
                    else if (PointList1[0].X == PointList1[1].X)
                    {
                        if (PointList1[1].Y > PointList2[0].Y && PointList1[1].Y < PointList2[PointList2.Length - 1].Y)
                        {
                            StartP = PointList1[1];
                            MiddleP = PointList2[PointList2.Length / 2];
                            EndP = PointList2[PointList2.Length - 1];
                        }
                        else
                        {
                            StartP = PointList1[0];
                            MiddleP = PointList2[PointList2.Length / 2];
                            EndP = PointList2[0];
                        }
                    }
                    else
                    {
                        StartP = PointList1[0];
                        MiddleP = PointList2[PointList2.Length / 2];
                        EndP = PointList2[0];
                    }
                }
                else
                {
                    throw new KellImageProcessException("达不到画弧要求");
                }
            }
            else if (PointList1.Length > 2)
            {
                if (PointList2.Length == 0)
                {
                    StartP = PointList1[PointList1.Length - 1];
                    MiddleP = PointList1[PointList1.Length / 2];
                    EndP = PointList1[0];
                }
                else if (PointList2.Length == 1)
                {
                    EndP = PointList2[0];
                    if (PointList1[0].X > PointList1[PointList1.Length - 1].X)
                    {
                        StartP = PointList1[PointList1.Length - 1];
                        MiddleP = PointList1[PointList1.Length / 2];
                    }
                    else
                    {
                        StartP = PointList1[0];
                        MiddleP = PointList1[PointList1.Length / 2];
                    }
                }
                else if (PointList2.Length >= 2)
                {
                    if (PointList1[0].X > PointList1[PointList1.Length - 1].X)
                    {
                        StartP = PointList1[PointList1.Length - 1];
                        MiddleP = PointList1[PointList1.Length / 2];
                        EndP = PointList2[PointList2.Length - 1];
                    }
                    else
                    {
                        StartP = PointList1[0];
                        MiddleP = PointList1[PointList1.Length / 2];
                        EndP = PointList2[0];
                    }
                }
            }
        }

        /// <summary>
        /// 从给定的位图中利用自适应算法获取最清晰的(不一定是最亮的)点集，注意：minRect>=step
        /// Bright效果比Gray好
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="step">步长>=1</param>
        /// <param name="minRect">检测框大小>=2</param>
        /// <returns></returns>
        public static List<Point> GetTheMaxTrenchantPsByBright(Bitmap bmp, int step, int minRect)
        {
            List<Point> tmpPs = new List<Point>();
            //const int step = 2;
            //const int minRect = 2;
            step = step < 1 ? 1 : step;//防止步长为0或反向（为负）
            step = step > minRect ? minRect : step;//不能大于检测框大小，以防遗漏待检区域
            minRect = minRect < 2 ? 2 : minRect;
            if (bmp != null && bmp.Width > minRect && bmp.Height > minRect)
            {
                byte autoFit = GetAutoFitTrenchantByBright(bmp);//获取清晰度自适应阀值
                for (int j = 0; j < bmp.Height - minRect; j += step)
                {
                    for (int i = 0; i < bmp.Width - minRect; i += step)
                    {
                        Rectangle rect = new Rectangle(i, j, minRect, minRect);
                        byte tmp = (byte)(GetTrenchantByBrightness(bmp, rect) * 255);
                        if (tmp > autoFit)
                        {
                            tmpPs.Add(new Point(i, j));
                        }
                    }
                }
            }
            return tmpPs;
        }

        /// <summary>
        /// 从给定的位图中利用自适应算法获取最清晰的(不一定是最亮的)点集，注意：minRect>=step
        /// Gray效果没有Bright好
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="step">步长>=1</param>
        /// <param name="minRect">检测框大小>=2</param>
        /// <returns></returns>
        public static List<Point> GetTheMaxTrenchantPsByGray(Bitmap bmp, int step, int minRect)
        {
            List<Point> tmpPs = new List<Point>();
            //const int step = 2;
            //const int minRect = 2;
            step = step < 1 ? 1 : step;//防止步长为0或反向（为负）
            step = step > minRect ? minRect : step;//不能大于检测框大小，以防遗漏待检区域
            minRect = minRect < 2 ? 2 : minRect;
            if (bmp != null && bmp.Width > minRect && bmp.Height > minRect)
            {
                byte autoFit = GetAutoFitTrenchantByGray(bmp);//获取清晰度自适应阀值
                for (int j = 0; j < bmp.Height - minRect; j += step)
                {
                    for (int i = 0; i < bmp.Width - minRect; i += step)
                    {
                        Rectangle rect = new Rectangle(i, j, minRect, minRect);
                        byte tmp = GetTrenchantByGray(bmp, rect);// * 255;
                        if (tmp > autoFit)
                        {
                            tmpPs.Add(new Point(i, j));
                        }
                    }
                }
            }
            return tmpPs;
        }

        /// <summary>
        /// 基色分量的自适应灰度值，范围[0,255]
        /// </summary>
        /// <param name="RGB">R=2,G=1,B=0,其它值为Gray</param>
        /// <returns></returns>
        private byte GetAutoFitTrenchant(int RGB)
        {
            byte autoFit = 0;
            if (bmp != null)
            {
                // 图像灰度化
                Bitmap bmp1 = Gray((Bitmap)bmp.Clone());
                // 建立直方图，并获取灰度统计信息
                Histogram histogram = new Histogram(bmp1);
                int[] GrayLevel = null;
                switch (RGB)
                {
                    case 0:
                        GrayLevel = histogram.Blue.Value;
                        break;
                    case 1:
                        GrayLevel = histogram.Green.Value;
                        break;
                    case 2:
                        GrayLevel = histogram.Red.Value;
                        break;
                    default:
                        GrayLevel = histogram.Gray.Value;
                        break;
                }
                int peak1, peak2, valley;
                int peak1Index, peak2Index, valleyIndex;
                // 取双峰
                peak1 = peak2 = GrayLevel[0];
                peak1Index = peak2Index = 0;
                for (int i = 1; i < 256; i++)
                {
                    // 如果产生新的高峰，则将第一峰退居第二峰，新的高峰升为第一峰
                    if (GrayLevel[i] > peak1)
                    {
                        peak2 = peak1;
                        peak2Index = peak1Index;
                        peak1 = GrayLevel[i];
                        peak1Index = i;
                    }
                } // i
                // 判断两个峰值索引
                int max = peak1Index;
                int min = peak2Index;
                if (max < min)
                {
                    int t = max;
                    max = min;
                    min = t;
                }
                // 取峰谷
                valley = GrayLevel[min];
                valleyIndex = min;
                for (int i = min; i < max; i++)
                {
                    if (GrayLevel[i] < valley)
                    {
                        valley = GrayLevel[i];
                        valleyIndex = i;
                    }
                } // i
                // 找到谷值
                autoFit = (byte)valleyIndex;
            }

            return autoFit;
        }

        /// <summary>
        /// 亮度自适应值
        /// </summary>
        /// <returns></returns>
        private static byte GetAutoFitTrenchantByBright(Bitmap bmp)
        {
            byte autoFit = 0;
            if (bmp != null)
            {
                // 图像灰度化
                //Bitmap bmp1 = Gray((Bitmap)bmp.Clone());
                Bitmap bmp1 = (Bitmap)bmp.Clone();//不必灰度化
                // 建立直方图，并获取灰度统计信息
                Histogram histogram = new Histogram(bmp1);
                bmp1.Dispose();
                int[] GrayLevel = histogram.Bright.Value;
                int peak1, peak2, valley;
                int peak1Index, peak2Index, valleyIndex;
                // 取双峰
                peak1 = peak2 = GrayLevel[0];
                peak1Index = peak2Index = 0;
                for (int i = 1; i < 256; i++)
                {
                    // 如果产生新的高峰，则将第一峰退居第二峰，新的高峰升为第一峰
                    if (GrayLevel[i] > peak1)
                    {
                        peak2 = peak1;
                        peak2Index = peak1Index;
                        peak1 = GrayLevel[i];
                        peak1Index = i;
                    }
                } // i
                // 判断两个峰值索引
                int max = peak1Index;
                int min = peak2Index;
                if (max < min)
                {
                    int t = max;
                    max = min;
                    min = t;
                }
                // 取峰谷
                valley = GrayLevel[min];
                valleyIndex = min;
                for (int i = min; i < max; i++)
                {
                    if (GrayLevel[i] < valley)
                    {
                        valley = GrayLevel[i];
                        valleyIndex = i;
                    }
                } // i
                // 找到谷值
                autoFit = (byte)valleyIndex;
            }

            return autoFit;
        }

        /// <summary>
        /// 灰度自适应值
        /// </summary>
        /// <returns></returns>
        private static byte GetAutoFitTrenchantByGray(Bitmap bmp)
        {
            byte autoFit = 0;
            if (bmp != null)
            {
                // 图像灰度化
                //Bitmap bmp1 = Gray((Bitmap)bmp.Clone());
                Bitmap bmp1 = (Bitmap)bmp.Clone();//不必灰度化
                // 建立直方图，并获取灰度统计信息
                Histogram histogram = new Histogram(bmp1);
                bmp1.Dispose();
                int[] GrayLevel = histogram.Gray.Value;
                int peak1, peak2, valley;
                int peak1Index, peak2Index, valleyIndex;
                // 取双峰
                peak1 = peak2 = GrayLevel[0];
                peak1Index = peak2Index = 0;
                for (int i = 1; i < 256; i++)
                {
                    // 如果产生新的高峰，则将第一峰退居第二峰，新的高峰升为第一峰
                    if (GrayLevel[i] > peak1)
                    {
                        peak2 = peak1;
                        peak2Index = peak1Index;
                        peak1 = GrayLevel[i];
                        peak1Index = i;
                    }
                } // i
                // 判断两个峰值索引
                int max = peak1Index;
                int min = peak2Index;
                if (max < min)
                {
                    int t = max;
                    max = min;
                    min = t;
                }
                // 取峰谷
                valley = GrayLevel[min];
                valleyIndex = min;
                for (int i = min; i < max; i++)
                {
                    if (GrayLevel[i] < valley)
                    {
                        valley = GrayLevel[i];
                        valleyIndex = i;
                    }
                } // i
                // 找到谷值
                autoFit = (byte)valleyIndex;
            }

            return autoFit;
        }

        /// <summary>
        /// 从给定的位图数组中获取最清晰的一张
        /// Bright效果比Gray好
        /// </summary>
        /// <param name="bitMaps">位图数组</param>
        /// <param name="size">清晰度检测正方形区域的宽度</param>
        /// <returns></returns>
        public static Bitmap GetTheMaxTrenchantBitmapByBright(List<Bitmap> bitMaps, int size)
        {
            if (bitMaps != null && bitMaps.Count > 0)
            {
                float max = 0;
                int maxIndex = 0;
                if (size < 5)
                    size = 5;
                for (int i = 0; i < bitMaps.Count; i++)
                {
                    Bitmap bmp = bitMaps[i];
                    int width = bmp.Width;
                    int height = bmp.Height;
                    int minSize = width < height ? width : height;
                    if (size > minSize)
                        size = minSize;
                    int midX = (int)0.5 * width;
                    int midY = (int)0.5 * height;
                    int halfSize = (int)0.5 * size;
                    Rectangle rect = new Rectangle(midX - halfSize, midY - halfSize, size, size);
                    float tmp = GetTrenchantByBrightness(bmp, rect);
                    if (max < tmp)
                    {
                        max = tmp;
                        maxIndex = i;
                    }
                }
                return bitMaps[maxIndex];
            }
            return null;
        }

        /// <summary>
        /// 从给定的位图数组中获取最清晰的一张
        /// Gray效果没有Bright好
        /// </summary>
        /// <param name="bitMaps">位图数组</param>
        /// <param name="size">清晰度检测正方形区域的宽度</param>
        /// <returns></returns>
        public static Bitmap GetTheMaxTrenchantBitmapByGray(List<Bitmap> bitMaps, int size)
        {
            if (bitMaps != null && bitMaps.Count > 0)
            {
                byte max = 0;
                int maxIndex = 0;
                if (size < 5)
                    size = 5;
                for (int i = 0; i < bitMaps.Count; i++)
                {
                    Bitmap bmp = bitMaps[i];
                    int width = bmp.Width;
                    int height = bmp.Height;
                    int minSize = width < height ? width : height;
                    if (size > minSize)
                        size = minSize;
                    int midX = (int)0.5 * width;
                    int midY = (int)0.5 * height;
                    int halfSize = (int)0.5 * size;
                    Rectangle rect = new Rectangle(midX - halfSize, midY - halfSize, size, size);
                    byte tmp = GetTrenchantByGray(bmp, rect);
                    if (max < tmp)
                    {
                        max = tmp;
                        maxIndex = i;
                    }
                }
                return bitMaps[maxIndex];
            }
            return null;
        }

        /// <summary>
        /// 获取位图中指定区域的清晰度[0.0,1.0]
        /// </summary>
        /// <param name="rect"></param>
        /// <returns></returns>
        public float GetTrenchant(Rectangle rect)
        {
            float max = 0;
            if (rect.X >= 0 && rect.Y >= 0 && rect.Width + rect.X <= bmp.Width && rect.Height + rect.Y <= bmp.Height)
            {
                float maxx1 = GetOneContrast(rect);
                int width = rect.Width;
                int height = rect.Height;
                Rectangle rect1 = new Rectangle(rect.X + (int)(0.25 * width), rect.Y + (int)(0.25 * height), (int)(0.5 * width), (int)(0.5 * height));
                float maxx2 = GetOneContrast(rect1);
                //Rectangle rect1 = new Rectangle(rect.X + (int)(0.165 * width), rect.Y + (int)(0.165 * height), (int)(0.66 * width), (int)(0.66 * height));
                //double maxx2 = GetOneContrast(rect1);
                //Rectangle rect2 = new Rectangle(rect.X + (int)(0.33 * width), rect.Y + (int)(0.33 * height), (int)(0.33 * width), (int)(0.33 * height));
                //double maxx3 = GetOneContrast(rect2);
                //max = (maxx1 + maxx2 + maxx3) / 3;//取全区域和各种子区域的对比反差的平均值，以提高准确度
                max = (maxx1 + maxx2) / 2;//取全区域和各种子区域的对比反差的平均值，以提高准确度
            }
            return max;
        }

        /// <summary>
        /// 获取位图中指定区域的清晰度[0.0,1.0]，最新算法
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="rect">指定区域</param>
        /// <returns></returns>
        public static double GetTrenchant(Bitmap bmp, Rectangle rect)
        {
            double Trenchant = 0;
            if (bmp != null && rect.X >= 0 && rect.Y >= 0 && rect.Width + rect.X <= bmp.Width && rect.Height + rect.Y <= bmp.Height && rect.Width >= 2 && rect.Height >= 2)
            {
                double factor = Math.Sqrt(2) / 2;
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                byte R, G, B;
                byte[] ZS = new byte[3];
                byte[] S = new byte[3];
                byte[] YS = new byte[3];
                byte[] Z = new byte[3];
                byte[] Y = new byte[3];
                byte[] ZX = new byte[3];
                byte[] X = new byte[3];
                byte[] YX = new byte[3];
                double sum = 0;
                BitmapData data = bmp.LockBits(rect, ImageLockMode.ReadOnly, bmp.PixelFormat);
                unsafe
                {
                    byte* p = (byte*)data.Scan0;
                    int stride = data.Stride;
                    int offset = stride - BPP * rect.Width;
                    p += BPP;
                    p += stride;
                    for (int y = 1; y < rect.Height - 1; y++)
                    {
                        for (int x = 1; x < rect.Width - 1; x++)
                        {
                            R = p[2];
                            G = p[1];
                            B = p[0];
                            double d = Math.Sqrt(R * R + G * G + B * B);
                            if (rect.Width > 2 && rect.Height > 2)
                            {
                                for (int i = 0; i < 3; i++)
                                {
                                    ZS[i] = p[i - stride - BPP];
                                    S[i] = p[i - stride];
                                    YS[i] = p[i - stride + BPP];
                                    Z[i] = p[i - BPP];
                                    Y[i] = p[i + BPP];
                                    ZX[i] = p[i + stride - BPP];
                                    X[i] = p[i + stride];
                                    YX[i] = p[i + stride + BPP];
                                }
                                double dZS = Math.Sqrt(ZS[2] * ZS[2] + ZS[1] * ZS[1] + ZS[0] * ZS[0]);
                                double dS = Math.Sqrt(Z[2] * Z[2] + Z[1] * Z[1] + Z[0] * Z[0]);
                                double dYS = Math.Sqrt(YS[2] * YS[2] + YS[1] * YS[1] + YS[0] * YS[0]);
                                double dZ = Math.Sqrt(Z[2] * Z[2] + Z[1] * Z[1] + Z[0] * Z[0]);
                                double dY = Math.Sqrt(Y[2] * Y[2] + Y[1] * Y[1] + Y[0] * Y[0]);
                                double dZX = Math.Sqrt(ZX[2] * ZX[2] + ZX[1] * ZX[1] + ZX[0] * ZX[0]);
                                double dX = Math.Sqrt(X[2] * X[2] + X[1] * X[1] + X[0] * X[0]);
                                double dYX = Math.Sqrt(YX[2] * YX[2] + YX[1] * YX[1] + YX[0] * YX[0]);
                                sum += Math.Abs(dZS - d) * factor + Math.Abs(dS - d) + Math.Abs(dYS - d) * factor + Math.Abs(dZ - d) + Math.Abs(dY - d) + Math.Abs(dZX - d) * factor + Math.Abs(dX - d) + Math.Abs(dYX - d) * factor;
                            }
                            else if (rect.Width == 2 || rect.Height == 2)
                            {
                                for (int i = 0; i < 3; i++)
                                {
                                    ZS[i] = p[i - stride - BPP];
                                    S[i] = p[i - stride];
                                    Z[i] = p[i - BPP];
                                }
                                double dZS = Math.Sqrt(ZS[2] * ZS[2] + ZS[1] * ZS[1] + ZS[0] * ZS[0]);
                                double dS = Math.Sqrt(Z[2] * Z[2] + Z[1] * Z[1] + Z[0] * Z[0]);
                                double dZ = Math.Sqrt(Z[2] * Z[2] + Z[1] * Z[1] + Z[0] * Z[0]);
                                sum += Math.Abs(dZS - d) * factor + Math.Abs(dS - d) + Math.Abs(dZ - d);
                            }
                            p += BPP;
                        }  // x
                        p += offset;
                    } // y
                }
                bmp.UnlockBits(data);
                if (rect.Width > 2 && rect.Height > 2)
                {
                    Trenchant = sum / ((rect.Width - 2) * (rect.Height - 2));
                }
                else if (rect.Width == 2 || rect.Height == 2)
                {
                    Trenchant = sum / ((rect.Width - 1) * (rect.Height - 1) * 4);
                }
            }
                /*
            else if (bmp != null && rect.X >= 0 && rect.Y >= 0 && rect.Width + rect.X <= bmp.Width && rect.Height + rect.Y <= bmp.Height && (rect.Width == 2 && rect.Height == 2))
            {
                //A, B
                //C, D
                Color A, B, C, D;
                A = bmp.GetPixel(rect.X, rect.Y);
                B = bmp.GetPixel(rect.X + 1, rect.Y);
                C = bmp.GetPixel(rect.X, rect.Y + 1);
                D = bmp.GetPixel(rect.X + 1, rect.Y + 1);
                double d = Math.Sqrt(D.R * D.R + D.G * D.G + D.B * D.B);
                double AD = Math.Sqrt(A.R * A.R + A.G * A.G + A.B * A.B) - d;
                double BD = Math.Sqrt(B.R * B.R + B.G * B.G + B.B * B.B) - d;
                double CD = Math.Sqrt(C.R * C.R + C.G * C.G + C.B * C.B) - d;
                Trenchant = Math.Abs(AD) + Math.Abs(BD) + Math.Abs(CD);
            }*/
            double tr = Trenchant / 256;
            return tr > 1 ? 1 : tr;
            //return Trenchant / 256;
        }

        /// <summary>
        /// 获取位图中指定区域的清晰度(改进算法)[0.0,1.0]，较省时
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="rect">指定区域</param>
        /// <returns></returns>
        public static float GetTrenchantByBrightness(Bitmap bmp, Rectangle rect)
        {
            float max = 0;
            if (rect.X >= 0 && rect.Y >= 0 && rect.Width + rect.X <= bmp.Width && rect.Height + rect.Y <= bmp.Height)
            {
                max = GetOneContrastByBrightness(bmp, rect);
            }
            return max;
        }

        /// <summary>
        /// 获取位图中指定区域的清晰度(改进算法1)[0.0,1.0]，较耗时
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="rect">指定区域</param>
        /// <returns></returns>
        public static float GetTrenchantByBrightness1(Bitmap bmp, Rectangle rect)
        {
            float max = 0;
            if (rect.X >= 0 && rect.Y >= 0 && rect.Width + rect.X <= bmp.Width && rect.Height + rect.Y <= bmp.Height)
            {
                max = GetOneContrastByBrightness1(bmp, rect);
            }
            return max;
        }

        /// <summary>
        /// 获取位图中指定区域的清晰度(改进算法)
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="rect">指定区域</param>
        /// <returns></returns>
        public static byte GetTrenchantByGray(Bitmap bmp, Rectangle rect)
        {
            byte max = 0;
            if (rect.X >= 0 && rect.Y >= 0 && rect.Width + rect.X <= bmp.Width && rect.Height + rect.Y <= bmp.Height)
            {
                max = GetOneContrastByGray(bmp, rect);
            }
            return max;
        }

        /// <summary>
        /// 在水平和垂直两个方向上取亮度的对比反差即清晰度（应该用减法而非除法）的最大值[0.0,1.0]
        /// </summary>
        /// <param name="rect"></param>
        /// <returns></returns>
        private float GetOneContrast(Rectangle rect)
        {
            //byte contrast = 0;
            float max1 = 0, max2 = 0, max11 = 0, max12 = 0, max21 = 0, max22 = 0;
            float min11 = 1F, min12 = 1F, min21 = 1F, min22 = 1F;
            byte R, G, B;//, gray = 128;
            float bright, bright1, bright2, bright3;
            byte R1, G1, B1;//, gray1 = 128;
            byte R2, G2, B2;//, gray2 = 128;
            byte R3, G3, B3;//, gray3 = 128;
            for (int x = rect.X; x < rect.Width + rect.X; x++)
            {
                R = bmp.GetPixel(x, rect.Y).R;
                G = bmp.GetPixel(x, rect.Y).G;
                B = bmp.GetPixel(x, rect.Y).B;
                bright = Color.FromArgb(R, G, B).GetBrightness();
                if (max11 < bright)
                    max11 = bright;
                if (min11 > bright)
                    min11 = bright;
                /*
                if (R != G || G != B)
                {
                    gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                    if (max11 < gray)
                        max11 = gray;
                    if (min11 > gray)
                        min11 = gray;
                }
                else
                {
                    if (max11 < B)
                        max11 = B;
                    if (min11 > B)
                        min11 = B;
                }*/
                R1 = bmp.GetPixel(x, rect.Height + rect.Y - 1).R;
                G1 = bmp.GetPixel(x, rect.Height + rect.Y - 1).G;
                B1 = bmp.GetPixel(x, rect.Height + rect.Y - 1).B;
                bright1 = Color.FromArgb(R1, G1, B1).GetBrightness();
                if (max12 < bright1)
                    max12 = bright1;
                if (min12 > bright1)
                    min12 = bright1;
                /*
                if (R1 != G1 || G1 != B1)
                {
                    gray1 = (byte)((19661 * R1 + 38666 * G1 + 7209 * B1) >> 16);
                    if (max12 < gray1)
                        max12 = gray1;
                    if (min12 > gray1)
                        min12 = gray1;
                }
                else
                {
                    if (max12 < B1)
                        max12 = B1;
                    if (min12 > B1)
                        min12 = B1;
                }*/
            }
            if (max11 > max12)
                max1 = (byte)(max11 - min12);
            else
                max1 = (byte)(max12 - min11);

            for (int y = rect.Y; y < rect.Height + rect.Y; y++)
            {
                R2 = bmp.GetPixel(rect.X, y).R;
                G2 = bmp.GetPixel(rect.X, y).G;
                B2 = bmp.GetPixel(rect.X, y).B;
                bright2 = Color.FromArgb(R2, G2, B2).GetBrightness();
                if (max21 < bright2)
                    max21 = bright2;
                if (min21 > bright2)
                    min21 = bright2;
                /*
                if (R2 != G2 || G2 != B2)
                {
                    gray2 = (byte)((19661 * R2 + 38666 * G2 + 7209 * B2) >> 16);
                    if (max21 < gray2)
                        max21 = gray2;
                    if (min21 > gray2)
                        min21 = gray2;
                }
                else
                {
                    if (max21 < B2)
                        max21 = B2;
                    if (min21 > B2)
                        min21 = B2;
                }*/
                R3 = bmp.GetPixel(rect.Width + rect.X - 1, y).R;
                G3 = bmp.GetPixel(rect.Width + rect.X - 1, y).G;
                B3 = bmp.GetPixel(rect.Width + rect.X - 1, y).B;
                bright3 = Color.FromArgb(R3, G3, B3).GetBrightness();
                if (max22 < bright3)
                    max22 = bright3;
                if (min22 > bright3)
                    min22 = bright3;
                /*
                if (R3 != G3 || G3 != B3)
                {
                    gray3 = (byte)((19661 * R3 + 38666 * G3 + 7209 * B3) >> 16);
                    if (max22 < gray3)
                        max22 = gray3;
                    if (min22 > gray3)
                        min22 = gray3;
                }
                else
                {
                    if (max22 < B3)
                        max22 = B3;
                    if (min22 > B3)
                        min22 = B3;
                }*/
            }
            if (max21 > max22)
                max2 = (byte)(max21 - min22);
            else
                max2 = (byte)(max22 - min21);

            //contrast = max1 > max2 ? max1 : max2;
            return max1 > max2 ? max1 : max2;
        }

        /// <summary>
        /// 改进的清晰度基算法[0.0,1.0]，较省时
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="rect">指定区域</param>
        /// <returns></returns>
        private static float GetOneContrastByBrightness(Bitmap bmp, Rectangle rect)
        {
            float max = 0;
            float min = 1F;
            byte R, G, B;
            float bright = 0;
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            BitmapData data = bmp.LockBits(rect, ImageLockMode.ReadOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - BPP * rect.Width;
                for (int y = 0; y < rect.Height; y++)
                {
                    for (int x = 0; x < rect.Width; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];
                        bright = Color.FromArgb(R, G, B).GetBrightness();
                        if (max < bright)
                            max = bright;
                        if (min > bright)
                            min = bright;
                        p += BPP;
                    }  // x
                    p += offset;
                } // y
            }
            bmp.UnlockBits(data);
            return (max - min);
        }

        /// <summary>
        /// 改进的清晰度基算法1[0.0,1.0]，较耗时
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="rect">指定区域</param>
        /// <returns></returns>
        private static float GetOneContrastByBrightness1(Bitmap bmp, Rectangle rect)
        {
            float sum = 0;
            byte R, G, B;
            float bright = 0;
            float b = GetBrightness(bmp, rect);
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            BitmapData data = bmp.LockBits(rect, ImageLockMode.ReadOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - BPP * rect.Width;
                for (int y = 0; y < rect.Height; y++)
                {
                    for (int x = 0; x < rect.Width; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];
                        bright = Color.FromArgb(R, G, B).GetBrightness();
                        sum += (bright - b) * (bright - b);
                        p += BPP;
                    }  // x
                    p += offset;
                } // y
            }
            bmp.UnlockBits(data);
            return (float)Math.Sqrt(sum / (rect.Width * rect.Height));
        }

        /// <summary>
        /// 改进的清晰度基算法[0.0,1.0]，较省时
        /// </summary>
        /// <returns></returns>
        private static float GetOneContrastByBrightness(Bitmap bmp)
        {
            float max = 0;
            float min = 1F;
            byte R, G, B;
            float bright = 0;

            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            BitmapData data = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - BPP * bmp.Width;
                for (int y = 0; y < bmp.Height; y++)
                {
                    for (int x = 0; x < bmp.Width; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];
                        bright = Color.FromArgb(R, G, B).GetBrightness();
                        if (max < bright)
                            max = bright;
                        if (min > bright)
                            min = bright;
                        p += BPP;
                    }  // x
                    p += offset;
                } // y
            }
            bmp.UnlockBits(data);
            return (max - min);
        }

        /// <summary>
        /// 改进的清晰度基算法1[0.0,1.0]，较耗时
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <returns></returns>
        private static float GetOneContrastByBrightness1(Bitmap bmp)
        {
            float sum = 0;
            byte R, G, B;
            float bright = 0;
            float b = GetBrightness(bmp);
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            BitmapData data = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - BPP * bmp.Width;
                for (int y = 0; y < bmp.Height; y++)
                {
                    for (int x = 0; x < bmp.Width; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];
                        bright = Color.FromArgb(R, G, B).GetBrightness();
                        sum += (bright - b) * (bright - b);
                        p += BPP;
                    }  // x
                    p += offset;
                } // y
            }
            bmp.UnlockBits(data);
            return (float)Math.Sqrt(sum / (bmp.Width * bmp.Height));
        }

        /// <summary>
        /// 改进的清晰度基算法(最新)
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="rect">指定区域</param>
        /// <returns></returns>
        private static byte GetOneContrastByGray(Bitmap bmp, Rectangle rect)
        {
            byte contrast = 0;
            byte R, G, B, gray = 128;
            int sum = 0;
            byte g = GetGray(bmp, rect);
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            BitmapData data = bmp.LockBits(rect, ImageLockMode.ReadOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - BPP * rect.Width;
                for (int y = 0; y < rect.Height; y++)
                {
                    for (int x = 0; x < rect.Width; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];

                        if (R != G || G != B)
                        {
                            gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                        }
                        else
                        {
                            gray = B;
                        }
                        sum += (gray - g) * (gray - g);
                        p += BPP;
                    }  // x
                    p += offset;
                } // y
            }
            bmp.UnlockBits(data);
            contrast = (byte)(Math.Sqrt(sum / (rect.Width * rect.Height)));
            return contrast;
        }

        /// <summary>
        /// 改进的清晰度基算法(最新)
        /// </summary>
        /// <returns></returns>
        private static byte GetOneContrastByGray(Bitmap bmp)
        {
            byte contrast = 0;
            byte R, G, B, gray = 128;
            int sum = 0;
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            byte g = GetGray(bmp);
            BitmapData data = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - BPP * bmp.Width;
                for (int y = 0; y < bmp.Height; y++)
                {
                    for (int x = 0; x < bmp.Width; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];

                        if (R != G || G != B)
                        {
                            gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                        }
                        else
                        {
                            gray = B;
                        }
                        sum += (gray - g) * (gray - g);
                        p += BPP;
                    }  // x
                    p += offset;
                } // y
            }
            bmp.UnlockBits(data);
            contrast = (byte)(Math.Sqrt(sum / (bmp.Width * bmp.Height)));
            return contrast;
        }

        /// <summary>
        /// 获取位图中指定区域的灰度
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="rect">指定区域</param>
        /// <returns></returns>
        public static byte GetGray(Bitmap bmp, Rectangle rect)
        {
            byte R, G, B, gray = 128;
            int sum = 0;
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            BitmapData data = bmp.LockBits(rect, ImageLockMode.ReadOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - BPP * rect.Width;
                for (int y = 0; y < rect.Height; y++)
                {
                    for (int x = 0; x < rect.Width; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];

                        if (R != G || G != B)
                        {
                            gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                        }
                        else
                        {
                            gray = B;
                        }
                        sum += gray;
                        p += BPP;
                    }  // x
                    p += offset;
                } // y
            }
            bmp.UnlockBits(data);
            return (byte)(sum / (rect.Width * rect.Height));
        }
        /// <summary>
        /// 获取位图的灰度
        /// </summary>
        /// <returns></returns>
        public byte GetGray()
        {
            byte R, G, B, gray = 128;
            int sum = 0;
            BitmapData data = bmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - BPP * width;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];

                        if (R != G || G != B)
                        {
                            gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                        }
                        else
                        {
                            gray = B;
                        }
                        sum += gray;
                        p += BPP;
                    }  // x
                    p += offset;
                } // y
            }
            bmp.UnlockBits(data);
            return (byte)(sum / (width * height));
        }
        /// <summary>
        /// 获取位图的灰度
        /// </summary>
        /// <param name="bmp"></param>
        /// <returns></returns>
        public static byte GetGray(Bitmap bmp)
        {
            byte R, G, B, gray = 128;
            int sum = 0;
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            BitmapData data = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - BPP * bmp.Width;
                for (int y = 0; y < bmp.Height; y++)
                {
                    for (int x = 0; x < bmp.Width; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];

                        if (R != G || G != B)
                        {
                            gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                        }
                        else
                        {
                            gray = B;
                        }
                        sum += gray;
                        p += BPP;
                    }  // x
                    p += offset;
                } // y
            }
            bmp.UnlockBits(data);
            return (byte)(sum / (bmp.Width * bmp.Height));
        }
        /// <summary>
        /// 获取位图中指定区域的灰度中值
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="rect"></param>
        /// <returns></returns>
        public static byte GetMidGray(Bitmap bmp, Rectangle rect)
        {
            Bitmap b = FastClipBitmap((Bitmap)bmp.Clone(), rect);
            Histogram hist = new Histogram(b);
            b.Dispose();
            return (byte)hist.Gray.Median;
        }
        /// <summary>
        /// 获取位图的灰度中值
        /// </summary>
        /// <param name="bmp"></param>
        /// <returns></returns>
        public static byte GetMidGray(Bitmap bmp)
        {
            Histogram hist = new Histogram((Bitmap)bmp.Clone());
            return (byte)hist.Gray.Median;
        }
        /// <summary>
        /// 获取位图的亮度
        /// </summary>
        /// <returns></returns>
        public float GetBrightness()
        {
            byte R, G, B;
            float sum = 0;
            BitmapData data = bmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - BPP * width;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];

                        sum += Color.FromArgb(R, G, B).GetBrightness();
                        p += BPP;
                    }  // x
                    p += offset;
                } // y
            }
            bmp.UnlockBits(data);
            return (byte)(sum / (width * height));
        }
        /// <summary>
        /// 获取位图的亮度
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <returns></returns>
        public static float GetBrightness(Bitmap bmp)
        {
            byte R, G, B;
            float sum = 0;
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            BitmapData data = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - BPP * bmp.Width;
                for (int y = 0; y < bmp.Height; y++)
                {
                    for (int x = 0; x < bmp.Width; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];

                        sum += Color.FromArgb(R, G, B).GetBrightness();
                        p += BPP;
                    }  // x
                    p += offset;
                } // y
            }
            bmp.UnlockBits(data);
            return sum / (bmp.Width * bmp.Height);
        }
        /// <summary>
        /// 获取位图中指定区域的亮度
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="rect"></param>
        /// <returns></returns>
        public static float GetBrightness(Bitmap bmp, Rectangle rect)
        {
            byte R, G, B;
            float sum = 0;
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            BitmapData data = bmp.LockBits(rect, ImageLockMode.ReadOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - BPP * rect.Width;
                for (int y = 0; y < rect.Height; y++)
                {
                    for (int x = 0; x < rect.Width; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];

                        sum += Color.FromArgb(R, G, B).GetBrightness();
                        p += BPP;
                    }  // x
                    p += offset;
                } // y
            }
            bmp.UnlockBits(data);
            return sum / (rect.Width * rect.Height);
        }
        /// <summary>
        /// 获取位图中指定区域的灰度中值
        /// </summary>
        /// <param name="rect"></param>
        /// <returns></returns>
        public byte GetMidBrightness(Rectangle rect)
        {
            Bitmap b = FastClipBitmap((Bitmap)bmp.Clone(), rect);
            Histogram hist = new Histogram(b);
            b.Dispose();
            return (byte)hist.Bright.Median;
        }
        /// <summary>
        /// 获取位图的灰度中值
        /// </summary>
        /// <returns></returns>
        public byte GetMidBrightness()
        {
            Histogram hist = new Histogram((Bitmap)bmp.Clone());
            return (byte)hist.Bright.Median;
        }
        /// <summary>
        /// 根据blkORwht, threshold是否找到边缘
        /// </summary>
        /// <param name="rect"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public bool FindEdge(Rectangle rect, int blkORwht, byte threshold)
        {
            bool find = false;
            if (bmp != null)
            {
                Rectangle rec = Rectangle.Intersect(rect, new Rectangle(0, 0, width, height));
                int X = rec.X;
                int Y = rec.Y;
                int Width = rec.Width + rec.X;
                int Height = rec.Height + rec.Y;
                byte R, G, B;
                byte gray = 128;
                unsafe
                {
                    BitmapData data = bmp.LockBits(new Rectangle(X, Y, rec.Width, rec.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
                    byte* p = (byte*)data.Scan0;
                    int stride = data.Stride;
                    int offset = stride - BPP * rec.Width;
                    for (int y = Y; y < Height; y++)
                    {
                        for (int x = X; x < Width; x++)
                        {
                            R = p[2];
                            G = p[1];
                            B = p[0];
                            if (R != G || G != B)
                            {
                                //尚未是灰度化，先计算灰度
                                // gray = 0.3*R + 0.59*G + 0.11*B
                                gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                                //bmp.SetPixel(x, y, Color.FromArgb(gray, gray, gray);
                                if (blkORwht == 0)//黑边
                                {
                                    if (gray < threshold)
                                    {
                                        find = true;
                                        break;
                                    }
                                }
                                else//白边
                                {
                                    if (gray > threshold)
                                    {
                                        find = true;
                                        break;
                                    }
                                }
                            }
                            else
                            {
                                if (blkORwht == 0)//黑边
                                {
                                    if (B < threshold)
                                    {
                                        find = true;
                                        break;
                                    }
                                }
                                else//白边
                                {
                                    if (B > threshold)
                                    {
                                        find = true;
                                        break;
                                    }
                                }
                            }
                            p += BPP;
                        }
                        if (find)
                            break;
                        p += offset;
                    }
                    bmp.UnlockBits(data);
                }
            }
            return find;
        }
        /// <summary>
        /// 根据color是否找到边缘
        /// </summary>
        /// <param name="rect"></param>
        /// <param name="color"></param>
        /// <returns></returns>
        public bool FindEdge(Rectangle rect, Color color)
        {
            bool find = false;
            if (bmp != null)
            {
                if (color == Color.Empty)
                    return false;
                Rectangle rec = Rectangle.Intersect(rect, new Rectangle(0, 0, width, height));
                int X = rec.X;
                int Y = rec.Y;
                int Width = rec.Width + rec.X;
                int Height = rec.Height + rec.Y;
                unsafe
                {
                    BitmapData data = bmp.LockBits(new Rectangle(X, Y, rec.Width, rec.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
                    byte* p = (byte*)data.Scan0;
                    int stride = data.Stride;
                    int offset = stride - BPP * rec.Width;
                    for (int y = Y; y < Height; y++)
                    {
                        for (int x = X; x < Width; x++)
                        {
                            if (p[0] == color.B && p[1] == color.G && p[2] == color.R && p[3] == color.A)
                            {
                                find = true;
                                break;
                            }
                            p += BPP;
                        }
                        if (find)
                            break;
                        p += offset;
                    }
                    bmp.UnlockBits(data);
                }
            }
            return find;
        }
        /// <summary>
        /// 根据点集生成二值位图(因为三重循环的缘故，速度太慢，留待以后改进)
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="blkORwht">点集为黑还是白，0为黑，非0则为白</param>
        /// <returns></returns>
        public static Bitmap GetBmpFromPs(Point[] ps, int blkORwht)
        {
            Bitmap bmp = null;
            if (ps != null)
            {
                bmp = new Bitmap(ps.GetLength(0), ps.GetLength(1), PixelFormat.Format24bppRgb);
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                BitmapData data = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
                unsafe
                {
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride - BPP * bmp.Width;
                    for (int y = 0; y < bmp.Height; y++)
                    {
                        for (int x = 0; x < bmp.Width; x++)
                        {
                            for (int i = 0; i < ps.Length; i++)
                            {
                                if (blkORwht == 0)
                                {
                                    if (ps[i].X == x && ps[i].Y == y)
                                    {
                                        p[0] = p[1] = p[2] = 0;
                                    }
                                    else
                                    {
                                        p[0] = p[1] = p[2] = 255;
                                    }
                                }
                                else
                                {
                                    if (ps[i].X == x && ps[i].Y == y)
                                    {
                                        p[0] = p[1] = p[2] = 255;
                                    }
                                    else
                                    {
                                        p[0] = p[1] = p[2] = 0;
                                    }
                                }
                            }
                            p += BPP;
                        }  // x
                        p += offset;
                    } // y
                }
                bmp.UnlockBits(data);
            }
            return bmp;
        }

        /// <summary>
        /// 根据点集生成彩色位图(因为三重循环的缘故，速度太慢，留待以后改进)
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="color">输出点集的颜色</param>
        /// <returns></returns>
        public static Bitmap GetBmpFromPs(Point[] ps, Color color)
        {
            Bitmap bmp = null;
            if (ps != null)
            {
                bmp = new Bitmap(ps.GetLength(0), ps.GetLength(1), PixelFormat.Format24bppRgb);
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                BitmapData data = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadWrite, bmp.PixelFormat);
                unsafe
                {
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride - BPP * bmp.Width;
                    for (int y = 0; y < bmp.Height; y++)
                    {
                        for (int x = 0; x < bmp.Width; x++)
                        {
                            for (int i = 0; i < ps.Length; i++)
                            {
                                if (ps[i].X == x && ps[i].Y == y)
                                {
                                    p[0] = color.B;
                                    p[1] = color.G;
                                    p[2] = color.R;
                                }
                            }
                            p += BPP;
                        }  // x
                        p += offset;
                    } // y
                }
                bmp.UnlockBits(data);
            }
            return bmp;
        }

        /// <summary>
        /// 移除点集中指定索引的一个点
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="index"></param>
        /// <returns></returns>
        public Point[] RemovePs(Point[] ps, int index)
        {
            Point[] ps1 = new Point[ps.Length - 1];
            int c = 0;
            for (int i = 0; i < ps.Length; i++)
            {
                if (i != index)
                {
                    ps1[c] = ps[i];
                    c++;
                }
            }
            return ps1;
        }
        /*
        /// <summary>
        /// 快速移除点集中指定索引的一个点
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public void RemovePs(int index)
        {
            tmpPs.RemoveAt(index);
        }
        */
        /// <summary>
        /// 灰度直方图的灰度值百分比数组(256个元素，元素值的取值范围[0, 100])
        /// </summary>
        /// <param name="bmp"></param>
        /// <returns></returns>
        public static uint[] GrayHistogram(Bitmap bmp)
        {
            if (bmp != null)
            {
                uint[] hist = new uint[256];
                uint maxScale = 1;
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                BitmapData data = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
                unsafe
                {
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride - BPP * bmp.Width;
                    for (int y = 0; y < bmp.Height; y++)
                    {
                        for (int x = 0; x < bmp.Width; x++)
                        {
                            byte B = p[0];
                            byte G = p[1];
                            byte R = p[2];
                            // 获取像素的亮度
                            //byte brightness = (byte)(0.3 * R + 0.59 * G + 0.11 * B);
                            byte brightness = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                            hist[brightness]++;
                            if (hist[brightness] > maxScale)
                            {
                                maxScale = hist[brightness];          //记住最大的数值
                            }
                            p += BPP;
                        }  // x
                        p += offset;
                    } // y
                }
                bmp.UnlockBits(data);
                int totalPixel = bmp.Width * bmp.Height;
                for (int i = 0; i < 256; i++)
                {
                    histogram[i] = (uint)(100.0 * hist[i] / maxScale);
                }
            }
            return histogram;
        }

        /// <summary>
        /// 获取给定图片的指定通道下的直方图，返回直方图的大小为[256, 256]
        /// </summary>
        /// <param name="img"></param>
        /// <param name="hst"></param>
        /// <returns></returns>
        public static Bitmap Histogram(Image img, HistogramType hst)
        {
            int[] scales = new int[256];                       //保存各阶亮度的统计
            byte R, G, B;
            int maxScale = 1;
            Bitmap bmp = (Bitmap)img;
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            BitmapData data = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - bmp.Width * 4;
                for (int y = 0; y < bmp.Height; y++)
                {
                    for (int x = 0; x < bmp.Width; x++)
                    {
                        R = p[2];
                        G = p[1];
                        B = p[0];
                        byte grayscale = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                        switch (hst)
                        {
                            case HistogramType.Red:
                                grayscale = R;
                                break;
                            case HistogramType.Green:
                                grayscale = G;
                                break;
                            case HistogramType.Blue:
                                grayscale = B;
                                break;
                            case HistogramType.Gray:
                            default:
                                break;
                        }

                        scales[grayscale]++;                       //该亮度的统计增一

                        if (scales[grayscale] > maxScale)
                        {
                            maxScale = scales[grayscale];          //记住最大的数值
                        }

                        p += 4;
                    }

                    p += offset;
                }
            }
            bmp.UnlockBits(data);

            for (int i = 0; i < scales.Length; i++)
            {
                scales[i] = scales[i] * 255 / maxScale;        //把亮度数组缩小到0~255区间，以便用图像直观表示出来
            }

            Bitmap result = new Bitmap(256, 256);              //准备一个直方图
            int HistogramHeight = 256;
            using (Graphics g = Graphics.FromImage(result))
            {
                for (int x = 0; x < result.Height; x++)
                {
                    switch (hst)
                    {
                        case HistogramType.Red:
                            g.DrawLine(Pens.Red, x, HistogramHeight, x, HistogramHeight - scales[x]);  //每个色阶画一条线，长度依据该色阶的统计数值
                            break;
                        case HistogramType.Green:
                            g.DrawLine(Pens.Green, x, HistogramHeight, x, HistogramHeight - scales[x]);  //每个色阶画一条线，长度依据该色阶的统计数值
                            break;
                        case HistogramType.Blue:
                            g.DrawLine(Pens.Blue, x, HistogramHeight, x, HistogramHeight - scales[x]);  //每个色阶画一条线，长度依据该色阶的统计数值
                            break;
                        case HistogramType.Gray:
                        default:
                            g.DrawLine(Pens.Gray, x, HistogramHeight, x, HistogramHeight - scales[x]);  //每个色阶画一条线，长度依据该色阶的统计数值
                            break;
                    }
                }
            }
            return result;
        }

        /// <summary>
        /// 透明化
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="transparent">透明度，范围[0, 255]</param>
        /// <returns></returns>
        public static Bitmap Opacitize(Bitmap bmp, byte transparent)
        {
            Bitmap b = new Bitmap(bmp.Width, bmp.Height, bmp.PixelFormat);
            Rectangle rect = new Rectangle(0, 0, bmp.Width, bmp.Height);
            int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
            BitmapData data = bmp.LockBits(rect, ImageLockMode.ReadOnly, bmp.PixelFormat);
            BitmapData d = b.LockBits(rect, ImageLockMode.WriteOnly, bmp.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                byte* p1 = (byte*)d.Scan0;
                int stride = data.Stride;
                int offset = stride - bmp.Width * 4;
                for (int y = 0; y < bmp.Height; y++)
                {
                    for (int x = 0; x < bmp.Width; x++)
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            p1[i] = p[i];
                        }
                        p1[3] = transparent;
                        p += 4;
                        p1 += 4;
                    }
                    p += offset;
                    p1 += offset;
                }
            }
            bmp.UnlockBits(data);
            b.UnlockBits(d);
            return b;
        }

        /// <summary>
        /// 初始化一个整型二维数组
        /// </summary>
        /// <param name="width">数组宽</param>
        /// <param name="height">数组高</param>
        /// <param name="init">初始值</param>
        /// <returns></returns>
        public static byte[,] InitArray(int width, int height, byte init)
        {
            byte[,] dst = new byte[width, height];

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    dst[x, y] = init;
                } // x
            } // y

            return dst;
        } // end of InitArray
        /// <summary>
        /// 初始化一个浮点二维数组
        /// </summary>
        /// <param name="width">数组宽</param>
        /// <param name="height">数组高</param>
        /// <param name="init">初始值</param>
        /// <returns></returns>
        public static float[,] InitArray(int width, int height, float init)
        {
            float[,] dst = new float[width, height];

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    dst[x, y] = init;
                } // x
            } // y

            return dst;
        } // end of InitArray

        /// <summary>
        /// 初始化一个颜色数组
        /// </summary>
        /// <param name="width">数组宽</param>
        /// <param name="height">数组高</param>
        /// <param name="init">初始值</param>
        /// <returns></returns>
        public static Color[,] InitArray(int width, int height, Color init)
        {
            Color[,] dst = new Color[width, height];

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    dst[x, y] = init;
                } // x
            } // y

            return dst;
        } // end of InitArray
        /// <summary>
        /// 连通
        /// </summary>
        /// <param name="inputPs"></param>
        /// <returns></returns>
        public static List<Point> Linking(List<Point> inputPs)
        {
            List<Point> tmpPs = new List<Point>();
            int minX = int.MaxValue, minY = int.MaxValue, maxX = 0, maxY = 0;
            byte[,] bMx;
            Point LocationPs;
            Size SizePs;
            foreach (Point p in inputPs)
            {
                if (minX > p.X)
                    minX = p.X;
                if (maxX < p.X)
                    maxX = p.X;
                if (minY > p.Y)
                    minY = p.Y;
                if (maxY < p.Y)
                    maxY = p.Y;
            }
            LocationPs = new Point(minX, minY);
            SizePs = new Size(maxX - minX, maxY - minY);
            bMx = InitArray((int)SizePs.Width + 1, (int)SizePs.Height + 1, 0);
            foreach (PointF p in inputPs)
            {
                bMx[(int)(p.X - LocationPs.X), (int)(p.Y - LocationPs.Y)] = 1;
            }
            for (int i = 0; i < bMx.GetLength(0); i++)
            {
                for (int j = 0; j < bMx.GetLength(1); j++)
                {
                    if (i - 2 < 0 && j - 2 >= 0 && i + 2 <= bMx.GetLength(0) - 1 && j + 2 <= bMx.GetLength(1) - 1)//左
                    {
                        if (bMx[i, j] == 1 && j > 1 && bMx[i, j - 2] == 1)//up
                        {
                            bMx[i, j - 1] = 1;
                        }
                        if (bMx[i, j] == 1 && i < bMx.GetLength(0) - 2 && bMx[i + 2, j] == 1)//right
                        {
                            bMx[i + 1, j] = 1;
                        }
                        if (bMx[i, j] == 1 && j < bMx.GetLength(1) - 2 && bMx[i, j + 2] == 1)//down
                        {
                            bMx[i, j + 1] = 1;
                        }
                    }
                    else if (j - 2 < 0 && i + 2 <= bMx.GetLength(0) - 1 && j + 2 <= bMx.GetLength(1) - 1)//上
                    {
                        if (bMx[i, j] == 1 && i > 1 && bMx[i - 2, j] == 1)//left
                        {
                            bMx[i - 1, j] = 1;
                        }
                        if (bMx[i, j] == 1 && i < bMx.GetLength(0) - 2 && bMx[i + 2, j] == 1)//right
                        {
                            bMx[i + 1, j] = 1;
                        }
                        if (bMx[i, j] == 1 && j < bMx.GetLength(1) - 2 && bMx[i, j + 2] == 1)//down
                        {
                            bMx[i, j + 1] = 1;
                        }
                    }
                    else if (i + 2 > bMx.GetLength(0) - 1 && j + 2 <= bMx.GetLength(1) - 1)//右
                    {
                        if (bMx[i, j] == 1 && i > 1 && bMx[i - 2, j] == 1)//left
                        {
                            bMx[i - 1, j] = 1;
                        }
                        if (bMx[i, j] == 1 && j > 1 && bMx[i, j - 2] == 1)//up
                        {
                            bMx[i, j - 1] = 1;
                        }
                        if (bMx[i, j] == 1 && j < bMx.GetLength(1) - 2 && bMx[i, j + 2] == 1)//down
                        {
                            bMx[i, j + 1] = 1;
                        }
                    }
                    else if (j + 2 > bMx.GetLength(1) - 1)//下
                    {
                        if (bMx[i, j] == 1 && i > 1 && bMx[i - 2, j] == 1)//left
                        {
                            bMx[i - 1, j] = 1;
                        }
                        if (bMx[i, j] == 1 && j > 1 && bMx[i, j - 2] == 1)//up
                        {
                            bMx[i, j - 1] = 1;
                        }
                        if (bMx[i, j] == 1 && i < bMx.GetLength(0) - 2 && bMx[i + 2, j] == 1)//right
                        {
                            bMx[i + 1, j] = 1;
                        }
                    }
                    else
                    {
                        if (bMx[i, j] == 1 && bMx[i, j - 2] == 1)//up
                        {
                            bMx[i, j - 1] = 1;
                        }
                        if (bMx[i, j] == 1 && bMx[i, j + 2] == 1)//down
                        {
                            bMx[i, j + 1] = 1;
                        }
                        if (bMx[i, j] == 1 && bMx[i - 2, j] == 1)//left
                        {
                            bMx[i - 1, j] = 1;
                        }
                        if (bMx[i, j] == 1 && bMx[i + 2, j] == 1)//right
                        {
                            bMx[i + 1, j] = 1;
                        }
                    }
                }
            }
            for (int i = 0; i < bMx.GetLength(0); i++)
            {
                for (int j = 0; j < bMx.GetLength(1); j++)
                {
                    if (bMx[i, j] == 1)
                        tmpPs.Add(new Point(LocationPs.X + i, LocationPs.Y + j));
                }
            }
            return tmpPs;
        }
        /// <summary>
        /// 连通
        /// </summary>
        /// <param name="inputPs"></param>
        /// <returns></returns>
        public static List<PointF> Linking(List<PointF> inputPs)
        {
            List<PointF> tmpPs1 = new List<PointF>();
            float minX = int.MaxValue, minY = int.MaxValue, maxX = 0, maxY = 0;
            float[,] bMx;
            PointF LocationPs;
            SizeF SizePs;
            foreach (PointF p in inputPs)
            {
                if (minX > p.X)
                    minX = p.X;
                if (maxX < p.X)
                    maxX = p.X;
                if (minY > p.Y)
                    minY = p.Y;
                if (maxY < p.Y)
                    maxY = p.Y;
            }
            LocationPs = new PointF(minX, minY);
            SizePs = new SizeF(maxX - minX, maxY - minY);
            bMx = InitArray((int)SizePs.Width + 1, (int)SizePs.Height + 1, 0.0F);
            foreach (PointF p in inputPs)
            {
                bMx[(int)(p.X - LocationPs.X), (int)(p.Y - LocationPs.Y)] = 1F;
            }
            for (int i = 0; i < bMx.GetLength(0); i++)
            {
                for (int j = 0; j < bMx.GetLength(1); j++)
                {
                    if (i - 2 < 0 && j - 2 >= 0 && i + 2 <= bMx.GetLength(0) - 1 && j + 2 <= bMx.GetLength(1) - 1)//左
                    {
                        if (bMx[i, j] == 1 && j > 1 && bMx[i, j - 2] == 1)//up
                        {
                            bMx[i, j - 1] = 1F;
                        }
                        if (bMx[i, j] == 1 && i < bMx.GetLength(0) - 2 && bMx[i + 2, j] == 1)//right
                        {
                            bMx[i + 1, j] = 1F;
                        }
                        if (bMx[i, j] == 1 && j < bMx.GetLength(1) - 2 && bMx[i, j + 2] == 1)//down
                        {
                            bMx[i, j + 1] = 1F;
                        }
                    }
                    else if (j - 2 < 0 && i + 2 <= bMx.GetLength(0) - 1 && j + 2 <= bMx.GetLength(1) - 1)//上
                    {
                        if (bMx[i, j] == 1 && i > 1 && bMx[i - 2, j] == 1)//left
                        {
                            bMx[i - 1, j] = 1F;
                        }
                        if (bMx[i, j] == 1 && i < bMx.GetLength(0) - 2 && bMx[i + 2, j] == 1)//right
                        {
                            bMx[i + 1, j] = 1F;
                        }
                        if (bMx[i, j] == 1 && j < bMx.GetLength(1) - 2 && bMx[i, j + 2] == 1)//down
                        {
                            bMx[i, j + 1] = 1F;
                        }
                    }
                    else if (i + 2 > bMx.GetLength(0) - 1 && j + 2 <= bMx.GetLength(1) - 1)//右
                    {
                        if (bMx[i, j] == 1 && i > 1 && bMx[i - 2, j] == 1)//left
                        {
                            bMx[i - 1, j] = 1F;
                        }
                        if (bMx[i, j] == 1 && j > 1 && bMx[i, j - 2] == 1)//up
                        {
                            bMx[i, j - 1] = 1F;
                        }
                        if (bMx[i, j] == 1 && j < bMx.GetLength(1) - 2 && bMx[i, j + 2] == 1)//down
                        {
                            bMx[i, j + 1] = 1F;
                        }
                    }
                    else if (j + 2 > bMx.GetLength(1) - 1)//下
                    {
                        if (bMx[i, j] == 1 && i > 1 && bMx[i - 2, j] == 1)//left
                        {
                            bMx[i - 1, j] = 1F;
                        }
                        if (bMx[i, j] == 1 && j > 1 && bMx[i, j - 2] == 1)//up
                        {
                            bMx[i, j - 1] = 1F;
                        }
                        if (bMx[i, j] == 1 && i < bMx.GetLength(0) - 2 && bMx[i + 2, j] == 1)//right
                        {
                            bMx[i + 1, j] = 1F;
                        }
                    }
                    else
                    {
                        if (bMx[i, j] == 1 && bMx[i, j - 2] == 1)//up
                        {
                            bMx[i, j - 1] = 1F;
                        }
                        if (bMx[i, j] == 1 && bMx[i, j + 2] == 1)//down
                        {
                            bMx[i, j + 1] = 1F;
                        }
                        if (bMx[i, j] == 1 && bMx[i - 2, j] == 1)//left
                        {
                            bMx[i - 1, j] = 1F;
                        }
                        if (bMx[i, j] == 1 && bMx[i + 2, j] == 1)//right
                        {
                            bMx[i + 1, j] = 1F;
                        }
                    }
                }
            }
            for (int i = 0; i < bMx.GetLength(0); i++)
            {
                for (int j = 0; j < bMx.GetLength(1); j++)
                {
                    if (bMx[i, j] == 1F)
                        tmpPs1.Add(new PointF(LocationPs.X + i, LocationPs.Y + j));
                }
            }
            return tmpPs1;
        }

        /// <summary>
        /// 取要细化点集的中点的脚标(1个或2个)
        /// </summary>
        /// <returns></returns>
        public static List<int> GetThinningMidPs(List<int> thinningPs)
        {
            List<int> midPs = new List<int>();
            int len = thinningPs.Count;
            if (len % 2 == 0)
            {
                midPs.Add(thinningPs[len / 2 - 1]);
                midPs.Add(thinningPs[len / 2]);
            }
            else
            {
                midPs.Add(thinningPs[len / 2]);
            }
            return midPs;
        }
        /// <summary>
        /// 取要细化点集的中点的X或Y值(1个)
        /// </summary>
        /// <returns></returns>
        public static int GetThinningMidPs1(List<int> thinningPs)
        {
            int midPs;
            int len = thinningPs.Count;
            if (len == 1)
                return thinningPs[0];
            int sum = 0;
            for (int i = 0; i < len; i++)
            {
                sum += thinningPs[i];
            }
            midPs = sum / len;
            return midPs;
        }
        /// <summary>
        /// 取要细化点集的中点的X或Y值(1个)
        /// </summary>
        /// <returns></returns>
        public static float GetThinningMidPs1(List<float> thinningPs)
        {
            float midPs;
            int len = thinningPs.Count;
            if (len == 1)
                return thinningPs[0];
            float sum = 0;
            for (int i = 0; i < len; i++)
            {
                sum += thinningPs[i];
            }
            midPs = sum / len;
            return midPs;
        }
        /// <summary>
        /// 细化 ---- 暂时无效
        /// </summary>
        /// <param name="inputPs"></param>
        /// <returns></returns>
        public List<Point> Thinning(Point[] inputPs)
        {
            List<Point> tmpPs = new List<Point>();
            int minX = int.MaxValue, minY = int.MaxValue, maxX = 0, maxY = 0;
            byte[,] bMx;
            Point LocationPs;
            Size SizePs;
            foreach (Point p in inputPs)
            {
                if (minX > p.X)
                    minX = p.X;
                if (maxX < p.X)
                    maxX = p.X;
                if (minY > p.Y)
                    minY = p.Y;
                if (maxY < p.Y)
                    maxY = p.Y;
            }
            LocationPs = new Point(minX, minY);
            SizePs = new Size(maxX - minX, maxY - minY);
            bMx = InitArray(SizePs.Width + 1, SizePs.Height + 1, 0);

            foreach (Point p in inputPs)
            {
                bMx[p.X - LocationPs.X, p.Y - LocationPs.Y] = 1;
            }

            byte[,] bMx1 = (byte[,])bMx.Clone();

            //水平方向取中点(一点或两点)
            for (int i = 0; i < bMx.GetLength(0); i++)
            {
                int have = 0;
                System.Collections.ArrayList al = new System.Collections.ArrayList();
                for (int j = 0; j < bMx.GetLength(1); j++)
                {
                    List<int> al1 = new List<int>();
                    if (bMx[i, j] == 1)
                    {
                        have++;
                        al1.Add(j);
                    }
                    else if (have > 0)
                    {
                        al.Add(al1);
                        have = 0;
                    }
                }
                int max = 0;
                int maxInd = 0;
                for (int i1 = 0; i1 < al.Count; i1++)
                {//选择长度最大的取其中点                    
                    int len = ((List<int>)al[i1]).Count;
                    if (max < len)
                    {
                        max = len;
                        maxInd = i1;
                    }
                }
                List<int> ps = new List<int>();
                for (int i2 = 0; i2 < max; i2++)
                {
                    ps.Add(((List<int>)al[maxInd])[i2]);
                }
                List<int> ps1 = null;
                if (ps.Count > 0)
                {
                    ps1 = GetThinningMidPs(ps);
                    for (int j1 = 0; j1 < bMx.GetLength(1); j1++)
                    {
                        for (int j2 = 0; j2 < ps1.Count; j2++)
                        {
                            if (j1 == ps1[j2])
                            {
                                bMx[i, j1] = 1;
                            }
                            else
                            {
                                bMx[i, j1] = 0;
                            }
                        }
                    }
                }
            }

            //垂直方向取中点(一点或两点)
            for (int j = 0; j < bMx1.GetLength(1); j++)
            {
                int have = 0;
                System.Collections.ArrayList al = new System.Collections.ArrayList();
                for (int i = 0; i < bMx1.GetLength(0); i++)
                {
                    List<int> al1 = new List<int>();
                    if (bMx1[i, j] == 1)
                    {
                        have++;
                        al1.Add(j);
                    }
                    else if (have > 0)
                    {
                        al.Add(al1);
                        have = 0;
                    }
                }
                int max = 0;
                int maxInd = 0;
                for (int i1 = 0; i1 < al.Count; i1++)
                {//选择长度最大的取其中点                    
                    int len = ((List<int>)al[i1]).Count;
                    if (max < len)
                    {
                        max = len;
                        maxInd = i1;
                    }
                }
                List<int> ps = new List<int>();
                for (int i2 = 0; i2 < max; i2++)
                {
                    ps.Add(((List<int>)al[maxInd])[i2]);
                }
                List<int> ps1 = null;
                if (ps.Count > 0)
                {
                    ps1 = GetThinningMidPs(ps);
                    for (int j1 = 0; j1 < bMx1.GetLength(0); j1++)
                    {
                        for (int j2 = 0; j2 < ps1.Count; j2++)
                        {
                            if (j1 == ps1[j2])
                            {
                                bMx1[j, j1] = 1;
                            }
                            else
                            {
                                bMx1[j, j1] = 0;
                            }
                        }
                    }
                }
            }

            //合并bMx和bMx1
            for (int j = 0; j < bMx.GetLength(1); j++)
            {
                for (int i = 0; i < bMx.GetLength(0); i++)
                {
                    if (bMx1[i, j] == 1)
                        bMx[i, j] = 1;
                }
            }

            for (int i = 0; i < bMx.GetLength(0); i++)
            {
                for (int j = 0; j < bMx.GetLength(1); j++)
                {
                    if (bMx[i, j] == 1)
                        tmpPs.Add(new Point(LocationPs.X + i, LocationPs.Y + j));
                }
            }
            return tmpPs;
        }
        /// <summary>
        /// 取点集的中间点（如弧线、直线的中点，与原点集的分布有关且与分布的疏密有关，有修正功能）
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static Point GetMidPointOfPs1(List<Point> ps)
        {
            //Point[] p = new Point[1];
            Point p;
            if (ps.Count < 1)
                throw new Exception("点数太少！ps.Count < 1");
            p = ps[ps.Count / 2];//初始化为中间值
            List<int> psXInd = new List<int>();
            List<int> psYInd = new List<int>();

            for (int i = 0; i < ps.Count; i++)
            {
                psXInd.Add(ps[i].X);
                psYInd.Add(ps[i].Y);
            }
            //先连通
            //List<Point> tmpPs = Linking(ps);
            //Point[] ps1 = tmpPs.ToArray();
            //后细化
            int psX = GetThinningMidPs1(psXInd);
            int psY = GetThinningMidPs1(psYInd);
            /*//返回点一定在原点集上
            foreach (Point pe in ps1)
            {
                if (pe.X == psX && pe.Y == psY)
                {
                    Point[] pp = new Point[1];
                    pp[0] = pe;
                    return pp;
                }
            }*/
            //返回点不一定在原点集上，即含有修正功能
            p = new Point(psX, psY);
            return p;
        }
        /// <summary>
        /// 取点集的中间点（如弧线、直线的中点，与原点集的分布有关且与分布的疏密有关，有修正功能）
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static PointF GetMidPointOfPs1(List<PointF> ps)
        {
            //Point[] p = new Point[1];
            PointF p;
            if (ps.Count < 1)
                throw new Exception("点数太少！ps.Count < 1");
            p = ps[ps.Count / 2];//初始化为中间值
            List<float> psXInd = new List<float>();
            List<float> psYInd = new List<float>();

            for (int i = 0; i < ps.Count; i++)
            {
                psXInd.Add(ps[i].X);
                psYInd.Add(ps[i].Y);
            }
            //先连通
            //List<PointF> tmpPs = Linking(ps);
            //PointF[] ps1 = tmpPs.ToArray();
            //后细化
            float psX = GetThinningMidPs1(psXInd);
            float psY = GetThinningMidPs1(psYInd);
            /*//返回点一定在原点集上
            foreach (Point pe in ps1)
            {
                if (pe.X == psX && pe.Y == psY)
                {
                    Point[] pp = new Point[1];
                    pp[0] = pe;
                    return pp;
                }
            }*/
            //返回点不一定在原点集上，即含有修正功能
            p = new PointF(psX, psY);
            return p;
        }
        /// <summary>
        /// 取点集的中间点（返回点一定在原点集上，与原点集的分布无关）
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static Point GetMidPointOfPs2(List<Point> ps)
        {
            //Point[] p = new Point[1];
            Point p;
            if (ps.Count < 1)
                throw new Exception("点数太少！ps.Count< 1");
            p = ps[ps.Count / 2];
            return p;
        }
        /// <summary>
        /// 取点集的中间点（返回点一定在原点集上，与原点集的分布无关）
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static PointF GetMidPointOfPs2(List<PointF> ps)
        {
            //Point[] p = new Point[1];
            PointF p;
            if (ps.Count < 1)
                throw new Exception("点数太少！ps.Count< 1");
            p = ps[ps.Count / 2];
            return p;
        }
        /// <summary>
        /// 取点集的中间点（与原点集的分布有关，但与分布的疏密无关）
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static Point GetMidPointOfPs3(List<Point> ps)
        {
            Rectangle re = RectangleOnMorePoint(ps);
            return new Point(re.X + re.Width / 2, re.Y + re.Height / 2);
        }
        /// <summary>
        /// 取点集的中间点（与原点集的分布有关，但与分布的疏密无关）
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static PointF GetMidPointOfPs3(List<PointF> ps)
        {
            RectangleF re = RectangleOnMorePoint(ps);
            return new PointF(re.X + re.Width / 2, re.Y + re.Height / 2);
        }
        /// <summary>
        /// 点集的外包正规矩形
        /// </summary>
        /// <param name="ps">点集</param>
        /// <returns></returns>
        public static Rectangle RectangleOnMorePoint(List<Point> ps)
        {
            Rectangle re = Rectangle.Empty;
            if (ps.Count > 0)
            {
                int maxx = ps[0].X, minx = ps[0].X, maxy = ps[0].Y, miny = ps[0].Y;
                for (int i = 1; i < ps.Count; i++)
                {
                    if (maxx < ps[i].X)
                        maxx = ps[i].X;
                    if (maxy < ps[i].Y)
                        maxy = ps[i].Y;
                    if (minx > ps[i].X)
                        minx = ps[i].X;
                    if (miny > ps[i].Y)
                        miny = ps[i].Y;
                }
                re = new Rectangle(minx, miny, maxx - minx, maxy - miny);
            }
            return re;
        }
        /// <summary>
        /// 点集的外包正规矩形
        /// </summary>
        /// <param name="ps">点集</param>
        /// <returns></returns>
        public static RectangleF RectangleOnMorePoint(List<PointF> ps)
        {
            RectangleF re = RectangleF.Empty;
            if (ps.Count > 0)
            {
                float maxx = ps[0].X, minx = ps[0].X, maxy = ps[0].Y, miny = ps[0].Y;
                for (int i = 1; i < ps.Count; i++)
                {
                    if (maxx < ps[i].X)
                        maxx = ps[i].X;
                    if (maxy < ps[i].Y)
                        maxy = ps[i].Y;
                    if (minx > ps[i].X)
                        minx = ps[i].X;
                    if (miny > ps[i].Y)
                        miny = ps[i].Y;
                }
                re = new RectangleF(minx, miny, maxx - minx, maxy - miny);
            }
            return re;
        }
        /// <summary>
        /// 取点集的指定范围点集（返回点集一定在原点集上，与原点集的分布无关）
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="index1"></param>
        /// <param name="index2"></param>
        /// <returns>返回点集，错误的范围索引则返回null</returns>
        public List<Point> GetRangePoints(List<Point> ps, int index1, int index2)
        {
            index1 = index1 < 0 ? 0 : index1;
            index1 = index1 > ps.Count - 1 ? ps.Count - 1 : index1;
            index2 = index2 < 0 ? 0 : index2;
            index2 = index2 > ps.Count - 1 ? ps.Count - 1 : index2;
            if (index1 > index2)
                return null;
            //Point[] range = new Point[index2 - index1 + 1];
            List<Point> range = new List<Point>();
            //int ii = 0;
            for (int i = index1; i < index2 + 1; i++)
            {
                //range[ii] = ps[i];
                range.Add(ps[i]);
                //ii++;
            }
            return range;
        }
        /// <summary>
        /// 取点集的指定范围点集（返回点集一定在原点集上，与原点集的分布无关）
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="index1"></param>
        /// <param name="index2"></param>
        /// <returns>返回点集，错误的范围索引则返回null</returns>
        public List<PointF> GetRangePoints(List<PointF> ps, int index1, int index2)
        {
            index1 = index1 < 0 ? 0 : index1;
            index1 = index1 > ps.Count - 1 ? ps.Count - 1 : index1;
            index2 = index2 < 0 ? 0 : index2;
            index2 = index2 > ps.Count - 1 ? ps.Count - 1 : index2;
            if (index1 > index2)
                return null;
            //Point[] range = new Point[index2 - index1 + 1];
            List<PointF> range = new List<PointF>();
            //int ii = 0;
            for (int i = index1; i < index2 + 1; i++)
            {
                //range[ii] = ps[i];
                range.Add(ps[i]);
                //ii++;
            }
            return range;
        }
        /// <summary>
        /// 取关键点
        /// </summary>
        /// <param name="bmp">要检测的位图</param>
        /// <param name="rect">检测区域</param>
        /// <param name="blkORwht">0为黑边，非0为白边</param>
        /// <param name="threshold">灰度阀值</param>
        /// <param name="needGrayExt"></param>
        /// <param name="needFilter"></param>
        /// <param name="needEdgeDetect"></param>
        /// <param name="edgeThinning"></param>
        /// <returns></returns>
        public static List<Point> GetKeyPs(Bitmap bmp, Rectangle rect, int blkORwht, byte threshold, bool needGrayExt, bool needFilter, bool needEdgeDetect, bool edgeThinning)
        {
            List<Point> tmpPs = new List<Point>();
            if (bmp != null)
            {
                if (rect.X > bmp.Width - 1 || rect.Y > bmp.Height - 1 || rect == Rectangle.Empty)
                    return null;
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                Rectangle rec = Rectangle.Intersect(rect, new Rectangle(0, 0, bmp.Width, bmp.Height));
                int X = rec.X;
                int Y = rec.Y;
                int Width = rec.Width + rec.X;
                int Height = rec.Height + rec.Y;
                Bitmap bm = null;
                if (needFilter)
                    bm = Filtering(bmp.Clone(rec, bmp.PixelFormat), blkORwht, threshold, 1);
                else
                    bm = Bitize(bmp.Clone(rec, bmp.PixelFormat), threshold, needGrayExt);
                if (needGrayExt)
                    bm = GrayExtend(bm, false);
                if (needEdgeDetect)
                {
                    bm = GetBitBmpByRoberts(bm, threshold, edgeThinning);
                }
                BitmapData data = bm.LockBits(new Rectangle(0, 0, rec.Width, rec.Height), ImageLockMode.ReadWrite, bm.PixelFormat);
                unsafe
                {
                    byte R, G, B;
                    byte gray = 128;
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride - BPP * (Width - X);
                    for (int y = Y; y < Height; y++)
                    {
                        for (int x = X; x < Width; x++)
                        {
                            R = p[2];
                            G = p[1];
                            B = p[0];
                            if (R != G || G != B)
                            {
                                //尚未是灰度化，先计算灰度
                                // gray = 0.3*R + 0.59*G + 0.11*B
                                gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                                //再提取关键点
                                if (blkORwht == 0 && !needEdgeDetect)//黑边及不取边缘
                                {
                                    if (gray < threshold)
                                        tmpPs.Add(new Point(x, y));
                                }
                                else//白边
                                {
                                    if (gray > threshold)
                                        tmpPs.Add(new Point(x, y));
                                }
                            }
                            else
                            {
                                if (blkORwht == 0 && !needEdgeDetect)//黑边及不取边缘
                                {
                                    if (B < threshold)
                                        tmpPs.Add(new Point(x, y));
                                }
                                else//白边
                                {
                                    if (B > threshold)
                                        tmpPs.Add(new Point(x, y));
                                }
                            }
                            p += BPP;
                        } // x
                        p += offset;
                    } // y
                }
                bm.UnlockBits(data);
            }
            return tmpPs;
        }
        /// <summary>
        /// 取关键点
        /// </summary>
        /// <param name="bmp">要检测的位图</param>
        /// <param name="path">检测路径</param>
        /// <param name="blkORwht">0为黑边，非0为白边</param>
        /// <param name="threshold">灰度阀值</param>
        /// <param name="needGrayExt"></param>
        /// <param name="needFilter"></param>
        /// <param name="needEdgeDetect"></param>
        /// <param name="edgeThinning"></param>
        /// <returns></returns>
        public static List<Point> GetKeyPs(Bitmap bmp, GraphicsPath path, int blkORwht, byte threshold, bool needGrayExt, bool needFilter, bool needEdgeDetect, bool edgeThinning)
        {
            List<Point> tmpPs = new List<Point>();
            if (bmp != null)
            {
                RectangleF rectF = path.GetBounds();
                Rectangle rect = Rectangle.Intersect(Rectangle.Truncate(rectF), new Rectangle(0, 0, bmp.Width, bmp.Height));
                if (rect.X > bmp.Width - 1 || rect.Y > bmp.Height - 1 || rect == Rectangle.Empty)
                    return null;
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                int X = rect.X;
                int Y = rect.Y;
                int Width = rect.Width;
                int Height = rect.Height;
                Bitmap bm = null;
                if (needFilter)
                    bm = Filtering(bmp.Clone(rect, bmp.PixelFormat), blkORwht, threshold, 1);
                else
                    bm = Bitize(bmp.Clone(rect, bmp.PixelFormat), threshold, needGrayExt);
                if (needGrayExt)
                    bm = GrayExtend(bm, true);
                if (needEdgeDetect)
                {
                    bm = GetBitBmpByRoberts(bm, threshold, edgeThinning);
                }
                BitmapData data = bm.LockBits(new Rectangle(0, 0, Width, Height), ImageLockMode.ReadWrite, bm.PixelFormat);
                unsafe
                {
                    byte R, G, B;
                    byte gray = 128;
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride - BPP * Width;
                    for (int y = Y; y < Height + Y; y++)
                    {
                        for (int x = X; x < Width + X; x++)
                        {
                            Point pp = new Point(x, y);
                            if (path.IsVisible(pp))
                            {
                                R = p[2];
                                G = p[1];
                                B = p[0];
                                if (R != G || G != B)
                                {
                                    //尚未是灰度化，先计算灰度
                                    gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                                    //再提取关键点
                                    if (blkORwht == 0 && !needEdgeDetect)//黑边及不取边缘
                                    {
                                        if (gray < threshold)
                                            tmpPs.Add(pp);
                                    }
                                    else//白边
                                    {
                                        if (gray > threshold)
                                            tmpPs.Add(pp);
                                    }
                                }
                                else
                                {
                                    if (blkORwht == 0 && !needEdgeDetect)//黑边及不取边缘
                                    {
                                        if (B < threshold)
                                            tmpPs.Add(pp);
                                    }
                                    else//白边
                                    {
                                        if (B > threshold)
                                            tmpPs.Add(pp);
                                    }
                                }
                            }
                            p += BPP;
                        } // x
                        p += offset;
                    } // y
                }
                bm.UnlockBits(data);
            }
            return tmpPs;
        }
        /// <summary>
        /// 取关键点
        /// </summary>
        /// <param name="bmp">要检测的位图</param>
        /// <param name="region">检测区域</param>
        /// <param name="blkORwht">0为黑边，非0为白边</param>
        /// <param name="threshold">灰度阀值</param>
        /// <param name="needGrayExt"></param>
        /// <param name="needFilter"></param>
        /// <param name="needEdgeDetect"></param>
        /// <param name="edgeThinning"></param>
        /// <returns></returns>
        public static List<Point> GetKeyPs(Bitmap bmp, Region region, int blkORwht, byte threshold, bool needGrayExt, bool needFilter, bool needEdgeDetect, bool edgeThinning)
        {
            List<Point> tmpPs = new List<Point>();
            if (bmp != null)
            {
                Graphics g = Graphics.FromImage(bmp);
                RectangleF rectF = region.GetBounds(g);
                Rectangle rect = Rectangle.Intersect(Rectangle.Truncate(rectF), new Rectangle(0, 0, bmp.Width, bmp.Height));
                if (rect.X > bmp.Width - 1 || rect.Y > bmp.Height - 1 || rect == Rectangle.Empty)
                    return null;
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                int X = rect.X;
                int Y = rect.Y;
                int Width = rect.Width;
                int Height = rect.Height;
                Bitmap bm = null;
                if (needFilter)
                    bm = Filtering(bmp.Clone(rect, bmp.PixelFormat), blkORwht, threshold, 1);
                else
                    bm = Bitize(bmp.Clone(rect, bmp.PixelFormat), threshold, needGrayExt);
                if (needGrayExt)
                    bm = GrayExtend(bm, true);
                if (needEdgeDetect)
                {
                    bm = GetBitBmpByRoberts(bm, threshold, edgeThinning);
                }
                BitmapData data = bm.LockBits(new Rectangle(0, 0, Width, Height), ImageLockMode.ReadWrite, bm.PixelFormat);
                unsafe
                {
                    byte R, G, B;
                    byte gray = 128;
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride - BPP * Width;
                    for (int y = Y; y < Height + Y; y++)
                    {
                        for (int x = X; x < Width + X; x++)
                        {
                            Point pp = new Point(x, y);
                            if (region.IsVisible(pp))
                            {
                                R = p[2];
                                G = p[1];
                                B = p[0];
                                if (R != G || G != B)
                                {
                                    //尚未是灰度化，先计算灰度
                                    gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                                    //再提取关键点
                                    if (blkORwht == 0 && !needEdgeDetect)//黑边及不取边缘
                                    {
                                        if (gray < threshold)
                                            tmpPs.Add(pp);
                                    }
                                    else//白边
                                    {
                                        if (gray > threshold)
                                            tmpPs.Add(pp);
                                    }
                                }
                                else
                                {
                                    if (blkORwht == 0 && !needEdgeDetect)//黑边及不取边缘
                                    {
                                        if (B < threshold)
                                            tmpPs.Add(pp);
                                    }
                                    else//白边
                                    {
                                        if (B > threshold)
                                            tmpPs.Add(pp);
                                    }
                                }
                            }
                            p += BPP;
                        } // x
                        p += offset;
                    } // y
                }
                bm.UnlockBits(data);
            }
            return tmpPs;
        }
        /// <summary>
        /// 快速剪裁指定路径的位图
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="path"></param>
        /// <returns></returns>
        public static Bitmap FastClipBitmap(Bitmap bmp, GraphicsPath path)
        {
            RectangleF rect = path.GetBounds();
            //Bitmap srcBmp = (Bitmap)bmp.Clone();
            PixelFormat pf = bmp.PixelFormat;
            int BPP = Image.GetPixelFormatSize(pf) / 8;
            GraphicsUnit gu = GraphicsUnit.Pixel;
            Rectangle rect1 = Rectangle.Truncate(RectangleF.Intersect(bmp.GetBounds(ref gu), rect));
            if (rect1 == Rectangle.Empty || rect1.Width == 0 || rect1.Height == 0)
                return null;
            Bitmap dstImage = new Bitmap(rect1.Width, rect1.Height, pf);
            BitmapData srcData = bmp.LockBits(new Rectangle(rect1.X, rect1.Y, rect1.Width, rect1.Height), ImageLockMode.ReadOnly, pf);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, rect1.Width, rect1.Height), ImageLockMode.WriteOnly, pf);
            int stride1 = srcData.Stride;
            int offset1 = stride1 - BPP * rect1.Width;
            int stride2 = dstData.Stride;
            int offset2 = stride2 - BPP * rect1.Width;
            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;
                for (int y = rect1.Y; y < rect1.Height + rect1.Y; y++)
                {
                    for (int x = rect1.X; x < rect1.Width + rect1.X; x++)
                    {
                        /*if (PathIsRect(path))//加上这个更慢！
                        {
                            for (int i = 0; i < 4; i++)
                            {
                                dst[i] = src[i];
                            } // i
                        }
                        else  
                        {*/
                        if (path.IsVisible(x, y))
                        {
                            for (int i = 0; i < BPP; i++)
                            {
                                dst[i] = src[i];
                            } // i
                        }
                        else//路径外的像素点全置为透明
                        {
                            for (int i = 0; i < 3; i++)
                            {
                                dst[i] = 255;
                            } // i
                            if (BPP > 3)
                                dst[3] = 0;
                        }
                        //}
                        src += BPP;
                        dst += BPP;
                    } // x
                    src += offset1;
                    dst += offset2;
                } // y
            }
            bmp.UnlockBits(srcData);
            dstImage.UnlockBits(dstData);
            return dstImage;
        }
        /// <summary>
        /// 获取指定黑边(或白边)以threshold阀值呈现出来的边缘点集(指定路径)
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="path"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public static List<Point> GetEdgeByPath(Bitmap bmp, GraphicsPath path, byte threshold)
        {
            List<Point> tmpPs = new List<Point>();

            if (bmp != null)
            {
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                RectangleF rectF = path.GetBounds();
                Rectangle rect = Rectangle.Intersect(Rectangle.Truncate(rectF), new Rectangle(0, 0, bmp.Width, bmp.Height));
                int X = rect.X;
                int Y = rect.Y;
                int Width = rect.Width + X;
                int Height = rect.Height + Y;
                /*
                Bitmap bm = GetBmpByRoberts(bmp);
                //#####阀值化###
                Bitmap edgeImage = Thresholding(bm, rect, threshold, needGrayExt);              
                bm.Dispose();
                */
                Bitmap edgeImage = GetBitBmpByRoberts(bmp, rect, threshold, false);
                unsafe
                {
                    BitmapData data = edgeImage.LockBits(new Rectangle(1, 1, rect.Width - 1, rect.Height - 1), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                    byte* p = (byte*)data.Scan0;
                    int stride = data.Stride;
                    int offset = stride - BPP * (rect.Width - 1);
                    for (int j = Y + 1; j < Height; j++)
                    {
                        for (int i = X + 1; i < Width; i++)
                        {
                            Point pp = new Point(i, j);
                            /*if (PathIsRect(path))//加上这个更慢！
                            {
                                if (blkORwht == 0)
                                {
                                    if (p[0] == 0)
                                        MergePs(pp);
                                }
                                else
                                {
                                    if (p[0] == 255)
                                        MergePs(pp);
                                }
                            }
                            else if (path.IsVisible(pp))
                            {*/
                            if (path.IsVisible(pp))
                            {
                                byte gray = (byte)((19661 * p[2] + 38666 * p[1] + 7209 * p[0]) >> 16);
                                if (gray > 127)
                                    tmpPs.Add(pp);
                            }
                            p += BPP;
                        }
                        p += offset;
                    }
                    edgeImage.UnlockBits(data);
                }
            }
            return tmpPs;
        }
        /// <summary>
        /// 获取指定黑边(或白边)以threshold阀值呈现出来的边缘点集(指定区域)
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="rect"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public static List<Point> GetEdgeByRect(Bitmap bmp, Rectangle rect, byte threshold)
        {
            List<Point> tmpPs = new List<Point>();
            if (bmp != null)
            {
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                Rectangle rec = Rectangle.Intersect(rect, new Rectangle(0, 0, bmp.Width, bmp.Height));
                int X = rec.X;
                int Y = rec.Y;
                int Width = rec.Width + rec.X;
                int Height = rec.Height + rec.Y;
                /*
                Bitmap bm = GetBmpByRoberts(bmp);
                //#####阀值化###
                Bitmap edgeImage = Thresholding(bm, rec, threshold, needGrayExt);
                
                bm.Dispose();
                */
                Bitmap edgeImage = GetBitBmpByRoberts(bmp, rec, threshold, false);
                unsafe
                {
                    BitmapData data = edgeImage.LockBits(new Rectangle(1, 1, rec.Width - 1, rec.Height - 1), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                    byte* p = (byte*)data.Scan0;
                    int stride = data.Stride;
                    int offset = stride - BPP * (rec.Width - 1);
                    for (int j = Y + 1; j < Height; j++)
                    {
                        for (int i = X + 1; i < Width; i++)
                        {
                            Point pp = new Point(i, j);
                            byte gray = (byte)((19661 * p[2] + 38666 * p[1] + 7209 * p[0]) >> 16);
                            if (gray > 127)
                                tmpPs.Add(pp);
                            p += BPP;
                        }
                        p += offset;
                    }
                    edgeImage.UnlockBits(data);
                }
            }
            return tmpPs;
        }
        /// <summary>
        /// 获取指定黑边(或白边)以threshold阀值呈现出来的边缘点集(全位图)
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public static List<Point> GetEdge(Bitmap bmp, byte threshold)
        {
            List<Point> tmpPs = new List<Point>();

            if (bmp != null)
            {
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                Rectangle rec = new Rectangle(0, 0, bmp.Width, bmp.Height);
                int X = rec.X;
                int Y = rec.Y;
                int Width = rec.Width + rec.X;
                int Height = rec.Height + rec.Y;
                /*
                Bitmap bm = GetBmpByRoberts(bmp);
                //#####阀值化###
                Bitmap edgeImage = Thresholding(bm, rec, threshold, needGrayExt);                
                bm.Dispose();
                */
                Bitmap edgeImage = GetBitBmpByRoberts(bmp, threshold, false);
                unsafe
                {
                    BitmapData data = edgeImage.LockBits(new Rectangle(1, 1, rec.Width - 1, rec.Height - 1), ImageLockMode.ReadOnly, edgeImage.PixelFormat);
                    byte* p = (byte*)data.Scan0;
                    int stride = data.Stride;
                    int offset = stride - BPP * (rec.Width - 1);
                    for (int j = Y + 1; j < Height; j++)
                    {
                        for (int i = X + 1; i < Width; i++)
                        {
                            Point pp = new Point(i, j);
                            byte gray = (byte)((19661 * p[2] + 38666 * p[1] + 7209 * p[0]) >> 16);
                            if (gray > 127)
                                tmpPs.Add(pp);
                            p += BPP;
                        }
                        p += offset;
                    }
                    edgeImage.UnlockBits(data);
                }
            }
            return tmpPs;
        }
        /// <summary>
        /// 获取指定黑边(或白边)以threshold阀值呈现出来的边缘位图(指定路径)
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="rect"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public static Bitmap GetEdgeBmpByRect(Bitmap bmp, Rectangle rect, byte threshold)
        {
            Bitmap resultBmp = null;
            if (bmp != null)
            {
                GraphicsUnit gu = GraphicsUnit.Pixel;
                Rectangle rec = Rectangle.Intersect(rect, Rectangle.Truncate(bmp.GetBounds(ref gu)));
                /*
                //#####阀值化###
                Bitmap b = GetBmpByRoberts(bmp); 

                resultBmp = Thresholding(b, rec, threshold, needGrayExt);

                b.Dispose();
                */
                Bitmap edgeImage = GetBitBmpByRoberts(bmp, rec, threshold, false);
            }
            return resultBmp;
        }

        /// <summary>
        /// 检测整张位图在指定黑边(或白边)以threshold为阀值的情况下，是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="needFilter"></param>
        /// <param name="needEdgeDetect"></param>
        /// <param name="edgeThinning"></param>
        /// <returns></returns>
        public static bool IsExistKeyPoint4(Bitmap bmp, int blkORwht, byte threshold, bool needGrayExt, bool needFilter, bool needEdgeDetect, bool edgeThinning)
        {
            lock (lockObj)
            {
                MultiThreading4.exist1 = false;
                MultiThreading4.exist2 = false;
                MultiThreading4.exist3 = false;
                MultiThreading4.exist4 = false;
                if (bmp != null)
                {
                    int X = 0;
                    int Y = 0;
                    int Width = bmp.Width;
                    int Height = bmp.Height;
                    int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;

                    Bitmap bm = null;
                    if (needFilter)
                        bm = Filtering((Bitmap)bmp.Clone(), blkORwht, threshold, 1);
                    else
                        bm = Bitize((Bitmap)bmp.Clone(), threshold, needGrayExt);
                    if (needGrayExt)
                        bm = GrayExtend(bm, true);
                    if (needEdgeDetect)
                    {
                        bm = GetBitBmpByRoberts(bm, threshold, edgeThinning);
                    }

                    Rectangle rec = new Rectangle(0, 0, bm.Width, bm.Height);
                    Rectangle rec1 = new Rectangle(0, 0, rec.Width / 2, rec.Height / 2);
                    Rectangle rec2 = new Rectangle(rec.Width / 2, 0, rec.Width - rec.Width / 2, rec.Height / 2);
                    Rectangle rec3 = new Rectangle(0, rec.Height / 2, rec.Width / 2, rec.Height - rec.Height / 2);
                    Rectangle rec4 = new Rectangle(rec.Width / 2, rec.Height / 2, rec.Width - rec.Width / 2, rec.Height - rec.Height / 2);
                    RectBitmap rectbmp1 = new RectBitmap();
                    rectbmp1.blkORwht = blkORwht;
                    rectbmp1.threshold = threshold;
                    rectbmp1.bmp = FastClipBitmap(bm, rec1);
                    rec1.X += X;
                    rec1.Y += Y;
                    rectbmp1.rect = rec1;
                    RectBitmap rectbmp2 = new RectBitmap();
                    rectbmp2.blkORwht = blkORwht;
                    rectbmp2.threshold = threshold;
                    rectbmp2.bmp = FastClipBitmap(bm, rec2);
                    rec2.X += X;
                    rec2.Y += Y;
                    rectbmp2.rect = rec2;
                    RectBitmap rectbmp3 = new RectBitmap();
                    rectbmp3.blkORwht = blkORwht;
                    rectbmp3.threshold = threshold;
                    rectbmp3.bmp = FastClipBitmap(bm, rec3);
                    rec3.X += X;
                    rec3.Y += Y;
                    rectbmp3.rect = rec3;
                    RectBitmap rectbmp4 = new RectBitmap();
                    rectbmp4.blkORwht = blkORwht;
                    rectbmp4.threshold = threshold;
                    rectbmp4.bmp = FastClipBitmap(bm, rec4);
                    rec4.X += X;
                    rec4.Y += Y;
                    rectbmp4.rect = rec4;
                    AsyncMethodCaller caller = new AsyncMethodCaller(MultiThreading4.CheckWholeRect);
                    IAsyncResult result1 = caller.BeginInvoke(rectbmp1, out MultiThreading4.exist1, new AsyncCallback(MultiThreading4.CallbackTask1), caller);
                    IAsyncResult result2 = caller.BeginInvoke(rectbmp2, out MultiThreading4.exist2, new AsyncCallback(MultiThreading4.CallbackTask2), caller);
                    IAsyncResult result3 = caller.BeginInvoke(rectbmp3, out MultiThreading4.exist3, new AsyncCallback(MultiThreading4.CallbackTask3), caller);
                    IAsyncResult result4 = caller.BeginInvoke(rectbmp4, out MultiThreading4.exist4, new AsyncCallback(MultiThreading4.CallbackTask4), caller);
                    while (!result1.IsCompleted || !result2.IsCompleted || !result3.IsCompleted || !result4.IsCompleted)
                    {
                        //Application.DoEvents();
                    }
                    if (MultiThreading4.exist1 || MultiThreading4.exist2 || MultiThreading4.exist3 || MultiThreading4.exist4)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// 检测指定黑边(或白边)以threshold阀值在指定的区域中，是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="rect"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="needFilter"></param>
        /// <param name="needEdgeDetect"></param>
        /// <param name="edgeThinning"></param>
        /// <returns></returns>
        public static bool IsExistKeyPointByRect4(Bitmap bmp, Rectangle rect, int blkORwht, byte threshold, bool needGrayExt, bool needFilter, bool needEdgeDetect, bool edgeThinning)
        {
            lock (lockObj)
            {
                MultiThreading4.existRect1 = false;
                MultiThreading4.existRect2 = false;
                MultiThreading4.existRect3 = false;
                MultiThreading4.existRect4 = false;
                if (bmp != null)
                {
                    Rectangle rec = Rectangle.Intersect(rect, new Rectangle(0, 0, bmp.Width, bmp.Height));
                    int X = rec.X;
                    int Y = rec.Y;
                    int Width = rec.Width + rec.X;
                    int Height = rec.Height + rec.Y;
                    int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;

                    Bitmap bm = null;
                    if (needFilter)
                        bm = Filtering(bmp.Clone(rec, bmp.PixelFormat), blkORwht, threshold, 1);
                    else
                        bm = Bitize(bmp.Clone(rec, bmp.PixelFormat), threshold, needGrayExt);
                    if (needGrayExt)
                        bm = GrayExtend(bm, true);
                    if (needEdgeDetect)
                    {
                        bm = GetBitBmpByRoberts(bm, threshold, edgeThinning);
                    }

                    Rectangle rec1 = new Rectangle(0, 0, rec.Width / 2, rec.Height / 2);
                    Rectangle rec2 = new Rectangle(rec.Width / 2, 0, rec.Width - rec.Width / 2, rec.Height / 2);
                    Rectangle rec3 = new Rectangle(0, rec.Height / 2, rec.Width / 2, rec.Height - rec.Height / 2);
                    Rectangle rec4 = new Rectangle(rec.Width / 2, rec.Height / 2, rec.Width - rec.Width / 2, rec.Height - rec.Height / 2);
                    RectBitmap rectbmp1 = new RectBitmap();
                    rectbmp1.blkORwht = blkORwht;
                    rectbmp1.threshold = threshold;
                    rectbmp1.bmp = FastClipBitmap(bm, rec1);
                    rec1.X += X;
                    rec1.Y += Y;
                    rectbmp1.rect = rec1;
                    RectBitmap rectbmp2 = new RectBitmap();
                    rectbmp2.blkORwht = blkORwht;
                    rectbmp2.threshold = threshold;
                    rectbmp2.bmp = FastClipBitmap(bm, rec2);
                    rec2.X += X;
                    rec2.Y += Y;
                    rectbmp2.rect = rec2;
                    RectBitmap rectbmp3 = new RectBitmap();
                    rectbmp3.blkORwht = blkORwht;
                    rectbmp3.threshold = threshold;
                    rectbmp3.bmp = FastClipBitmap(bm, rec3);
                    rec3.X += X;
                    rec3.Y += Y;
                    rectbmp3.rect = rec3;
                    RectBitmap rectbmp4 = new RectBitmap();
                    rectbmp4.blkORwht = blkORwht;
                    rectbmp4.threshold = threshold;
                    rectbmp4.bmp = FastClipBitmap(bm, rec4);
                    rec4.X += X;
                    rec4.Y += Y;
                    rectbmp4.rect = rec4;
                    AsyncMethodCaller caller = new AsyncMethodCaller(MultiThreading4.CheckSomeRect);
                    IAsyncResult result1 = caller.BeginInvoke(rectbmp1, out MultiThreading4.existRect1, new AsyncCallback(MultiThreading4.CallbackTaskWithRect1), caller);
                    IAsyncResult result2 = caller.BeginInvoke(rectbmp2, out MultiThreading4.existRect2, new AsyncCallback(MultiThreading4.CallbackTaskWithRect2), caller);
                    IAsyncResult result3 = caller.BeginInvoke(rectbmp3, out MultiThreading4.existRect3, new AsyncCallback(MultiThreading4.CallbackTaskWithRect3), caller);
                    IAsyncResult result4 = caller.BeginInvoke(rectbmp4, out MultiThreading4.existRect4, new AsyncCallback(MultiThreading4.CallbackTaskWithRect4), caller);
                    while (!result1.IsCompleted || !result2.IsCompleted || !result3.IsCompleted || !result4.IsCompleted)
                    {
                        //Application.DoEvents();
                    }
                    if (MultiThreading4.existRect1 || MultiThreading4.existRect2 || MultiThreading4.existRect3 || MultiThreading4.existRect4)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// 检测指定黑边(或白边)以threshold阀值在指定的路径下，是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="path"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="needFilter"></param>
        /// <param name="needEdgeDetect"></param>
        /// <param name="edgeThinning"></param>
        /// <returns></returns>
        public static bool IsExistKeyPointByPath4(Bitmap bmp, GraphicsPath path, int blkORwht, byte threshold, bool needGrayExt, bool needFilter, bool needEdgeDetect, bool edgeThinning)
        {
            lock (lockObj)
            {
                MultiThreading4.existPath1 = false;
                MultiThreading4.existPath2 = false;
                MultiThreading4.existPath3 = false;
                MultiThreading4.existPath4 = false;
                if (bmp != null)
                {
                    RectangleF rect = path.GetBounds();
                    Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(rect), new Rectangle(0, 0, bmp.Width, bmp.Height));
                    int X = rec.X;
                    int Y = rec.Y;

                    Bitmap bm = null;
                    if (needFilter)
                        bm = Filtering(bmp.Clone(rec, bmp.PixelFormat), blkORwht, threshold, 1);
                    else
                        bm = Bitize(bmp.Clone(rec, bmp.PixelFormat), threshold, needGrayExt);
                    if (needGrayExt)
                        bm = GrayExtend(bm, true);
                    if (needEdgeDetect)
                    {
                        bm = GetBitBmpByRoberts(bm, threshold, edgeThinning);
                    }

                    Rectangle rec1 = new Rectangle(0, 0, rec.Width / 2, rec.Height / 2);
                    Rectangle rec2 = new Rectangle(rec.Width / 2, 0, rec.Width - rec.Width / 2, rec.Height / 2);
                    Rectangle rec3 = new Rectangle(0, rec.Height / 2, rec.Width / 2, rec.Height - rec.Height / 2);
                    Rectangle rec4 = new Rectangle(rec.Width / 2, rec.Height / 2, rec.Width - rec.Width / 2, rec.Height - rec.Height / 2);
                    RectBitmapWithPath rectbmp1 = new RectBitmapWithPath();
                    rectbmp1.blkORwht = blkORwht;
                    rectbmp1.threshold = threshold;
                    rectbmp1.path = (GraphicsPath)path.Clone();
                    rectbmp1.bmp = FastClipBitmap(bm, rec1);
                    rec1.X += X;
                    rec1.Y += Y;
                    rectbmp1.rect = rec1;
                    RectBitmapWithPath rectbmp2 = new RectBitmapWithPath();
                    rectbmp2.blkORwht = blkORwht;
                    rectbmp2.threshold = threshold;
                    rectbmp2.path = (GraphicsPath)path.Clone();
                    rectbmp2.bmp = FastClipBitmap(bm, rec2);
                    rec2.X += X;
                    rec2.Y += Y;
                    rectbmp2.rect = rec2;
                    RectBitmapWithPath rectbmp3 = new RectBitmapWithPath();
                    rectbmp3.blkORwht = blkORwht;
                    rectbmp3.threshold = threshold;
                    rectbmp3.path = (GraphicsPath)path.Clone();
                    rectbmp3.bmp = FastClipBitmap(bm, rec3);
                    rec3.X += X;
                    rec3.Y += Y;
                    rectbmp3.rect = rec3;
                    RectBitmapWithPath rectbmp4 = new RectBitmapWithPath();
                    rectbmp4.blkORwht = blkORwht;
                    rectbmp4.threshold = threshold;
                    rectbmp4.path = (GraphicsPath)path.Clone();
                    rectbmp4.bmp = FastClipBitmap(bm, rec4);
                    rec4.X += X;
                    rec4.Y += Y;
                    rectbmp4.rect = rec4;
                    AsyncMethodCallerWithPath caller = new AsyncMethodCallerWithPath(MultiThreading4.CheckSomeRectWithPath);
                    IAsyncResult result1 = caller.BeginInvoke(rectbmp1, out MultiThreading4.existPath1, new AsyncCallback(MultiThreading4.CallbackTaskWithPath1), caller);
                    IAsyncResult result2 = caller.BeginInvoke(rectbmp2, out MultiThreading4.existPath2, new AsyncCallback(MultiThreading4.CallbackTaskWithPath2), caller);
                    IAsyncResult result3 = caller.BeginInvoke(rectbmp3, out MultiThreading4.existPath3, new AsyncCallback(MultiThreading4.CallbackTaskWithPath3), caller);
                    IAsyncResult result4 = caller.BeginInvoke(rectbmp4, out MultiThreading4.existPath4, new AsyncCallback(MultiThreading4.CallbackTaskWithPath4), caller);
                    while (!result1.IsCompleted || !result2.IsCompleted || !result3.IsCompleted || !result4.IsCompleted)
                    {
                        //Application.DoEvents();
                    }
                    if (MultiThreading4.existPath1 || MultiThreading4.existPath2 || MultiThreading4.existPath3 || MultiThreading4.existPath4)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// 检测指定黑边(或白边)以threshold阀值在指定的区域中，是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="region"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="needFilter"></param>
        /// <param name="needEdgeDetect"></param>
        /// <param name="edgeThinning"></param>
        /// <returns></returns>
        public static bool IsExistKeyPointByRegion4(Bitmap bmp, Region region, int blkORwht, byte threshold, bool needGrayExt, bool needFilter, bool needEdgeDetect, bool edgeThinning)
        {
            lock (lockObj)
            {
                MultiThreading4.existRegion1 = false;
                MultiThreading4.existRegion2 = false;
                MultiThreading4.existRegion3 = false;
                MultiThreading4.existRegion4 = false;
                if (bmp != null)
                {
                    Region reg = region.Clone();
                    reg.Intersect(new Rectangle(0, 0, bmp.Width, bmp.Height));
                    Graphics g = Graphics.FromImage(bmp);
                    Rectangle rec = Rectangle.Truncate(reg.GetBounds(g));
                    int X = rec.X;
                    int Y = rec.Y;

                    Bitmap bm = null;
                    if (needFilter)
                        bm = Filtering(bmp.Clone(rec, bmp.PixelFormat), blkORwht, threshold, 1);
                    else
                        bm = Bitize(bmp.Clone(rec, bmp.PixelFormat), threshold, needGrayExt);
                    if (needGrayExt)
                        bm = GrayExtend(bm, true);
                    if (needEdgeDetect)
                    {
                        bm = GetBitBmpByRoberts(bm, threshold, edgeThinning);
                    }

                    Rectangle rec1 = new Rectangle(0, 0, rec.Width / 2, rec.Height / 2);
                    Rectangle rec2 = new Rectangle(rec.Width / 2, 0, rec.Width - rec.Width / 2, rec.Height / 2);
                    Rectangle rec3 = new Rectangle(0, rec.Height / 2, rec.Width / 2, rec.Height - rec.Height / 2);
                    Rectangle rec4 = new Rectangle(rec.Width / 2, rec.Height / 2, rec.Width - rec.Width / 2, rec.Height - rec.Height / 2);
                    RectBitmapWithRegion rectbmp1 = new RectBitmapWithRegion();
                    rectbmp1.blkORwht = blkORwht;
                    rectbmp1.threshold = threshold;
                    rectbmp1.region = reg.Clone();
                    rectbmp1.bmp = FastClipBitmap(bm, rec1);
                    rec1.X += X;
                    rec1.Y += Y;
                    rectbmp1.rect = rec1;
                    RectBitmapWithRegion rectbmp2 = new RectBitmapWithRegion();
                    rectbmp2.blkORwht = blkORwht;
                    rectbmp2.threshold = threshold;
                    rectbmp2.region = reg.Clone();
                    rectbmp2.bmp = FastClipBitmap(bm, rec2);
                    rec2.X += X;
                    rec2.Y += Y;
                    rectbmp2.rect = rec2;
                    RectBitmapWithRegion rectbmp3 = new RectBitmapWithRegion();
                    rectbmp3.blkORwht = blkORwht;
                    rectbmp3.threshold = threshold;
                    rectbmp3.region = reg.Clone();
                    rectbmp3.bmp = FastClipBitmap(bm, rec3);
                    rec3.X += X;
                    rec3.Y += Y;
                    rectbmp3.rect = rec3;
                    RectBitmapWithRegion rectbmp4 = new RectBitmapWithRegion();
                    rectbmp4.blkORwht = blkORwht;
                    rectbmp4.threshold = threshold;
                    rectbmp4.region = reg.Clone();
                    rectbmp4.bmp = FastClipBitmap(bm, rec4);
                    rec4.X += X;
                    rec4.Y += Y;
                    rectbmp4.rect = rec4;
                    AsyncMethodCallerWithRegion caller = new AsyncMethodCallerWithRegion(MultiThreading4.CheckSomeRectWithRegion);
                    IAsyncResult result1 = caller.BeginInvoke(rectbmp1, out MultiThreading4.existRegion1, new AsyncCallback(MultiThreading4.CallbackTaskWithRegion1), caller);
                    IAsyncResult result2 = caller.BeginInvoke(rectbmp2, out MultiThreading4.existRegion2, new AsyncCallback(MultiThreading4.CallbackTaskWithRegion2), caller);
                    IAsyncResult result3 = caller.BeginInvoke(rectbmp3, out MultiThreading4.existRegion3, new AsyncCallback(MultiThreading4.CallbackTaskWithRegion3), caller);
                    IAsyncResult result4 = caller.BeginInvoke(rectbmp4, out MultiThreading4.existRegion4, new AsyncCallback(MultiThreading4.CallbackTaskWithRegion4), caller);
                    reg.Dispose();
                    while (!result1.IsCompleted || !result2.IsCompleted || !result3.IsCompleted || !result4.IsCompleted)
                    {
                        //Application.DoEvents();
                    }
                    if (MultiThreading4.existRegion1 || MultiThreading4.existRegion2 || MultiThreading4.existRegion3 || MultiThreading4.existRegion4)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// 检测整张位图在指定黑边(或白边)以threshold为阀值的情况下，是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="needFilter"></param>
        /// <param name="needEdgeDetect"></param>
        /// <param name="edgeThinning"></param>
        /// <returns></returns>
        public static bool IsExistKeyPoint9(Bitmap bmp, int blkORwht, byte threshold, bool needGrayExt, bool needFilter, bool needEdgeDetect, bool edgeThinning)
        {
            lock (lockObj)
            {
                MultiThreading9.exist1 = false;
                MultiThreading9.exist2 = false;
                MultiThreading9.exist3 = false;
                MultiThreading9.exist4 = false;
                MultiThreading9.exist5 = false;
                MultiThreading9.exist6 = false;
                MultiThreading9.exist7 = false;
                MultiThreading9.exist8 = false;
                MultiThreading9.exist9 = false;
                if (bmp != null)
                {
                    int X = 0;
                    int Y = 0;
                    int Width = bmp.Width;
                    int Height = bmp.Height;
                    int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;

                    Bitmap bm = null;
                    if (needFilter)
                        bm = Filtering((Bitmap)bmp.Clone(), blkORwht, threshold, 1);
                    else
                        bm = Bitize((Bitmap)bmp.Clone(), threshold, needGrayExt);
                    if (needGrayExt)
                        bm = GrayExtend(bm, true);
                    if (needEdgeDetect)
                    {
                        bm = GetBitBmpByRoberts(bm, threshold, edgeThinning);
                    }

                    Rectangle rec = new Rectangle(0, 0, bm.Width, bm.Height);
                    Rectangle rec1 = new Rectangle(0, 0, rec.Width / 3, rec.Height / 3);
                    Rectangle rec2 = new Rectangle(rec.Width / 3, 0, rec.Width / 3, rec.Height / 3);
                    Rectangle rec3 = new Rectangle(rec.Width * 2 / 3, 0, rec.Width - rec.Width * 2 / 3, rec.Height / 3);
                    Rectangle rec4 = new Rectangle(0, rec.Height / 3, rec.Width / 3, rec.Height / 3);
                    Rectangle rec5 = new Rectangle(rec.Width / 3, rec.Height / 3, rec.Width / 3, rec.Height / 3);
                    Rectangle rec6 = new Rectangle(rec.Width * 2 / 3, rec.Height / 3, rec.Width - rec.Width * 2 / 3, rec.Height / 3);
                    Rectangle rec7 = new Rectangle(0, rec.Height * 2 / 3, rec.Width / 3, rec.Height - rec.Height * 2 / 3);
                    Rectangle rec8 = new Rectangle(rec.Width / 3, rec.Height * 2 / 3, rec.Width / 3, rec.Height - rec.Height * 2 / 3);
                    Rectangle rec9 = new Rectangle(rec.Width * 2 / 3, rec.Height * 2 / 3, rec.Width - rec.Width * 2 / 3, rec.Height - rec.Height * 2 / 3);
                    RectBitmap rectbmp1 = new RectBitmap();
                    rectbmp1.blkORwht = blkORwht;
                    rectbmp1.threshold = threshold;
                    rectbmp1.bmp = FastClipBitmap(bm, rec1);
                    rec1.X += X;
                    rec1.Y += Y;
                    rectbmp1.rect = rec1;
                    RectBitmap rectbmp2 = new RectBitmap();
                    rectbmp2.blkORwht = blkORwht;
                    rectbmp2.threshold = threshold;
                    rectbmp2.bmp = FastClipBitmap(bm, rec2);
                    rec2.X += X;
                    rec2.Y += Y;
                    rectbmp2.rect = rec2;
                    RectBitmap rectbmp3 = new RectBitmap();
                    rectbmp3.blkORwht = blkORwht;
                    rectbmp3.threshold = threshold;
                    rectbmp3.bmp = FastClipBitmap(bm, rec3);
                    rec3.X += X;
                    rec3.Y += Y;
                    rectbmp3.rect = rec3;
                    RectBitmap rectbmp4 = new RectBitmap();
                    rectbmp4.blkORwht = blkORwht;
                    rectbmp4.threshold = threshold;
                    rectbmp4.bmp = FastClipBitmap(bm, rec4);
                    rec4.X += X;
                    rec4.Y += Y;
                    rectbmp4.rect = rec4;
                    RectBitmap rectbmp5 = new RectBitmap();
                    rectbmp5.blkORwht = blkORwht;
                    rectbmp5.threshold = threshold;
                    rectbmp5.bmp = FastClipBitmap(bm, rec5);
                    rec5.X += X;
                    rec5.Y += Y;
                    rectbmp5.rect = rec5;
                    RectBitmap rectbmp6 = new RectBitmap();
                    rectbmp6.blkORwht = blkORwht;
                    rectbmp6.threshold = threshold;
                    rectbmp6.bmp = FastClipBitmap(bm, rec6);
                    rec6.X += X;
                    rec6.Y += Y;
                    rectbmp6.rect = rec6;
                    RectBitmap rectbmp7 = new RectBitmap();
                    rectbmp7.blkORwht = blkORwht;
                    rectbmp7.threshold = threshold;
                    rectbmp7.bmp = FastClipBitmap(bm, rec7);
                    rec7.X += X;
                    rec7.Y += Y;
                    rectbmp7.rect = rec7;
                    RectBitmap rectbmp8 = new RectBitmap();
                    rectbmp8.blkORwht = blkORwht;
                    rectbmp8.threshold = threshold;
                    rectbmp8.bmp = FastClipBitmap(bm, rec8);
                    rec8.X += X;
                    rec8.Y += Y;
                    rectbmp8.rect = rec8;
                    RectBitmap rectbmp9 = new RectBitmap();
                    rectbmp9.blkORwht = blkORwht;
                    rectbmp9.threshold = threshold;
                    rectbmp9.bmp = FastClipBitmap(bm, rec9);
                    rec9.X += X;
                    rec9.Y += Y;
                    rectbmp9.rect = rec9;
                    AsyncMethodCaller caller = new AsyncMethodCaller(MultiThreading9.CheckWholeRect);
                    IAsyncResult result1 = caller.BeginInvoke(rectbmp1, out MultiThreading9.exist1, new AsyncCallback(MultiThreading9.CallbackTask1), caller);
                    IAsyncResult result2 = caller.BeginInvoke(rectbmp2, out MultiThreading9.exist2, new AsyncCallback(MultiThreading9.CallbackTask2), caller);
                    IAsyncResult result3 = caller.BeginInvoke(rectbmp3, out MultiThreading9.exist3, new AsyncCallback(MultiThreading9.CallbackTask3), caller);
                    IAsyncResult result4 = caller.BeginInvoke(rectbmp4, out MultiThreading9.exist4, new AsyncCallback(MultiThreading9.CallbackTask4), caller);
                    IAsyncResult result5 = caller.BeginInvoke(rectbmp5, out MultiThreading9.exist5, new AsyncCallback(MultiThreading9.CallbackTask5), caller);
                    IAsyncResult result6 = caller.BeginInvoke(rectbmp6, out MultiThreading9.exist6, new AsyncCallback(MultiThreading9.CallbackTask6), caller);
                    IAsyncResult result7 = caller.BeginInvoke(rectbmp7, out MultiThreading9.exist7, new AsyncCallback(MultiThreading9.CallbackTask7), caller);
                    IAsyncResult result8 = caller.BeginInvoke(rectbmp8, out MultiThreading9.exist8, new AsyncCallback(MultiThreading9.CallbackTask8), caller);
                    IAsyncResult result9 = caller.BeginInvoke(rectbmp9, out MultiThreading9.exist9, new AsyncCallback(MultiThreading9.CallbackTask9), caller);
                    while (!result1.IsCompleted || !result2.IsCompleted || !result3.IsCompleted || !result4.IsCompleted || !result5.IsCompleted || !result6.IsCompleted || !result7.IsCompleted || !result8.IsCompleted || !result9.IsCompleted)
                    {
                        //Application.DoEvents();
                    }
                    if (MultiThreading9.exist1 || MultiThreading9.exist2 || MultiThreading9.exist3 || MultiThreading9.exist4 || MultiThreading9.exist5 || MultiThreading9.exist6 || MultiThreading9.exist7 || MultiThreading9.exist8 || MultiThreading9.exist9)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// 检测指定黑边(或白边)以threshold阀值在指定的区域中，是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="rect"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="needFilter"></param>
        /// <param name="needEdgeDetect"></param>
        /// <param name="edgeThinning"></param>
        /// <returns></returns>
        public static bool IsExistKeyPointByRect9(Bitmap bmp, Rectangle rect, int blkORwht, byte threshold, bool needGrayExt, bool needFilter, bool needEdgeDetect, bool edgeThinning)
        {
            lock (lockObj)
            {
                MultiThreading9.existRect1 = false;
                MultiThreading9.existRect2 = false;
                MultiThreading9.existRect3 = false;
                MultiThreading9.existRect4 = false;
                MultiThreading9.existRect5 = false;
                MultiThreading9.existRect6 = false;
                MultiThreading9.existRect7 = false;
                MultiThreading9.existRect8 = false;
                MultiThreading9.existRect9 = false;
                if (bmp != null)
                {
                    Rectangle rec = Rectangle.Intersect(rect, new Rectangle(0, 0, bmp.Width, bmp.Height));
                    int X = rec.X;
                    int Y = rec.Y;
                    int Width = rec.Width + rec.X;
                    int Height = rec.Height + rec.Y;
                    int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;

                    Bitmap bm = null;
                    if (needFilter)
                        bm = Filtering(bmp.Clone(rec, bmp.PixelFormat), blkORwht, threshold, 1);
                    else
                        bm = Bitize(bmp.Clone(rec, bmp.PixelFormat), threshold, needGrayExt);
                    if (needGrayExt)
                        bm = GrayExtend(bm, true);
                    if (needEdgeDetect)
                    {
                        bm = GetBitBmpByRoberts(bm, threshold, edgeThinning);
                    }

                    Rectangle rec1 = new Rectangle(0, 0, rec.Width / 3, rec.Height / 3);
                    Rectangle rec2 = new Rectangle(rec.Width / 3, 0, rec.Width / 3, rec.Height / 3);
                    Rectangle rec3 = new Rectangle(rec.Width * 2 / 3, 0, rec.Width - rec.Width * 2 / 3, rec.Height / 3);
                    Rectangle rec4 = new Rectangle(0, rec.Height / 3, rec.Width / 3, rec.Height / 3);
                    Rectangle rec5 = new Rectangle(rec.Width / 3, rec.Height / 3, rec.Width / 3, rec.Height / 3);
                    Rectangle rec6 = new Rectangle(rec.Width * 2 / 3, rec.Height / 3, rec.Width - rec.Width * 2 / 3, rec.Height / 3);
                    Rectangle rec7 = new Rectangle(0, rec.Height * 2 / 3, rec.Width / 3, rec.Height - rec.Height * 2 / 3);
                    Rectangle rec8 = new Rectangle(rec.Width / 3, rec.Height * 2 / 3, rec.Width / 3, rec.Height - rec.Height * 2 / 3);
                    Rectangle rec9 = new Rectangle(rec.Width * 2 / 3, rec.Height * 2 / 3, rec.Width - rec.Width * 2 / 3, rec.Height - rec.Height * 2 / 3);
                    RectBitmap rectbmp1 = new RectBitmap();
                    rectbmp1.blkORwht = blkORwht;
                    rectbmp1.threshold = threshold;
                    rectbmp1.bmp = FastClipBitmap(bm, rec1);
                    rec1.X += X;
                    rec1.Y += Y;
                    rectbmp1.rect = rec1;
                    RectBitmap rectbmp2 = new RectBitmap();
                    rectbmp2.blkORwht = blkORwht;
                    rectbmp2.threshold = threshold;
                    rectbmp2.bmp = FastClipBitmap(bm, rec2);
                    rec2.X += X;
                    rec2.Y += Y;
                    rectbmp2.rect = rec2;
                    RectBitmap rectbmp3 = new RectBitmap();
                    rectbmp3.blkORwht = blkORwht;
                    rectbmp3.threshold = threshold;
                    rectbmp3.bmp = FastClipBitmap(bm, rec3);
                    rec3.X += X;
                    rec3.Y += Y;
                    rectbmp3.rect = rec3;
                    RectBitmap rectbmp4 = new RectBitmap();
                    rectbmp4.blkORwht = blkORwht;
                    rectbmp4.threshold = threshold;
                    rectbmp4.bmp = FastClipBitmap(bm, rec4);
                    rec4.X += X;
                    rec4.Y += Y;
                    rectbmp4.rect = rec4;
                    RectBitmap rectbmp5 = new RectBitmap();
                    rectbmp5.blkORwht = blkORwht;
                    rectbmp5.threshold = threshold;
                    rectbmp5.bmp = FastClipBitmap(bm, rec5);
                    rec5.X += X;
                    rec5.Y += Y;
                    rectbmp5.rect = rec5;
                    RectBitmap rectbmp6 = new RectBitmap();
                    rectbmp6.blkORwht = blkORwht;
                    rectbmp6.threshold = threshold;
                    rectbmp6.bmp = FastClipBitmap(bm, rec6);
                    rec6.X += X;
                    rec6.Y += Y;
                    rectbmp6.rect = rec6;
                    RectBitmap rectbmp7 = new RectBitmap();
                    rectbmp7.blkORwht = blkORwht;
                    rectbmp7.threshold = threshold;
                    rectbmp7.bmp = FastClipBitmap(bm, rec7);
                    rec7.X += X;
                    rec7.Y += Y;
                    rectbmp7.rect = rec7;
                    RectBitmap rectbmp8 = new RectBitmap();
                    rectbmp8.blkORwht = blkORwht;
                    rectbmp8.threshold = threshold;
                    rectbmp8.bmp = FastClipBitmap(bm, rec8);
                    rec8.X += X;
                    rec8.Y += Y;
                    rectbmp8.rect = rec8;
                    RectBitmap rectbmp9 = new RectBitmap();
                    rectbmp9.blkORwht = blkORwht;
                    rectbmp9.threshold = threshold;
                    rectbmp9.bmp = FastClipBitmap(bm, rec9);
                    rec9.X += X;
                    rec9.Y += Y;
                    rectbmp9.rect = rec9;
                    AsyncMethodCaller caller = new AsyncMethodCaller(MultiThreading9.CheckSomeRect);
                    IAsyncResult result1 = caller.BeginInvoke(rectbmp1, out MultiThreading9.existRect1, new AsyncCallback(MultiThreading9.CallbackTaskWithRect1), caller);
                    IAsyncResult result2 = caller.BeginInvoke(rectbmp2, out MultiThreading9.existRect2, new AsyncCallback(MultiThreading9.CallbackTaskWithRect2), caller);
                    IAsyncResult result3 = caller.BeginInvoke(rectbmp3, out MultiThreading9.existRect3, new AsyncCallback(MultiThreading9.CallbackTaskWithRect3), caller);
                    IAsyncResult result4 = caller.BeginInvoke(rectbmp4, out MultiThreading9.existRect4, new AsyncCallback(MultiThreading9.CallbackTaskWithRect4), caller);
                    IAsyncResult result5 = caller.BeginInvoke(rectbmp5, out MultiThreading9.existRect5, new AsyncCallback(MultiThreading9.CallbackTaskWithRect5), caller);
                    IAsyncResult result6 = caller.BeginInvoke(rectbmp6, out MultiThreading9.existRect6, new AsyncCallback(MultiThreading9.CallbackTaskWithRect6), caller);
                    IAsyncResult result7 = caller.BeginInvoke(rectbmp7, out MultiThreading9.existRect7, new AsyncCallback(MultiThreading9.CallbackTaskWithRect7), caller);
                    IAsyncResult result8 = caller.BeginInvoke(rectbmp8, out MultiThreading9.existRect8, new AsyncCallback(MultiThreading9.CallbackTaskWithRect8), caller);
                    IAsyncResult result9 = caller.BeginInvoke(rectbmp9, out MultiThreading9.existRect9, new AsyncCallback(MultiThreading9.CallbackTaskWithRect9), caller);
                    while (!result1.IsCompleted || !result2.IsCompleted || !result3.IsCompleted || !result4.IsCompleted || !result5.IsCompleted || !result6.IsCompleted || !result7.IsCompleted || !result8.IsCompleted || !result9.IsCompleted)
                    {
                        //Application.DoEvents();
                    }
                    if (MultiThreading9.existRect1 || MultiThreading9.existRect2 || MultiThreading9.existRect3 || MultiThreading9.existRect4 || MultiThreading9.existRect5 || MultiThreading9.existRect6 || MultiThreading9.existRect7 || MultiThreading9.existRect8 || MultiThreading9.existRect9)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// 检测指定黑边(或白边)以threshold阀值在指定的路径下，是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="path"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="needFilter"></param>
        /// <param name="needEdgeDetect"></param>
        /// <param name="edgeThinning"></param>
        /// <returns></returns>
        public static bool IsExistKeyPointByPath9(Bitmap bmp, GraphicsPath path, int blkORwht, byte threshold, bool needGrayExt, bool needFilter, bool needEdgeDetect, bool edgeThinning)
        {
            lock (lockObj)
            {
                MultiThreading9.existPath1 = false;
                MultiThreading9.existPath2 = false;
                MultiThreading9.existPath3 = false;
                MultiThreading9.existPath4 = false;
                MultiThreading9.existPath5 = false;
                MultiThreading9.existPath6 = false;
                MultiThreading9.existPath7 = false;
                MultiThreading9.existPath8 = false;
                MultiThreading9.existPath9 = false;
                if (bmp != null)
                {
                    RectangleF rect = path.GetBounds();
                    Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(rect), new Rectangle(0, 0, bmp.Width, bmp.Height));
                    int X = rec.X;
                    int Y = rec.Y;

                    Bitmap bm = null;
                    if (needFilter)
                        bm = Filtering(bmp.Clone(rec, bmp.PixelFormat), blkORwht, threshold, 1);
                    else
                        bm = Bitize(bmp.Clone(rec, bmp.PixelFormat), threshold, needGrayExt);
                    if (needGrayExt)
                        bm = GrayExtend(bm, true);
                    if (needEdgeDetect)
                    {
                        bm = GetBitBmpByRoberts(bm, threshold, edgeThinning);
                    }

                    Rectangle rec1 = new Rectangle(0, 0, rec.Width / 3, rec.Height / 3);
                    Rectangle rec2 = new Rectangle(rec.Width / 3, 0, rec.Width / 3, rec.Height / 3);
                    Rectangle rec3 = new Rectangle(rec.Width * 2 / 3, 0, rec.Width - rec.Width * 2 / 3, rec.Height / 3);
                    Rectangle rec4 = new Rectangle(0, rec.Height / 3, rec.Width / 3, rec.Height / 3);
                    Rectangle rec5 = new Rectangle(rec.Width / 3, rec.Height / 3, rec.Width / 3, rec.Height / 3);
                    Rectangle rec6 = new Rectangle(rec.Width * 2 / 3, rec.Height / 3, rec.Width - rec.Width * 2 / 3, rec.Height / 3);
                    Rectangle rec7 = new Rectangle(0, rec.Height * 2 / 3, rec.Width / 3, rec.Height - rec.Height * 2 / 3);
                    Rectangle rec8 = new Rectangle(rec.Width / 3, rec.Height * 2 / 3, rec.Width / 3, rec.Height - rec.Height * 2 / 3);
                    Rectangle rec9 = new Rectangle(rec.Width * 2 / 3, rec.Height * 2 / 3, rec.Width - rec.Width * 2 / 3, rec.Height - rec.Height * 2 / 3);
                    RectBitmapWithPath rectbmp1 = new RectBitmapWithPath();
                    rectbmp1.blkORwht = blkORwht;
                    rectbmp1.threshold = threshold;
                    rectbmp1.path = (GraphicsPath)path.Clone();
                    rectbmp1.bmp = FastClipBitmap(bm, rec1);
                    rec1.X += X;
                    rec1.Y += Y;
                    rectbmp1.rect = rec1;
                    RectBitmapWithPath rectbmp2 = new RectBitmapWithPath();
                    rectbmp2.blkORwht = blkORwht;
                    rectbmp2.threshold = threshold;
                    rectbmp2.path = (GraphicsPath)path.Clone();
                    rectbmp2.bmp = FastClipBitmap(bm, rec2);
                    rec2.X += X;
                    rec2.Y += Y;
                    rectbmp2.rect = rec2;
                    RectBitmapWithPath rectbmp3 = new RectBitmapWithPath();
                    rectbmp3.blkORwht = blkORwht;
                    rectbmp3.threshold = threshold;
                    rectbmp3.path = (GraphicsPath)path.Clone();
                    rectbmp3.bmp = FastClipBitmap(bm, rec3);
                    rec3.X += X;
                    rec3.Y += Y;
                    rectbmp3.rect = rec3;
                    RectBitmapWithPath rectbmp4 = new RectBitmapWithPath();
                    rectbmp4.blkORwht = blkORwht;
                    rectbmp4.threshold = threshold;
                    rectbmp4.path = (GraphicsPath)path.Clone();
                    rectbmp4.bmp = FastClipBitmap(bm, rec4);
                    rec4.X += X;
                    rec4.Y += Y;
                    rectbmp4.rect = rec4;
                    RectBitmapWithPath rectbmp5 = new RectBitmapWithPath();
                    rectbmp5.blkORwht = blkORwht;
                    rectbmp5.threshold = threshold;
                    rectbmp5.path = (GraphicsPath)path.Clone();
                    rectbmp5.bmp = FastClipBitmap(bm, rec5);
                    rec5.X += X;
                    rec5.Y += Y;
                    rectbmp5.rect = rec5;
                    RectBitmapWithPath rectbmp6 = new RectBitmapWithPath();
                    rectbmp6.blkORwht = blkORwht;
                    rectbmp6.threshold = threshold;
                    rectbmp6.path = (GraphicsPath)path.Clone();
                    rectbmp6.bmp = FastClipBitmap(bm, rec6);
                    rec6.X += X;
                    rec6.Y += Y;
                    rectbmp6.rect = rec6;
                    RectBitmapWithPath rectbmp7 = new RectBitmapWithPath();
                    rectbmp7.blkORwht = blkORwht;
                    rectbmp7.threshold = threshold;
                    rectbmp7.path = (GraphicsPath)path.Clone();
                    rectbmp7.bmp = FastClipBitmap(bm, rec7);
                    rec7.X += X;
                    rec7.Y += Y;
                    rectbmp7.rect = rec7;
                    RectBitmapWithPath rectbmp8 = new RectBitmapWithPath();
                    rectbmp8.blkORwht = blkORwht;
                    rectbmp8.threshold = threshold;
                    rectbmp8.path = (GraphicsPath)path.Clone();
                    rectbmp8.bmp = FastClipBitmap(bm, rec8);
                    rec8.X += X;
                    rec8.Y += Y;
                    rectbmp8.rect = rec8;
                    RectBitmapWithPath rectbmp9 = new RectBitmapWithPath();
                    rectbmp9.blkORwht = blkORwht;
                    rectbmp9.threshold = threshold;
                    rectbmp9.path = (GraphicsPath)path.Clone();
                    rectbmp9.bmp = FastClipBitmap(bm, rec9);
                    rec9.X += X;
                    rec9.Y += Y;
                    rectbmp9.rect = rec9;
                    AsyncMethodCallerWithPath caller = new AsyncMethodCallerWithPath(MultiThreading9.CheckSomeRectWithPath);
                    IAsyncResult result1 = caller.BeginInvoke(rectbmp1, out MultiThreading9.existPath1, new AsyncCallback(MultiThreading9.CallbackTaskWithPath1), caller);
                    IAsyncResult result2 = caller.BeginInvoke(rectbmp2, out MultiThreading9.existPath2, new AsyncCallback(MultiThreading9.CallbackTaskWithPath2), caller);
                    IAsyncResult result3 = caller.BeginInvoke(rectbmp3, out MultiThreading9.existPath3, new AsyncCallback(MultiThreading9.CallbackTaskWithPath3), caller);
                    IAsyncResult result4 = caller.BeginInvoke(rectbmp4, out MultiThreading9.existPath4, new AsyncCallback(MultiThreading9.CallbackTaskWithPath4), caller);
                    IAsyncResult result5 = caller.BeginInvoke(rectbmp5, out MultiThreading9.existPath5, new AsyncCallback(MultiThreading9.CallbackTaskWithPath5), caller);
                    IAsyncResult result6 = caller.BeginInvoke(rectbmp6, out MultiThreading9.existPath6, new AsyncCallback(MultiThreading9.CallbackTaskWithPath6), caller);
                    IAsyncResult result7 = caller.BeginInvoke(rectbmp7, out MultiThreading9.existPath7, new AsyncCallback(MultiThreading9.CallbackTaskWithPath7), caller);
                    IAsyncResult result8 = caller.BeginInvoke(rectbmp8, out MultiThreading9.existPath8, new AsyncCallback(MultiThreading9.CallbackTaskWithPath8), caller);
                    IAsyncResult result9 = caller.BeginInvoke(rectbmp9, out MultiThreading9.existPath9, new AsyncCallback(MultiThreading9.CallbackTaskWithPath9), caller);
                    while (!result1.IsCompleted || !result2.IsCompleted || !result3.IsCompleted || !result4.IsCompleted || !result5.IsCompleted || !result6.IsCompleted || !result7.IsCompleted || !result8.IsCompleted || !result9.IsCompleted)
                    {
                        //Application.DoEvents();
                    }
                    if (MultiThreading9.existPath1 || MultiThreading9.existPath2 || MultiThreading9.existPath3 || MultiThreading9.existPath4 || MultiThreading9.existPath5 || MultiThreading9.existPath6 || MultiThreading9.existPath7 || MultiThreading9.existPath8 || MultiThreading9.existPath9)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// 检测指定黑边(或白边)以threshold阀值在指定的区域中，是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="region"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="needFilter"></param>
        /// <param name="needEdgeDetect"></param>
        /// <param name="edgeThinning"></param>
        /// <returns></returns>
        public static bool IsExistKeyPointByRegion9(Bitmap bmp, Region region, int blkORwht, byte threshold, bool needGrayExt, bool needFilter, bool needEdgeDetect, bool edgeThinning)
        {
            lock (lockObj)
            {
                MultiThreading9.existRegion1 = false;
                MultiThreading9.existRegion2 = false;
                MultiThreading9.existRegion3 = false;
                MultiThreading9.existRegion4 = false;
                MultiThreading9.existRegion5 = false;
                MultiThreading9.existRegion6 = false;
                MultiThreading9.existRegion7 = false;
                MultiThreading9.existRegion8 = false;
                MultiThreading9.existRegion9 = false;
                if (bmp != null)
                {
                    Region reg = region.Clone();
                    reg.Intersect(new Rectangle(0, 0, bmp.Width, bmp.Height));
                    Graphics g = Graphics.FromImage(bmp);
                    Rectangle rec = Rectangle.Truncate(reg.GetBounds(g));
                    int X = rec.X;
                    int Y = rec.Y;

                    Bitmap bm = null;
                    if (needFilter)
                        bm = Filtering(bmp.Clone(rec, bmp.PixelFormat), blkORwht, threshold, 1);
                    else
                        bm = Bitize(bmp.Clone(rec, bmp.PixelFormat), threshold, needGrayExt);
                    if (needGrayExt)
                        bm = GrayExtend(bm, true);
                    if (needEdgeDetect)
                    {
                        bm = GetBitBmpByRoberts(bm, threshold, edgeThinning);
                    }

                    Rectangle rec1 = new Rectangle(0, 0, rec.Width / 3, rec.Height / 3);
                    Rectangle rec2 = new Rectangle(rec.Width / 3, 0, rec.Width / 3, rec.Height / 3);
                    Rectangle rec3 = new Rectangle(rec.Width * 2 / 3, 0, rec.Width - rec.Width * 2 / 3, rec.Height / 3);
                    Rectangle rec4 = new Rectangle(0, rec.Height / 3, rec.Width / 3, rec.Height / 3);
                    Rectangle rec5 = new Rectangle(rec.Width / 3, rec.Height / 3, rec.Width / 3, rec.Height / 3);
                    Rectangle rec6 = new Rectangle(rec.Width * 2 / 3, rec.Height / 3, rec.Width - rec.Width * 2 / 3, rec.Height / 3);
                    Rectangle rec7 = new Rectangle(0, rec.Height * 2 / 3, rec.Width / 3, rec.Height - rec.Height * 2 / 3);
                    Rectangle rec8 = new Rectangle(rec.Width / 3, rec.Height * 2 / 3, rec.Width / 3, rec.Height - rec.Height * 2 / 3);
                    Rectangle rec9 = new Rectangle(rec.Width * 2 / 3, rec.Height * 2 / 3, rec.Width - rec.Width * 2 / 3, rec.Height - rec.Height * 2 / 3);
                    RectBitmapWithRegion rectbmp1 = new RectBitmapWithRegion();
                    rectbmp1.blkORwht = blkORwht;
                    rectbmp1.threshold = threshold;
                    rectbmp1.region = reg.Clone();
                    rectbmp1.bmp = FastClipBitmap(bm, rec1);
                    rec1.X += X;
                    rec1.Y += Y;
                    rectbmp1.rect = rec1;
                    RectBitmapWithRegion rectbmp2 = new RectBitmapWithRegion();
                    rectbmp2.blkORwht = blkORwht;
                    rectbmp2.threshold = threshold;
                    rectbmp2.region = reg.Clone();
                    rectbmp2.bmp = FastClipBitmap(bm, rec2);
                    rec2.X += X;
                    rec2.Y += Y;
                    rectbmp2.rect = rec2;
                    RectBitmapWithRegion rectbmp3 = new RectBitmapWithRegion();
                    rectbmp3.blkORwht = blkORwht;
                    rectbmp3.threshold = threshold;
                    rectbmp3.region = reg.Clone();
                    rectbmp3.bmp = FastClipBitmap(bm, rec3);
                    rec3.X += X;
                    rec3.Y += Y;
                    rectbmp3.rect = rec3;
                    RectBitmapWithRegion rectbmp4 = new RectBitmapWithRegion();
                    rectbmp4.blkORwht = blkORwht;
                    rectbmp4.threshold = threshold;
                    rectbmp4.region = reg.Clone();
                    rectbmp4.bmp = FastClipBitmap(bm, rec4);
                    rec4.X += X;
                    rec4.Y += Y;
                    rectbmp4.rect = rec4;
                    RectBitmapWithRegion rectbmp5 = new RectBitmapWithRegion();
                    rectbmp5.blkORwht = blkORwht;
                    rectbmp5.threshold = threshold;
                    rectbmp5.region = reg.Clone();
                    rectbmp5.bmp = FastClipBitmap(bm, rec5);
                    rec5.X += X;
                    rec5.Y += Y;
                    rectbmp5.rect = rec5;
                    RectBitmapWithRegion rectbmp6 = new RectBitmapWithRegion();
                    rectbmp6.blkORwht = blkORwht;
                    rectbmp6.threshold = threshold;
                    rectbmp6.region = reg.Clone();
                    rectbmp6.bmp = FastClipBitmap(bm, rec6);
                    rec6.X += X;
                    rec6.Y += Y;
                    rectbmp6.rect = rec6;
                    RectBitmapWithRegion rectbmp7 = new RectBitmapWithRegion();
                    rectbmp7.blkORwht = blkORwht;
                    rectbmp7.threshold = threshold;
                    rectbmp7.region = reg.Clone();
                    rectbmp7.bmp = FastClipBitmap(bm, rec7);
                    rec7.X += X;
                    rec7.Y += Y;
                    rectbmp7.rect = rec7;
                    RectBitmapWithRegion rectbmp8 = new RectBitmapWithRegion();
                    rectbmp8.blkORwht = blkORwht;
                    rectbmp8.threshold = threshold;
                    rectbmp8.region = reg.Clone();
                    rectbmp8.bmp = FastClipBitmap(bm, rec8);
                    rec8.X += X;
                    rec8.Y += Y;
                    rectbmp8.rect = rec8;
                    RectBitmapWithRegion rectbmp9 = new RectBitmapWithRegion();
                    rectbmp9.blkORwht = blkORwht;
                    rectbmp9.threshold = threshold;
                    rectbmp9.region = reg.Clone();
                    rectbmp9.bmp = FastClipBitmap(bm, rec9);
                    rec9.X += X;
                    rec9.Y += Y;
                    rectbmp4.rect = rec9;
                    AsyncMethodCallerWithRegion caller = new AsyncMethodCallerWithRegion(MultiThreading9.CheckSomeRectWithRegion);
                    IAsyncResult result1 = caller.BeginInvoke(rectbmp1, out MultiThreading9.existRegion1, new AsyncCallback(MultiThreading9.CallbackTaskWithRegion1), caller);
                    IAsyncResult result2 = caller.BeginInvoke(rectbmp2, out MultiThreading9.existRegion2, new AsyncCallback(MultiThreading9.CallbackTaskWithRegion2), caller);
                    IAsyncResult result3 = caller.BeginInvoke(rectbmp3, out MultiThreading9.existRegion3, new AsyncCallback(MultiThreading9.CallbackTaskWithRegion3), caller);
                    IAsyncResult result4 = caller.BeginInvoke(rectbmp4, out MultiThreading9.existRegion4, new AsyncCallback(MultiThreading9.CallbackTaskWithRegion4), caller);
                    IAsyncResult result5 = caller.BeginInvoke(rectbmp5, out MultiThreading9.existRegion5, new AsyncCallback(MultiThreading9.CallbackTaskWithRegion5), caller);
                    IAsyncResult result6 = caller.BeginInvoke(rectbmp6, out MultiThreading9.existRegion6, new AsyncCallback(MultiThreading9.CallbackTaskWithRegion6), caller);
                    IAsyncResult result7 = caller.BeginInvoke(rectbmp7, out MultiThreading9.existRegion7, new AsyncCallback(MultiThreading9.CallbackTaskWithRegion7), caller);
                    IAsyncResult result8 = caller.BeginInvoke(rectbmp8, out MultiThreading9.existRegion8, new AsyncCallback(MultiThreading9.CallbackTaskWithRegion8), caller);
                    IAsyncResult result9 = caller.BeginInvoke(rectbmp9, out MultiThreading9.existRegion9, new AsyncCallback(MultiThreading9.CallbackTaskWithRegion9), caller);
                    reg.Dispose();
                    while (!result1.IsCompleted || !result2.IsCompleted || !result3.IsCompleted || !result4.IsCompleted || !result5.IsCompleted || !result6.IsCompleted || !result7.IsCompleted || !result8.IsCompleted || !result9.IsCompleted)
                    {
                        //Application.DoEvents();
                    }
                    if (MultiThreading9.existRegion1 || MultiThreading9.existRegion2 || MultiThreading9.existRegion3 || MultiThreading9.existRegion4 || MultiThreading9.existRegion5 || MultiThreading9.existRegion6 || MultiThreading9.existRegion7 || MultiThreading9.existRegion8 || MultiThreading9.existRegion9)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// 检测指定黑边(或白边)以threshold阀值在指定的路径下，是否存在连通路径，并返回连通的路径
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="path"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="isCZ"></param>
        /// <param name="linkedRoad"></param>
        /// <returns></returns>
        public static bool IsExistLinkedRoadByPath(Bitmap bmp, GraphicsPath path, int blkORwht, byte threshold, bool needGrayExt, bool isCZ, out List<Point> linkedRoad)
        {
            ThreadCommon.linked1 = false;
            ThreadCommon.linked2 = false;
            ThreadCommon.linked3 = false;
            linkedRoad = null;
            if (bmp != null)
            {
                RectangleF rect = path.GetBounds();
                Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(rect), new Rectangle(0, 0, bmp.Width, bmp.Height));
                int X = rec.X;
                int Y = rec.Y;
                Bitmap bm = bmp.Clone(rec, bmp.PixelFormat);
                if (needGrayExt)
                    bm = GrayExtend(bm, false);//灰度拉伸

                RectBitmapWithPath rectbmp1 = new RectBitmapWithPath();
                rectbmp1.blkORwht = blkORwht;
                rectbmp1.threshold = threshold;
                rectbmp1.path = (GraphicsPath)path.Clone();
                rectbmp1.bmp = (Bitmap)bm.Clone();
                rectbmp1.rect = rec;
                RectBitmapWithPath rectbmp2 = new RectBitmapWithPath();
                rectbmp2.blkORwht = blkORwht;
                rectbmp2.threshold = threshold;
                rectbmp2.path = (GraphicsPath)path.Clone();
                rectbmp2.bmp = (Bitmap)bm.Clone();
                rectbmp2.rect = rec;
                RectBitmapWithPath rectbmp3 = new RectBitmapWithPath();
                rectbmp3.blkORwht = blkORwht;
                rectbmp3.threshold = threshold;
                rectbmp3.path = (GraphicsPath)path.Clone();
                rectbmp3.bmp = (Bitmap)bm.Clone();
                rectbmp3.rect = rec;
                AsyncMethodCallerLinkedWithPath caller = new AsyncMethodCallerLinkedWithPath(ThreadCommon.CheckSomeRoadWithPath);
                IAsyncResult result2 = caller.BeginInvoke(rectbmp2, isCZ, 2, out ThreadCommon.linked2, new AsyncCallback(ThreadCommon.CallbackRoadWithPath2), caller);
                IAsyncResult result1 = caller.BeginInvoke(rectbmp1, isCZ, 1, out ThreadCommon.linked1, new AsyncCallback(ThreadCommon.CallbackRoadWithPath1), caller);
                IAsyncResult result3 = caller.BeginInvoke(rectbmp3, isCZ, 3, out ThreadCommon.linked3, new AsyncCallback(ThreadCommon.CallbackRoadWithPath3), caller);
                while (!result1.IsCompleted || !result2.IsCompleted || !result3.IsCompleted)
                {
                    //Application.DoEvents();
                }
                if (ThreadCommon.linked1 || ThreadCommon.linked2 || ThreadCommon.linked3)
                {
                    if (ThreadCommon.lps2.Linked)
                    {
                        linkedRoad = ThreadCommon.lps2.PointList;
                    }
                    else if (ThreadCommon.lps1.Linked)
                    {
                        linkedRoad = ThreadCommon.lps1.PointList;
                    }
                    else
                    {
                        linkedRoad = ThreadCommon.lps3.PointList;
                    }
                    return true;
                }
                else
                {
                    return false;
                }
            }
            return false;
        }

        /// <summary>
        /// 用轮廓提取的方法寻找所有合适路径
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="path"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="isCZ"></param>
        /// <param name="linkedRoad"></param>
        /// <returns></returns>
        public static bool IsExistLinkedWayByPath(Bitmap bmp, GraphicsPath path, int blkORwht, byte threshold, bool needGrayExt, bool isCZ, out List<Point> linkedRoad)
        {
            linkedRoad = null;
            RectangleF rect = path.GetBounds();
            Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(rect), new Rectangle(0, 0, bmp.Width, bmp.Height));
            Point recLoc = rec.Location;
            bool link = CheckWayLinked(bmp, path, blkORwht, threshold, needGrayExt, recLoc, isCZ, out linkedRoad);

            /*int X = rec.X;
            int Y = rec.Y;
            Bitmap bm = bmp.Clone(rec, bmp.PixelFormat);
            if (needGrayExt)
                bm = GrayExtend(bm, false);//灰度拉伸
            bool[,] ways = Image2Array(bm, blkORwht, threshold);
            tmpPoints.Clear();
            isLinkCZ = isCZ;
            
            MAX = 0;
            for (int j = 0; j < ways.GetLength(1); j++)
            {
                for (int i = 0; i < ways.GetLength(0); i++)
                {
                    if (ways[i, j])
                        MAX++;
                }
            }
            int enumber, cost = 0;
            Graph G = new Graph();
            G.vnum = MAX;
            G.e = new EDGE[MAX * (MAX - 1) / 2];
            enumber = CreatGraph(ways, ref G);
            List<EDGE> tree;
            cost = Kruskal(G, enumber, out tree);
            foreach (EDGE e in tree)
            {
                tmpPoints.Add(e.v1);
                tmpPoints.Add(e.v2);
            }*/
            /*
            if (isCZ)
            {
                for (int i = 0; i < rec.Width; i++)
                {
                    if (ways[i, 0])
                        FindWay1(ways, i, 0);
                }
            }
            else
            {
                for (int j = 0; j < rec.Height; j++)
                {
                    if (ways[0, j])
                        FindWay1(ways, 0, j);
                }
            }*/
            return link;
        }

        /// <summary>
        /// 确定水平或垂直方向的关键点路径跟踪
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="path"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExtend"></param>
        /// <param name="rec"></param>
        /// <param name="isCZ"></param>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static bool CheckWayLinked(Bitmap bmp, GraphicsPath path, int blkORwht, byte threshold, bool needGrayExtend, Point rec, bool isCZ, out List<Point> ps)
        {
            ps = new List<Point>();
            bool[,] b = Image2Array(FastClipBitmap(bmp, path), blkORwht, threshold);
            int width = b.GetLength(0);
            int height = b.GetLength(1);
            // 不考虑图像最外围一圈
            int topRect = 1;
            int bottomRect = height - 1;
            int leftRect = 1;
            int rightRect = width - 1;
            for (int y = topRect; y < bottomRect; y++)
            {
                for (int x = leftRect; x < rightRect; x++)
                {
                    if (!b[x, y])
                    {
                        bool sum =
                          b[x - 1, y - 1] && b[x, y - 1] && b[x + 1, y - 1] &&
                          b[x - 1, y] && b[x + 1, y] &&
                          b[x - 1, y + 1] && b[x, y + 1] && b[x + 1, y + 1];

                        // 如果周围八点全为关键点，则取当前点
                        if (sum)
                            ps.Add(new Point(x + rec.X, y + rec.Y));
                    }
                    else
                    {
                        ps.Add(new Point(x + rec.X, y + rec.Y));
                    }
                } // x
            } // y
            bool link = false;
            ushort[,] sign = ImageSign(bmp, path, blkORwht, threshold, needGrayExtend);

            // 找出最大标记号，即找出区域数
            int max = 0;
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    if (sign[x, y] > max)
                        max = sign[x, y];
                }
            }
            if (max > 0)
            {
                List<int> access = new List<int>();
                int maxC = 0;
                int pc = 0;
                for (int i = 1; i < max + 1; i++)
                {
                    pc = 0;
                    access.Clear();
                    for (int y = topRect; y < bottomRect; y++)
                    {
                        for (int x = leftRect; x < rightRect; x++)
                        {
                            Point pp = new Point(x, y);
                            ushort sig = sign[x, y];
                            if (sig == i)
                            {
                                if (isCZ)
                                {
                                    if (!access.Contains(y))
                                    {
                                        pc++;
                                        access.Add(y);
                                    }
                                }
                                else
                                {
                                    if (!access.Contains(x))
                                    {
                                        pc++;
                                        access.Add(x);
                                    }
                                }
                            }
                        }
                    }
                    if (maxC < pc)
                        maxC = pc;
                }
                if (isCZ)
                {
                    if (maxC > height - 3)
                    {
                        link = true;
                    }
                }
                else
                {
                    if (maxC > width - 3)
                    {
                        link = true;
                    }
                }
            }
            return link;
        } // end of CheckWayLinked

        /*
        struct EDGE
        {
            public Point v1, v2;          //每条边的两个顶点
            public int weight;         //边的权值
        }

        struct Graph
        {
            public int vnum;          //图的顶点数目
            public EDGE[] e;    //图中的边
        }

        struct Alist
        {
            public Point v;                 //用链表存储同一个连通分量的顶点
            public Point next;
        }

        static void HeapAdjust(EDGE[] data, int s, int m)
        { //将元素序列data[s..m]调整为小根堆,堆顶元素为data[s]
            int j;
            EDGE t;
            t = data[s];            //备份元素data[s],为其找到适当位置后插入
            for (j = 2 * s + 1; j <= m; j = j * 2 + 1)        //沿值较小的孩子结点向下筛选
            {
                if (j < m && data[j].weight > data[j + 1].weight) j++;
                if (!(t.weight > data[j].weight)) break;
                data[s] = data[j]; s = j;        //用s记录待插入元素的位置
            }
            data[s] = t;             //将备份元素插入由s所指出的位置
        }

        static int CreatGraph(bool[,] way, ref Graph p)  //输入图的顶点数和边数,创建一个图并返回图的边数
        {
            int m = p.e.Length, w = 1, k = 0;
            Point v1, v2;
            //printf("输入图的顶点数:");
            //scanf("%d",&n);
            if (p.vnum < 1) return 0;
            //printf("输入图的边数: ");
            //scanf("%d",&m);
            while (k < m)
            { //printf("输入第 %d 条边的两个顶点和权值:",k+1);
                //scanf("%d%d%d",&v1,&v2,&w);
                //if(v1>=0 && v1<n && v2>=0 && v2<n)

                v1 = new Point();
                v2 = new Point();
                p.e[k].v1 = v1;
                p.e[k].v2 = v2;
                p.e[k].weight = w;
                k++;
            }
            return m;
        }
        static Alist GetTheAlistInArray(Alist[] a, Point p)
        {
            foreach (Alist ae in a)
            {
                if (ae.v == p)
                    return ae;
            }
            return new Alist();
        }

        static void SetTheAlistInArray(ref Alist[] a, Point p, Alist b)
        {
            for (int i = 0; i < a.Length; i++)
            {
                if (a[i].v == p)
                {
                    a[i] = b;
                    return;
                }
            }
        }

        static int Kruskal(Graph G, int enumber, out List<EDGE> tree) //用kruskal算法求图G的最小生成树,返回其代价
        {
            tree = new List<EDGE>();
            int i, k, m, cost = 0;
            Point v1, v2;
            Alist p, q = new Alist();
            Alist[] a = new Alist[MAX];
            for (i = 0; i < G.vnum; i++)     //将每个连通分量的顶点存放在一个单链表中
            {
                a[i].v = G.e[i].v1;
                a[i].next = Point.Empty;
            }
            for (i = enumber - 1; i >= 0; i--)    //按边上的权值建立小根堆
                HeapAdjust(G.e, i, enumber - 1);
            k = G.vnum;                  //k用于计算图中连通分量的数目
            m = enumber - 1;
            i = 0;
            do
            {
                v1 = G.e[0].v1;
                v2 = G.e[0].v2;
                p = GetTheAlistInArray(a, v1);
                while (p.v != Point.Empty && p.v != v2)    //判断当前所选择边的顶点是否在同一个连通分量中
                {
                    q = p;
                    p.v = p.next;
                }
                if (p.v == Point.Empty)                  //如果不在同一个连通分量中
                {
                    p = q;
                    p.next = GetTheAlistInArray(a, G.e[0].v2).v;
                    p = GetTheAlistInArray(a, G.e[0].v1);       //加入边(v1,v2),将两个连通分量合并为一个
                    while (p.v != Point.Empty)
                    {
                        SetTheAlistInArray(ref a, p.v, GetTheAlistInArray(a, G.e[0].v1));
                        p.v = p.next;
                    }
                    k--;                  //连通分量数目减少一个
                    EDGE newTree = new EDGE();
                    newTree.v1 = v1;
                    newTree.v2 = v2;
                    newTree.weight = G.e[0].weight;
                    tree.Add(newTree);        //加入最小生成树中
                    cost += G.e[0].weight;
                    i++;
                }
                G.e[0] = G.e[m];
                m--;
                HeapAdjust(G.e, 0, m); //找下一条权值最小的边
            } while (k > 1);  //当所有的顶点不在同一个连通分量时,继续循环
            return cost;
        }

        static void FindWay1(bool[,] p, int x, int y)
        {
            int m = p.GetLength(0);
            int n = p.GetLength(1);
            int i, tx, ty;
            int[,] d = new int[2, 2] { { 0, 1 }, { 1, 0 } };
            Stack<int> vs = new Stack<int>();
            v.Push(x);
            v.Push(y);
            vs.Push(0);
            while (v.Count > 0)
            {
                int top = v.Pop();
                x = v.Peek();
                y = top;
                v.Push(top);
                if (x == m - 1 && y == n - 1)
                {
                    OutPut(p);
                    i = 2;
                }
                else
                {
                    for (i = vs.Peek(); i < 2; i++)
                    {
                        tx = x + d[i, 0];
                        ty = y + d[i, 1];
                        if (0 <= tx && tx < m && 0 <= ty && ty < n && !p[tx, ty])
                        {
                            p[tx, ty] = true;
                            vs.Pop();
                            vs.Push(i);
                            vs.Push(0);
                            v.Push(tx);
                            v.Push(ty);
                            break;
                        }
                    }
                }
                if (i >= 2)
                {
                    p[x, y] = false;
                    v.Pop();
                    v.Pop();
                    vs.Pop();
                    if (vs.Count > 0)
                    {
                        int peek = vs.Pop();
                        peek++;
                        vs.Push(peek);
                    }
                }
            }
        }
        */
        static void OutPut(bool[,] p)
        {
            //这里可以修改行列数,p[i][j]为1代表该点在路径上
            List<Point> ps = new List<Point>();
            int i = 0, j = 0;
            int m = p.GetLength(0) - 1;
            int n = p.GetLength(1) - 1;
            if (isLinkCZ)
            {
                while (j < n)
                {
                    if (i == m)
                    {
                        j++;
                    }
                    else if (j == n)
                    {
                        break;
                    }
                    else
                    {
                        if (p[i, j + 1])
                            j++;
                        else
                        {
                            if (p[i + 1, j])
                                i++;
                            else if (p[i - 1, j])
                                i--;
                        }
                    }
                    Point pp = new Point(i, j);
                    if (!tmpPoints.Contains(pp))
                        tmpPoints.Add(pp);
                }
            }
            else
            {
                while (i < m)
                {
                    if (i == m)
                    {
                        break;
                    }
                    else if (j == n)
                    {
                        i++;
                    }
                    else
                    {
                        if (p[i + 1, j])
                            i++;
                        else
                        {
                            if (p[i, j + 1])
                                j++;
                            else if (p[i, j - 1])
                                j--;
                        }
                    }
                    Point pp = new Point(i, j);
                    if (!tmpPoints.Contains(pp))
                        tmpPoints.Add(pp);
                }
            }
        }

        /// <summary>
        /// 用次方法之前要先置空tmpPoints
        /// </summary>
        /// <param name="p"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        static void FindWay(bool[,] p, int i, int j)//递归函数,寻找所有合适路径   
        {
            int m = p.GetLength(0) - 1;
            int n = p.GetLength(1) - 1;
            if ((i == m - 1) && (j == n - 1))
            {
                OutPut(p);
            }
            else if (i == m - 1)
            {
                p[i, j + 1] = true;
                p[i + 1, j] = false;
                FindWay(p, i, j + 1);
            }     //到底了,只能向右   
            else if (j == n - 1)
            {
                p[i + 1, j] = true;
                p[i, j + 1] = false;
                FindWay(p, i + 1, j);
            }     //最右了,只能向下   
            else
            {
                p[i, j + 1] = true;
                p[i + 1, j] = false;
                FindWay(p, i, j + 1);
                p[i + 1, j] = true;
                p[i, j + 1] = false;
                FindWay(p, i + 1, j);
            }
        }
        
        /// <summary>
        /// 判断给定的路径是否为正矩形
        /// </summary>
        /// <param name="path">给定的路径</param>
        /// <returns></returns>
        public static bool PathIsRect(GraphicsPath path)
        {
            if (path.PointCount == 5 && (path.PathPoints[0].X == path.PathPoints[1].X && path.PathPoints[1].Y == path.PathPoints[2].Y && path.PathPoints[2].X == path.PathPoints[3].X && path.PathPoints[3].Y == path.PathPoints[0].Y) || (path.PathPoints[0].Y == path.PathPoints[1].Y && path.PathPoints[1].X == path.PathPoints[2].X && path.PathPoints[2].Y == path.PathPoints[3].Y && path.PathPoints[3].X == path.PathPoints[0].X))
                return true;
            else
                return false;
        }
        /// <summary>
        /// 平面中两点的距离
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static double dist(Point2 p1, Point2 p2)
        {
            return Math.Sqrt((p2.X - p1.X) * (p2.X - p1.X) + (p2.Y - p1.Y) * (p2.Y - p1.Y));
        }
        /// <summary>
        /// 返回两个矢量l1和l2的夹角的余弦(-1 --- 1)注意：如果想从余弦求夹角的话，注意反余弦函数的定义域是从 0到pi
        /// </summary>
        /// <param name="l1"></param>
        /// <param name="l2"></param>
        /// <returns></returns>
        public static double cosine(Line2 l1, Line2 l2)
        {
            double dis1 = dist(l1.Point2, l1.Point1);
            double dis2 = dist(l2.Point2, l2.Point1);
            if (dis1 * dis2 > Common.Eps)
                return (((l1.Point2.X - l1.Point1.X) * (l2.Point2.X - l2.Point1.X) + (l1.Point2.Y - l1.Point1.Y) * (l2.Point2.Y - l2.Point1.Y)) / (dis1 * dis2));
            else
                return double.MaxValue;
        }
        //r=Multiply(sp,ep,op),得到(sp-op)*(ep-op)的叉积
        //r>0:ep在向量op-sp的逆时针方向；
        //r=0:op-sp-ep三点共线；
        //r<0:ep在向量op-sp的顺时针方向
        /// <summary>
        /// 求平面向量的叉积，并可以根据返回值判断三点的位置关系
        /// r》0:ep在向量op-sp的逆时针方向；
        /// r==0:op-sp-ep三点共线；
        /// r《0:ep在向量op-sp的顺时针方向
        /// </summary>
        /// <param name="sp"></param>
        /// <param name="ep"></param>
        /// <param name="op"></param>
        /// <returns></returns>
        public static double Multiply(Point2 sp, Point2 ep, Point2 op)
        {
            return ((sp.X - op.X) * (ep.Y - op.Y) - (ep.X - op.X) * (sp.Y - op.Y));
        }
        /// <summary>
        /// 获取点集的Graham扫描凸包壳，返回凸包上的点集按照逆时针方向排列（没问题！只是会忽略一个点）
        /// </summary>
        /// <param name="PointSet"></param>
        public static Point2[] Graham_Scan(Point2[] PointSet)
        {
            if (PointSet == null || PointSet.Length == 0)
                return null;
            List<Point2> ch = new List<Point2>();
            int i, j, k = 0, top = 2;
            int n = PointSet.Length;
            Point2 tmp;
            // 选取PointSet中y坐标最小的点PointSet[k]，如果这样的点有多个，则取最左边的一个

            for (i = 1; i < n; i++)
                if (PointSet[i].Y < PointSet[k].Y || (PointSet[i].Y == PointSet[k].Y) && (PointSet[i].X < PointSet[k].X))
                {
                    k = i;
                }
            tmp = PointSet[0];
            PointSet[0] = PointSet[k];
            PointSet[k] = tmp; // 现在PointSet中y坐标最小的点在PointSet[0]
            for (i = 1; i < n - 1; i++) /* 对顶点按照相对PointSet[0]的极角从小到大进行排序，极角相
同
的按照距离PointSet[0]从近到远进行排序 */
            {
                k = i;
                for (j = i + 1; j < n; j++)
                {
                    if (Multiply(PointSet[j], PointSet[k], PointSet[0]) > 0 ||  // 极角更小   
                          (Multiply(PointSet[j], PointSet[k], PointSet[0]) == 0) && /* 极角相等，距离更短 */
                                    dist(PointSet[0], PointSet[j]) < dist(PointSet[0], PointSet[k]))
                    {
                        k = j;
                    }
                }
                tmp = PointSet[i];
                PointSet[i] = PointSet[k];
                PointSet[k] = tmp;
            }
            ch.Add(PointSet[0]);
            ch.Add(PointSet[1]);
            ch.Add(PointSet[2]);
            int len = 2;
            for (i = 3; i < n; i++)
            {
                while (top > 0 && top < ch.Count && Multiply(PointSet[i], ch[top], ch[top - 1]) >= 0)
                {
                    top--;
                }
                ch.Add(PointSet[i]);
                len++;
            }
            ch.RemoveAt(len);
            return ch.ToArray();
        }
        /// <summary>
        /// 卷包裹法求点集凸壳（有问题...）
        /// </summary>
        /// <param name="PointSet"></param>
        public static Point2[] ConvexClosure(Point2[] PointSet)
        {
            if (PointSet == null || PointSet.Length == 0)
                return null;
            List<Point2> ch = new List<Point2>();
            int top = 0, i, index, first;
            int n = PointSet.Length;
            double curmax, curcos, curdis;
            Point2 tmp;
            Line2 l1 = new Line2(), l2 = new Line2();
            bool[] use = new bool[n];
            tmp = PointSet[0];
            index = 0;
            // 选取y最小点，如果多于一个，则选取最左点
            for (i = 1; i < n; i++)
            {
                if (PointSet[i].Y < tmp.Y || PointSet[i].Y == tmp.Y && PointSet[i].X < tmp.X)
                {
                    index = i;
                }
                use[i] = false;
            }
            tmp = PointSet[index];
            first = index;
            use[index] = true;

            index = -1;
            ch.Add(tmp);
            top++;
            tmp.X -= 100;
            l1.Point1 = tmp;
            l1.Point2 = ch[0];
            l2.Point1 = ch[0];

            while (index != first)
            {
                curmax = -100;
                curdis = 0;
                // 选取与最后一条确定边夹角最小的点，即余弦值最大者
                for (i = 0; i < n; i++)
                {
                    if (use[i])
                        continue;
                    l2.Point2 = PointSet[i];
                    curcos = cosine(l1, l2); // 根据cos值求夹角余弦，范围在 （-1 -- 1 ）
                    if (curcos > curmax || Math.Abs(curcos - curmax) < Common.Eps && dist(l2.Point1, l2.Point2) > curdis)
                    {
                        curmax = curcos;
                        index = i;
                        curdis = dist(l2.Point1, l2.Point2);
                    }

                }
                use[first] = false;//清空第first个顶点标志，使最后能形成封闭的hull

                use[index] = true;
                ch.Add(PointSet[index]);
                top++;
                l1.Point1 = ch[top - 2];
                l1.Point2 = ch[top - 1];
                l2.Point1 = ch[top - 1];
            }
            ch.RemoveAt(top - 1);
            return ch.ToArray();
        }

        /// <summary>
        /// 寻找点集中Y极小、X极大的点
        /// </summary>
        /// <param name="ps">点集</param>
        /// <returns></returns>
        public Point FindTheMinYPoint(List<Point> ps)
        {
            Point minY = ps[0];
            int miny = ps[0].Y;
            int maxx = ps[0].X;
            for (int i = 1; i < ps.Count; i++)
            {
                if (miny > ps[i].Y)
                {
                    maxx = ps[i].X;
                    miny = ps[i].Y;
                    minY = ps[i];
                }
                else if (miny == ps[i].Y)
                {
                    if (maxx < ps[i].X)
                    {
                        maxx = ps[i].X;
                        minY = ps[i];
                    }
                }
            }
            return minY;
        }

        /// <summary>
        /// 寻找点集中Y极小、X极大的点
        /// </summary>
        /// <param name="ps">点集</param>
        /// <returns></returns>
        public PointF FindTheMinYPoint(List<PointF> ps)
        {
            PointF minY = PointF.Empty;
            float miny = ps[0].Y;
            float maxx = ps[0].X;
            for (int i = 1; i < ps.Count; i++)
            {
                if (miny > ps[i].Y)
                {
                    maxx = ps[i].X;
                    miny = ps[i].Y;
                    minY = ps[i];
                }
                else if (miny == ps[i].Y)
                {
                    if (maxx < ps[i].X)
                    {
                        maxx = ps[i].X;
                        minY = ps[i];
                    }
                }
            }
            return minY;
        }

        /// <summary>
        /// 返回顶角在o点，起始边为os，终止边为oe的夹角(单位：弧度)
        /// </summary>
        /// <param name="o"></param>
        /// <param name="s"></param>
        /// <param name="e"></param>
        /// <returns></returns>
        public double GetAngleValue(Point o, Point s, Point e)
        {
            double cosfi, fi, norm;
            double dsx = s.X - o.X;
            double dsy = s.Y - o.Y;
            double dex = e.X - o.X;
            double dey = e.Y - o.Y;

            cosfi = dsx * dex + dsy * dey;
            norm = (dsx * dsx + dey * dey) * (dex * dex + dey * dey);
            cosfi /= Math.Sqrt(norm);

            if (cosfi >= 1.0) return 0;
            if (cosfi <= -1.0) return -Math.PI;

            fi = Math.Acos(cosfi);
            if (dsx * dey - dsy * dex > 0) return fi;// 说明矢量os 在矢量 oe的顺时针方向
            return -fi;
        }

        /// <summary>
        /// 返回顶角在o点，起始边为os，终止边为oe的夹角(单位：弧度)
        /// </summary>
        /// <param name="o"></param>
        /// <param name="s"></param>
        /// <param name="e"></param>
        /// <returns></returns>
        public double GetAngleValue(PointF o, PointF s, PointF e)
        {
            double cosfi, fi, norm;
            double dsx = s.X - o.X;
            double dsy = s.Y - o.Y;
            double dex = e.X - o.X;
            double dey = e.Y - o.Y;

            cosfi = dsx * dex + dsy * dey;
            norm = (dsx * dsx + dey * dey) * (dex * dex + dey * dey);
            cosfi /= Math.Sqrt(norm);

            if (cosfi >= 1.0) return 0;
            if (cosfi <= -1.0) return -Math.PI;

            fi = Math.Acos(cosfi);
            if (dsx * dey - dsy * dex > 0) return fi;// 说明矢量os 在矢量 oe的顺时针方向
            return -fi;
        }

        private List<Point> GetOutSideOfPoints(List<Point> ps, bool realOutSide)
        {
            List<Point> result = new List<Point>();
            Point minY = FindTheMinYPoint(ps);
            result.Add(minY);
            Point OriginX = new Point(minY.X + 10, minY.Y);
            List<LAD> lads = new List<LAD>();
            for (int i = 0; i < ps.Count; i++)
            {
                if (ps[i].X != minY.X || ps[i].Y != minY.Y)
                {
                    LAD lad;
                    lad.Location = ps[i];
                    lad.Angle = GetAngleValue(minY, OriginX, ps[i]);
                    lad.Distance = Distance2D(minY, ps[i]);
                    lads.Add(lad);
                }
            }
            //先做角度排序
            lads.Sort(new AngleSort());
            List<int> start = new List<int>();
            List<int> end = new List<int>();
            double lastAngle = lads[0].Angle;
            bool hasStart = false;
            for (int i = 1; i < lads.Count; i++)
            {
                //当多点同角时，比较距离
                if (!hasStart && lastAngle == lads[i].Angle)
                {
                    start.Add(i - 1);
                    hasStart = true;
                }
                else if (hasStart && lastAngle != lads[i].Angle)
                {
                    end.Add(i - 1);
                    hasStart = false;
                }
                lastAngle = lads[i].Angle;
            }
            for (int i = 0; i < start.Count; i++)
            {
                int count = end[i] - start[i];
                List<LAD> lads1 = new List<LAD>();
                LAD[] tmp = new LAD[count];
                lads.CopyTo(start[i], tmp, 0, count);
                for (int j = 0; j < tmp.Length; j++)
                {
                    lads1.Add(tmp[j]);
                }
                //再做距离排序
                lads1.Sort(new DistanceSort());
                lads.RemoveRange(start[i], count);
                lads.Insert(start[i], lads1[0]);
            }
            for (int i = 0; i < lads.Count; i++)
            {
                result.Add(lads[i].Location);
            }
            if (realOutSide)
            {
                RemoveAOP(ref result);
            }
            return result;
        }

        /// <summary>
        /// 移除凹点
        /// </summary>
        /// <param name="result"></param>
        private void RemoveAOP(ref List<Point> result)
        {
            if (result.Count > 2)
            {
                bool haveAOP = false;
                for (int i = 0; i < result.Count - 2; i++)
                {
                    double tmp = Multiply(result[i], result[(i + 2) % result.Count], result[(i + 1) % result.Count]);
                    if (tmp <= 0)
                    {
                        //i+1为凹点
                        result.RemoveAt((i + 1) % result.Count);
                    }
                }
                for (int i = 0; i < result.Count - 2; i++)
                {
                    double tmp = Multiply(result[i], result[(i + 2) % result.Count], result[(i + 1) % result.Count]);
                    if (tmp <= 0)
                    {
                        haveAOP = true;
                        break;
                    }
                }
                if (haveAOP)
                {
                    RemoveAOP(ref result);
                }
            }
        }

        /// <summary>
        /// 点角度距离结构
        /// </summary>
        public struct LAD
        {
            /// <summary>
            /// 点
            /// </summary>
            public Point Location;
            /// <summary>
            /// 角度
            /// </summary>
            public double Angle;
            /// <summary>
            /// 距离
            /// </summary>
            public double Distance;
            /// <summary>
            /// 构造函数
            /// </summary>
            /// <param name="Location"></param>
            /// <param name="Angle"></param>
            /// <param name="Distance"></param>
            public LAD(Point Location, double Angle, double Distance)
            {
                this.Location = Location;
                this.Angle = Angle;
                this.Distance = Distance;
            }
        }

        /// <summary>
        /// 角度排序类
        /// </summary>
        public class AngleSort : IComparer<LAD>
        {
            int IComparer<LAD>.Compare(LAD x, LAD y)
            {
                CaseInsensitiveComparer cic = new CaseInsensitiveComparer();
                return cic.Compare(y.Angle, x.Angle);
            }
        }
        /// <summary>
        /// 距离排序类
        /// </summary>
        public class DistanceSort : IComparer<LAD>
        {
            int IComparer<LAD>.Compare(LAD x, LAD y)
            {
                CaseInsensitiveComparer cic = new CaseInsensitiveComparer();
                return cic.Compare(y.Distance, x.Distance);
            }
        }

        /// <summary>
        /// r《0:ep在向量op-sp的逆时针方向，r==0:op-sp-ep三点共线，r》0:ep在向量op-sp的顺时针方向，注意方向是以小于180度的哪个方向为准，op为旋转中心
        /// </summary>
        /// <param name="sp"></param>
        /// <param name="ep"></param>
        /// <param name="op"></param>
        /// <returns></returns>
        public static double Multiply(Point sp, Point ep, Point op)
        {
            return ((sp.X - op.X) * (ep.Y - op.Y) - (ep.X - op.X) * (sp.Y - op.Y));
        }

        // 如果无特别说明，输入多边形顶点要求按逆时针排列
        /*
        返回值：输入的多边形是简单多边形，返回true
        要 求：输入顶点序列按逆时针排序
        说 明：简单多边形定义：
        1：循环排序中相邻线段对的交是他们之间共有的单个点
        2：不相邻的线段不相交
        本程序默认第一个条件已经满足
        */
        /// <summary>
        /// 输入的多边形是否为简单多边形
        /// </summary>
        /// <param name="polygon"></param>
        /// <returns></returns>
        public static bool IsSimple(Point[] polygon)
        {
            int cn;
            Line2D l1, l2;
            int vcount = polygon.Length;
            for (int i = 0; i < vcount; i++)
            {
                l1.Point1 = polygon[i];
                l1.Point2 = polygon[(i + 1) % vcount];
                cn = vcount - 3;
                while (cn > 0)
                {
                    l2.Point1 = polygon[(i + 2) % vcount];
                    l2.Point2 = polygon[(i + 3) % vcount];
                    if (intersect(l1, l2))
                        break;
                    cn--;
                }
                if (cn > 0)
                    return false;
            }
            return true;
        }

        /// <summary>
        /// 二维线结构
        /// </summary>
        public struct Line2D
        {
            /// <summary>
            /// 线上的一点
            /// </summary>
            public Point Point1;
            /// <summary>
            /// 线上的另一点
            /// </summary>
            public Point Point2;
        }

        /// <summary>
        /// 如果线段u和v相交(包括相交在端点处)时，返回true
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns> 
        public static bool intersect(Line2D u, Line2D v)
        {
            return ((Math.Max(u.Point1.X, u.Point2.X) >= Math.Min(v.Point1.X, v.Point2.X)) &&                     //排斥实验
                    (Math.Max(v.Point1.X, v.Point2.X) >= Math.Min(u.Point1.X, u.Point2.X)) &&
                    (Math.Max(u.Point1.Y, u.Point2.Y) >= Math.Min(v.Point1.Y, v.Point2.Y)) &&
                    (Math.Max(v.Point1.Y, v.Point2.Y) >= Math.Min(u.Point1.Y, u.Point2.Y)) &&
                    (Multiply(v.Point1, u.Point2, u.Point1) * Multiply(u.Point2, v.Point2, u.Point1) >= 0) &&         //跨立实验

                    (Multiply(u.Point1, v.Point2, v.Point1) * Multiply(v.Point2, u.Point2, v.Point1) >= 0));
        }

        /// <summary>
        /// 返回值：按输入顺序返回多边形顶点的凸凹性判断，bc[i]=true,iff:第i个顶点是凸顶点
        /// </summary>
        /// <param name="polygon"></param>
        /// <param name="bc"></param>
        private static void CheckConvex(Point[] polygon, out bool[] bc)
        {
            int index = 0;
            Point tp = polygon[0];
            int vcount = polygon.Length;
            for (int i = 1; i < vcount; i++) // 寻找凸顶点
            {
                if (polygon[i].Y < tp.Y || (polygon[i].Y == tp.Y && polygon[i].X < tp.X))
                {
                    tp = polygon[i];
                    index = i;
                }
            }
            int count = vcount;// - 1;
            bc = new bool[count];
            bc[index] = true;
            while (count > 0) // 判断凸凹性
            {
                if (index < polygon.Length)
                {
                    if (Multiply(polygon[(index + 1) % vcount], polygon[(index + 2) % vcount], polygon[index]) >= 0)
                        bc[(index + 1) % vcount] = true;
                    else
                        bc[(index + 1) % vcount] = false;
                    index++;
                    count--;
                }
            }
        }

        /// <summary>
        /// 返回值：多边形polygon是凸多边形时，返回true 
        /// </summary>
        /// <param name="polygon"></param>
        /// <returns></returns>
        public static bool IsConvex(Point[] polygon)
        {
            int vcount = polygon.Length;
            bool[] bc = new bool[vcount];
            CheckConvex(polygon, out bc);
            for (int i = 0; i < vcount; i++) // 逐一检查顶点，是否全部是凸顶点
                if (!bc[i])
                    return false;
            return true;
        }

        /// <summary>
        /// 返回多边形面积(signed)；输入顶点按逆时针排列时，返回正值；否则返回负值
        /// </summary>
        /// <param name="polygon"></param>
        /// <returns></returns>
        public static double Area_of_polygon(Point[] polygon)
        {
            int vcount = polygon.Length;
            double s;
            if (vcount < 3) return 0;
            s = polygon[0].Y * (polygon[vcount - 1].X - polygon[1].X);
            for (int i = 1; i < vcount; i++)
                s += polygon[i].Y * (polygon[(i - 1)].X - polygon[(i + 1) % vcount].X);
            return s / 2;//Math.Abs(s / 2);
        }

        /// <summary>
        /// 如果输入顶点按逆时针排列，返回true
        /// </summary>
        /// <param name="polygon"></param>
        /// <returns></returns>
        public static bool IsConterClock(Point[] polygon)
        {
            return Area_of_polygon(polygon) < 0;
        }

        /// <summary>
        /// 另一种判断多边形顶点排列方向的方法，逆时针为true
        /// </summary>
        /// <param name="polygon"></param>
        /// <returns></returns>
        public static bool IsCCwize(Point[] polygon)
        {
            int index;
            Point a, b, v;
            v = polygon[0];
            index = 0;
            int vcount = polygon.Length;
            for (int i = 1; i < vcount; i++) // 找到最低且最左顶点，肯定是凸顶点
            {
                if (polygon[i].Y < v.Y || polygon[i].Y == v.Y && polygon[i].X < v.X)
                {
                    v = polygon[i];
                    index = i;
                    break;
                }
            }
            a = polygon[(index - 1 + vcount) % vcount]; // 顶点v的前一顶点
            b = polygon[(index + 1) % vcount]; // 顶点v的后一顶点
            return Multiply(v, b, a) < 0;//逆时针为正
        }

        /// <summary>
        /// 检测指定黑边(或白边)以threshold阀值在指定的区域下的 凸包 所在Rectangle的中点(不处理透明/半透明像素)，并返回指定路径中去噪后的小位图
        /// </summary>
        /// <param name="rect"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public Point GetCenterPointByRectTB(Rectangle rect, int blkORwht, byte threshold)
        {
            int minX = int.MaxValue, minY = int.MaxValue, maxX = 0, maxY = 0;
            //List<Point2> psArr = new List<Point2>();
            List<Point> psArr = new List<Point>();
            if (bmp != null)
            {
                GraphicsUnit gu = GraphicsUnit.Pixel;
                Rectangle rec = Rectangle.Intersect(rect, Rectangle.Truncate(bmp.GetBounds(ref gu)));
                int X = rec.X;
                int Y = rec.Y;
                int Width = rec.Width + rec.X;
                int Height = rec.Height + rec.Y;
                //先滤波、去噪
                Bitmap b = Filtering(bmp.Clone(rec, bmp.PixelFormat), blkORwht, threshold, 1);

                if (b != null)
                {
                    unsafe
                    {
                        byte gray, R, G, B, A;
                        BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadOnly, b.PixelFormat);
                        byte* p = (byte*)data.Scan0;
                        int stride = data.Stride;
                        int offset = stride - BPP * b.Width;
                        //注意：由于滤波的原因造成最外一圈的不正常，因此不处理最外一圈
                        p += stride;//不处理最上面一行
                        for (int j = Y + 1; j < Height - 1; j++)
                        {
                            p += BPP;//不处理最左一列
                            for (int i = X + 1; i < Width - 1; i++)
                            {
                                A = p[3];
                                R = p[2];
                                G = p[1];
                                B = p[0];
                                if (A < 255)//不处理透明（半透明）像素，避免由于Color.Transparent[0, 255, 255, 255]透明色导致gray = B时的误判
                                    continue;
                                gray = B;//因为在去噪的时候已经转化为灰度图了
                                if (blkORwht == 0)//黑边
                                {
                                    if (gray < threshold)
                                    {
                                        /*if (minX > i)
                                            minX = i;
                                        if (maxX < i)
                                            maxX = i;
                                        if (minY > j)
                                            minY = j;
                                        if (maxY < j)
                                            maxY = j;*/
                                        psArr.Add(new Point(i, j));
                                    }
                                }
                                else//白边
                                {
                                    if (gray > threshold)
                                    {
                                        /*if (minX > i)
                                            minX = i;
                                        if (maxX < i)
                                            maxX = i;
                                        if (minY > j)
                                            minY = j;
                                        if (maxY < j)
                                            maxY = j;*/
                                        psArr.Add(new Point(i, j));
                                    }
                                }
                                p += BPP;
                            }
                            p += BPP;//不处理最右一列
                            p += offset;
                        }
                        //不处理最下面一行（丢弃）
                        b.UnlockBits(data);
                    }
                }
            }
            //Point2[] ps1 = Graham_Scan(psArr.ToArray()); //ConvexClosure(psArr.ToArray()); //Graham_Scan太慢
            List<Point> ps1 = GetOutSideOfPoints(psArr, true);
            if (ps1 != null)
            {
                foreach (Point p1 in ps1)
                {
                    if (minX > p1.X)
                        minX = (int)p1.X;
                    if (maxX < p1.X)
                        maxX = (int)p1.X;
                    if (minY > p1.Y)
                        minY = (int)p1.Y;
                    if (maxY < p1.Y)
                        maxY = (int)p1.Y;
                }
            }
            Point center = Point.Empty;
            if (minX != int.MaxValue && maxX != 0 && minY != int.MaxValue && maxY != 0)
            {
                //center = new Point((maxX + minX) / 2, (maxY + minY) / 2);
                Rectangle re = Rectangle.FromLTRB(minX, minY, maxX, maxY);
                center = new Point(re.X + re.Width / 2, re.Y + re.Height / 2);
            }
            return center;
        }

        /// <summary>
        /// 检测指定黑边(或白边)以threshold阀值在指定的路径下关键点集的中点(不处理透明/半透明像素)
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="path"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public static Point GetCenterPointByPath(Bitmap bmp, GraphicsPath path, int blkORwht, byte threshold)
        {
            int minX = int.MaxValue, minY = int.MaxValue, maxX = 0, maxY = 0;
            if (bmp != null)
            {
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                RectangleF rectF = path.GetBounds();
                Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(rectF), new Rectangle(0, 0, bmp.Width, bmp.Height));
                int X = rec.X;
                int Y = rec.Y;
                int Width = rec.Width + rec.X;
                int Height = rec.Height + rec.Y;

                //先滤波、去噪
                Bitmap b = Filtering(bmp.Clone(rec, bmp.PixelFormat), blkORwht, threshold, 1);

                if (b != null)
                {
                    unsafe
                    {
                        byte gray, R, G, B, A;
                        BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
                        byte* p = (byte*)data.Scan0;
                        int stride = data.Stride;
                        int offset = stride - BPP * b.Width;
                        //注意：由于滤波的原因造成最外一圈的不正常，因此不处理最外一圈
                        p += stride;//不处理最上面一行
                        for (int j = Y + 1; j < Height - 1; j++)
                        {
                            p += BPP;//不处理最左一列
                            for (int i = X + 1; i < Width - 1; i++)
                            {
                                /*if (PathIsRect(path))
                                {
                                    A = p[3];
                                    R = p[2];
                                    G = p[1];
                                    B = p[0];
                                    if (A < 255)//不处理透明（半透明）像素，避免由于Color.Transparent[0, 255, 255, 255]透明色导致gray = B时的误判
                                        continue;
                                    gray = B;//因为在去噪的时候已经转化为灰度图了
                                    if (blkORwht == 0)//黑边
                                    {
                                        if (gray < threshold)
                                        {
                                            if (minX > i)
                                                minX = i;
                                            if (maxX < i)
                                                maxX = i;
                                            if (minY > j)
                                                minY = j;
                                            if (maxY < j)
                                                maxY = j;
                                        }
                                    }
                                    else//白边
                                    {
                                        if (gray > threshold)
                                        {
                                            if (minX > i)
                                                minX = i;
                                            if (maxX < i)
                                                maxX = i;
                                            if (minY > j)
                                                minY = j;
                                            if (maxY < j)
                                                maxY = j;
                                        }
                                    }
                                }
                                else if (path.IsVisible(i, j))
                                {*/
                                if (path.IsVisible(i, j))
                                {
                                    A = p[3];
                                    R = p[2];
                                    G = p[1];
                                    B = p[0];
                                    if (A < 255)//不处理透明（半透明）像素，避免由于Color.Transparent[0, 255, 255, 255]透明色导致gray = B时的误判
                                        continue;
                                    gray = B;//因为在去噪的时候已经转化为灰度图了
                                    if (blkORwht == 0)//黑边
                                    {
                                        if (gray < threshold)
                                        {
                                            if (minX > i)
                                                minX = i;
                                            if (maxX < i)
                                                maxX = i;
                                            if (minY > j)
                                                minY = j;
                                            if (maxY < j)
                                                maxY = j;
                                        }
                                    }
                                    else//白边
                                    {
                                        if (gray > threshold)
                                        {
                                            if (minX > i)
                                                minX = i;
                                            if (maxX < i)
                                                maxX = i;
                                            if (minY > j)
                                                minY = j;
                                            if (maxY < j)
                                                maxY = j;
                                        }
                                    }
                                }
                                p += BPP;
                            }
                            p += BPP;//不处理最右一列
                            p += offset;
                        }
                        //不处理最下面一行（丢弃）
                        b.UnlockBits(data);
                    }
                }
            }
            Point center = Point.Empty;
            if (minX != int.MaxValue && maxX != 0 && minY != int.MaxValue && maxY != 0)
                center = new Point((maxX + minX) / 2, (maxY + minY) / 2);
            return center;
        }
        /// <summary>
        /// 检测指定黑边(或白边)以threshold阀值在指定的区域下关键点集的中点(不处理透明/半透明像素)，并返回指定路径中去噪后的小位图
        /// </summary>
        /// <param name="rect"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public Point GetCenterPointByRect(Rectangle rect, int blkORwht, byte threshold)
        {
            int minX = int.MaxValue, minY = int.MaxValue, maxX = 0, maxY = 0;
            if (bmp != null)
            {
                Rectangle rec = Rectangle.Intersect(rect, new Rectangle(0, 0, width, height));
                int X = rec.X;
                int Y = rec.Y;
                int Width = rec.Width + rec.X;
                int Height = rec.Height + rec.Y;
                //先滤波、去噪
                Bitmap b = Filtering(bmp.Clone(rec, bmp.PixelFormat), blkORwht, threshold, 1);

                if (b != null)
                {
                    unsafe
                    {
                        byte gray, R, G, B, A;
                        BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
                        byte* p = (byte*)data.Scan0;
                        int stride = data.Stride;
                        int offset = stride - BPP * b.Width;
                        //注意：由于滤波的原因造成最外一圈的不正常，因此不处理最外一圈
                        p += stride;//不处理最上面一行
                        for (int j = Y + 1; j < Height - 1; j++)
                        {
                            p += BPP;//不处理最左一列
                            for (int i = X + 1; i < Width - 1; i++)
                            {
                                A = p[3];
                                R = p[2];
                                G = p[1];
                                B = p[0];
                                if (A < 255)//不处理透明（半透明）像素，避免由于Color.Transparent[0, 255, 255, 255]透明色导致gray = B时的误判
                                    continue;
                                gray = B;//因为在去噪的时候已经转化为灰度图了
                                if (blkORwht == 0)//黑边
                                {
                                    if (gray < threshold)
                                    {
                                        if (minX > i)
                                            minX = i;
                                        if (maxX < i)
                                            maxX = i;
                                        if (minY > j)
                                            minY = j;
                                        if (maxY < j)
                                            maxY = j;
                                    }
                                }
                                else//白边
                                {
                                    if (gray > threshold)
                                    {
                                        if (minX > i)
                                            minX = i;
                                        if (maxX < i)
                                            maxX = i;
                                        if (minY > j)
                                            minY = j;
                                        if (maxY < j)
                                            maxY = j;
                                    }
                                }
                                p += BPP;
                            }
                            p += BPP;//不处理最右一列
                            p += offset;
                        }
                        //不处理最下面一行（丢弃）
                        b.UnlockBits(data);
                    }
                }
            }
            Point center = Point.Empty;
            if (minX != int.MaxValue && maxX != 0 && minY != int.MaxValue && maxY != 0)
                center = new Point((maxX + minX) / 2, (maxY + minY) / 2);
            return center;
        }

        /// <summary>
        /// 检测指定黑边(或白边)以threshold阀值在指定的路径下的 凸包 所在Rectangle的中点(不处理透明/半透明像素)，并返回指定路径中去噪后的小位图
        /// </summary>
        /// <param name="path"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public Point GetCenterPointByPathTB(GraphicsPath path, int blkORwht, byte threshold)
        {
            int minX = int.MaxValue, minY = int.MaxValue, maxX = 0, maxY = 0;
            //List<Point2> psArr = new List<Point2>();
            List<Point> psArr = new List<Point>();
            if (bmp != null)
            {
                RectangleF rectF = path.GetBounds();
                Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(rectF), new Rectangle(0, 0, width, height));
                int X = rec.X;
                int Y = rec.Y;
                int Width = rec.Width + rec.X;
                int Height = rec.Height + rec.Y;
                //先滤波、去噪
                Bitmap b = Filtering(bmp.Clone(rec, bmp.PixelFormat), blkORwht, threshold, 1);

                if (b != null)
                {
                    unsafe
                    {
                        byte gray, R, G, B, A;
                        BitmapData data = b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
                        byte* p = (byte*)data.Scan0;
                        int stride = data.Stride;
                        int offset = stride - BPP * b.Width;
                        //注意：由于滤波的原因造成最外一圈的不正常，因此不处理最外一圈
                        p += stride;//不处理最上面一行
                        for (int j = Y + 1; j < Height - 1; j++)
                        {
                            p += BPP;//不处理最左一列
                            for (int i = X + 1; i < Width - 1; i++)
                            {
                                /*if (PathIsRect(path))
                                {
                                    A = p[3];
                                    R = p[2];
                                    G = p[1];
                                    B = p[0];
                                    if (A < 255)//不处理透明（半透明）像素，避免由于Color.Transparent[0, 255, 255, 255]透明色导致gray = B时的误判
                                        continue;
                                    gray = B;//因为在去噪的时候已经转化为灰度图了
                                    if (blkORwht == 0)//黑边
                                    {
                                        if (gray < threshold)
                                        {
                                            psArr.Add(new Point2(i, j));
                                        }
                                    }
                                    else//白边
                                    {
                                        if (gray > threshold)
                                        {
                                            psArr.Add(new Point2(i, j));
                                        }
                                    }
                                }
                                else if (path.IsVisible(i, j))
                                {*/
                                if (path.IsVisible(i, j))
                                {
                                    A = p[3];
                                    R = p[2];
                                    G = p[1];
                                    B = p[0];
                                    if (A < 255)//不处理透明（半透明）像素，避免由于Color.Transparent[0, 255, 255, 255]透明色导致gray = B时的误判
                                        continue;
                                    gray = B;//因为在去噪的时候已经转化为灰度图了
                                    if (blkORwht == 0)//黑边
                                    {
                                        if (gray < threshold)
                                        {
                                            psArr.Add(new Point(i, j));
                                        }
                                    }
                                    else//白边
                                    {
                                        if (gray > threshold)
                                        {
                                            psArr.Add(new Point(i, j));
                                        }
                                    }
                                }
                                p += BPP;
                            }
                            p += BPP;//不处理最右一列
                            p += offset;
                        }
                        //不处理最下面一行（丢弃）
                        b.UnlockBits(data);
                    }
                }
            }
            //Point2[] ps1 = Graham_Scan(psArr.ToArray()); //ConvexClosure(psArr.ToArray());//Graham_Scan太慢
            List<Point> ps1 = GetOutSideOfPoints(psArr, true);
            if (ps1 != null)
            {
                foreach (Point p1 in ps1)
                {
                    if (minX > p1.X)
                        minX = (int)p1.X;
                    if (maxX < p1.X)
                        maxX = (int)p1.X;
                    if (minY > p1.Y)
                        minY = (int)p1.Y;
                    if (maxY < p1.Y)
                        maxY = (int)p1.Y;
                }
            }
            Point center = Point.Empty;
            if (minX != int.MaxValue && maxX != 0 && minY != int.MaxValue && maxY != 0)
            {
                //center = new Point((maxX + minX) / 2, (maxY + minY) / 2);
                Rectangle re = Rectangle.FromLTRB(minX, minY, maxX, maxY);
                center = new Point(re.X + re.Width / 2, re.Y + re.Height / 2);
            }
            return center;
        }
        /// <summary>
        /// 无方向的线检测（返回直线的两点），有问题...
        /// </summary>
        /// <param name="rgn"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public bool GetLineNoDirection(Region rgn, int blkORwht, byte threshold, bool needGrayExt, out PointF p1, out PointF p2)
        {
            List<PointF> tmp = new List<PointF>();
            p1 = PointF.Empty;
            p2 = PointF.Empty;
            if (bmp != null)
            {
                Graphics g = Graphics.FromImage((Image)bmp.Clone());
                RectangleF rect = rgn.GetBounds(g);
                Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(rect), new Rectangle(0, 0, width, height));
                int X = rec.X;
                int Y = rec.Y;
                int Width = rec.Width + rec.X;
                int Height = rec.Height + rec.Y;

                Bitmap b = Thresholding(bmp, rec, threshold, needGrayExt);
                Bitmap edgeImage = GetBmpByRoberts(b);
                unsafe
                {
                    BitmapData data = edgeImage.LockBits(new Rectangle(0, 0, rec.Width, rec.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
                    byte* p = (byte*)data.Scan0;
                    int stride = data.Stride;
                    int offset = stride - BPP * rec.Width;
                    for (int j = Y; j < Height; j++)
                    {
                        for (int i = X; i < Width; i++)
                        {
                            PointF pp = new PointF(i, j);
                            if (rgn.IsVisible(pp))
                            {
                                if (blkORwht == 0)
                                {
                                    if (p[0] == 0)
                                        tmp.Add(pp);
                                }
                                else
                                {
                                    if (p[0] == 255)
                                        tmp.Add(pp);
                                }
                            }
                            p += BPP;
                        }
                        p += offset;
                    }
                    edgeImage.UnlockBits(data);
                }
                //PointF[] ps = tmp.ToArray();
                //if (ps != null && ps.Length > 1)
                if (tmp.Count > 0)
                {
                    GetLine(tmp, out p1, out p2);
                    return true;
                }
            }
            return false;
        }
        /// <summary>
        /// 有方向的线检测（返回直线的两点）
        /// </summary>
        /// <param name="rgn"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="needGrayExt"></param>
        /// <param name="angle"></param>
        /// <param name="step">步长必须大于0</param>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="linePs"></param>
        /// <returns></returns>
        public bool GetLineByDirection(Region rgn, int blkORwht, byte threshold, bool needGrayExt, double angle, byte step, out PointF p1, out PointF p2, out List<PointF> linePs)
        {
            p1 = PointF.Empty;
            p2 = PointF.Empty;
            linePs = null;
            if (step < 1)
                return false;
            if (bmp != null)
            {
                List<PointF> tmp = new List<PointF>();
                Graphics g = Graphics.FromImage((Image)bmp.Clone());
                RectangleF rect = rgn.GetBounds(g);
                Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(rect), new Rectangle(0, 0, width, height));
                int X = rec.X;
                int Y = rec.Y;
                int Width = rec.Width + rec.X;
                int Height = rec.Height + rec.Y;

                Bitmap b = Thresholding(bmp, rec, threshold, needGrayExt);
                b = GetBmpByRoberts((Bitmap)b.Clone());
                System.Collections.ArrayList al = GetPointListsByDir((Bitmap)b.Clone(), threshold, rgn, angle, step);
                if (al != null)
                {
                    linePs = new List<PointF>();

                    if (al.Count == 1)
                    {
                        List<PointF> s1 = (List<PointF>)al[0];
                        foreach (PointF p in s1)
                        {
                            tmp.Add(p);
                            linePs.Add(p);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < al.Count; i++)
                        {
                            List<PointF> ps = (List<PointF>)al[i];
                            if (ps != null && ps.Count > 0)
                            {
                                tmp.Add(ps[0]);//取第一个扫描到的点！
                                for (int j = 0; j < ps.Count; j++)
                                {
                                    linePs.Add(ps[j]);
                                }
                            }
                        }
                    }
                }
                /*string s="";
                foreach (Point p in tmp)
                {
                    s += p.ToString() + "\n";
                }
                //MessageBox.Show(s);*/
                b.Dispose();
                //PointF[] ps1 = tmp.ToArray();
                //if (ps1 != null && ps1.Length > 1)
                if (tmp.Count > 0)
                {
                    GetLine(tmp, out p1, out p2);
                    return true;
                }
            }
            return false;
        }
        /// <summary>
        /// 获得点集序列数组ArrayList(List{PointF}>,List{PointF},List{PointF},...)，其中每一个List{PointF}是该方向上的一条扫描线所经过的关键点集，直线方向在[135，180]或[315，360]的范围内 有问题...
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="blkORwht">取白点1还是黑点0</param>
        /// <param name="rgn">区域</param>
        /// <param name="angle">方向角度，单位：弧度</param>
        /// <param name="step">步长，必须大于0</param>
        private System.Collections.ArrayList GetPointListsByDir(Bitmap bmp, int blkORwht, Region rgn, double angle, byte step)
        {
            if (step < 1)
                return null;
            System.Collections.ArrayList LinshiList = null;
            if (bmp != null)
            {
                LinshiList = new System.Collections.ArrayList();

                bool order1 = true;//角度标识
                bool order2 = true;//角度标识
                bool order3 = true;//角度标识
                bool order4 = true;//角度标识
                bool order5 = true;//角度标识
                bool order6 = true;//角度标识
                angle = angle % (2 * Math.PI);
                if (angle < 0)
                    angle += 2 * Math.PI;
                if (angle > Math.PI && angle <= 5 * Math.PI / 4)
                {
                    order1 = false;//在(180,225]
                }
                else if (angle > 5 * Math.PI / 4 && angle < 3 * Math.PI / 2)
                {
                    order2 = false;//在(225,270)
                }
                else if (angle > 3 * Math.PI / 2 && angle <= 7 * Math.PI / 4)
                {
                    order3 = false;//在(270,315]
                }
                else if (angle > 7 * Math.PI / 4 && angle < Math.PI * 2)
                {
                    order4 = false;//在(315,360)
                }
                else if (angle == 3 * Math.PI / 2)
                {
                    order5 = false;//为270
                }
                else if (angle == Math.PI)
                {
                    order6 = false;//为180
                }

                if (angle > Math.PI)
                    angle = angle - Math.PI;
                Graphics g = Graphics.FromImage((Image)bmp.Clone());
                RectangleF rect = rgn.GetBounds(g);
                float offsetX = rect.X;
                float offsetY = rect.Y;
                rgn.Transform(new Matrix(1, 0, 0, 1, -offsetX, -offsetY));//.Translate(offsetX, offsetY);
                //区域的左上角
                int TopX = 0;
                int TopY = 0;
                //区域的高度
                int heigh = bmp.Height;
                //区域的宽度
                int width = bmp.Width;

                BitmapData srcData = bmp.LockBits(new Rectangle(TopX, TopY, width, heigh), ImageLockMode.ReadOnly, PixelFormat.Format32bppRgb);
                unsafe
                {
                    byte* p = (byte*)srcData.Scan0;
                    int stride = srcData.Stride;
                    if (angle > 0 && angle <= Math.PI / 4)
                    {
                        p += stride * heigh;

                        for (int j = TopY + heigh - 1; j > TopY + step - 1; j -= step)
                        {
                            List<PointF> EdgePoints = new List<PointF>();
                            for (int i = TopX; i < TopX + width; i++)
                            {
                                PointF temp = new PointF(i, j);

                                if (rgn.IsVisible(temp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p[0] == 0)
                                            EdgePoints.Add(new PointF(i + offsetX, j + offsetY));
                                    }
                                    else
                                    {
                                        if (p[0] == 255)
                                            EdgePoints.Add(new PointF(i + offsetX, j + offsetY));
                                    }
                                }
                                p += BPP;
                            }
                            if (EdgePoints.Count > 0)
                            {
                                if (!order1)
                                    EdgePoints = GetReverseList(EdgePoints);
                                LinshiList.Add(EdgePoints);
                            }
                            p -= stride * step + BPP * width;
                        }
                    }
                    else if (angle > Math.PI / 4 && angle < Math.PI / 2)
                    {
                        p += stride * heigh;

                        for (int i = TopX; i < TopX + width - step + 1; i += step)
                        {
                            List<PointF> EdgePoints = new List<PointF>();
                            for (int j = TopY + heigh - 1; j > TopY; j--)
                            {
                                PointF temp = new PointF(i, j);
                                if (rgn.IsVisible(temp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p[0] == 0)
                                            EdgePoints.Add(new PointF(i + offsetX, j + offsetY));
                                    }
                                    else
                                    {
                                        if (p[0] == 255)
                                            EdgePoints.Add(new PointF(i + offsetX, j + offsetY));
                                    }
                                }
                                p -= stride;
                            }
                            if (EdgePoints.Count > 0)
                            {
                                if (!order2)
                                    EdgePoints = GetReverseList(EdgePoints);
                                LinshiList.Add(EdgePoints);
                            }
                            p += stride * heigh + BPP * step;
                        }
                    }
                    else if (angle > Math.PI / 2 && angle <= 3 * Math.PI / 4)
                    {
                        p += BPP * width + stride * heigh;

                        for (int i = TopX + width - 1; i > TopX + step - 1; i -= step)
                        {
                            List<PointF> EdgePoints = new List<PointF>();
                            for (int j = TopY + heigh - 1; j > TopY; j--)
                            {
                                PointF temp = new PointF(i, j);
                                if (rgn.IsVisible(temp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p[0] == 0)
                                            EdgePoints.Add(new PointF(i + offsetX, j + offsetY));
                                    }
                                    else
                                    {
                                        if (p[0] == 255)
                                            EdgePoints.Add(new PointF(i + offsetX, j + offsetY));
                                    }
                                }
                                p -= stride;
                            }
                            if (EdgePoints.Count > 0)
                            {
                                if (!order3)
                                    EdgePoints = GetReverseList(EdgePoints);
                                LinshiList.Add(EdgePoints);
                            }
                            p += stride * heigh - BPP * step;
                        }
                    }
                    else if (angle > 3 * Math.PI / 4 && angle < Math.PI)//非常奇怪，在这个范围内取到的点集都有问题！
                    {
                        p += BPP * width + stride * heigh;

                        for (int j = TopY + heigh - 1; j > TopY + step - 1; j -= step)
                        {
                            List<PointF> EdgePoints = new List<PointF>();
                            for (int i = TopX + width - 1; i > TopX; i--)
                            {
                                PointF temp = new PointF(i, j);
                                if (rgn.IsVisible(temp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p[0] == 0)
                                            EdgePoints.Add(new PointF(i + offsetX, j + offsetY));
                                    }
                                    else
                                    {
                                        if (p[0] == 255)
                                            EdgePoints.Add(new PointF(i + offsetX, j + offsetY));
                                    }
                                }
                                p -= BPP;
                            }
                            if (EdgePoints.Count > 0)
                            {
                                if (!order4)
                                    EdgePoints = GetReverseList(EdgePoints);
                                LinshiList.Add(EdgePoints);
                            }
                            p += BPP * width - stride * step;
                        }
                    }
                    else if (angle == Math.PI / 2)
                    {
                        p += stride * heigh;

                        for (int i = TopX; i < TopX + width - step + 1; i += step)
                        {
                            List<PointF> EdgePoints = new List<PointF>();
                            for (int j = TopY + heigh - 1; j > TopY; j--)
                            {
                                PointF temp = new PointF(i, j);
                                if (rgn.IsVisible(temp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p[0] == 0)
                                            EdgePoints.Add(new PointF(i + offsetX, j + offsetY));
                                    }
                                    else
                                    {
                                        if (p[0] == 255)
                                            EdgePoints.Add(new PointF(i + offsetX, j + offsetY));
                                    }
                                }
                                p -= stride;
                            }
                            if (EdgePoints.Count > 0)
                            {
                                if (!order5)
                                    EdgePoints = GetReverseList(EdgePoints);
                                LinshiList.Add(EdgePoints);
                            }
                            p += stride * heigh + BPP * step;
                        }
                    }
                    else if (angle == 0)
                    {
                        for (int j = TopY; j < TopY + heigh - step + 1; j += step)
                        {
                            List<PointF> EdgePoints = new List<PointF>();
                            for (int i = TopX; i < TopX + width; i++)
                            {
                                PointF temp = new PointF(i, j);
                                if (rgn.IsVisible(temp))
                                {
                                    if (blkORwht == 0)
                                    {
                                        if (p[0] == 0)
                                            EdgePoints.Add(new PointF(i + offsetX, j + offsetY));
                                    }
                                    else
                                    {
                                        if (p[0] == 255)
                                            EdgePoints.Add(new PointF(i + offsetX, j + offsetY));
                                    }
                                }
                                p += BPP;
                            }
                            if (EdgePoints.Count > 0)
                            {
                                if (!order6)
                                    EdgePoints = GetReverseList(EdgePoints);
                                LinshiList.Add(EdgePoints);
                            }
                            p += stride * step - BPP * width;
                        }
                    }
                }
                bmp.UnlockBits(srcData);
            }
            return LinshiList;
        }
        /// <summary>
        /// 获得点集序列数组ArrayList(List{PointF},List{PointF},List{PointF},...)，其中每一个List{PointF}是该方向上的一条扫描线所经过的关键点集，有问题...
        /// </summary>
        /// <param name="bmp">位图</param>
        /// <param name="blkORwht">取白点1还是黑点0</param>
        /// <param name="threshold">阈值</param>
        /// <param name="rgn">区域</param>
        /// <param name="angle">方向角度，单位：弧度</param>
        public static System.Collections.ArrayList GetPointListsByDirection(Bitmap bmp, int blkORwht, byte threshold, Region rgn, double angle)
        {
            System.Collections.ArrayList LinshiList = null;
            if (bmp != null)
            {
                //距离
                double d = 0;
                LinshiList = new System.Collections.ArrayList();
                //角度标识，是否大于180度。
                bool order = true;
                angle = angle % (2 * Math.PI);
                if (angle >= Math.PI)
                {
                    angle = angle - Math.PI;
                    order = false;
                }
                Graphics g = Graphics.FromImage((Image)bmp.Clone());
                RectangleF rect = rgn.GetBounds(g);
                //区域的高度
                int heigh = (int)rect.Height;
                //区域的宽度
                int width = (int)rect.Width;
                //区域的左上角
                int TopX = (int)rect.X;
                int TopY = (int)rect.Y;
                double h = 0;

                double YuZhi = 0;
                double b = 0;
                double k = 0;
                int h1 = 0, w1 = 0;
                //求域值
                if (angle > Math.PI / 2)
                {
                    h = Math.Cos(angle - Math.PI / 2);
                    YuZhi = h / 2;

                    k = Math.Tan(angle);
                    b = -TopY - k * (TopX + width);
                    h1 = (int)(TopY + (k * TopX + b));
                    w1 = (int)((-b - TopY - heigh) / k) - TopX;
                }
                else if (angle == Math.PI / 2)
                {

                }
                else if (angle == 0)
                {

                }
                else
                {
                    h = Math.Cos(angle);
                    YuZhi = h / 2;

                    k = Math.Tan(angle);
                    b = -TopY - k * TopX;
                    h1 = (int)(TopY + (k * (TopX + width) + b));
                    w1 = TopX - (int)((-b - TopY - heigh) / k);
                }

                //直线方程Y = k*X+b; 

                //求高度
                int Heigh = heigh + h1;
                //求宽度
                int Width = width + w1;
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                BitmapData srcData = bmp.LockBits(new Rectangle(TopX, TopY, width, heigh), ImageLockMode.ReadOnly, bmp.PixelFormat);
                unsafe
                {
                    byte* pIn = (byte*)srcData.Scan0.ToPointer();
                    int stride = srcData.Stride;
                    byte* t, temp9, temp8;
                    byte value1;
                    byte max = 0, value2, value3;
                    temp9 = pIn;
                    //循环求边缘点 
                    if (angle > 0 && angle <= Math.PI / 4)
                    {
                        for (int i = TopY; i < TopY + Heigh; i = i + 2)
                        {
                            b = -i - k * TopX;   //直线方程为：y = kx+b；    
                            t = temp9 + stride + BPP;
                            List<PointF> EdgePoints = new List<PointF>();
                            for (int q = TopY + 1; q < TopY + heigh - 1; q++)
                            {
                                for (int p = TopX + 1; p < TopX + width - 1; p++)
                                {
                                    PointF temp = new PointF(p, q);

                                    if (rgn.IsVisible(temp))
                                    {
                                        d = GetPLDistance(k, b, new Point(p, -q));
                                        if (d < YuZhi)
                                        {
                                            #region//符合边缘条件的点，就加入到临时数组中去。
                                            pIn = t;//(byte)((19661 * R + 38666 * G + 7209 * B) >> 16)
                                            value1 = (byte)((19661 * pIn[2] + 38666 * pIn[1] + 7209 * pIn[0]) >> 16);//转化成灰度
                                            //正上
                                            temp8 = pIn - stride;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            max = value3;

                                            //右上
                                            temp8 = pIn - stride + BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //右侧
                                            temp8 = pIn + BPP;//(byte)((19661 * R + 38666 * G + 7209 * B) >> 16)
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //右下
                                            temp8 = pIn + stride + BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //正下
                                            temp8 = pIn + stride;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //左下
                                            temp8 = pIn + stride - BPP;//(byte)((19661 * R + 38666 * G + 7209 * B) >> 16)
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //左侧
                                            temp8 = pIn - BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //左上
                                            temp8 = pIn + stride - BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            if (blkORwht == 0)
                                            {
                                                if (max < threshold)
                                                {
                                                    EdgePoints.Add(temp);
                                                }
                                            }
                                            else
                                            {
                                                if (max > threshold)
                                                {
                                                    EdgePoints.Add(temp);
                                                }
                                            }
                                            #endregion
                                        }
                                    }
                                    t = t + BPP;
                                }

                                t += stride - width * BPP + 2 * BPP;
                            }
                            if (EdgePoints.Count > 0)
                            {
                                if (!order)
                                {
                                    EdgePoints = GetReverseList(EdgePoints);
                                }
                                LinshiList.Add(EdgePoints);
                            }
                        }
                    }
                    else if (angle > Math.PI / 4 && angle < Math.PI / 2)
                    {
                        for (int i = TopX; i < TopX + Width; i = i + 2)
                        {
                            b = -TopY - k * i;//直线方程为：y = kx+b； 
                            int N9 = 0;
                            t = temp9 + stride + BPP;
                            List<PointF> EdgePoints = new List<PointF>();
                            for (int p = TopX + 1; p < TopX + width - 1; p++)
                            {
                                for (int q = TopY + 1; q < TopY + heigh - 1; q++)
                                {
                                    PointF temp = new PointF(p, q);
                                    if (rgn.IsVisible(temp))
                                    {
                                        d = GetPLDistance(k, b, new Point(p, -q));
                                        if (d < YuZhi)
                                        {
                                            #region//符合边缘条件的点，就加入到临时数组中去。
                                            pIn = t;//(byte)((19661 * R + 38666 * G + 7209 * B) >> 16)
                                            value1 = (byte)((19661 * pIn[2] + 38666 * pIn[1] + 7209 * pIn[0]) >> 16);//转化成灰度
                                            //正上
                                            temp8 = pIn - stride;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            max = value3;

                                            //右上
                                            temp8 = pIn - stride + BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //右侧
                                            temp8 = pIn + BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //右下
                                            temp8 = pIn + stride + BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //正下
                                            temp8 = pIn + stride;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //左下
                                            temp8 = pIn + stride - BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //左侧
                                            temp8 = pIn - BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //左上
                                            temp8 = pIn + stride - BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }

                                            if (blkORwht == 0)
                                            {
                                                if (max < threshold)
                                                {
                                                    EdgePoints.Add(temp);
                                                }
                                            }
                                            else
                                            {
                                                if (max > threshold)
                                                {
                                                    EdgePoints.Add(temp);
                                                }
                                            }
                                            #endregion
                                        }
                                    }
                                    t = t + stride;
                                }
                                N9++;
                                t = temp9 + stride + BPP + N9 * BPP;
                            }
                            if (EdgePoints.Count > 0)
                            {
                                if (!order)
                                {
                                    EdgePoints = GetReverseList(EdgePoints);
                                }
                                LinshiList.Add(EdgePoints);
                            }
                        }
                    }
                    else if (angle > Math.PI / 2 && angle <= 3 * Math.PI / 4)
                    {//(0,45)
                        for (int i = TopX + width; i > TopX - w1; i = i - 2)
                        {
                            b = -TopY - k * i;//直线方程为：y = kx+b；     

                            int N9 = 1;
                            t = temp9 + BPP + stride;
                            List<PointF> EdgePoints = new List<PointF>();
                            for (int p = TopX + 1; p < TopX + width - 1; p++)
                            {
                                for (int q = TopY + 1; q < TopY + heigh - 1; q++)
                                {
                                    PointF temp = new PointF(p, q);
                                    if (rgn.IsVisible(temp))
                                    {
                                        d = GetPLDistance(k, b, new Point(p, -q));
                                        if (d < YuZhi)
                                        {
                                            #region//符合边缘条件的点，就加入到临时数组中去。
                                            pIn = t;//(byte)((19661 * R + 38666 * G + 7209 * B) >> 16)
                                            value1 = (byte)((19661 * pIn[2] + 38666 * pIn[1] + 7209 * pIn[0]) >> 16);//转化成灰度
                                            //正上
                                            temp8 = pIn - stride;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            max = value3;

                                            //右上
                                            temp8 = pIn - stride + BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //右侧
                                            temp8 = pIn + BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //右下
                                            temp8 = pIn + stride + BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //正下
                                            temp8 = pIn + stride;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //左下
                                            temp8 = pIn + stride - BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //左侧
                                            temp8 = pIn - BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //左上
                                            temp8 = pIn + stride - BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            if (blkORwht == 0)
                                            {
                                                if (max < threshold)
                                                {
                                                    EdgePoints.Add(temp);
                                                }
                                            }
                                            else
                                            {
                                                if (max > threshold)
                                                {
                                                    EdgePoints.Add(temp);
                                                }
                                            }
                                            #endregion
                                        }
                                    }
                                    t = t + stride;
                                }
                                N9++;
                                t = temp9 + stride + N9 * BPP;
                            }
                            if (EdgePoints.Count > 0)
                            {
                                if (!order)
                                {
                                    EdgePoints = GetReverseList(EdgePoints);
                                }
                                LinshiList.Add(EdgePoints);
                            }
                        }
                    }
                    else if (angle > 3 * Math.PI / 4 && angle < Math.PI)
                    {//(45,90)
                        for (int i = TopY; i < TopY + Heigh; i = i + 2)
                        {
                            b = -i - k * (TopX + width);   //直线方程为：y = kx+b；    
                            t = temp9 + stride + BPP;
                            List<PointF> EdgePoints = new List<PointF>();
                            for (int q = TopY + 1; q < TopY + heigh - 1; q++)
                            {
                                //for (int q = Math.Max(0,(int)((h1-i)/k)); q < Math.Min( width,(int)(heigh-i)/k); q++)
                                for (int p = TopX + 1; p < TopX + width - 1; p++)
                                {
                                    PointF temp = new PointF(p, q);
                                    if (rgn.IsVisible(temp))
                                    {
                                        d = GetPLDistance(k, b, new Point(p, -q));
                                        if (d < YuZhi)
                                        {
                                            #region//符合边缘条件的点，就加入到临时数组中去。
                                            pIn = t;//(byte)((19661 * R + 38666 * G + 7209 * B) >> 16)
                                            value1 = (byte)((19661 * pIn[2] + 38666 * pIn[1] + 7209 * pIn[0]) >> 16);//转化成灰度
                                            //正上
                                            temp8 = pIn - stride;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            max = value3;

                                            //右上
                                            temp8 = pIn - stride + BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //右侧
                                            temp8 = pIn + BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //右下
                                            temp8 = pIn + stride + BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //正下
                                            temp8 = pIn + stride;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //左下
                                            temp8 = pIn + stride - BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //左侧
                                            temp8 = pIn - BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }
                                            //左上
                                            temp8 = pIn + stride - BPP;
                                            value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                            value3 = (byte)Math.Abs(value1 - value2);
                                            if (max < value3)
                                            {
                                                max = value3;
                                            }

                                            if (blkORwht == 0)
                                            {
                                                if (max < threshold)
                                                {
                                                    EdgePoints.Add(temp);
                                                }
                                            }
                                            else
                                            {
                                                if (max > threshold)
                                                {
                                                    EdgePoints.Add(temp);
                                                }
                                            }
                                            #endregion
                                        }
                                    }
                                    t = t + BPP;
                                }
                                t += stride - width * BPP + 2 * BPP;
                            }
                            if (EdgePoints.Count > 0)
                            {
                                if (!order)
                                {
                                    EdgePoints = GetReverseList(EdgePoints);
                                }
                                LinshiList.Add(EdgePoints);
                            }
                        }
                    }
                    else if (angle == 0)
                    {
                        t = temp9 + stride + BPP;
                        for (int q = TopY + 1; q < TopY + heigh - 1; q++)
                        {
                            List<PointF> EdgePoints = new List<PointF>();
                            for (int p = TopX + 1; p < TopX + width - 1; p++)
                            {
                                PointF temp = new PointF(p, q);
                                if (rgn.IsVisible(temp))
                                {
                                    #region//符合边缘条件的点，就加入到临时数组中去。
                                    pIn = t;//(byte)((19661 * R + 38666 * G + 7209 * B) >> 16)
                                    value1 = (byte)((19661 * pIn[2] + 38666 * pIn[1] + 7209 * pIn[0]) >> 16);//转化成灰度
                                    //正上
                                    temp8 = pIn - stride;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    max = value3;

                                    //右上
                                    temp8 = pIn - stride + BPP;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    if (max < value3)
                                    {
                                        max = value3;
                                    }
                                    //右侧
                                    temp8 = pIn + BPP;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    if (max < value3)
                                    {
                                        max = value3;
                                    }
                                    //右下
                                    temp8 = pIn + stride + BPP;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    if (max < value3)
                                    {
                                        max = value3;
                                    }
                                    //正下
                                    temp8 = pIn + stride;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    if (max < value3)
                                    {
                                        max = value3;
                                    }
                                    //左下
                                    temp8 = pIn + stride - BPP;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    if (max < value3)
                                    {
                                        max = value3;
                                    }
                                    //左侧
                                    temp8 = pIn - BPP;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    if (max < value3)
                                    {
                                        max = value3;
                                    }
                                    //左上
                                    temp8 = pIn + stride - BPP;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    if (max < value3)
                                    {
                                        max = value3;
                                    }
                                    if (blkORwht == 0)
                                    {
                                        if (max < threshold)
                                        {
                                            EdgePoints.Add(temp);
                                        }
                                    }
                                    else
                                    {
                                        if (max > threshold)
                                        {
                                            EdgePoints.Add(temp);
                                        }
                                    }
                                    #endregion
                                }
                                t = t + BPP;
                                if (EdgePoints.Count > 0)
                                {
                                    if (!order)
                                    {
                                        EdgePoints = GetReverseList(EdgePoints);
                                    }
                                    LinshiList.Add(EdgePoints);
                                }
                            }
                            t += stride - width * BPP + 2 * BPP;
                        }
                    }
                    else if (angle == Math.PI / 2)
                    {
                        int N9 = 1;
                        t = temp9 + stride + BPP;
                        for (int p = TopX + 1; p < TopX + width - 1; p++)
                        {
                            List<PointF> EdgePoints = new List<PointF>();
                            for (int q = TopY + 1; q < TopY + heigh - 1; q++)
                            {
                                PointF temp = new PointF(p, q);
                                if (rgn.IsVisible(temp))
                                {
                                    #region//符合边缘条件的点，就加入到临时数组中去。
                                    pIn = t;//(byte)((19661 * R + 38666 * G + 7209 * B) >> 16)
                                    value1 = (byte)((19661 * pIn[2] + 38666 * pIn[1] + 7209 * pIn[0]) >> 16);//转化成灰度
                                    //正上
                                    temp8 = pIn - stride;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    max = value3;

                                    //右上
                                    temp8 = pIn - stride + BPP;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    if (max < value3)
                                    {
                                        max = value3;
                                    }
                                    //右侧
                                    temp8 = pIn + BPP;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    if (max < value3)
                                    {
                                        max = value3;
                                    }
                                    //右下
                                    temp8 = pIn + stride + BPP;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    if (max < value3)
                                    {
                                        max = value3;
                                    }
                                    //正下
                                    temp8 = pIn + stride;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    if (max < value3)
                                    {
                                        max = value3;
                                    }
                                    //左下
                                    temp8 = pIn + stride - BPP;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    if (max < value3)
                                    {
                                        max = value3;
                                    }
                                    //左侧
                                    temp8 = pIn - BPP;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    if (max < value3)
                                    {
                                        max = value3;
                                    }
                                    //左上
                                    temp8 = pIn + stride - BPP;
                                    value2 = (byte)((19661 * temp8[2] + 38666 * temp8[1] + 7209 * temp8[0]) >> 16);//转化成灰度
                                    value3 = (byte)Math.Abs(value1 - value2);
                                    if (max < value3)
                                    {
                                        max = value3;
                                    }
                                    if (blkORwht == 0)
                                    {
                                        if (max < threshold)
                                        {
                                            EdgePoints.Add(temp);
                                        }
                                    }
                                    else
                                    {
                                        if (max > threshold)
                                        {
                                            EdgePoints.Add(temp);
                                        }
                                    }
                                    #endregion
                                }
                                t = t + stride;
                            }

                            if (EdgePoints.Count > 0)
                            {
                                if (!order)
                                {
                                    EdgePoints = GetReverseList(EdgePoints);
                                }
                                LinshiList.Add(EdgePoints);
                            }

                            N9++;
                            t = temp9 + stride + N9 * BPP;
                        }

                    }

                }
                bmp.UnlockBits(srcData);
            }
            return LinshiList;
        }
        /// <summary>
        /// 获得点p到直线的距离。直线为：k*x-y+b= 0
        /// </summary>
        /// <param name="k">直线斜率</param>
        /// <param name="b">参数b</param>
        /// <param name="p">点</param>
        /// <returns></returns>
        private static double GetPLDistance(double k, double b, Point p)
        {
            return Math.Abs(k * p.X - p.Y + b) / Math.Sqrt(1 + k * k);
        }
        /// <summary>
        /// 获取逆序点数组
        /// </summary>
        /// <param name="al">临时点数组</param>
        /// <returns>返回获取的点数组</returns>
        private static List<PointF> GetReverseList(List<PointF> al)
        {
            List<PointF> points = new List<PointF>();
            for (int i = 0; i < al.Count; i++)
            {
                points.Add(al[al.Count - 1 - i]);
            }
            return points;
        }
        /// <summary>
        /// 改进的多点取线
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="startP"></param>
        /// <param name="endP"></param>
        public static void GetLine1(List<PointF> ps, out PointF startP, out PointF endP)
        {
            double y = 0.0, x = 0.0, yx = 0.0, xx = 0.0;
            ////y=b+kx
            for (int j = 0; j < ps.Count; j++)
            {
                x += ps[j].X;
                y += ps[j].Y;
                yx += ps[j].X * ps[j].Y;
                xx += ps[j].X * ps[j].X;
            }
            double k = (y * x - yx * ps.Count) / (x * x - xx * ps.Count);
            double b = (y - k * x) / ps.Count;
            PointF[] p4 = RectangleOnMorePoint1(ps);
            PointF p1 = new PointF();
            p1.X = 0;
            p1.Y = (float)b;
            PointF p2 = new PointF();
            p2.X = 10;
            p2.Y = (float)(b + k * 10);
            PointF[] pp1 = new PointF[2];
            pp1[0] = Pvint(p4[0], p1, p2);
            pp1[1] = Pvint(p4[2], p1, p2);
            PointF[] pp2 = new PointF[2];
            pp2[0] = Pvint(p4[1], p1, p2);
            pp2[1] = Pvint(p4[3], p1, p2);
            double t1 = Tpoint(pp1[0], pp1[1]);
            double t2 = Tpoint(pp2[0], pp2[1]);
            if (t1 > t2)
            {
                startP = pp1[0];
                endP = pp1[1];
            }
            else
            {
                startP = pp2[0];
                endP = pp2[1];
            }
        }
        /// <summary>
        /// 改进的多点取线
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="startP"></param>
        /// <param name="endP"></param>
        public static void GetLine1(List<Point> ps, out PointF startP, out PointF endP)
        {
            double y = 0.0, x = 0.0, yx = 0.0, xx = 0.0;
            ////y=b+kx
            for (int j = 0; j < ps.Count; j++)
            {
                x += ps[j].X;
                y += ps[j].Y;
                yx += ps[j].X * ps[j].Y;
                xx += ps[j].X * ps[j].X;
            }
            double k = (y * x - yx * ps.Count) / (x * x - xx * ps.Count);
            double b = (y - k * x) / ps.Count;
            Point[] p4 = RectangleOnMorePoint1(ps);
            PointF p1 = new PointF();
            p1.X = 0;
            p1.Y = (float)b;
            PointF p2 = new PointF();
            p2.X = 10;
            p2.Y = (float)(b + k * 10);
            PointF[] pp1 = new PointF[2];
            pp1[0] = Pvint(p4[0], p1, p2);
            pp1[1] = Pvint(p4[2], p1, p2);
            PointF[] pp2 = new PointF[2];
            pp2[0] = Pvint(p4[1], p1, p2);
            pp2[1] = Pvint(p4[3], p1, p2);
            double t1 = Tpoint(pp1[0], pp1[1]);
            double t2 = Tpoint(pp2[0], pp2[1]);
            if (t1 > t2)
            {
                startP = pp1[0];
                endP = pp1[1];
            }
            else
            {
                startP = pp2[0];
                endP = pp2[1];
            }
        }
        /// <summary>
        ///求两点的矩离
        /// </summary>
        /// <param name="a">第一点</param>
        /// <param name="b">第二点</param>
        /// <returns>返回两点的矩离</returns>
        public static double Tpoint(PointF a, PointF b)
        {
            double y = 0.0;
            y = (a.X - b.X) * (a.X - b.X) + (a.Y - b.Y) * (a.Y - b.Y);
            y = Math.Sqrt(y);
            return y;
        }
        /// <summary>
        /// 点集的外包矩阵
        /// </summary>
        /// <param name="pfs">点集</param>
        /// <returns>正规矩阵的四个顶点,顺时钟</returns>
        public static Point[] RectangleOnMorePoint1(List<Point> pfs)
        {
            Point[] re = new Point[4];
            if (pfs.Count > 0)
            {
                int maxx = 0, minx = 0, maxy = 0, miny = 0;
                minx = pfs[0].X;
                miny = pfs[0].Y;
                for (int i = 1; i < pfs.Count; i++)
                {
                    if (maxx < pfs[i].X)
                        maxx = pfs[i].X;
                    if (maxy < pfs[i].Y)
                        maxy = pfs[i].Y;
                    if (minx > pfs[i].X)
                        minx = pfs[i].X;
                    if (miny > pfs[i].Y)
                        miny = pfs[i].Y;
                }
                re[0].X = minx;
                re[0].Y = miny;
                re[1].X = maxx;
                re[1].Y = miny;
                re[2].X = maxx;
                re[2].Y = maxy;
                re[3].X = minx;
                re[3].Y = maxy;
            }
            return re;
        }
        /// <summary>
        /// 点集的外包矩阵
        /// </summary>
        /// <param name="pfs">点集</param>
        /// <returns>正规矩阵的四个顶点,顺时钟</returns>
        public static PointF[] RectangleOnMorePoint1(List<PointF> pfs)
        {
            PointF[] re = new PointF[4];
            if (pfs.Count > 0)
            {
                float maxx = 0, minx = 0, maxy = 0, miny = 0;
                minx = pfs[0].X;
                miny = pfs[0].Y;
                for (int i = 1; i < pfs.Count; i++)
                {
                    if (maxx < pfs[i].X)
                        maxx = pfs[i].X;
                    if (maxy < pfs[i].Y)
                        maxy = pfs[i].Y;
                    if (minx > pfs[i].X)
                        minx = pfs[i].X;
                    if (miny > pfs[i].Y)
                        miny = pfs[i].Y;
                }
                re[0].X = minx;
                re[0].Y = miny;
                re[1].X = maxx;
                re[1].Y = miny;
                re[2].X = maxx;
                re[2].Y = maxy;
                re[3].X = minx;
                re[3].Y = maxy;
            }
            return re;
        }
        /// <summary>
        /// 过直线外一点的在直线上的垂足
        /// </summary>
        /// <param name="pt">直线外的一点</param>
        /// <param name="lt1">直线上的一点</param>
        /// <param name="lt2">直线上的另一点</param>
        /// <returns></returns>
        public static PointF Pvint(PointF pt, PointF lt1, PointF lt2)
        {

            PointF t = new PointF();
            double x = 0.0;
            double y = 0.0;
            if (lt1.X != lt2.X)
            {
                if (lt1.Y != lt2.Y)
                {


                    double k1 = 0.0;
                    double k2 = 0.0;
                    k1 = (lt2.Y - lt1.Y) / ((double)lt2.X - lt1.X);
                    k2 = -1 / k1;//垂直线的斜率

                    x = ((lt1.Y * lt2.X - lt2.Y * lt1.X) / (lt2.X - lt1.X) + k2 * pt.X - pt.Y) / (k2 - k1);
                    y = k2 * x - k2 * pt.X + pt.Y;
                }
                else
                {
                    x = pt.X;//与X平行
                    y = lt1.Y;

                }
            }
            else
            {
                x = lt1.X;//与Y平行
                y = pt.Y;
            }
            t.X = (float)x;
            t.Y = (float)y;
            return t;
        }
        /// <summary>
        /// 直线拟合算法，多点点分布定线
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="startP"></param>
        /// <param name="endP"></param>
        public static void GetLine(List<PointF> ps, out PointF startP, out PointF endP)
        {
            Algorithm NiheLine = new Algorithm();
            NiheLine.NiHeLine(ps, out startP, out endP);
            /*//屏蔽掉，因为GetBoundFromPs有问题...
            startP = PointF.Empty;
            endP = PointF.Empty;
            if (ps.Length < 2)
            {
                throw new Exception("检测到的点数少于两点，无法取线！");
            }
            if (ps.Length == 2)
            {
                startP = ps[0];
                endP = ps[1];
                return;
            }
            PointF[] pp = new PointF[2];
            RectangleF bound = RectangleF.Empty;
            try
            {
                bound = GetBoundFromPs(ps, out pp);
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.Message);
            }
            if (bound == RectangleF.Empty)
                return;
            
            if (startP == endP)
            {
                startP = ps[0];
                endP = ps[ps.Length - 1];
            }
            else
            {
                startP = pp[0];
                endP = pp[1];
            }
            Vector.Vector2 v = new Vector.Vector2(new Vector.Point2(startP.X, startP.Y), new Vector.Point2(endP.X, endP.Y));
            Vector.Point2 startP1, endP1;
            double lenS = 0;
            double lenE = 0;
            if (Distance(startP, new PointF(bound.X, bound.Y)) < Distance(endP, new PointF(bound.X, bound.Y)))
            {
                lenS = Distance(startP, new PointF(bound.X, bound.Y));
                lenE = Distance(endP, new PointF(bound.X + bound.Width, bound.Y + bound.Height)) + Vector.Vector2.GetLength(v);
                startP1 = Vector.Vector2.GetPointOfLengthFormStartP(v, -lenS);
                endP1 = Vector.Vector2.GetPointOfLengthFormStartP(v, lenE);
            }
            else
            {
                lenS = Distance(startP, new PointF(bound.X + bound.Width, bound.Y + bound.Height)) + Vector.Vector2.GetLength(v);
                lenE = Distance(endP, new PointF(bound.X, bound.Y));
                startP1 = Vector.Vector2.GetPointOfLengthFormStartP(v, lenS);
                endP1 = Vector.Vector2.GetPointOfLengthFormStartP(v, -lenE);
            }

            startP = new PointF(startP1.X, startP1.Y);
            endP = new PointF(endP1.X, endP1.Y);
            */
        }
        /// <summary>
        /// 两点之间的距离
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static double Distance2D(Point p1, Point p2)
        {
            return Math.Sqrt((p1.X - p2.X) * (p1.X - p2.X) + (p1.Y - p2.Y) * (p1.Y - p2.Y));
        }
        /// <summary>
        /// 两点之间的距离
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static double Distance2D(PointF p1, PointF p2)
        {
            return Math.Sqrt((p1.X - p2.X) * (p1.X - p2.X) + (p1.Y - p2.Y) * (p1.Y - p2.Y));
        }
        /// <summary>
        /// 两点之间的距离
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static double Distance3D(Point3 p1, Point3 p2)
        {
            return Math.Sqrt((p1.X - p2.X) * (p1.X - p2.X) + (p1.Y - p2.Y) * (p1.Y - p2.Y) + (p1.Z - p2.Z) * (p1.Z - p2.Z));
        }
        /// <summary>
        /// 求矩形中心
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <param name="p4"></param>
        /// <param name="center"></param>
        public static void GetRectangleCenter(PointF p1, PointF p2, PointF p3, PointF p4, out PointF center)
        {
            Point3 pp1 = new Point3(p1.X, p1.Y, 0);
            Point3 pp2 = new Point3(p2.X, p2.Y, 0);
            Point3 pp3 = new Point3(p3.X, p3.Y, 0);
            Point3 pp4 = new Point3(p4.X, p4.Y, 0);
            Vector3 v1 = new Vector3(pp1, pp3, false);
            Vector3 v2 = new Vector3(pp2, pp4, false);
            Point3 c = Vector3.GetCross(v1, v2);
            if (c == null)
                throw new Exception("找不到矩形中点！");
            center = new PointF(c.X, c.Y);
        }
        /// <summary>
        /// 获取有序点集构造封闭曲线(凸多边形)的周长，注意：有序
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static double GetPsLength(Point[] ps)//封闭曲线(凸多边形)的长度
        {
            if (ps == null || ps.Length == 0)
            {
                throw new Exception("点集为空！");
            }
            double lg = 0.0;
            for (int i = 0; i < ps.Length - 1; i++)
            {
                lg += Distance2D(ps[i], ps[i + 1]);
            }
            lg += Distance2D(ps[ps.Length - 1], ps[0]);
            //lg = ps.Length;
            return lg;
        }
        /// <summary>
        /// 偶数为正，奇数为负
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        public static int IsNeg(int i)
        {
            int neg = 1;
            if (i % 2 == 0)
                neg = 1;
            else
                neg = -1;
            return neg;
        }

        /// <summary>
        /// 获取矩阵D在i行j列的子矩阵
        /// </summary>
        /// <param name="D"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        /// <returns></returns>
        public static double[,] Met(double[,] D, int i, int j)
        {
            int xR = D.GetLength(0);
            int yR = D.GetLength(1);
            int minR = xR < yR ? xR : yR;
            double[,] met;
            if (minR == 1)
            {
                met = new double[1, 1];
                met[0, 0] = D[0, 0];
            }
            else
            {
                met = new double[xR - 1, yR - 1];
                int i1cnt = 0, j1cnt = 0;
                for (int i1 = 0; i1 < xR; i1++)
                {
                    if (i1 == i)
                        continue;
                    j1cnt = 0;
                    for (int j1 = 0; j1 < yR; j1++)
                    {
                        if (j1 == j)
                            continue;
                        met[i1cnt, j1cnt] = D[i1, j1];
                        j1cnt++;
                    }
                    i1cnt++;
                }
            }
            return met;
        }

        /// <summary>
        /// 获取矩阵的行列式值
        /// </summary>
        /// <param name="D"></param>
        /// <returns></returns>
        public static double Det(double[,] D)
        {
            double det = 0;
            int xR = D.GetLength(0);
            int yR = D.GetLength(1);
            int minR = xR < yR ? xR : yR;
            if (minR == 1)
            {
                det = D[0, 0];
            }
            else
            {
                for (int j = 0; j < yR; j++)
                {
                    det += D[0, j] * IsNeg(0 + j) * Det(Met(D, 0, j));
                }
            }
            return det;
        }
        /// <summary>
        /// 凸多边形(封闭曲线)的面积
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static double GetPsMJ(Point[] ps)//凸多边形的面积(封闭曲线) 通过行列式计算，已经没问题。
        {
            if (ps == null || ps.Length == 0)
            {
                throw new Exception("点集为空！");
            }
            if (ps != null && ps.Length < 3)
            {
                throw new Exception("点集少于3点！");
            }
            double mj = 0;
            double[,] area = new double[3, 3];
            area[0, 0] = 1;
            area[0, 1] = 1;
            area[0, 2] = 1;
            //定点
            area[1, 0] = ps[0].X;
            area[2, 0] = ps[0].Y;
            for (int i = 1; i < ps.Length - 1; i++)
            {
                //area[1, 0] = ps[0].X;
                area[1, 1] = ps[i].X;
                area[1, 2] = ps[i + 1].X;
                //area[2, 0] = ps[0].Y;
                area[2, 1] = ps[i].Y;
                area[2, 2] = ps[i + 1].Y;
                mj += Det(area);
            }
            //要记得闭合多边形：
            //area[1, 0] = ps[0].X;
            area[1, 1] = ps[ps.Length - 1].X;
            area[1, 2] = ps[0].X;
            //area[2, 0] = ps[0].Y;
            area[2, 1] = ps[ps.Length - 1].Y;
            area[2, 2] = ps[0].Y;
            mj += Det(area);
            /*for (int i = 1; i < ps.Length - 1; i++)
            {
                mj += ps[i].X * ps[i + 1].Y - ps[i + 1].X * ps[i].Y;
            }*/
            mj = Math.Abs(0.5 * mj);
            return mj;
        }
        /// <summary>
        /// 已知平行四边形的三个顶点(a,b,c)，计算第四个顶点d的坐标. 注意：已知的三个顶点可以是有序的，b为中间点
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <returns></returns>
        public static Point Rect4thP(Point a, Point b, Point c)
        {
            Point d = new Point();
            /*
            if (Math.Abs(Dotmultiply(a, b, c)) < Common.Eps) // 说明c点是直角拐角处
            {
                d.X = a.X + b.X - c.X;
                d.Y = a.Y + b.Y - c.Y;
            }
            if (Math.Abs(Dotmultiply(a, c, b)) < Common.Eps) // 说明b点是直角拐角处
            {
                d.X = a.X + c.X - b.X;
                d.Y = a.Y + c.Y - b.Y;
            }
            if (Math.Abs(Dotmultiply(c, b, a)) < Common.Eps) // 说明a点是直角拐角处
            {
                d.X = c.X + b.X - a.X;
                d.Y = c.Y + b.Y - a.Y;
            }*/
            Vector2 v1 = new Vector2(new Point2(b.X, b.Y), new Point2(a.X, a.Y));
            Vector2 v2 = new Vector2(new Point2(b.X, b.Y), new Point2(c.X, c.Y));
            Vector2 v = v1 + v2;
            d.X = (int)v.End.X;
            d.Y = (int)v.End.Y;
            return d;
        }
        /// <summary>
        /// 已知平行四边形的三个顶点(a,b,c)，计算第四个顶点d的坐标. 注意：已知的三个顶点可以是有序的，b为中间点
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <returns></returns>
        public static PointF Rect4thP(PointF a, PointF b, PointF c)
        {
            PointF d = new PointF();
            /*
            if (Math.Abs(Dotmultiply(a, b, c)) < Common.Eps) // 说明c点是直角拐角处
            {
                d.X = a.X + b.X - c.X;
                d.Y = a.Y + b.Y - c.Y;
            }
            if (Math.Abs(Dotmultiply(a, c, b)) < Common.Eps) // 说明b点是直角拐角处
            {
                d.X = a.X + c.X - b.X;
                d.Y = a.Y + c.Y - b.Y;
            }
            if (Math.Abs(Dotmultiply(c, b, a)) < Common.Eps) // 说明a点是直角拐角处
            {
                d.X = c.X + b.X - a.X;
                d.Y = c.Y + b.Y - a.Y;
            }*/
            Vector2 v1 = new Vector2(new Point2(b.X, b.Y), new Point2(a.X, a.Y));
            Vector2 v2 = new Vector2(new Point2(b.X, b.Y), new Point2(c.X, c.Y));
            Vector2 v = v1 + v2;
            d.X = v.End.X;
            d.Y = v.End.Y;
            return d;
        }
        /// <summary>
        /// 返回点集决定的矩形区域，以及点集中距离最远的两点的方向
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="pp"></param>
        /// <returns></returns>
        public static Rectangle GetBoundFromPs(Point[] ps, out Point[] pp)
        {
            if (ps == null || ps.Length < 2)
            {
                throw new Exception("points is null or ps.Length < 2 !");
            }
            Point startP = Point.Empty;
            Point endP = Point.Empty;
            pp = new Point[2];

            if (ps.Length > 2)
            {
                pp = LinkMidPoint(ps);
            }
            else
            {
                pp[0] = ps[0];
                pp[1] = ps[1];
            }

            Point[] minmax1 = GetLineMinMaxX(ps);
            Point[] minmax2 = GetLineMinMaxY(ps);
            float xl = 0;
            if (Math.Abs(pp[1].X - pp[0].X) < Common.Eps)
                xl = float.MaxValue;
            else
                xl = (pp[0].Y - pp[1].Y) / (pp[1].X - pp[0].X);
            if (xl >= -1 && xl < 1)
            {
                startP = minmax1[0];
                endP = minmax1[1];
            }
            else if (xl < -1 || xl >= 1)
            {
                startP = minmax2[0];
                endP = minmax2[1];
            }
            /*if (startP.X > endP.X && startP.Y > endP.Y)
            {
                Point tmp = startP;
                startP = endP;
                endP = tmp;
            }*/
            int width = Math.Abs(endP.X - startP.X);
            int height = Math.Abs(endP.Y - startP.Y);
            return new Rectangle(startP.X, startP.Y, width, height);
        }
        /// <summary>
        /// 返回点集决定的矩形区域，以及点集中距离最远的两点的方向
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="pp"></param>
        /// <returns></returns>
        public static RectangleF GetBoundFromPs(PointF[] ps, out PointF[] pp)
        {
            if (ps == null || ps.Length < 2)
            {
                throw new Exception("points is null or ps.Length < 2 !");
            }
            PointF startP = PointF.Empty;
            PointF endP = PointF.Empty;
            pp = new PointF[2];

            if (ps.Length > 2)
            {
                pp = LinkMidPoint(ps);
            }
            else
            {
                pp[0] = ps[0];
                pp[1] = ps[1];
            }

            PointF[] minmax1 = GetLineMinMaxX(ps);
            PointF[] minmax2 = GetLineMinMaxY(ps);
            float xl = 0;
            if (Math.Abs(pp[1].X - pp[0].X) < Common.Eps)
                xl = float.MaxValue;
            else
                xl = (pp[0].Y - pp[1].Y) / (pp[1].X - pp[0].X);
            if (xl >= -1 && xl < 1)
            {
                startP = minmax1[0];
                endP = minmax1[1];
            }
            else if (xl < -1 || xl >= 1)
            {
                startP = minmax2[0];
                endP = minmax2[1];
            }
            /*if (startP.X > endP.X && startP.Y > endP.Y)
            {
                Point tmp = startP;
                startP = endP;
                endP = tmp;
            }*/
            float width = Math.Abs(endP.X - startP.X);
            float height = Math.Abs(endP.Y - startP.Y);
            return new RectangleF(startP.X, startP.Y, width, height);
        }
        /// <summary>
        /// 获取点集在X方向的最大最小值
        /// </summary>
        /// <param name="l"></param>
        /// <returns></returns>
        public static Point[] GetLineMinMaxX(Point[] l)
        {
            Point[] dd = new Point[2];
            if (l != null && l.Length > 0)
            {
                Point max = l[0];
                Point min = l[0];
                foreach (Point p in l)
                {
                    if (p.X < min.X)
                    {
                        min = p;
                    }
                    if (p.X > max.X)
                    {
                        max = p;
                    }
                }
                dd[0] = min;
                dd[1] = max;
            }
            return dd;
        }
        /// <summary>
        /// 获取点集在X方向的最大最小值
        /// </summary>
        /// <param name="l"></param>
        /// <returns></returns>
        public static PointF[] GetLineMinMaxX(PointF[] l)
        {
            PointF[] dd = new PointF[2];
            if (l != null && l.Length > 0)
            {
                PointF max = l[0];
                PointF min = l[0];
                foreach (PointF p in l)
                {
                    if (p.X < min.X)
                    {
                        min = p;
                    }
                    if (p.X > max.X)
                    {
                        max = p;
                    }
                }
                dd[0] = min;
                dd[1] = max;
            }
            return dd;
        }
        /// <summary>
        /// 获取点集在Y方向的最大最小值
        /// </summary>
        /// <param name="l"></param>
        /// <returns></returns>
        public static Point[] GetLineMinMaxY(Point[] l)
        {
            Point[] dd = new Point[2];
            if (l != null && l.Length > 0)
            {
                Point max = l[0];
                Point min = l[0];
                foreach (Point p in l)
                {
                    if (p.Y < min.Y)
                    {
                        min = p;
                    }
                    if (p.Y > max.Y)
                    {
                        max = p;
                    }
                }
                dd[0] = min;
                dd[1] = max;
            }
            return dd;
        }
        /// <summary>
        /// 获取点集在Y方向的最大最小值
        /// </summary>
        /// <param name="l"></param>
        /// <returns></returns>
        public static PointF[] GetLineMinMaxY(PointF[] l)
        {
            PointF[] dd = new PointF[2];
            if (l != null && l.Length > 0)
            {
                PointF max = l[0];
                PointF min = l[0];
                foreach (PointF p in l)
                {
                    if (p.Y < min.Y)
                    {
                        min = p;
                    }
                    if (p.Y > max.Y)
                    {
                        max = p;
                    }
                }
                dd[0] = min;
                dd[1] = max;
            }
            return dd;
        }
        /// <summary>
        /// 连中点
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static Point[] LinkMidPoint(Point[] ps)//连中点
        {
            if (ps.Length > 10000)
                throw new Exception("The points count > 10000, it's too more points !");
            int midPoints = ps.Length - 1;
            Point[] tmpMidPs = new Point[midPoints];
            for (int i = 0; i < midPoints; i++)
            {
                tmpMidPs[i] = GetMidPoint(ps[i], ps[i + 1]);
            }
            if (midPoints > 2)
            {
                tmpMidPs = LinkMidPoint(tmpMidPs);
            }
            return tmpMidPs;
        }
        /// <summary>
        /// 连中点
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static PointF[] LinkMidPoint(PointF[] ps)//连中点
        {
            if (ps.Length > 10000)
                throw new Exception("The points count > 10000, it's too more points !");
            int midPoints = ps.Length - 1;
            PointF[] tmpMidPs = new PointF[midPoints];
            for (int i = 0; i < midPoints; i++)
            {
                tmpMidPs[i] = GetMidPoint(ps[i], ps[i + 1]);
            }
            if (midPoints > 2)
            {
                tmpMidPs = LinkMidPoint(tmpMidPs);
            }
            return tmpMidPs;
        }
        /// <summary>
        /// 获取两点的中点
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static Point GetMidPoint(Point p1, Point p2)
        {
            Point mid = Point.Empty;
            try//避免堆栈溢出造成程序终止
            {
                mid = new Point((p1.X + p2.X) / 2, (p1.Y + p2.Y) / 2);
            }
            catch
            { }
            return mid;
        }
        /// <summary>
        /// 获取两点的中点
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static PointF GetMidPoint(PointF p1, PointF p2)
        {
            PointF mid = PointF.Empty;
            try//避免堆栈溢出造成程序终止
            {
                mid = new PointF((p1.X + p2.X) / 2, (p1.Y + p2.Y) / 2);
            }
            catch
            { }
            return mid;
        }
        /// <summary>
        /// 带区域的Roberts算子边缘检测
        /// </summary>
        /// <param name="b"></param>
        /// <param name="rect"></param>
        /// <returns></returns>
        public static Bitmap GetBmpByRoberts(Bitmap b, Rectangle rect)
        {
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            Bitmap dstImage = new Bitmap(rect.Width, rect.Height, b.PixelFormat);
            BitmapData srcData = b.LockBits(rect, ImageLockMode.ReadOnly, b.PixelFormat);
            BitmapData dstData = dstImage.LockBits(rect, ImageLockMode.WriteOnly, b.PixelFormat);
            int stride1 = srcData.Stride;
            int offset1 = stride1 - BPP * rect.Width;
            int stride2 = dstData.Stride;
            int offset2 = stride2 - BPP * rect.Width;
            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;
                int A, B;   // A(x-1, y-1)    B(x, y-1)
                int C, D;   // C(x-1,   y)    D(x,   y)
                // 不处理最上边
                src += stride1;
                dst += stride2;
                for (int y = 1; y < rect.Height; y++)
                {
                    //不处理最左边
                    src += BPP;
                    dst += BPP;
                    for (int x = 1; x < rect.Width; x++)
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            A = src[i - stride1 - BPP];
                            B = src[i - stride1];
                            C = src[i - BPP];
                            D = src[i];
                            //int tmp = (int)(Math.Sqrt((A - D) * (A - D) + (B - C) * (B - C))) + 128;
                            //tmp = tmp > 255 ? 255 : tmp;
                            dst[i] = (byte)(Math.Sqrt((A - D) * (A - D) + (B - C) * (B - C)));
                        } // i
                        if (BPP > 3)
                            dst[3] = src[3];
                        src += BPP;
                        dst += BPP;
                    } // x
                    src += offset1;
                    dst += offset2;
                } // y
            }
            b.UnlockBits(srcData);
            dstImage.UnlockBits(dstData);

            //dstImage = Effect.SharpenMore(dstImage);

            return dstImage;
        }
        /// <summary>
        /// 带区域的Roberts算子二值边缘检测
        /// </summary>
        /// <param name="b"></param>
        /// <param name="rect"></param>
        /// <returns></returns>
        public static Bitmap GetBitBmpByRoberts(Bitmap b, Rectangle rect)
        {
            byte threshold = OTSUThresholdValue(b.Clone(rect, b.PixelFormat));
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            Bitmap dstImage = new Bitmap(rect.Width, rect.Height, b.PixelFormat);
            BitmapData srcData = b.LockBits(rect, ImageLockMode.ReadOnly, b.PixelFormat);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, rect.Width, rect.Height), ImageLockMode.WriteOnly, b.PixelFormat);
            int stride1 = srcData.Stride;
            int offset1 = stride1 - BPP * rect.Width;
            int stride2 = dstData.Stride;
            int offset2 = stride2 - BPP * rect.Width;
            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;
                int A, B;   // A(x-1, y-1)    B(x, y-1)
                int C, D;   // C(x-1,   y)    D(x,   y)
                // 不处理最上边
                src += stride1;
                dst += stride2;
                for (int y = 1; y < rect.Height; y++)
                {
                    //不处理最左边
                    src += BPP;
                    dst += BPP;
                    for (int x = 1; x < rect.Width; x++)
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            A = src[i - stride1 - BPP];
                            B = src[i - stride1];
                            C = src[i - BPP];
                            D = src[i];
                            int tmp = (int)(Math.Sqrt((A - D) * (A - D) + (B - C) * (B - C)));
                            tmp = tmp >= threshold ? 255 : 0;
                            dst[i] = (byte)tmp;
                        } // i
                        if (BPP > 3)
                            dst[3] = src[3];
                        src += BPP;
                        dst += BPP;
                    } // x
                    src += offset1;
                    dst += offset2;
                } // y
            }
            b.UnlockBits(srcData);
            dstImage.UnlockBits(dstData);

            //dstImage = Effect.SharpenMore(dstImage);

            return dstImage;
        }
        /// <summary>
        /// 带区域的Roberts算子二值边缘检测
        /// </summary>
        /// <param name="b"></param>
        /// <param name="rect"></param>
        /// <param name="threshold"></param>
        /// <param name="needThinning"></param>
        /// <returns></returns>
        public static Bitmap GetBitBmpByRoberts(Bitmap b, Rectangle rect, byte threshold, bool needThinning)
        {
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            Bitmap dstImage = new Bitmap(rect.Width, rect.Height, b.PixelFormat);
            BitmapData srcData = b.LockBits(rect, ImageLockMode.ReadOnly, b.PixelFormat);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, rect.Width, rect.Height), ImageLockMode.WriteOnly, b.PixelFormat);
            int stride1 = srcData.Stride;
            int offset1 = stride1 - BPP * rect.Width;
            int stride2 = dstData.Stride;
            int offset2 = stride2 - BPP * rect.Width;
            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;
                int A, B;   // A(x-1, y-1)    B(x, y-1)
                int C, D;   // C(x-1,   y)    D(x,   y)
                // 不处理最上边
                src += stride1;
                dst += stride2;
                for (int y = 1; y < rect.Height; y++)
                {
                    //不处理最左边
                    src += BPP;
                    dst += BPP;
                    for (int x = 1; x < rect.Width; x++)
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            A = src[i - stride1 - BPP];
                            B = src[i - stride1];
                            C = src[i - BPP];
                            D = src[i];
                            int tmp = (int)(Math.Sqrt((A - D) * (A - D) + (B - C) * (B - C)));
                            tmp = tmp >= threshold ? 255 : 0;
                            dst[i] = (byte)tmp;
                        } // i
                        if (BPP > 3)
                            dst[3] = src[3];
                        src += BPP;
                        dst += BPP;
                    } // x
                    src += offset1;
                    dst += offset2;
                } // y
            }
            b.UnlockBits(srcData);
            dstImage.UnlockBits(dstData);

            if (needThinning)
                dstImage = ErosionCross(dstImage, 1, threshold);
            //dstImage = Effect.SharpenMore(dstImage);

            return dstImage;
        }
        /// <summary>
        /// Roberts算子边缘检测
        /// </summary>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Bitmap GetBmpByRoberts(Bitmap b)
        {
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            Bitmap dstImage = new Bitmap(width, height, b.PixelFormat);
            BitmapData srcData = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, b.PixelFormat);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, b.PixelFormat);
            int stride = srcData.Stride;
            int offset = stride - width * BPP;
            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;
                int A, B;   // A(x-1, y-1)    B(x, y-1)
                int C, D;   // C(x-1,   y)    D(x,   y)
                // 指向第一行
                src += stride;
                dst += stride;
                // 不处理最上边和最左边
                for (int y = 1; y < height; y++)
                {
                    // 指向每行第一列
                    src += BPP;
                    dst += BPP;
                    for (int x = 1; x < width; x++)
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            A = src[i - stride - BPP];
                            B = src[i - stride];
                            C = src[i - BPP];
                            D = src[i];
                            //int tmp = (int)(Math.Sqrt((A - D) * (A - D) + (B - C) * (B - C))) + 128;
                            //tmp = tmp > 255 ? 255 : tmp;
                            dst[i] = (byte)(Math.Sqrt((A - D) * (A - D) + (B - C) * (B - C)));
                        } // i
                        if (BPP > 3)
                            dst[3] = src[3];
                        src += BPP;
                        dst += BPP;
                    } // x
                    src += offset;
                    dst += offset;
                } // y
            }
            b.UnlockBits(srcData);
            dstImage.UnlockBits(dstData);
            return dstImage;
        }
        /// <summary>
        /// Roberts算子二值边缘检测
        /// </summary>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Bitmap GetBitBmpByRoberts(Bitmap b)
        {
            byte threshold = OTSUThresholdValue(b);
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            Bitmap dstImage = new Bitmap(width, height, b.PixelFormat);
            BitmapData srcData = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, b.PixelFormat);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, b.PixelFormat);
            int stride1 = srcData.Stride;
            int offset1 = stride1 - width * BPP;
            int stride2 = dstData.Stride;
            int offset2 = stride2 - width * BPP;
            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;
                int A, B;   // A(x-1, y-1)    B(x, y-1)
                int C, D;   // C(x-1,   y)    D(x,   y)
                // 指向第一行
                src += stride1;
                dst += stride2;
                // 不处理最上边和最左边
                for (int y = 1; y < height; y++)
                {
                    // 指向每行第一列
                    src += BPP;
                    dst += BPP;
                    for (int x = 1; x < width; x++)
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            A = src[i - stride1 - BPP];
                            B = src[i - stride1];
                            C = src[i - BPP];
                            D = src[i];
                            int tmp = (int)(Math.Sqrt((A - D) * (A - D) + (B - C) * (B - C)));
                            tmp = tmp >= threshold ? 255 : 0;
                            dst[i] = (byte)tmp;
                        } // i
                        if (BPP > 3)
                            dst[3] = src[3];
                        src += BPP;
                        dst += BPP;
                    } // x
                    src += offset1;
                    dst += offset2;
                } // y
            }
            b.UnlockBits(srcData);
            dstImage.UnlockBits(dstData);
            return dstImage;
        }
        /// <summary>
        /// Roberts算子二值边缘检测
        /// </summary>
        /// <param name="b"></param>
        /// <param name="threshold"></param>
        /// <param name="needThinning"></param>
        /// <returns></returns>
        public static Bitmap GetBitBmpByRoberts(Bitmap b, byte threshold, bool needThinning)
        {
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            Bitmap dstImage = new Bitmap(width, height, b.PixelFormat);
            BitmapData srcData = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, b.PixelFormat);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, b.PixelFormat);
            int stride1 = srcData.Stride;
            int offset1 = stride1 - width * BPP;
            int stride2 = dstData.Stride;
            int offset2 = stride2 - width * BPP;
            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;
                int A, B;   // A(x-1, y-1)    B(x, y-1)
                int C, D;   // C(x-1,   y)    D(x,   y)
                // 指向第一行
                src += stride1;
                dst += stride2;
                // 不处理最上边和最左边
                for (int y = 1; y < height; y++)
                {
                    // 指向每行第一列
                    src += BPP;
                    dst += BPP;
                    for (int x = 1; x < width; x++)
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            A = src[i - stride1 - BPP];
                            B = src[i - stride1];
                            C = src[i - BPP];
                            D = src[i];
                            int tmp = (int)(Math.Sqrt((A - D) * (A - D) + (B - C) * (B - C)));
                            tmp = tmp >= threshold ? 255 : 0;
                            dst[i] = (byte)tmp;
                        } // i
                        if (BPP > 3)
                            dst[3] = src[3];
                        src += BPP;
                        dst += BPP;
                    } // x
                    src += offset1;
                    dst += offset2;
                } // y
            }
            b.UnlockBits(srcData);
            dstImage.UnlockBits(dstData);
            if (needThinning)
                dstImage = ErosionCross(dstImage, 1, threshold);
            return dstImage;
        }
        /// <summary>
        /// 亮化
        /// </summary>
        /// <param name="color"></param>
        public void Brighten(Color color)
        {
            if (bmp != null)
            {
                int[] co = new int[4];
                co[0] = color.B;
                co[1] = color.G;
                co[2] = color.R;
                co[3] = color.A;
                BitmapData data = bmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, bmp.PixelFormat);
                unsafe
                {
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride - BPP * width;
                    int pixel = 0;
                    for (int y = 0; y < height; y++)
                    {
                        for (int x = 0; x < width; x++)
                        {
                            // 处理像素 B, G, R, A 分量
                            for (int i = 0; i < 4; i++)
                            {
                                pixel = p[i] + co[i];
                                if (pixel < 0) pixel = 0;
                                if (pixel > 255) pixel = 255;
                                p[i] = (byte)pixel;
                            } // i
                            p += BPP;
                        }  // x
                        p += offset;
                    } // y
                }
                bmp.UnlockBits(data);
            }
        }
        /// <summary>
        /// 暗化
        /// </summary>
        /// <param name="color"></param>
        public void Darken(Color color)
        {
            if (bmp != null)
            {
                int[] co = new int[4];
                co[0] = color.B;
                co[1] = color.G;
                co[2] = color.R;
                co[3] = color.A;
                BitmapData data = bmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, bmp.PixelFormat);
                unsafe
                {
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride - BPP * width;
                    int pixel = 0;
                    for (int y = 0; y < height; y++)
                    {
                        for (int x = 0; x < width; x++)
                        {
                            // 处理像素 B, G, R, A 分量
                            for (int i = 0; i < 4; i++)
                            {
                                pixel = p[i] - co[i];
                                if (pixel < 0) pixel = 0;
                                if (pixel > 255) pixel = 255;
                                p[i] = (byte)pixel;
                            } // i
                            p += BPP;
                        }  // x
                        p += offset;
                    } // y
                }
                bmp.UnlockBits(data);
            }
        }
        /*
        public Color GetPixel(int x, int y)
        {
            Color color = Color.Empty;
            if (bmp != null)
            {
                int width = bmp.Width;
                int height = bmp.Height;
                BitmapData data = bmp.LockBits(new Rectangle(x, y, 1, 1), ImageLockMode.ReadWrite, PixelFormat.Format32bppArgb);
                unsafe
                {
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride * y + BPP * x;
                    p += offset;
                    byte aa = 0, rr = 0, gg = 0, bb = 0;
                    // 处理像素 B, G, R ,A 三分量
                    bb = p[0];
                    if (bb < 0) bb = 0;
                    if (bb > 255) bb = 255;
                    gg = p[1];
                    if (gg < 0) gg = 0;
                    if (gg > 255) gg = 255;
                    rr = p[2];
                    if (rr < 0) rr = 0;
                    if (rr > 255) rr = 255;
                    aa = p[3];
                    if (aa < 0) aa = 0;
                    if (aa > 255) aa = 255;
                    color = Color.FromArgb(aa, rr, gg, bb);
                }
                bmp.UnlockBits(data);
            }
            return color;
        }
        public void SetPixel(int x, int y, Color color)
        {
            if (bmp != null)
            {
                int width = bmp.Width;
                int height = bmp.Height;
                BitmapData data = bmp.LockBits(new Rectangle(x, y, 1, 1), ImageLockMode.ReadWrite, PixelFormat.Format32bppArgb);
                unsafe
                {
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride * y + BPP * x;
                    p += offset;
                    // 处理像素 B, G, R ,A 三分量
                    p[0] = color.B;
                    if (p[0] < 0) p[0] = 0;
                    if (p[0] > 255) p[0] = 255;
                    p[1] = color.G;
                    if (p[1] < 0) p[1] = 0;
                    if (p[1] > 255) p[1] = 255;
                    p[2] = color.R;
                    if (p[2] < 0) p[2] = 0;
                    if (p[2] > 255) p[2] = 255;
                    p[3] = color.A;
                    if (p[3] < 0) p[3] = 0;
                    if (p[3] > 255) p[3] = 255;
                }
                bmp.UnlockBits(data);
            }
        }
        */
        /// <summary>
        /// 灰度化
        /// </summary>
        /// <param name="bmp"></param>
        /// <returns></returns>
        public static Bitmap Gray(Bitmap bmp)
        {
            Bitmap dstBmp = null;
            if (bmp != null)
            {
                dstBmp = (Bitmap)bmp.Clone();
                int width = dstBmp.Width;
                int height = dstBmp.Height;
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                BitmapData data = dstBmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, bmp.PixelFormat);
                unsafe
                {
                    byte R, G, B, gray;
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride - BPP * width;
                    for (int y = 0; y < height; y++)
                    {
                        for (int x = 0; x < width; x++)
                        {
                            R = p[2];
                            G = p[1];
                            B = p[0];
                            // gray = 0.3*R + 0.59*G + 0.11*B
                            gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                            p[0] = p[1] = p[2] = gray;
                            p += BPP;
                        } // x
                        p += offset;
                    } // y
                }
                dstBmp.UnlockBits(data);
            }
            return dstBmp;
        }
        /// <summary>
        /// 灰度化
        /// </summary>
        /// <returns></returns>
        public Bitmap Gray()
        {
            Bitmap dstBmp = null;
            if (bmp != null)
            {
                dstBmp = (Bitmap)bmp.Clone();
                int width = dstBmp.Width;
                int height = dstBmp.Height;
                int BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                BitmapData data = dstBmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, bmp.PixelFormat);
                unsafe
                {
                    byte R, G, B, gray;
                    byte* p = (byte*)data.Scan0;
                    int offset = data.Stride - BPP * width;
                    for (int y = 0; y < height; y++)
                    {
                        for (int x = 0; x < width; x++)
                        {
                            R = p[2];
                            G = p[1];
                            B = p[0];
                            // gray = 0.3*R + 0.59*G + 0.11*B
                            gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                            p[0] = p[1] = p[2] = gray;
                            p += BPP;
                        } // x
                        p += offset;
                    } // y
                }
                dstBmp.UnlockBits(data);
            }
            return dstBmp;
        }
        /// <summary>
        /// 将二值数组转换为双色位图流
        /// </summary>
        /// <param name="GrayArray">二值数组</param>
        /// <param name="pf">位图的像素格式</param>
        /// <param name="bgColor">背景色</param>
        /// <param name="fgColor">前景色</param>
        /// <returns></returns>
        public static Bitmap BinaryImage(byte[,] GrayArray, PixelFormat pf, Color bgColor, Color fgColor)
        {
            // 获取二维数组的宽高
            int width = GrayArray.GetLength(0);
            int height = GrayArray.GetLength(1);
            // 背景色
            byte bgAlpha = (byte)bgColor.A;
            byte bgRed = (byte)bgColor.R;
            byte bgGreen = (byte)bgColor.G;
            byte bgBlue = (byte)bgColor.B;
            // 前景色
            byte fgAlpha = (byte)fgColor.A;
            byte fgRed = (byte)fgColor.R;
            byte fgGreen = (byte)fgColor.G;
            byte fgBlue = (byte)fgColor.B;
            Bitmap b = new Bitmap(width, height, pf);
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        // 设置图像前景、背景颜色
                        if (GrayArray[x, y] > 128)
                        {
                            p[3] = fgAlpha;
                            p[2] = fgRed;
                            p[1] = fgGreen;
                            p[0] = fgBlue;
                        }
                        else
                        {
                            p[3] = bgAlpha;
                            p[2] = bgRed;
                            p[1] = bgGreen;
                            p[0] = bgBlue;
                        }
                        p += BPP;
                    } //  x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            return b;
        } // end of BinaryImage
        /// <summary>
        /// 将位图流转换为二维布尔数组
        /// </summary>
        /// <param name="b"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public static bool[,] Image2Array(Bitmap b, int blkORwht, byte threshold)
        {
            if (b == null)
                return null;
            int width = b.Width;
            int height = b.Height;
            // 申请一个二维数组
            bool[,] GrayArray = new bool[width, height];
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        bool val = false;
                        byte gray = 0;
                        if (p[0] == p[1] && p[1] == p[2])
                            gray = p[0];
                        else
                            gray = (byte)((19661 * p[2] + 38666 * p[1] + 7209 * p[0]) >> 16);
                        if (blkORwht == 0)
                        {
                            if (gray < threshold)
                                val = true;
                        }
                        else
                        {
                            if (gray > threshold)
                                val = true;
                        }
                        GrayArray[x, y] = val;
                        p += BPP;
                    } //  x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            return GrayArray;
        } // end of Image2Array
        /// <summary>
        /// 将位图流转换为二维数组
        /// </summary>
        /// <param name="b">位图流</param>
        /// <returns></returns>
        public static byte[,] Image2Array(Bitmap b)
        {
            if (b == null)
                return null;
            int width = b.Width;
            int height = b.Height;
            // 申请一个二维数组
            byte[,] GrayArray = new byte[width, height];
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        byte gray = 0;
                        if (p[0] == p[1] && p[1] == p[2])
                            gray = p[0];
                        else
                            gray = (byte)((19661 * p[2] + 38666 * p[1] + 7209 * p[0]) >> 16);
                        GrayArray[x, y] = gray;// p[0];
                        p += BPP;
                    } //  x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            return GrayArray;
        } // end of Image2Array
        /// <summary>
        /// 将位图流转换为二值数组
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="threshold">阈值</param>
        /// <param name="needGrayExtend">是否需要灰度拉伸</param>
        /// <returns></returns>
        public static byte[,] BinaryArray(Bitmap b, byte threshold, bool needGrayExtend)
        {
            int width = b.Width;
            int height = b.Height;
            //先将位图灰度化
            Bitmap bmp = null;

            /*if (needGray)
            {
                if (needGrayExtend)
                    bmp = GrayExtend(b, false);//灰度拉伸
                else
                    bmp = Gray(b);//灰度化
            }
            else
            {*/
            if (needGrayExtend)
                bmp = GrayExtend(b, true);//灰度拉伸
            else
                bmp = (Bitmap)b.Clone();
            //}
            // 将灰度图转化为灰度数组
            byte[,] GrayArray = Image2Array(bmp);
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    // 当前像素颜色灰度值与指定阈值相比较
                    if (GrayArray[x, y] >= threshold)
                        GrayArray[x, y] = 255;
                    else
                        GrayArray[x, y] = 0;
                } //  x
            } // y
            return GrayArray;
        } // end of BinaryArray
        /// <summary>
        /// 二值化
        /// </summary>
        /// <param name="bmp">位图流</param>
        /// <param name="threshold">阈值</param>
        /// <param name="needGrayExtend">是否需要灰度拉伸</param>
        /// <returns></returns>
        public static Bitmap Bitize(Bitmap bmp, byte threshold, bool needGrayExtend)
        {
            Bitmap dstImage = null;
            if (bmp != null)
            {
                byte[,] GrayArray = BinaryArray(bmp, threshold, needGrayExtend);
                dstImage = BinaryImage(GrayArray, bmp.PixelFormat, Color.Black, Color.White);
            }
            return dstImage;
        }

        /// <summary>
        /// 自适应二值化，但效果不是很理想
        /// </summary>
        /// <param name="b"></param>
        /// <param name="isTwoPeek"></param>
        /// <returns></returns>
        public Bitmap AutoFitBitize(Bitmap b, bool isTwoPeek)
        {
            if (isTwoPeek)
                currentThreshold = TwoPeakThresholdValue(b);
            else
                currentThreshold = OTSUThresholdValue(b);
            // 根据找到的谷值对图像进行二值化
            return Bitize(bmp, currentThreshold, false);
        } // end of AutoFitBitize

        private byte TwoPeakThresholdValue(Bitmap bitmap)
        {
            // 图像灰度化
            bmp = Gray(bitmap);
            // 建立直方图，并获取灰度统计信息
            Histogram histogram = new Histogram(bmp);
            int[] GrayLevel = histogram.Red.Value;
            int peak1, peak2, valley;
            int peak1Index, peak2Index, valleyIndex;
            // 取双峰
            peak1 = peak2 = GrayLevel[0];
            peak1Index = peak2Index = 0;
            for (int i = 1; i < 256; i++)
            {
                // 如果产生新的高峰，则将第一峰退居第二峰，新的高峰升为第一峰
                if (GrayLevel[i] > peak1)
                {
                    peak2 = peak1;
                    peak2Index = peak1Index;
                    peak1 = GrayLevel[i];
                    peak1Index = i;
                }
            } // i
            // 判断两个峰值索引
            int max = peak1Index;
            int min = peak2Index;
            if (max < min)
            {
                int t = max;
                max = min;
                min = t;
            }
            // 取峰谷
            valley = GrayLevel[min];
            valleyIndex = min;
            for (int i = min; i < max; i++)
            {
                if (GrayLevel[i] < valley)
                {
                    valley = GrayLevel[i];
                    valleyIndex = i;
                }
            } // i
            return (byte)valleyIndex;
        }
        /// <summary>
        /// 对图像进行平滑处理
        /// </summary>
        /// <param name="b">位图流</param>
        /// <returns></returns>
        public static Bitmap Smooth(Bitmap b)
        {
            //    1  1  1
            //    1  1  1
            //    1  1  1 / 9
            Matrix3x3 m = new Matrix3x3();
            m.Init(1);
            m.Scale = 9;
            return m.Convolute(b);
        } // end of Smooth

    } //end of AccessPixel
}