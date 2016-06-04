using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;

namespace KellImageProcess
{
    /// <summary>
    /// 数值分析算法
    /// </summary>
    public class Algorithm
    {
        //const double Common.Eps = 0.00001;

        /// <summary>
        /// 根据线段的始点、角度和长度，确定线段的正向终点
        /// </summary>
        /// <param name="startP"></param>
        /// <param name="angle">弧度为单位</param>
        /// <param name="length"></param>
        /// <returns></returns>
        public static PointF GetEndPofLineK(PointF startP, double angle, double length)
        {
            float k = (float)Math.Tan(angle);
            double d2 = length * length;
            float x = startP.X + (float)Math.Sqrt(d2 / (1 + k * k));
            float y = startP.Y - (float)Math.Sqrt(d2 - (d2 / (1 + k * k))) * Math.Sign(k);
            return new PointF(x, y);
        }
        /// <summary>
        /// 根据线段的始点、角度和长度，确定线段的反向终点
        /// </summary>
        /// <param name="startP"></param>
        /// <param name="angle">弧度为单位</param>
        /// <param name="length"></param>
        /// <returns></returns>
        public static PointF GetEndPofLine_K(PointF startP, double angle, double length)
        {
            float k = (float)Math.Tan(angle);
            double d2 = length * length;
            float x = startP.X - (float)Math.Sqrt(d2 / (1 + k * k));
            float y = startP.Y + (float)Math.Sqrt(d2 - (d2 / (1 + k * k))) * Math.Sign(k);
            return new PointF(x, y);
        }

        /// <summary>
        /// 求组合数
        /// </summary>
        /// <param name="m">取m个</param>
        /// <param name="n">总共n个</param>
        /// <returns></returns>
        private long nCr(int m, int n)
        {
            long v = 1, u = 1;
            for (int i = m; i >= 1; i--)
            {
                v = v * n;
                n--;
                u = u * i;
            }
            return v / u;
        }

        /// <summary>
        /// 求两点的方向角，范围[0, 360)
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public double GetFXJ(Point p1, Point p2)
        {
            double fxj = Math.Atan2(p1.Y - p2.Y, p2.X - p1.X);
            if (fxj < 0)
                fxj += 2 * Math.PI;
            return fxj * 180 / Math.PI;
        }
        /// <summary>
        /// 求有序三点的方向角差(单位：度)
        /// </summary>
        /// <param name="p1">始点</param>
        /// <param name="p2">中间点</param>
        /// <param name="p3">终点</param>
        /// <returns></returns>
        public double GetDiffFXJ(Point p1, Point p2, Point p3)
        {
            double fxj1 = GetFXJ(p1, p2);
            double fxj2 = GetFXJ(p2, p3);
            return fxj2 - fxj1;
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
                //MessageBox.Show(y.ToString() + "\n" + ys[0].ToString() + "\n" + s[4].ToString());
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
        /// 求转弯点
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public bool GetCornerPoint(List<Point> ps, out Point p)
        {
            p = Point.Empty;
            if (ps != null && ps.Count > 2)
            {
                p = ps[1];
                double maxJJ = Math.Abs(GetDiffFXJ(ps[0], ps[1], ps[2]));
                for (int i = 1; i < ps.Count - 2; i++)
                {
                    double current = Math.Abs(GetDiffFXJ(ps[i], ps[i + 1], ps[i + 2]));
                    if (current > maxJJ)
                    {
                        maxJJ = current;
                        p = ps[i + 1];
                    }
                }
                return true;
            }
            else
            {
                throw new KellImageProcessException("点集为空或点数少于3，无法求转弯点！");
            }
        }

        /// <summary>
        /// 求拐点
        /// </summary>
        /// <param name="ps">点集的个数要不能少于点集粒度(badDegree)的两倍</param>
        /// <param name="badDegree">点集粒度，常用[3,5,7,9,11]</param>
        /// <param name="gd"></param>
        /// <returns></returns>
        public bool GetGD(List<Point> ps, int badDegree, out List<Point> gd)
        {
            gd = new List<Point>();
            if (ps != null && ps.Count >= badDegree * 2)
            {
                int points = ps.Count / badDegree;
                List<Point> tmp = new List<Point>();
                for (int j = 0; j < badDegree; j++)
                {
                    tmp.Add(ps[j]);
                }
                double last = GetDSN(tmp, tmp[tmp.Count / 2], 2); ;
                for (int i = 1; i < points; i++)
                {
                    tmp.Clear();
                    for (int j = 0; j < badDegree; j++)
                    {
                        tmp.Add(ps[badDegree * i + j]);
                    }
                    double current = GetDSN(tmp, tmp[tmp.Count / 2], 2);
                    if (current * last < 0)
                    {
                        gd.Add(tmp[tmp.Count / 2]);
                        last = current;
                    }
                }
                return true;
            }
            else
            {
                throw new KellImageProcessException("点集为空或点数少于点集粒度的两倍" + Convert.ToString(badDegree * 2) + "，无法求拐点！");
            }
        }

        /// <summary>
        /// 计算点pt到直线ab的距离
        /// </summary>
        /// <param name="pt">直线外的一点</param>
        /// <param name="a">直线上的一点</param>
        /// <param name="b">直线上的另一点</param>
        /// <returns>pt点到直线上的距离</returns>
        public double Dlong(Point pt, Point a, Point b)
        {
            double A = a.Y - b.Y;//为x的系数
            double B = b.X - a.X;//为y的系数
            double C = a.X * b.Y - a.Y * b.X;//常数

            double d = 0.0;//距离
            d = Math.Abs(A * pt.X + B * pt.Y + C) / Math.Sqrt(A * A + B * B);
            return d;

        }

        /// <summary>
        /// 计算点pt到直线ab的距离
        /// </summary>
        /// <param name="pt">直线外的一点</param>
        /// <param name="a">直线上的一点</param>
        /// <param name="b">直线上的另一点</param>
        /// <returns>pt点到直线上的距离</returns>
        public double Dlong(PointF pt, PointF a, PointF b)
        {
            double A = a.Y - b.Y;//为x的系数
            double B = b.X - a.X;//为y的系数
            double C = a.X * b.Y - a.Y * b.X;//常数

            double d = 0.0;//距离
            d = Math.Abs(A * pt.X + B * pt.Y + C) / Math.Sqrt(A * A + B * B);
            return d;

        }
        /// <summary>
        /// 判断点是否在弧上
        /// </summary>
        /// <param name="p">一点</param>
        /// <param name="arcCenter">弧心</param>
        /// <param name="r">半径</param>
        /// <param name="startangle">开始角</param>
        /// <param name="sweepangle">扇过角</param>
        /// <returns>在弧上返回</returns>
        public bool BoolSpotOnArc(PointF p, PointF arcCenter, double r, double startangle, double sweepangle)
        {
            float d = (float)LineLength(p, arcCenter);
            double tep = d - (float)r;
            if (Math.Abs(tep) < 0.1)
            {
                PointF pt = new PointF();
                pt.X = p.X - arcCenter.X;
                pt.Y = p.Y - arcCenter.Y;
                double a = Math.Atan2(-pt.Y, pt.X);
                if (a > 0)
                    a = 360 - a * 180 / Math.PI;
                else if (a < 0)
                    a = -a * 180 / Math.PI;
                if ((a > startangle && a <= (startangle + sweepangle)) || (a <= startangle + sweepangle - 360))
                    return true;
                else
                    return false;
            }
            else
                return false;
        }

        /// <summary>
        /// 计算凸多边形的面积
        /// </summary>
        /// <param name="ps">多边形的各个顶点</param>
        /// <returns>多边形围成的面积</returns>
        public double PolygonDimension(List<PointF> ps)
        {
            //将点组成矩阵
            double sum = 0.0;
            for (int i = 0; i < ps.Count; i++)
            {
                double[,] ary = new double[3, 3];
                ary[0, 0] = 0;
                ary[0, 1] = ps[i].X;
                ary[0, 2] = ps[(i + 1) % ps.Count].X;
                ary[1, 0] = 0;
                ary[1, 1] = ps[i].Y;
                ary[1, 2] = ps[(i + 1) % ps.Count].Y;
                ary[2, 0] = 1;
                ary[2, 1] = 1;
                ary[2, 2] = 1;
                sum += MatrixDeterminant(ary) * 0.5;
            }
            return Math.Abs(sum);
        }

        /// <summary>
        /// 计算凸多边形的面积
        /// </summary>
        /// <param name="ps">多边形的各个顶点</param>
        /// <returns>多边形围成的面积</returns>
        public double PolygonDimension(List<Point> ps)
        {
            //将点组成矩阵
            double sum = 0.0;
            for (int i = 0; i < ps.Count; i++)
            {
                double[,] ary = new double[3, 3];
                ary[0, 0] = 0;
                ary[0, 1] = ps[i].X;
                ary[0, 2] = ps[(i + 1) % ps.Count].X;
                ary[1, 0] = 0;
                ary[1, 1] = ps[i].Y;
                ary[1, 2] = ps[(i + 1) % ps.Count].Y;
                ary[2, 0] = 1;
                ary[2, 1] = 1;
                ary[2, 2] = 1;
                sum += MatrixDeterminant(ary) * 0.5;
            }
            return Math.Abs(sum);
        }

        /// <summary>
        /// 计算矩阵的行列式值
        /// </summary>
        /// <param name="ary">矩阵</param>
        /// <returns></returns>
        private double MatrixDeterminant(double[,] ary)
        {
            int nRow = ary.GetLength(0);
            double sum = 0.0;
            if (nRow == 1)
            {
                sum = ary[0, 0];

            }
            else if (nRow == 2)
            {
                sum = ary[0, 0] * ary[1, 1] - ary[0, 1] * ary[1, 0];
            }
            else
            {
                for (int i = 0; i < nRow; i++)
                {
                    double[,] newary = new double[nRow - 1, nRow - 1];
                    for (int k = 1; k < nRow; k++)
                    {
                        bool fag = false;
                        for (int l = 0; l < nRow; l++)
                        {
                            if (l == i)
                                fag = true;
                            else if (fag && l != i)
                                newary[k - 1, l - 1] = ary[k, l];
                            else
                                newary[k - 1, l] = ary[k, l];
                        }
                    }
                    //将重组的计算行列式
                    sum = sum + (Math.Pow(-1.0, i) * ary[0, i] * MatrixDeterminant(newary));
                }
            }
            return sum;
        }
        /// <summary>
        /// 得到线段上的所有的点
        /// </summary>
        /// <param name="lt1">线段上的一点</param>
        /// <param name="lt2">线段上的另一点</param>
        /// <returns>点集</returns>
        public PointF[] AllLinePoints(PointF lt1, PointF lt2)
        {
            List<PointF> lf = new List<PointF>();
            float dx = lt2.X - lt1.X;
            float dy = lt2.Y - lt1.Y;
            PointF p = new PointF();
            p.X = lt1.X < lt2.X ? lt1.X : lt2.X;
            p.Y = lt1.Y < lt2.Y ? lt1.Y : lt2.Y;
            for (int i = 0; i < Math.Abs(dx); i++)
            {
                for (int j = 0; j < Math.Abs(dy); j++)
                {
                    PointF p1 = new PointF();
                    p1.X = p.X + i;
                    p1.Y = p.Y + j;
                    double tep = Dlong(p1, lt1, lt2);
                    if (tep < Math.Sqrt(2) / 2.0)
                        lf.Add(p1);
                }
            }
            PointF[] pfs = new PointF[lf.Count];
            pfs = lf.ToArray();
            return pfs;
        }

        /// <summary>
        /// 得到在圆上的所有点
        /// </summary>
        /// <param name="c">圆心</param>
        /// <param name="r">半径</param>
        /// <returns>点集</returns>
        public PointF[] AllPointPoints(PointF c, double r)
        {
            List<PointF> lf = new List<PointF>();
            PointF p = c;
            p.X = (float)(p.X - r);
            p.Y = (float)(p.Y - r);
            for (int i = 0; i < 2 * r; i++)
            {
                for (int j = 0; j < 2 * r; j++)
                {
                    PointF tep = new PointF();
                    tep.X = p.X + j;
                    tep.Y = p.Y + i;
                    double d1 = LineLength(tep, c);
                    if (Math.Abs(d1 - r) < Math.Sqrt(2) / 2.0)
                    {
                        lf.Add(tep);
                    }
                }
            }
            return lf.ToArray();
        }


        /// <summary>
        /// 得到在弧上的所有点
        /// </summary>
        /// <param name="p">弧心</param>
        /// <param name="r">弧半径</param>
        /// <param name="startangle">开始角</param>
        /// <param name="sweepangle">扇过角</param>
        /// <returns>点集</returns>
        public PointF[] AllArcPoints(PointF p, double r, double startangle, double sweepangle)
        {
            List<PointF> lf = new List<PointF>();
            PointF pp = p;
            pp.X = (float)(pp.X - r);
            pp.Y = (float)(pp.Y - r);
            for (int i = 0; i < 2 * r; i++)
            {
                for (int j = 0; j < 2 * r; j++)
                {
                    PointF tep = new PointF();
                    tep.X = pp.X + j;
                    tep.Y = pp.Y + i;
                    double d1 = LineLength(tep, p);
                    if (Math.Abs(d1 - r) < Math.Sqrt(2) / 2.0)
                    {
                        PointF pk = new PointF();
                        pk.X = p.X - tep.X;
                        pk.Y = tep.Y - p.Y;
                        double b3 = Math.Atan2(pk.Y, pk.X);
                        b3 = 180 - b3 * 180 / Math.PI;
                        if ((startangle + sweepangle) > 360)
                        {
                            if (b3 >= startangle || b3 <= (startangle + sweepangle) - 360)
                            {
                                lf.Add(tep);
                            }
                        }
                        else
                        {
                            if (b3 >= startangle && b3 <= (sweepangle + sweepangle))
                            {
                                lf.Add(tep);
                            }

                        }
                    }
                }
            }
            return lf.ToArray();

        }
        /// <summary>
        /// 是否是凸多边形
        /// </summary>
        /// <param name="pf">点集</param>
        /// <returns>是true</returns>
        public bool BoolRaisedPolygon(PointF[] pf)
        {
            for (int i = 0; i < pf.Length; i++)
            {
                if (i == pf.Length - 1)
                {
                    List<PointF> tep = new List<PointF>();
                    for (int j = 1; j < pf.Length - 1; j++)
                    {
                        PointF d = Pvint(pf[j], pf[i], pf[0]);
                        d.X = pf[j].X - d.X;
                        d.Y = pf[j].Y - d.Y;
                        tep.Add(d);
                    }
                    //组成的向量两两比较
                    for (int w = 0; w < tep.Count - 1; w++)
                    {
                        if (tep[w].X * tep[w + 1].X < 0)
                        {
                            return false;
                        }
                    }
                }
                else
                {
                    List<PointF> tep = new List<PointF>();
                    for (int j = 0; j < pf.Length; j++)
                    {
                        bool fag = (j != i && j != i + 1);
                        if (fag)
                        {
                            PointF d = Pvint(pf[j], pf[i], pf[i + 1]);
                            d.X = pf[j].X - d.X;
                            d.Y = pf[j].Y - d.Y;
                            tep.Add(d);
                        }
                    }
                    for (int w = 0; w < tep.Count - 1; w++)
                    {
                        if (tep[w].X * tep[w + 1].X < 0)
                        {
                            return false;
                        }
                    }

                }
            }
            return true;
        }
        /// <summary>
        /// 所有的点顺时钟所在的直线
        /// </summary>
        /// <param name="a">点集</param>
        /// <returns></returns>
        public List<Linepoint> ToLine(List<PointF> a)
        {
            List<Linepoint> linepoint = new List<Linepoint>();
            for (int i = 0; i < a.Count - 1; i++)
            {
                Linepoint lp;
                lp.lt1 = a[i];
                lp.lt2 = a[i + 1];
                linepoint.Add(lp);
            }
            Linepoint lp1;
            lp1.lt1 = a[a.Count - 1];
            lp1.lt2 = a[0];
            linepoint.Add(lp1);
            return linepoint;
        }
        /// <summary>
        /// 计算凸多边形的内包圆
        /// </summary>
        /// <param name="pf"></param>
        /// <param name="c"></param>
        /// <param name="r"></param>
        public void RaisePolygonInCircle(List<PointF> pf, out PointF c, out double r)
        {
            double sumx = 0.0, sumy = 0.0;
            for (int i = 0; i < pf.Count; i++)
            {
                sumx += pf[i].X;
                sumy += pf[i].Y;
            }
            sumx = sumx / pf.Count;
            sumy = sumy / pf.Count;

            PointF p = new PointF();
            p.X = (float)(sumx);
            p.Y = (float)(sumy);
            double mind = Dlong(p, pf[0], pf[1]);
            for (int i = 0; i < pf.Count; i++)
            {
                if (i == pf.Count - 1)
                {
                    double tepd = Dlong(p, pf[i], pf[0]);
                    if (tepd < mind)
                        mind = tepd;
                }
                else
                {
                    double tepd = Dlong(p, pf[i], pf[i + 1]);
                    if (tepd < mind)
                        mind = tepd;
                }
            }
            List<Linepoint> linepoint = ToLine(pf);
            int s = 0;

            mind += Common.Eps;
            while (true)
            {
                bool fag = true;
                for (int j = 0; j < linepoint.Count; j++)
                {
                    double temp = Dlong(p, linepoint[j].lt1, linepoint[j].lt2);
                    if (temp < mind)
                    {
                        PointF p1 = Pvint(p, linepoint[j].lt1, linepoint[j].lt2);
                        PointF pk = new PointF();
                        pk.X = (float)((p.X - p1.X) / temp);
                        pk.Y = (float)((p.Y - p1.Y) / temp);
                        p.X = (float)(p.X + pk.X * 0.05);
                        p.Y = (float)(p.Y + pk.Y * 0.05);
                        fag = false;
                        s++;
                        break;
                    }
                }
                if (fag)
                {
                    mind = mind + 0.05;
                    s = 0;
                }
                if (s > 1000)
                {
                    //找一个最小的为半径
                    mind = Dlong(p, linepoint[0].lt1, linepoint[0].lt2);
                    for (int x = 0; x < linepoint.Count; x++)
                    {
                        double tep = Dlong(p, linepoint[x].lt1, linepoint[x].lt2);
                        if (tep < mind)
                            mind = tep;
                    }
                    break;
                }
            }
            c = p;
            r = mind;
        }
        /// <summary>
        /// 点集的外包正规矩形
        /// </summary>
        /// <param name="ps">点集</param>
        /// <returns></returns>
        public RectangleF RectangleOnMorePoint(List<PointF> ps)
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
            RectangleF re = new RectangleF(minx, miny, maxx - minx, maxy - miny);

            return re;
        }
        /// <summary>
        /// 点集的外包正规矩形
        /// </summary>
        /// <param name="ps">点集</param>
        /// <returns></returns>
        public Rectangle RectangleOnMorePoint(List<Point> ps)
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
            Rectangle re = new Rectangle(minx, miny, maxx - minx, maxy - miny);

            return re;
        }

        /// <summary>
        /// 两点间的距离
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public double LineLength(PointF p1, PointF p2)
        {
            return Math.Sqrt((p1.X - p2.X) * (p1.X - p2.X)) + ((p1.Y - p2.Y) * (p1.Y - p2.Y));
        }

        /// <summary>
        /// 两点间的距离
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public double LineLength(Point p1, Point p2)
        {
            return Math.Sqrt((p1.X - p2.X) * (p1.X - p2.X)) + ((p1.Y - p2.Y) * (p1.Y - p2.Y));
        }

        /// <summary>
        /// 两线段的交点
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <param name="p4"></param>
        /// <returns></returns>
        public PointF GetCrossPointOfTwoLineSegment(PointF p1, PointF p2, PointF p3, PointF p4)
        {
            double a1 = p2.Y - p1.Y;
            double a2 = p4.Y - p3.Y;
            double b1 = p1.X - p2.X;
            double b2 = p3.X - p4.X;
            double c1 = (p2.X * p1.Y) - (p1.X * p2.Y);
            double c2 = (p4.X * p3.Y) - (p3.X * p4.Y);
            if ((a1 * b2) == (a2 * b1))
            {
                return new PointF(float.MaxValue, float.MaxValue);
            }
            double x = ((c2 * b1) - (c1 * b2)) / ((a1 * b2) - (a2 * b1));
            double y = ((c1 * a2) - (c2 * a1)) / ((a1 * b2) - (b1 * a2));
            if ((((x <= Math.Max(p1.X, p2.X)) && (x >= Math.Min(p1.X, p2.X))) && (x <= Math.Max(p3.X, p4.X))) && (x >= Math.Min(p3.X, p4.X)))
            {
                return new PointF((float)x, (float)y);
            }
            return PointF.Empty;
        }

        /// <summary>
        /// 两线段的交点
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <param name="p4"></param>
        /// <returns></returns>
        public PointF GetCrossPointOfTwoLineSegment(Point p1, Point p2, Point p3, Point p4)
        {
            double a1 = p2.Y - p1.Y;
            double a2 = p4.Y - p3.Y;
            double b1 = p1.X - p2.X;
            double b2 = p3.X - p4.X;
            double c1 = (p2.X * p1.Y) - (p1.X * p2.Y);
            double c2 = (p4.X * p3.Y) - (p3.X * p4.Y);
            if ((a1 * b2) == (a2 * b1))
            {
                return new PointF(float.MaxValue, float.MaxValue);
            }
            double x = ((c2 * b1) - (c1 * b2)) / ((a1 * b2) - (a2 * b1));
            double y = ((c1 * a2) - (c2 * a1)) / ((a1 * b2) - (b1 * a2));
            if ((((x <= Math.Max(p1.X, p2.X)) && (x >= Math.Min(p1.X, p2.X))) && (x <= Math.Max(p3.X, p4.X))) && (x >= Math.Min(p3.X, p4.X)))
            {
                return new PointF((float)x, (float)y);
            }
            return PointF.Empty;
        }

        /// <summary>
        /// 两直线的交点
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <param name="p4"></param>
        /// <returns></returns>
        public PointF GetCrossPointOfTwoLine(Point p1, Point p2, Point p3, Point p4)
        {
            double a1 = p2.Y - p1.Y;
            double a2 = p4.Y - p3.Y;
            double b1 = p1.X - p2.X;
            double b2 = p3.X - p4.X;
            double c1 = (p2.X * p1.Y) - (p1.X * p2.Y);
            double c2 = (p4.X * p3.Y) - (p3.X * p4.Y);
            if ((a1 * b2) == (a2 * b1))
            {
                return PointF.Empty;
            }
            double x = ((c2 * b1) - (c1 * b2)) / ((a1 * b2) - (a2 * b1));
            double y = ((a1 * c2) - (a2 * c1)) / ((a2 * b1) - (a1 * b2));
            return new PointF((float)x, (float)y);
        }

        /// <summary>
        /// 两直线的交点
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <param name="p4"></param>
        /// <returns></returns>
        public PointF GetCrossPointOfTwoLine(PointF p1, PointF p2, PointF p3, PointF p4)
        {
            double a1 = p2.Y - p1.Y;
            double a2 = p4.Y - p3.Y;
            double b1 = p1.X - p2.X;
            double b2 = p3.X - p4.X;
            double c1 = (p2.X * p1.Y) - (p1.X * p2.Y);
            double c2 = (p4.X * p3.Y) - (p3.X * p4.Y);
            if (a1 * b2 == a2 * b1)//两直线平行
            {
                return PointF.Empty;
            }
            double x = (c2 * b1 - c1 * b2) / (a1 * b2 - a2 * b1);
            double y = (a1 * c2 - a2 * c1) / (a2 * b1 - a1 * b2);
            return new PointF((float)x, (float)y);
        }

        /// <summary>
        /// 直线与线段的交点
        /// </summary>
        /// <param name="p1">直线的一点</param>
        /// <param name="p2">直线的另一点</param>
        /// <param name="p3">线段的一点</param>
        /// <param name="p4">线段的另一点</param>
        /// <returns></returns>
        public PointF GetCrossPointOfLineToLineSegment(PointF p1, PointF p2, PointF p3, PointF p4)
        {
            PointF cross = GetCrossPointOfTwoLine(p1, p2, p3, p4);
            if (OnLineSEG(p3, p4, cross))
                return cross;
            else
                return PointF.Empty;
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
        public static double Multiply(Point sp, Point ep, Point op)
        {
            return ((sp.X - op.X) * (ep.Y - op.Y) - (ep.X - op.X) * (sp.Y - op.Y));
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
        public static double Multiply(PointF sp, PointF ep, PointF op)
        {
            return ((sp.X - op.X) * (ep.Y - op.Y) - (ep.X - op.X) * (sp.Y - op.Y));
        }

        /// <summary>
        /// 判断点p是否在线段(p1,p2)上，条件：(p在线段l所在的直线上)而且(点p在线段l的正规矩形内)
        /// </summary>
        /// <param name="p1">线段的一点</param>
        /// <param name="p2">线段的另一点</param>
        /// <param name="p"></param>
        /// <returns></returns>
        public static bool OnLineSEG(PointF p1, PointF p2, PointF p)
        {
            return (OnLine(p1, p2, p)) && (((p.X - p1.X) * (p.X - p2.X) <= 0) && ((p.Y - p1.Y) * (p.Y - p2.Y) <= 0));
        }
        /// <summary>
        /// 判断点p是否在直线(p1,p2)上
        /// </summary>
        /// <param name="p1">直线的一点</param>
        /// <param name="p2">直线的另一点</param>
        /// <param name="p"></param>
        /// <returns></returns>
        public static bool OnLine(PointF p1, PointF p2, PointF p)
        {
            return Multiply(p2, p, p1) == 0;
        }

        /// <summary>
        /// 直线与多边形的交点
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public List<PointF> CrossPointOfLineAndPolygon(PointF p1, PointF p2, PointF[] p)
        {
            int n = p.Length;
            List<PointF> list = new List<PointF>();
            for (int i = 0; i < n; i++)
            {
                PointF px = GetCrossPointOfTwoLine(p1, p2, p[i], p[(i + 1) % n]);
                if ((px != PointF.Empty) && (px != new PointF(float.MaxValue, float.MaxValue)))
                {
                    list.Add(px);
                }
            }
            return list;
        }

        /// <summary>
        /// 直线与矩形的交点
        /// </summary>
        /// <param name="pa"></param>
        /// <param name="pb"></param>
        /// <param name="rect"></param>
        /// <returns></returns>
        public List<PointF> CrossPointOfLineAndRect(PointF pa, PointF pb, RectangleF rect)
        {
            PointF ru = new PointF(rect.X + rect.Width, rect.Y);
            PointF rd = new PointF(rect.X + rect.Width, rect.Y + rect.Height);
            PointF ld = new PointF(rect.X, rect.Y + rect.Height);
            List<PointF> list = new List<PointF>();
            PointF p1 = GetCrossPointOfLineToLineSegment(pa, pb, rect.Location, ru);
            PointF p2 = GetCrossPointOfLineToLineSegment(pa, pb, ru, rd);
            PointF p3 = GetCrossPointOfLineToLineSegment(pa, pb, rd, ld);
            PointF p4 = GetCrossPointOfLineToLineSegment(pa, pb, ld, rect.Location);
            if ((p1 != PointF.Empty) && (p1 != new PointF(float.MaxValue, float.MaxValue)))
            {
                list.Add(p1);
            }
            if ((p2 != PointF.Empty) && (p2 != new PointF(float.MaxValue, float.MaxValue)))
            {
                list.Add(p2);
            }
            if ((p3 != PointF.Empty) && (p3 != new PointF(float.MaxValue, float.MaxValue)))
            {
                list.Add(p3);
            }
            if ((p4 != PointF.Empty) && (p4 != new PointF(float.MaxValue, float.MaxValue)))
            {
                list.Add(p4);
            }
            return list;
        }

        /// <summary>
        /// 直线与矩形的交点
        /// </summary>
        /// <param name="pa"></param>
        /// <param name="pb"></param>
        /// <param name="rect"></param>
        /// <returns></returns>
        public List<PointF> CrossPointOfLineAndRect(Point pa, Point pb, Rectangle rect)
        {
            Point ru = new Point(rect.X + rect.Width, rect.Y);
            Point rd = new Point(rect.X + rect.Width, rect.Y + rect.Height);
            Point ld = new Point(rect.X, rect.Y + rect.Height);
            List<PointF> list = new List<PointF>();
            PointF p1 = GetCrossPointOfLineToLineSegment(pa, pb, rect.Location, ru);
            PointF p2 = GetCrossPointOfLineToLineSegment(pa, pb, ru, rd);
            PointF p3 = GetCrossPointOfLineToLineSegment(pa, pb, rd, ld);
            PointF p4 = GetCrossPointOfLineToLineSegment(pa, pb, ld, rect.Location);
            if ((p1 != PointF.Empty) && (p1 != new PointF(float.MaxValue, float.MaxValue)))
            {
                list.Add(p1);
            }
            if ((p2 != PointF.Empty) && (p2 != new PointF(float.MaxValue, float.MaxValue)))
            {
                list.Add(p2);
            }
            if ((p3 != PointF.Empty) && (p3 != new PointF(float.MaxValue, float.MaxValue)))
            {
                list.Add(p3);
            }
            if ((p4 != PointF.Empty) && (p4 != new PointF(float.MaxValue, float.MaxValue)))
            {
                list.Add(p4);
            }
            return list;
        }

        /// <summary>
        /// 拟合直线算法.拟合直线形如：y = a*x+b;
        /// </summary>
        /// <param name="pf">拟合点集合</param>
        /// <param name="a">系数a</param>
        /// <param name="b">参数b</param>
        public void NiHeLine(List<PointF> pf, out float a, out float b)
        {
            int n = pf.Count;
            float[] pfx = new float[pf.Count];
            float[] pfy = new float[pf.Count];
            for (int i = 0; i < n; i++)
            {
                pfx[i] = pf[i].X;
                pfy[i] = pf[i].Y;
            }
            float sumx, sumy;
            float sumxy;
            float sumxx, sumyy;
            sumx = Sum(pfx);
            sumy = Sum(pfy);
            sumxx = PingFangHe(pfx);
            sumyy = PingFangHe(pfy);
            sumxy = Multiply(pf);

            a = (n * sumxy - sumx * sumy) / (n * sumxx - sumx * sumx);
            b = (1 / n) * sumy - (a / n) * sumx;
        }
        /// <summary>
        /// 拟合直线算法.拟合直线形如：y = a*x+b;
        /// </summary>
        /// <param name="ps">拟合点集合</param>
        /// <param name="a">系数a</param>
        /// <param name="b">参数b</param>
        public void NiHeLine(List<PointF> ps, out double a, out double b)
        {
            int n = ps.Count;
            double sumx = 0, sumy = 0;
            double sumxy = 0;
            double sumxx = 0, sumyy = 0;
            float minX = ps[0].X, maxY = ps[0].Y;
            for (int i = 0; i < n; i++)
            {
                sumx = sumx + ps[i].X;
                sumy = sumy - ps[i].Y;

                sumxx = sumxx + ps[i].X * ps[i].X;
                sumyy = sumyy + ps[i].Y * ps[i].Y;
                sumxy = sumxy - ps[i].X * ps[i].Y;
            }

            if (n * sumxx == sumx * sumx)
            {
                a = 0;
                b = 0;
            }
            else
            {
                a = (n * sumxy - sumx * sumy) / (n * sumxx - sumx * sumx);
                b = (sumxx * sumy - sumx * sumxy) / (n * sumxx - sumx * sumx);

            }

        }
        /// <summary>
        /// 拟合直线算法.拟合直线形如：y = a*x+b;
        /// </summary>
        /// <param name="ps">拟合点集合</param>
        /// <param name="a">系数a</param>
        /// <param name="b">参数b</param>
        public void NiHeLine(List<Point> ps, out double a, out double b)
        {
            int n = ps.Count;
            double sumx = 0, sumy = 0;
            double sumxy = 0;
            double sumxx = 0, sumyy = 0;
            float minX = ps[0].X, maxY = ps[0].Y;
            for (int i = 0; i < n; i++)
            {
                sumx = sumx + ps[i].X;
                sumy = sumy - ps[i].Y;

                sumxx = sumxx + ps[i].X * ps[i].X;
                sumyy = sumyy + ps[i].Y * ps[i].Y;
                sumxy = sumxy - ps[i].X * ps[i].Y;
            }

            if (n * sumxx == sumx * sumx)
            {
                a = 0;
                b = 0;
            }
            else
            {
                a = (n * sumxy - sumx * sumy) / (n * sumxx - sumx * sumx);
                b = (sumxx * sumy - sumx * sumxy) / (n * sumxx - sumx * sumx);

            }

        }

        /// <summary>
        /// 过直线外一点的在直线上的垂足
        /// </summary>
        /// <param name="pt">直线外的一点</param>
        /// <param name="lt1">直线上的一点</param>
        /// <param name="lt2">直线上的另一点</param>
        /// <returns></returns>
        public Point Pvint(Point pt, Point lt1, Point lt2)
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
        /// 过直线外一点的在直线上的垂足
        /// </summary>
        /// <param name="pt">直线外的一点</param>
        /// <param name="lt1">直线上的一点</param>
        /// <param name="lt2">直线上的另一点</param>
        /// <returns></returns>
        public PointF Pvint(PointF pt, PointF lt1, PointF lt2)
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
        /// 拟合直线算法.拟合直线形如：y = a*x+b;
        /// </summary>
        /// <param name="ps">拟合点集合</param>
        /// <param name="PointA">拟合的线段端点A</param>
        /// <param name="PointB">拟合的线段端点B</param>
        public void NiHeLine(List<PointF> ps, out PointF PointA, out PointF PointB)
        {
            /*//师
            int n = ps.Count;
            double sumx = 0, sumy = 0;
            double sumxy = 0;
            double sumxx = 0, sumyy = 0;
            float MinY = ps[0].Y, MaxY = ps[0].Y;
            int minY = 0, maxY = 0;

            double a, b;
            float MinX = ps[0].X;
            float MaxX = ps[0].X;
            int minX = 0, maxX = 0;
            for (int i = 0; i < n; i++)
            {
                sumx = sumx + ps[i].X;
                sumy = sumy - ps[i].Y;

                sumxx = sumxx + ps[i].X * ps[i].X;
                sumyy = sumyy + ps[i].Y * ps[i].Y;
                sumxy = sumxy - ps[i].X * ps[i].Y;

                if (ps[i].Y < minY)
                {
                    MinY = ps[i].Y;
                    minY = i;
                }
                if (ps[i].Y > maxY)
                {
                    MaxY = ps[i].Y;
                    maxY = i;
                }
                if (ps[i].X < MinX)
                {
                    MinX = ps[i].X;
                    minX = i;
                }
                if (ps[i].X > MaxX)
                {
                    MaxX = ps[i].X;
                    maxX = i;
                }
            }

            //a = (n * sumxy - sumx * sumy) / (n * sumxx - sumx * sumx);
            //b = (1 / n) * sumy - (a / n) * sumx;
            if (Math.Abs(n * sumxx - sumx * sumx) < 0.005)
            {
                PointA = ps[maxY];
                PointB = ps[minY];
            }
            else
            {
                a = (n * sumxy - sumx * sumy) / (n * sumxx - sumx * sumx);
                b = (sumxx * sumy - sumx * sumxy) / (n * sumxx - sumx * sumx);
                if (Math.Abs(a) < 0.01)
                {
                    PointA = ps[maxX];
                    PointB = ps[minX];
                }
                else if (Math.Abs(a) > 250)
                {
                    PointA = ps[maxY];
                    PointB = ps[minY];
                }
                else
                {
                    int Min = (int)(-(MinY + b) / a);
                    int Max = (int)(-(MaxY + b) / a);
                    PointA = new PointF(Min, MinY);
                    PointB = new PointF(Max, MaxY);
                }
            }*/
            PointA = PointF.Empty;
            PointB = PointF.Empty;
            if (ps.Count < 2)
            {
                //MessageBox.Show("点数少于2点，无法定线！");
                return;
            }
            int i;
            double x = 0.0;
            double y = 0.0;
            double xy = 0.0;
            double x2 = 0.0;
            for (i = 0; i < ps.Count; i++)
            {
                x += ps[i].X;
                y += ps[i].Y;
                xy += ps[i].X * ps[i].Y;
                x2 += ps[i].X * ps[i].X;
            }
            double a = ((ps.Count * xy) - (x * y)) / ((ps.Count * x2) - (x * x));
            double b = (y - (a * x)) / ps.Count;
            PointF p1 = new PointF(0, (float)b);
            PointF p2 = new PointF(1, (float)(a + b));
            RectangleF rect = RectangleOnMorePoint(ps);
            //rect = new RectangleF(rect.X - 1, rect.Y - 1, rect.Width + 2, rect.Height + 2);
            //List<PointF> cross = CrossPointOfLineAndRect(p1, p2, rect);
            //PointA = cross[0];
            //PointB = cross[1];
            List<PointF> p4 = new List<PointF>();
            PointF ru = new PointF(rect.X + rect.Width, rect.Y);
            PointF rd = new PointF(rect.X + rect.Width, rect.Y + rect.Height);
            PointF ld = new PointF(rect.X, rect.Y + rect.Height);
            p4.Add(rect.Location);
            p4.Add(ru);
            p4.Add(rd);
            p4.Add(ld);
            PointF[] pp1 = new PointF[2];
            pp1[0] = Pvint(p4[0], p1, p2);
            pp1[1] = Pvint(p4[2], p1, p2);
            PointF[] pp2 = new PointF[2];
            pp2[0] = Pvint(p4[1], p1, p2);
            pp2[1] = Pvint(p4[3], p1, p2);
            double t1 = LineLength(pp1[0], pp1[1]);
            double t2 = LineLength(pp2[0], pp2[1]);
            if (t1 > t2)
            {
                PointA = pp1[0];
                PointB = pp1[1];
            }
            else
            {
                PointA = pp2[0];
                PointB = pp2[1];
            }
        }

        /// <summary>
        /// 拟合直线算法.拟合直线形如：y = a*x+b;
        /// </summary>
        /// <param name="ps">拟合点集合</param>
        /// <param name="PointA">拟合的线段端点A</param>
        /// <param name="PointB">拟合的线段端点B</param>
        public void NiHeLine(List<Point> ps, out PointF PointA, out PointF PointB)
        {
            /*//师
            int n = ps.Count;
            double sumx = 0, sumy = 0;
            double sumxy = 0;
            double sumxx = 0, sumyy = 0;
            float MinY = ps[0].Y, MaxY = ps[0].Y;
            int minY = 0, maxY = 0;

            double a, b;
            float MinX = ps[0].X;
            float MaxX = ps[0].X;
            int minX = 0, maxX = 0;
            for (int i = 0; i < n; i++)
            {
                sumx = sumx + ps[i].X;
                sumy = sumy - ps[i].Y;

                sumxx = sumxx + ps[i].X * ps[i].X;
                sumyy = sumyy + ps[i].Y * ps[i].Y;
                sumxy = sumxy - ps[i].X * ps[i].Y;

                if (ps[i].Y < minY)
                {
                    MinY = ps[i].Y;
                    minY = i;
                }
                if (ps[i].Y > maxY)
                {
                    MaxY = ps[i].Y;
                    maxY = i;
                }
                if (ps[i].X < MinX)
                {
                    MinX = ps[i].X;
                    minX = i;
                }
                if (ps[i].X > MaxX)
                {
                    MaxX = ps[i].X;
                    maxX = i;
                }
            }

            //a = (n * sumxy - sumx * sumy) / (n * sumxx - sumx * sumx);
            //b = (1 / n) * sumy - (a / n) * sumx;
            if (Math.Abs(n * sumxx - sumx * sumx) < 0.005)
            {
                PointA = ps[maxY];
                PointB = ps[minY];
            }
            else
            {
                a = (n * sumxy - sumx * sumy) / (n * sumxx - sumx * sumx);
                b = (sumxx * sumy - sumx * sumxy) / (n * sumxx - sumx * sumx);
                if (Math.Abs(a) < 0.01)
                {
                    PointA = ps[maxX];
                    PointB = ps[minX];
                }
                else if (Math.Abs(a) > 250)
                {
                    PointA = ps[maxY];
                    PointB = ps[minY];
                }
                else
                {
                    int Min = (int)(-(MinY + b) / a);
                    int Max = (int)(-(MaxY + b) / a);
                    PointA = new PointF(Min, MinY);
                    PointB = new PointF(Max, MaxY);
                }
            }
            */
            PointA = PointF.Empty;
            PointB = PointF.Empty;
            if (ps.Count < 2)
            {
                //MessageBox.Show("点数少于2点，无法定线！");
                return;
            }
            int i;
            double x = 0.0;
            double y = 0.0;
            double xy = 0.0;
            double x2 = 0.0;
            for (i = 0; i < ps.Count; i++)
            {
                x += ps[i].X;
                y += ps[i].Y;
                xy += ps[i].X * ps[i].Y;
                x2 += ps[i].X * ps[i].X;
            }
            double a = ((ps.Count * xy) - (x * y)) / ((ps.Count * x2) - (x * x));
            double b = (y - (a * x)) / ps.Count;
            PointF p1 = new PointF(0, (float)b);
            PointF p2 = new PointF(1, (float)(a + b));
            Rectangle rect = RectangleOnMorePoint(ps);
            //rect = new RectangleF(rect.X - 1, rect.Y - 1, rect.Width + 2, rect.Height + 2);
            //List<PointF> cross = CrossPointOfLineAndRect(p1, p2, rect);
            //PointA = cross[0];
            //PointB = cross[1];
            List<PointF> p4 = new List<PointF>();
            PointF ru = new PointF(rect.X + rect.Width, rect.Y);
            PointF rd = new PointF(rect.X + rect.Width, rect.Y + rect.Height);
            PointF ld = new PointF(rect.X, rect.Y + rect.Height);
            p4.Add(rect.Location);
            p4.Add(ru);
            p4.Add(rd);
            p4.Add(ld);
            PointF[] pp1 = new PointF[2];
            pp1[0] = Pvint(p4[0], p1, p2);
            pp1[1] = Pvint(p4[2], p1, p2);
            PointF[] pp2 = new PointF[2];
            pp2[0] = Pvint(p4[1], p1, p2);
            pp2[1] = Pvint(p4[3], p1, p2);
            double t1 = LineLength(pp1[0], pp1[1]);
            double t2 = LineLength(pp2[0], pp2[1]);
            if (t1 > t2)
            {
                PointA = pp1[0];
                PointB = pp1[1];
            }
            else
            {
                PointA = pp2[0];
                PointB = pp2[1];
            }
        }

        /// <summary>
        /// 求和算法
        /// </summary>
        /// <param name="pf">集合</param>
        /// <returns>返回集合的总和</returns>
        private float Sum(float[] pf)
        {
            int n = pf.Length;
            float sum = 0;
            for (int i = 0; i < n; i++)
            {
                sum += pf[i];
            }
            return sum;
        }
        /// <summary>
        /// 求集合点Xi*Yi的总和
        /// </summary>
        /// <param name="pf">点的集合</param>
        /// <returns>返回总和</returns>
        private float Multiply(List<PointF> pf)
        {
            int n = pf.Count;
            float Mul = 0;
            for (int i = 0; i < n; i++)
            {
                Mul += pf[i].X * pf[i].Y;
            }
            return Mul;
        }
        /// <summary>
        /// 求数组平方和
        /// </summary>
        /// <param name="pf">数组</param>
        /// <returns>返回总和</returns>
        private float PingFangHe(float[] pf)
        {
            int n = pf.Length;
            float sum = 0;
            for (int i = 0; i < n; i++)
            {
                sum += pf[i] * pf[i];
            }
            return sum;
        }
        /// <summary>
        /// 拟合圆算法
        /// </summary>
        /// <param name="ps">圆集合</param>
        /// <param name="center">圆心坐标</param>
        /// <param name="Rand">圆半径</param>
        public void NiHeCircle(List<PointF> ps, out PointF center, out double Rand)
        {
            int m_nNum = ps.Count;
            if (m_nNum < 3)
            {
                center = PointF.Empty;
                Rand = 0;
            }
            int i = 0;

            double X1 = 0;
            double Y1 = 0;
            double X2 = 0;
            double Y2 = 0;
            double X3 = 0;
            double Y3 = 0;
            double X1Y1 = 0;
            double X1Y2 = 0;
            double X2Y1 = 0;

            for (i = 0; i < m_nNum; i++)
            {
                X1 = X1 + ps[i].X;
                Y1 = Y1 + ps[i].Y;
                X2 = X2 + ps[i].X * ps[i].X;
                Y2 = Y2 + ps[i].Y * ps[i].Y;
                X3 = X3 + ps[i].X * ps[i].X * ps[i].X;
                Y3 = Y3 + ps[i].Y * ps[i].Y * ps[i].Y;
                X1Y1 = X1Y1 + ps[i].X * ps[i].Y;
                X1Y2 = X1Y2 + ps[i].X * ps[i].Y * ps[i].Y;
                X2Y1 = X2Y1 + ps[i].X * ps[i].X * ps[i].Y;
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

            center = new PointF(Convert.ToSingle(a / (-2)), Convert.ToSingle(b / (-2)));
            Rand = Math.Sqrt(a * a + b * b - 4 * c) / 2;
        }

        /// <summary>
        /// 拟合圆算法
        /// </summary>
        /// <param name="ps">圆集合</param>
        /// <param name="center">圆心坐标</param>
        /// <param name="Rand">圆半径</param>
        public void NiHeCircle(List<Point> ps, out PointF center, out double Rand)
        {
            int m_nNum = ps.Count;
            if (m_nNum < 3)
            {
                center = PointF.Empty;
                Rand = 0;
            }
            int i = 0;

            double X1 = 0;
            double Y1 = 0;
            double X2 = 0;
            double Y2 = 0;
            double X3 = 0;
            double Y3 = 0;
            double X1Y1 = 0;
            double X1Y2 = 0;
            double X2Y1 = 0;

            for (i = 0; i < m_nNum; i++)
            {
                X1 = X1 + ps[i].X;
                Y1 = Y1 + ps[i].Y;
                X2 = X2 + ps[i].X * ps[i].X;
                Y2 = Y2 + ps[i].Y * ps[i].Y;
                X3 = X3 + ps[i].X * ps[i].X * ps[i].X;
                Y3 = Y3 + ps[i].Y * ps[i].Y * ps[i].Y;
                X1Y1 = X1Y1 + ps[i].X * ps[i].Y;
                X1Y2 = X1Y2 + ps[i].X * ps[i].Y * ps[i].Y;
                X2Y1 = X2Y1 + ps[i].X * ps[i].X * ps[i].Y;
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

            center = new PointF(Convert.ToSingle(a / (-2)), Convert.ToSingle(b / (-2)));
            Rand = Math.Sqrt(a * a + b * b - 4 * c) / 2;
        }

        /// <summary>
        /// 拟合弧
        /// </summary>
        /// <param name="pf">弧上点的数组。其中第一点和最后一点一定要是起始点。</param>
        /// <param name="center">弧心坐标</param>
        /// <param name="rand">弧半径</param>       
        public void NiHeArc(List<PointF> pf, out PointF center, out double rand)
        {
            //m_dX,m_dY分别存储样本点坐标，TWE为点个数
            int TWE = pf.Count;
            double x0 = (double)(pf[0].X);
            double y0 = (double)(pf[0].Y);
            double x1 = (double)(pf[TWE - 1].X);
            double y1 = (double)(pf[TWE - 1].Y);

            //获取中间点集的中点
            double x2 = 0, y2 = 0;
            for (int i = 1; i < TWE - 1; i++)
            {
                x2 += (double)(pf[i].X);
                y2 += (double)(pf[i].X);
            }
            x2 /= (TWE - 2); y2 /= (TWE - 2);

            //分别求两个弦的中点和垂直平分线 
            double X1 = (x0 + x1) / 2.0;
            double Y1 = (y0 + y1) / 2.0;
            double K1 = -(x1 - x0) / (y1 - y0);

            double X2 = (x1 + x2) / 2.0;
            double Y2 = (y1 + y2) / 2.0;
            double K2 = -(x2 - x1) / (y2 - y1);

            //求取圆心和半径
            double m_X = (Y2 - Y1 + K1 * X1 - K2 * X2) / (K1 - K2);
            double m_Y = Y1 + K1 * (m_X - X1);
            double r = Math.Sqrt((x0 - m_X) * (x0 - m_X) + (y0 - m_Y) * (y0 - m_Y));


            //求圆心到大弦中点的距离
            double H = Math.Sqrt((m_X - X1) * (m_X - X1) + (m_Y - Y1) * (m_Y - Y1));
            //求误差平方和
            double tmp = 0.0;
            for (int i = 1; i < TWE - 1; i++)
            {
                double len = r - Math.Sqrt((pf[i].X - m_X) * (pf[i].X - m_X) + (pf[i].Y - m_Y) * (pf[i].Y - m_Y));
                tmp += len * len;
            }
            tmp /= (TWE - 2);

            //采用最小二乘法搜索最合理圆心位置
            double realX = m_X;
            double realY = m_Y;
            //虚拟圆心坐标 virX,virY virr vireps;
            double virX, virY, virr, vireps;

            for (int j = -100; j <= 100; j++)
            {
                virX = realX + 0.01 * H * j * 1.0 / (Math.Sqrt(K1 * K1 + 1));//虚拟圆心x坐标
                virY = realY + 0.01 * H * j * K1 / (Math.Sqrt(K1 * K1 + 1));//虚拟圆心y坐标
                virr = Math.Sqrt((x0 - virX) * (x0 - virX) + (y0 - virY) * (y0 - virY));//虚拟半径

                vireps = 0.0;//虚拟精度
                for (int i = 1; i < TWE - 1; i++)
                {
                    double len = virr - Math.Sqrt((pf[i].X - virX) * (pf[i].X - virX) + (pf[i].Y - virY) * (pf[i].Y - virY));
                    vireps += len * len;
                }
                vireps /= (TWE - 2);

                if (vireps < tmp)
                {
                    tmp = vireps;
                    m_X = virX;
                    m_Y = virY;
                    r = virr;
                }
            }
            //求相关参数
            //最长弧
            double a = Math.Sqrt((pf[0].X - pf[TWE - 1].X) * (pf[0].X - pf[TWE - 1].X) + (pf[0].Y - pf[TWE - 1].Y) * (pf[0].Y - pf[TWE - 1].Y));
            //弧心角
            double a1 = 2.0 * Math.Asin(a / 2.0 / r) * 180.0 / Math.PI;
            //弧长
            double c = r * 2.0 * Math.Asin(a / 2.0 / r);
            center = new PointF((float)m_X, (float)m_Y);
            rand = r;
        }

        /// <summary>
        /// 拟合弧
        /// </summary>
        /// <param name="pf">弧上点的数组。其中第一点和最后一点一定要是起始点。</param>
        /// <param name="center">弧心坐标</param>
        /// <param name="rand">弧半径</param>       
        public void NiHeArc(List<Point> pf, out PointF center, out double rand)
        {
            //m_dX,m_dY分别存储样本点坐标，TWE为点个数
            int TWE = pf.Count;
            double x0 = (double)(pf[0].X);
            double y0 = (double)(pf[0].Y);
            double x1 = (double)(pf[TWE - 1].X);
            double y1 = (double)(pf[TWE - 1].Y);

            //获取中间点集的中点
            double x2 = 0, y2 = 0;
            for (int i = 1; i < TWE - 1; i++)
            {
                x2 += (double)(pf[i].X);
                y2 += (double)(pf[i].X);
            }
            x2 /= (TWE - 2); y2 /= (TWE - 2);

            //分别求两个弦的中点和垂直平分线 
            double X1 = (x0 + x1) / 2.0;
            double Y1 = (y0 + y1) / 2.0;
            double K1 = -(x1 - x0) / (y1 - y0);

            double X2 = (x1 + x2) / 2.0;
            double Y2 = (y1 + y2) / 2.0;
            double K2 = -(x2 - x1) / (y2 - y1);

            //求取圆心和半径
            double m_X = (Y2 - Y1 + K1 * X1 - K2 * X2) / (K1 - K2);
            double m_Y = Y1 + K1 * (m_X - X1);
            double r = Math.Sqrt((x0 - m_X) * (x0 - m_X) + (y0 - m_Y) * (y0 - m_Y));


            //求圆心到大弦中点的距离
            double H = Math.Sqrt((m_X - X1) * (m_X - X1) + (m_Y - Y1) * (m_Y - Y1));
            //求误差平方和
            double tmp = 0.0;
            for (int i = 1; i < TWE - 1; i++)
            {
                double len = r - Math.Sqrt((pf[i].X - m_X) * (pf[i].X - m_X) + (pf[i].Y - m_Y) * (pf[i].Y - m_Y));
                tmp += len * len;
            }
            tmp /= (TWE - 2);

            //采用最小二乘法搜索最合理圆心位置
            double realX = m_X;
            double realY = m_Y;
            //虚拟圆心坐标 virX,virY virr vireps;
            double virX, virY, virr, vireps;

            for (int j = -100; j <= 100; j++)
            {
                virX = realX + 0.01 * H * j * 1.0 / (Math.Sqrt(K1 * K1 + 1));//虚拟圆心x坐标
                virY = realY + 0.01 * H * j * K1 / (Math.Sqrt(K1 * K1 + 1));//虚拟圆心y坐标
                virr = Math.Sqrt((x0 - virX) * (x0 - virX) + (y0 - virY) * (y0 - virY));//虚拟半径

                vireps = 0.0;//虚拟精度
                for (int i = 1; i < TWE - 1; i++)
                {
                    double len = virr - Math.Sqrt((pf[i].X - virX) * (pf[i].X - virX) + (pf[i].Y - virY) * (pf[i].Y - virY));
                    vireps += len * len;
                }
                vireps /= (TWE - 2);

                if (vireps < tmp)
                {
                    tmp = vireps;
                    m_X = virX;
                    m_Y = virY;
                    r = virr;
                }
            }
            //求相关参数
            //最长弧
            double a = Math.Sqrt((pf[0].X - pf[TWE - 1].X) * (pf[0].X - pf[TWE - 1].X) + (pf[0].Y - pf[TWE - 1].Y) * (pf[0].Y - pf[TWE - 1].Y));
            //弧心角
            double a1 = 2.0 * Math.Asin(a / 2.0 / r) * 180.0 / Math.PI;
            //弧长
            double c = r * 2.0 * Math.Asin(a / 2.0 / r);
            center = new PointF((float)m_X, (float)m_Y);
            rand = r;
        }

        //r=Dotmultiply(p1,p2,p0),得到向量(p1-p0)和(p2-p0)的点积，如果两个向量都非零矢量
        //r<0:两向量夹角为锐角；r=0：两向量夹角为直角；r>0:两向量夹角为钝角
        /// <summary>
        /// 求平面向量的点积，并可以根据返回值判断向量的夹角是锐是直还是钝
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p0"></param>
        /// <returns></returns>
        public static double Dotmultiply(Point2 p1, Point2 p2, Point2 p0)
        {
            return ((p1.X - p0.X) * (p2.X - p0.X) + (p1.Y - p0.Y) * (p2.Y - p0.Y));
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
        /// 判断点p是否在线段l上，条件：(p在线段l所在的直线上)而且(点p在以线段l为对角线的矩形内)
        /// </summary>
        /// <param name="l"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public static bool OnLineSEG(Line2 l, Point2 p)
        {
            return ((OnLine(l, p)) && (((p.X - l.Point1.X) * (p.X - l.Point2.X) <= 0) && ((p.Y - l.Point1.Y) * (p.Y - l.Point2.Y) <= 0)));
        }
        /// <summary>
        /// 判断点p是否在直线l上
        /// </summary>
        /// <param name="l"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public static bool OnLine(Line2 l, Point2 p)
        {
            return Multiply(l.Point2, p, l.Point1) == 0;
        }

        /// <summary>
        /// 返回点p以点o为圆心逆时针旋转alpha(单位：弧度)后所在的位置
        /// </summary>
        /// <param name="o"></param>
        /// <param name="alpha"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public static Point2 Rotate(Point2 o, double alpha, Point2 p)
        {
            Point2 tp = new Point2();
            p.X -= o.X;
            p.Y -= o.Y;
            tp.X = p.X * (float)Math.Cos(alpha) - p.Y * (float)Math.Sin(alpha) + o.X;
            tp.Y = p.Y * (float)Math.Cos(alpha) + p.X * (float)Math.Sin(alpha) + o.Y;
            return tp;
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
        /// 求点C到线段AB所在直线的垂足 P
        /// </summary>
        /// <param name="p"></param>
        /// <param name="l"></param>
        /// <returns></returns> 
        public static Point2 Perpendicular(Point2 p, Line2 l)
        {
            float r = (float)Relation(p, l);
            Point2 tp = new Point2(l.Point1.X + r * (l.Point2.X - l.Point1.X), l.Point1.Y + r * (l.Point2.Y - l.Point1.Y));
            return tp;
        }

        /// <summary>
        /// 求点p到线段l的最短距离,并返回线段上距该点最近的点np
        /// 注意：np是线段l上到点p最近的点，不一定是垂足
        /// </summary>
        /// <param name="p"></param>
        /// <param name="l"></param>
        /// <param name="np"></param>
        /// <returns></returns>
        public static double PtoLinesegDist(Point2 p, Line2 l, out Point2 np)
        {
            double r = Relation(p, l);
            if (r < 0)
            {
                np = l.Point1;
                return dist(p, l.Point1);
            }
            if (r > 1)
            {
                np = l.Point2;
                return dist(p, l.Point2);
            }
            np = Perpendicular(p, l);
            return dist(p, np);
        }

        /// <summary>
        /// 点与线段的关系
        /// </summary>
        /// <param name="p"></param>
        /// <param name="l"></param>
        /// <returns></returns>
        public static double Relation(Point2 p, Line2 l)
        {
            Line2 tl = new Line2(l.Point1, p);
            return Dotmultiply(tl.Point2, l.Point2, l.Point1) / (dist(l.Point1, l.Point2) * dist(l.Point1, l.Point2));
        }

        /// <summary>
        /// 求点p到线段l所在直线的距离(即点p到垂足的距离),请注意本函数与PtoLinesegDist函数的区别
        /// </summary>
        /// <param name="p"></param>
        /// <param name="l"></param>
        /// <returns></returns>  
        public static double PtoLdist(Point2 p, Line2 l)
        {
            return Math.Abs(Multiply(p, l.Point2, l.Point1)) / dist(l.Point1, l.Point2);
        }

        /// <summary>
        /// 计算点到折线集的最近距离,并返回最近点.
        ///注意：调用的是ptolineseg()函数
        /// </summary>
        /// <param name="pointset"></param>
        /// <param name="p"></param>
        /// <param name="q"></param>
        /// <returns></returns>
        public static double PtoPointSet(Point2[] pointset, Point2 p, out Point2 q)
        {
            double cd = double.PositiveInfinity, td;
            Line2 l;
            Point2 tq, cq = null;
            int vcount = pointset.Length;
            for (int i = 0; i < vcount - 1; i++)
            {
                l = new Line2(pointset[i], pointset[i + 1]);
                td = PtoLinesegDist(p, l, out tq);
                if (td < cd)
                {
                    cd = td;
                    cq = tq;
                }
            }
            q = cq;
            return cd;
        }

        /// <summary>
        /// 判断圆是否在多边形内
        /// ptolineseg()函数的应用2
        /// </summary>
        /// <param name="center"></param>
        /// <param name="radius"></param>
        /// <param name="polygon"></param>
        /// <returns></returns>
        public static bool CircleInsidePolygon(Point2 center, double radius, Point2[] polygon)
        {
            Point2 q = Point2.ZERO;
            double d = PtoPointSet(polygon, center, out q);
            if (d < radius || Math.Abs(d - radius) < Common.Eps)
                return true;
            else
                return false;
        }

        /// <summary>
        /// 返回两个矢量l1和l2的夹角的余弦(-1 --- 1)注意：如果想从余弦求夹角的话，注意反余弦函数的定义域是从 0到pi
        /// </summary>
        /// <param name="l1"></param>
        /// <param name="l2"></param>
        /// <returns></returns>
        public static double cosine(Line2 l1, Line2 l2)
        {
            return (((l1.Point2.X - l1.Point1.X) * (l2.Point2.X - l2.Point1.X) + (l1.Point2.Y - l1.Point1.Y) * (l2.Point2.Y - l2.Point1.Y)) / (dist(l1.Point2, l1.Point1) * dist(l2.Point2, l2.Point1)));
        }
        /// <summary>
        /// 四点定角，p2为顶点
        /// </summary>
        /// <param name="ps"></param>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        public static void GetAngle(Point3[] ps, out Point3 p1, out Point3 p2, out Point3 p3)
        {
            p1 = null;
            p2 = null;
            p3 = null;
            Point3[] pps = new Point3[4];
            if (ps.Length < 4)
                return;
            for (int i = 0; i < 4; i++)
            {
                pps[i] = new Point3(ps[i].X, ps[i].Y, ps[i].Z);
            }
            Vector3 v1 = new Vector3(pps[0], pps[1], false);
            Vector3 v2 = new Vector3(pps[2], pps[3], false);
            //p2 = new Point3(pps[1].X, pps[1].Y, pps[1].Z);//初始化为第二个点
            Point3 cross = Vector3.GetCross(v1, v2);
            p2 = cross;
            double max = 0;
            int maxIndex = 0;
            if (cross == null)
                return;
            for (int i = 0; i < 2; i++)
            {
                if (pps[i] != null)
                {
                    if (max < PPdistance.GetPPdis3(cross.X, cross.Y, cross.Z, pps[i].X, pps[i].Y, pps[i].Z))
                    {
                        max = PPdistance.GetPPdis3(cross.X, cross.Y, cross.Z, pps[i].X, pps[i].Y, pps[i].Z);
                        maxIndex = i;
                    }
                }
            }
            max = 0;
            p1 = new Point3(pps[maxIndex].X, pps[maxIndex].Y, pps[maxIndex].Z);
            for (int i = 2; i < 4; i++)
            {
                if (max < PPdistance.GetPPdis3(cross.X, cross.Y, cross.Z, pps[i].X, pps[i].Y, pps[i].Z))
                {
                    max = PPdistance.GetPPdis3(cross.X, cross.Y, cross.Z, pps[i].X, pps[i].Y, pps[i].Z);
                    maxIndex = i;
                }
            }
            p3 = new Point3(pps[maxIndex].X, pps[maxIndex].Y, pps[maxIndex].Z);
        }

        /* 返回顶角在o点，起始边为os，终止边为oe的夹角(单位：弧度)
角度小于pi，返回正值
角度大于pi，返回负值
可以用于求平面中角、线段、平面之间的夹角
*/
        /// <summary>
        /// 返回顶角在o点，起始边为os，终止边为oe的夹角(单位：弧度)
        /// </summary>
        /// <param name="o"></param>
        /// <param name="s"></param>
        /// <param name="e"></param>
        /// <returns></returns>
        public static double GetAngleValue(Point2 o, Point2 s, Point2 e)
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
        /// 返回线段l1与l2之间的夹角 单位：弧度 范围(-pi，pi)
        /// </summary>
        /// <param name="l1"></param>
        /// <param name="l2"></param>
        /// <returns></returns>
        public static double lsangle(Line2 l1, Line2 l2)
        {
            Point2 o = new Point2(), s = new Point2(), e = new Point2();
            o.X = o.Y = 0;
            s.X = l1.Point2.X - l1.Point1.X;
            s.Y = l1.Point2.Y - l1.Point1.Y;
            e.X = l2.Point2.X - l2.Point1.X;
            e.Y = l2.Point2.Y - l2.Point1.Y;
            return GetAngleValue(o, s, e);
        }

        /// <summary>
        /// 如果线段u和v相交(包括相交在端点处)时，返回true
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns> 
        public static bool intersect(Line2 u, Line2 v)
        {
            return ((Math.Max(u.Point1.X, u.Point2.X) >= Math.Min(v.Point1.X, v.Point2.X)) &&                     //排斥实验
                    (Math.Max(v.Point1.X, v.Point2.X) >= Math.Min(u.Point1.X, u.Point2.X)) &&
                    (Math.Max(u.Point1.Y, u.Point2.Y) >= Math.Min(v.Point1.Y, v.Point2.Y)) &&
                    (Math.Max(v.Point1.Y, v.Point2.Y) >= Math.Min(u.Point1.Y, u.Point2.Y)) &&
                    (Multiply(v.Point1, u.Point2, u.Point1) * Multiply(u.Point2, v.Point2, u.Point1) >= 0) &&         //跨立实验

                    (Multiply(u.Point1, v.Point2, v.Point1) * Multiply(v.Point2, u.Point2, v.Point1) >= 0));
        }

        /// <summary>
        /// (线段u和v相交)而且(交点不是双方的端点) 时返回true
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        public static bool intersect_A(Line2 u, Line2 v)
        {
            return ((intersect(u, v)) &&
                   (!OnLine(u, v.Point1)) &&
                   (!OnLine(u, v.Point2)) &&
                   (!OnLine(v, u.Point2)) &&
                   (!OnLine(v, u.Point1)));
        }


        /// <summary>
        /// 线段v所在直线与线段u相交时返回true；方法：判断线段u是否跨立线段v
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        public static bool intersect_l(Line2 u, Line2 v)
        {
            return Multiply(u.Point1, v.Point2, v.Point1) * Multiply(v.Point2, u.Point2, v.Point1) >= 0;
        }

        /// <summary>
        /// 根据已知两点坐标，求过这两点的直线解析方程： a*x+b*y+c = 0  (a >= 0)
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static Line2 makeline(Point2 p1, Point2 p2)
        {
            Line2 tl = new Line2();
            int sign = 1;
            float a = p2.Y - p1.Y;
            if (tl.A < 0)
            {
                sign = -1;
                a = sign * a;
            }
            float b = sign * (p1.X - p2.X);
            float c = sign * (p1.Y * p2.X - p1.X * p2.Y);
            tl = new Line2(a, b, c, 1);
            return tl;
        }

        /// <summary>
        /// 求点p关于直线l的对称点
        /// </summary>
        /// <param name="l"></param>
        /// <param name="p"></param>
        /// <returns></returns>  
        public static Point2 Symmetry(Line2 l, Point2 p)
        {
            float x = ((l.B * l.B - l.A * l.A) * p.X - 2 * l.A * l.B * p.Y - 2 * l.A * l.C) / (l.A * l.A + l.B * l.B);
            float y = ((l.A * l.A - l.B * l.B) * p.Y - 2 * l.A * l.B * p.X - 2 * l.B * l.C) / (l.A * l.A + l.B * l.B);
            Point2 tp = new Point2(x, y);
            return tp;
        }

        /// <summary>
        /// 如果两条直线 l1(a1*x+b1*y+c1 = 0), l2(a2*x+b2*y+c2 = 0)相交，返回true，且返回交点p 
        /// </summary>
        /// <param name="l1"></param>
        /// <param name="l2"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public static bool LineIntersect(Line2 l1, Line2 l2, out Point2 p) // 是 L1，L2
        {
            p = null;
            double d = l1.A * l2.B - l2.A * l1.B;
            if (Math.Abs(d) < Common.Eps) // 不相交
                return false;
            float x = (l2.C * l1.B - l1.C * l2.B) / (float)d;
            float y = (l2.A * l1.C - l1.A * l2.C) / (float)d;
            p = new Point2(x, y);
            return true;
        }

        /// <summary>
        /// 如果线段l1和l2相交，返回true，且返回交点inter
        /// </summary>
        /// <param name="l1"></param>
        /// <param name="l2"></param>
        /// <param name="inter"></param>
        /// <returns></returns>
        public static bool LineSEGIntersection(Line2 l1, Line2 l2, out Point2 inter)
        {
            Line2 ll1, ll2;
            ll1 = makeline(l1.Point1, l1.Point2);
            ll2 = makeline(l2.Point1, l2.Point2);
            if (LineIntersect(ll1, ll2, out inter))
            {
                return OnLine(l1, inter);
            }
            else
                return false;
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
        public static bool IsSimple(Point2[] polygon)
        {
            int cn;
            Line2 l1 = new Line2(), l2 = new Line2();
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
        /// 返回值：按输入顺序返回多边形顶点的凸凹性判断，bc[i]=true,iff:第i个顶点是凸顶点
        /// </summary>
        /// <param name="polygon"></param>
        /// <param name="bc"></param>
        public static void CheckConvex(Point2[] polygon, out bool[] bc)
        {
            int index = 0;
            Point2 tp = polygon[0];
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
                if (Multiply(polygon[(index + 1) % vcount], polygon[(index + 2) % vcount], polygon[index]) >= 0)
                    bc[(index + 1) % vcount] = true;
                else
                    bc[(index + 1) % vcount] = false;
                index++;
                count--;
            }
        }

        /// <summary>
        /// 返回值：多边形polygon是凸多边形时，返回true 
        /// </summary>
        /// <param name="polygon"></param>
        /// <returns></returns>
        public static bool IsConvex(Point2[] polygon)
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
        public static double Area_of_polygon(Point2[] polygon)
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
        public static bool IsConterClock(Point2[] polygon)
        {
            return Area_of_polygon(polygon) > 0;
        }

        /// <summary>
        /// 另一种判断多边形顶点排列方向的方法，逆时针为true
        /// </summary>
        /// <param name="polygon"></param>
        /// <returns></returns>
        public static bool IsCCwize(Point2[] polygon)
        {
            int index;
            Point2 a, b, v;
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
            return Multiply(v, b, a) > 0;//逆时针为正
        }

        /// <summary>
        /// 射线法判断点q与多边形polygon的位置关系
        /// 要求polygon为简单多边形，顶点逆时针排列
        ///如果点在多边形内：   返回0
        ///如果点在多边形边上： 返回1
        ///如果点在多边形外： 返回2
        /// </summary>
        /// <param name="Polygon"></param>
        /// <param name="q"></param>
        /// <returns></returns>
        public static int InsidePolygon(Point2[] Polygon, Point2 q)
        {
            int c = 0, i, n;
            Line2 l1 = new Line2(), l2 = new Line2();
            int vcount = Polygon.Length;
            bool bintersect_a, bonline1, bonline2, bonline3;
            double r1, r2;

            l1.Point1 = q;
            l1.Point2 = q;
            l1.Point2.X = float.PositiveInfinity;
            n = vcount;
            for (i = 0; i < vcount; i++)
            {
                l2.Point1 = Polygon[i];
                l2.Point2 = Polygon[(i + 1) % n];
                if (OnLine(l2, q)) return 1; // 如果点在边上，返回1
                if ((bintersect_a = intersect_A(l1, l2)) || // 相交且不在端点
                (
                (bonline1 = OnLine(l1, Polygon[(i + 1) % n])) && // 第二个端点在射线上
                (
                (!(bonline2 = OnLine(l1, Polygon[(i + 2) % n]))) && /* 前一个端点和后一个
端点在射线两侧 */
                ((r1 = Multiply(Polygon[i], Polygon[(i + 1) % n], l1.Point1) * Multiply(Polygon[(i + 1) % n], Polygon[(i + 2) % n], l1.Point1)) > 0) ||
                (bonline3 = OnLine(l1, Polygon[(i + 2) % n])) && /* 下一条边是水平线，
前一个端点和后一个端点在射线两侧  */
                      ((r2 = Multiply(Polygon[i], Polygon[(i + 2) % n], l1.Point1) * Multiply(Polygon[(i + 2) % n], Polygon[(i + 3) % n], l1.Point1)) > 0)
                        )
                      )
                   ) c++;
            }
            if (c % 2 == 1)
                return 0;
            else
                return 2;
        }

        /// <summary>
        /// 点q是凸多边形polygon内时，返回true；注意：多边形polygon一定要是凸多边形
        /// </summary>
        /// <param name="polygon"></param>
        /// <param name="q"></param>
        /// <returns></returns>
        public static bool InsideConvexPolygon(Point2[] polygon, Point2 q) // 可用于三角形！
        {
            Point2 p = Point2.ZERO;
            Line2 l = new Line2();
            int vcount = polygon.Length;
            int i;
            for (i = 0; i < vcount; i++) // 寻找一个肯定在多边形polygon内的点p：多边形顶点平均值
            {
                p.X += polygon[i].X;
                p.Y += polygon[i].Y;
            }
            p.X /= vcount;
            p.Y /= vcount;

            for (i = 0; i < vcount; i++)
            {
                l.Point1 = polygon[i];
                l.Point2 = polygon[(i + 1) % vcount];
                if (Multiply(p, l.Point2, l.Point1) * Multiply(q, l.Point2, l.Point1) < 0) /* 点p和点q在边l的两侧，说明
点q肯定在多边形外 */
                    break;
            }
            return (i == vcount);
        }

        /* 判断线段是否在简单多边形内(注意：如果多边形是凸多边形，下面的算法可以化简)

原理：
必要条件一：线段的两个端点都在多边形内；
必要条件二：线段和多边形的所有边都不内交；
用途：1. 判断折线是否在简单多边形内
  2. 判断简单多边形是否在另一个简单多边形内
*/
        /// <summary>
        /// 判断线段是否在简单多边形内
        /// </summary>
        /// <param name="polygon"></param>
        /// <param name="l"></param>
        /// <returns></returns>
        public static bool LinesegInsidePolygon(Point2[] polygon, Line2 l)
        {
            // 判断线段l的端点是否不都在多边形内
            if (InsidePolygon(polygon, l.Point1) != 0 || InsidePolygon(polygon, l.Point2) != 0)
                return false;
            int top = 0, i, j;
            int vcount = polygon.Length;
            Point2[] PointSet = new Point2[vcount];
            Point2 tmp = new Point2();
            Line2 s = new Line2();

            for (i = 0; i < vcount; i++)
            {
                s.Point1 = polygon[i];
                s.Point2 = polygon[(i + 1) % vcount];
                if (OnLine(s, l.Point1)) //线段l的起始端点在线段s上
                    PointSet[top++] = l.Point1;
                else if (OnLine(s, l.Point2)) //线段l的终止端点在线段s上
                    PointSet[top++] = l.Point2;
                else
                {
                    if (OnLine(l, s.Point1)) //线段s的起始端点在线段l上
                        PointSet[top++] = s.Point1;
                    else if (OnLine(l, s.Point2)) // 线段s的终止端点在线段l上
                        PointSet[top++] = s.Point2;
                    else
                    {
                        if (intersect(l, s)) // 这个时候如果相交，肯定是内交，返回false
                            return false;
                    }
                }
            }

            for (i = 0; i < top - 1; i++) /* 冒泡排序，x坐标小的排在前面；x坐标相同者，
y坐标小的排在前面 */
            {
                for (j = i + 1; j < top; j++)
                {
                    if (PointSet[i].X > PointSet[j].X || Math.Abs(PointSet[i].X - PointSet[j].X) < Common.Eps && PointSet[i].Y > PointSet[j].Y)
                    {
                        tmp = PointSet[i];
                        PointSet[i] = PointSet[j];
                        PointSet[j] = tmp;
                    }
                }
            }

            for (i = 0; i < top - 1; i++)
            {
                tmp.X = (PointSet[i].X + PointSet[i + 1].X) / 2; //得到两个相邻交点的中点
                tmp.Y = (PointSet[i].Y + PointSet[i + 1].Y) / 2;
                if (InsidePolygon(polygon, tmp) != 0)
                    return false;
            }
            return true;
        }

        /* 求任意简单多边形polygon的重心
需要调用下面几个函数：
void AddPosPart(); 增加右边区域的面积
void AddNegPart(); 增加左边区域的面积
void AddRegion(); 增加区域面积
在使用该程序时，如果把xtr,ytr,wtr,xtl,ytl,wtl设成全局变量就可以使这些函数的形式
得到化简,但要注意函数的声明和调用要做相应变化
*/
        /// <summary>
        /// 增加右边区域的面积
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="w"></param>
        /// <param name="xtr"></param>
        /// <param name="ytr"></param>
        /// <param name="wtr"></param>
        public static void AddPosPart(double x, double y, double w, ref double xtr, ref double ytr, ref double wtr)
        {
            if (Math.Abs(wtr + w) < Common.Eps) return; // detect zero regions
            xtr = (wtr * xtr + w * x) / (wtr + w);
            ytr = (wtr * ytr + w * y) / (wtr + w);
            wtr = w + wtr;
            return;
        }
        /// <summary>
        /// 增加左边区域的面积
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="w"></param>
        /// <param name="xtl"></param>
        /// <param name="ytl"></param>
        /// <param name="wtl"></param>
        public static void AddNegPart(double x, double y, double w, ref double xtl, ref double ytl, ref double wtl)
        {
            if (Math.Abs(wtl + w) < Common.Eps)
                return; // detect zero regions

            xtl = (wtl * xtl + w * x) / (wtl + w);
            ytl = (wtl * ytl + w * y) / (wtl + w);
            wtl = w + wtl;
            return;
        }
        /// <summary>
        /// 增加区域面积
        /// </summary>
        /// <param name="x1"></param>
        /// <param name="y1"></param>
        /// <param name="x2"></param>
        /// <param name="y2"></param>
        /// <param name="xtr"></param>
        /// <param name="ytr"></param>
        /// <param name="wtr"></param>
        /// <param name="xtl"></param>
        /// <param name="ytl"></param>
        /// <param name="wtl"></param>
        public static void AddRegion(double x1, double y1, double x2, double y2, ref double xtr, ref double ytr, ref double wtr, ref double xtl, ref double ytl, ref double wtl)
        {
            if (Math.Abs(x1 - x2) < Common.Eps)
                return;

            if (x2 > x1)
            {
                AddPosPart((x2 + x1) / 2, y1 / 2, (x2 - x1) * y1, ref xtr, ref ytr, ref wtr); /* rectangle 全局
变量变化处 */
                AddPosPart((x1 + x2 + x2) / 3, (y1 + y1 + y2) / 3, (x2 - x1) * (y2 - y1) / 2, ref xtr, ref ytr, ref wtr);
                // triangle 全局变量变化处
            }
            else
            {
                AddNegPart((x2 + x1) / 2, y1 / 2, (x2 - x1) * y1, ref xtl, ref ytl, ref wtl);
                // rectangle 全局变量变化处
                AddNegPart((x1 + x2 + x2) / 3, (y1 + y1 + y2) / 3, (x2 - x1) * (y2 - y1) / 2, ref xtl, ref ytl, ref wtl);
                // triangle  全局变量变化处
            }
        }

        /// <summary>
        /// 求多边形的重心
        /// </summary>
        /// <param name="polygon"></param>
        /// <returns></returns>
        public static Point2 cg_simple(Point2[] polygon)
        {
            int vcount = polygon.Length;
            double xtr, ytr, wtr, xtl, ytl, wtl;
            //注意： 如果把xtr,ytr,wtr,xtl,ytl,wtl改成全局变量后这里要删去
            Point2 p1, p2, tp = new Point2();
            xtr = ytr = wtr = 0.0;
            xtl = ytl = wtl = 0.0;
            for (int i = 0; i < vcount; i++)
            {
                p1 = polygon[i];
                p2 = polygon[(i + 1) % vcount];
                AddRegion(p1.X, p1.Y, p2.X, p2.Y, ref xtr, ref ytr, ref wtr, ref xtl, ref ytl, ref wtl); //全局变量变化处
            }
            tp.X = (float)((wtr * xtr + wtl * xtl) / (wtr + wtl));
            tp.Y = (float)((wtr * ytr + wtl * ytl) / (wtr + wtl));
            return tp;
        }

        /// <summary>
        /// 求凸多边形的重心,要求输入多边形按逆时针排序
        /// </summary>
        /// <param name="polygon"></param>
        /// <returns></returns>
        public static Point2 GravityCenter(Point2[] polygon)
        {
            int vcount = polygon.Length;
            Point2 tp = new Point2();
            double x, y, s, x0, y0, cs, k;
            x = 0; y = 0; s = 0;
            for (int i = 1; i < vcount - 1; i++)
            {
                x0 = (polygon[0].X + polygon[i].X + polygon[i + 1].X) / 3;
                y0 = (polygon[0].Y + polygon[i].Y + polygon[i + 1].Y) / 3; //求当前三角形的重心
                cs = Multiply(polygon[i], polygon[i + 1], polygon[0]) / 2;
                //三角形面积可以直接利用该公式求解
                if (Math.Abs(s) < Common.Eps)
                {
                    x = x0; y = y0; s += cs; continue;
                }
                k = cs / s; //求面积比例
                x = (x + k * x0) / (1 + k);
                y = (y + k * y0) / (1 + k);
                s += cs;
            }
            tp.X = (float)x;
            tp.Y = (float)y;
            return tp;
        }

        /* 给定一简单多边形，找出一个肯定在该多边形内的点
        定理1：每个多边形至少有一个凸顶点
        定理2：顶点数>=4的简单多边形至少有一条对角线
        结论： x坐标最大，最小的点肯定是凸顶点
               y坐标最大，最小的点肯定是凸顶点           
        */
        /// <summary>
        /// 给定一简单多边形，找出一个肯定在该多边形内的点
        /// </summary>
        /// <param name="polygon"></param>
        /// <returns></returns>
        public static Point2 A_point_insidepoly(Point2[] polygon)
        {
            int vcount = polygon.Length;
            Point2 v, a, b, r = new Point2();
            int i, index;
            v = polygon[0];
            index = 0;
            for (i = 1; i < vcount; i++) //寻找一个凸顶点
            {
                if (polygon[i].Y < v.Y)
                {
                    v = polygon[i];
                    index = i;
                }
            }
            a = polygon[(index - 1 + vcount) % vcount]; //得到v的前一个顶点
            b = polygon[(index + 1) % vcount]; //得到v的后一个顶点
            Point2[] tri = new Point2[3];
            Point2 q = new Point2();
            tri[0] = a; tri[1] = v; tri[2] = b;
            double md = double.PositiveInfinity;
            int in1 = index;
            bool bin = false;
            for (i = 0; i < vcount; i++) //寻找在三角形avb内且离顶点v最近的顶点q
            {
                if (i == index)
                    continue;
                if (i == (index - 1 + vcount) % vcount)
                    continue;
                if (i == (index + 1) % vcount)
                    continue;
                if (!InsideConvexPolygon(tri, polygon[i]))
                    continue;
                bin = true;
                if (dist(v, polygon[i]) < md)
                {
                    q = polygon[i];
                    md = dist(v, q);
                }
            }
            if (!bin) //没有顶点在三角形avb内，返回线段ab中点
            {
                r.X = (a.X + b.X) / 2;
                r.Y = (a.Y + b.Y) / 2;
                return r;
            }
            r.X = (v.X + q.X) / 2; //返回线段vq的中点
            r.Y = (v.Y + q.Y) / 2;
            return r;
        }

        /* 求从多边形外一点p出发到一个简单多边形的切线,如果存在返回切点,其中rp点是右切
        点,lp是左切点
        注意：p点一定要在多边形外
        输入顶点序列是逆时针排列
        原 理:如果点在多边形内肯定无切线;凸多边形有唯一的两个切点,凹多边形就可能有多于
        两个的切点;
        如果polygon是凸多边形，切点只有两个只要找到就可以,可以化简此算法
        如果是凹多边形还有一种算法可以求解:先求凹多边形的凸包,然后求凸包的切线
        */
        /// <summary>
        /// 求从多边形外一点p出发到一个简单多边形的切线,如果存在返回切点,
        /// 其中rp点是右切点,lp是左切点
        /// </summary>
        /// <param name="polygon"></param>
        /// <param name="p"></param>
        /// <param name="rp"></param>
        /// <param name="lp"></param>
        public static void PointTangentPoly(Point2[] polygon, Point2 p, out Point2 rp, out Point2 lp)
        {
            int vcount = polygon.Length;
            Line2 ep = new Line2(), en = new Line2();
            bool blp, bln;
            rp = polygon[0];
            lp = polygon[0];
            for (int i = 1; i < vcount; i++)
            {
                ep.Point1 = polygon[(i + vcount - 1) % vcount];
                ep.Point2 = polygon[i];
                en.Point1 = polygon[i];
                en.Point2 = polygon[(i + 1) % vcount];
                blp = Multiply(ep.Point2, p, ep.Point1) >= 0;                // p is to the left of pre edge

                bln = Multiply(en.Point2, p, en.Point1) >= 0;                // p is to the left of next edge

                if (!blp && bln)
                {
                    if (Multiply(polygon[i], rp, p) > 0)           // polygon[i] is above rp
                        rp = polygon[i];
                }
                if (blp && !bln)
                {
                    if (Multiply(lp, polygon[i], p) > 0)           // polygon[i] is below lp
                        lp = polygon[i];
                }
            }
        }

        /// <summary>
        /// 如果多边形polygon的核存在，返回true，返回核上的一点p.顶点按逆时针方向输入
        /// </summary>
        /// <param name="polygon"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public static bool Core_exist(Point2[] polygon, out Point2 p)
        {
            p = null;
            int vcount = polygon.Length;
            int i, j, k;
            Line2 l = new Line2();
            Line2[] lineset = new Line2[vcount];
            for (i = 0; i < vcount; i++)
            {
                lineset[i] = makeline(polygon[i], polygon[(i + 1) % vcount]);
            }
            for (i = 0; i < vcount; i++)
            {
                for (j = 0; j < vcount; j++)
                {
                    if (i == j) continue;
                    if (LineIntersect(lineset[i], lineset[j], out p))
                    {
                        for (k = 0; k < vcount; k++)
                        {
                            l.Point1 = polygon[k];
                            l.Point2 = polygon[(k + 1) % vcount];
                            if (Multiply(p, l.Point2, l.Point1) > 0)
                                //多边形顶点按逆时针方向排列，核肯定在每条边的左侧或边上
                                break;
                        }
                        if (k == vcount)             //找到了一个核上的点
                            break;
                    }
                }
                if (j < vcount)
                    break;
            }
            if (i < vcount)
                return true;
            else
                return false;
        }

        /// <summary>
        /// 返回值： 点p在圆内(包括边界)时，返回true
        ///用途： 因为圆为凸集，所以判断点集，折线，多边形是否在圆内时，
        ///只需要逐一判断点是否在圆内即可。
        /// </summary>
        /// <param name="o"></param>
        /// <param name="r"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public static bool Point_in_circle(Point2 o, double r, Point2 p)
        {
            double d2 = (p.X - o.X) * (p.X - o.X) + (p.Y - o.Y) * (p.Y - o.Y);
            double r2 = r * r;
            return d2 < r2 || Math.Abs(d2 - r2) < Common.Eps;
        }


        /// <summary>
        /// 用 途：求不共线的三点确定一个圆
        /// 返回值：如果三点共线，返回false；反之，返回true。圆心由q返回，半径由r返回
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <param name="q"></param>
        /// <param name="r"></param>
        /// <returns></returns>
        public static bool Cocircle(Point2 p1, Point2 p2, Point2 p3, out Point2 q, out double r)
        {
            q = null;
            r = 0;
            float x12 = p2.X - p1.X;
            float y12 = p2.Y - p1.Y;
            float x13 = p3.X - p1.X;
            float y13 = p3.Y - p1.Y;
            float z2 = x12 * (p1.X + p2.X) + y12 * (p1.Y + p2.Y);
            float z3 = x13 * (p1.X + p3.X) + y13 * (p1.Y + p3.Y);
            float d = 2.0F * (x12 * (p3.Y - p2.Y) - y12 * (p3.X - p2.X));
            if (Math.Abs(d) < Common.Eps) //共线，圆不存在
                return false;
            q = new Point2();
            q.X = (y13 * z2 - y12 * z3) / d;
            q.Y = (x12 * z3 - x13 * z2) / d;
            r = dist(p1, q);
            return true;
        }

        /*
说明：因为矩形的特殊性，常用算法可以化简：
1.判断矩形是否包含点
只要判断该点的横坐标和纵坐标是否夹在矩形的左右边和上下边之间。
2.判断线段、折线、多边形是否在矩形中
因为矩形是个凸集，所以只要判断所有端点是否都在矩形中就可以了。
3.判断圆是否在矩形中
圆在矩形中的充要条件是：圆心在矩形中且圆的半径小于等于圆心到矩形四边的距离的最
小值。
*/
        /// <summary>
        /// 已知矩形的三个顶点(a,b,c)，计算第四个顶点d的坐标. 注意：已知的三个顶点可以是无序的
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <returns></returns>
        public static Point2 Rect4thP(Point2 a, Point2 b, Point2 c)
        {
            Point2 d = new Point2();
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
            }
            return d;
        }

        //两圆关系：

        /* 两圆：
        相离： return 1；
        外切： return 2；
        相交： return 3；
        内切： return 4；
        内含： return 5；
        */
        /// <summary>
        /// 两圆关系
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="r1"></param>
        /// <param name="p2"></param>
        /// <param name="r2"></param>
        /// <returns></returns>
        public static int CircleRelation(Point2 p1, double r1, Point2 p2, double r2)
        {
            double d = Math.Sqrt((p1.X - p2.X) * (p1.X - p2.X) + (p1.Y - p2.Y) * (p1.Y - p2.Y));


            if (Math.Abs(d - r1 - r2) < Common.Eps) // 必须保证前两个if先被判定！
                return 2;
            if (Math.Abs(d - Math.Abs(r1 - r2)) < Common.Eps)
                return 4;
            if (d > r1 + r2)
                return 1;
            if (d < Math.Abs(r1 - r2))
                return 5;
            if (Math.Abs(r1 - r2) < d && d < r1 + r2)
                return 3;
            return 0; // indicate an error!
        }


        //判断圆是否在矩形内：

        /// <summary>
        /// 判定圆是否在矩形内，是就返回true（设矩形水平，且其四个顶点由左上开始按顺时针排列）
        /// </summary>
        /// <param name="pc"></param>
        /// <param name="r"></param>
        /// <param name="pr1"></param>
        /// <param name="pr2"></param>
        /// <param name="pr3"></param>
        /// <param name="pr4"></param>
        /// <returns></returns>
        public static bool CircleRecRelation(Point2 pc, double r, Point2 pr1, Point2 pr2, Point2 pr3, Point2 pr4)
        {
            if (pr1.X < pc.X && pc.X < pr2.X && pr3.Y < pc.Y && pc.Y < pr2.Y)
            {
                Line2 line1 = new Line2(pr1, pr2);
                Line2 line2 = new Line2(pr2, pr3);
                Line2 line3 = new Line2(pr3, pr4);
                Line2 line4 = new Line2(pr4, pr1);
                if (r < PtoLdist(pc, line1) && r < PtoLdist(pc, line2) &&
                r < PtoLdist(pc, line3) && r < PtoLdist(pc, line4))
                    return true;
            }
            return false;
        }


        //点到平面的距离：

        /// <summary>
        /// 点到平面的距离,平面用一般式表示ax+by+cz+d=0
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <param name="d"></param>
        /// <returns></returns>
        public static double P2planeDist(double x, double y, double z, double a, double b, double c,
         double d)
        {
            return Math.Abs(a * x + b * y + c * z + d) / Math.Sqrt(a * a + b * b + c * c);
        }


        //点是否在直线同侧：

        /// <summary>
        /// 两个点是否在直线同侧，是则返回true
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="line"></param>
        /// <returns></returns>
        public static bool SameSide(Point2 p1, Point2 p2, Line2 line)
        {
            return (line.A * p1.X + line.B * p1.Y + line.C) *
            (line.A * p2.X + line.B * p2.Y + line.C) > 0;
        }


        //镜面反射线：

        // 已知入射线、镜面，求反射线。
        // a1,b1,c1为镜面直线方程(a1 x + b1 y + c1 = 0 ,下同)系数; 
        //a2,b2,c2为入射光直线方程系数; 
        //a,b,c为反射光直线方程系数.
        // 光是有方向的，使用时注意：入射光向量:<-b2,a2>；反射光向量:<b,-a>.
        // 不要忘记结果中可能会有"negative zeros"
        /// <summary>
        /// 已知入射线、镜面，求反射线。
        /// </summary>
        /// <param name="a1"></param>
        /// <param name="b1"></param>
        /// <param name="c1"></param>
        /// <param name="a2"></param>
        /// <param name="b2"></param>
        /// <param name="c2"></param>
        /// <param name="mirrorReflectLine"></param>
        public static void Reflect(double a1, double b1, double c1, double a2, double b2, double c2, out Line2 mirrorReflectLine)
        {
            mirrorReflectLine = null;
            float a = 0, b = 0, c = 0;
            float n, m;
            float tpb, tpa;
            tpb = (float)(b1 * b2 + a1 * a2);
            tpa = (float)(a2 * b1 - a1 * b2);
            m = (float)((tpb * b1 + tpa * a1) / (b1 * b1 + a1 * a1));
            n = (float)((tpa * b1 - tpb * a1) / (b1 * b1 + a1 * a1));
            if (Math.Abs(a1 * b2 - a2 * b1) < 1e-20)
            {
                a = (float)a2; b = (float)b2; c = (float)c2;
                return;
            }
            double xx, yy; //(xx,yy)是入射线与镜面的交点。
            xx = (b1 * c2 - b2 * c1) / (a1 * b2 - a2 * b1);
            yy = (a2 * c1 - a1 * c2) / (a1 * b2 - a2 * b1);
            a = n;
            b = -m;
            c = m * (float)yy - (float)xx * n;
            mirrorReflectLine = new Line2(a, b, c, 1);
        }


        //矩形包含：
        /// <summary>
        /// 矩形2（C，D）是否在1（A，B）内
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="C"></param>
        /// <param name="D"></param>
        /// <returns></returns>
        public static bool r2inr1(double A, double B, double C, double D)
        {
            double X, Y, L, K, DMax;
            if (A < B)
            {
                double tmp = A;
                A = B;
                B = tmp;
            }
            if (C < D)
            {
                double tmp = C;
                C = D;
                D = tmp;
            }

            if (A > C && B > D)                 // trivial case 
                return true;
            else
                if (D >= B)
                    return false;
                else
                {
                    X = Math.Sqrt(A * A + B * B);         // outer rectangle's diagonal 
                    Y = Math.Sqrt(C * C + D * D);         // inner rectangle's diagonal 
                    if (Y < B) // check for marginal conditions
                        return true; // the inner rectangle can freely rotate inside
                    else
                        if (Y > X)
                            return false;
                        else
                        {
                            L = (B - Math.Sqrt(Y * Y - A * A)) / 2;
                            K = (A - Math.Sqrt(Y * Y - B * B)) / 2;
                            DMax = Math.Sqrt(L * L + K * K);
                            if (D >= DMax)
                                return false;
                            else
                                return true;
                        }
                }
        }


        //两圆交点：

        // 两圆已经相交（相切）
        /// <summary>
        /// 两圆交点
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="r1"></param>
        /// <param name="p2"></param>
        /// <param name="r2"></param>
        /// <param name="rp1"></param>
        /// <param name="rp2"></param>
        public static void c2point(Point2 p1, double r1, Point2 p2, double r2, out Point2 rp1, out Point2 rp2)
        {
            rp1 = new Point2();
            rp2 = new Point2();
            float a, b, r;
            a = p2.X - p1.X;
            b = p2.Y - p1.Y;
            r = (float)(a * a + b * b + r1 * r1 - r2 * r2) / 2;
            if (a == 0 && b != 0)
            {
                rp1.Y = rp2.Y = r / b;
                rp1.X = (float)Math.Sqrt(r1 * r1 - rp1.Y * rp1.Y);
                rp2.X = -rp1.X;
            }
            else if (a != 0 && b == 0)
            {
                rp1.X = rp2.X = r / a;
                rp1.Y = (float)Math.Sqrt(r1 * r1 - rp1.X * rp2.X);
                rp2.Y = -rp1.Y;
            }
            else if (a != 0 && b != 0)
            {
                double delta;
                delta = b * b * r * r - (a * a + b * b) * (r * r - r1 * r1 * a * a);
                rp1.Y = (b * r + (float)Math.Sqrt(delta)) / (a * a + b * b);
                rp2.Y = (b * r - (float)Math.Sqrt(delta)) / (a * a + b * b);
                rp1.X = (r - b * rp1.Y) / a;
                rp2.X = (r - b * rp2.Y) / a;
            }

            rp1.X += p1.X;
            rp1.Y += p1.Y;
            rp2.X += p1.X;
            rp2.Y += p1.Y;
        }

        /// <summary>
        /// 交换两点
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        public static void SwapPoints(ref Point2 p1, ref Point2 p2)
        {
            Point2 tmp = p1;
            p1 = p2;
            p2 = p1;
        }
        /// <summary>
        /// 交换两个浮点数值
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        public static void SwapDoubles(ref double p1, ref double p2)
        {
            double tmp = p1;
            p1 = p2;
            p2 = p1;
        }

        //两圆公共面积：

        // 必须保证相交
        /// <summary>
        /// 两圆公共面积
        /// 必须保证相交
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="r1"></param>
        /// <param name="p2"></param>
        /// <param name="r2"></param>
        /// <returns></returns>
        public static double c2area(Point2 p1, double r1, Point2 p2, double r2)
        {
            Point2 rp1, rp2;
            c2point(p1, r1, p2, r2, out rp1, out rp2);

            if (r1 > r2) //保证r2>r1
            {
                SwapPoints(ref p1, ref p2);
                SwapDoubles(ref r1, ref r2);
            }
            double a, b, rr;
            a = p1.X - p2.X;
            b = p1.Y - p2.Y;
            rr = Math.Sqrt(a * a + b * b);

            double dx1, dy1, dx2, dy2;
            double sita1, sita2;
            dx1 = rp1.X - p1.X;
            dy1 = rp1.Y - p1.Y;
            dx2 = rp2.X - p1.X;
            dy2 = rp2.Y - p1.Y;
            sita1 = Math.Acos((dx1 * dx2 + dy1 * dy2) / r1 / r1);

            dx1 = rp1.X - p2.X;
            dy1 = rp1.Y - p2.Y;
            dx2 = rp2.X - p2.X;
            dy2 = rp2.Y - p2.Y;
            sita2 = Math.Acos((dx1 * dx2 + dy1 * dy2) / r2 / r2);
            double s = 0;
            if (rr < r2) //相交弧为优弧
                s = r1 * r1 * (Math.PI - sita1 / 2 + Math.Sin(sita1) / 2) + r2 * r2 * (sita2 - Math.Sin(sita2)) / 2;
            else //相交弧为劣弧
                s = (r1 * r1 * (sita1 - Math.Sin(sita1)) + r2 * r2 * (sita2 - Math.Sin(sita2))) / 2;

            return s;
        }

        //圆和直线关系：

        //0----相离 1----相切 2----相交
        /// <summary>
        /// 圆和直线关系
        /// 0----相离 1----相切 2----相交
        /// </summary>
        /// <param name="p"></param>
        /// <param name="r"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <param name="rp1"></param>
        /// <param name="rp2"></param>
        /// <returns></returns>
        public static int clpoint(Point2 p, double r, double a, double b, double c, out Point2 rp1, out Point2 rp2)
        {
            rp1 = new Point2();
            rp2 = new Point2();
            int res = 0;
            c = c + a * p.X + b * p.Y;
            float tmp;
            if (a == 0 && b != 0)
            {
                tmp = (float)(-c / b);
                if (r * r < tmp * tmp)
                    res = 0;
                else if (r * r == tmp * tmp)
                {
                    res = 1;
                    rp1.Y = tmp;
                    rp1.X = 0;
                }
                else
                {
                    res = 2;
                    rp1.Y = rp2.Y = tmp;
                    rp1.X = (float)Math.Sqrt(r * r - tmp * tmp);
                    rp2.X = -rp1.X;
                }
            }
            else if (a != 0 && b == 0)
            {
                tmp = (float)(-c / a);
                if (r * r < tmp * tmp)
                    res = 0;
                else if (r * r == tmp * tmp)
                {
                    res = 1;
                    rp1.X = tmp;
                    rp1.Y = 0;
                }
                else
                {
                    res = 2;
                    rp1.X = rp2.X = tmp;
                    rp1.Y = (float)Math.Sqrt(r * r - tmp * tmp);
                    rp2.Y = -rp1.Y;
                }
            }
            else if (a != 0 && b != 0)
            {
                double delta;
                delta = b * b * c * c - (a * a + b * b) * (c * c - a * a * r * r);
                if (delta < 0)
                    res = 0;
                else if (delta == 0)
                {
                    res = 1;
                    rp1.Y = (float)(-b * c / (a * a + b * b));
                    rp1.X = (float)((-c - b * rp1.Y) / a);
                }
                else
                {
                    res = 2;
                    rp1.Y = (float)((-b * c + Math.Sqrt(delta)) / (a * a + b * b));
                    rp2.Y = (float)((-b * c - Math.Sqrt(delta)) / (a * a + b * b));
                    rp1.X = (float)((-c - b * rp1.Y) / a);
                    rp2.X = (float)((-c - b * rp2.Y) / a);
                }
            }
            rp1.X += p.X;
            rp1.Y += p.Y;
            rp2.X += p.X;
            rp2.Y += p.Y;
            return res;
        }

        /// <summary>
        /// 圆和直线关系
        /// 0----相离 1----相切 2----相交
        /// rp1, rp2为交点
        /// </summary>
        /// <param name="p"></param>
        /// <param name="r"></param>
        /// <param name="xl"></param>
        /// <param name="pp"></param>
        /// <param name="rp1"></param>
        /// <param name="rp2"></param>
        /// <returns></returns>
        public static int clpoint(Point2 p, double r, double xl, Point2 pp, out Point2 rp1, out Point2 rp2)
        {
            rp1 = new Point2();
            rp2 = new Point2();
            int res = 0;
            double a = xl, b = -1, c = pp.Y - xl * pp.X;
            c = c + a * p.X + b * p.Y;
            float tmp;
            if (a == 0 && b != 0)
            {
                tmp = (float)(-c / b);
                if (r * r < tmp * tmp)
                    res = 0;
                else if (r * r == tmp * tmp)
                {
                    res = 1;
                    rp1.Y = tmp;
                    rp1.X = 0;
                }
                else
                {
                    res = 2;
                    rp1.Y = rp2.Y = tmp;
                    rp1.X = (float)Math.Sqrt(r * r - tmp * tmp);
                    rp2.X = -rp1.X;
                }
            }
            else if (a != 0 && b == 0)
            {
                tmp = (float)(-c / a);
                if (r * r < tmp * tmp)
                    res = 0;
                else if (r * r == tmp * tmp)
                {
                    res = 1;
                    rp1.X = tmp;
                    rp1.Y = 0;
                }
                else
                {
                    res = 2;
                    rp1.X = rp2.X = tmp;
                    rp1.Y = (float)Math.Sqrt(r * r - tmp * tmp);
                    rp2.Y = -rp1.Y;
                }
            }
            else if (a != 0 && b != 0)
            {
                double delta;
                delta = b * b * c * c - (a * a + b * b) * (c * c - a * a * r * r);
                if (delta < 0)
                    res = 0;
                else if (delta == 0)
                {
                    res = 1;
                    rp1.Y = (float)(-b * c / (a * a + b * b));
                    rp1.X = (float)((-c - b * rp1.Y) / a);
                }
                else
                {
                    res = 2;
                    rp1.Y = (float)((-b * c + Math.Sqrt(delta)) / (a * a + b * b));
                    rp2.Y = (float)((-b * c - Math.Sqrt(delta)) / (a * a + b * b));
                    rp1.X = (float)((-c - b * rp1.Y) / a);
                    rp2.X = (float)((-c - b * rp2.Y) / a);
                }
            }
            rp1.X += p.X;
            rp1.Y += p.Y;
            rp2.X += p.X;
            rp2.Y += p.Y;
            return res;
        }

        //内切圆：
        /// <summary>
        /// 求内切圆
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <param name="rp"></param>
        /// <param name="r"></param>
        public static void incircle(Point2 p1, Point2 p2, Point2 p3, out Point2 rp, out double r)
        {
            rp = new Point2();
            double dx31, dy31, dx21, dy21, d31, d21, a1, b1, c1;
            dx31 = p3.X - p1.X;
            dy31 = p3.Y - p1.Y;
            dx21 = p2.X - p1.X;
            dy21 = p2.Y - p1.Y;

            d31 = Math.Sqrt(dx31 * dx31 + dy31 * dy31);
            d21 = Math.Sqrt(dx21 * dx21 + dy21 * dy21);
            a1 = dx31 * d21 - dx21 * d31;
            b1 = dy31 * d21 - dy21 * d31;
            c1 = a1 * p1.X + b1 * p1.Y;

            double dx32, dy32, dx12, dy12, d32, d12, a2, b2, c2;
            dx32 = p3.X - p2.X;
            dy32 = p3.Y - p2.Y;
            dx12 = -dx21;
            dy12 = -dy21;

            d32 = Math.Sqrt(dx32 * dx32 + dy32 * dy32);
            d12 = d21;
            a2 = dx12 * d32 - dx32 * d12;
            b2 = dy12 * d32 - dy32 * d12;
            c2 = a2 * p2.X + b2 * p2.Y;

            rp.X = (float)((c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1));
            rp.Y = (float)((c2 * a1 - c1 * a2) / (a1 * b2 - a2 * b1));
            r = Math.Abs(dy21 * rp.X - dx21 * rp.Y + dx21 * p1.Y - dy21 * p1.X) / d21;
        }


        //求切点：

        // p---圆心坐标， r---圆半径， sp---圆外一点， rp1,rp2---切点坐标
        /// <summary>
        /// 求切点
        /// p---圆心坐标， r---圆半径， sp---圆外一点， rp1,rp2---切点坐标
        /// </summary>
        /// <param name="p"></param>
        /// <param name="r"></param>
        /// <param name="sp"></param>
        /// <param name="rp1"></param>
        /// <param name="rp2"></param>
        public static void cutpoint(Point2 p, double r, Point2 sp, out Point2 rp1, out Point2 rp2)
        {
            Point2 p2 = new Point2();
            p2.X = (p.X + sp.X) / 2;
            p2.Y = (p.Y + sp.Y) / 2;

            double dx2, dy2, r2;
            dx2 = p2.X - p.X;
            dy2 = p2.Y - p.Y;
            r2 = Math.Sqrt(dx2 * dx2 + dy2 * dy2);
            c2point(p, r, p2, r2, out rp1, out rp2);
        }


        //线段的左右旋：

        /* l2在l1的左/右方向（l1为基准线）
        返回 0： 重合；
        返回 1： 右旋；
        返回 –1： 左旋；
        */
        /// <summary>
        /// l2在l1的左/右方向（l1为基准线），并且l1、l2不能为空
        ///返回 0： 重合；
        ///返回 1： 右旋；
        ///返回 –1： 左旋
        /// </summary>
        /// <param name="l1"></param>
        /// <param name="l2"></param>
        /// <returns></returns>
        public static int rotate(Line2 l1, Line2 l2)
        {
            double dx1, dx2, dy1, dy2;
            dx1 = l1.Point1.X - l1.Point2.X;
            dy1 = l1.Point1.Y - l1.Point2.Y;
            dx2 = l2.Point1.X - l2.Point2.X;
            dy2 = l2.Point1.Y - l2.Point2.Y;

            double d;
            d = dx1 * dy2 - dx2 * dy1;
            if (d == 0)
                return 0;
            else if (d > 0)
                return -1;
            else
                return 1;
        }

        /////////////////////////测试3点共线函数///////////////
        /// <summary>
        /// 测试3点共线函数
        /// </summary>
        /// <param name="x1"></param>
        /// <param name="y1"></param>
        /// <param name="x2"></param>
        /// <param name="y2"></param>
        /// <param name="x3"></param>
        /// <param name="y3"></param>
        /// <returns></returns>
        public static int Check3PInLine(float x1, float y1, float x2, float y2, float x3, float y3)
        {
            float k1, k2;
            k1 = 0;
            k2 = 0;
            if (x1 == x2)
            {
                if (x2 == x3)
                {
                    return 1;//在一条直线上。
                }
                else
                {
                    if (y1 == y2)
                    {
                        return 7;
                    }
                    else
                    {
                        if (y2 == y3)
                        {
                            return 8;
                        }
                        else
                        {
                            return 2;
                        }
                    }
                }
            }
            else
            {
                if (x2 == x3)
                {
                    if (y2 == y3)
                    {
                        return 7;
                    }
                    else
                    {
                        if (y1 == y2)
                        {
                            return 9;
                        }
                        else
                        {
                            return 3;
                        }
                    }
                }
                else
                {
                    if (y1 == y2)
                    {
                        if (y2 == y3)
                        {
                            return 1;
                        }
                        else
                        {
                            return 4;
                        }
                    }
                    else
                    {
                        if (y2 == y3)
                        {
                            return 5;
                        }
                        else
                        {
                            k1 = (y2 - y1) / (x2 - x1);
                            k2 = (y3 - y2) / (x3 - x2);
                            if (k1 == k2)
                            {
                                return 1;
                            }
                            else
                            {
                                return 6;
                            }
                        }
                    }
                }
            }
        }
    }
}
