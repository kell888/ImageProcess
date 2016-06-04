using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;

namespace KellImageProcess
{
    /// <summary>
    /// 多线程公用类
    /// </summary>
    public static class ThreadCommon
    {
        /// <summary>
        /// 第1方向是否连通
        /// </summary>
        public static bool linked1;
        /// <summary>
        /// 第2方向是否连通
        /// </summary>
        public static bool linked2;
        /// <summary>
        /// 第3方向是否连通
        /// </summary>
        public static bool linked3;
        /// <summary>
        /// 第1方向的相邻点集
        /// </summary>
        public static LinkedPointList lps1;
        /// <summary>
        /// 第2方向的相邻点集
        /// </summary>
        public static LinkedPointList lps2;
        /// <summary>
        /// 第3方向的相邻点集
        /// </summary>
        public static LinkedPointList lps3;
        /// <summary>
        /// 第1方向连通路径回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackRoadWithPath1(IAsyncResult ar)
        {
            AsyncMethodCallerLinkedWithPath caller = (AsyncMethodCallerLinkedWithPath)ar.AsyncState;
            caller.EndInvoke(out linked1, ar);
        }
        /// <summary>
        /// 第2方向连通路径回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackRoadWithPath2(IAsyncResult ar)
        {
            AsyncMethodCallerLinkedWithPath caller = (AsyncMethodCallerLinkedWithPath)ar.AsyncState;
            caller.EndInvoke(out linked2, ar);
        }
        /// <summary>
        /// 第3方向连通路径回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackRoadWithPath3(IAsyncResult ar)
        {
            AsyncMethodCallerLinkedWithPath caller = (AsyncMethodCallerLinkedWithPath)ar.AsyncState;
            caller.EndInvoke(out linked3, ar);
        }
        /// <summary>
        /// 判断指定的点是否为内存位图中的关键点
        /// </summary>
        /// <param name="pp"></param>
        /// <param name="data"></param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        public static bool CheckThePointIsKeyPoint(Point pp, BitmapData data, int blkORwht, byte threshold)
        {
            if (pp.X < 0 || pp.X > data.Width - 1 || pp.Y < 0 || pp.Y > data.Height - 1)
                return false;
            int BPP = Image.GetPixelFormatSize(data.PixelFormat) / 8;
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                p = (byte*)data.Scan0 + BPP * pp.X + stride * pp.Y;
                byte R = p[2];
                byte G = p[1];
                byte B = p[0];
                if (R != G || G != B)
                {
                    //尚未是灰度化，先计算灰度
                    // gray = 0.3*R + 0.59*G + 0.11*B
                    byte gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                    //再提取关键点
                    if (blkORwht == 0)//黑边
                    {
                        if (gray < threshold)
                        {
                            return true;
                        }
                    }
                    else
                    {
                        if (gray > threshold)
                        {
                            return true;
                        }
                    }
                }
                else
                {
                    if (blkORwht == 0)//黑边
                    {
                        if (B < threshold)
                        {
                            return true;
                        }
                    }
                    else
                    {
                        if (B > threshold)
                        {
                            return true;
                        }
                    }
                }
            }
            return false;
        }
        /// <summary>
        /// 判断指定的两点是否相临(八临域)
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static bool LinkedPoints(Point p1, Point p2)
        {
            bool linked = false;
            if ((p1.X == p2.X) && (p1.Y == p2.Y - 1 || p1.Y == p2.Y + 1) ||
                (p1.Y == p2.Y) && (p1.X == p2.X - 1 || p1.X == p2.X + 1) ||
                (p1.X == p2.X - 1 || p1.X == p2.X + 1) && (p1.Y == p2.Y - 1 || p1.Y == p2.Y + 1))
                linked = true;
            return linked;
        }
        /// <summary>
        /// 判断指定的两点是否相临(有临近系数near控制)
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="near"></param>
        /// <returns></returns>
        public static bool LinkedPoints(Point p1, Point p2, double near)
        {
            bool linked = false;
            if (Math.Sqrt((p1.X - p2.X) * (p1.X - p2.X) + (p1.Y - p2.Y) * (p1.Y - p2.Y)) <= near)
                linked = true;
            return linked;
        }
        /// <summary>
        /// 检测指定位图内是否存在连通点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="isCZ"></param>
        /// <param name="direction"></param>
        /// <param name="linked"></param>
        public static void CheckSomeRoadWithPath(RectBitmapWithPath bmp, bool isCZ, int direction, out bool linked)
        {
            linked = false;
            RectBitmapWithPath rectbmp = (RectBitmapWithPath)bmp;
            Bitmap bm = rectbmp.bmp;
            Rectangle rec = rectbmp.rect;
            int blkORwht = rectbmp.blkORwht;
            byte threshold = rectbmp.threshold;
            GraphicsPath path = rectbmp.path;
            int BPP = Image.GetPixelFormatSize(bm.PixelFormat) / 8;
            List<Point> ps = new List<Point>();
            List<Point> errPs = new List<Point>();
            unsafe
            {
                byte gray, R, G, B;
                BitmapData data = bm.LockBits(new Rectangle(0, 0, rec.Width, rec.Height), ImageLockMode.ReadOnly, bm.PixelFormat);
                byte* p1, p2, p3;
                int stride = data.Stride;
                if (isCZ)//垂直方向的连通判断
                {
                    switch (direction)
                    {
                        case 1:
                            #region dir=1
                            int i = 0;
                            int j = 0;
                            while (i >= 0 && i < rec.Width && j >= 0 && j < rec.Height)
                            {
                                #region 判断逻辑：西南方向
                                if (linked1 || linked2 || linked3)
                                {
                                    break;
                                }
                                Point pp = new Point(i + rec.X, j + rec.Y);
                                if (path.IsVisible(pp) && !errPs.Contains(pp))
                                {
                                    //g.DrawEllipse(Pens.Blue, pp.X, pp.Y, 1, 1);
                                    p1 = (byte*)data.Scan0 + BPP * i + stride * j;
                                    R = p1[2];
                                    G = p1[1];
                                    B = p1[0];
                                    if (R != G || G != B)
                                    {
                                        //尚未是灰度化，先计算灰度
                                        // gray = 0.3*R + 0.59*G + 0.11*B
                                        gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                                        //再提取关键点
                                        if (blkORwht == 0)//黑边
                                        {
                                            if (gray < threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (j == rec.Height - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i = ps[ps.Count - 1].X - 1 - rec.X;
                                                    if (i < 0)
                                                        i = 0;
                                                    j = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i = ps[ps.Count - 1].X - rec.X;
                                                        j = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i, j + 1), data, blkORwht, threshold))
                                                        {
                                                            j = j + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i + 1, j + 1), data, blkORwht, threshold))
                                                        {
                                                            i = i + 1;
                                                            j = j + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i = ps[ps.Count - 1].X - rec.X;
                                                                j = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (j == rec.Height - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i = errPs[errPs.Count - 1].X + 1 - rec.X;
                                                        j = errPs[errPs.Count - 1].Y - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                        else//白边
                                        {
                                            if (gray > threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (j == rec.Height - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i = ps[ps.Count - 1].X - 1 - rec.X;
                                                    if (i < 0)
                                                        i = 0;
                                                    j = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i = ps[ps.Count - 1].X - rec.X;
                                                        j = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i, j + 1), data, blkORwht, threshold))
                                                        {
                                                            j = j + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i + 1, j + 1), data, blkORwht, threshold))
                                                        {
                                                            i = i + 1;
                                                            j = j + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i = ps[ps.Count - 1].X - rec.X;
                                                                j = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (j == rec.Height - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i = errPs[errPs.Count - 1].X + 1 - rec.X;
                                                        j = errPs[errPs.Count - 1].Y - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if (blkORwht == 0)//黑边
                                        {
                                            if (B < threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (j == rec.Height - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i = ps[ps.Count - 1].X - 1 - rec.X;
                                                    if (i < 0)
                                                        i = 0;
                                                    j = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i = ps[ps.Count - 1].X - rec.X;
                                                        j = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i, j + 1), data, blkORwht, threshold))
                                                        {
                                                            j = j + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i + 1, j + 1), data, blkORwht, threshold))
                                                        {
                                                            i = i + 1;
                                                            j = j + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i = ps[ps.Count - 1].X - rec.X;
                                                                j = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (j == rec.Height - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i = errPs[errPs.Count - 1].X + 1 - rec.X;
                                                        j = errPs[errPs.Count - 1].Y - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                        else//白边
                                        {
                                            if (B > threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (j == rec.Height - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i = ps[ps.Count - 1].X - 1 - rec.X;
                                                    if (i < 0)
                                                        i = 0;
                                                    j = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i = ps[ps.Count - 1].X - rec.X;
                                                        j = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i, j + 1), data, blkORwht, threshold))
                                                        {
                                                            j = j + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i + 1, j + 1), data, blkORwht, threshold))
                                                        {
                                                            i = i + 1;
                                                            j = j + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i = ps[ps.Count - 1].X - rec.X;
                                                                j = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (j == rec.Height - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i = errPs[errPs.Count - 1].X + 1 - rec.X;
                                                        j = errPs[errPs.Count - 1].Y - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    break;
                                }
                                #endregion
                            }
                            lps1.Linked = linked;
                            lps1.PointList = ps;
                            #endregion
                            break;
                        case 2:
                            #region dir=2
                            int i1 = 0;
                            int j1 = 0;
                            while (i1 >= 0 && i1 < rec.Width && j1 >= 0 && j1 < rec.Height)
                            {
                                #region 判断逻辑：正南方向
                                if (linked1 || linked2 || linked3)
                                {
                                    break;
                                }
                                Point pp = new Point(i1 + rec.X, j1 + rec.Y);
                                if (path.IsVisible(pp) && !errPs.Contains(pp))
                                {
                                    //g.DrawEllipse(Pens.Blue, pp.X, pp.Y, 1, 1);
                                    p2 = (byte*)data.Scan0 + BPP * i1 + stride * j1;
                                    R = p2[2];
                                    G = p2[1];
                                    B = p2[0];
                                    if (R != G || G != B)
                                    {
                                        //尚未是灰度化，先计算灰度
                                        // gray = 0.3*R + 0.59*G + 0.11*B
                                        gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                                        //再提取关键点
                                        if (blkORwht == 0)//黑边
                                        {
                                            if (gray < threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (j1 == rec.Height - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i1 = ps[ps.Count - 1].X - rec.X;
                                                    j1 = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i1 = ps[ps.Count - 1].X - rec.X;
                                                        j1 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i1 - 1, j1 + 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 - 1;
                                                            j1 = j1 + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i1 + 1, j1 + 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 + 1;
                                                            j1 = j1 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i1 = ps[ps.Count - 1].X - rec.X;
                                                                j1 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (j1 == rec.Height - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i1 = errPs[errPs.Count - 1].X + 1 - rec.X;
                                                        j1 = errPs[errPs.Count - 1].Y - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                        else//白边
                                        {
                                            if (gray > threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (j1 == rec.Height - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i1 = ps[ps.Count - 1].X - rec.X;
                                                    j1 = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i1 = ps[ps.Count - 1].X - rec.X;
                                                        j1 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i1 - 1, j1 + 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 - 1;
                                                            j1 = j1 + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i1 + 1, j1 + 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 + 1;
                                                            j1 = j1 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i1 = ps[ps.Count - 1].X - rec.X;
                                                                j1 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (j1 == rec.Height - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i1 = errPs[errPs.Count - 1].X + 1 - rec.X;
                                                        j1 = errPs[errPs.Count - 1].Y - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if (blkORwht == 0)//黑边
                                        {
                                            if (B < threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (j1 == rec.Height - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i1 = ps[ps.Count - 1].X - rec.X;
                                                    j1 = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i1 = ps[ps.Count - 1].X - rec.X;
                                                        j1 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i1 - 1, j1 + 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 - 1;
                                                            j1 = j1 + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i1 + 1, j1 + 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 + 1;
                                                            j1 = j1 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i1 = ps[ps.Count - 1].X - rec.X;
                                                                j1 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (j1 == rec.Height - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i1 = errPs[errPs.Count - 1].X + 1 - rec.X;
                                                        j1 = errPs[errPs.Count - 1].Y - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                        else//白边
                                        {
                                            if (B > threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (j1 == rec.Height - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i1 = ps[ps.Count - 1].X - rec.X;
                                                    j1 = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i1 = ps[ps.Count - 1].X - rec.X;
                                                        j1 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i1 - 1, j1 + 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 - 1;
                                                            j1 = j1 + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i1 + 1, j1 + 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 + 1;
                                                            j1 = j1 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i1 = ps[ps.Count - 1].X - rec.X;
                                                                j1 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (j1 == rec.Height - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i1 = errPs[errPs.Count - 1].X + 1 - rec.X;
                                                        j1 = errPs[errPs.Count - 1].Y - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    break;
                                }
                                #endregion
                            }
                            lps2.Linked = linked;
                            lps2.PointList = ps;
                            #endregion
                            break;
                        case 3:
                            #region dir=3
                            int i2 = 0;
                            int j2 = 0;
                            while (i2 >= 0 && i2 < rec.Width && j2 >= 0 && j2 < rec.Height)
                            {
                                #region 判断逻辑：东南方向
                                if (linked1 || linked2 || linked3)
                                {
                                    break;
                                }
                                Point pp = new Point(i2 + rec.X, j2 + rec.Y);
                                if (path.IsVisible(pp) && !errPs.Contains(pp))
                                {
                                    //g.DrawEllipse(Pens.Blue, pp.X, pp.Y, 1, 1);
                                    p3 = (byte*)data.Scan0 + BPP * i2 + stride * j2;
                                    R = p3[2];
                                    G = p3[1];
                                    B = p3[0];
                                    if (R != G || G != B)
                                    {
                                        //尚未是灰度化，先计算灰度
                                        // gray = 0.3*R + 0.59*G + 0.11*B
                                        gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                                        //再提取关键点
                                        if (blkORwht == 0)//黑边
                                        {
                                            if (gray < threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (j2 == rec.Height - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i2 = ps[ps.Count - 1].X + 1 - rec.X;
                                                    if (i2 > rec.Width - 1)
                                                        i2 = rec.Width - 1;
                                                    j2 = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i2 = ps[ps.Count - 1].X - rec.X;
                                                        j2 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i2 - 1, j2 + 1), data, blkORwht, threshold))
                                                        {
                                                            i2 = i2 - 1;
                                                            j2 = j2 + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i2, j2 + 1), data, blkORwht, threshold))
                                                        {
                                                            j2 = j2 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i2 = ps[ps.Count - 1].X - rec.X;
                                                                j2 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (j2 == rec.Height - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i2 = errPs[errPs.Count - 1].X + 1 - rec.X;
                                                        j2 = errPs[errPs.Count - 1].Y - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                        else//白边
                                        {
                                            if (gray > threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (j2 == rec.Height - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i2 = ps[ps.Count - 1].X + 1 - rec.X;
                                                    if (i2 > rec.Width - 1)
                                                        i2 = rec.Width - 1;
                                                    j2 = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i2 = ps[ps.Count - 1].X - rec.X;
                                                        j2 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i2 - 1, j2 + 1), data, blkORwht, threshold))
                                                        {
                                                            i2 = i2 - 1;
                                                            j2 = j2 + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i2, j2 + 1), data, blkORwht, threshold))
                                                        {
                                                            j2 = j2 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i2 = ps[ps.Count - 1].X - rec.X;
                                                                j2 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (j2 == rec.Height - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i2 = errPs[errPs.Count - 1].X + 1 - rec.X;
                                                        j2 = errPs[errPs.Count - 1].Y - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if (blkORwht == 0)//黑边
                                        {
                                            if (B < threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (j2 == rec.Height - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i2 = ps[ps.Count - 1].X + 1 - rec.X;
                                                    if (i2 > rec.Width - 1)
                                                        i2 = rec.Width - 1;
                                                    j2 = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i2 = ps[ps.Count - 1].X - rec.X;
                                                        j2 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i2 - 1, j2 + 1), data, blkORwht, threshold))
                                                        {
                                                            i2 = i2 - 1;
                                                            j2 = j2 + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i2, j2 + 1), data, blkORwht, threshold))
                                                        {
                                                            j2 = j2 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i2 = ps[ps.Count - 1].X - rec.X;
                                                                j2 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (j2 == rec.Height - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i2 = errPs[errPs.Count - 1].X + 1 - rec.X;
                                                        j2 = errPs[errPs.Count - 1].Y - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                        else//白边
                                        {
                                            if (B > threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (j2 == rec.Height - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i2 = ps[ps.Count - 1].X + 1 - rec.X;
                                                    if (i2 > rec.Width - 1)
                                                        i2 = rec.Width - 1;
                                                    j2 = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i2 = ps[ps.Count - 1].X - rec.X;
                                                        j2 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i2 - 1, j2 + 1), data, blkORwht, threshold))
                                                        {
                                                            i2 = i2 - 1;
                                                            j2 = j2 + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i2, j2 + 1), data, blkORwht, threshold))
                                                        {
                                                            j2 = j2 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i2 = ps[ps.Count - 1].X - rec.X;
                                                                j2 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (j2 == rec.Height - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i2 = errPs[errPs.Count - 1].X + 1 - rec.X;
                                                        j2 = errPs[errPs.Count - 1].Y - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    break;
                                }
                                #endregion
                            }
                            lps3.Linked = linked;
                            lps3.PointList = ps;
                            #endregion
                            break;
                    }
                }
                else//水平方向的连通判断
                {
                    switch (direction)
                    {
                        case 1:
                            #region dir=1
                            int i = 0;
                            int j = 0;
                            while (i >= 0 && i < rec.Width && j >= 0 && j < rec.Height)
                            {
                                #region 判断逻辑：东北方向
                                if (linked1 || linked2 || linked3)
                                {
                                    break;
                                }
                                Point pp = new Point(i + rec.X, j + rec.Y);
                                if (path.IsVisible(pp) && !errPs.Contains(pp))
                                {
                                    //g.DrawEllipse(Pens.Blue, pp.X, pp.Y, 1, 1);
                                    p1 = (byte*)data.Scan0 + BPP * i + stride * j;
                                    R = p1[2];
                                    G = p1[1];
                                    B = p1[0];
                                    if (R != G || G != B)
                                    {
                                        //尚未是灰度化，先计算灰度
                                        // gray = 0.3*R + 0.59*G + 0.11*B
                                        gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                                        //再提取关键点
                                        if (blkORwht == 0)//黑边
                                        {
                                            if (gray < threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (i == rec.Width - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i = ps[ps.Count - 1].X + 1 - rec.X;
                                                    j = ps[ps.Count - 1].Y - 1 - rec.Y;
                                                    if (j < 0)
                                                        j = 0;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i = ps[ps.Count - 1].X - rec.X;
                                                        j = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i + 1, j), data, blkORwht, threshold))
                                                        {
                                                            i = i + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i + 1, j + 1), data, blkORwht, threshold))
                                                        {
                                                            i = i + 1;
                                                            j = j + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i = ps[ps.Count - 1].X - rec.X;
                                                                j = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (i == rec.Width - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i = errPs[errPs.Count - 1].X - rec.X;
                                                        j = errPs[errPs.Count - 1].Y + 1 - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                        else//白边
                                        {
                                            if (gray > threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (i == rec.Width - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i = ps[ps.Count - 1].X + 1 - rec.X;
                                                    j = ps[ps.Count - 1].Y - 1 - rec.Y;
                                                    if (j < 0)
                                                        j = 0;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i = ps[ps.Count - 1].X - rec.X;
                                                        j = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i + 1, j), data, blkORwht, threshold))
                                                        {
                                                            i = i + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i + 1, j + 1), data, blkORwht, threshold))
                                                        {
                                                            i = i + 1;
                                                            j = j + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i = ps[ps.Count - 1].X - rec.X;
                                                                j = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (i == rec.Width - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i = errPs[errPs.Count - 1].X - rec.X;
                                                        j = errPs[errPs.Count - 1].Y + 1 - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if (blkORwht == 0)//黑边
                                        {
                                            if (B < threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (i == rec.Width - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i = ps[ps.Count - 1].X + 1 - rec.X;
                                                    j = ps[ps.Count - 1].Y - 1 - rec.Y;
                                                    if (j < 0)
                                                        j = 0;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i = ps[ps.Count - 1].X - rec.X;
                                                        j = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i + 1, j), data, blkORwht, threshold))
                                                        {
                                                            i = i + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i + 1, j + 1), data, blkORwht, threshold))
                                                        {
                                                            i = i + 1;
                                                            j = j + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i = ps[ps.Count - 1].X - rec.X;
                                                                j = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (i == rec.Width - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i = errPs[errPs.Count - 1].X - rec.X;
                                                        j = errPs[errPs.Count - 1].Y + 1 - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                        else//白边
                                        {
                                            if (B > threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (i == rec.Width - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i = ps[ps.Count - 1].X + 1 - rec.X;
                                                    j = ps[ps.Count - 1].Y - 1 - rec.Y;
                                                    if (j < 0)
                                                        j = 0;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i = ps[ps.Count - 1].X - rec.X;
                                                        j = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i + 1, j), data, blkORwht, threshold))
                                                        {
                                                            i = i + 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i + 1, j + 1), data, blkORwht, threshold))
                                                        {
                                                            i = i + 1;
                                                            j = j + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i = ps[ps.Count - 1].X - rec.X;
                                                                j = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (i == rec.Width - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i = errPs[errPs.Count - 1].X - rec.X;
                                                        j = errPs[errPs.Count - 1].Y + 1 - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    break;
                                }
                                #endregion
                            }
                            lps1.Linked = linked;
                            lps1.PointList = ps;
                            #endregion
                            break;
                        case 2:
                            #region dir=2
                            int i1 = 0;
                            int j1 = 0;
                            while (i1 >= 0 && i1 < rec.Width && j1 >= 0 && j1 < rec.Height)
                            {
                                #region 判断逻辑：正东方向
                                if (linked1 || linked2 || linked3)
                                {
                                    break;
                                }
                                Point pp = new Point(i1 + rec.X, j1 + rec.Y);
                                if (path.IsVisible(pp) && !errPs.Contains(pp))
                                {
                                    //g.DrawEllipse(Pens.Blue, pp.X, pp.Y, 1, 1);
                                    p2 = (byte*)data.Scan0 + BPP * i1 + stride * j1;
                                    R = p2[2];
                                    G = p2[1];
                                    B = p2[0];
                                    if (R != G || G != B)
                                    {
                                        //尚未是灰度化，先计算灰度
                                        // gray = 0.3*R + 0.59*G + 0.11*B
                                        gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                                        //再提取关键点
                                        if (blkORwht == 0)//黑边
                                        {
                                            if (gray < threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (i1 == rec.Width - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i1 = ps[ps.Count - 1].X + 1 - rec.X;
                                                    j1 = ps[ps.Count - 1].Y - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i1 = ps[ps.Count - 1].X - rec.X;
                                                        j1 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i1 + 1, j1 - 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 + 1;
                                                            j1 = j1 - 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i1 + 1, j1 + 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 + 1;
                                                            j1 = j1 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i1 = ps[ps.Count - 1].X - rec.X;
                                                                j1 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (i1 == rec.Width - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i1 = errPs[errPs.Count - 1].X - rec.X;
                                                        j1 = errPs[errPs.Count - 1].Y + 1 - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                        else//白边
                                        {
                                            if (gray > threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (i1 == rec.Width - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i1 = ps[ps.Count - 1].X + 1 - rec.X;
                                                    j1 = ps[ps.Count - 1].Y - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i1 = ps[ps.Count - 1].X - rec.X;
                                                        j1 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i1 + 1, j1 - 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 + 1;
                                                            j1 = j1 - 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i1 + 1, j1 + 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 + 1;
                                                            j1 = j1 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i1 = ps[ps.Count - 1].X - rec.X;
                                                                j1 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (i1 == rec.Width - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i1 = errPs[errPs.Count - 1].X - rec.X;
                                                        j1 = errPs[errPs.Count - 1].Y + 1 - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if (blkORwht == 0)//黑边
                                        {
                                            if (B < threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (i1 == rec.Width - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i1 = ps[ps.Count - 1].X + 1 - rec.X;
                                                    j1 = ps[ps.Count - 1].Y - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i1 = ps[ps.Count - 1].X - rec.X;
                                                        j1 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i1 + 1, j1 - 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 + 1;
                                                            j1 = j1 - 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i1 + 1, j1 + 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 + 1;
                                                            j1 = j1 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i1 = ps[ps.Count - 1].X - rec.X;
                                                                j1 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (i1 == rec.Width - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i1 = errPs[errPs.Count - 1].X - rec.X;
                                                        j1 = errPs[errPs.Count - 1].Y + 1 - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                        else//白边
                                        {
                                            if (B > threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (i1 == rec.Width - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i1 = ps[ps.Count - 1].X + 1 - rec.X;
                                                    j1 = ps[ps.Count - 1].Y - rec.Y;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i1 = ps[ps.Count - 1].X - rec.X;
                                                        j1 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i1 + 1, j1 - 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 + 1;
                                                            j1 = j1 - 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i1 + 1, j1 + 1), data, blkORwht, threshold))
                                                        {
                                                            i1 = i1 + 1;
                                                            j1 = j1 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i1 = ps[ps.Count - 1].X - rec.X;
                                                                j1 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (i1 == rec.Width - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i1 = errPs[errPs.Count - 1].X - rec.X;
                                                        j1 = errPs[errPs.Count - 1].Y + 1 - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    break;
                                }
                                #endregion
                            }
                            lps2.Linked = linked;
                            lps2.PointList = ps;
                            #endregion
                            break;
                        case 3:
                            #region dir=3
                            int i2 = 0;
                            int j2 = 0;
                            while (i2 >= 0 && i2 < rec.Width && j2 >= 0 && j2 < rec.Height)
                            {
                                #region 判断逻辑：东南方向
                                if (linked1 || linked2 || linked3)
                                {
                                    break;
                                }
                                Point pp = new Point(i2 + rec.X, j2 + rec.Y);
                                if (path.IsVisible(pp) && !errPs.Contains(pp))
                                {
                                    //g.DrawEllipse(Pens.Blue, pp.X, pp.Y, 1, 1);
                                    p3 = (byte*)data.Scan0 + BPP * i2 + stride * j2;
                                    R = p3[2];
                                    G = p3[1];
                                    B = p3[0];
                                    if (R != G || G != B)
                                    {
                                        //尚未是灰度化，先计算灰度
                                        // gray = 0.3*R + 0.59*G + 0.11*B
                                        gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                                        //再提取关键点
                                        if (blkORwht == 0)//黑边
                                        {
                                            if (gray < threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (i2 == rec.Width - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i2 = ps[ps.Count - 1].X + 1 - rec.X;
                                                    j2 = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                    if (j2 > rec.Height - 1)
                                                        j2 = rec.Height - 1;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i2 = ps[ps.Count - 1].X - rec.X;
                                                        j2 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i2 + 1, j2 - 1), data, blkORwht, threshold))
                                                        {
                                                            i2 = i2 + 1;
                                                            j2 = j2 - 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i2 + 1, j2), data, blkORwht, threshold))
                                                        {
                                                            i2 = i2 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i2 = ps[ps.Count - 1].X - rec.X;
                                                                j2 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (i2 == rec.Width - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i2 = errPs[errPs.Count - 1].X - rec.X;
                                                        j2 = errPs[errPs.Count - 1].Y + 1 - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                        else//白边
                                        {
                                            if (gray > threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (i2 == rec.Width - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i2 = ps[ps.Count - 1].X + 1 - rec.X;
                                                    j2 = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                    if (j2 > rec.Height - 1)
                                                        j2 = rec.Height - 1;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i2 = ps[ps.Count - 1].X - rec.X;
                                                        j2 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i2 + 1, j2 - 1), data, blkORwht, threshold))
                                                        {
                                                            i2 = i2 + 1;
                                                            j2 = j2 - 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i2 + 1, j2), data, blkORwht, threshold))
                                                        {
                                                            i2 = i2 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i2 = ps[ps.Count - 1].X - rec.X;
                                                                j2 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (i2 == rec.Width - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i2 = errPs[errPs.Count - 1].X - rec.X;
                                                        j2 = errPs[errPs.Count - 1].Y + 1 - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if (blkORwht == 0)//黑边
                                        {
                                            if (B < threshold && LinkedPoints(pp, ps[ps.Count - 1]))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (i2 == rec.Width)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i2 = ps[ps.Count - 1].X + 1 - rec.X;
                                                    j2 = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                    if (j2 > rec.Height - 1)
                                                        j2 = rec.Height - 1;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i2 = ps[ps.Count - 1].X - rec.X;
                                                        j2 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i2 + 1, j2 - 1), data, blkORwht, threshold))
                                                        {
                                                            i2 = i2 + 1;
                                                            j2 = j2 - 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i2 + 1, j2), data, blkORwht, threshold))
                                                        {
                                                            i2 = i2 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i2 = ps[ps.Count - 1].X - rec.X;
                                                                j2 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (i2 == rec.Width - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i2 = errPs[errPs.Count - 1].X - rec.X;
                                                        j2 = errPs[errPs.Count - 1].Y + 1 - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                        else//白边
                                        {
                                            if (B > threshold && (ps.Count == 0 ? true : LinkedPoints(pp, ps[ps.Count - 1])))
                                            {
                                                if (!ps.Contains(pp))
                                                    ps.Add(pp);
                                                if (i2 == rec.Width - 1)
                                                {
                                                    linked = true;
                                                    break;
                                                }
                                                else
                                                {
                                                    //找到后一个要检查的点
                                                    i2 = ps[ps.Count - 1].X + 1 - rec.X;
                                                    j2 = ps[ps.Count - 1].Y + 1 - rec.Y;
                                                    if (j2 > rec.Height - 1)
                                                        j2 = rec.Height - 1;
                                                }
                                            }
                                            else
                                            {
                                                if (ps.Count > 0)
                                                {
                                                    if (!errPs.Contains(ps[ps.Count - 1]))
                                                        errPs.Add(ps[ps.Count - 1]);
                                                    ps.RemoveAt(ps.Count - 1);
                                                    //返回前一个连通点
                                                    if (ps.Count > 0)
                                                    {
                                                        i2 = ps[ps.Count - 1].X - rec.X;
                                                        j2 = ps[ps.Count - 1].Y - rec.Y;
                                                        if (CheckThePointIsKeyPoint(new Point(i2 + 1, j2 - 1), data, blkORwht, threshold))
                                                        {
                                                            i2 = i2 + 1;
                                                            j2 = j2 - 1;
                                                        }
                                                        else if (CheckThePointIsKeyPoint(new Point(i2 + 1, j2), data, blkORwht, threshold))
                                                        {
                                                            i2 = i2 + 1;
                                                        }
                                                        else
                                                        {
                                                            if (!errPs.Contains(ps[ps.Count - 1]))
                                                                errPs.Add(ps[ps.Count - 1]);
                                                            ps.RemoveAt(ps.Count - 1);
                                                            if (ps.Count > 0)
                                                            {
                                                                i2 = ps[ps.Count - 1].X - rec.X;
                                                                j2 = ps[ps.Count - 1].Y - rec.Y;
                                                            }
                                                        }
                                                        if (i2 == rec.Width - 1)
                                                        {
                                                            linked = true;
                                                            break;
                                                        }
                                                    }
                                                    else
                                                    {
                                                        i2 = errPs[errPs.Count - 1].X - rec.X;
                                                        j2 = errPs[errPs.Count - 1].Y + 1 - rec.Y;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    break;
                                }
                                #endregion
                            }
                            lps3.Linked = linked;
                            lps3.PointList = ps;
                            #endregion
                            break;
                    }
                }
                bm.UnlockBits(data);
            }
        }
    }
    /// <summary>
    /// 4线程类
    /// </summary>
    public static class MultiThreading4
    {
        /// <summary>
        /// 第1区域的是否已存在关键点
        /// </summary>
        public static bool exist1;
        /// <summary>
        /// 第2区域的是否已存在关键点
        /// </summary>
        public static bool exist2;
        /// <summary>
        /// 第3区域的是否已存在关键点
        /// </summary>
        public static bool exist3;
        /// <summary>
        /// 第4区域的是否已存在关键点
        /// </summary>
        public static bool exist4;
        /// <summary>
        /// 第1区域的指定矩形范围内是否已存在关键点
        /// </summary>
        public static bool existRect1;
        /// <summary>
        /// 第2区域的指定矩形范围内是否已存在关键点
        /// </summary>
        public static bool existRect2;
        /// <summary>
        /// 第3区域的指定矩形范围内是否已存在关键点
        /// </summary>
        public static bool existRect3;
        /// <summary>
        /// 第4区域的指定矩形范围内是否已存在关键点
        /// </summary>
        public static bool existRect4;
        /// <summary>
        /// 第1区域的指定路径范围内是否已存在关键点
        /// </summary>
        public static bool existPath1;
        /// <summary>
        /// 第2区域的指定路径范围内是否已存在关键点
        /// </summary>
        public static bool existPath2;
        /// <summary>
        /// 第3区域的指定路径范围内是否已存在关键点
        /// </summary>
        public static bool existPath3;
        /// <summary>
        /// 第4区域的指定路径范围内是否已存在关键点
        /// </summary>
        public static bool existPath4;
        /// <summary>
        /// 第1区域的指定区域范围内是否已存在关键点
        /// </summary>
        public static bool existRegion1;
        /// <summary>
        /// 第2区域的指定区域范围内是否已存在关键点
        /// </summary>
        public static bool existRegion2;
        /// <summary>
        /// 第3区域的指定区域范围内是否已存在关键点
        /// </summary>
        public static bool existRegion3;
        /// <summary>
        /// 第4区域的指定区域范围内是否已存在关键点
        /// </summary>
        public static bool existRegion4;
        /// <summary>
        /// 第1区域是否存在关键点回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTask1(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out exist1, ar);
        }
        /// <summary>
        /// 第2区域是否存在关键点回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTask2(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out exist2, ar);
        }
        /// <summary>
        /// 第3区域是否存在关键点回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTask3(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out exist3, ar);
        }
        /// <summary>
        /// 第4区域是否存在关键点回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTask4(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out exist4, ar);
        }
        /// <summary>
        /// 第1区域是否存在关键点回调函数(位图中某矩形范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRect1(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out existRect1, ar);
        }
        /// <summary>
        /// 第2区域是否存在关键点回调函数(位图中某矩形范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRect2(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out existRect2, ar);
        }
        /// <summary>
        /// 第3区域是否存在关键点回调函数(位图中某矩形范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRect3(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out existRect3, ar);
        }
        /// <summary>
        /// 第4区域是否存在关键点回调函数(位图中某矩形范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRect4(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out existRect4, ar);
        }
        /// <summary>
        /// 第1区域是否存在关键点回调函数(位图中某路径范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithPath1(IAsyncResult ar)
        {
            AsyncMethodCallerWithPath caller = (AsyncMethodCallerWithPath)ar.AsyncState;
            caller.EndInvoke(out existPath1, ar);
        }
        /// <summary>
        /// 第2区域是否存在关键点回调函数(位图中某路径范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithPath2(IAsyncResult ar)
        {
            AsyncMethodCallerWithPath caller = (AsyncMethodCallerWithPath)ar.AsyncState;
            caller.EndInvoke(out existPath2, ar);
        }
        /// <summary>
        /// 第3区域是否存在关键点回调函数(位图中某路径范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithPath3(IAsyncResult ar)
        {
            AsyncMethodCallerWithPath caller = (AsyncMethodCallerWithPath)ar.AsyncState;
            caller.EndInvoke(out existPath3, ar);
        }
        /// <summary>
        /// 第4区域是否存在关键点回调函数(位图中某路径范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithPath4(IAsyncResult ar)
        {
            AsyncMethodCallerWithPath caller = (AsyncMethodCallerWithPath)ar.AsyncState;
            caller.EndInvoke(out existPath4, ar);
        }
        /// <summary>
        /// 第1区域是否存在关键点回调函数(位图中某区域内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRegion1(IAsyncResult ar)
        {
            AsyncMethodCallerWithRegion caller = (AsyncMethodCallerWithRegion)ar.AsyncState;
            caller.EndInvoke(out existRegion1, ar);
        }
        /// <summary>
        /// 第2区域是否存在关键点回调函数(位图中某区域内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRegion2(IAsyncResult ar)
        {
            AsyncMethodCallerWithRegion caller = (AsyncMethodCallerWithRegion)ar.AsyncState;
            caller.EndInvoke(out existRegion2, ar);
        }
        /// <summary>
        /// 第3区域是否存在关键点回调函数(位图中某区域内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRegion3(IAsyncResult ar)
        {
            AsyncMethodCallerWithRegion caller = (AsyncMethodCallerWithRegion)ar.AsyncState;
            caller.EndInvoke(out existRegion3, ar);
        }
        /// <summary>
        /// 第4区域是否存在关键点回调函数(位图中某区域内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRegion4(IAsyncResult ar)
        {
            AsyncMethodCallerWithRegion caller = (AsyncMethodCallerWithRegion)ar.AsyncState;
            caller.EndInvoke(out existRegion4, ar);
        }
        /// <summary>
        /// 检测指定路径范围内是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="exist"></param>
        public static void CheckSomeRectWithPath(RectBitmapWithPath bmp, out bool exist)
        {
            exist = false;
            RectBitmapWithPath rectbmp = (RectBitmapWithPath)bmp;
            Bitmap bm = rectbmp.bmp;
            Rectangle rec = rectbmp.rect;
            int blkORwht = rectbmp.blkORwht;
            byte threshold = rectbmp.threshold;
            GraphicsPath path = rectbmp.path;
            int BPP = Image.GetPixelFormatSize(bm.PixelFormat) / 8;
            unsafe
            {
                byte gray, R, G, B;
                BitmapData data = bm.LockBits(new Rectangle(0, 0, rec.Width, rec.Height), ImageLockMode.ReadOnly, bm.PixelFormat);
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * rec.Width;
                for (int j = rec.Y; j < rec.Height + rec.Y; j++)
                {
                    for (int i = rec.X; i < rec.Width + rec.X; i++)
                    {
                        if (existPath1 || existPath2 || existPath3 || existPath4)
                        {
                            exist = true;
                            break;
                        }
                        if (path.IsVisible(i, j))
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
                                if (blkORwht == 0)//黑边
                                {
                                    if (gray < threshold)
                                    {
                                        exist = true;
                                        break;
                                    }
                                }
                                else//白边
                                {
                                    if (gray > threshold)
                                    {
                                        exist = true;
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
                                        exist = true;
                                        break;
                                    }
                                }
                                else//白边
                                {
                                    if (B > threshold)
                                    {
                                        exist = true;
                                        break;
                                    }
                                }
                            }
                        }
                        p += BPP;
                    }
                    if (exist)
                        break;
                    p += offset;
                }
                bm.UnlockBits(data);
            }
        }
        /// <summary>
        /// 检测指定区域范围内是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="exist"></param>
        public static void CheckSomeRectWithRegion(RectBitmapWithRegion bmp, out bool exist)
        {
            exist = false;
            RectBitmapWithRegion rectbmp = (RectBitmapWithRegion)bmp;
            Bitmap bm = rectbmp.bmp;
            Rectangle rec = rectbmp.rect;
            int blkORwht = rectbmp.blkORwht;
            byte threshold = rectbmp.threshold;
            Region reg = rectbmp.region;
            int BPP = Image.GetPixelFormatSize(bm.PixelFormat) / 8;
            unsafe
            {
                byte gray, R, G, B;
                BitmapData data = bm.LockBits(new Rectangle(0, 0, rec.Width, rec.Height), ImageLockMode.ReadOnly, bm.PixelFormat);
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * rec.Width;
                for (int j = rec.Y; j < rec.Height + rec.Y; j++)
                {
                    for (int i = rec.X; i < rec.Width + rec.X; i++)
                    {
                        if (existRegion1 || existRegion2 || existRegion3 || existRegion4)
                        {
                            exist = true;
                            break;
                        }
                        if (reg.IsVisible(i, j))
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
                                if (blkORwht == 0)//黑边
                                {
                                    if (gray < threshold)
                                    {
                                        exist = true;
                                        break;
                                    }
                                }
                                else//白边
                                {
                                    if (gray > threshold)
                                    {
                                        exist = true;
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
                                        exist = true;
                                        break;
                                    }
                                }
                                else//白边
                                {
                                    if (B > threshold)
                                    {
                                        exist = true;
                                        break;
                                    }
                                }
                            }
                        }
                        p += BPP;
                    }
                    if (exist)
                        break;
                    p += offset;
                }
                bm.UnlockBits(data);
            }
        }
        /// <summary>
        /// 检测整个位图范围内是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="exist"></param>
        public static void CheckWholeRect(RectBitmap bmp, out bool exist)
        {
            exist = false;
            RectBitmap rectbmp = (RectBitmap)bmp;
            Bitmap bm = rectbmp.bmp;
            Rectangle rec = rectbmp.rect;
            int blkORwht = rectbmp.blkORwht;
            byte threshold = rectbmp.threshold;
            int BPP = Image.GetPixelFormatSize(bm.PixelFormat) / 8;
            unsafe
            {
                byte gray, R, G, B;
                BitmapData data = bm.LockBits(new Rectangle(0, 0, rec.Width, rec.Height), ImageLockMode.ReadOnly, bm.PixelFormat);
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * rec.Width;
                for (int j = rec.Y; j < rec.Height + rec.Y; j++)
                {
                    for (int i = rec.X; i < rec.Width + rec.X; i++)
                    {
                        if (exist1 || exist2 || exist3 || exist4)
                        {
                            exist = true;
                            break;
                        }
                        R = p[2];
                        G = p[1];
                        B = p[0];
                        if (R != G || G != B)
                        {
                            //尚未是灰度化，先计算灰度
                            // gray = 0.3*R + 0.59*G + 0.11*B
                            gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                            //再提取关键点
                            if (blkORwht == 0)//黑边
                            {
                                if (gray < threshold)
                                {
                                    exist = true;
                                    break;
                                }
                            }
                            else//白边
                            {
                                if (gray > threshold)
                                {
                                    exist = true;
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
                                    exist = true;
                                    break;
                                }
                            }
                            else//白边
                            {
                                if (B > threshold)
                                {
                                    exist = true;
                                    break;
                                }
                            }
                        }
                        p += BPP;
                    }
                    if (exist)
                        break;
                    p += offset;
                }
                bm.UnlockBits(data);
            }
        }
        /// <summary>
        /// 检测指定矩形范围内是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="exist"></param>
        public static void CheckSomeRect(RectBitmap bmp, out bool exist)
        {
            exist = false;
            RectBitmap rectbmp = (RectBitmap)bmp;
            Bitmap bm = rectbmp.bmp;
            Rectangle rec = rectbmp.rect;
            int blkORwht = rectbmp.blkORwht;
            byte threshold = rectbmp.threshold;
            int BPP = Image.GetPixelFormatSize(bm.PixelFormat) / 8;
            unsafe
            {
                byte gray, R, G, B;
                BitmapData data = bm.LockBits(new Rectangle(0, 0, rec.Width, rec.Height), ImageLockMode.ReadOnly, bm.PixelFormat);
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * rec.Width;
                for (int j = rec.Y; j < rec.Height + rec.Y; j++)
                {
                    for (int i = rec.X; i < rec.Width + rec.X; i++)
                    {
                        if (existRect1 || existRect2 || existRect3 || existRect4)
                        {
                            exist = true;
                            break;
                        }
                        R = p[2];
                        G = p[1];
                        B = p[0];
                        if (R != G || G != B)
                        {
                            //尚未是灰度化，先计算灰度
                            // gray = 0.3*R + 0.59*G + 0.11*B
                            gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                            //再提取关键点
                            if (blkORwht == 0)//黑边
                            {
                                if (gray < threshold)
                                {
                                    exist = true;
                                    break;
                                }
                            }
                            else//白边
                            {
                                if (gray > threshold)
                                {
                                    exist = true;
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
                                    exist = true;
                                    break;
                                }
                            }
                            else//白边
                            {
                                if (B > threshold)
                                {
                                    exist = true;
                                    break;
                                }
                            }
                        }
                        p += BPP;
                    }
                    if (exist)
                        break;
                    p += offset;
                }
                bm.UnlockBits(data);
            }
        }
    }
    /// <summary>
    /// 9线程类
    /// </summary>
    public static class MultiThreading9
    {
        /// <summary>
        /// 第1区域的是否已存在关键点
        /// </summary>
        public static bool exist1;
        /// <summary>
        /// 第2区域的是否已存在关键点
        /// </summary>
        public static bool exist2;
        /// <summary>
        /// 第3区域的是否已存在关键点
        /// </summary>
        public static bool exist3;
        /// <summary>
        /// 第4区域的是否已存在关键点
        /// </summary>
        public static bool exist4;
        /// <summary>
        /// 第5区域的是否已存在关键点
        /// </summary>
        public static bool exist5;
        /// <summary>
        /// 第6区域的是否已存在关键点
        /// </summary>
        public static bool exist6;
        /// <summary>
        /// 第7区域的是否已存在关键点
        /// </summary>
        public static bool exist7;
        /// <summary>
        /// 第8区域的是否已存在关键点
        /// </summary>
        public static bool exist8;
        /// <summary>
        /// 第9区域的是否已存在关键点
        /// </summary>
        public static bool exist9;
        /// <summary>
        /// 第1区域的指定矩形范围内是否已存在关键点
        /// </summary>
        public static bool existRect1;
        /// <summary>
        /// 第2区域的指定矩形范围内是否已存在关键点
        /// </summary>
        public static bool existRect2;
        /// <summary>
        /// 第3区域的指定矩形范围内是否已存在关键点
        /// </summary>
        public static bool existRect3;
        /// <summary>
        /// 第4区域的指定矩形范围内是否已存在关键点
        /// </summary>
        public static bool existRect4;
        /// <summary>
        /// 第5区域的指定矩形范围内是否已存在关键点
        /// </summary>
        public static bool existRect5;
        /// <summary>
        /// 第6区域的指定矩形范围内是否已存在关键点
        /// </summary>
        public static bool existRect6;
        /// <summary>
        /// 第7区域的指定矩形范围内是否已存在关键点
        /// </summary>
        public static bool existRect7;
        /// <summary>
        /// 第8区域的指定矩形范围内是否已存在关键点
        /// </summary>
        public static bool existRect8;
        /// <summary>
        /// 第9区域的指定矩形范围内是否已存在关键点
        /// </summary>
        public static bool existRect9;
        /// <summary>
        /// 第1区域的指定路径范围内是否已存在关键点
        /// </summary>
        public static bool existPath1;
        /// <summary>
        /// 第2区域的指定路径范围内是否已存在关键点
        /// </summary>
        public static bool existPath2;
        /// <summary>
        /// 第3区域的指定路径范围内是否已存在关键点
        /// </summary>
        public static bool existPath3;
        /// <summary>
        /// 第4区域的指定路径范围内是否已存在关键点
        /// </summary>
        public static bool existPath4;
        /// <summary>
        /// 第5区域的指定路径范围内是否已存在关键点
        /// </summary>
        public static bool existPath5;
        /// <summary>
        /// 第6区域的指定路径范围内是否已存在关键点
        /// </summary>
        public static bool existPath6;
        /// <summary>
        /// 第7区域的指定路径范围内是否已存在关键点
        /// </summary>
        public static bool existPath7;
        /// <summary>
        /// 第8区域的指定路径范围内是否已存在关键点
        /// </summary>
        public static bool existPath8;
        /// <summary>
        /// 第9区域的指定路径范围内是否已存在关键点
        /// </summary>
        public static bool existPath9;
        /// <summary>
        /// 第1区域的指定区域范围内是否已存在关键点
        /// </summary>
        public static bool existRegion1;
        /// <summary>
        /// 第2区域的指定区域范围内是否已存在关键点
        /// </summary>
        public static bool existRegion2;
        /// <summary>
        /// 第3区域的指定区域范围内是否已存在关键点
        /// </summary>
        public static bool existRegion3;
        /// <summary>
        /// 第4区域的指定区域范围内是否已存在关键点
        /// </summary>
        public static bool existRegion4;
        /// <summary>
        /// 第5区域的指定区域范围内是否已存在关键点
        /// </summary>
        public static bool existRegion5;
        /// <summary>
        /// 第6区域的指定区域范围内是否已存在关键点
        /// </summary>
        public static bool existRegion6;
        /// <summary>
        /// 第7区域的指定区域范围内是否已存在关键点
        /// </summary>
        public static bool existRegion7;
        /// <summary>
        /// 第8区域的指定区域范围内是否已存在关键点
        /// </summary>
        public static bool existRegion8;
        /// <summary>
        /// 第9区域的指定区域范围内是否已存在关键点
        /// </summary>
        public static bool existRegion9;
        /// <summary>
        /// 第1区域是否存在关键点回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTask1(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out exist1, ar);
        }
        /// <summary>
        /// 第2区域是否存在关键点回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTask2(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out exist2, ar);
        }
        /// <summary>
        /// 第3区域是否存在关键点回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTask3(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out exist3, ar);
        }
        /// <summary>
        /// 第4区域是否存在关键点回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTask4(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out exist4, ar);
        }
        /// <summary>
        /// 第5区域是否存在关键点回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTask5(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out exist5, ar);
        }
        /// <summary>
        /// 第6区域是否存在关键点回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTask6(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out exist6, ar);
        }
        /// <summary>
        /// 第7区域是否存在关键点回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTask7(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out exist7, ar);
        }
        /// <summary>
        /// 第8区域是否存在关键点回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTask8(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out exist8, ar);
        }
        /// <summary>
        /// 第9区域是否存在关键点回调函数
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTask9(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out exist9, ar);
        }
        /// <summary>
        /// 第1区域是否存在关键点回调函数(位图中某矩形范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRect1(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out existRect1, ar);
        }
        /// <summary>
        /// 第2区域是否存在关键点回调函数(位图中某矩形范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRect2(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out existRect2, ar);
        }
        /// <summary>
        /// 第3区域是否存在关键点回调函数(位图中某矩形范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRect3(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out existRect3, ar);
        }
        /// <summary>
        /// 第4区域是否存在关键点回调函数(位图中某矩形范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRect4(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out existRect4, ar);
        }
        /// <summary>
        /// 第5区域是否存在关键点回调函数(位图中某矩形范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRect5(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out existRect5, ar);
        }
        /// <summary>
        /// 第6区域是否存在关键点回调函数(位图中某矩形范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRect6(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out existRect6, ar);
        }
        /// <summary>
        /// 第7区域是否存在关键点回调函数(位图中某矩形范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRect7(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out existRect7, ar);
        }
        /// <summary>
        /// 第8区域是否存在关键点回调函数(位图中某矩形范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRect8(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out existRect8, ar);
        }
        /// <summary>
        /// 第9区域是否存在关键点回调函数(位图中某矩形范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRect9(IAsyncResult ar)
        {
            AsyncMethodCaller caller = (AsyncMethodCaller)ar.AsyncState;
            caller.EndInvoke(out existRect9, ar);
        }
        /// <summary>
        /// 第1区域是否存在关键点回调函数(位图中某路径范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithPath1(IAsyncResult ar)
        {
            AsyncMethodCallerWithPath caller = (AsyncMethodCallerWithPath)ar.AsyncState;
            caller.EndInvoke(out existPath1, ar);
        }
        /// <summary>
        /// 第2区域是否存在关键点回调函数(位图中某路径范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithPath2(IAsyncResult ar)
        {
            AsyncMethodCallerWithPath caller = (AsyncMethodCallerWithPath)ar.AsyncState;
            caller.EndInvoke(out existPath2, ar);
        }
        /// <summary>
        /// 第3区域是否存在关键点回调函数(位图中某路径范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithPath3(IAsyncResult ar)
        {
            AsyncMethodCallerWithPath caller = (AsyncMethodCallerWithPath)ar.AsyncState;
            caller.EndInvoke(out existPath3, ar);
        }
        /// <summary>
        /// 第4区域是否存在关键点回调函数(位图中某路径范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithPath4(IAsyncResult ar)
        {
            AsyncMethodCallerWithPath caller = (AsyncMethodCallerWithPath)ar.AsyncState;
            caller.EndInvoke(out existPath4, ar);
        }
        /// <summary>
        /// 第5区域是否存在关键点回调函数(位图中某路径范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithPath5(IAsyncResult ar)
        {
            AsyncMethodCallerWithPath caller = (AsyncMethodCallerWithPath)ar.AsyncState;
            caller.EndInvoke(out existPath5, ar);
        }
        /// <summary>
        /// 第6区域是否存在关键点回调函数(位图中某路径范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithPath6(IAsyncResult ar)
        {
            AsyncMethodCallerWithPath caller = (AsyncMethodCallerWithPath)ar.AsyncState;
            caller.EndInvoke(out existPath6, ar);
        }
        /// <summary>
        /// 第7区域是否存在关键点回调函数(位图中某路径范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithPath7(IAsyncResult ar)
        {
            AsyncMethodCallerWithPath caller = (AsyncMethodCallerWithPath)ar.AsyncState;
            caller.EndInvoke(out existPath7, ar);
        }
        /// <summary>
        /// 第8区域是否存在关键点回调函数(位图中某路径范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithPath8(IAsyncResult ar)
        {
            AsyncMethodCallerWithPath caller = (AsyncMethodCallerWithPath)ar.AsyncState;
            caller.EndInvoke(out existPath8, ar);
        }
        /// <summary>
        /// 第9区域是否存在关键点回调函数(位图中某路径范围内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithPath9(IAsyncResult ar)
        {
            AsyncMethodCallerWithPath caller = (AsyncMethodCallerWithPath)ar.AsyncState;
            caller.EndInvoke(out existPath9, ar);
        }
        /// <summary>
        /// 第1区域是否存在关键点回调函数(位图中某区域内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRegion1(IAsyncResult ar)
        {
            AsyncMethodCallerWithRegion caller = (AsyncMethodCallerWithRegion)ar.AsyncState;
            caller.EndInvoke(out existRegion1, ar);
        }
        /// <summary>
        /// 第2区域是否存在关键点回调函数(位图中某区域内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRegion2(IAsyncResult ar)
        {
            AsyncMethodCallerWithRegion caller = (AsyncMethodCallerWithRegion)ar.AsyncState;
            caller.EndInvoke(out existRegion2, ar);
        }
        /// <summary>
        /// 第3区域是否存在关键点回调函数(位图中某区域内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRegion3(IAsyncResult ar)
        {
            AsyncMethodCallerWithRegion caller = (AsyncMethodCallerWithRegion)ar.AsyncState;
            caller.EndInvoke(out existRegion3, ar);
        }
        /// <summary>
        /// 第4区域是否存在关键点回调函数(位图中某区域内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRegion4(IAsyncResult ar)
        {
            AsyncMethodCallerWithRegion caller = (AsyncMethodCallerWithRegion)ar.AsyncState;
            caller.EndInvoke(out existRegion4, ar);
        }
        /// <summary>
        /// 第5区域是否存在关键点回调函数(位图中某区域内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRegion5(IAsyncResult ar)
        {
            AsyncMethodCallerWithRegion caller = (AsyncMethodCallerWithRegion)ar.AsyncState;
            caller.EndInvoke(out existRegion5, ar);
        }
        /// <summary>
        /// 第6区域是否存在关键点回调函数(位图中某区域内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRegion6(IAsyncResult ar)
        {
            AsyncMethodCallerWithRegion caller = (AsyncMethodCallerWithRegion)ar.AsyncState;
            caller.EndInvoke(out existRegion6, ar);
        }
        /// <summary>
        /// 第7区域是否存在关键点回调函数(位图中某区域内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRegion7(IAsyncResult ar)
        {
            AsyncMethodCallerWithRegion caller = (AsyncMethodCallerWithRegion)ar.AsyncState;
            caller.EndInvoke(out existRegion7, ar);
        }
        /// <summary>
        /// 第8区域是否存在关键点回调函数(位图中某区域内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRegion8(IAsyncResult ar)
        {
            AsyncMethodCallerWithRegion caller = (AsyncMethodCallerWithRegion)ar.AsyncState;
            caller.EndInvoke(out existRegion8, ar);
        }
        /// <summary>
        /// 第9区域是否存在关键点回调函数(位图中某区域内)
        /// </summary>
        /// <param name="ar"></param>
        public static void CallbackTaskWithRegion9(IAsyncResult ar)
        {
            AsyncMethodCallerWithRegion caller = (AsyncMethodCallerWithRegion)ar.AsyncState;
            caller.EndInvoke(out existRegion9, ar);
        }
        /// <summary>
        /// 检测指定路径范围内是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="exist"></param>
        public static void CheckSomeRectWithPath(RectBitmapWithPath bmp, out bool exist)
        {
            exist = false;
            RectBitmapWithPath rectbmp = (RectBitmapWithPath)bmp;
            Bitmap bm = rectbmp.bmp;
            Rectangle rec = rectbmp.rect;
            int blkORwht = rectbmp.blkORwht;
            byte threshold = rectbmp.threshold;
            GraphicsPath path = rectbmp.path;
            int BPP = Image.GetPixelFormatSize(bm.PixelFormat) / 8;
            unsafe
            {
                byte gray, R, G, B;
                BitmapData data = bm.LockBits(new Rectangle(0, 0, rec.Width, rec.Height), ImageLockMode.ReadOnly, bm.PixelFormat);
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * rec.Width;
                for (int j = rec.Y; j < rec.Height + rec.Y; j++)
                {
                    for (int i = rec.X; i < rec.Width + rec.X; i++)
                    {
                        if (existPath1 || existPath2 || existPath3 || existPath4 || existPath5 || existPath6 || existPath7 || existPath8 || existPath9)
                        {
                            exist = true;
                            break;
                        }
                        if (path.IsVisible(i, j))
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
                                if (blkORwht == 0)//黑边
                                {
                                    if (gray < threshold)
                                    {
                                        exist = true;
                                        break;
                                    }
                                }
                                else//白边
                                {
                                    if (gray > threshold)
                                    {
                                        exist = true;
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
                                        exist = true;
                                        break;
                                    }
                                }
                                else//白边
                                {
                                    if (B > threshold)
                                    {
                                        exist = true;
                                        break;
                                    }
                                }
                            }
                        }
                        p += BPP;
                    }
                    if (exist)
                        break;
                    p += offset;
                }
                bm.UnlockBits(data);
            }
        }
        /// <summary>
        /// 检测指定区域范围内是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="exist"></param>
        public static void CheckSomeRectWithRegion(RectBitmapWithRegion bmp, out bool exist)
        {
            exist = false;
            RectBitmapWithRegion rectbmp = (RectBitmapWithRegion)bmp;
            Bitmap bm = rectbmp.bmp;
            Rectangle rec = rectbmp.rect;
            int blkORwht = rectbmp.blkORwht;
            byte threshold = rectbmp.threshold;
            Region reg = rectbmp.region;
            int BPP = Image.GetPixelFormatSize(bm.PixelFormat) / 8;
            unsafe
            {
                byte gray, R, G, B;
                BitmapData data = bm.LockBits(new Rectangle(0, 0, rec.Width, rec.Height), ImageLockMode.ReadOnly, bm.PixelFormat);
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * rec.Width;
                for (int j = rec.Y; j < rec.Height + rec.Y; j++)
                {
                    for (int i = rec.X; i < rec.Width + rec.X; i++)
                    {
                        if (existRegion1 || existRegion2 || existRegion3 || existRegion4 || existRegion5 || existRegion6 || existRegion7 || existRegion8 || existRegion9)
                        {
                            exist = true;
                            break;
                        }
                        if (reg.IsVisible(i, j))
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
                                if (blkORwht == 0)//黑边
                                {
                                    if (gray < threshold)
                                    {
                                        exist = true;
                                        break;
                                    }
                                }
                                else//白边
                                {
                                    if (gray > threshold)
                                    {
                                        exist = true;
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
                                        exist = true;
                                        break;
                                    }
                                }
                                else//白边
                                {
                                    if (B > threshold)
                                    {
                                        exist = true;
                                        break;
                                    }
                                }
                            }
                        }
                        p += BPP;
                    }
                    if (exist)
                        break;
                    p += offset;
                }
                bm.UnlockBits(data);
            }
        }
        /// <summary>
        /// 检测整个位图范围内是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="exist"></param>
        public static void CheckWholeRect(RectBitmap bmp, out bool exist)
        {
            exist = false;
            RectBitmap rectbmp = (RectBitmap)bmp;
            Bitmap bm = rectbmp.bmp;
            Rectangle rec = rectbmp.rect;
            int blkORwht = rectbmp.blkORwht;
            byte threshold = rectbmp.threshold;
            int BPP = Image.GetPixelFormatSize(bm.PixelFormat) / 8;
            unsafe
            {
                byte gray, R, G, B;
                BitmapData data = bm.LockBits(new Rectangle(0, 0, rec.Width, rec.Height), ImageLockMode.ReadOnly, bm.PixelFormat);
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * rec.Width;
                for (int j = rec.Y; j < rec.Height + rec.Y; j++)
                {
                    for (int i = rec.X; i < rec.Width + rec.X; i++)
                    {
                        if (exist1 || exist2 || exist3 || exist4 || exist5 || exist6 || exist7 || exist8 || exist9)
                        {
                            exist = true;
                            break;
                        }
                        R = p[2];
                        G = p[1];
                        B = p[0];
                        if (R != G || G != B)
                        {
                            //尚未是灰度化，先计算灰度
                            // gray = 0.3*R + 0.59*G + 0.11*B
                            gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                            //再提取关键点
                            if (blkORwht == 0)//黑边
                            {
                                if (gray < threshold)
                                {
                                    exist = true;
                                    break;
                                }
                            }
                            else//白边
                            {
                                if (gray > threshold)
                                {
                                    exist = true;
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
                                    exist = true;
                                    break;
                                }
                            }
                            else//白边
                            {
                                if (B > threshold)
                                {
                                    exist = true;
                                    break;
                                }
                            }
                        }
                        p += BPP;
                    }
                    if (exist)
                        break;
                    p += offset;
                }
                bm.UnlockBits(data);
            }
        }
        /// <summary>
        /// 检测指定矩形范围内是否存在关键点
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="exist"></param>
        public static void CheckSomeRect(RectBitmap bmp, out bool exist)
        {
            exist = false;
            RectBitmap rectbmp = (RectBitmap)bmp;
            Bitmap bm = rectbmp.bmp;
            Rectangle rec = rectbmp.rect;
            int blkORwht = rectbmp.blkORwht;
            byte threshold = rectbmp.threshold;
            int BPP = Image.GetPixelFormatSize(bm.PixelFormat) / 8;
            unsafe
            {
                byte gray, R, G, B;
                BitmapData data = bm.LockBits(new Rectangle(0, 0, rec.Width, rec.Height), ImageLockMode.ReadOnly, bm.PixelFormat);
                byte* p = (byte*)data.Scan0;
                int stride = data.Stride;
                int offset = stride - BPP * rec.Width;
                for (int j = rec.Y; j < rec.Height + rec.Y; j++)
                {
                    for (int i = rec.X; i < rec.Width + rec.X; i++)
                    {
                        if (existRect1 || existRect2 || existRect3 || existRect4 || existRect5 || existRect6 || existRect7 || existRect8 || existRect9)
                        {
                            exist = true;
                            break;
                        }
                        R = p[2];
                        G = p[1];
                        B = p[0];
                        if (R != G || G != B)
                        {
                            //尚未是灰度化，先计算灰度
                            // gray = 0.3*R + 0.59*G + 0.11*B
                            gray = (byte)((19661 * R + 38666 * G + 7209 * B) >> 16);
                            //再提取关键点
                            if (blkORwht == 0)//黑边
                            {
                                if (gray < threshold)
                                {
                                    exist = true;
                                    break;
                                }
                            }
                            else//白边
                            {
                                if (gray > threshold)
                                {
                                    exist = true;
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
                                    exist = true;
                                    break;
                                }
                            }
                            else//白边
                            {
                                if (B > threshold)
                                {
                                    exist = true;
                                    break;
                                }
                            }
                        }
                        p += BPP;
                    }
                    if (exist)
                        break;
                    p += offset;
                }
                bm.UnlockBits(data);
            }
        }
    }
}
