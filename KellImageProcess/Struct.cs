using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Drawing.Drawing2D;

namespace KellImageProcess
{
    delegate void AsyncMethodCallerLinkedWithPath(RectBitmapWithPath bitmapWithPath, bool isCZ, int direction, out bool linked);
    delegate void AsyncMethodCaller(RectBitmap callParam, out bool exist);
    delegate void AsyncMethodCallerWithPath(RectBitmapWithPath callParam, out bool exist);
    delegate void AsyncMethodCallerWithRegion(RectBitmapWithRegion callParam, out bool exist);

    /// <summary>
    /// 相邻点集结构
    /// </summary>
    public struct LinkedPointList
    {
        /// <summary>
        /// 是否连通
        /// </summary>
        public bool Linked;
        /// <summary>
        /// 点集
        /// </summary>
        public List<Point> PointList;
    }
    /// <summary>
    /// 直方图类型
    /// </summary>
    public enum HistogramType
    {
        /// <summary>
        /// 灰度直方图
        /// </summary>
        Gray,
        /// <summary>
        /// 红色通道直方图
        /// </summary>
        Red,
        /// <summary>
        /// 绿色通道直方图
        /// </summary>
        Green,
        /// <summary>
        /// 蓝色通道直方图
        /// </summary>
        Blue
    }
    /// <summary>
    /// 指定区域的位图结构
    /// </summary>
    public struct RectBitmap
    {
        /// <summary>
        /// 整体位图
        /// </summary>
        public Bitmap bmp;
        /// <summary>
        /// 指定区域
        /// </summary>
        public Rectangle rect;
        /// <summary>
        /// 0为黑检测，其他为白检测
        /// </summary>
        public int blkORwht;
        /// <summary>
        /// 灰度阀值
        /// </summary>
        public byte threshold;
    }
    /// <summary>
    /// 指定路径的位图结构
    /// </summary>
    public struct RectBitmapWithPath
    {
        /// <summary>
        /// 整体位图
        /// </summary>
        public Bitmap bmp;
        /// <summary>
        /// 指定区域
        /// </summary>
        public Rectangle rect;
        /// <summary>
        /// 0为黑检测，其他为白检测
        /// </summary>
        public int blkORwht;
        /// <summary>
        /// 灰度阀值
        /// </summary>
        public byte threshold;
        /// <summary>
        /// 外部大位图的路径
        /// </summary>
        public GraphicsPath path;
    }
    /// <summary>
    /// 指定路径的位图结构
    /// </summary>
    public struct RectBitmapWithRegion
    {
        /// <summary>
        /// 整体位图
        /// </summary>
        public Bitmap bmp;
        /// <summary>
        /// 指定区域
        /// </summary>
        public Rectangle rect;
        /// <summary>
        /// 0为黑检测，其他为白检测
        /// </summary>
        public int blkORwht;
        /// <summary>
        /// 灰度阀值
        /// </summary>
        public byte threshold;
        /// <summary>
        /// 外部大位图的区域
        /// </summary>
        public Region region;
    }
    /// <summary>
    /// 带位置信息的像素
    /// </summary>
    public struct PointColor
    {
        /// <summary>
        /// 位置
        /// </summary>
        public Point Location;
        /// <summary>
        /// 颜色
        /// </summary>
        public Color Color;
    }
    /// <summary>
    /// 直线
    /// </summary>
    public struct Linepoint
    {
        /// <summary>
        /// 在直线上的一点
        /// </summary>
        public PointF lt1;
        /// <summary>
        /// 在直线上的另一点
        /// </summary>
        public PointF lt2;
    }
}
