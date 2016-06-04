using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Drawing.Imaging;
using System.Drawing.Drawing2D;

namespace KellImageProcess
{
    #region 基础类
    /// <summary>
    /// 3×3 转换矩阵
    /// </summary>
    public class Matrix3x3
    {
        /// <summary>
        /// 每个像素占用4个字节
        /// </summary>
        public const int BPP = 4;
        int topLeft = 0, topMid = 0, topRight = 0;
        int midLeft = 0, center = 1, midRight = 0;
        int bottomLeft = 0, bottomMid = 0, bottomRight = 0;
        int scale = 1;
        int kernelOffset = 0;
        /// <summary>
        /// 获取或设置左上点权值
        /// </summary>
        public int TopLeft
        {
            get
            {
                return topLeft;
            }
            set
            {
                topLeft = value;
            }
        }
        /// <summary>
        /// 获取或设置正上点权值
        /// </summary>
        public int TopMid
        {
            get
            {
                return topMid;
            }
            set
            {
                topMid = value;
            }
        }
        /// <summary>
        /// 获取或设置右上点权值
        /// </summary>
        public int TopRight
        {
            get
            {
                return topRight;
            }
            set
            {
                topRight = value;
            }
        }
        /// <summary>
        /// 获取或设置左点权值
        /// </summary>
        public int MidLeft
        {
            get
            {
                return midLeft;
            }
            set
            {
                midLeft = value;
            }
        }
        /// <summary>
        /// 获取或设置中心点权值
        /// </summary>
        public int Center
        {
            get
            {
                return center;
            }
            set
            {
                center = value;
            }
        }
        /// <summary>
        /// 获取或设置右点权值
        /// </summary>
        public int MidRight
        {
            get
            {
                return midRight;
            }
            set
            {
                midRight = value;
            }
        }
        /// <summary>
        /// 获取或设置左下点权值
        /// </summary>
        public int BottomLeft
        {
            get
            {
                return bottomLeft;
            }
            set
            {
                bottomLeft = value;
            }
        }
        /// <summary>
        /// 获取或设置正下点权值
        /// </summary>
        public int BottomMid
        {
            get
            {
                return bottomMid;
            }
            set
            {
                bottomMid = value;
            }
        }
        /// <summary>
        /// 获取或设置右下点权值
        /// </summary>
        public int BottomRight
        {
            get
            {
                return bottomRight;
            }
            set
            {
                bottomRight = value;
            }
        }
        /// <summary>
        /// 获取或设置缩放比例
        /// </summary>
        public int Scale
        {
            get
            {
                return scale;
            }
            set
            {
                scale = value;
            }
        }
        /// <summary>
        /// 获取或设置偏移量
        /// </summary>
        public int Offset
        {
            get
            {
                return kernelOffset;
            }
            set
            {
                kernelOffset = value;
            }
        }
        /// <summary>
        /// 初始化窗口所有点为同一权值
        /// </summary>
        /// <param name="degree">权值</param>
        public void Init(int degree)
        {
            topLeft = topMid = topRight =
            midLeft = center = midRight =
            bottomLeft = bottomMid = bottomRight = degree;
        } // end of Init
        /// <summary>
        /// 将图像按 3X3 窗口进行卷积转换
        /// </summary>
        /// <param name="srcImage">位图流</param>
        /// <returns></returns>
        public Bitmap Convolute(Bitmap srcImage)
        {
            // 避免被零除
            if (scale == 0) scale = 1;
            int width = srcImage.Width;
            int height = srcImage.Height;
            Bitmap dstImage = (Bitmap)srcImage.Clone();
            BitmapData srcData = srcImage.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, srcImage.PixelFormat);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, dstImage.PixelFormat);
            // 图像实际处理区域
            // 图像最外围一圈不处理
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
                // 移向第一行
                src += stride;
                dst += stride;
                for (int y = rectTop; y < rectBottom; y++)
                {
                    // 移向每行第一列
                    src += BPP;
                    dst += BPP;
                    for (int x = rectLeft; x < rectRight; x++)
                    {
                        // 如果当前像素为透明色，则跳过不处理
                        if (src[3] > 0)
                        {
                            // 处理 B, G, R 三分量
                            for (int i = 0; i < 3; i++)
                            {
                                pixel =
                                  src[i - stride - BPP] * topLeft +
                                  src[i - stride] * topMid +
                                  src[i - stride + BPP] * topRight +
                                  src[i - BPP] * midLeft +
                                  src[i] * center +
                                  src[i + BPP] * midRight +
                                  src[i + stride - BPP] * bottomLeft +
                                  src[i + stride] * bottomMid +
                                  src[i + stride + BPP] * bottomRight;
                                pixel = pixel / scale + kernelOffset;
                                if (pixel < 0) pixel = 0;
                                if (pixel > 255) pixel = 255;
                                dst[i] = (byte)pixel;
                            } // i
                        }
                        // 向后移一像素
                        src += BPP;
                        dst += BPP;
                    } // x
                    // 移向下一行
                    // 这里得注意要多移一列，因最右列不处理
                    src += (offset + BPP);
                    dst += (offset + BPP);
                } // y
            }
            srcImage.UnlockBits(srcData);
            dstImage.UnlockBits(dstData);
            return dstImage;
        } // end of Convolute
    }//end of Matrix3x3
    /// <summary>
    /// 数据统计类
    /// </summary>
    public class Statistics
    {
        private int[] Sequence;
        private int length;
        private int min, max;
        private int minIndex, maxIndex;
        /// <summary>
        /// 获取数组序列
        /// </summary>
        public int[] Value
        {
            get
            {
                return Sequence;
            }
        }
        /// <summary>
        /// 获取经均衡化后的数组序列
        /// </summary>
        public byte[] Equalizer
        {
            get
            {
                // 先计算概率
                double[] Probability = this.Probability;
                // S 为亮度级的定积分，即离散图像的亮度变换函数
                double[] S = new double[256];
                // L 数组用于记录均衡化后的新亮度值
                byte[] L = new byte[256];
                // 进行均衡化处理
                for (int i = 0; i < 256; i++)
                {
                    if (i == 0)
                    {
                        S[0] = Probability[0];
                    }
                    else
                    {
                        S[i] = S[i - 1] + Probability[i];
                    }
                    L[i] = (byte)(255 * S[i] + 0.5);
                } // i
                return L;
            }
        }
        /// <summary>
        /// 获取概率
        /// </summary>
        public double[] Probability
        {
            get
            {
                double total = (double)this.Sum;
                double[] probability = new double[256];
                // 计算各亮度级概率密度
                for (int i = 0; i < 256; i++)
                {
                    probability[i] = Sequence[i] / total;
                } // i
                return probability;
            }
        }
        /// <summary>
        /// 获取序列累加和
        /// </summary>
        public int Sum
        {
            get
            {
                int total = 0;
                for (int i = 0; i < length; i++)
                {
                    total += Sequence[i];
                } // i
                return total;
            }
        }
        /// <summary>
        /// 获取加权平均数
        /// </summary>
        public double Mean
        {
            get
            {
                int mean = 0;
                for (int i = 0; i < length; i++)
                {
                    mean += i * Sequence[i];
                } // i
                return (double)mean / (double)this.Sum;
            }
        }
        /// <summary>
        /// 获取标准偏差
        /// </summary>
        public double StdDev
        {
            get
            {
                double mean = this.Mean;
                int total = this.Sum;
                double stddev = 0;
                for (int i = 0; i < length; i++)
                {
                    double t = (double)i - mean;
                    stddev += t * t * Sequence[i];
                } // i
                return Math.Sqrt(stddev / total);
            }
        }
        /// <summary>
        /// 获取中值
        /// </summary>
        public int Median
        {
            get
            {
                int halfTotal = this.Sum / 2;

                // 查找中值
                int total = 0;
                int median = 0;
                while (total < halfTotal)
                {
                    total += Sequence[median];
                    median++;
                } // while
                return median - 1;
            }
        }
        /// <summary>
        /// 获取最大值
        /// </summary>
        public int Maximum
        {
            get
            {
                return max;
            }
        }
        /// <summary>
        /// 获取最小值
        /// </summary>
        public int Minimum
        {
            get
            {
                return min;
            }
        }
        /// <summary>
        /// 获取最大值索引
        /// </summary>
        public int MaxIndex
        {
            get
            {
                return maxIndex;
            }
        }
        /// <summary>
        /// 获取最小值索引
        /// </summary>
        public int MinIndex
        {
            get
            {
                return minIndex;
            }
        }
        /// <summary>
        /// 建立数据统计资料
        /// </summary>
        /// <param name="sequence">数组序列</param>
        public Statistics(int[] sequence)
        {
            this.Sequence = sequence;
            length = Sequence.Length;

            MaxMin();
        }
        /// <summary>
        /// 求最大值、最小值
        /// </summary>
        private void MaxMin()
        {
            max = min = Sequence[0];
            maxIndex = minIndex = 0;
            // 计算最大、最小值
            for (int i = 1; i < length; i++)
            {
                int t = Sequence[i];
                if (t > max)
                {
                    max = t;
                    maxIndex = i;
                }
                if (t < min)
                {
                    min = t;
                    minIndex = i;
                }
            } // i
        } // end of MaxMin
    } // end of Statistics
    /// <summary>
    /// HSL 色彩空间结构体
    /// </summary>
    public struct HSL
    {
        float h, s, l;
        /// <summary>
        /// 获取或设置色调[0, 360]
        /// </summary>
        public float Hue
        {
            get
            {
                return h;
            }
            set
            {
                h = (float)((int)(value) % 360);
            }
        }
        /// <summary>
        /// 获取或设置饱和度[0, 1]
        /// </summary>
        public float Saturation
        {
            get
            {
                return s;
            }
            set
            {
                s = value;
                if (s < 0.0) s = 0.0f;
                if (s > 1.0) s = 1.0f;
            }
        }
        /// <summary>
        /// 获取或设置亮度[0, 1]
        /// </summary>
        public float Luminance
        {
            get
            {
                return l;
            }
            set
            {
                l = value;
                if (l < 0.0f) l = 0.0f;
                if (l > 1.0f) l = 1.0f;
            }
        }
        /// <summary>
        /// 判断 HSL 结构体是否相等
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj)
        {
            return this == (HSL)obj;
        }
        /// <summary>
        /// 返回哈希代码
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode()
        {
            return base.GetHashCode();
        }
        /// <summary>
        /// 判断 HSL 结构体是否相等
        /// </summary>
        /// <param name="lHsl">HSL 结构体 1</param>
        /// <param name="rHsl">HSL 结构体 2</param>
        /// <returns></returns>
        public static bool operator ==(HSL lHsl, HSL rHsl)
        {
            if ((lHsl.Hue == rHsl.Hue) &&
                (lHsl.Saturation == rHsl.Saturation) &&
                (lHsl.Luminance == rHsl.Luminance))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        /// <summary>
        /// 判断 HSL 结构体是否不相等
        /// </summary>
        /// <param name="lHsl">HSL 结构体 1</param>
        /// <param name="rHsl">HSL 结构体 2</param>
        /// <returns></returns>
        public static bool operator !=(HSL lHsl, HSL rHsl)
        {
            return !(lHsl == rHsl);
        }
        /// <summary>
        /// 根据 (hue, saturation, luminance) 分量建立 PhotoSprite.ColorSpace.HSL 结构体
        /// </summary>
        /// <param name="hue">色调[0, 360]</param>
        /// <param name="saturation">饱和度[0, 1]</param>
        /// <param name="luminance">亮度[0, 1]</param>
        /// <returns></returns>
        public static HSL FromHsl(float hue, float saturation, float luminance)
        {
            HSL hsl = new HSL();
            hsl.Hue = hue;
            hsl.Saturation = saturation;
            hsl.Luminance = luminance;
            return hsl;
        } // end of FromHsl
        /// <summary>
        /// 根据 Color 结构体建立 PhotoSprite.ColorSpace.HSL 结构体
        /// </summary>
        /// <param name="color">RGB 颜色结构体</param>
        /// <returns></returns>
        public static HSL FromColor(Color color)
        {
            HSL hsl = new HSL();
            hsl.Hue = color.GetHue();
            hsl.Saturation = color.GetSaturation();
            hsl.Luminance = color.GetBrightness();

            return hsl;
        } // end of FromColor
        /// <summary>
        /// 根据 (red, green, blue) 颜色分量建立 PhotoSprite.ColorSpace.HSL 结构体
        /// </summary>
        /// <param name="red">red 分量</param>
        /// <param name="green">green 分量</param>
        /// <param name="blue">blue 分量</param>
        /// <returns></returns>
        public static HSL FromRgb(byte red, byte green, byte blue)
        {
            return FromColor(Color.FromArgb(red, green, blue));
        } // end of FromRgb
        /// <summary>
        /// 获取 RGB 颜色值
        /// </summary>
        /// <returns></returns>
        public Color ToRgb()
        {
            double r = 0, g = 0, b = 0;
            double t1, t2;
            double normalisedH = h / 360.0;
            if (l == 0)
            {
                r = g = b = 0;
            }
            else
            {
                if (s == 0)
                {
                    r = g = b = l;
                }
                else
                {
                    t2 = ((l <= 0.5) ? l * (1.0 + s) : l + s - l * s);
                    t1 = 2.0 * l - t2;
                    double[] T3 = new double[] { 
              normalisedH + 1.0 / 3.0, 
              normalisedH, 
              normalisedH - 1.0 / 3.0 };
                    double[] C = new double[3];
                    for (int i = 0; i < 3; i++)
                    {
                        if (T3[i] < 0) T3[i] += 1.0;
                        if (T3[i] > 1) T3[i] -= 1.0;

                        if (6.0 * T3[i] < 1.0)
                            C[i] = t1 + (t2 - t1) * T3[i] * 6.0;
                        else if (2.0 * T3[i] < 1.0)
                            C[i] = t2;
                        else if (3.0 * T3[i] < 2.0)
                            C[i] = (t1 + (t2 - t1) * ((2.0 / 3.0) - T3[i]) * 6.0);
                        else
                            C[i] = t1;
                    } // i
                    r = C[0];
                    g = C[1];
                    b = C[2];
                } // s==0
            } // l==0
            return Color.FromArgb((int)(255 * r), (int)(255 * g), (int)(255 * b));
        } // end of ToRgb
        /// <summary>
        /// 获取 RGB 结构中 red 分量值
        /// </summary>
        /// <returns></returns>
        public byte GetRed()
        {
            return (byte)ToRgb().R;
        }
        /// <summary>
        /// 获取 RGB 结构中 green 分量值
        /// </summary>
        /// <returns></returns>
        public byte GetGreen()
        {
            return (byte)ToRgb().G;
        }
        /// <summary>
        /// 获取 RGB 结构中 blue 分量值
        /// </summary>
        /// <returns></returns>
        public byte GetBlue()
        {
            return (byte)ToRgb().B;
        }
    } // end of class HSL

    /// <summary>
    /// 图像调整类
    /// </summary>
    public class Adjustment
    {
        //public const int BPP = 4;
        /************************************************************
         * 
         * 色彩平衡、亮度、对比度、HSL 调整、Gamma 矫正
         * 
         ************************************************************/
        /// <summary>
        /// 图像色彩平衡
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="red">红色分量[-255, 255]</param>
        /// <param name="green">绿色分量[-255, 255]</param>
        /// <param name="blue">蓝色分量[-255, 255]</param>
        /// <returns></returns>
        public static Bitmap ColorBalance(Bitmap b, int red, int green, int blue)
        {
            if (red < -255) red = -255;
            if (red > 255) red = 255;
            if (green < -255) green = -255;
            if (green > 255) green = 255;
            if (blue < -255) blue = -255;
            if (blue > 255) blue = 255;
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                int pixel;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        pixel = p[2] + red;
                        if (pixel < 0) pixel = 0;
                        if (pixel > 255) pixel = 255;
                        p[2] = (byte)pixel;
                        pixel = p[1] + green;
                        if (pixel < 0) pixel = 0;
                        if (pixel > 255) pixel = 255;
                        p[1] = (byte)pixel;
                        pixel = p[0] + blue;
                        if (pixel < 0) pixel = 0;
                        if (pixel > 255) pixel = 255;
                        p[0] = (byte)pixel;
                        p += BPP;
                    } // x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            return b;
        } // end of ColorBalance
        /// <summary>
        /// 图像亮度调整
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="degree">亮度值[-255, 255]</param>
        /// <returns></returns>
        public static Bitmap Brightness(Bitmap b, int degree)
        {
            if (degree < -255) degree = -255;
            if (degree > 255) degree = 255;
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                int pixel = 0;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        // 处理像素 B, G, R 亮度三分量
                        for (int i = 0; i < 3; i++)
                        {
                            pixel = p[i] + degree;
                            if (pixel < 0) pixel = 0;
                            if (pixel > 255) pixel = 255;
                            p[i] = (byte)pixel;
                        } // i
                        p += BPP;
                    }  // x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            return b;
        } // end of Brightness

        /// <summary>
        /// 图像对比度调整
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="degree">对比度[-100, 100]</param>
        /// <returns></returns>
        public static Bitmap Contrast(Bitmap b, int degree)
        {
            if (degree < -100) degree = -100;
            if (degree > 100) degree = 100;
            double pixel = 0;
            double contrast = (100.0 + degree) / 100.0;
            contrast *= contrast;
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        // 处理指定位置像素的对比度
                        for (int i = 0; i < 3; i++)
                        {
                            pixel = ((p[i] / 255.0 - 0.5) * contrast + 0.5) * 255;
                            if (pixel < 0) pixel = 0;
                            if (pixel > 255) pixel = 255;
                            p[i] = (byte)pixel;
                        } // i
                        p += BPP;
                    } // x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            return b;
        } // end of Contrast

        /// <summary>
        /// 按指定的色调、饱和度、亮度对图像进行调整
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="hue">色调[-180, 180]</param>
        /// <param name="saturation">饱和度[-1, 1]</param>
        /// <param name="luminance">亮度[-1, 1]</param>
        /// <returns></returns>
        public static Bitmap AdjustHsl(Bitmap b, float hue, float saturation, float luminance)
        {
            if (hue < -180.0f) hue = -180.0f;
            if (hue > 180.0f) hue = 180f;
            if (saturation < -1.0f) saturation = -1.0f;
            if (saturation > 1.0f) saturation = 1.0f;
            if (luminance < -1.0f) luminance = -1.0f;
            if (luminance > 1.0f) luminance = 1.0f;
            int width = b.Width;
            int height = b.Height;
            Bitmap dstImage = new Bitmap(width, height, b.PixelFormat);
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        HSL hsl = HSL.FromRgb(p[2], p[1], p[0]);
                        hsl.Hue += hue;
                        hsl.Saturation += saturation;
                        hsl.Luminance += luminance;
                        p[0] = hsl.GetBlue();
                        p[1] = hsl.GetGreen();
                        p[2] = hsl.GetRed();
                        p += BPP;
                    } // x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            return b;
        } // end of AdjustHsl
        /// <summary>
        /// 图像 Gamma 矫正
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="degree">Gamma 矫正量[0.1, 5.0]</param>
        /// <returns></returns>
        public static Bitmap GammaCorrect(Bitmap b, double degree)
        {
            if (degree < 0.1) degree = 0.1;
            if (degree > 5.0) degree = 5.0;
            byte[] Gamma = new byte[256];
            double g = 1 / degree;
            // 建立 Gamma 矫正映射表
            for (int i = 0; i < 256; i++)
            {
                int pixel = (int)((255.0 * Math.Pow(i / 255.0, g)) + 0.5);
                Gamma[i] = (byte)((pixel > 255) ? 255 : pixel);
            } // i
            // 根据 Gamma 矫正映射表对图像色彩进行映射处理
            return Mapping(b, Gamma, ChannelMode.White);
        } // end of GammaCorrect
        /************************************************************
         * 
         * 负像、交错负像、伪彩色
         * 
         ************************************************************/
        /// <summary>
        /// 反转，负像
        /// </summary>
        /// <param name="b">位图流</param>
        /// <returns></returns>
        public static Bitmap Invert(Bitmap b)
        {
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        // 与255进行异或运算，相当于按位非
                        p[2] ^= 0xFF;
                        p[1] ^= 0xFF;
                        p[0] ^= 0xFF;
                        p += BPP;
                    } // x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            return b;
        } // end of Invert

        /// <summary>
        /// 交叉反转
        /// </summary>
        /// <param name="b">位图流</param>
        /// <returns></returns>
        public static Bitmap Interleaving(Bitmap b)
        {
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        // 将指定位置的像素颜色反转
                        if ((x + y) % 2 == 0)
                        {
                            p[0] ^= 0xFF;
                            p[1] ^= 0xFF;
                            p[2] ^= 0xFF;
                        }
                        p += BPP;
                    } // x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            return b;
        } // end of Interleaving
        /// <summary>
        /// 按函数曲线形式映射伪彩色
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="graied">已经灰度化</param>
        /// <returns></returns>
        public static Bitmap PseudoColor(Bitmap b, bool graied)
        {
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                int R, G, B, gray;
                R = G = B = gray = 0;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        if (graied)
                            gray = p[0];
                        else
                            gray = (byte)((19661 * p[2] + 38666 * p[1] + 7209 * p[0]) >> 16);
                        // 伪彩色处理
                        switch (gray / 64)
                        {
                            case 0:
                                R = 0;
                                G = 4 * gray;
                                B = 255;
                                break;
                            case 1:
                                R = 0;
                                G = 255;
                                B = 511 - 4 * gray;
                                break;
                            case 2:
                                R = 4 * gray - 511;
                                G = 255;
                                B = 0;
                                break;
                            case 3:
                                R = 255;
                                G = 1023 - 4 * gray;
                                B = 0;
                                break;
                        }
                        p[0] = (byte)B;
                        p[1] = (byte)G;
                        p[2] = (byte)R;
                        p += BPP;
                    } // x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            return b;
        } // end of PseudoColor
        /// <summary>
        /// 按色彩表形式映射伪彩色
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="colorTable">色彩映射表</param>
        /// <param name="graied">已经灰度化</param>
        /// <returns></returns>
        public static Bitmap PseudoColor(Bitmap b, Color[] colorTable, bool graied)
        {
            int width = b.Width;
            int height = b.Height;
            // 颜色对照表
            int lenColorTable = colorTable.Length;
            uint[] ColorTable = new uint[lenColorTable];
            // 将颜色转换为数字
            for (int i = 0; i < lenColorTable; i++)
            {
                ColorTable[i] = (uint)((colorTable[i].A << 24) |
                  (colorTable[i].R << 16) |
                  (colorTable[i].G << 8) |
                  (colorTable[i].B << 0));
            } // i
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                int gray = 0;
                uint color = 0;
                // 避免在下面的循环中，出现色彩表索引越界
                lenColorTable--;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        if (graied)
                            gray = p[0];
                        else
                            gray = (byte)((19661 * p[2] + 38666 * p[1] + 7209 * p[0]) >> 16);
                        // 转化灰度级，并映射到色彩表
                        color = ColorTable[lenColorTable * gray / 255];
                        p[3] = (byte)(color >> 24); // A
                        p[2] = (byte)(color >> 16); // R
                        p[1] = (byte)(color >> 8);  // G
                        p[0] = (byte)(color);       // B
                        p += BPP;
                    } // x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            return b;
        } // end of PseudoColor
        /************************************************************
         * 
         * 轮换通道、提取通道、过滤通道
         * 
         ************************************************************/
        /// <summary>
        /// 轮换通道
        /// </summary>
        /// <param name="b">位图流</param>
        /// <returns></returns>
        public static Bitmap RotateChannel(Bitmap b)
        {
            int width = b.Width;
            int height = b.Height;
            Bitmap dstImage = (Bitmap)b.Clone();
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData srcData = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, b.PixelFormat);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, b.PixelFormat);
            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;
                int offset = srcData.Stride - width * BPP;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        dst[2] = src[1]; // R <- G
                        dst[1] = src[0]; // G <- B
                        dst[0] = src[2]; // B <- R
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
        } // end of RotateChannel
        /// <summary>
        /// 通道模式
        /// </summary>
        public enum ChannelMode : int
        {
            /// <summary>
            /// 蓝色通道
            /// </summary>
            Blue = 1,
            /// <summary>
            /// 绿色通道
            /// </summary>
            Green = 2,
            /// <summary>
            /// 红色通道
            /// </summary>
            Red = 4,
            /// <summary>
            /// Alpha 通道
            /// </summary>
            Alpha = 8,
            /// <summary>
            /// 青色 = 绿色 + 蓝色
            /// </summary>
            Cyan = 3,
            /// <summary>
            /// 品红 = 红色 + 蓝色
            /// </summary>
            Megenta = 5,
            /// <summary>
            /// 黄色 = 红色 + 绿色
            /// </summary>
            Yellow = 6,
            /// <summary>
            /// 白色 = 红色 + 绿色 + 蓝色
            /// </summary>
            White = 7
        }
        /// <summary>
        /// 提取通道
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="channelMode">通道模式[A, R, G, B]</param>
        /// <returns></returns>
        public static Bitmap ExtractChannel(Bitmap b, ChannelMode channelMode)
        {
            int channel = (int)Math.Log((double)channelMode, 2.0);
            int width = b.Width;
            int height = b.Height;
            Bitmap dstImage = (Bitmap)b.Clone();
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData srcData = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, b.PixelFormat);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.WriteOnly, b.PixelFormat);
            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;
                int offset = srcData.Stride - width * BPP;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        dst[0] = dst[1] = dst[2] = src[channel];
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
        } // end of ExtractChannel
        /// <summary>
        /// 过滤通道
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="channelMode">通道模式</param>
        /// <returns></returns>
        public static Bitmap FilterChannel(Bitmap b, ChannelMode channelMode)
        {
            int channel = (int)channelMode;
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                for (int i = 0; i < 3; i++)
                {
                    if (((int)Math.Pow(2, i) & channel) > 0)
                        continue;
                    p = (byte*)data.Scan0;
                    for (int y = 0; y < height; y++)
                    {
                        for (int x = 0; x < width; x++)
                        {
                            p[i] = 0;
                            p += BPP;
                        } // x
                        p += offset;
                    } // y
                } // i
            }
            b.UnlockBits(data);
            return b;
        } // end of FilterChannel
        /************************************************************
         * 
         * 映射
         * 
         ************************************************************/
        /// <summary>
        /// 图像色彩映射
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="Map">映射表</param>
        /// <param name="channelMode">通道模式</param>
        /// <returns></returns>
        public static Bitmap Mapping(Bitmap b, byte[] Map, ChannelMode channelMode)
        {
            int channel = (int)channelMode;
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                for (int i = 0; i < 3; i++)
                {
                    if (((int)Math.Pow(2, i) & channel) == 0)
                        continue;
                    p = (byte*)data.Scan0;
                    for (int y = 0; y < height; y++)
                    {
                        for (int x = 0; x < width; x++)
                        {
                            p[i] = Map[p[i]];
                            p += BPP;
                        } // x
                        p += offset;
                    } // y
                } // i
            }
            b.UnlockBits(data);
            return b;
        } // end of Mapping
    }

    /// <summary>
    /// 直方图类
    /// </summary>
    public class Histogram
    {
        /// <summary>
        /// 每个像素占用4个字节
        /// </summary>
        public const int BPP = 4;
        /// <summary>
        /// 支持绘制的色彩模式
        /// </summary>
        public enum ColorMode
        {
            /// <summary>
            /// 红色
            /// </summary>
            Red,
            /// <summary>
            /// 绿色
            /// </summary>
            Green,
            /// <summary>
            /// 蓝色
            /// </summary>
            Blue,
            /// <summary>
            /// 亮度
            /// </summary>
            Brightness,
            /// <summary>
            /// 灰度
            /// </summary>
            Gray
        }
        private Statistics red;
        private Statistics green;
        private Statistics blue;
        private Statistics bright;
        private Statistics gray;
        private Bitmap b;
        private int width = 0;
        private int height = 0;
        /// <summary>
        /// 获取红色分量统计资料
        /// </summary>
        public Statistics Red
        {
            get
            {
                return red;
            }
        }
        /// <summary>
        /// 获取绿色分量统计资料
        /// </summary>
        public Statistics Green
        {
            get
            {
                return green;
            }
        }
        /// <summary>
        /// 获取蓝色分量统计资料
        /// </summary>
        public Statistics Blue
        {
            get
            {
                return blue;
            }
        }
        /// <summary>
        /// 获取亮度统计资料
        /// </summary>
        public Statistics Bright
        {
            get
            {
                return bright;
            }
        }
        /// <summary>
        /// 获取灰度统计资料
        /// </summary>
        public Statistics Gray
        {
            get
            {
                return gray;
            }
        }
        /// <summary>
        /// 建立图像直方图
        /// </summary>
        /// <param name="b">位图流</param>
        public Histogram(Bitmap b)
        {
            if (b != null)
            {
                this.b = b;
                width = b.Width;
                height = b.Height;
                CountRgb();
                CountBright();
                CountGray();
            }
        } // end of Histogram
        /// <summary>
        /// 统计色彩 R、G、B 三分量频率
        /// </summary>
        private void CountRgb()
        {
            int[] Red = new int[256];
            int[] Green = new int[256];
            int[] Blue = new int[256];
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        Red[p[2]]++;
                        Green[p[1]]++;
                        Blue[p[0]]++;
                        p += BPP;
                    } // x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            // 生成 R、G、B 三分量统计类
            this.red = new Statistics(Red);
            this.green = new Statistics(Green);
            this.blue = new Statistics(Blue);
        } // end of CountRgb
        /// <summary>
        /// 统计亮度频率
        /// </summary>
        private void CountBright()
        {
            int[] Bright = new int[256];
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        Bright[(int)(Color.FromArgb(p[2], p[1], p[0]).GetBrightness() * 255)]++;
                        p += BPP;
                    } // x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            // 生成 Bright 统计类
            this.bright = new Statistics(Bright);
        } // end of CountBright
        /// <summary>
        /// 统计亮度频率
        /// </summary>
        private void CountGray()
        {
            int[] Gray = new int[256];
            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadOnly, b.PixelFormat);
            unsafe
            {
                byte* p = (byte*)data.Scan0;
                int offset = data.Stride - width * BPP;
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        Gray[(19661 * p[2] + 38666 * p[1] + 7209 * p[0]) >> 16]++;
                        p += BPP;
                    } // x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            // 生成 Bright 统计类
            this.gray = new Statistics(Gray);
        } // end of CountGray
        /// <summary>
        /// 绘制直方图
        /// </summary>
        /// <param name="diagramHeight">图表高度</param>
        /// <param name="viewByLog">按 Log 函数绘制直方图</param>
        /// <param name="colorMode">色彩模式</param>
        /// <returns></returns>
        public Bitmap DrawDiagram(int diagramHeight, bool viewByLog, ColorMode colorMode)
        {
            // 获取亮度概率
            double[] Probability = this.Red.Probability;
            // 亮度概率最大值
            double maxProbability = Probability[this.Red.MaxIndex];
            // 用于绘制直方图的颜色
            Color color = Color.Red;
            switch (colorMode)
            {
                case ColorMode.Red:
                    Probability = this.Red.Probability;
                    maxProbability = Probability[this.Red.MaxIndex];
                    color = Color.Red;
                    break;
                case ColorMode.Green:
                    Probability = this.Green.Probability;
                    maxProbability = Probability[this.Green.MaxIndex];
                    color = Color.Green;
                    break;
                case ColorMode.Blue:
                    Probability = this.Blue.Probability;
                    maxProbability = Probability[this.Blue.MaxIndex];
                    color = Color.Blue;
                    break;
                case ColorMode.Brightness:
                    Probability = this.Bright.Probability;
                    maxProbability = Probability[this.Bright.MaxIndex];
                    color = Color.Black;
                    break;
                case ColorMode.Gray:
                    Probability = this.Gray.Probability;
                    maxProbability = Probability[this.Gray.MaxIndex];
                    color = Color.Gray;
                    break;
            } // switch
            Pen pen = new Pen(color, 1);
            Bitmap dstImage = new Bitmap(256, diagramHeight);
            System.Drawing.Graphics g = System.Drawing.Graphics.FromImage(dstImage);
            int y = 0;
            // 绘制直方图
            for (int i = 0; i < 256; i++)
            {
                // 当前亮度级与最大亮度级的比率
                double percent = Probability[i] / maxProbability;
                if (viewByLog)
                {
                    // 因为不同色彩亮度级概率密度分布不同，分布曲线可能平缓，也可能陡峭，
                    // 鉴于这种差异，我们对色彩概率进行对数处理，以拉开分布情况，便于查阅
                    y = (int)(diagramHeight * (1 - Math.Log(100 * percent + 1, 101)));
                }
                else
                {
                    y = (int)(diagramHeight * (1 - percent));
                }
                g.DrawLine(pen, i, y, i, diagramHeight);
            } // i
            g.Save();
            g.Dispose();
            return dstImage;
        } // end of DrawDiagram
        /// <summary>
        /// 绘制经过直方图均衡化处理的图像
        /// </summary>
        /// <returns></returns>
        public Bitmap Equalizer()
        {
            // 图像各通道均衡化后映射表
            byte[] RedMap = this.Red.Equalizer;
            byte[] GreenMap = this.Green.Equalizer;
            byte[] BlueMap = this.Blue.Equalizer;
            Bitmap dstImage = (Bitmap)b.Clone();
            // 通道映射
            dstImage = Adjustment.Mapping(dstImage, RedMap, Adjustment.ChannelMode.Red);
            dstImage = Adjustment.Mapping(dstImage, GreenMap, Adjustment.ChannelMode.Green);
            dstImage = Adjustment.Mapping(dstImage, BlueMap, Adjustment.ChannelMode.Blue);

            return dstImage;
        } // end of Equalizer
    }//end of Histogram

    /// <summary>
    /// 逻辑运算类
    /// </summary>
    public class Logic
    {
        /// <summary>
        /// 逻辑运算方法
        /// </summary>
        public enum LogicMethod
        {
            /// <summary>
            /// 与运算
            /// </summary>
            And,

            /// <summary>
            /// 或运算
            /// </summary>
            Or,

            /// <summary>
            /// 异或运算
            /// </summary>
            Xor
        }

        /// <summary>
        /// 获取或设置背景区域
        /// </summary>
        public static Region BackgroundRegion
        {
            get
            {
                return bgRegion;
            }
            set
            {
                bgRegion = value;
            }
        }
        private static Region bgRegion = new Region();

        /// <summary>
        /// 获取或设置前景区域
        /// </summary>
        public static Region ForegroundRegion
        {
            get
            {
                return fgRegion;
            }
            set
            {
                fgRegion = value;
            }
        }
        private static Region fgRegion = new Region();


        /// <summary>
        /// 图像逻辑运算
        /// </summary>
        /// <param name="bgImage">二值背景</param>
        /// <param name="fgImage">二值前景</param>
        /// <param name="logicMethod">逻辑运算方法</param>
        /// <returns></returns>
        public static Bitmap LogicOperate(Bitmap bgImage, Bitmap fgImage, LogicMethod logicMethod)
        {
            Bitmap dstImage = (Bitmap)bgImage.Clone();
            Graphics g = System.Drawing.Graphics.FromImage(dstImage);

            // 计算有效区域
            Region validRegion = bgRegion;
            validRegion.Intersect(fgRegion);
            RectangleF validRect = validRegion.GetBounds(g);
            RectangleF fgRect = fgRegion.GetBounds(g);

            RegionClip bgRegionClip = new RegionClip(validRegion);
            Bitmap background = RegionClip.Hold((Bitmap)bgImage.Clone());

            RegionClip fgRegionClip = new RegionClip(validRegion);
            validRegion.Translate(-fgRect.X, -fgRect.Y);
            Bitmap foreground = RegionClip.Hold((Bitmap)fgImage.Clone());
            validRegion.Translate(fgRect.X, fgRect.Y);

            // 先将原始二值图转化为二维数组
            byte[,] bgGray = AccessPixel.Image2Array(background);
            byte[,] fgGray = AccessPixel.Image2Array(foreground);

            // 进行逻辑运算处理后的灰度二维数组
            byte[,] dstGray = null;

            // 进行逻辑运算
            switch (logicMethod)
            {
                case LogicMethod.And:
                    dstGray = LogicAnd(bgGray, fgGray);
                    break;

                case LogicMethod.Or:
                    dstGray = LogicOr(bgGray, fgGray);
                    break;

                case LogicMethod.Xor:
                    dstGray = LogicXor(bgGray, fgGray);
                    break;
            }

            // 将二值数组转化为灰度图
            Bitmap validImage = AccessPixel.Array2Image(dstGray, fgImage.PixelFormat);

            g.DrawImage(validImage, validRect,
              new Rectangle(0, 0, (int)validRect.Width, (int)validRect.Height), GraphicsUnit.Pixel);

            bgImage.Dispose();
            fgImage.Dispose();

            return dstImage;
        } // end of LogicOperate


        /// <summary>
        /// 图像逻辑与运算
        /// </summary>
        /// <param name="bg">背景二值化数组</param>
        /// <param name="fg">前景二值化数组</param>
        /// <returns></returns>
        private static byte[,] LogicAnd(byte[,] bg, byte[,] fg)
        {
            int width = bg.GetLength(0);
            int height = bg.GetLength(1);

            // 初始化目标数组为 255，即白色
            byte[,] dst = AccessPixel.InitArray(width, height, 255);

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    dst[x, y] = (byte)(bg[x, y] & fg[x, y]);
                } // x
            } // y

            return dst;
        } // end of LogicAnd


        /// <summary>
        /// 图像逻辑或运算
        /// </summary>
        /// <param name="bg">背景二值化数组</param>
        /// <param name="fg">前景二值化数组</param>
        /// <returns></returns>
        private static byte[,] LogicOr(byte[,] bg, byte[,] fg)
        {
            int width = bg.GetLength(0);
            int height = bg.GetLength(1);

            // 初始化目标数组为 255，即白色
            byte[,] dst = AccessPixel.InitArray(width, height, 255);

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    dst[x, y] = (byte)(bg[x, y] | fg[x, y]);
                } // x
            } // y

            return dst;
        } // end of LogicOr


        /// <summary>
        /// 逻辑非运算
        /// </summary>
        /// <param name="b">二值位图流</param>
        /// <returns></returns>
        public static Bitmap LogicNot(Bitmap b)
        {
            // 对二值图像直接进行负像处理
            b = Adjustment.Invert(b);

            return b;
        } // end of LogicNot


        /// <summary>
        /// 图像逻辑异或运算
        /// </summary>
        /// <param name="bg">背景二值化数组</param>
        /// <param name="fg">前景二值化数组</param>
        /// <returns></returns>
        private static byte[,] LogicXor(byte[,] bg, byte[,] fg)
        {
            int width = bg.GetLength(0);
            int height = bg.GetLength(1);

            // 初始化目标数组为 255，即白色
            byte[,] dst = AccessPixel.InitArray(width, height, 255);

            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    dst[x, y] = (byte)(bg[x, y] ^ fg[x, y]);
                } // x
            } // y

            return dst;
        } // end of LogicXor
    }// end of Logic

    /// <summary>
    /// 区域修整类
    /// </summary>
    public class RegionClip
    {
        private static Region region = null;

        /// <summary>
        /// 获取或设置区域
        /// </summary>
        public static Region SelectedRegion
        {
            get
            {
                return region;
            }
            set
            {
                region = value;
            }
        }


        /// <summary>
        /// 初始化区域修整类
        /// </summary>
        /// <param name="region">修整区域</param>
        public RegionClip(Region region)
        {
            RegionClip.region = region;
        }


        /// <summary>
        /// 去除掉图像中指定区域
        /// </summary>
        /// <param name="b">位图</param>
        /// <returns></returns>
        public static Bitmap Remove(Bitmap b)
        {
            int width = b.Width;
            int height = b.Height;

            Bitmap dstImage = new Bitmap(width, height);
            Graphics g = System.Drawing.Graphics.FromImage(dstImage);
            Region all = new Region(new Rectangle(0, 0, width, height));
            Region validRegion = (Region)region.Clone();
            validRegion.Complement(all);
            g.SetClip(validRegion, System.Drawing.Drawing2D.CombineMode.Replace);
            g.DrawImage(b, new Rectangle(0, 0, width, height), new Rectangle(0, 0, width, height), GraphicsUnit.Pixel);

            return dstImage;
        } // end of Remove


        /// <summary>
        /// 保留下图像中指定区域
        /// </summary>
        /// <param name="b">位图</param>
        /// <returns></returns>
        public static Bitmap Hold(Bitmap b)
        {
            return ImageTransform.Crop(b, (Region)region.Clone());
        } // end of Hold


        /// <summary>
        /// 将背景中指定区域用前景的相应区域替换
        /// </summary>
        /// <param name="bgImage">背景</param>
        /// <param name="fgImage">前景</param>
        /// <returns></returns>
        public static Bitmap Replace(Bitmap bgImage, Bitmap fgImage)
        {
            int width = bgImage.Width;
            int height = bgImage.Height;

            Graphics g = System.Drawing.Graphics.FromImage(bgImage);
            g.SetClip(region, System.Drawing.Drawing2D.CombineMode.Replace);
            g.DrawImage(fgImage, new Rectangle(0, 0, width, height), new Rectangle(0, 0, width, height), GraphicsUnit.Pixel);

            fgImage.Dispose();
            return bgImage;
        } // end of Replace
    } // end of RegionClip

    /// <summary>
    /// 图像变换类
    /// </summary>
    public class ImageTransform
    {
        /// <summary>
        /// 未知区域设置模式
        /// </summary>
        public enum AreasMode
        {
            /// <summary>
            /// 透明
            /// </summary>
            Transparent,

            /// <summary>
            /// 重复边缘像素
            /// </summary>
            RepeatEdgePixels,

            /// <summary>
            /// 四周环绕
            /// </summary>
            WrapAround
        }

        /// <summary>
        /// 修整模式
        /// </summary>
        public enum TrimMode : int
        {
            /// <summary>
            /// 不修整任何边
            /// </summary>
            None = 0,

            /// <summary>
            /// 修整上边
            /// </summary>
            Top = 1,

            /// <summary>
            /// 修整下边
            /// </summary>
            Bottom = 2,

            /// <summary>
            /// 修整左边
            /// </summary>
            Left = 4,

            /// <summary>
            /// 修整右边
            /// </summary>
            Right = 8
        }


        /************************************************************
         * 
         * 平移、尺寸、裁剪
         * 旋转、翻转、转置、倾斜、修整
         * 
         ************************************************************/

        /// <summary>
        /// 图像平移
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="horizontal">水平偏移量</param>
        /// <param name="vertical">垂直偏移量</param>
        /// <param name="areaMode">平移后图像留下的未知区域设置方式</param>
        /// <returns></returns>
        public static Bitmap Translate(Bitmap b, int horizontal, int vertical, AreasMode areaMode)
        {
            int width = b.Width;
            int height = b.Height;

            horizontal %= width;
            vertical %= height;

            Bitmap dstImage = new Bitmap(width, height);
            System.Drawing.Graphics g = System.Drawing.Graphics.FromImage(dstImage);

            // 图像偏移部分
            Point diagonal = new Point(horizontal, vertical);
            g.DrawImage(b, diagonal);

            switch (areaMode)
            {
                case AreasMode.WrapAround:
                    if (horizontal < 0)
                        diagonal.X += width;
                    else
                        diagonal.X -= width;
                    // 水平移出部分
                    g.DrawImage(b, diagonal);

                    if (vertical < 0)
                        diagonal.Y += height;
                    else
                        diagonal.Y -= height;
                    // 对角移出部分
                    g.DrawImage(b, diagonal);

                    // 垂直移出部分
                    diagonal.X = horizontal;
                    g.DrawImage(b, diagonal);
                    break;

                case AreasMode.RepeatEdgePixels:
                    // 水平移出部分
                    diagonal.X += (horizontal < 0 ? width : -width);
                    g.DrawImage(b,
                      new Rectangle(diagonal.X, diagonal.Y, width, height),
                      new Rectangle(horizontal < 0 ? width - 1 : 0, 0, 1, height),
                      GraphicsUnit.Pixel
                      );

                    // 对角移出部分
                    diagonal.Y += (vertical < 0 ? height : -height);
                    g.DrawImage(b,
                      new Rectangle(diagonal.X, diagonal.Y, width, height),
                      new Rectangle(horizontal < 0 ? width - 1 : 0, vertical < 0 ? height - 1 : 0, 1, 1),
                      GraphicsUnit.Pixel
                      );

                    // 垂直移出部分
                    diagonal.X = horizontal;
                    g.DrawImage(b,
                     new Rectangle(diagonal.X, diagonal.Y, width, height),
                     new Rectangle(0, vertical < 0 ? height - 1 : 0, width, 1),
                     GraphicsUnit.Pixel
                     );
                    break;

                case AreasMode.Transparent:
                    break;
            } // switch


            g.Save();
            g.Dispose();

            return dstImage;
        } // end of Translate


        /// <summary>
        /// 图像尺寸调节
        /// </summary>
        /// <param name="b">原始图像</param>
        /// <param name="dstWidth">目标宽度</param>
        /// <param name="dstHeight">目标高度</param>
        /// <returns></returns>
        public static Bitmap Resize(Bitmap b, int dstWidth, int dstHeight)
        {
            Bitmap dstImage = new Bitmap(dstWidth, dstHeight);
            System.Drawing.Graphics g = System.Drawing.Graphics.FromImage(dstImage);

            // 设置插值模式
            g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.Bilinear;

            // 设置平滑模式
            g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.HighQuality;

            g.DrawImage(b,
              new System.Drawing.Rectangle(0, 0, dstImage.Width, dstImage.Height),
              new System.Drawing.Rectangle(0, 0, b.Width, b.Height),
              System.Drawing.GraphicsUnit.Pixel);
            g.Save();
            g.Dispose();

            return dstImage;
        } // end of Resize


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

            Bitmap dstImage = new Bitmap(width, height);
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
        /// 以逆时针方向为正方向对图像进行旋转
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="angle">旋转角度[0, 360]</param>
        /// <returns></returns>
        public static Bitmap Rotate(Bitmap b, int angle)
        {
            angle = angle % 360;

            // 弧度转化
            double radian = angle * Math.PI / 180.0;
            double cos = Math.Cos(radian);
            double sin = Math.Sin(radian);

            // 原图宽高
            int w = b.Width;
            int h = b.Height;

            // 新图的宽高
            int W = (int)(Math.Max(Math.Abs(w * cos - h * sin), Math.Abs(w * cos + h * sin)));
            int H = (int)(Math.Max(Math.Abs(w * sin - h * cos), Math.Abs(w * sin + h * cos)));

            // 目标位图，即旋转后的图
            Bitmap dstImage = new Bitmap(W, H);
            System.Drawing.Graphics g = System.Drawing.Graphics.FromImage(dstImage);
            g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.Bilinear;
            g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.HighQuality;

            // 偏移量
            Point offset = new Point((W - w) / 2, (H - h) / 2);

            // 构造图像显示区域：让图像的中心点与窗口的中心点一致
            Rectangle rect = new Rectangle(offset.X, offset.Y, w, h);
            Point center = new Point(rect.X + rect.Width / 2, rect.Y + rect.Height / 2);

            // 以图像的中心点旋转
            g.TranslateTransform(center.X, center.Y);
            g.RotateTransform(360 - angle);

            // 恢复图像在水平和垂直方向的平移
            g.TranslateTransform(-center.X, -center.Y);

            // 绘制旋转后的结果图
            g.DrawImage(b, rect);

            // 重置绘图的所有变换
            g.ResetTransform();

            g.Save();

            return dstImage;
        } // end of Rotate


        /// <summary>
        /// 对图像进行翻转变换，即镜像
        /// </summary>
        /// <param name="b">原始图像</param>
        /// <param name="isHorz">是否按水平方向进行翻转</param>
        /// <returns></returns>
        public static Bitmap Flip(Bitmap b, bool isHorz)
        {
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            Bitmap dstImage = new Bitmap(width, height);

            BitmapData srcData = b.LockBits(new Rectangle(0, 0, width, height),
              ImageLockMode.ReadOnly, b.PixelFormat);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, width, height),
              ImageLockMode.WriteOnly, dstImage.PixelFormat);

            int stride = srcData.Stride;
            System.IntPtr srcScan0 = srcData.Scan0;
            System.IntPtr dstScan0 = dstData.Scan0;
            int offset = stride - width * BPP;

            unsafe
            {
                byte* src = (byte*)srcScan0;
                byte* dst = (byte*)dstScan0;

                if (isHorz)
                {
                    // 水平翻转
                    for (int y = 0; y < height; y++)
                    {
                        for (int x = 0; x < width; x++)
                        {
                            dst = (byte*)dstScan0 + (y * stride) + ((width - x - 1) * BPP);

                            dst[0] = src[0]; // B
                            dst[1] = src[1]; // G
                            dst[2] = src[2]; // R
                            //dst[3] = src[3]; // A

                            src += BPP;
                        } // x

                        src += offset;
                    } // y
                }
                else
                {
                    // 垂直翻转
                    for (int y = 0; y < height; y++)
                    {
                        for (int x = 0; x < width; x++)
                        {
                            dst = (byte*)dstScan0 + ((height - y - 1) * stride) + (x * BPP);

                            dst[0] = src[0]; // B
                            dst[1] = src[1]; // G
                            dst[2] = src[2]; // R
                            //dst[3] = src[3]; // A

                            src += BPP;
                        } // x

                        src += offset;
                    } // y
                } // isHorz
            }

            b.UnlockBits(srcData);
            dstImage.UnlockBits(dstData);

            return dstImage;
        } // end of Flip


        /// <summary>
        /// 对图像进行转置变换
        /// </summary>
        /// <param name="b">原始图像</param>
        /// <returns></returns>
        public static Bitmap Transpose(Bitmap b)
        {
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            Bitmap dstImage = new Bitmap(height, width);

            BitmapData srcData = b.LockBits(new Rectangle(0, 0, width, height),
              ImageLockMode.ReadOnly, b.PixelFormat);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, height, width),
              ImageLockMode.WriteOnly, b.PixelFormat);

            int srcStride = srcData.Stride;
            int dstStride = dstData.Stride;
            System.IntPtr srcScan0 = srcData.Scan0;
            System.IntPtr dstScan0 = dstData.Scan0;
            int offset = srcData.Stride - width * BPP;

            unsafe
            {
                byte* src = (byte*)srcScan0;
                byte* dst = (byte*)dstScan0;

                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        dst = (byte*)dstScan0 + (x * dstStride) + (y * BPP);

                        dst[0] = src[0];
                        dst[1] = src[1];
                        dst[2] = src[2];
                        //dst[3] = src[3];

                        src += BPP;
                    } // x

                    src += offset;
                } // y
            }

            dstImage.UnlockBits(dstData);
            b.UnlockBits(srcData);

            return dstImage;
        } // end of Transpose


        /// <summary>
        /// 对图像进行水平方向倾斜变换，基准点为平行四边形左上顶点
        /// </summary>
        /// <param name="b">原始图像</param>
        /// <param name="horz">水平方向倾斜量</param>
        /// <returns></returns>
        public static Bitmap SlantHorz(Bitmap b, int horz)
        {
            int width = b.Width;
            int height = b.Height;

            Bitmap dstImage = new Bitmap(Math.Abs(horz) + width, height);

            System.Drawing.Graphics g = System.Drawing.Graphics.FromImage(dstImage);

            g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.Bilinear;
            g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.HighQuality;

            Point[] dstPoints = new Point[3];

            // 分别计算平行四边形的左上、右上和左下三个顶点
            if (horz != 0)
            {
                if (horz > 0)
                {
                    dstPoints[0] = new Point(horz, 0);
                    dstPoints[1] = new Point(width + horz, 0);
                    dstPoints[2] = new Point(0, height);
                }
                else
                {
                    dstPoints[0] = new Point(0, 0);
                    dstPoints[1] = new Point(width, 0);
                    dstPoints[2] = new Point(-horz, height);
                }
            }

            // 绘水平倾斜图
            g.DrawImage(b, dstPoints);

            g.Save();
            g.Dispose();

            return dstImage;
        } // end of SlantHorz


        /// <summary>
        /// 对图像进行垂直方向倾斜变换，基准点为平行四边形左上顶点
        /// </summary>
        /// <param name="b">原始图像</param>
        /// <param name="vert">垂直方向倾斜量</param>
        /// <returns></returns>
        public static Bitmap SlantVert(Bitmap b, int vert)
        {
            int width = b.Width;
            int height = b.Height;

            Bitmap dstImage = new Bitmap(width, Math.Abs(vert) + height);

            System.Drawing.Graphics g = System.Drawing.Graphics.FromImage(dstImage);

            g.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.Bilinear;
            g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.HighQuality;

            Point[] dstPoints = new Point[3];

            // 分别计算平行四边形的左上、右上和左下三个顶点
            if (vert != 0)
            {
                if (vert > 0)
                {
                    dstPoints[0] = new Point(0, vert);
                    dstPoints[1] = new Point(width, 0);
                    dstPoints[2] = new Point(0, height + vert);
                }
                else
                {
                    dstPoints[0] = new Point(0, 0);
                    dstPoints[1] = new Point(width, -vert);
                    dstPoints[2] = new Point(0, height);
                }

            }

            // 绘水平倾斜图
            g.DrawImage(b, dstPoints);

            g.Save();
            g.Dispose();

            return dstImage;
        } // end of SlantVert


        /// <summary>
        /// 修整
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="trimAway">修整范围</param>
        /// <returns></returns>
        public static Bitmap Trim(Bitmap b, TrimMode trimAway)
        {
            int width = b.Width;
            int height = b.Height;
            int BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            // 原始图像大小
            int area = width * height;

            BitmapData data = b.LockBits(new Rectangle(0, 0, width, height),
              ImageLockMode.ReadWrite, b.PixelFormat);

            int stride = data.Stride;
            System.IntPtr scan0 = data.Scan0;
            int offset = stride - width * BPP;

            // 修整后的有效区域
            int rectTop = 0;
            int rectBottom = height;
            int rectLeft = 0;
            int rectRight = width;

            unsafe
            {
                byte* p = (byte*)scan0;
                int i = 1;

                // 上
                if ((trimAway & TrimMode.Top) == TrimMode.Top)
                {
                    p = (byte*)scan0;
                    i = 1;
                    while (i < area)
                    {
                        if (p[3] != 0)
                        {
                            // 获取上边界
                            rectTop = i / width;
                            break;
                        }

                        p += BPP;
                        if (i % width == 0)
                            p += offset;
                        i++;
                    } // while
                }

                // 下
                if ((trimAway & TrimMode.Bottom) == TrimMode.Bottom)
                {
                    p = (byte*)scan0 + (height - 1) * stride + (width - 1) * BPP;
                    i = 1;
                    while (i < area)
                    {
                        if (p[3] != 0)
                        {
                            // 获取下边界
                            rectBottom = height - i / width;
                            break;
                        }

                        p -= BPP;
                        if (i % width == 0)
                            p -= offset;
                        i++;
                    } // while
                }

                // 左
                if ((trimAway & TrimMode.Left) == TrimMode.Left)
                {
                    p = (byte*)scan0;
                    i = 1;
                    while (i < area)
                    {
                        if (p[3] != 0)
                        {
                            // 获取左边界
                            rectLeft = i / height;
                            break;
                        }

                        p += stride;
                        if (i % height == 0)
                        {
                            p = (byte*)scan0;
                            p += (i / height) * BPP;
                        }
                        i++;
                    } // while
                }

                // 右
                if ((trimAway & TrimMode.Right) == TrimMode.Right)
                {
                    p = (byte*)scan0 + (width - 1) * BPP;
                    i = 1;
                    while (i < area)
                    {
                        if (p[3] != 0)
                        {
                            // 获取右边界
                            rectRight = width - i / height;
                            break;
                        }

                        p += stride;
                        if (i % height == 0)
                        {
                            p = (byte*)scan0 + (width - 1) * BPP;
                            p -= (i / height) * BPP;
                        }
                        i++;
                    } // while
                }
            }

            b.UnlockBits(data);


            // 修整有效区域后的图像
            Bitmap dst = new Bitmap(rectRight - rectLeft, rectBottom - rectTop);
            System.Drawing.Graphics g = System.Drawing.Graphics.FromImage(dst);
            g.DrawImage(b,
              new Rectangle(0, 0, dst.Width, dst.Height),
              new Rectangle(rectLeft, rectTop, dst.Width, dst.Height),
              GraphicsUnit.Pixel
              );
            g.Save();

            return dst;
        } // end of Trim
    }// end of ImageTransform

    /// <summary>
    /// 逻辑算法类
    /// </summary>
    public class MathLogic
    {
        /// <summary>
        /// 利用快速排序查找中值
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        public static byte SearchMid(List<byte> sequence)
        {
            int len = sequence.Count;
            /*
            bool isMovable = true;
            byte t = 0;
            //先排序（升序）
            while (isMovable)
            {
                isMovable = false;
                for (int i = 1; i < len; i++)
                {
                    if (sequence[i - 1] > sequence[i])
                    {
                        // 交换 sequence[i-1], sequence[i]
                        Swap(sequence[i - 1], sequence[i]);

                        isMovable = true;
                    }
                } // i
            } // isMovable*/
            List<byte> newSequence = sequence;
            if (sequence.Count > 1)
                newSequence = BinSort(sequence, 0, sequence.Count - 1);
            // 取中值
            if (newSequence.Count > len / 2)
                return newSequence[len / 2];
            else
                return 0;
        }
        /*
        private static void Swap(byte a, byte b)
        {
            byte t = a;
            a = b;
            b = t;
        }*/

        /// <summary>
        /// 二分快速排序
        /// </summary>
        /// <param name="a"></param>
        /// <param name="l"></param>
        /// <param name="r"></param>
        /// <returns></returns>
        public static List<byte> BinSort(List<byte> a, int l, int r)
        {
            int mi = (l + r) / 2;
            int mid = a[mi]; //将当前序列在中间位置的数定义为中间数
            int i = l, j = r;
            do
            {
                while (a[i] < mid && i < r) i++; //在左半部分寻找比中间数大的数
                while (a[j] > mid && j > l) j--; //在右半部分寻找比中间数小的数
                if (i <= j)  //若找到一组与排序目标不一致的数对则交换它们
                {
                    byte tmp = a[i];
                    a[i] = a[j];
                    a[j] = tmp;
                    i++;
                    j--;//继续找
                }
            } while (i <= j);
            if (l < j) BinSort(a, l, j); //若未到两个数的边界，则递归搜索左右区间
            if (r > i) BinSort(a, i, r);
            return a;
        }

        /// <summary>
        /// 二分(对半)查找
        /// </summary>
        /// <param name="a">要查找的数组</param>
        /// <param name="k">要查找的对象</param>
        /// <returns></returns>
        public static int BinSearch(List<int> a, int k)
        {
            int low, hig, mid;
            low = 1;
            hig = a.Count;
            mid = (low + hig) / 2;
            while (a[mid] != k && low <= hig)
            {
                if (a[mid] > k)
                    hig = mid - 1;
                else
                    low = mid + 1;
                mid = (low + hig) / 2;
            }
            if (low > hig)
                mid = 0;
            return mid;
        }
    }
    #endregion
}
