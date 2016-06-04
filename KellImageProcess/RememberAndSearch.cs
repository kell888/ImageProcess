using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;

namespace KellImageProcess
{
    /// <summary>
    /// 可序列化的记忆对象结构
    /// </summary>
    [Serializable]
    public struct RememberObject
    {
        /// <summary>
        /// 实际位置
        /// </summary>
        public Point RealPosition;
        /// <summary>
        /// 未经过二值化的记忆区域中的位图
        /// </summary>
        public Bitmap BaseBmp;
        /// <summary>
        /// 当前记忆对象的灰度阀值
        /// </summary>
        public byte Threshold;
        /// <summary>
        /// 表示一个空RememberObject对象
        /// </summary>
        public static readonly RememberObject Empty;
        /// <summary>
        /// 操作符重载，相等
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <returns></returns>
        public static bool operator ==(RememberObject left, RememberObject right)
        {
            return ((left.RealPosition == right.RealPosition) && (left.BaseBmp == right.BaseBmp) && (left.Threshold == right.Threshold));
        }
        /// <summary>
        /// 操作符重载，不相等
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <returns></returns>
        public static bool operator !=(RememberObject left, RememberObject right)
        {
            return !(left == right);
        }
        /// <summary>
        /// 重载Equals
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj)
        {
            if (!(obj is RememberObject))
            {
                return false;
            }
            RememberObject rem = (RememberObject)obj;
            return ((rem.RealPosition == this.RealPosition) && (rem.BaseBmp == this.BaseBmp) && (rem.Threshold == this.Threshold));
        }
        /// <summary>
        /// 重载GetHashCode
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode()
        {
            return base.GetHashCode();
        }
    }

    /// <summary>
    /// 记忆和匹配类
    /// </summary>
    public class RememberAndSearch
    {
        static int BPP = 4;
        static int blkORwht = 1;
        Rectangle rememberRect;
        Bitmap rememberImage;
        Bitmap offsetImage;
        Bitmap baseImage;
        Bitmap tmpImage;
        bool stop;
        int matchCount;
        Size pointSize = new Size(10, 10);
        static int SplitWidth;
        static int SplitHeight;
        int matchRate = 10;
        static byte threshold = 200;
        bool needGrayExtend;
        int reservedCount = 3;
        List<Suited> list;
        Rectangle effectRect;
        Rectangle checkRect;
        Rectangle offsetCheckRect;
        RememberObject rememberObject = RememberObject.Empty;
        bool showInTime;
        Point rememberPosition;
        Point objectBasePosition;
        Point objectOffsetPosition;
        List<int> baseString;
        Size effectSize;
        Point _offsetPosition;
        bool needFilter = true;
        static PixelFormat pf = PixelFormat.Format32bppArgb;
        //System.Timers.Timer timer = new System.Timers.Timer(10);
        //bool finish1, finish2;
        int maxMatch = 0;//匹配率
        List<Suited> resultList = new List<Suited>();

        /// <summary>
        /// 构造函数
        /// </summary>
        /// <param name="size">要记忆的基准位图的大小</param>
        public RememberAndSearch(Size size)
        {
            rememberImage = new Bitmap(size.Width, size.Height);
            baseImage = new Bitmap(size.Width, size.Height);
            tmpImage = new Bitmap(size.Width, size.Height);
            //timer.Elapsed += new System.Timers.ElapsedEventHandler(timer_Elapsed);
        }
        /*
        private void timer_Elapsed(object sender, System.Timers.ElapsedEventArgs e)
        {
            if (finish1 && finish2)
                timer.Stop();
        }
        */
        /// <summary>
        /// 是否需要除噪
        /// </summary>
        public bool NeedFilter
        {
            get
            {
                return needFilter;
            }
            set
            {
                needFilter = value;
            }
        }

        /// <summary>
        /// 获取基准位图中目标的实际位置
        /// </summary>
        public Point ObjectBasePosition
        {
            get
            {
                return objectBasePosition;
            }
        }

        /// <summary>
        /// 获取偏移位图中目标的实际位置
        /// </summary>
        public Point ObjectOffsetPosition
        {
            get
            {
                return objectOffsetPosition;
            }
        }

        /// <summary>
        /// 获取记忆的基准位图的实际位置
        /// </summary>
        public Point RememberPosition
        {
            get
            {
                return rememberPosition;
            }
        }

        /// <summary>
        /// 获取偏移对象需要匹配的有效区域
        /// </summary>
        public Rectangle OffsetEffectRect
        {
            get
            {
                return offsetCheckRect;
            }
        }

        /// <summary>
        /// 获取或设置搜索匹配时是否要实时显示
        /// </summary>
        public bool ShowInTime
        {
            get
            {
                return showInTime;
            }
            set
            {
                showInTime = value;
            }
        }

        /// <summary>
        /// 获取或设置二值化阀值
        /// </summary>
        public byte Threshold
        {
            get
            {
                return threshold;
            }
            set
            {
                threshold = value;
            }
        }

        /// <summary>
        /// 获取或设置二值化时是否需要灰度拉伸
        /// </summary>
        public bool NeedGrayExtend
        {
            get
            {
                return needGrayExtend;
            }
            set
            {
                needGrayExtend = value;
            }
        }

        /// <summary>
        /// 获取当前的记忆对象
        /// </summary>
        public RememberObject TheRememberObject
        {
            get
            {
                return rememberObject;
            }
        }

        /// <summary>
        /// 获取当前状态是否曾经记忆基准对象
        /// </summary>
        public bool HasRemember
        {
            get
            {
                return rememberObject.BaseBmp != null;
            }
        }

        /// <summary>
        /// 清空记忆对象
        /// </summary>
        public void ClearRemember()
        {
            rememberObject = RememberObject.Empty;
        }

        /// <summary>
        /// 获取基准对象的临时位图
        /// </summary>
        public Bitmap BaseImage
        {
            get
            {
                try
                {
                    return (Bitmap)baseImage.Clone();
                }
                catch { return null; }
            }
        }

        /// <summary>
        /// 获取偏移对象的临时匹配位图
        /// </summary>
        public Bitmap TempImage
        {
            get
            {
                try
                {
                    return (Bitmap)tmpImage.Clone();
                }
                catch { return null; }
            }
        }

        /// <summary>
        /// 获取当前搜索匹配的位置
        /// </summary>
        public Point CurrentSearchPosition
        {
            get
            {
                return checkRect.Location;
            }
        }

        /// <summary>
        /// 获取基准对象需要匹配的有效区域
        /// </summary>
        public Rectangle BaseEffectRectangle
        {
            get
            {
                return effectRect;
            }
        }

        /// <summary>
        /// 获取或设置要保留的最匹配的位置对象个数
        /// </summary>
        public int ReservedCount
        {
            get
            {
                return reservedCount;
            }
            set
            {
                reservedCount = value;
            }
        }

        /// <summary>
        /// 获取当前搜索中最匹配的位置对象列表(默认最多reservedCount个)
        /// </summary>
        public List<Suited> SuitedPositionList
        {
            get
            {
                return list;
            }
        }

        /// <summary>
        /// 获取或设置可接受的匹配率(范围：(0,100]，默认为10)
        /// </summary>
        public int MatchRate
        {
            get
            {
                return matchRate;
            }
            set
            {
                matchRate = value;
            }
        }

        /// <summary>
        /// 获取或设置匹配的精度(值越小越精确)，范围：(2, 2)到(rememberImage.Width / 2, rememberImage.Height / 2)，默认为(10, 10)
        /// </summary>
        public Size PointSize
        {
            get
            {
                return pointSize;
            }
            set
            {
                Size tmp = value;
                if (value.Width < 2 || value.Height < 2)
                {
                    tmp = new Size(2, 2);
                }
                else if (value.Width > rememberImage.Width / 2 || value.Height > rememberImage.Height / 2)
                {
                    int min = Math.Min(rememberImage.Width / 2, rememberImage.Height / 2);
                    tmp = new Size(min, min);
                }
                pointSize = tmp;
                SplitWidth = rememberImage.Width / pointSize.Width;
                SplitHeight = rememberImage.Height / pointSize.Height;
            }
        }

        /// <summary>
        /// 获取当前匹配的次数
        /// </summary>
        public int MatchCount
        {
            get
            {
                return matchCount;
            }
        }

        /// <summary>
        /// 停止匹配搜索
        /// </summary>
        public void StopSearch()
        {
            stop = true;
        }

        /// <summary>
        /// 获取或设置白匹配还是黑匹配
        /// 0为黑匹配，其他值为白匹配(默认)
        /// </summary>
        public int BlkORwht
        {
            get
            {
                return blkORwht;
            }
            set
            {
                blkORwht = value;
            }
        }

        /// <summary>
        /// 获取记忆中的基准对象的二值位图
        /// </summary>
        public Bitmap RememberImage
        {
            get
            {
                return rememberImage;
            }
        }

        /// <summary>
        /// 获取要匹配偏移对象的二值位图
        /// </summary>
        public Bitmap OffsetImage
        {
            get
            {
                return offsetImage;
            }
        }

        /// <summary>
        /// 将记忆对象保存到文件
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns></returns>
        public bool SaveRememberObjectToFile(string filePath)
        {
            if (rememberObject == RememberObject.Empty)
                return false;
            try
            {
                IFormatter formatter = new BinaryFormatter();
                using (Stream stream = new FileStream(filePath, FileMode.Create, FileAccess.Write, FileShare.None))
                {
                    formatter.Serialize(stream, rememberObject);
                    return true;
                }
            }
            catch
            {
                return false;
            }
        }

        /// <summary>
        /// 从文件中载入记忆对象
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns></returns>
        public bool LoadRememberObjectFromFile(string filePath)
        {
            try
            {
                IFormatter formatter = new BinaryFormatter();
                using (Stream stream = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    rememberObject = (RememberObject)formatter.Deserialize(stream);
                    rememberImage = Bitize((Bitmap)rememberObject.BaseBmp.Clone(), threshold, needGrayExtend, needFilter);
                    if (rememberImage.Width < 20 && pointSize.Width >= 10)
                        pointSize.Width = rememberImage.Width / 2;
                    if (rememberImage.Height < 20 && pointSize.Height >= 10)
                        pointSize.Height = rememberImage.Height / 2;
                    SplitWidth = rememberImage.Width / pointSize.Width;
                    SplitHeight = rememberImage.Height / pointSize.Height;
                    rememberRect = new Rectangle(rememberObject.RealPosition, new Size(rememberImage.Width, rememberImage.Height));
                    rememberPosition = rememberObject.RealPosition;
                    threshold = rememberObject.Threshold;
                    this.tmpImage = new Bitmap(rememberImage.Width, rememberImage.Height);
                    this.baseImage = new Bitmap(rememberImage.Width, rememberImage.Height);
                    return true;
                }
            }
            catch
            {
                return false;
            }
        }

        /// <summary>
        /// 记忆基准位图(搜索匹配Search()之前要先调用此方法)，搜索过一次以后若想要记忆另外的对象则要调用此方法(载入新的记忆对象)
        /// </summary>
        /// <param name="baseImage">记忆区域中的基准位图</param>
        /// <param name="realPosition">基准位图的实际位置</param>
        /// <param name="thresholdValue">二值化阀值</param>
        /// <param name="needGrayExtend">是否需要灰度拉伸</param>
        public void RememberBaseImage(Bitmap baseImage, Point realPosition, byte thresholdValue, bool needGrayExtend)
        {
            if (baseImage == null)
                return;
            rememberObject.BaseBmp = baseImage;
            rememberObject.RealPosition = realPosition;
            rememberObject.Threshold = thresholdValue;

            rememberImage = Bitize((Bitmap)baseImage.Clone(), threshold, needGrayExtend, needFilter);

            if (rememberImage.Width < 20 && pointSize.Width >= 10)
                pointSize.Width = rememberImage.Width / 2;
            if (rememberImage.Height < 20 && pointSize.Height >= 10)
                pointSize.Height = rememberImage.Height / 2;
            SplitWidth = rememberImage.Width / pointSize.Width;
            SplitHeight = rememberImage.Height / pointSize.Height;
            rememberRect = new Rectangle(realPosition, new Size(rememberImage.Width, rememberImage.Height));
            rememberPosition = realPosition;
            threshold = thresholdValue;
            this.needGrayExtend = needGrayExtend;
            this.tmpImage = new Bitmap(rememberImage.Width, rememberImage.Height);
            this.baseImage = new Bitmap(rememberImage.Width, rememberImage.Height);

            baseString = GraspRawData(rememberImage, out effectRect);

            objectBasePosition = new Point(realPosition.X + effectRect.X, realPosition.Y + effectRect.Y);

            effectSize = new Size(effectRect.Width, effectRect.Height);
        }

        /// <summary>
        /// 以指定的颜色刷新画布
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="clearColor"></param>
        public void ClearGraphics(Bitmap bmp, Color clearColor)
        {
            if (bmp != null)
            {
                Graphics g = Graphics.FromImage(bmp);
                g.FillRectangle(new SolidBrush(clearColor), 0, 0, bmp.Width, bmp.Height);
                g.Dispose();
            }
        }

        /// <summary>
        /// 记忆新的基准对象并开始搜索匹配，并返回最匹配的位置列表(最多ReservedCount个)，且一旦记忆过基准位图就不要用此方法，而应该用Search(Bitmap, Point)，否则会浪费时间
        /// </summary>
        /// <param name="baseImage">基准位图</param>
        /// <param name="realPosition">新的基准位图的实际位置</param>
        /// <param name="threshold">二值化阀值</param>
        /// <param name="needGrayExtend">是否需要灰度拉伸</param>
        /// <param name="offsetBmp">偏移对象位图</param>
        /// <param name="offsetPosition">偏移位图的位置</param>
        /// <returns></returns>
        public List<Suited> Search(Bitmap baseImage, Point realPosition, byte threshold, bool needGrayExtend, Bitmap offsetBmp, Point offsetPosition)
        {
            if (baseImage != null)
            {
                RememberBaseImage(baseImage, realPosition, threshold, needGrayExtend);

                if (effectRect == Rectangle.Empty || effectRect.IsEmpty || effectRect.Width == 0 || effectRect.Height == 0)
                {
                    return null;
                }

                offsetImage = Bitize((Bitmap)offsetBmp.Clone(), threshold, needGrayExtend, needFilter);

                Rectangle rec = GetEffectRect(offsetImage);
                offsetCheckRect = new Rectangle(offsetPosition.X + rec.X, offsetPosition.Y + rec.Y, rec.Width, rec.Height);

                if (offsetCheckRect == Rectangle.Empty || offsetCheckRect.IsEmpty || offsetCheckRect.Width == 0 || offsetCheckRect.Height == 0)
                {
                    return null;
                }

                if (effectRect.Width > offsetCheckRect.Width || effectRect.Height > offsetCheckRect.Height)
                {
                    return null;
                }

                _offsetPosition = offsetPosition;

                list = RecognizeObject(rec);

                return list;
            }
            else
            {
                return null;
            }
        }

        /// <summary>
        /// 在已有的记忆中开始搜索匹配，并返回最匹配的位置列表(最多ReservedCount个)
        /// </summary>
        /// <param name="offsetBmp">偏移对象位图</param>
        /// <param name="offsetPosition">偏移位图的位置</param>
        /// <returns></returns>
        public List<Suited> Search(Bitmap offsetBmp, Point offsetPosition)
        {
            if (rememberObject == RememberObject.Empty)
                return null;
            if (rememberImage != null)
            {
                if (effectRect == Rectangle.Empty || effectRect.IsEmpty || effectRect.Width == 0 || effectRect.Height == 0)
                {
                    return null;
                }
                
                offsetImage = Bitize((Bitmap)offsetBmp.Clone(), threshold, needGrayExtend, needFilter);

                Rectangle rec = GetEffectRect(offsetImage);
                offsetCheckRect = new Rectangle(offsetPosition.X + rec.X, offsetPosition.Y + rec.Y, rec.Width, rec.Height);

                if (offsetCheckRect == Rectangle.Empty || offsetCheckRect.IsEmpty || offsetCheckRect.Width == 0 || offsetCheckRect.Height == 0)
                {
                    return null;
                }
                
                if (effectRect.Width > offsetCheckRect.Width || effectRect.Height > offsetCheckRect.Height)
                {
                    return null;
                }
                
                _offsetPosition = offsetPosition;

                list = RecognizeObject(rec);

                return list;
            }
            else
            {
                return null;
            }
        }

        /// <summary>
        /// 获取指定位图的有效匹配区域
        /// </summary>
        /// <param name="bmp">二值位图</param>
        /// <returns></returns>
        public static Rectangle GetEffectRect(Bitmap bmp)
        {
            if (bmp != null)
            {
                BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                bool bool1stScan = true;
                int ax, ay, bx, by;
                ax = ay = bx = by = 0;
                //分割手写体字符的图片
                BitmapData dat = bmp.LockBits(new Rectangle(0, 0, bmp.Width, bmp.Height), ImageLockMode.ReadOnly, bmp.PixelFormat);
                unsafe
                {
                    byte* t = (byte*)dat.Scan0;
                    int strid = dat.Stride;
                    int ro = (bmp.Width / SplitWidth) / 2;
                    int co = (bmp.Height / SplitHeight) / 2;

                    for (int j = 0; j < bmp.Height; j += co)
                    {
                        t = (byte*)dat.Scan0 + strid * j;
                        for (int i = 0; i < bmp.Width; i += ro)
                        {
                            byte keyValue = 255;
                            if (blkORwht == 0)//黑边
                                keyValue = 0;
                            if (t[0] == keyValue)
                            {
                                if (!bool1stScan)
                                {
                                    if (i < ax) ax = i;
                                    if (i >= bx) bx = i;
                                    if (j < ay) ay = j;
                                    if (j >= by) by = j;
                                }
                                else
                                {
                                    bool1stScan = false;
                                    ax = i;
                                    bx = i;
                                    ay = j;
                                    by = j;
                                }
                            }
                            t += BPP * ro;
                        }
                    }
                }
                bmp.UnlockBits(dat);

                return new Rectangle(ax, ay, bx - ax, by - ay);
            }
            else
            {
                return Rectangle.Empty;
            }
        }

        /// <summary>
        /// 获取最相似的识别对象数组
        /// </summary>
        /// <param name="offsetCheckRect">偏移位图的有效匹配区域</param>
        /// <returns></returns>
        private List<Suited> RecognizeObject(Rectangle offsetCheckRect)
        {
            if (baseString == null || baseString.Count == 0)
                return null;

            resultList.Clear();
            matchCount = 0;
            /*
            finish1 = false;
            finish2 = false;

            timer.Start();
            System.Threading.Thread thr1 = new System.Threading.Thread(new System.Threading.ThreadStart(GetMatchString1));
            System.Threading.Thread thr2 = new System.Threading.Thread(new System.Threading.ThreadStart(GetMatchString2));
            thr1.Start();
            thr2.Start();*/
            /*while (!finish1 || !finish2)
            {
                System.Threading.Thread.Sleep(10);
            }*/
            for (int j = offsetCheckRect.Y; j < offsetCheckRect.Y + offsetCheckRect.Height; j += pointSize.Height)
            {
                for (int i = offsetCheckRect.X; i < offsetCheckRect.X + offsetCheckRect.Width; i += pointSize.Width)
                {
                    if (stop)
                        break;
                    if (i <= offsetCheckRect.Width + offsetCheckRect.X - effectSize.Width && j <= offsetCheckRect.Height + offsetCheckRect.Y - effectSize.Height)
                    {
                        Point p = new Point(i, j);
                        checkRect = new Rectangle(p, effectSize);
                        List<int> arrCurrentEigenvalue = CheckGraspRawData(FastClipBitmap(offsetImage, checkRect));
                        maxMatch = Match(baseString, arrCurrentEigenvalue);
                        //保留的匹配率,建议10%
                        if (maxMatch >= matchRate)
                        {
                            Suited suited = new Suited(new Point(p.X + _offsetPosition.X, p.Y + _offsetPosition.Y), maxMatch);//位置
                            resultList.Add(suited);

                            objectOffsetPosition = new Point(_offsetPosition.X + checkRect.X, _offsetPosition.Y + checkRect.Y);
                        }
                        matchCount++;
                    }
                }
                if (stop)
                {
                    stop = false;
                    break;
                }
            }

            /*if (stop)
                stop = false;*/

            //找出匹配率最高的
            resultList.Sort(0, resultList.Count, new MaxMatchSort());

            //分离结果
            //保存匹配率最高的reservedCount个位置
            while (resultList.Count > reservedCount)
            {
                resultList.RemoveAt(resultList.Count - 1);
            }
            return resultList;
        }
        /*
        private void GetMatchString1()
        {
            for (int j = offsetCheckRect.Y; j < offsetCheckRect.Y + offsetCheckRect.Height / 2; j += pointSize.Height)
            {
                for (int i = offsetCheckRect.X; i < offsetCheckRect.X + offsetCheckRect.Width / 2; i += pointSize.Width)
                {
                    if (stop)
                    {
                        finish1 = true;
                        return;
                    }
                    if (i <= offsetCheckRect.Width + offsetCheckRect.X - effectSize.Width && j <= offsetCheckRect.Height + offsetCheckRect.Y - effectSize.Height)
                    {
                        Point p = new Point(i, j);
                        checkRect = new Rectangle(p, effectSize);
                        List<string> arrCurrentEigenvalue = CheckGraspRawData(FastClipBitmap(offsetImage, checkRect));
                        maxMatch = Match(baseString, arrCurrentEigenvalue);
                        //保留的匹配率,建议10%
                        if (maxMatch >= matchRate)
                        {
                            Suited suited = new Suited(p, maxMatch);//位置
                            resultList.Add(suited);

                            objectOffsetPosition = new Point(_offsetPosition.X + checkRect.X, _offsetPosition.Y + checkRect.Y);
                        }
                        matchCount++;
                    }
                }
            }
            finish1 = true;
        }

        private void GetMatchString2()
        {
            for (int j = offsetCheckRect.Y + offsetCheckRect.Height / 2; j < offsetCheckRect.Y + offsetCheckRect.Height; j += pointSize.Height)
            {
                for (int i = offsetCheckRect.X + offsetCheckRect.Width / 2; i < offsetCheckRect.X + offsetCheckRect.Width; i += pointSize.Width)
                {
                    if (stop)
                    {
                        finish1 = true;
                        return;
                    }
                    if (i <= offsetCheckRect.Width + offsetCheckRect.X - effectSize.Width && j <= offsetCheckRect.Height + offsetCheckRect.Y - effectSize.Height)
                    {
                        Point p = new Point(i, j);
                        checkRect = new Rectangle(p, effectSize);
                        List<string> arrCurrentEigenvalue = CheckGraspRawData(FastClipBitmap(offsetImage, checkRect));
                        maxMatch = Match(baseString, arrCurrentEigenvalue);
                        //保留的匹配率,建议10%
                        if (maxMatch >= matchRate)
                        {
                            Suited suited = new Suited(p, maxMatch);//位置
                            resultList.Add(suited);

                            objectOffsetPosition = new Point(_offsetPosition.X + checkRect.X, _offsetPosition.Y + checkRect.Y);
                        }
                        matchCount++;
                    }
                }
            }
            finish1 = true;
        }
        */
        /// <summary>
        /// 根据二进制编码字符串获取当前位置的匹配率
        /// </summary>
        /// <param name="Eigenvalue">匹配的基准对象二进制编码字符串</param>
        /// <param name="arrCurrentEigenvalue">需要匹配的偏移对象二进制编码字符串</param>
        /// <returns>返回匹配的百分点[0, 100]</returns>
        private int Match(List<int> Eigenvalue, List<int> arrCurrentEigenvalue)
        {
            if (Eigenvalue == null || Eigenvalue.Count == 0 || arrCurrentEigenvalue == null || arrCurrentEigenvalue.Count == 0)
                return 0;
            int maxMatch = 0;
            int c = 0;
            int e = 0;
            int d = 0;
            int maxLen = SplitWidth * SplitHeight;
            if (Eigenvalue.Count < maxLen)
                return 0;

            for (int i = 0; i < SplitWidth; i++)
            {
                for (int j = 0; j < SplitHeight; j++)
                {
                    if (d >= arrCurrentEigenvalue.Count)
                        return Convert.ToInt32(Math.Round(100.0 * maxMatch / maxLen));
                    c = Eigenvalue[d];
                    e = arrCurrentEigenvalue[d];
                    if (c == e)
                        maxMatch++;
                    /*
                    if (c == 1)//注意：此逻辑不会因为blkORwht的变化而变化
                    {
                        if (e == 1)
                            maxMatch++;
                    }*/
                    d++;
                }
            }
            return (int)(100.0 * maxMatch / maxLen);//maxMatch;//
        }

        /// <summary>
        /// 根据给定的位图获取基准对象的二进制编码结构，最为关键的代码
        /// </summary>
        /// <param name="bmp"></param>
        /// <param name="effectRect"></param>
        /// <returns></returns>
        private List<int> GraspRawData(Bitmap bmp, out Rectangle effectRect)
        {
            effectRect = Rectangle.Empty;
            if (bmp == null)
                return null;
            List<int> arrCurrentEigenvalue = new List<int>();
            if (blkORwht == 0)
                this.ClearGraphics(baseImage, Color.White);
            else
                this.ClearGraphics(baseImage, Color.Black);

            effectRect = GetEffectRect(bmp);

            //分割
            if (effectRect.Width != 0 && effectRect.Height != 0 && ((effectRect.Width > (bmp.Width / SplitWidth) / 2) || (effectRect.Height > (bmp.Height / SplitHeight) / 2)))
            {
                Bitmap memBmp = new Bitmap(tmpImage.Width, tmpImage.Height, System.Drawing.Imaging.PixelFormat.Format16bppRgb555);
                BPP = 2;
                Graphics memGp = Graphics.FromImage(memBmp);
                memGp.DrawImage(bmp, new Rectangle(0, 0, memBmp.Width - 1, memBmp.Height - 1), effectRect, GraphicsUnit.Pixel);
                int pw = (baseImage.Width / SplitWidth);
                int ph = (baseImage.Height / SplitHeight);
                BitmapData data = memBmp.LockBits(new Rectangle(0, 0, memBmp.Width, memBmp.Height), ImageLockMode.ReadOnly, memBmp.PixelFormat);
                unsafe
                {
                    byte* p = (byte*)data.Scan0;
                    int stride = data.Stride;
                    int row = memBmp.Width / SplitWidth;
                    int col = memBmp.Height / SplitHeight;

                    for (int j = 0; j < memBmp.Height; j += col)
                    {
                        p = (byte*)data.Scan0 + stride * j;
                        for (int i = 0; i < memBmp.Width; i += row)
                        {
                            byte keyValue = 255;
                            if (blkORwht == 0)//黑边
                                keyValue = 0;
                            if (p[0] == keyValue)
                            {
                                if (showInTime)
                                {
                                    using (Graphics g = Graphics.FromImage(baseImage))
                                    {
                                        g.FillRectangle(new SolidBrush(Color.Red), i - (pw / 2), j - (ph / 2), pw, ph);
                                        g.DrawRectangle(new Pen(Color.Yellow, 1f), i - (pw / 2), j - (ph / 2), pw, ph);
                                    }
                                }
                                arrCurrentEigenvalue.Add(1);
                            }
                            else
                                arrCurrentEigenvalue.Add(0);
                            p += BPP * row;
                        }
                    }
                }
                memBmp.UnlockBits(data);
                memGp.Dispose();
                memBmp.Dispose();
            }
            return arrCurrentEigenvalue;
        }

        /// <summary>
        /// 根据给定的位图的指定区域获取偏移对象的二进制编码结构，最为关键的代码
        /// </summary>
        /// <param name="bmp"></param>
        /// <returns></returns>
        private List<int> CheckGraspRawData(Bitmap bmp)
        {
            if (bmp == null)
                return null;
            Rectangle checkRect = new Rectangle(0, 0, bmp.Width, bmp.Height);
            List<int> arrCurrentEigenvalue = new List<int>();
            if (blkORwht == 0)
                this.ClearGraphics(tmpImage, Color.White);
            else
                this.ClearGraphics(tmpImage, Color.Black);
            //分割
            if (checkRect.Width != 0 && checkRect.Height != 0 && ((checkRect.Width > (bmp.Width / SplitWidth) / 2) || (checkRect.Height > (bmp.Height / SplitHeight) / 2)))
            {
                int width = tmpImage.Width;
                int height = tmpImage.Height;
                Bitmap memBmp = new Bitmap(width, height, System.Drawing.Imaging.PixelFormat.Format16bppRgb555);
                BPP = 2;
                Graphics memGp = Graphics.FromImage(memBmp);
                //MessageBox.Show(effectRect.ToString());
                memGp.DrawImage(bmp, new Rectangle(0, 0, memBmp.Width - 1, memBmp.Height - 1), checkRect, GraphicsUnit.Pixel);
                int pw = (width / SplitWidth);
                int ph = (height / SplitHeight);
                BitmapData data = memBmp.LockBits(new Rectangle(0, 0, memBmp.Width, memBmp.Height), ImageLockMode.ReadOnly, memBmp.PixelFormat);
                unsafe
                {
                    byte* p = (byte*)data.Scan0;
                    int stride = data.Stride;
                    int offset = stride - BPP * memBmp.Width;
                    int row = memBmp.Width / SplitWidth;
                    int col = memBmp.Height / SplitHeight;
                    for (int j = 0; j < memBmp.Height; j += col)
                    {
                        p = (byte*)data.Scan0 + stride * j;
                        for (int i = 0; i < memBmp.Width; i += row)
                        {
                            byte keyValue = 255;
                            if (blkORwht == 0)//黑边
                                keyValue = 0;
                            if (p[0] == keyValue)
                            {
                                if (showInTime)
                                {
                                    using (Graphics g = Graphics.FromImage(tmpImage))
                                    {
                                        g.FillRectangle(new SolidBrush(Color.Red), i - (pw / 2), j - (ph / 2), pw, ph);
                                        g.DrawRectangle(new Pen(Color.Yellow, 1f), i - (pw / 2), j - (ph / 2), pw, ph);
                                    }
                                }
                                arrCurrentEigenvalue.Add(1);
                            }
                            else
                                arrCurrentEigenvalue.Add(0);

                            p += BPP * row;
                        }
                    }
                }
                memBmp.UnlockBits(data);
                memGp.Dispose();
                memBmp.Dispose();
            }
            return arrCurrentEigenvalue;
        }

        /// <summary>
        /// 快速剪裁指定区域的位图
        /// </summary>
        /// <param name="srcBmp"></param>
        /// <param name="rect"></param>
        /// <returns></returns>
        public static Bitmap FastClipBitmap(Bitmap srcBmp, Rectangle rect)
        {
            BPP = Image.GetPixelFormatSize(srcBmp.PixelFormat) / 8;
            GraphicsUnit gu = GraphicsUnit.Pixel;
            RectangleF rectF = srcBmp.GetBounds(ref gu);
            Rectangle rec = Rectangle.Intersect(Rectangle.Truncate(rectF), rect);
            if (rec.Width == 0 || rec.Height == 0)
                return null;
            int X = rec.X;
            int Y = rec.Y;
            int Width = rec.Width + rec.X;
            int Height = rec.Height + rec.Y;
            Bitmap dstImage = new Bitmap(rec.Width, rec.Height, srcBmp.PixelFormat);
            BitmapData srcData = srcBmp.LockBits(new Rectangle(X, Y, rec.Width, rec.Height), ImageLockMode.ReadOnly, srcBmp.PixelFormat);
            BitmapData dstData = dstImage.LockBits(new Rectangle(0, 0, rec.Width, rec.Height), ImageLockMode.WriteOnly, dstImage.PixelFormat);
            int stride1 = srcData.Stride;
            int offset1 = stride1 - BPP * rec.Width;
            int stride2 = dstData.Stride;
            int offset2 = stride2 - BPP * rec.Width;
            unsafe
            {
                byte* src = (byte*)srcData.Scan0;
                byte* dst = (byte*)dstData.Scan0;
                for (int y = Y; y < Height; y++)
                {
                    for (int x = X; x < Width; x++)
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
            srcBmp.UnlockBits(srcData);
            dstImage.UnlockBits(dstData);
            return dstImage;
        }

        /// <summary>
        /// 初始化一个整型二维数组
        /// </summary>
        /// <param name="width">数组宽</param>
        /// <param name="height">数组高</param>
        /// <param name="init">初始值</param>
        /// <returns></returns>
        private static byte[,] InitArray(int width, int height, byte init)
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
        /// 将二维数组转换为灰度位图流
        /// </summary>
        /// <param name="GrayArray">灰度数组</param>
        /// <returns></returns>
        private static Bitmap Array2Image(byte[,] GrayArray)
        {
            int width = GrayArray.GetLength(0);
            int height = GrayArray.GetLength(1);
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
                        p[0] = p[1] = p[2] = GrayArray[x, y];
                        if (BPP > 3)
                            p[3] = 255;

                        p += BPP;
                    } //  x
                    p += offset;
                } // y
            }
            b.UnlockBits(data);
            return b;
        } // end of Array2Image

        /// <summary>
        /// 十字型腐蚀
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="graied">已经灰度化</param>
        /// <returns></returns>
        private static Bitmap ErosionCross(Bitmap b, int blkORwht, byte threshold, bool graied)
        {
            // 先将原始二值图转化为二维数组
            byte[,] srcGray = Image2Array(b, graied);

            // 进行十字型腐蚀处理
            byte[,] dstGray = ErosionCross(srcGray, blkORwht, threshold);

            // 转换为灰度图像
            return Array2Image(dstGray);
        } // end of ErosionCross

        /// <summary>
        /// 十字型膨胀
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="blkORwht"></param>
        /// <param name="threshold"></param>
        /// <param name="graied">已经灰度化</param>
        /// <returns></returns>
        private static Bitmap DilationCross(Bitmap b, int blkORwht, byte threshold, bool graied)
        {
            // 先将原始二值图转化为二维数组
            byte[,] srcGray = Image2Array(b, graied);

            // 进行十字型膨胀处理
            byte[,] dstGray = DilationCross(srcGray, blkORwht, threshold);

            // 转换为灰度图像
            return Array2Image(dstGray);
        } // end of DilationCross

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
            b = ErosionCross(b, blkORwht, threshold, graied);
            b = DilationCross(b, blkORwht, threshold, true);

            return b;
        } // end of Opening

        /// <summary>
        /// 利用迭代开运算的方法进行图像滤波，返回二值位图，可是有一个致命的缺点：有边框！
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
            return bmp;
        }

        /// <summary>
        /// 二值化
        /// </summary>
        /// <param name="bmp">位图流</param>
        /// <param name="threshold">阈值</param>
        /// <param name="needGrayExtend">是否需要灰度拉伸</param>
        /// <param name="needFilter">是否需要除噪</param>
        /// <returns></returns>
        public static Bitmap Bitize(Bitmap bmp, byte threshold, bool needGrayExtend, bool needFilter)
        {
            Bitmap dstImage = null;
            if (bmp != null)
            {
                Bitmap bmp1 = null;
                if (needGrayExtend)
                    bmp1 = GrayExtend(bmp);//灰度拉伸
                else
                    bmp1 = Gray(bmp);//灰度化
                if (needFilter)
                    bmp1 = Filtering(bmp1, blkORwht, threshold, 1);
                byte[,] GrayArray = BinaryArray(bmp1, threshold);
                dstImage = BinaryImage(GrayArray, Color.White, Color.Black);
                bmp1.Dispose();
            }
            return dstImage;
        }

        /// <summary>
        /// 获取经过灰度拉伸后的位图(灰度图)
        /// </summary>
        /// <param name="srcBmp"></param>
        /// <returns></returns>
        public static Bitmap GrayExtend(Bitmap srcBmp)
        {
            Bitmap dstBmp = null;
            if (srcBmp != null)
            {
                //先灰度化
                BPP = Image.GetPixelFormatSize(srcBmp.PixelFormat) / 8;
                dstBmp = Gray((Bitmap)srcBmp.Clone());
                int width = dstBmp.Width;
                int height = dstBmp.Height;
                BitmapData data = dstBmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, srcBmp.PixelFormat);
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
        /// 将位图流转换为二维数组
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="graied">已经灰度化</param>
        /// <returns></returns>
        private static byte[,] Image2Array(Bitmap b, bool graied)
        {
            pf = b.PixelFormat;
            BPP = Image.GetPixelFormatSize(b.PixelFormat) / 8;
            int width = b.Width;
            int height = b.Height;
            // 申请一个二维数组
            byte[,] GrayArray = new byte[width, height];
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
                        if (graied)
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
        /// 灰度化
        /// </summary>
        /// <param name="bmp"></param>
        /// <returns></returns>
        public static Bitmap Gray(Bitmap bmp)
        {
            Bitmap dstBmp = null;
            if (bmp != null)
            {
                BPP = Image.GetPixelFormatSize(bmp.PixelFormat) / 8;
                dstBmp = (Bitmap)bmp.Clone();
                int width = dstBmp.Width;
                int height = dstBmp.Height;
                BitmapData data = dstBmp.LockBits(new Rectangle(0, 0, width, height), ImageLockMode.ReadWrite, dstBmp.PixelFormat);
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
        /// 将灰度化的位图流转换为二值数组
        /// </summary>
        /// <param name="b">位图流</param>
        /// <param name="threshold">阈值</param>
        /// <returns></returns>
        private static byte[,] BinaryArray(Bitmap b, byte threshold)
        {
            int width = b.Width;
            int height = b.Height;
            // 将灰度图转化为灰度数组
            byte[,] GrayArray = Image2Array(b, true);
            b.Dispose();
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
        /// 将二值数组转换为双色位图流
        /// </summary>
        /// <param name="GrayArray">二值数组</param>
        /// <param name="bgColor">背景色</param>
        /// <param name="fgColor">前景色</param>
        /// <returns></returns>
        private static Bitmap BinaryImage(byte[,] GrayArray, Color bgColor, Color fgColor)
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
            BPP = Image.GetPixelFormatSize(pf) / 8;
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
                        if (GrayArray[x, y] < threshold)
                        {
                            if (BPP > 3)
                                p[3] = fgAlpha;
                            p[2] = fgRed;
                            p[1] = fgGreen;
                            p[0] = fgBlue;
                        }
                        else
                        {
                            if (BPP > 3)
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

    }

    /// <summary>
    /// 位置匹配结构
    /// </summary>
    public struct Suited
    {
        /// <summary>
        /// 位置
        /// </summary>
        public Point Location;
        /// <summary>
        /// 匹配率
        /// </summary>
        public int MaxMatch;
        /// <summary>
        /// 构造函数
        /// </summary>
        /// <param name="Location"></param>
        /// <param name="MaxMatch"></param>
        public Suited(Point Location, int MaxMatch)
        {
            this.Location = Location;
            this.MaxMatch = MaxMatch;
        }
    }

    /// <summary>
    /// 匹配率排序类
    /// </summary>
    public class MaxMatchSort : IComparer<Suited>
    {
        int IComparer<Suited>.Compare(Suited x, Suited y)
        {
            CaseInsensitiveComparer cic = new CaseInsensitiveComparer();
            return cic.Compare(y.MaxMatch, x.MaxMatch);
        }
    }
}