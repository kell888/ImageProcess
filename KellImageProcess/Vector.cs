using System;
using System.Collections.Generic;
using System.Text;

namespace KellImageProcess
{
    /// <summary>
    /// 数字测试类
    /// </summary>
    public class TestNumeric
    {
        /// <summary>
        /// 是否为数字组合
        /// </summary>
        /// <param name="str"></param>
        /// <returns></returns>
        public static bool IsNumeric(string str)
        {
            if (str == null || str.Length == 0)
                return false;
            foreach (char c in str)
            {
                if (!Char.IsNumber(c))
                {
                    return false;
                }
            }
            return true;
        }
    /// <summary>
    /// 是否为整数
    /// </summary>
    /// <param name="str"></param>
    /// <returns></returns>
        public static bool IsInt(string str)
        {
            System.Text.RegularExpressions.Regex reg1
                = new System.Text.RegularExpressions.Regex(@"^[-]?\d+\d*$");
            return reg1.IsMatch(str);
        }

        /// <summary>
        /// 是否为正整数，不同于IsNumeric()
        /// </summary>
        /// <param name="str"></param>
        /// <returns></returns>
        public static bool IsPosInt(string str)
        {
            System.Text.RegularExpressions.Regex reg1
                = new System.Text.RegularExpressions.Regex(@"^\d*$");
            return reg1.IsMatch(str);
        }
        /// <summary>
        /// 是否为实数
        /// </summary>
        /// <param name="str"></param>
        /// <returns></returns>
        public static bool IsReal(string str)
        {
            System.Text.RegularExpressions.Regex reg1
                = new System.Text.RegularExpressions.Regex(@"^[-]?\d+[.]?\d*$");
            return reg1.IsMatch(str);
        }
        /// <summary>
        /// 是否为正实数
        /// </summary>
        /// <param name="str"></param>
        /// <returns></returns>
        public static bool IsPosReal(string str)
        {
            System.Text.RegularExpressions.Regex reg1
                = new System.Text.RegularExpressions.Regex(@"^\d+[.]?\d*$");
            return reg1.IsMatch(str);
        }
    }
    /// <summary>
    /// 二维向量
    /// </summary>
    public class Vector2
    {
        /// <summary>
        /// 判断在给定的误差范围内两个数值是否相等
        /// </summary>
        /// <param name="comp1"></param>
        /// <param name="comp2"></param>
        /// <returns></returns>
        public static bool CompareEps(double comp1, double comp2)
        {
            return Math.Abs(comp1 - comp2) <= Common.Eps;
        }
        double dir = 0;
        Point2 start = new Point2();
        Point2 end = new Point2();
        /// <summary>
        /// 构造函数
        /// </summary>
        public Vector2()
        {

        }
        /// <summary>
        /// 始点为原点，终点为p的向量
        /// </summary>
        /// <param name="p"></param>
        public Vector2(Point2 p)
        {
            start = new Point2(0, 0);
            end = p;
        }
        /// <summary>
        /// 始点为startP，终点为endP的向量
        /// </summary>
        /// <param name="startP"></param>
        /// <param name="endP"></param>
        public Vector2(Point2 startP, Point2 endP)
        {
            start = startP;
            end = endP;
        }
        /// <summary>
        /// 加法操作符
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        /// <returns></returns>
        public static Vector2 operator +(Vector2 vector1, Vector2 vector2)
        {
            Vector2 vector = new Vector2(vector1.Start, new Point2(vector2.End.X + vector1.End.X - vector1.Start.X, vector2.End.Y + vector1.End.Y - vector1.Start.Y));
            return vector;
        }
        /// <summary>
        /// 减法操作符
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        /// <returns></returns>
        public static Vector2 operator -(Vector2 vector1, Vector2 vector2)
        {
            Vector2 vector = new Vector2(vector2.End, vector1.End);
            return vector;
        }
        /// <summary>
        /// 规范化
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static Vector2 Formatting(Vector2 vector)
        {
            Point2 offset = vector.start;
            Vector2 vv = new Vector2(new Point2(vector.end.X - offset.X, vector.end.Y - offset.Y));
            return vv;
        }
        /// <summary>
        /// 是否相等
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        /// <returns></returns>
        public static bool IsEqual(Vector2 vector1, Vector2 vector2)
        {
            Vector2 v1 = Vector2.Formatting(vector1);
            Vector2 v2 = Vector2.Formatting(vector2);
            return Point2.IsEqual(v1.start, v2.start) && Point2.IsEqual(v1.end, v2.end);
        }
        /// <summary>
        /// 方向
        /// </summary>
        public double Direction
        {
            get
            {
                return dir = VectorDirection.GetDirection2(this);
            }
        }
        /// <summary>
        /// 获取或设置向量的始点
        /// </summary>
        public Point2 Start
        {
            get
            {
                return start;
            }
            set
            {
                start = value;
            }
        }
        /// <summary>
        /// 获取或设置向量的终点
        /// </summary>
        public Point2 End
        {
            get
            {
                return end;
            }
            set
            {
                end = value;
            }
        }
        /// <summary>
        /// 空向量
        /// </summary>
        public static Vector2 Empty
        {
            get
            {
                return new Vector2(new Point2(0, 0));
            }
        }
        /// <summary>
        /// 将向量转化为可视化字符串
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return start.X.ToString() + "," + start.Y.ToString() + " ----> " + end.X.ToString() + "," + end.Y.ToString();
        }
        /// <summary>
        /// 长度(模)
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static double GetLength(Vector2 vector)
        {
            Point2 startP = vector.Start;
            Point2 endP = vector.End;
            return Math.Sqrt((endP.X - startP.X) * (endP.X - startP.X) + (endP.Y - startP.Y) * (endP.Y - startP.Y));
        }
        /// <summary>
        /// 获取直线上距离开始点len的点
        /// </summary>
        /// <param name="vector">直线的向量</param>
        /// <param name="len">距离</param>
        /// <returns></returns>
        public static Point2 GetPointOfLengthFormStartP(Vector2 vector, double len)
        {
            Point2 startP = vector.Start;
            double k = Math.Tan(vector.Direction);
            int sign = 1;
            if (!double.IsNaN(k))
                sign = Math.Sign(k);
            //else
                //sign = -1;
            double d2 = len * len;
            float x = startP.X - (float)Math.Sqrt(d2 / (1 + k * k));
            float y = startP.Y + (float)Math.Sqrt(d2 - (d2 / (1 + k * k))) * sign;
            return new Point2(x, y);
            /*if (len < 0)
            {
                float x = 0;
                float y = 0;
                if (vector.End.X - vector.Start.X != 0)
                {
                    double angle = Math.Atan2(vector.End.Y - vector.Start.Y, vector.End.X - vector.Start.X);
                    x = (float)Math.Abs(len * Math.Cos(angle));//保证为正！
                    y = (float)vector.Direction * x;
                }
                else
                {
                    double angle = Math.Atan2(vector.End.Y - vector.Start.Y, vector.End.X - vector.Start.X);
                    y = (float)Math.Abs(len * Math.Sin(angle));//保证为正！
                    x = y / (float)vector.Direction;
                }
                return new Point2(vector.Start.X - x, vector.Start.Y - y);
            }
            else if (len == 0)
            {
                return vector.Start;
            }
            else
            {
                float x = 0;
                float y = 0;
                if (vector.End.X - vector.Start.X != 0)
                {
                    double angle = Math.Atan2(vector.End.Y - vector.Start.Y, vector.End.X - vector.Start.X);
                    x = (float)Math.Abs(len * Math.Cos(angle));//保证为正！
                    y = (float)vector.Direction * x;
                }
                else
                {
                    double angle = Math.Atan2(vector.End.Y - vector.Start.Y, vector.End.X - vector.Start.X);
                    y = (float)Math.Abs(len * Math.Sin(angle));//保证为正！
                    x = y / (float)vector.Direction;
                }
                return new Point2(vector.Start.X + x, vector.Start.Y + y);
            }*/
        }
        /// <summary>
        /// 缩放
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="scalingFactor"></param>
        /// <returns></returns>
        public static Vector2 Scale(Vector2 vector, double scalingFactor)
        {
            double len=Vector2.GetLength(vector);
            len *= scalingFactor;
            Point2 newEnd = GetPointOfLengthFormStartP(vector, len);
            return new Vector2(vector.Start, newEnd);
        }
        /// <summary>
        /// 单位向量
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static Vector2 GetUnitVector(Vector2 vector)
        {
            vector = Vector2.Formatting(vector);
            float x = vector.End.X;
            float y = vector.End.Y;
            float length = (float)Math.Sqrt(x * x + y * y);
            x = x / length;
            y = y / length;
            Vector2 unit = new Vector2(new Point2(x, y));
            return unit;
        }
        /// <summary>
        /// 根据有序三点(中间点为顶点)获取角度(0,180)
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <returns></returns>
        public double GetJJ(Point2 p1, Point2 p2, Point2 p3)
        {
            double angle;
            double jj1, jj2;
            jj1 = Math.Atan2(p2.Y - p1.Y, p1.X - p2.X);
            jj2 = Math.Atan2(p2.Y - p3.Y, p3.X - p2.X);
            if (jj1 < 0)
                jj1 += 2 * Math.PI;
            if (jj2 < 0)
                jj2 += 2 * Math.PI;
            double mina = 0, maxa = 0;
            mina = jj1 < jj2 ? jj1 : jj2;
            maxa = jj1 > jj2 ? jj1 : jj2;
            angle = (maxa - mina) * 180 / Math.PI;
            if (angle > 180)
                angle = 360 - angle;
            return angle;
        }
        /// <summary>
        /// 两向量的夹角[0,180)
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static double GetJJ(Vector2 v1, Vector2 v2)
        {
            double angle;
            double jj1 = Math.Atan2(v1.Start.Y - v1.end.Y, v1.end.X - v1.Start.X);
            double jj2 = Math.Atan2(v2.Start.Y - v2.end.Y, v2.end.X - v2.Start.X);
            if (jj1 < 0)
                jj1 += 2 * Math.PI;
            if (jj2 < 0)
                jj2 += 2 * Math.PI;
            double mina = 0, maxa = 0;
            mina = jj1 < jj2 ? jj1 : jj2;
            maxa = jj1 > jj2 ? jj1 : jj2;
            angle = (maxa - mina) * 180 / Math.PI;
            if (angle > 180)
                angle = 360 - angle;
            return angle;
        }
        /// <summary>
        /// 获取两直线的交点
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static Point2 GetCross(Vector2 v1, Vector2 v2)
        {
            Line2 l1 = Line2.MakeLine(v1.Start, v1.End);
            Line2 l2 = Line2.MakeLine(v2.Start, v2.End);
            float a1 = l1.A;
            float b1 = l1.B;
            float c1 = l1.C;
            float a2 = l2.A;
            float b2 = l2.B;
            float c2 = l2.C;

            float x = (c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1);
            float y = (a1 * c2 - a2 * c1) / (a1 * b2 - a2 * b1);

            return new Point2(x, y);
        }
        /// <summary>
        /// 获取两直线的交点
        /// </summary>
        /// <param name="line1"></param>
        /// <param name="line2"></param>
        /// <returns></returns>
        public static Point2 GetCross(Line2 line1, Line2 line2)
        {
            Line2 l1 = Line2.MakeLine(line1.Point1, line1.Point2);
            Line2 l2 = Line2.MakeLine(line2.Point1, line2.Point2);
            float a1 = l1.A;
            float b1 = l1.B;
            float c1 = l1.C;
            float a2 = l2.A;
            float b2 = l2.B;
            float c2 = l2.C;
            //System.Windows.Forms.MessageBox.Show("a1=" + a1.ToString() + "\na2=" + a2.ToString());
            float x = -(c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1);
            float y = -(a1 * c2 - a2 * c1) / (a1 * b2 - a2 * b1);

            return new Point2(x, y);
        }
        /// <summary>
        /// 检测画角的4点条件是否符合。
        /// </summary>    
        /// <param name="x1">画角第一点X轴坐标</param>
        /// <param name="y1">画角第一点Y轴坐标</param>
        /// <param name="x2">画角第二点X轴坐标</param>
        /// <param name="y2">画角第二点Y轴坐标</param>
        /// <param name="x3">画角第三点X轴坐标</param>
        /// <param name="y3">画角第三点Y轴坐标</param>
        /// <param name="x4">画角第四点Y轴坐标</param>
        /// <param name="y4">画角第四点Y轴坐标</param>
        /// <returns>返回画角类型值</returns>
        public static int CeShi4DianHuaJiao(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4)
        {
            float k1, k2;
            k1 = 0;
            k2 = 0;
            if (x1 == x2)
            {
                if (x3 == x4)
                {
                    if (y1 == y2)
                    {
                        return 0;
                    }
                    else
                    {
                        if (y3 == y4)
                        {
                            return 0;
                        }
                        else
                        {
                            return 8;
                        }
                    }
                }
                else
                {
                    if (y1 == y2)
                    {
                        return 0;
                    }
                    else
                    {
                        if (y3 == y4)
                        {
                            return 1;
                        }
                        else
                        {
                            return 2;
                        }
                    }
                }
            }//x1!=x2
            else
            {
                if (x3 == x4)
                {
                    if (y3 == y4)
                    {
                        return 0;
                    }
                    else
                    {
                        if (y1 == y2)
                        {
                            return 3;
                        }
                        else
                        {
                            return 4;
                        }
                    }
                }
                else
                {
                    if (y1 == y2)
                    {
                        if (y3 == y4)
                        {
                            return 8;
                        }
                        else
                        {
                            return 5;
                        }

                    }
                    else
                    {
                        if (y3 == y4)
                        {
                            return 6;
                        }
                        else
                        {
                            k1 = (y1 - y2) / (x1 - x2);
                            k2 = (y3 - y4) / (x3 - x4);
                            if (k1 == k2)
                            {
                                return 8;
                            }
                            else
                            {
                                return 7;
                            }
                        }
                    }
                }
            }
        }
        /// <summary>
        /// 返回呼出交点的数组
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <param name="pf"></param>
        /// <returns></returns>
        public static bool HuChuJiaoDian(Vector2 v1, Vector2 v2, out Point2 pf)
        {
            pf = null;

            float x1 = v1.Start.X, y1 = v1.Start.Y, x2 = v1.End.X, y2 = v1.End.Y, x3 = v2.Start.X, y3 = v2.Start.Y, x4 = v2.End.X, y4 = v2.End.Y;
            float X1, Y1;
            float b1, b2, k1, k2;
            switch (CeShi4DianHuaJiao(x1, y1, x2, y2, x3, y3, x4, y4))
            {
                case 1:
                    X1 = x1;
                    Y1 = y3;
                    pf = new Point2(X1, Y1);
                    break;

                case 2:
                    k2 = (y3 - y4) / (x3 - x4);
                    b2 = y3 - k2 * x3;

                    X1 = x1;
                    Y1 = k2 * x1 + b2;
                    pf = new Point2(X1, Y1);
                    break;

                case 3:
                    X1 = x3;
                    Y1 = y1;
                    pf = new Point2(X1, Y1);

                    break;

                case 4:
                    k1 = (y1 - y2) / (x1 - x2);
                    b1 = y1 - k1 * x1;

                    X1 = x3;
                    Y1 = k1 * x3 + b1;
                    pf = new Point2(X1, Y1);

                    break;

                case 5:
                    Y1 = y1;
                    k2 = (y3 - y4) / (x3 - x4);
                    b2 = y3 - k2 * x3;
                    X1 = (Y1 - b2) / k2;
                    pf = new Point2(X1, Y1);

                    break;

                case 6:
                    Y1 = y3;
                    k1 = (y1 - y2) / (x1 - x2);
                    b1 = y1 - k1 * x1;
                    X1 = (Y1 - b1) / k1;
                    pf = new Point2(X1, Y1);
                    break;

                case 7:
                    k1 = (y1 - y2) / (x1 - x2);
                    k2 = (y3 - y4) / (x3 - x4);
                    b1 = y1 - k1 * x1;
                    b2 = y3 - k2 * x3;

                    X1 = (b2 - b1) / (k1 - k2);
                    Y1 = k1 * X1 + b1;
                    pf = new Point2(X1, Y1);
                    break;

                case 0:
                case 8:
                default:
                    break;
            }
            if (pf == null)
                return false;
            else
                return true;
        }
        /// <summary>
        /// 获取两直线的交点
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <param name="cross"></param>
        /// <returns></returns>
        public static bool GetCross(Vector2 v1, Vector2 v2, out Point2 cross)
        {
            if (HuChuJiaoDian(v1, v2, out cross))
                return true;
            else
                return false;
            /*
            cross = null;
            Point2 startP1 = v1.Start, startP2 = v2.Start;
            Line2 l1 = new Line2(v1.Start, v1.end);
            Line2 l2 = new Line2(v2.Start, v2.end);
            float xl1 = XL.Getxl2(l1), xl2 = XL.Getxl2(l2);
            //if (GetJJ(v1, v2) != 0)//夹角不能为零
            if (!Parallel(v1, v2))//两向量不能平行
            {
                float x = 0, y = 0;
                if (xl1 - xl2 != 0)//分母不能为零
                {
                    //当 y 相等时：
                    x = (xl1 * startP1.X - xl2 * startP2.X + startP2.Y - startP1.Y) / (xl1 - xl2);
                    //当 x 相等时：
                    y = (xl2 * startP1.Y - xl1 * startP2.Y + xl1 * xl2 * (startP2.X - startP1.X)) / (xl2 - xl1);
                }
                cross = new Point2(x, y);
                return true;
            }
            return false;*/
        }
        /// <summary>
        /// 获取两线段的交点，有问题，因为该算法是延长线的交点
        /// </summary>
        /// <param name="l1"></param>
        /// <param name="l2"></param>
        /// <param name="cross"></param>
        /// <returns></returns>
        public static bool GetCross(Line2 l1, Line2 l2, out Point2 cross)
        {
            cross = null;
            Point2 startP1 = l1.Point1, startP2 = l2.Point1;
            float xl1 = XL.Getxl2(l1), xl2 = XL.Getxl2(l2);
            Vector2 v1 = Vector2.LineToVector(l1);
            Vector2 v2 = Vector2.LineToVector(l2);
            //if (GetJJ(v1, v2) != 0)//夹角不能为零
            if (!Parallel(v1, v2))//两向量不能平行
            {
                float x = 0, y = 0;
                if (xl1 - xl2 != 0)//分母不能为零
                {
                    //当 y 相等时：
                    x = (xl1 * startP1.X - xl2 * startP2.X + startP2.Y - startP1.Y) / (xl1 - xl2);
                    //当 x 相等时：
                    y = (xl2 * startP1.Y - xl1 * startP2.Y + xl1 * xl2 * (startP2.X - startP1.X)) / (xl2 - xl1);
                }
                cross = new Point2(x, y);
                return true;
            }
            return false;
        }
        /// <summary>
        /// 两向量是否平行
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static bool Parallel(Vector2 v1, Vector2 v2)
        {
            Vector2 vv1 = Formatting(v1);
            Vector2 vv2 = Formatting(v2);
            if (vv2.End.X == 0 || vv2.End.Y == 0)
                return CompareEps(v1.Direction, v2.Direction);
            else
                return CompareEps(vv1.End.X / vv2.End.X, vv1.End.Y / vv2.End.Y);
        }
        /// <summary>
        /// 两向量是否垂直
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static bool Perpendicular(Vector2 v1, Vector2 v2)
        {
            Vector2 vv1 = Formatting(v1);
            Vector2 vv2 = Formatting(v2);
            return CompareEps(vv1.End.X * vv2.End.X + vv1.End.Y * vv2.End.Y, 0);
        }
        /// <summary>
        /// 平移(正方向)
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="offset"></param>
        /// <returns></returns>
        public static Vector2 Transform(Vector2 vector, Point2 offset)
        {
            Vector2 vv = new Vector2(new Point2(vector.start.X + offset.X, vector.start.Y + offset.Y), new Point2(vector.end.X + offset.X, vector.end.Y + offset.Y));
            return vv;
        }
        /// <summary>
        /// 旋转(顺时针)
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="angle"></param>
        /// <returns></returns>
        public static Vector2 Transform(Vector2 vector, double angle)
        {
            Point2 p1 = vector.Start;
            vector = Vector2.Formatting(vector);
            Point2 rotationP = Rotation(vector.End, angle);
            Point2 p2 = new Point2(rotationP.X + p1.X, rotationP.Y + p1.Y);
            Vector2 vv = new Vector2(p1, p2);
            return vv;
        }
        /// <summary>
        /// 求点p经过旋转弧度a之后的新点坐标
        /// </summary>
        /// <param name="p"></param>
        /// <param name="a"></param>
        /// <returns></returns>
        public static Point2 Rotation(Point2 p, double a)
        {
            Point2 pp = new Point2();
            double[,] A = new double[1, 2];
            double[,] B = new double[2, 2];
            B[0, 0] = Math.Cos(a);
            B[0, 1] = Math.Sin(a);
            B[1, 0] = 0 - Math.Sin(a);
            B[1, 1] = Math.Cos(a);
            A[0, 0] = p.X;
            A[0, 1] = p.Y;
            A = Mul.MatrixMul(A, B);
            pp.X = (float)A[0, 0];
            pp.Y = (float)A[0, 1];
            return pp;
        }
        /// <summary>
        /// 综合转换
        /// </summary>
        /// <param name="srcPt"></param>
        /// <param name="tX"></param>
        /// <param name="tY"></param>
        /// <param name="sX"></param>
        /// <param name="sY"></param>
        /// <param name="Theta">以弧度计量的角度</param>
        /// <returns></returns>
        public static Point2 GetTransformation(Point2 srcPt, double tX, double tY, double sX, double sY, double Theta)
        {
            //(tX,tY) describes the translation
            //(sX,sY) describes the scale
            //Theta describes the rotation
            //SrcPt is the point to be transformed
            //[RETURN] is the transformed point

            //   The general formulas:
            //   X' = xSxCosq - ySySinq + dx
            //   Y' = xSxSinq + ySyCosq + dy

            double Cosq = Math.Cos(Theta);
            double Sinq = Math.Sin(Theta);
            Point2 retPt = new Point2();
            retPt.X = (float)((srcPt.X * sX * Cosq) - (srcPt.Y * sY * Sinq) + tX);
            retPt.Y = (float)((srcPt.X * sX * Sinq) + (srcPt.Y * sY * Cosq) + tY);
            return retPt;
        }
        /// <summary>
        /// 基于指定基矩阵的向量
        /// </summary>
        /// <param name="vector">原向量</param>
        /// <param name="baseMatrix">基矩阵</param>
        /// <returns></returns>
        public static Vector2 BaseToSystem(Vector2 vector, double[,] baseMatrix)
        {
            Vector2 v = null;
            if (baseMatrix.GetLength(0) == 2 && baseMatrix.GetLength(1) == 2)
            {
                v = Mul.Vector2MulMatrix(vector, baseMatrix);
            }
            return v;
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
        /// 求逆阵
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static double[,] GetInverseMatrix(double[,] a)
        {
            double[,] b = new double[2, 2];
            double t = 0;
            //b[0] = new double[2];
            //b[1] = new double[2];
            //b[2] = new double[2];
            //求伴随阵
            for (int ci = 0; ci < 2; ci++)
            {
                for (int cj = 0; cj < 2; cj++)
                {
                    //求余子式
                    for (int i = 0; i < 2; i++)
                    {
                        for (int j = 0; j < 2; j++)
                        {
                            if (i != ci && j != cj)
                            {
                                t = a[i, j];
                            }
                        }
                    }
                    b[cj, ci] = IsNeg(ci + cj) * t;
                }
            }
            //求矩阵的行列式值
            double rValue = 0;
            for (int i = 0; i < 2; i++)
            {
                rValue = rValue + a[i, 0] * b[0, i];
            }
            if (rValue != 0) //行列式值不为0，才存在逆阵 
            {
                //数乘
                rValue = 1 / rValue;
                for (int i = 0; i < 2; i++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        b[i, j] = b[i, j] * rValue;
                    }
                }
            }
            return b;
        }
        /// <summary>
        /// 线转化为向量
        /// </summary>
        /// <param name="line"></param>
        /// <returns></returns>
        public static Vector2 LineToVector(Line2 line)
        {
            Vector2 vector;
            vector = new Vector2(new Point2(line.Point2.X - line.Point1.X, line.Point2.Y - line.Point1.Y));
            return vector;
        }
        /// <summary>
        /// 向量转化为线
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static Line2 VectorToLine(Vector2 vector)
        {
            Line2 line;
            line = new Line2(vector.Start, vector.End);
            return line;
        }
        /// <summary>
        /// 求垂直平分线，第一个点Point1即为垂足，且长度为1
        /// </summary>
        /// <param name="line"></param>
        /// <returns></returns>
        public static Line2 GetCZAVGLine(Line2 line)
        {
            Line2 l = null;
            float x1 = line.Point1.X;
            float x2 = line.Point2.X;
            float y1 = line.Point1.Y;
            float y2 = line.Point2.Y;
            float midX1 = (x1 + x2) / 2;
            float midY1 = (y1 + y2) / 2;
            Point2 p1 = new Point2(midX1, midY1);
            //bool inLine;
            //p1为垂足，不能作为线外的点
            //l = GetCZLine( , line, out inLine);


            /*
            float xl = XL.Getxl2(line);
            float xl1 = -1 / xl;
            float x = (float)Math.Sqrt(1 / (1 + xl1 * xl1));
            float y = xl1 * x;
            Vector2 v = new Vector2(new Point2(x, y));
            v = Vector2.Transform(v, p1);
            Point2 p2 = GetPointOfLengthFormStartP(v, 1);
            l = new Line2(p1, p2);
            */
            return l;
        }
        /// <summary>
        /// 根据线段的始点、角度和长度，确定线段的正向终点
        /// </summary>
        /// <param name="startP"></param>
        /// <param name="angle">弧度为单位</param>
        /// <param name="length"></param>
        /// <returns></returns>
        public static Point2 GetEndPofLineK(Point2 startP, double angle, double length)
        {
            float k = (float)Math.Tan(angle);
            double d2 = length * length;
            float x = startP.X + (float)Math.Sqrt(d2 / (1 + k * k));
            float y = startP.Y - (float)Math.Sqrt(d2 - (d2 / (1 + k * k))) * Math.Sign(k);
            return new Point2(x, y);
        }
        /// <summary>
        /// 根据线段的始点、角度和长度，确定线段的反向终点
        /// </summary>
        /// <param name="startP"></param>
        /// <param name="angle">弧度为单位</param>
        /// <param name="length"></param>
        /// <returns></returns>
        public static Point2 GetEndPofLine_K(Point2 startP, double angle, double length)
        {
            float k = (float)Math.Tan(angle);
            double d2 = length * length;
            float x = startP.X - (float)Math.Sqrt(d2 / (1 + k * k));
            float y = startP.Y + (float)Math.Sqrt(d2 - (d2 / (1 + k * k))) * Math.Sign(k);
            return new Point2(x, y);
        }
        /// <summary>
        /// 求过点p垂直于line的垂直线，第一个点Point1即为垂足，长度为1 注意：这个函数非常重要！很多地方都用上！
        /// </summary>
        /// <param name="p"></param>
        /// <param name="line"></param>
        /// <param name="ifInLine"></param>
        /// <returns></returns>
        public static Line2 GetCZLine(Point2 p, Line2 line, out bool ifInLine)
        {
            float xl = -1/(line.K);
            Line2 l = new Line2(p, xl, 1);
            Point2 cross;
            LLcross.GetLLcrs2(line, l, out cross);
            double angle = Vector2.LineToVector(line).Direction + 90;
            Point2 q = GetEndPofLineK(line.Point1, angle, 1);
            ifInLine = InLine(line.Point1, line.Point2, p);
            return new Line2(cross, XL.Getxl2(l), 1);
            /*Line2 l = null;
            ifInLine = false;
            float xl = XL.Getxl2(line);
            l = new Line2(p, xl, 1);
            Point2 cross;
            LLcross.GetLLcrs2(line, l, out cross);
            if (cross != null)
                ifInLine = InLine(line.Point1, line.Point2, cross);
            return l;*/
        }
        /// <summary>
        /// X是否在AB所在的直线上
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="X"></param>
        /// <returns></returns>
        private static bool InLine(Point2 A, Point2 B, Point2 X)
        {
            if ((A.X - X.X) * (X.X - B.X) >= 0 && (A.Y - X.Y) * (X.Y - B.Y) >= 0)
                return true;
            return false;
        }
        /// <summary>
        /// 求角平分线(中分线)，v1,v2不能单位化，否则有误。
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <param name="avgLine"></param>
        public static void GetAngleAvgLine(Vector2 v1, Vector2 v2, out Vector2 avgLine)
        {
            //根据角度相等有
            //(v.x * v1.x + v.y * v1.y) / |v1| = (v.x * v2.x + v.y * v2.y) / |v2|
            //于是有
            //v.x = -|v1|v2.y + |v2|v1.y
            //v.y = |v1|v2.x - |v2|v1.x
            Point2 cross;
            Vector2.GetCross(v1, v2, out cross);
            cross = new Point2(-cross.X, -cross.Y);
            v1 = Vector2.Transform(v1, cross);
            v1 = Vector2.Transform(v2, cross);
            float x = -(float)Vector2.GetLength(v1) * v2.End.Y + (float)Vector2.GetLength(v2) * v1.End.Y;
            float y = (float)Vector2.GetLength(v1) * v2.End.X - (float)Vector2.GetLength(v2) * v1.End.X;
            avgLine = new Vector2(new Point2(x, y));
            cross = new Point2(-cross.X, -cross.Y);
            avgLine = Vector2.Transform(avgLine, cross);
        }
        /// <summary>
        /// 求角平分线(中分线)
        /// </summary>
        /// <param name="l1"></param>
        /// <param name="l2"></param>
        /// <param name="avgLine"></param>
        public static void GetAngleAvgLine(Line2 l1, Line2 l2, out Line2 avgLine)
        {
            //根据角度相等有(v.x * v1.x + v.y * v1.y) / |v1| = (v.x * v2.x + v.y * v2.y) / |v2|
            //于是有
            //v.x = -|v1|v2.y + |v2|v1.y
            //v.y = |v1|v2.x - |v2|v1.x
            Vector2 v1 = Vector2.LineToVector(l1);
            Vector2 v2 = Vector2.LineToVector(l2);
            Point2 cross;
            Vector2.GetCross(v1, v2, out cross);
            cross = new Point2(-cross.X, -cross.Y);
            v1 = Vector2.Transform(v1, cross);
            v1 = Vector2.Transform(v2, cross);
            float x = -(float)Vector2.GetLength(v1) * v2.End.Y + (float)Vector2.GetLength(v2) * v1.End.Y;
            float y = (float)Vector2.GetLength(v1) * v2.End.X - (float)Vector2.GetLength(v2) * v1.End.X;
            Vector2 avgVector = new Vector2(new Point2(x, y));
            cross = new Point2(-cross.X, -cross.Y);
            avgVector = Vector2.Transform(avgVector, cross);
            avgLine = Vector2.VectorToLine(avgVector);
        }
    }
    /// <summary>
    /// 三维向量
    /// </summary>
    public class Vector3
    {
        /// <summary>
        /// 判断在给定的误差范围内两个数值是否相等
        /// </summary>
        /// <param name="comp1"></param>
        /// <param name="comp2"></param>
        /// <returns></returns>
        public static bool CompareEps(double comp1, double comp2)
        {
            return Math.Abs(comp1 - comp2) <= Common.Eps;
        }
        double[] dir = new double[2];
        Point3 start = new Point3();
        Point3 end = new Point3();
        Point3 x0point = new Point3();
        /// <summary>
        /// 构造函数
        /// </summary>
        public Vector3()
        {

        }
        /// <summary>
        /// 始点为原点，终点为p的向量
        /// </summary>
        /// <param name="p"></param>
        public Vector3(Point3 p)
        {
            start = new Point3(0, 0, 0);
            end = p;
        }
        /// <summary>
        /// 如果pointDirection为true，则start = new Point3(0, 0, 0);end = startP;x0point = endP;
        /// </summary>
        /// <param name="startP"></param>
        /// <param name="endP"></param>
        /// <param name="pointDirection"></param>
        public Vector3(Point3 startP, Point3 endP, bool pointDirection)
        {
            if (pointDirection)
            {
                start = new Point3(0, 0, 0);
                end = startP;
                x0point = endP;
            }
            else
            {
                start = startP;
                end = endP;
            }
        }
        /// <summary>
        /// 加法操作符
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        /// <returns></returns>
        public static Vector3 operator +(Vector3 vector1, Vector3 vector2)
        {
            Vector3 vector = new Vector3(vector1.Start, new Point3(vector2.End.X + vector1.End.X - vector1.Start.X, vector2.End.Y + vector1.End.Y - vector1.Start.Y, vector2.End.Z + vector1.End.Z - vector1.Start.Z), false);
            return vector;
        }
        /// <summary>
        /// 减法操作符
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        /// <returns></returns>
        public static Vector3 operator -(Vector3 vector1, Vector3 vector2)
        {
            Vector3 vector = new Vector3(vector2.End, vector1.End, false);
            return vector;
        }
        /// <summary>
        /// 规范化
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static Vector3 Formatting(Vector3 vector)
        {
            Point3 offset = vector.start;
            Vector3 vv = new Vector3(new Point3(vector.end.X - offset.X, vector.end.Y - offset.Y, vector.end.Z - offset.Z));
            return vv;
        }
        /// <summary>
        /// 是否相等
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        /// <returns></returns>
        public static bool IsEqual(Vector3 vector1, Vector3 vector2)
        {
            Vector3 v1 = Vector3.Formatting(vector1);
            Vector3 v2 = Vector3.Formatting(vector2);
            return Point3.IsEqual(v1.start, v2.start) && Point3.IsEqual(v1.end, v2.end);
        }
        /// <summary>
        /// 方向
        /// </summary>
        public double[] Direction
        {
            get
            {
                dir = VectorDirection.GetDirection3(this);
                return dir;
            }
        }
        /// <summary>
        /// 获取或设置向量的始点
        /// </summary>
        public Point3 Start
        {
            get
            {
                return start;
            }
            set
            {
                start = value;
            }
        }
        /// <summary>
        /// 获取或设置向量的终点
        /// </summary>
        public Point3 End
        {
            get
            {
                return end;
            }
            set
            {
                end = value;
            }
        }
        /// <summary>
        /// 向量上的某一点
        /// </summary>
        public Point3 X0point
        {
            get
            {
                return x0point;
            }
            set
            {
                x0point = value;
            }
        }
        /// <summary>
        /// 空向量
        /// </summary>
        public static Vector3 Empty
        {
            get
            {
                return new Vector3(new Point3(0, 0, 0));
            }
        }
        /// <summary>
        /// 将向量转化为可视化字符串
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return start.X.ToString() + "," + start.Y.ToString() + "," + start.Z.ToString() + " ----> " + end.X.ToString() + "," + end.Y.ToString() + "," + end.Z.ToString();
        }
        /// <summary>
        /// 长度(模)
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static double GetLength(Vector3 vector)
        {
            Point3 startP = vector.Start;
            Point3 endP = vector.End;
            return Math.Sqrt((endP.X - startP.X) * (endP.X - startP.X) + (endP.Y - startP.Y) * (endP.Y - startP.Y) + (endP.Z - startP.Z) * (endP.Z - startP.Z));
        }
        /// <summary>
        /// 获取直线上距离开始点len的点，可能有问题...
        /// </summary>
        /// <param name="vector">直线的向量</param>
        /// <param name="len">距离</param>
        /// <returns></returns>
        public static Point3 GetPointOfLengthFormStartP(Vector3 vector, double len)
        {
            Point3 startP = vector.Start;
            float k1 = (float)Math.Tan(vector.Direction[0]);
            float k2 = (float)Math.Tan(vector.Direction[1]);
            double d2 = len * len;
            float x = startP.X - (float)Math.Sqrt(d2 / (1 + k1 * k1));
            float y = startP.Y + (float)Math.Sqrt(d2 - (d2 / (1 + k1 * k1))) * Math.Sign(k1);
            float z = startP.Z + (float)Math.Sqrt(d2 - (d2 / (1 + k2 * k2))) * Math.Sign(k2);
            return new Point3(x, y, z);
            /*
            double length = Vector3.GetLength(vector);
            if (length > 0)
            {
                if (len < 0)//当len < 0 时，即在向量vector的反方向。
                {
                    float x = 0;
                    float y = 0;
                    if (vector.End.X - vector.Start.X == 0)
                    {
                        y = (float)(Math.Abs(vector.End.Y - vector.Start.Y) * (-len) / length);//在Y轴上的投影长，也可以用Vector3.GetProjL2L();//保证为正！
                        if (vector.Direction[0] != 0)
                            x = y / (float)vector.Direction[0];
                        else
                            x = float.MaxValue;
                    }
                    else
                    {
                        x = (float)(Math.Abs(vector.End.X - vector.Start.X) * (-len) / length);//在X轴上的投影长，也可以用Vector3.GetProjL2L();//保证为正！
                        y = (float)vector.Direction[0] * x;
                    }
                    float z = (float)vector.Direction[1] * x;
                    return new Point3(vector.Start.X - x, vector.Start.Y - y, vector.Start.Z - z);
                }
                else if (len == 0)
                {
                    return vector.Start;
                }
                else
                {
                    float x = 0;
                    float y = 0;
                    if (vector.End.X - vector.Start.X == 0)
                    {
                        y = (float)(Math.Abs(vector.End.Y - vector.Start.Y) * len / length);//在Y轴上的投影长，也可以用Vector3.GetProjL2L();//保证为正！
                        if (vector.Direction[0] != 0)
                            x = y / (float)vector.Direction[0];
                        else
                            x = float.MaxValue;
                    }
                    else
                    {
                        x = (float)(Math.Abs(vector.End.X - vector.Start.X) * len / length);//在X轴上的投影长，也可以用Vector3.GetProjL2L();//保证为正！
                        y = (float)vector.Direction[0] * x;
                    }
                    float z = (float)vector.Direction[1] * x;
                    return new Point3(vector.Start.X + x, vector.Start.Y + y, vector.Start.Z + z);
                }
            }
            else
            {
                return vector.Start;
            }*/
        }
        /// <summary>
        /// 缩放
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="scalingFactor"></param>
        /// <returns></returns>
        public static Vector3 Scale(Vector3 vector, double scalingFactor)
        {
            double len = Vector3.GetLength(vector);
            len *= scalingFactor;
            Point3 newEnd = GetPointOfLengthFormStartP(vector, len);
            return new Vector3(vector.Start, newEnd, false);
        }
        /// <summary>
        /// 单位向量
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static Vector3 GetUnitVector(Vector3 vector)
        {
            vector = Vector3.Formatting(vector);
            float x = vector.End.X;
            float y = vector.End.Y;
            float z = vector.End.Z;
            double length = Math.Sqrt(x * x + y * y + z * z);
            x = (float)(x / length);
            y = (float)(y / length);
            z = (float)(z / length);
            Vector3 unit = new Vector3(new Point3(x, y, z));
            return unit;
        }
        /// <summary>
        /// 获取当前向量与X轴的夹角
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static double GetAngleX(Vector3 vector)
        {
            Vector3 unit = Vector3.GetUnitVector(vector);
            return Math.Acos(unit.End.X);
        }
        /// <summary>
        /// 获取当前向量与Y轴的夹角
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static double GetAngleY(Vector3 vector)
        {
            Vector3 unit = Vector3.GetUnitVector(vector);
            return Math.Acos(unit.End.Y);
        }
        /// <summary>
        /// 获取当前向量与Z轴的夹角
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static double GetAngleZ(Vector3 vector)
        {
            Vector3 unit = Vector3.GetUnitVector(vector);
            return Math.Acos(unit.End.Z);
        }
        /// <summary>
        /// 两向量的夹角
        /// </summary>
        /// <param name="vector1"></param>
        /// <param name="vector2"></param>
        /// <returns></returns>
        public static double GetVectorsAngle(Vector3 vector1, Vector3 vector2)
        {
            Vector3 unit1 = Vector3.GetUnitVector(vector1);
            Vector3 unit2 = Vector3.GetUnitVector(vector2);
            return Math.Acos(unit1.End.X * unit2.End.X + unit1.End.Y * unit2.End.Y + unit1.End.Z * unit2.End.Z);
        }
        /// <summary>
        /// 3维转化为2维
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static Vector2 Vector3ToVector2(Vector3 vector)
        {
            Vector3 n = new Vector3(new Point3(0, 0, 1));
            float k = Mul.GetPmul(vector, n);
            Vector3 smV = new Vector3(new Point3(k * n.End.X, k * n.End.Y, k * n.End.Z));
            vector = vector - smV;
            Vector2 vector2 = new Vector2(new Point2(vector.End.X, vector.End.Y));
            return vector2;
        }/*
        /// <summary>
        /// 两向量的夹角[0,180)
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static double GetJJ(Vector3 v1, Vector3 v2)
        {
            Vector3 vv1 = Formatting(v1);
            Vector3 vv2 = Formatting(v2);
            CSharpAlgorithm.Algorithm.Matrix x1 = Mul.VectorToMatrix(vv1);
            CSharpAlgorithm.Algorithm.Matrix x2 = Mul.VectorToMatrix(vv2);
            CSharpAlgorithm.Algorithm.Matrix xMatrix = Mul.GetXmul(x1, x2);
            Vector3 baseVector = Mul.MatrixToVector(xMatrix);
            Vector2 vvv1 = Vector3ToVector2(vv1, baseVector);
            Vector2 vvv2 = Vector3ToVector2(vv2, baseVector);
            return Vector2.GetJJ(vvv1, vvv2);
        }*/
        /// <summary>
        /// 两向量的夹角[0,180)
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static double GetJJ(Vector3 v1, Vector3 v2)
        {
            Vector3 vv1 = GetUnitVector(v1);
            Vector3 vv2 = GetUnitVector(v2);
            return Math.Acos(Mul.GetPmul(vv1, vv2) / (GetLength(vv1) * GetLength(vv2))) * 180 / Math.PI;
        }
        /// <summary>
        /// 判断两线段是否相交，v1,v2不能格式化，否则有误
        /// </summary>
        /// <param name="l1"></param>
        /// <param name="l2"></param>
        /// <returns></returns>
        public static bool IsLineCross(Line3 l1, Line3 l2)
        {
            double x1 = l1.Point1.X;
            double y1 = l1.Point1.Y;
            double z1 = l1.Point1.Z;
            double x2 = l1.Point2.X;
            double y2 = l1.Point2.Y;
            double z2 = l1.Point2.Z;
            double x3 = l2.Point1.X;
            double y3 = l2.Point1.Y;
            double z3 = l2.Point1.Z;
            double x4 = l2.Point2.X;
            double y4 = l2.Point2.Y;
            double z4 = l2.Point2.Z;
            //t1,t2对应为l1,l2方程的参数， t1,t2范围必须为[0, 1] 且满足 z1+(z2-z1)t1=z3+(z4-z3)t2，否则线段不相交
            double t2 = (x2 * y3 - x1 * y3 - x2 * y1 - x3 * y2 + x1 * y2 + x3 * y1) / (x4 * y2 - x3 * y2 - x4 * y1 + x3 * y1 - x2 * y4 + x1 * y4 + x2 * y3 - x1 * y3);
            double t1 = (x3 + (x4 - x3) * t2 - x1) / (x2 - x1);
            bool inRange1 = Between(t1, 0, 1);
            bool inRange2 = Between(t2, 0, 1);
            if (inRange1 && inRange2)
            {
                if (z1 + (z2 - z1) * t1 == z3 + (z4 - z3) * t2)
                {
                    return true;
                }
            }
            return false;
        }
        /// <summary>
        /// 判断两直线是否相交，v1,v2不能格式化，否则有误
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static bool IsVectorCross(Vector3 v1, Vector3 v2)
        {
            double x1 = v1.Start.X;
            double y1 = v1.Start.Y;
            double z1 = v1.Start.Z;
            double x2 = v1.End.X;
            double y2 = v1.End.Y;
            double z2 = v1.End.Z;
            double x3 = v2.Start.X;
            double y3 = v2.Start.Y;
            double z3 = v2.Start.Z;
            double x4 = v2.End.X;
            double y4 = v2.End.Y;
            double z4 = v2.End.Z;
            //t1,t2对应为v1,v2方程的参数， 满足 z1+(z2-z1)t1=z3+(z4-z3)t2，否则直线不相交
            double t2 = (x2 * y3 - x1 * y3 - x2 * y1 - x3 * y2 + x1 * y2 + x3 * y1) / (x4 * y2 - x3 * y2 - x4 * y1 + x3 * y1 - x2 * y4 + x1 * y4 + x2 * y3 - x1 * y3);
            double t1 = (x3 + (x4 - x3) * t2 - x1) / (x2 - x1);
            if (z1 + (z2 - z1) * t1 == z3 + (z4 - z3) * t2)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// 范围
        /// </summary>
        /// <param name="val"></param>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <returns></returns>
        public static bool Between(double val, double min, double max)
        {
            return val >= min && val <= max;
        }
        /// <summary>
        /// 由直线的端点获取直线的方程组系数矩阵以及常数矩阵，注：coef[4,3], cons[4]
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <param name="coef"></param>
        /// <param name="cons"></param>
        public static void GetCoefAndConsFrom2Vectors(Vector3 v1, Vector3 v2, ref double[,] coef, ref double[] cons)
        {
            double x10 = v1.Start.X;
            double x20 = v2.Start.X;
            double y10 = v1.Start.Y;
            double y20 = v2.Start.Y;
            double z10 = v1.Start.Z;
            double z20 = v2.Start.Z;
            double x11 = v1.End.X;
            double x21 = v2.End.X;
            double y11 = v1.End.Y;
            double y21 = v2.End.Y;
            double z11 = v1.End.Z;
            double z21 = v2.End.Z;
            coef = new double[4, 3];
            cons = new double[4];
            //coef={{y11-y10,x10-x11,0},{z11-z10,0,x11-x10},{y21-y20,x20-x21,0},{z21-z20,0,x21-x20}};
            //cons={(y11-y10)*x10-(x11-x10)*y10,(z11-z10)*x10-(x11-x10)*z10,(y21-y20)*x20-(x21-x20)*y20,(z21-z20)*x20-(x21-x20)*z20};
            coef[0, 0] = y11 - y10;
            coef[0, 1] = x10 - x11;
            coef[0, 2] = 0;
            coef[1, 0] = z11 - z10;
            coef[1, 1] = 0;
            coef[1, 2] = x11 - x10;
            coef[2, 0] = y21 - y20;
            coef[2, 1] = x20 - x21;
            coef[2, 2] = 0;
            coef[3, 0] = z21 - z20;
            coef[3, 1] = 0;
            coef[3, 2] = x21 - x20;
            cons[0] = (y11 - y10) * x10 - (x11 - x10) * y10;
            cons[1] = (z11 - z10) * x10 - (x11 - x10) * z10;
            cons[2] = (y21 - y20) * x20 - (x21 - x20) * y20;
            cons[3] = (z21 - z20) * x20 - (x21 - x20) * z20;
        }
        /// <summary>
        /// 获取两直线的交点，v1,v2不能格式化，否则有误
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static Point3 GetCross(Vector3 v1, Vector3 v2)
        {
            Point3 cross = null;
            if (!Parallel(v1, v2))//两向量不能平行
            {
                if (IsVectorCross(v1, v2))
                {
                    /*double x10 = v1.Start.X;
                    double x20 = v2.Start.X;
                    double y10 = v1.Start.Y;
                    double y20 = v2.Start.Y;
                    double z10 = v1.Start.Z;
                    double z20 = v2.Start.Z;
                    double x11 = v1.End.X;
                    double x21 = v2.End.X;
                    double y11 = v1.End.Y;
                    double y21 = v2.End.Y;
                    double z11 = v1.End.Z;
                    double z21 = v2.End.Z;
                    double[] coef ={ };
                    double[] cons ={ };
                    CSharpAlgorithm.Algorithm.Matrix mtxCoef = new CSharpAlgorithm.Algorithm.Matrix(3, coef);
                    CSharpAlgorithm.Algorithm.Matrix mtxCons = new CSharpAlgorithm.Algorithm.Matrix(3, 1, cons);
                    CSharpAlgorithm.Algorithm.Matrix root = new CSharpAlgorithm.Algorithm.Matrix();//(3, 1);
                    CSharpAlgorithm.Algorithm.LEquations le = new CSharpAlgorithm.Algorithm.LEquations();
                    le.Init(mtxCoef, mtxCons);
                    if (le.GetRootsetGauss(root))
                    {};*/
                    //coef, cons 待定。
                    double[,] coef = null;// = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
                    double[] cons = null;// = { 1, 2, 3 };
                    GetCoefAndConsFrom2Vectors(v1, v2, ref coef, ref cons);//即如何由直线的两个端点求取其解析方程式，以便求直线的交点（根）
                    //CSharpAlgorithm.Algorithm.Matrix mtx = new CSharpAlgorithm.Algorithm.Matrix(coef);
                    //CSharpAlgorithm.Algorithm.Matrix mtx1 = new CSharpAlgorithm.Algorithm.Matrix(4, 1, cons);
                    //System.Windows.Forms.MessageBox.Show(mtx.ToString() + "\n" + mtx1.ToString());
                    //double[] root = GetRootsOfCramer3D(coef, cons);
                    double[] root = GetRootsByCramer(coef, cons);
                    float x = (float)root[0];
                    float y = (float)root[1];
                    float z = (float)root[2];
                    cross = new Point3(x, y, z);
                }
            }
            return cross;
        }
        /// <summary>
        /// 获取两线段的交点
        /// </summary>
        /// <param name="l1"></param>
        /// <param name="l2"></param>
        /// <returns></returns>
        public static Point3 GetCross(Line3 l1, Line3 l2)
        {
            Point3 cross = null;
            Vector3 v1 = Vector3.LineToVector(l1);
            Vector3 v2 = Vector3.LineToVector(l2);
            if (!Parallel(v1, v2))//两向量不能平行
            {
                if (IsVectorCross(v1, v2))
                {
                    /*double x10 = v1.Start.X;
                    double x20 = v2.Start.X;
                    double y10 = v1.Start.Y;
                    double y20 = v2.Start.Y;
                    double z10 = v1.Start.Z;
                    double z20 = v2.Start.Z;
                    double x11 = v1.End.X;
                    double x21 = v2.End.X;
                    double y11 = v1.End.Y;
                    double y21 = v2.End.Y;
                    double z11 = v1.End.Z;
                    double z21 = v2.End.Z;
                    double[] coef ={ };
                    double[] cons ={ };
                    CSharpAlgorithm.Algorithm.Matrix mtxCoef = new CSharpAlgorithm.Algorithm.Matrix(3, coef);
                    CSharpAlgorithm.Algorithm.Matrix mtxCons = new CSharpAlgorithm.Algorithm.Matrix(3, 1, cons);
                    CSharpAlgorithm.Algorithm.Matrix root = new CSharpAlgorithm.Algorithm.Matrix();//(3, 1);
                    CSharpAlgorithm.Algorithm.LEquations le = new CSharpAlgorithm.Algorithm.LEquations();
                    le.Init(mtxCoef, mtxCons);
                    if (le.GetRootsetGauss(root))
                    {};*/
                    //coef, cons 待定。
                    double[,] coef = null;// = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
                    double[] cons = null;// = { 1, 2, 3 };
                    GetCoefAndConsFrom2Vectors(v1, v2, ref coef, ref cons);//即如何由直线的两个端点求取其解析方程式，以便求直线的交点（根）
                    //CSharpAlgorithm.Algorithm.Matrix mtx = new CSharpAlgorithm.Algorithm.Matrix(coef);
                    //CSharpAlgorithm.Algorithm.Matrix mtx1 = new CSharpAlgorithm.Algorithm.Matrix(4, 1, cons);
                    //System.Windows.Forms.MessageBox.Show(mtx.ToString() + "\n" + mtx1.ToString());
                    //double[] root = GetRootsOfCramer3D(coef, cons);
                    double[] root = GetRootsByCramer(coef, cons);
                    float x = (float)root[0];
                    float y = (float)root[1];
                    float z = (float)root[2];
                    cross = new Point3(x, y, z);
                }
            }
            return cross;
        }
        /// <summary>
        /// 两向量是否平行
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static bool Parallel(Vector3 v1, Vector3 v2)
        {
            Vector3 vv1 = Formatting(v1);
            Vector3 vv2 = Formatting(v2);
            if (vv2.End.X == 0 || vv2.End.Y == 0 || vv2.End.Z == 0)
                return CompareEps(v1.Direction[0], v2.Direction[0]) && CompareEps(v1.Direction[1], v2.Direction[1]);
            else
                return CompareEps(vv1.End.X / vv2.End.X, vv1.End.Y / vv2.End.Y) && CompareEps(vv1.End.Y / vv2.End.Y, vv1.End.Z / vv2.End.Z);
        }
        /// <summary>
        /// 两向量是否垂直
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static bool Perpendicular(Vector3 v1, Vector3 v2)
        {
            Vector3 vv1 = Formatting(v1);
            Vector3 vv2 = Formatting(v2);
            return CompareEps(vv1.End.X * vv2.End.X + vv1.End.Y * vv2.End.Y + vv1.End.Z * vv2.End.Z, 0);
        }
        /// <summary>
        /// 平移(正方向)
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="offset"></param>
        /// <returns></returns>
        public static Vector3 Transform(Vector3 vector, Point3 offset)
        {
            Vector3 vv = new Vector3(new Point3(vector.Start.X + offset.X, vector.Start.Y + offset.Y, vector.Start.Z + offset.Z), new Point3(vector.End.X + offset.X, vector.End.Y + offset.Y, vector.End.Z + offset.Z), false);
            return vv;
        }
        /// <summary>
        /// 旋转(顺时针)
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="baseVector"></param>
        /// <returns></returns>
        public static Vector3 Transform(Vector3 vector, Vector3 baseVector)
        {
            Vector3 N = Vector3.GetUnitVector(baseVector);
            Vector3 U = new Vector3();
            U.Start = new Point3(0, 0, 0);
            double invLength = 0;
            if (Math.Abs(vector.End.X) >= Math.Abs(vector.End.Y))
            {
                invLength = 1 / Math.Sqrt(N.End.X * N.End.X + N.End.Z * N.End.Z);
                U.End.X = (float)(N.End.Z * invLength);
                U.End.Y = 0;
                U.End.Z = -(float)(N.End.X * invLength);
            }
            else
            {
                invLength = 1 / Math.Sqrt(N.End.Y * N.End.Y + N.End.Z * N.End.Z);
                U.End.X = 0;
                U.End.Y = (float)(N.End.Z * invLength);
                U.End.Z = -(float)(N.End.Y * invLength);
            }
            Vector3 V = Mul.GetXmul(N, U);
            Point3 startP = vector.Start;
            vector = Vector3.Formatting(vector);
            float x = Mul.GetPmul(vector, U);
            float y = Mul.GetPmul(vector, V);
            float z = Mul.GetPmul(vector, N);
            Vector3 vv = new Vector3(startP, new Point3(x + startP.X, y + startP.Y, z + startP.Z), false);
            return vv;
        }
        /// <summary>
        /// 两面交线
        /// </summary>
        /// <param name="A1"></param>
        /// <param name="B1"></param>
        /// <param name="C1"></param>
        /// <param name="D1"></param>
        /// <param name="A2"></param>
        /// <param name="B2"></param>
        /// <param name="C2"></param>
        /// <param name="D2"></param>
        /// <returns></returns>
        public static Vector3 GetVectorFrom2Faces(double A1, double B1, double C1, double D1, double A2, double B2, double C2, double D2)
        {
            Vector3 direction = Mul.GetXmul(new Vector3(new Point3((float)A1, (float)B1, (float)C1)), new Vector3(new Point3((float)A2, (float)B2, (float)C2)));
            float y = (float)((B1 * C1 * D1 - B1 * C2 * D1) / (B1 * B1 * C2 - B1 * B2 * C1));
            float z = (float)((B1 * D2 - B2 - D1) / (B2 * C1 - B1 * C2));
            Vector3 v = new Vector3();
            if (Vector3.Parallel(direction, direction))
            {
                Point3 xpoint = new Point3(0, y, z);
                v.End = direction.End;
                v.X0point = xpoint;
            }
            else
            {
                v.Start = new Point3(0, 0, 0);
                v.End = direction.End;
            }
            return v;
        }
        /// <summary>
        /// 线转化为向量
        /// </summary>
        /// <param name="line"></param>
        /// <returns></returns>
        public static Vector3 LineToVector(Line3 line)
        {
            Vector3 vector;
            vector = new Vector3(new Point3(line.Point2.X - line.Point1.X, line.Point2.Y - line.Point1.Y, line.Point2.Z - line.Point1.Z));
            return vector;
        }
        /// <summary>
        /// 向量转化为线
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static Line3 VectorToLine(Vector3 vector)
        {
            Line3 line;
            line = new Line3(vector.Start, vector.End);
            return line;
        }
        /// <summary>
        /// 点在线上的投影
        /// </summary>
        /// <param name="point"></param>
        /// <param name="line"></param>
        /// <returns></returns>
        public static Point3 GetProjP2L(Point3 point, Line3 line)
        {
            float x0 = point.X;
            float y0 = point.Y;
            float z0 = point.Z;
            float x1 = line.Point1.X;
            float y1 = line.Point1.Y;
            float z1 = line.Point1.Z;
            float m = line.Point2.X - line.Point1.X;
            float n = line.Point2.Y - line.Point1.Y;
            float p = line.Point2.Z - line.Point1.Z;
            float t = ((x0 - x1) * m + (y0 - y1) * n + (z0 - z1) * p) / (m * m + n * n + p * p);
            float x2 = x0 + m * t;
            float y2 = x0 + n * t;
            float z2 = x0 + p * t;
            return new Point3(x2, y2, z2);
        }
        /// <summary>
        /// 点在线上的投影
        /// </summary>
        /// <param name="point"></param>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static Point3 GetProjP2L(Point3 point, Vector3 vector)
        {
            float x0 = point.X;
            float y0 = point.Y;
            float z0 = point.Z;
            float x1 = vector.Start.X;
            float y1 = vector.Start.Y;
            float z1 = vector.Start.Z;
            float m = vector.End.X - vector.Start.X;
            float n = vector.End.Y - vector.Start.Y;
            float p = vector.End.Z - vector.Start.Z;
            float t = ((x0 - x1) * m + (y0 - y1) * n + (z0 - z1) * p) / (m * m + n * n + p * p);
            float x2 = x0 + m * t;
            float y2 = x0 + n * t;
            float z2 = x0 + p * t;
            return new Point3(x2, y2, z2);
        }
        /// <summary>
        /// 点在面上的投影
        /// </summary>
        /// <param name="point"></param>
        /// <param name="face"></param>
        /// <returns></returns>
        public static Point3 GetProjP2F(Point3 point, Face face)
        {
            float x0 = point.X;
            float y0 = point.Y;
            float z0 = point.Z;
            float x1 = face.P.X;
            float y1 = face.P.Y;
            float z1 = face.P.Z;
            float m = face.N.End.X;
            float n = face.N.End.Y;
            float p = face.N.End.Z;
            float t = ((x0 - x1) * m + (y0 - y1) * n + (z0 - z1) * p) / (m * m + n * n + p * p);
            float x2 = x0 + m * t;
            float y2 = x0 + n * t;
            float z2 = x0 + p * t;
            return new Point3(x2, y2, z2);
        }
        /// <summary>
        /// 线在线上的投影
        /// </summary>
        /// <param name="line1"></param>
        /// <param name="line2"></param>
        /// <returns></returns>
        public static Line3 GetProjL2L(Line3 line1, Line3 line2)
        {
            Point3 p1 = GetProjP2L(line1.Point1, line2);
            Point3 p2 = GetProjP2L(line1.Point2, line2);
            return new Line3(p1, p2);
        }
        /// <summary>
        /// 线在面上的投影
        /// </summary>
        /// <param name="line"></param>
        /// <param name="face"></param>
        /// <returns></returns>
        public static Line3 GetProjL2F(Line3 line, Face face)
        {
            Point3 p1 = GetProjP2F(line.Point1, face);
            Point3 p2 = GetProjP2F(line.Point2, face);
            return new Line3(p1, p2);
        }
        /// <summary>
        /// 线线夹角
        /// </summary>
        /// <param name="line1"></param>
        /// <param name="line2"></param>
        /// <returns></returns>
        public static double GetJJL2L(Line3 line1, Line3 line2)
        {
            Vector3 v1 = LineToVector(line1);
            Vector3 v2 = LineToVector(line2);
            return GetJJ(v1, v2);
        }
        /// <summary>
        /// 线面夹角
        /// </summary>
        /// <param name="line"></param>
        /// <param name="face"></param>
        /// <returns></returns>
        public static double GetJJL2F(Line3 line, Face face)
        {
            Vector3 v1 = LineToVector(line);
            double jj = GetJJ(v1, face.N);
            if (jj > 90)
                jj = 180 - jj;
            return jj;
        }
        /// <summary>
        /// 面面夹角
        /// </summary>
        /// <param name="face1"></param>
        /// <param name="face2"></param>
        /// <returns></returns>
        public static double GetJJF2F(Face face1, Face face2)
        {
            double jj = GetJJ(face1.N, face2.N);
            if (jj > 90)
                jj = 180 - jj;
            return jj;
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
        /// 求三元一次（线性）方程组的根
        /// </summary>
        /// <param name="coef"></param>
        /// <param name="cons"></param>
        /// <returns></returns>
        public static double[] GetRootsByCramer3D(double[,] coef, double[] cons)
        {
            if (coef.GetLength(0) == 3 && coef.GetLength(0) == coef.GetLength(1) && coef.GetLength(0) == cons.GetLength(0))
            {
                double d = Det(coef);
                if (d != 0)
                {
                    double[,] coef1 = new double[3, 3];
                    for (int i = 0; i < 3; i++)
                    {
                        coef1[i, 0] = cons[i];
                        coef1[i, 1] = coef[i, 1];
                        coef1[i, 2] = coef[i, 2];
                    }
                    double[,] coef2 = new double[3, 3];
                    for (int i = 0; i < 3; i++)
                    {
                        coef2[i, 1] = cons[i];
                        coef2[i, 0] = coef[i, 0];
                        coef2[i, 2] = coef[i, 2];
                    }
                    double[,] coef3 = new double[3, 3];
                    for (int i = 0; i < 3; i++)
                    {
                        coef3[i, 2] = cons[i];
                        coef3[i, 0] = coef[i, 0];
                        coef3[i, 1] = coef[i, 1];
                    }
                    double d1 = Det(coef1);
                    double d2 = Det(coef2);
                    double d3 = Det(coef3);
                    double[] root = new double[3];
                    root[0] = d1 / d;
                    root[1] = d2 / d;
                    root[2] = d3 / d;
                    return root;
                }
            }
            return null;
        }
        /// <summary>
        /// 求多元一次（线性）方程组的根
        /// </summary>
        /// <param name="coef"></param>
        /// <param name="cons"></param>
        /// <returns></returns>
        public static double[] GetRootsByCramer(double[,] coef, double[] cons)
        {
            int rows = coef.GetLength(0);
            int cols = coef.GetLength(1);
            if (rows == cons.Length)
            {
                double d = Det(coef);
                if (d != 0)
                {
                    double[][,] coefs = new double[cols][,];
                    for (int i = 0; i < cols; i++)
                    {
                        coefs[i] = new double[rows, cols];
                        for (int j = 0; j < rows; j++)
                        {
                            for (int k = 0; k < cols; k++)
                            {
                                if (i == k)
                                {
                                    coefs[i][j, k] = cons[j];
                                }
                                else
                                {
                                    coefs[i][j, k] = coef[j, k];
                                }
                            }
                        }
                    }

                    double[] ds = new double[cols];
                    for (int i = 0; i < cols; i++)
                    {
                        //CSharpAlgorithm.Algorithm.Matrix mtx = new CSharpAlgorithm.Algorithm.Matrix(coefs[i]);
                        //System.Windows.Forms.MessageBox.Show(mtx.ToString());
                        ds[i] = Det(coefs[i]);
                    }
                    double[] root = new double[cols];
                    for (int i = 0; i < cols; i++)
                    {
                        root[i] = ds[i] / d;
                    }
                    return root;
                }
            }
            return null;
        }
        /// <summary>
        /// 创建三维转换矩阵，注：transMtx = new double[4, 3]
        /// </summary>
        /// <param name="Rx"></param>
        /// <param name="Ry"></param>
        /// <param name="Rz"></param>
        /// <param name="Sx"></param>
        /// <param name="Sy"></param>
        /// <param name="Sz"></param>
        /// <param name="Tx"></param>
        /// <param name="Ty"></param>
        /// <param name="Tz"></param>
        /// <returns></returns>
        public static double[,] CreateTransMatrix(double Rx, double Ry, double Rz, double Sx, double Sy, double Sz, double Tx, double Ty, double Tz)
        {
            double CosRx, CosRy, CosRz;
            double SinRx, SinRy, SinRz;
            CosRx = Math.Cos(Rx); //Used 6x
            CosRy = Math.Cos(Ry); //Used 4x
            CosRz = Math.Cos(Rz); //Used 4x
            SinRx = Math.Sin(Rx); //Used 5x
            SinRy = Math.Sin(Ry); //Used 5x
            SinRz = Math.Sin(Rz); //Used 5x
            //total of 29 trig functions
            //23 trig functions cancelled out by
            //this optimisation; hence the 2.6x speed increase.

            double[,] transMtx = new double[4, 3];//4];
            transMtx[0, 0] = Sx * CosRy * CosRz;
            transMtx[0, 1] = Sx * CosRy * SinRz;
            transMtx[0, 2] = -(Sx * SinRy);

            transMtx[1, 0] = -(Sy * CosRx * SinRz) + Sy * SinRx * SinRy * CosRz;
            transMtx[1, 1] = Sy * CosRx * CosRz + Sy * SinRx * SinRy * SinRz;
            transMtx[1, 2] = Sy * SinRx * CosRy;

            transMtx[2, 0] = Sz * SinRx * SinRz + Sz * CosRx * SinRy * CosRz;
            transMtx[2, 1] = -(Sz * SinRx * CosRx) + Sz * CosRx * SinRy * SinRz;
            transMtx[2, 2] = Sz * CosRx * CosRy;

            transMtx[3, 0] = Tx;
            transMtx[3, 1] = Ty;
            transMtx[3, 2] = Tz;
            //transMtx[3, 3] = 1;
            return transMtx;
        }
        /// <summary>
        /// 综合转换
        /// </summary>
        /// <param name="srcPt"></param>
        /// <param name="Rx">以弧度计量的角度(x分量)</param>
        /// <param name="Ry">以弧度计量的角度(y分量)</param>
        /// <param name="Rz">以弧度计量的角度(z分量)</param>
        /// <param name="Sx"></param>
        /// <param name="Sy"></param>
        /// <param name="Sz"></param>
        /// <param name="Tx"></param>
        /// <param name="Ty"></param>
        /// <param name="Tz"></param>
        /// <returns></returns>
        public static Point3 GetTransformation(Point3 srcPt, double Rx, double Ry, double Rz, double Sx, double Sy, double Sz, double Tx, double Ty, double Tz)
        {
            double[,] trans = CreateTransMatrix(Rx, Ry, Rz, Sx, Sy, Sz, Tx, Ty, Tz);
            Point3 retPt = new Point3();
            double[,] A = new double[1, 4];
            A[0, 0] = srcPt.X;
            A[0, 1] = srcPt.Y;
            A[0, 2] = srcPt.Z;
            A[0, 3] = 1;
            double[,] B = new double[1, 3];
            B = Mul.MatrixMul(A, trans);
            retPt.X = (float)B[0, 0];
            retPt.Y = (float)B[0, 1];
            retPt.Z = (float)B[0, 2];
            return retPt;
        }
        /// <summary>
        /// 基于指定基矩阵的向量
        /// </summary>
        /// <param name="vector">原向量</param>
        /// <param name="baseMatrix">基矩阵</param>
        /// <returns></returns>
        public static Vector3 BaseToSystem(Vector3 vector, double[,] baseMatrix)
        {
            Vector3 v = null;
            if (baseMatrix.GetLength(0) == 3 && baseMatrix.GetLength(1) == 3)
            {
                v = Mul.Vector3MulMatrix(vector, baseMatrix);
            }
            return v;
        }
        /// <summary>
        /// 求逆阵
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static double[,] GetInverseMatrix(double[,] a)
        {
            double[,] b = new double[3, 3];
            double[] t = new double[4];
            //求伴随阵
            for (int ci = 0; ci < 3; ci++)
            {
                for (int cj = 0; cj < 3; cj++)
                {
                    //求余子式
                    int n = 0;
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            if (i != ci && j != cj)
                            {
                                t[n] = a[i, j];
                                n++;
                            }
                        }
                    }
                    b[cj, ci] = IsNeg(ci + cj) * (t[0] * t[3] - t[1] * t[2]);
                }
            }
            //求矩阵的行列式值
            double rValue = 0;
            for (int i = 0; i < 3; i++)
            {
                rValue = rValue + a[i, 0] * b[0, i];
            }
            if (rValue != 0) //行列式值不为0，才存在逆阵 
            {
                //数乘
                rValue = 1 / rValue;
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        b[i, j] = b[i, j] * rValue;
                    }
                }
            }
            return b;
        }
        /// <summary>
        /// 求垂直平分线，第一个点Point1即为垂足，且长度为1
        /// </summary>
        /// <param name="line"></param>
        /// <returns></returns>
        public static Line3 GetCZAVGLine(Line3 line)
        {
            Line3 l = null;
            float x1 = line.Point1.X;
            float x2 = line.Point2.X;
            float y1 = line.Point1.Y;
            float y2 = line.Point2.Y;
            float z1 = line.Point1.Z;
            float z2 = line.Point2.Z;
            float midX1 = (x1 + x2) / 2;
            float midY1 = (y1 + y2) / 2;
            float midZ1 = (z1 + z2) / 2;
            Point3 p1 = new Point3(midX1, midY1, midZ1);
            //bool inLine;
            ////p1为垂足，不能作为线外的点
            //l = GetCZLine( , line, out inLine);


            /*
            float xl = 0;
            if (x2 - x1 != 0)
                xl = (y2 - y1) / (x2 - x1);
            else
                xl = float.MaxValue;
            float xl1 = -1 / xl;            
            float xll = 0;
            if (x2 - x1 != 0)
                xll = (z2 - z1) / (x2 - x1);
            else
                xll = float.MaxValue;
            float xll1 = -1 / xll;
            float x = (float)Math.Sqrt(1 / (1 + xl1 * xl1 + xll1 * xll1));
            float y = xl1 * x;
            float z = xll1 * x;
            Vector3 v = new Vector3(new Point3(x, y, z));
            v = Vector3.Transform(v, p1);
            Point3 p2 = GetPointOfLengthFormStartP(v, 1);            
            l = new Line3(p1, p2);
            */
            return l;
        }
        /// <summary>
        /// 求过点p垂直于line的垂直线，第一个点Point1即为垂足
        /// </summary>
        /// <param name="p"></param>
        /// <param name="line"></param>
        /// <param name="ifInLine"></param>
        /// <returns></returns>
        public static Line3 GetCZLine(Point3 p, Line3 line, out bool ifInLine)
        {
            Line3 l = null;
            /*
            float x2, y2, z2;//垂足
            float x0 = p.X;
            float y0 = p.Y;
            float z0 = p.Z;
            float x1 = line.Point1.X;
            float y1 = line.Point1.Y;
            float z1 = line.Point1.Z;
            float m = line.Point2.X - line.Point1.X;
            float n = line.Point2.Y - line.Point1.Y;
            float k = line.Point2.Z - line.Point1.Z;
            float t = ((x0 - x1) * m + (y0 - y1) * n + (z0 - z1) * k) / (m * m + n * n + k * k);
            x2 = x0 + m * t;
            y2 = x0 + n * t;
            z2 = x0 + k * t;
            Point3 crossP = new Point3(x2, y2, z2);
            ifInLine = InLine(line.Point1, line.Point2, crossP);
            l = new Line3(crossP, p);
            */
            //求面：
            Vector3 n = new Vector3(line.Point1, line.Point2, false);
            Face face = new Face(n, p);
            //求线面交点
            Point3 cross = LFcross.GetCross(line, face);
            ifInLine = InLine(line.Point1, line.Point2, cross);
            l = new Line3(cross, p);
            return l;
        }
        /// <summary>
        /// X是否在AB所在的直线上
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="X"></param>
        /// <returns></returns>
        private static bool InLine(Point3 A, Point3 B, Point3 X)
        {
            if ((A.X - X.X) * (X.X - B.X) >= 0 && (A.Y - X.Y) * (X.Y - B.Y) >= 0 && (A.Z - X.Z) * (X.Z - B.Z) >= 0)
                return true;
            return false;
        }
        /// <summary>
        /// 返回与指定向量vector垂直且经过vector上一定点的平面
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="point"></param>
        /// <returns></returns>
        public static Face GetPerpendicular(Vector3 vector, Point3 point)
        {
            vector = Vector3.GetUnitVector(vector);
            Face face = new Face(vector, point);
            return face;
        }
    }
    /// <summary>
    /// 乘法类
    /// </summary>
    public class Mul
    {
        /// <summary>
        /// 矩阵的乘法中乘积和
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        private static double Sum(double[,] A, double[,] B, int x, int y)
        {
            double sum = 0;
            int R = B.GetLength(0);
            for (int i = 0; i < R; i++)
            {
                sum += A[x, i] * B[i, y];
            }
            return sum;
        }
        /// <summary>
        /// 矩阵的乘法
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static double[,] MatrixMul(double[,] A, double[,] B)
        {
            double[,] C = null;
            if (A.GetLength(1) == B.GetLength(0))
            {
                int R1 = A.GetLength(0);
                int R2 = B.GetLength(1);
                C = new double[R1, R2];
                for (int x = 0; x < R1; x++)
                {
                    for (int y = 0; y < R2; y++)
                    {
                        C[x, y] = Sum(A, B, x, y);
                    }
                }
            }
            else
            {
                //System.Windows.Forms.MessageBox.Show("Can not Process the MatrixMul, because A.ColumnCount <> B.RowCount!", "MessageBox");
            }
            return C;
        }
        /// <summary>
        /// 向量转化为矩阵
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static CSharpAlgorithm.Algorithm.Matrix VectorToMatrix(Vector3 vector)
        {
            CSharpAlgorithm.Algorithm.Matrix matrix = new CSharpAlgorithm.Algorithm.Matrix(1, 3);
            vector = Vector3.Formatting(vector);
            matrix.SetElement(0, 0, vector.End.X);
            matrix.SetElement(0, 0, vector.End.Y);
            matrix.SetElement(0, 0, vector.End.Z);
            return matrix;
        }
        /// <summary>
        /// 矩阵转化为向量
        /// </summary>
        /// <param name="matrix"></param>
        /// <returns></returns>
        public static Vector3 MatrixToVector(CSharpAlgorithm.Algorithm.Matrix matrix)
        {
            Vector3 vector = null;
            if (matrix.Columns == 3 && matrix.Rows == 1)
            {
                vector = new Vector3(new Point3((float)matrix.GetElement(0, 0), (float)matrix.GetElement(0, 1), (float)matrix.GetElement(0, 2)));
            }
            return vector;
        }
        /*
        /// <summary>
        /// 获取同时垂直于X1,X2的向量
        /// </summary>
        /// <param name="X1"></param>
        /// <param name="X2"></param>
        /// <returns></returns>
        public static CSharpAlgorithm.Algorithm.Matrix GetXmul(CSharpAlgorithm.Algorithm.Matrix X1, CSharpAlgorithm.Algorithm.Matrix X2)
        {
            CSharpAlgorithm.Algorithm.Matrix X = new CSharpAlgorithm.Algorithm.Matrix(1, 3);
            CSharpAlgorithm.Algorithm.Matrix X3 = new CSharpAlgorithm.Algorithm.Matrix(3, 3);
            X3.SetElement(0, 0, 0);
            X3.SetElement(1, 1, 0);
            X3.SetElement(2, 2, 0);
            X3.SetElement(0, 1, -X2.GetElement(0, 2));
            X3.SetElement(0, 2, X2.GetElement(0, 1));
            X3.SetElement(1, 0, X2.GetElement(0, 2));
            X3.SetElement(1, 2, -X2.GetElement(0, 0));
            X3.SetElement(2, 0, -X2.GetElement(0, 1));
            X3.SetElement(0, 1, X2.GetElement(0, 0));
            X.SetValue(X1.Multiply(X3));
            return X;
        }*/
        /// <summary>
        /// 叉乘，两向量平行的充要条件是叉乘为 (0, 0, 0)
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static Vector3 GetXmul(Vector3 v1, Vector3 v2)
        {
            v1 = Vector3.Formatting(v1);
            v2 = Vector3.Formatting(v2);
            Point3 p = new Point3(v1.End.Y * v2.End.Z - v1.End.Z * v2.End.Y, v1.End.Z * v2.End.X - v1.End.X * v2.End.Z, v1.End.X * v2.End.Y - v1.End.Y * v2.End.X);
            Vector3 v = new Vector3(p);
            return v;
        }
        /// <summary>
        /// 列向量转化为行向量
        /// </summary>
        /// <param name="X2"></param>
        /// <returns></returns>
        public static CSharpAlgorithm.Algorithm.Matrix ColumeToRow(CSharpAlgorithm.Algorithm.Matrix X2)
        {
            CSharpAlgorithm.Algorithm.Matrix X3 = null;
            if (X2.Columns == 1 && X2.Rows == 3)
            {
                X3 = new CSharpAlgorithm.Algorithm.Matrix(3, 1);
                X3.SetElement(0, 0, X2.GetElement(0, 0));
                X3.SetElement(1, 0, X2.GetElement(0, 1));
                X3.SetElement(2, 0, X2.GetElement(0, 2));
            }
            return X3;
        }
        /// <summary>
        /// 行向量转化为列向量
        /// </summary>
        /// <param name="X2"></param>
        /// <returns></returns>
        public static CSharpAlgorithm.Algorithm.Matrix RowToColume(CSharpAlgorithm.Algorithm.Matrix X2)
        {
            CSharpAlgorithm.Algorithm.Matrix X3 = null;
            if (X2.Columns == 3 && X2.Rows == 1)
            {
                X3 = new CSharpAlgorithm.Algorithm.Matrix(3, 1);
                X3.SetElement(0, 0, X2.GetElement(0, 0));
                X3.SetElement(1, 0, X2.GetElement(0, 1));
                X3.SetElement(2, 0, X2.GetElement(0, 2));
            }
            return X3;
        }
        /*
        /// <summary>
        /// 两向量垂直的充要条件是点积为 0
        /// </summary>
        /// <param name="X1"></param>
        /// <param name="X2"></param>
        /// <returns></returns>
        public static double GetPmul(CSharpAlgorithm.Algorithm.Matrix X1, CSharpAlgorithm.Algorithm.Matrix X2)
        {
            double p = 0;
            CSharpAlgorithm.Algorithm.Matrix X = new CSharpAlgorithm.Algorithm.Matrix();
            CSharpAlgorithm.Algorithm.Matrix X3 = RowToColume(X2);
            X.SetValue(X1.Multiply(X3));
            p = X[0, 0];
            return p;
        }*/
        /// <summary>
        /// 点乘，两向量垂直的充要条件是点乘为 0
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static float GetPmul(Vector3 v1, Vector3 v2)
        {
            v1 = Vector3.Formatting(v1);
            v2 = Vector3.Formatting(v2);
            return v1.End.X * v2.End.X + v1.End.Y * v2.End.Y + v1.End.Z * v2.End.Z;
        }
        /*
        /// <summary>
        /// 三向量共面的充要条件是混合积为 0
        /// </summary>
        /// <param name="X1"></param>
        /// <param name="X2"></param>
        /// <param name="X3"></param>
        /// <returns></returns>
        public static double MixMul(CSharpAlgorithm.Algorithm.Matrix X1, CSharpAlgorithm.Algorithm.Matrix X2, CSharpAlgorithm.Algorithm.Matrix X3)
        {
            return GetPmul(GetXmul(X1, X2), X3);
        }*/
        /// <summary>
        /// 混合积，三向量共面的充要条件是混合积为 0
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <param name="v3"></param>
        /// <returns></returns>
        public static double MixMul(Vector3 v1, Vector3 v2, Vector3 v3)
        {
            return GetPmul(GetXmul(v1, v2), v3);
        }
/// <summary>
/// 二维向量与矩阵相乘
/// </summary>
/// <param name="vector"></param>
/// <param name="baseMatrix"></param>
/// <returns></returns>
        public static Vector2 Vector2MulMatrix(Vector2 vector, double[,] baseMatrix)
        {
            vector = Vector2.Formatting(vector);
            double[,] point = new double[1, 2];
            point[0, 0] = vector.End.X;
            point[0, 1] = vector.End.Y;
            double[,] vs = MatrixMul(point, baseMatrix);
            Vector2 v = new Vector2(new Point2((float)vs[0, 0], (float)vs[0, 1]));
            return v;
        }
/// <summary>
        /// 三维向量与矩阵相乘
/// </summary>
/// <param name="vector"></param>
/// <param name="baseMatrix"></param>
/// <returns></returns>
        public static Vector3 Vector3MulMatrix(Vector3 vector, double[,] baseMatrix)
        {
            vector = Vector3.Formatting(vector);
            double[,] point = new double[1, 3];
            point[0, 0] = vector.End.X;
            point[0, 1] = vector.End.Y;
            point[0, 2] = vector.End.Z;
            double[,] vs = MatrixMul(point, baseMatrix);
            Vector3 v = new Vector3(new Point3((float)vs[0, 0], (float)vs[0, 1], (float)vs[0, 2]));
            return v;
        }
    }
    /// <summary>
    /// 二维点
    /// </summary>
    public class Point2
    {
        float xx = 0;
        float yy = 0;
        /// <summary>
        /// 构造函数
        /// </summary>
        public Point2()
        {

        }
        /// <summary>
        /// 初始化点
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        public Point2(float x, float y)
        {
            xx = x;
            yy = y;
        }
        /// <summary>
        /// 点的X坐标
        /// </summary>
        public float X
        {
            get
            {
                return xx;
            }
            set
            {
                xx = value;
            }
        }
        /// <summary>
        /// 点的Y坐标
        /// </summary>
        public float Y
        {
            get
            {
                return yy;
            }
            set
            {
                yy = value;
            }
        }
        /// <summary>
        /// 零点
        /// </summary>
        public static Point2 ZERO
        {
            get
            {
                return new Point2(0, 0);
            }
        }
        /// <summary>
        /// 是否相等
        /// </summary>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <returns></returns>
        public static bool IsEqual(Point2 point1, Point2 point2)
        {
            return (Math.Abs(point1.xx - point2.xx) < Common.Eps) && (Math.Abs(point1.yy - point2.yy) < Common.Eps);
        }
        /// <summary>
        /// 将二维坐标点转化为可视化字符串
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return xx.ToString() + "," + yy.ToString();
        }
        /// <summary>
        /// 字符串转化为点
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        public static Point2 Parse(string s)
        {
            Point2 p = null;
            string[] tmp = s.Split(',');
            if (tmp.Length == 2)
            {
                string s1 = tmp[0];
                string s2 = tmp[1];
                if (TestNumeric.IsReal(s1) && TestNumeric.IsReal(s2))
                    p = new Point2(float.Parse(s1), float.Parse(s2));
            }
            if (p == null)
                throw new Exception("Wrong Point2 Formation!");
            else
                return p;
        }
    }
/// <summary>
/// 三维点
/// </summary>
    public class Point3
    {
        float xx = 0;
        float yy = 0;
        float zz = 0;
        /// <summary>
        /// 构造函数
        /// </summary>
        public Point3()
        {

        }
        /// <summary>
        /// 初始化点
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        public Point3(float x, float y, float z)
        {
            xx = x;
            yy = y;
            zz = z;
        }
        /// <summary>
        /// 点的X坐标
        /// </summary>
        public float X
        {
            get
            {
                return xx;
            }
            set
            {
                xx = value;
            }
        }
        /// <summary>
        /// 点的Y坐标
        /// </summary>
        public float Y
        {
            get
            {
                return yy;
            }
            set
            {
                yy = value;
            }
        }
        /// <summary>
        /// 点的Z坐标
        /// </summary>
        public float Z
        {
            get
            {
                return zz;
            }
            set
            {
                zz = value;
            }
        }
        /// <summary>
        /// 零点
        /// </summary>
        public static Point3 ZERO
        {
            get
            {
                return new Point3(0, 0, 0);
            }
        }
        /// <summary>
        /// 是否相等
        /// </summary>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <returns></returns>
        public static bool IsEqual(Point3 point1, Point3 point2)
        {
            return (Math.Abs(point1.xx - point2.xx) < Common.Eps) && (Math.Abs(point1.yy - point2.yy) < Common.Eps) && (Math.Abs(point1.zz - point2.zz) < Common.Eps);
        }
        /// <summary>
        /// 将三维坐标点转化为可视化字符串
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return xx.ToString() + "," + yy.ToString() + "," + zz.ToString();
        }
        /// <summary>
        /// 字符串转化为点
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        public static Point3 Parse(string s)
        {
            Point3 p = null;
            string[] tmp = s.Split(',');
            if (tmp.Length == 3)
            {
                string s1 = tmp[0];
                string s2 = tmp[1];
                string s3 = tmp[2];
                if (TestNumeric.IsReal(s1) && TestNumeric.IsReal(s2) && TestNumeric.IsReal(s3))
                    p = new Point3(float.Parse(s1), float.Parse(s2), float.Parse(s3));
            }
            if (p == null)
                throw new Exception("Wrong Point3 Formation!");
            else
                return p;
        }
    }
    /// <summary>
    /// 面
    /// </summary>
    public class Face
    {
        Vector3 n = new Vector3();
        Point3 p = new Point3();
        double aa, bb, cc, dd;
        /// <summary>
        /// 构造函数
        /// </summary>
        public Face()
        {

        }
        /// <summary>
        /// 法向量与点构造面
        /// </summary>
        /// <param name="vector"></param>
        /// <param name="point"></param>
        public Face(Vector3 vector, Point3 point)
        {
            n = Vector3.GetUnitVector(vector);
            p = point;
            aa = n.End.X;
            bb = n.End.Y;
            cc = n.End.Z;
            dd = -(aa * point.X + bb * point.Y + cc * point.Z);
        }
        /// <summary>
        /// 解析式构造面
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <param name="d"></param>
        public Face(double a, double b, double c, double d)
        {
            aa = a;
            bb = b;
            cc = c;
            dd = d;
            n = Vector3.GetUnitVector(new Vector3(new Point3((float)a, (float)b, (float)c)));
            //p = new Point3();
        }
        /// <summary>
        /// 三点定面
        /// </summary>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <param name="point3"></param>
        public Face(Point3 point1, Point3 point2, Point3 point3)
        {
            n = Mul.GetXmul(new Vector3(point1, point2, false), new Vector3(point1, point3, false));
            n = Vector3.GetUnitVector(n);
            p = point1;
            aa = n.End.X;
            bb = n.End.Y;
            cc = n.End.Z;
            dd = -(aa * point1.X + bb * point1.Y + cc * point1.Z);
        }
        /// <summary>
        /// 法向量
        /// </summary>
        public Vector3 N
        {
            get
            {
                return n;
            }
            set
            {
                n = Vector3.GetUnitVector(value);
            }
        }
        /// <summary>
        /// 面上的某一点
        /// </summary>
        public Point3 P
        {
            get
            {
                return p;
            }
            set
            {
                p = value;
            }
        }
        /// <summary>
        /// 解析式中的A
        /// </summary>
        public double A
        {
            get
            {
                return aa;
            }
            set
            {
                aa = value;
            }
        }
        /// <summary>
        /// 解析式中的B
        /// </summary>
        public double B
        {
            get
            {
                return bb;
            }
            set
            {
                bb = value;
            }
        }
        /// <summary>
        /// 解析式中的C
        /// </summary>
        public double C
        {
            get
            {
                return cc;
            }
            set
            {
                cc = value;
            }
        }
        /// <summary>
        /// 解析式中的D
        /// </summary>
        public double D
        {
            get
            {
                return dd;
            }
            set
            {
                dd = value;
            }
        }
        /// <summary>
        /// 转化为可视化的字符串
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return "P = " + p.ToString() + ", N = " + n.End.ToString();
        }
        /// <summary>
        /// 旋转
        /// </summary>
        /// <param name="face"></param>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static Face Transform(Face face, Vector3 vector)
        {
            face.N = vector;
            return face;
        }
        /// <summary>
        /// 平移
        /// </summary>
        /// <param name="face"></param>
        /// <param name="point"></param>
        /// <returns></returns>
        public static Face Transform(Face face, Point3 point)
        {
            face.P = point;
            return face;
        }
        /// <summary>
        /// 先旋转，后平移
        /// </summary>
        /// <param name="face"></param>
        /// <param name="vector"></param>
        /// <param name="point"></param>
        /// <returns></returns>
        public static Face Transform(Face face, Vector3 vector, Point3 point)
        {
            face.N = vector;
            face.P = point;
            return face;
        }
        /// <summary>
        /// 先平移，后旋转
        /// </summary>
        /// <param name="face"></param>
        /// <param name="point"></param>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static Face Transform(Face face, Point3 point, Vector3 vector)
        {
            face.P = point;
            face.N = vector;
            return face;
        }
        /// <summary>
        /// 是否相等
        /// </summary>
        /// <param name="face1"></param>
        /// <param name="face2"></param>
        /// <returns></returns>
        public static bool IsEqual(Face face1, Face face2)
        {
            Vector3 v1 = face1.N;
            Vector3 v2 = face2.N;
            if (Vector3.IsEqual(v1, v2) && Point3.IsEqual(face1.P, face2.P))
                return true;
            else
                return false;
        }
        /// <summary>
        /// 空间中的两线段是否共面
        /// </summary>
        /// <param name="l1"></param>
        /// <param name="l2"></param>
        /// <returns></returns>
        public static bool ShareFace(Line3 l1, Line3 l2)
        {
            Vector3 v1 = new Vector3(l1.Point1, l1.Point2, false);
            Vector3 v2 = new Vector3(l2.Point1, l2.Point2, false);
            //System.Windows.Forms.MessageBox.Show(Mul.GetXmul(v1, v2).ToString());
            Vector3 ZERO = Vector3.Empty;
            return Vector3.IsEqual(Mul.GetXmul(v1, v2), ZERO) || Vector3.IsVectorCross(v1, v2);
        }
        /// <summary>
        /// 空间中的两向量是否共面
        /// </summary>
        /// <param name="v1"></param>
        /// <param name="v2"></param>
        /// <returns></returns>
        public static bool ShareFace(Vector3 v1, Vector3 v2)
        {
            //System.Windows.Forms.MessageBox.Show(Mul.GetXmul(v1, v2).ToString());
            Vector3 ZERO = Vector3.Empty;
            return Vector3.IsEqual(Mul.GetXmul(v1, v2), ZERO) || Vector3.IsVectorCross(v1, v2);
        }
    }
    /// <summary>
    /// 二维线
    /// </summary>
    public class Line2
    {
        Point2 point1 = new Point2();
        Point2 point2 = new Point2();
        float a, b, c;
        float k;
        /// <summary>
        /// 构造函数
        /// </summary>
        public Line2()
        {

        }
        /// <summary>
        /// 两点定线
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        public Line2(Point2 p1, Point2 p2)
        {
            point1 = p1;
            point2 = p2;
            if (Math.Abs(p2.X - p1.X) > Common.Eps)
            {
                k = (p2.Y - p1.Y) / (p2.X - p1.X);
            }
            else
            {
                k = float.MaxValue;
            }
            a = k;
            b = -1;
            c = p1.Y - k * p1.X;
        }
        /// <summary>
        /// 返回长度为length，斜率为xl，过p且以p为始点的线段
        /// </summary>
        /// <param name="p"></param>
        /// <param name="xl"></param>
        /// <param name="length"></param>
        public Line2(Point2 p, float xl, double length)
        {
            point1 = p;
            k = xl;
            float x = p.X + (float)(length * Math.Sqrt(1 / (1 + k * k)));
            float y = p.Y + (float)(length * Math.Sqrt(k / (1 + k * k)));
            point2 = new Point2(x, y);
            a = xl;
            b = -1;
            c = p.Y - xl * p.X;
        }
        /// <summary>
        /// 解析式构造线
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="C"></param>
        /// <param name="length"></param>
        public Line2(float A, float B, float C, double length)
        {
            if (A < 0)
            {
                a = -A;
                b = -B;
                c = -C;
            }
            else
            {
                a = A;
                b = B;
                c = C;
            }
            k = -a / b;
            point1 = new Point2(-c / b, 0);
            float x = point1.X + (float)(length * Math.Sqrt(1 / (1 + k * k)));
            float y = point1.Y + (float)(length * Math.Sqrt(k / (1 + k * k)));
            point2 = new Point2(x, y);
        }
        /// <summary>
        /// 根据已知两点坐标，求过这两点的直线解析方程： a*x+b*y+c = 0  (a >= 0) 
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static Line2 MakeLine(Point2 p1, Point2 p2)
        {
            /*Line2 tl = new Line2();
            int sign = 1;
            tl.a = p2.Y - p1.Y;
            if (tl.a < 0)
            {
                sign = -1;
                tl.a = sign * tl.a;
            }
            tl.b = sign * (p1.X - p2.X);
            tl.c = sign * (p1.Y * p2.X - p1.X * p2.Y);
            tl.k = -tl.a / tl.b;
            tl.point1 = p1;
            tl.point2 = p2;*/
            return new Line2(p1, p2);
        }
        /// <summary>
        /// 始点
        /// </summary>
        public Point2 Point1
        {
            get
            {
                return point1;
            }
            set
            {
                point1 = value;
            }
        }
        /// <summary>
        /// 终点
        /// </summary>
        public Point2 Point2
        {
            get
            {
                return point2;
            }
            set
            {
                point2 = value;
            }
        }
        /// <summary>
        /// 斜率
        /// </summary>
        public float K
        {
            get
            {
                return k;
            }
        }
        /// <summary>
        /// 解析式中的A
        /// </summary>
        public float A
        {
            get
            {
                return a;
            }
        }
        /// <summary>
        /// 解析式中的B
        /// </summary>
        public float B
        {
            get
            {
                return b;
            }
        }
        /// <summary>
        /// 解析式中的C
        /// </summary>
        public float C
        {
            get
            {
                return c;
            }
        }
        /// <summary>
        /// 转化为可视化的字符串
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return point1.X.ToString() + "," + point1.Y.ToString() + " ----> " + point2.X.ToString() + "," + point2.Y.ToString();
        }

    }
    /// <summary>
    /// 三维线
    /// </summary>
    public class Line3
    {
        Point3 point1 = new Point3();
        Point3 point2 = new Point3();
        float a, b, c;
        float kyx = 0;
        float kzx = 0;
        /// <summary>
        /// 构造函数
        /// </summary>
        public Line3()
        {

        }
        /// <summary>
        /// 两点定线
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        public Line3(Point3 p1, Point3 p2)
        {
            point1 = p1;
            point2 = p2;
            if (Math.Abs(p2.X - p1.X) > Common.Eps)
            {
                kyx = (p2.Y - p1.Y) / (p2.X - p1.X);
            }
            else
            {
                kyx = float.MaxValue;
            }
            a = kyx;
            b = -1;
            c = p1.Y - kyx * p1.X;
        }
        /// <summary>
        /// 返回长度为length，斜率为xl，过p且以p为始点的线段
        /// </summary>
        /// <param name="p"></param>
        /// <param name="xl"></param>
        /// <param name="length"></param>
        public Line3(Point3 p, float xl, double length)
        {
            point1 = p;
            kyx = xl;
            float x = p.X + (float)(length * Math.Sqrt(1 / (1 + xl * xl)));
            float y = p.Y + (float)(length * Math.Sqrt(xl / (1 + xl * xl))); 
            float z = 0;
            point2 = new Point3(x, y, z);
            a = xl;
            b = -1;
            c = p.Y - xl * p.X;
        }
        /// <summary>
        /// 解析式构造线
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="C"></param>
        /// <param name="length"></param>
        public Line3(float A, float B, float C, double length)
        {
            if (A < 0)
            {
                a = -A;
                b = -B;
                c = -C;
            }
            else
            {
                a = A;
                b = B;
                c = C;
            }
            kyx = -a / b;
            float z = 0;
            point1 = new Point3(-c / b, 0, z);
            float x = point1.X + (float)(length * Math.Sqrt(1 / (1 + kyx * kyx)));
            float y = point1.Y + (float)(length * Math.Sqrt(kyx / (1 + kyx * kyx)));
            point2 = new Point3(x, y, z);
        }
        /// <summary>
        /// 根据已知两点坐标，求过这两点的直线解析方程： a*x+b*y+c = 0  (a >= 0) 
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <returns></returns>
        public static Line3 MakeLine(Point3 p1, Point3 p2)
        {
            /*Line3 tl = new Line3();
            int sign = 1;
            tl.a = p2.Y - p1.Y;
            if (tl.a < 0)
            {
                sign = -1;
                tl.a = sign * tl.a;
            }
            tl.b = sign * (p1.X - p2.X);
            tl.c = sign * (p1.Y * p2.X - p1.X * p2.Y);
            tl.kyx = -tl.a / tl.b;
            tl.point1 = p1;
            tl.point2 = p2;*/
            return new Line3(p1, p2);
        }
        /// <summary>
        /// 始点
        /// </summary>
        public Point3 Point1
        {
            get
            {
                return point1;
            }
            set
            {
                point1 = value;
            }
        }
        /// <summary>
        /// 终点
        /// </summary>
        public Point3 Point2
        {
            get
            {
                return point2;
            }
            set
            {
                point2 = value;
            }
        }
        /// <summary>
        /// Y对X的斜率
        /// </summary>
        public float Kyx
        {
            get
            {
                return kyx;
            }
        }
        /// <summary>
        /// Z对X的斜率
        /// </summary>
        public float Kzx
        {
            get
            {
                return kzx;
            }
        }
        /// <summary>
        /// 解析式中的A
        /// </summary>
        public float A
        {
            get
            {
                return a;
            }
        }
        /// <summary>
        /// 解析式中的B
        /// </summary>
        public float B
        {
            get
            {
                return b;
            }
        }
        /// <summary>
        /// 解析式中的C
        /// </summary>
        public float C
        {
            get
            {
                return c;
            }
        }
        /// <summary>
        /// 转化为可视化的字符串
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return point1.X.ToString() + "," + point1.Y.ToString() + "," + point1.Z.ToString() + " ----> " + point2.X.ToString() + "," + point2.Y.ToString() + "," + point2.Z.ToString();
        }
    }
    /// <summary>
    /// 面面距离
    /// </summary>
    public class FFdistance
    {
        /// <summary>
        /// 求面面距离
        /// </summary>
        /// <param name="face1"></param>
        /// <param name="face2"></param>
        /// <returns></returns>
        public static double GetDistance(Face face1, Face face2)
        {
            double dis = 0;
            Vector3 f1N = Vector3.GetUnitVector(face1.N);
            Vector3 f2N = Vector3.GetUnitVector(face2.N);
            if (Vector3.IsEqual(f1N, f2N))
            {
                double d1 = face1.D;
                double d2 = face2.D;
                dis = Math.Abs(d1 - d2) / (Math.Sqrt(face1.A * face1.A + face1.B * face1.B + face1.C * face1.C));
            }
            return dis;
        }
    }
    /// <summary>
    /// 线面交点
    /// </summary>
    public class LFcross
    {
        /// <summary>
        /// 求线面交点
        /// </summary>
        /// <param name="line"></param>
        /// <param name="face"></param>
        /// <returns></returns>
        public static Point3 GetCross(Line3 line, Face face)
        {
            float t = 0;
            Vector3 v = Vector3.LineToVector(line);
            Point3 point = new Point3();
            float deno = Mul.GetPmul(v, face.N);
            if (Math.Abs(deno) < Common.Eps)
            {
                if (Math.Abs(line.Point1.X * face.A + line.Point1.Y * face.B + line.Point1.Z * face.C + face.D) < Common.Eps)
                {
                    t = 0;
                }
            }
            t = -(float)(face.A * line.Point1.X + face.B * line.Point1.Y + face.C * line.Point1.Z + face.D);
            t = t / deno;
            point = new Point3(line.Point1.X + v.End.X * t, line.Point1.Y + v.End.Y * t, line.Point1.Z + v.End.Z * t);
            return point;
        }
    }
    /// <summary>
    /// 线线交点
    /// </summary>
    public class LLcross
    {
        /// <summary>
        /// 求二维线线交点
        /// </summary>
        /// <param name="l1"></param>
        /// <param name="l2"></param>
        /// <param name="cross"></param>
        /// <returns></returns>
        public static bool GetLLcrs2(Line2 l1, Line2 l2, out Point2 cross)
        {
            Vector2 v1 = new Vector2(l1.Point1, l1.Point2);
            Vector2 v2 = new Vector2(l2.Point1, l2.Point2);
            return Vector2.GetCross(v1, v2, out cross);
        }
/// <summary>
        /// 求三维线线交点
/// </summary>
/// <param name="l1"></param>
/// <param name="l2"></param>
/// <returns></returns>
        public static Point3 GetLLcrs3(Line3 l1, Line3 l2)
        {
            Point3 cross = null;
            Vector3 v1 = new Vector3(l1.Point1, l1.Point2, false);
            Vector3 v2 = new Vector3(l2.Point1, l2.Point2, false);
            if (Vector3.IsVectorCross(v1, v2))//两线必须相交
            {
                //cross = Vector3.GetCross(v1, v2);
            }
            return cross;
        }
    }
    /// <summary>
    /// 线线距离
    /// </summary>
    public class LLdistance
    {
        /// <summary>
        ///  求二维线线距
        /// </summary>
        /// <param name="line1"></param>
        /// <param name="line2"></param>
        /// <returns></returns>
        public static double GetDistance2(Line2 line1, Line2 line2)
        {//先求垂直平分线，再求两线的交点，最后是点点距离。
            double dis;
            Vector2 v1=Vector2.LineToVector(line1);
            Vector2 v2=Vector2.LineToVector(line2);
            Line2 l=null;
            //if (Vector2.Parallel(v1,v2))
            if (Vector2.GetLength(v1) < Vector2.GetLength(v2))
            {
                l = Vector2.GetCZAVGLine(line1);
                Vector2 v = Vector2.LineToVector(l);
                Point2 cross;
                Vector2.GetCross(v, v2, out cross);
                dis = PPdistance.GetPPdis2(l.Point1.X, l.Point1.Y,cross.X,cross.Y);
            }
            else
            {
                l = Vector2.GetCZAVGLine(line2);
                Vector2 v = Vector2.LineToVector(l);
                Point2 cross;
                Vector2.GetCross(v, v1, out cross);
                dis = PPdistance.GetPPdis2(l.Point1.X, l.Point1.Y, cross.X, cross.Y);
            }
            return dis;
        }
        /// <summary>
        ///  求三维线线距
        /// </summary>
        /// <param name="line1"></param>
        /// <param name="line2"></param>
        /// <returns></returns>
        public static double GetDistance3(Line3 line1, Line3 line2)
        {
            //两条直线
            Vector3 v1 = Vector3.LineToVector(line1);
            Vector3 v2 = Vector3.LineToVector(line2);
            Vector3 v3 = Mul.GetXmul(v1, v2);
            //距离abs(dot(v3,p1)-dot(v3,p2))
            Point3 p1 = Vector3.GetProjP2L(line1.Point1, v3);
            Point3 p2 = Vector3.GetProjP2L(line2.Point1, v3);
            return PPdistance.GetPPdis3(p1.X, p1.Y, p1.Z, p2.X, p2.Y, p2.Z);
            /*
            Vector3 u = new Vector3(line1.Point1, line2.Point1, false);
            Vector3 v1 = new Vector3(line1.Point1, line1.Point2, false);
            Vector3 v2 = new Vector3(line2.Point1, line2.Point2, false);
            //Vector3 v = new Vector3();
            v1 = Vector3.GetUnitVector(v1);
            v2 = Vector3.GetUnitVector(v2);
            double a = Mul.GetPmul(v1, v1);
            double b = Mul.GetPmul(v1, v2);
            double c = Mul.GetPmul(v2, v2);
            double d = Mul.GetPmul(v1, u);
            double e = Mul.GetPmul(v2, u);
            double f = Mul.GetPmul(u, u);
            double det = a * c - b * b;
            double s = 0;
            double t = 0;
            if (det < Common.Eps)//平行
            {
                if (b > c)
                    t = d / b;
                else
                    t = e / c;
                return d * s + f;
            }
            else//不平行
            {
                double invDet = 1 / det;
                s = (b * e - c * d) * invDet;
                t = (a * e - b * d) * invDet;
                return s * (a * s + b * t + 2 * d) + t * (b * s - c * t + 2 * e) + f;
            }*/
        }
    }
    /// <summary>
    /// 点面距离
    /// </summary>
    public class PFdistance
    {
        /// <summary>
        /// 求点面距离，无返回垂足
        /// </summary>
        /// <param name="point"></param>
        /// <param name="face"></param>
        /// <returns></returns>
        public static double GetDistance(Point3 point, Face face)
        {
            return Math.Abs(face.A * point.X + face.B * point.Y + face.C * point.Z + face.D) / (Math.Sqrt(face.A * face.A + face.B * face.B + face.C * face.C));
        }
        /// <summary>
        /// 求点面距离，有返回垂足
        /// </summary>
        /// <param name="point"></param>
        /// <param name="face"></param>
        /// <param name="crossP"></param>
        /// <returns></returns>
        public static double GetDistance(Point3 point, Face face, out Point3 crossP)
        {
            float x0 = point.X;
            float y0 = point.Y;
            float z0 = point.Z;
            float x1 = face.P.X;
            float y1 = face.P.Y;
            float z1 = face.P.Z;
            float m = face.N.End.X;
            float n = face.N.End.Y;
            float p = face.N.End.Z;
            float t = ((x0 - x1) * m + (y0 - y1) * n + (z0 - z1) * p) / (m * m + n * n + p * p);
            float x2 = x0 + m * t;
            float y2 = x0 + n * t;
            float z2 = x0 + p * t;
            crossP = new Point3(x2, y2, z2);
            return PPdistance.GetPPdis3(x0, y0, z0, x2, y2, z2);
        }
    }
    /// <summary>
    /// 点线距离
    /// </summary>
    public class PLdistance
    {
        /// <summary>
        /// 求二维点线距
        /// </summary>
        /// <param name="point"></param>
        /// <param name="line"></param>
        /// <param name="crossP"></param>
        /// <returns></returns>
        public static double GetDistance2(Point2 point, Line2 line, out Point2 crossP)
        {//先求垂直平分线，再求两线的交点，最后是点点距离。
            Vector2 v = new Vector2(line.Point1, line.Point2);
            Line2 l = Vector2.GetCZAVGLine(line);
            Vector2 v1 = Vector2.Transform(Vector2.LineToVector(l), point);
            Point2 cross;
            Vector2.GetCross(v1, v, out cross);
            crossP = cross;
            return PPdistance.GetPPdis2(point.X, point.Y, cross.X, cross.Y);
        }
        /// <summary>
        /// 求三维点线距
        /// </summary>
        /// <param name="point"></param>
        /// <param name="line"></param>
        /// <param name="crossP"></param>
        /// <returns></returns>
        public static double GetDistance3(Point3 point, Line3 line, out Point3 crossP)
        {
            float x0 = point.X;
            float y0 = point.Y;
            float z0 = point.Z;
            float x1 = line.Point1.X;
            float y1 = line.Point1.Y;
            float z1 = line.Point1.Z;
            float m = line.Point2.X - line.Point1.X;
            float n = line.Point2.Y - line.Point1.Y;
            float p = line.Point2.Z - line.Point1.Z;
            float t = ((x0 - x1) * m + (y0 - y1) * n + (z0 - z1) * p) / (m * m + n * n + p * p);
            float x2 = x0 + m * t;
            float y2 = x0 + n * t;
            float z2 = x0 + p * t;
            crossP = new Point3(x2, y2, z2);
            return PPdistance.GetPPdis3(x0, y0, z0, x2, y2, z2);
        }
        /// <summary>
        /// 求三维点线的两点的距离
        /// </summary>
        /// <param name="point"></param>
        /// <param name="point1"></param>
        /// <param name="point2"></param>
        /// <param name="crossP"></param>
        /// <returns></returns>
        public static double GetDistance3(Point3 point, Point3 point1, Point3 point2, out Point3 crossP)
        {
            Line3 line12 = new Line3(point1, point2);
            bool inLine;
            Line3 line = Vector3.GetCZLine(point, new Line3(point1, point2), out inLine);
            
            float x0 = line.Point1.X;
            float y0 = line.Point1.Y;
            float z0 = line.Point1.Z;

            //(line12.Point2.X - line12.Point1.X) * (line.Point2.X - line.Point1.X) + (line12.Point2.Y - line12.Point1.Y) * (line.Point2.Y - line.Point1.Y) + (line12.Point2.Z - line12.Point1.Z) * (line.Point2.Z - line.Point1.Z) == 0;

            crossP = new Point3(x0, y0, z0);
            return PPdistance.GetPPdis3(x0, y0, z0, point.X, point.Y, point.Z);
        }
    }
    /// <summary>
    /// 点点距离
    /// </summary>
    public class PPdistance
    {
        /// <summary>
        /// 求二维点点距
        /// </summary>
        /// <param name="x1"></param>
        /// <param name="y1"></param>
        /// <param name="x2"></param>
        /// <param name="y2"></param>
        /// <returns></returns>
        public static double GetPPdis2(float x1, float y1, float x2, float y2)
        {
            return Math.Sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
        }
/// <summary>
        /// 求三维点点距
/// </summary>
/// <param name="x1"></param>
/// <param name="y1"></param>
/// <param name="z1"></param>
/// <param name="x2"></param>
/// <param name="y2"></param>
/// <param name="z2"></param>
/// <returns></returns>
        public static double GetPPdis3(float x1, float y1, float z1, float x2, float y2, float z2)
        {
            return Math.Sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
        }
    }
    /// <summary>
    /// 斜率类
    /// </summary>
    public class XL
    {
        /// <summary>
        /// 求二维斜率
        /// </summary>
        /// <param name="line"></param>
        /// <returns></returns>
        public static float Getxl2(Line2 line)
        {
            Point2 min = line.Point1;
            Point2 max = line.Point2;
            float xl = 0;
            if (Math.Abs(max.X - min.X) > Common.Eps)
            {
                xl = (max.Y - min.Y) / (max.X - min.X);
            }
            else
            {
                xl = float.MaxValue;
            }
            return xl;
        }
        /// <summary>
        /// 求三维斜率
        /// </summary>
        /// <param name="line"></param>
        /// <returns></returns>
        public static float[] Getxl3(Line3 line)
        {
            Point3 min = line.Point1;
            Point3 max = line.Point2;
            float[] xl = new float[2];
            float xl1 = 0, xl2 = 0;
            if (Math.Abs(max.X - min.X) > Common.Eps)
            {
                xl1 = (max.Y - min.Y) / (max.X - min.X);
                xl2 = (max.Z - min.Z) / (max.X - min.X);
            }
            else
            {
                xl1 = float.MaxValue;
                xl2 = float.MaxValue;
            }
            xl[0] = xl1;
            xl[1] = xl2;
            return xl;
        }
    }
    /// <summary>
    /// 向量方向类
    /// </summary>
    public class VectorDirection
    {
        /// <summary>
        /// 求二维直线的方向角[0,2*PI)，弧度为单位
        /// </summary>
        /// <param name="line"></param>
        /// <returns></returns>
        public static double GetDirection2(Line2 line)
        {
            Point2 min = line.Point1;
            Point2 max = line.Point2;
            double angle = 0;
            if (Math.Abs(max.X - min.X) > Common.Eps)
            {
                angle = Math.Atan2(min.Y - max.Y, max.X - min.X);
                angle = angle < 0 ? angle + 2 * Math.PI : angle;
            }
            else
            {
                if (min.Y - max.Y < 0)
                    angle = 3 * Math.PI / 2;
                else
                    angle = Math.PI / 2;
            }
            return angle;
        }
        /// <summary>
        /// 求三维直线的方向角[0,2*PI)，弧度为单位
        /// </summary>
        /// <param name="line"></param>
        /// <returns></returns>
        public static double[] GetDirection3(Line3 line)
        {
            Point3 min = line.Point1;
            Point3 max = line.Point2;
            double[] angle = new double[2];
            double angle1 = 0, angle2 = 0;
            if (Math.Abs(max.X - min.X) > Common.Eps)
            {
                angle1 = Math.Atan2(min.Y - max.Y, max.X - min.X);
                angle1 = angle1 < 0 ? angle1 + 2 * Math.PI : angle1;
                angle2 = Math.Atan2(min.Z - max.Z, max.X - min.X);
                angle2 = angle2 < 0 ? angle2 + 2 * Math.PI : angle2;
            }
            else
            {
                if (min.Y - max.Y < 0)
                    angle1 = 3 * Math.PI / 2;
                else
                    angle1 = Math.PI / 2;
                if (min.Z - max.Z < 0)
                    angle2 = 3 * Math.PI / 2;
                else
                    angle2 = Math.PI / 2;
            }
            angle[0] = angle1;
            angle[1] = angle2;
            return angle;
        }
        /// <summary>
        /// 求二维向量的方向角[0,2*PI)，弧度为单位
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static double GetDirection2(Vector2 vector)
        {
            Point2 min = vector.Start;
            Point2 max = vector.End;
            double angle = 0;
            if (Math.Abs(max.X - min.X) > Common.Eps)
            {
                angle = Math.Atan2(min.Y - max.Y, max.X - min.X);
                angle = angle < 0 ? angle + 2 * Math.PI : angle;
            }
            else
            {
                if (min.Y - max.Y < 0)
                    angle = 3 * Math.PI / 2;
                else
                    angle = Math.PI / 2;
            }
            return angle;
        }
        /// <summary>
        /// 求三维向量的方向角[0,2*PI)，弧度为单位
        /// </summary>
        /// <param name="vector"></param>
        /// <returns></returns>
        public static double[] GetDirection3(Vector3 vector)
        {
            Point3 min = vector.Start;
            Point3 max = vector.End;
            double[] angle = new double[2];
            double angle1 = 0, angle2 = 0;
            if (Math.Abs(max.X - min.X) > Common.Eps)
            {
                angle1 = Math.Atan2(min.Y - max.Y, max.X - min.X);
                angle1 = angle1 < 0 ? angle1 + 2 * Math.PI : angle1;
                angle2 = Math.Atan2(min.Z - max.Z, max.X - min.X);
                angle2 = angle2 < 0 ? angle2 + 2 * Math.PI : angle2;
            }
            else
            {
                if (min.Y - max.Y < 0)
                    angle1 = 3 * Math.PI / 2;
                else
                    angle1 = Math.PI / 2;
                if (min.Z - max.Z < 0)
                    angle2 = 3 * Math.PI / 2;
                else
                    angle2 = Math.PI / 2;
            }
            angle[0] = angle1;
            angle[1] = angle2;
            return angle;
        }
    }
}    