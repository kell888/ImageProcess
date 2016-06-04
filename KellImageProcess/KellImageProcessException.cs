using System;
using System.Collections.Generic;
using System.Text;

namespace KellImageProcess
{
    /// <summary>
    /// 由本程序集引发的异常类
    /// </summary>
    public class KellImageProcessException : Exception
    {
        string _msg;
        /// <summary>
        /// 异常信息
        /// </summary>
        public string Msg
        {
            get { return _msg; }
        }
        /// <summary>
        /// 构造函数
        /// </summary>
        /// <param name="msg"></param>
        /// <param name="inner"></param>
        public KellImageProcessException(string msg, Exception inner = null)
            : base(msg, inner)
        {
            this._msg = msg;
        }
    }
}
