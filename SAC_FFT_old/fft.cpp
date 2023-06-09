/*****************************************************************************
 * FileName:       fft.c
 * Author:         gotozw@yeah.net
 * Description:
 * Version:
 * Function List:
 *                 1.
 * History:
 *     <author>    <time>    <version >        <desc>
 *        ZW     2021-12-18     1.00      Created initial version
 *****************************************************************************/


/* Includes ------------------------------------------------------------------*/
#include "fft.h"
#include "stdio.h"
#include <string.h>
//#include <string.h>
//#include <iostream>

/*********************功能函数定义*************************/
//复数乘法函数
/*
void cplx_mul(cplx_t* x, cplx_t* y, cplx_t* r)
{
    r->real = x->real * y->real - x->image * y->image;
    r->image = x->real * y->image + x->image * y->real;
}
*/
//计算e的复数次方
/*
void cplx_exp(cplx_t *x, cplx_t *r)
{
    double expx = exp(x->real);
    r->real = expx * cos(x->image);
    r->image = expx * sin(x->image);
}
*/


//cos函数的选择函数，选择cos[l][k],复制到r函数里面去
void W_exp_32(u16 l , u16 k , cplx_t *r)
{
    switch(l)
    {
    case 1 :
        r->real = cos1;
        r->image = sin1;
        break;

    case 2:
        r->real = cos2[k];
        r->image = sin2[k];
        break;

    case 4:
        r->real = cos4[k];
        r->image = sin4[k];
        break;

    case 8:
        r->real = cos8[k];
        r->image = sin8[k];
        break;

    case 16:
        r->real = cos16[k];
        r->image = sin16[k];
        break;
    case 32:
        r->real = cos32[k];
        r->image = sin32[k];
        break;
    case 64:
        r->real = cos64[k];
        r->image = sin64[k];
        break;
    case 128:
        r->real = cos128[k];
        r->image = sin128[k];
        break;
    }
}

/**
  * @brief W_exp computer
  * @param l: 级数，必须为2的n次。如1,2,4,8,16.....
  * @param k: 0...l/2
  * @param r: W_exp   (out)
  * @retval None
  */

void W_exp(u16 l , u16 k , cplx_t *r)
{

    r->real = cos(k * C_PI / l); //  cos (k*2*pi/l) k=0..l/2
    r->image = -sin(k * C_PI / l); //  sin (k*2*pi/l) k=0..l/2

}



/*
void bit_reverse(cplx_t *x, u16 N)
{
    double t;
    cplx_t tmp;
    u16 i = 0, j = 0, k = 0;
    for(i = 0; i < N; i++)
    {
        k = i;
        j = 0;
        t = log(0.0 + N) / log(2.0);
        while((t--) > 0)
        {
            j <<= 1;
            j |= k & 1;
            k >>= 1;
        }
        if(j > i)
        {
            tmp = x[i];
            x[i] = x[j];
            x[j] = tmp;
        }
    }
}
*/


/*********************功能函数定义结束*************************/


//下面是类的函数
arduinoFFT::arduinoFFT(double samplingFrequency)
{// Constructor
//    this->_vReal = vReal;
//    this->_vImag = vImag;
//    this->_samples = samples;
//    this->_samplingFrequency = samplingFrequency;
//    this->_power = Exponent(samples);
    //this->g_ComputerFFT = input_ComputerFFT;
    this->FS_FFT = samplingFrequency;
}

arduinoFFT::~arduinoFFT(void)
{
    // Destructor
}


/**
  * @brief split two channel
  * @param x: fft_sam_all (in)
  * @param x: fft_out_all (out)
  * @param len: SAM_FREQ   (in)
  * @retval None
  */

void arduinoFFT::fft2( cplx_t *x, cplx_t *y, u16 len)
{
    u16 ex = 0, t = len;
    u16 i, j, k;
    cplx_t  w;
    float tr, ti, rr, ri, yr, yi;

    // if(NOT2POW(len)) return NULL;        //如果失败，返回空指针
    for(; !(t & 1); t >>= 1) ex++;  //len应该等于2的ex次方

    //y=(my_complex*)malloc(len*sizeof(my_complex));
    //if(!y) return NULL;

    //变址计算，库里-图基算法
    for(i = 0; i < len; i++)
    {
        k = i;
        j = 0;
        t = ex;
        while((t--) > 0)
        {
            j <<= 1;
            j |= k & 1;
            k >>= 1;
        }
        if(j >= i)
        {
            y[i] = x[j];
            y[j] = x[i];
        }
    }

    //用变址后的y向量进行计算
    for(i = 0; i < ex; i++)
    {
        t = 1 << i;
        for(j = 0; j < len; j += t << 1)
        {
            for(k = 0; k < t; k++)
            {
                //ti=-M_PI*k/t;
                //rr=cos(ti);
                //ri=sin(ti);
                W_exp( t , k , &w);
                rr = w.real;
                ri = w.image;

                tr = y[j+k+t].real;
                ti = y[j+k+t].image;

                yr = rr * tr - ri * ti;
                yi = rr * ti + ri * tr;

                tr = y[j+k].real;
                ti = y[j+k].image;

                y[j+k].real = tr + yr;
                y[j+k].image = ti + yi;
                y[j+k+t].real = tr - yr;
                y[j+k+t].image = ti - yi;
            }
        }
    }

    // return y;
}

/**
  * @brief split two channel
  * @param a: fft_out_all (in)
  * @param x: fft_real    (out)
  * @param y: fft_image   (out)
  * @param N: SAM_FREQ    (in)
  * @retval None
  */


void arduinoFFT::fft_split(cplx_t *a, cplx_t *x, cplx_t *y, u16 N)
{
    u16 i;
    x[0].real = a[0].real;
    x[0].image = 0;

    y[0].real = a[0].image;
    y[0].image = 0;
    for (i = 1; i < N; i++)
    {
        x[i].real = (a[i].real + a[N-i].real) / 2.0;
        x[i].image = (a[i].image - a[N-i].image) / 2.0;

        y[i].real = (a[i].image + a[N-i].image) / 2.0;
        y[i].image = (a[N-i].real - a[i].real) / 2.0;
    }
}

/**
  * @brief FFT test.
  * @param None
  *
  * @retval None
  */

void arduinoFFT::FFT_Computer(void)
{
    u16 j, i;
    //u16 count=pl20_sample; // 先保存指针
    static float faAngle;
    //g_Ad_using=FALSE;

    memcpy(&(this->SysPara), &DefaultSysPara, sizeof(this->SysPara));  //DefaultSysPara拷贝到SysPara


    for (i = 0; i < 1; i++)
    {
#if 0
        for (j = 0; j < SAMFREQ_FFT; j++)
        {
            this->g_ComputerFFT.faSam[j].real = nMeaBuf[i][j];         //电压
            this->g_ComputerFFT.faSam[j].image = nMeaBuf[i + 3][j];    //电流
        }
#endif
        //fft(faSam,SAMFREQ_FFT);
        /* 同时计算2个通道的FFT */
        fft2(this->g_ComputerFFT.faSam, this->g_ComputerFFT.faSam_out, SAMFREQ_FFT);

        /* 对计算结果按通道排序输出 */
        fft_split(this->g_ComputerFFT.faSam_out, this->g_ComputerFFT.faSam, this->g_ComputerFFT.faSam1, SAMFREQ_FFT);

        /* 计算各次谐波幅值 */

        for (j = 0; j < (SAMFREQ_FFT / 2); j++)
        {
            this->g_ComputerFFT.faValue_U[i][j] = this->SysPara.fUFFTCoef[i] * F_AMP(this->g_ComputerFFT.faSam[j]) / (SAMFREQ_FFT / 2);
            this->g_ComputerFFT.faValue_I[i][j] = this->SysPara.fIFFTCoef[i] * F_AMP(this->g_ComputerFFT.faSam1[j]) / (SAMFREQ_FFT / 2);
 

            if ((i == 0) && (j == 1))
            {
                faAngle = F_DEG(this->g_ComputerFFT.faSam[j]); // 计算基准相位
                //faAngle = 0; // 模拟时强制为0
            }
            this->g_ComputerFFT.faAngle_U[i][j] = F_DEG(this->g_ComputerFFT.faSam[j]) - faAngle;
            this->g_ComputerFFT.faAngle_I[i][j] = F_DEG(this->g_ComputerFFT.faSam1[j]) - faAngle;

        }
        /* 计算实部表示的函数直流幅值 */
        this->g_ComputerFFT.faValue_U[i][0] /= 2 * this->SysPara.fUFFTCoef[i]; //直流就是偏置电压，因此不需要乘以系数
        /* 计算虚部表示的函数直流幅值 */
        this->g_ComputerFFT.faValue_I[i][0] /= 2 * this->SysPara.fIFFTCoef[i];



        //printf("%d      %f      %f      %f      %f\n", j, this->g_ComputerFFT.faValue_U[i][j], this->g_ComputerFFT.faAngle_U[i][j], this->g_ComputerFFT.faValue_I[i][j], this->g_ComputerFFT.faAngle_I[i][j]);


        this->g_ComputerFFT.faThd_U[i] = this->g_ComputerFFT.faThd_I[i] = 0;
        for (j = 2; j < (SAMFREQ_FFT / 2); j++)
        {
            this->g_ComputerFFT.faThd_U[i] = this->g_ComputerFFT.faValue_U[i][j] * this->g_ComputerFFT.faValue_U[i][j] + this->g_ComputerFFT.faThd_U[i];
            this->g_ComputerFFT.faThd_I[i] = this->g_ComputerFFT.faValue_I[i][j] * this->g_ComputerFFT.faValue_I[i][j] + this->g_ComputerFFT.faThd_I[i];
        }
        /* 计算实部表示的函数均方根值 */
        this->g_ComputerFFT.faRms_U[i] = sqrt(this->g_ComputerFFT.faValue_U[i][1] * this->g_ComputerFFT.faValue_U[i][1] + this->g_ComputerFFT.faThd_U[i]);
        /* 计算虚部表示的函数均方根值 */
        this->g_ComputerFFT.faRms_I[i] = sqrt(this->g_ComputerFFT.faValue_I[i][1] * this->g_ComputerFFT.faValue_I[i][1] + this->g_ComputerFFT.faThd_I[i]);

        

        /* 计算实部表示的函数谐波含有率*/
        if (this->g_ComputerFFT.faValue_U[i][1] < 0.001)
            this->g_ComputerFFT.faThd_U[i] = 0;
        else
            this->g_ComputerFFT.faThd_U[i] = sqrt(this->g_ComputerFFT.faThd_U[i]) / this->g_ComputerFFT.faValue_U[i][1];
        /* 计算虚部表示的函数谐波含有率*/
        if (this->g_ComputerFFT.faValue_I[i][1] < 0.001)
            this->g_ComputerFFT.faThd_I[i] = 0;
        else
            this->g_ComputerFFT.faThd_I[i] = sqrt(this->g_ComputerFFT.faThd_I[i]) / this->g_ComputerFFT.faValue_I[i][1];

        printf("-----------------------------i=%d-------------------------------\n", i);
    }
    /*
    for (i = 0; i < 3; i++)
    {
        //直流谐波分量
        this->HWPercent[i][0] = this->HWPercent[i][1] = this->g_ComputerFFT.faThd_U[i];
        this->HWPercent[i + 3][0] = this->HWPercent[i + 3][1] = this->g_ComputerFFT.faThd_I[i];
    }
    for (i = 0; i < 3; i++)
    {
        if (this->g_ComputerFFT.faValue_U[i][1] < 0.001)
        {
            for (j = 2; j < (SAMFREQ_FFT / 2); j++)
            {
                this->HWPercent[i][j] = 0;
            }
        }
        else
        {
            for (j = 2; j < (SAMFREQ_FFT / 2); j++)
            {
                this->HWPercent[i][j] = this->g_ComputerFFT.faValue_U[i][j] / this->g_ComputerFFT.faValue_U[i][1];
            }
        }

        if (this->g_ComputerFFT.faValue_I[i][1] < 0.001)
        {
            for (j = 2; j < (SAMFREQ_FFT / 2); j++)
            {
                this->HWPercent[i + 3][j] = 0;
            }
        }
        else
        {
            for (j = 2; j < (SAMFREQ_FFT / 2); j++)
            {
                this->HWPercent[i + 3][j] = this->g_ComputerFFT.faValue_I[i][j] / this->g_ComputerFFT.faValue_I[i][1];
            }
        }
    }
    */
    //g_Ad_using=TRUE;
}



/******************* (C) COPYRIGHT 2010  *****END OF FILE****/
