/*****************************************************************************
 * FileName:       mytypes.h
 * Author:         gotozw@yeah.net
 * Description:
 * Version:
 * Function List:
 *                 1.
 * History:
 *     <author>    <time>    <version >        <desc>
 *        ZW     2021-12-18     1.00      Created initial version
 *****************************************************************************/


#ifndef _MYTYPES_H
#define _MYTYPES_H

#ifdef __cplusplus
extern "C" {
#endif



//typedef enum {FALSE = 0, TRUE = !FALSE} BOOL;
//typedef enum {FAILED = 0, PASSED = !FAILED} TestStatus;


//typedef  uint8_t   			BYTE;
//typedef  uint16_t		  	WORD;
//typedef  uint32_t   		DWORD;
//typedef  BOOL           bool;

//typedef uint8_t             uint8;                   /* defined for unsigned 8-bits integer variable 	无符号8位整型变量  */
//typedef int8_t              int8;                    /* defined for signed 8-bits integer variable		有符号8位整型变量  */
//typedef uint16_t            uint16;                  /* defined for unsigned 16-bits integer variable 	无符号16位整型变量 */
//typedef int16_t             int16;                   /* defined for signed 16-bits integer variable 		有符号16位整型变量 */
//typedef uint32_t            uint32;                  /* defined for unsigned 32-bits integer variable 	无符号32位整型变量 */
//typedef int32_t             int32;                   /* defined for signed 32-bits integer variable 		有符号32位整型变量 */
//typedef float               fp32;                    /* single precision floating point variable (32bits) 单精度浮点数（32位长度） */
//typedef double              fp64;                    /* double precision floating point variable (64bits) 双精度浮点数（64位长度） */

//typedef uint8_t             u8;                   /* defined for unsigned 8-bits integer variable 	无符号8位整型变量  */
//typedef uint16_t            u16;                  /* defined for unsigned 16-bits integer variable 	无符号16位整型变量 */
//typedef uint32_t            u32;                  /* defined for unsigned 32-bits integer variable 	无符号32位整型变量 */


/*
* typedef struct
{
	float real;
	float image;
}COMPLEX_T;
*/


//#define REAL(x) ((x).real)
//#define IMAGE(x) ((x).image)

//#define FLOAT2INT(x)     ((*((int *)(&x))))
//#define FLOAT2UINT(x)   ((*((unsigned int *)(&x))))

//#define C_PI	3.141592653589793238462643
//#define C_DPR	180.0/C_PI	/*degree per radius*/
//#define C_RPD	C_PI/180.0	/*radius per degree*/
//#define PI2 	6.28318530717959

//求复数的矢量角
#define F_DEG(xVal)	(atan2(xVal.image, xVal.real)*C_DPR)

//#define F_MIN(a,b)	(((a) < (b))?(a):(b))
//#define F_MAX(a,b)	(((a) > (b))?(a):(b))

//求复数的幅值
#define F_AMP(xVal)	(sqrt((xVal.real)*(xVal.real)+(xVal.image)*(xVal.image)))
//#define F_AMP2(xVal)	((xVal.real)*(xVal.real)+(xVal.image)*(xVal.image))

//#define SAMFREQ 32
//#define FS  (SAMFREQ*50)



//#define  countof(a) (sizeof(a) / sizeof(*(a)))
//#define  LOBYTE(w)     ((BYTE)(w))             /* 0x1234----->0x34 */
//#define  HIBYTE(w)     ((BYTE)((WORD)(w) >> 8))   /* 0x1234----->0x12 */
//#define  LOWORD(l)     ((WORD)(l))                /* 0x12345678----->0x5678 */
//#define  HIWORD(l)     ((WORD)((DWORD)(l) >> 16))   /* 0x12345678----->0x1234 */
//#define  MAKEWORD(low,high)  ((WORD)((BYTE)(low)|(((WORD)((BYTE)(high)))<< 8)))    			/* high:12;  low:34----->0x1234 */
//#define  MAKEDWORD(low,high) ((DWORD)(((WORD)(low))|(((DWORD)((WORD)(high)))<< 16)))   /* high:1234;  low:5678----->0x12345678 */
//#define  ADJUSTWORD(x)   (MAKEWORD(HIBYTE(x),LOBYTE(x)))                                     /* 0x1234----->0x3412 */
//#define  ADJUSTDWORD(x)  (MAKEDWORD(ADJUSTWORD(HIWORD(x)),ADJUSTWORD(LOWORD(x))))         /*0x12345678----->0x78563412 */


//#define U8_TO_U16(ucH8, ucL8)   (uint16_t)((uint8_t)(ucH8)*0x100U+(uint8_t)(ucL8))

/*
    typedef union
    {
     struct {	BYTE low_byte;
                BYTE mlow_byte;
                BYTE mhigh_byte;
                BYTE high_byte;
             }float_byte;

     struct {	WORD low_word;
                WORD high_word;
             }float_word;
           float  	value;
           DWORD    value_u32;
    }FLOAT;
*/



#ifdef __cplusplus
}
#endif


#endif /*_MYTYPES_H*/
