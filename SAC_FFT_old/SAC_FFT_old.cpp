// ConsoleApplication1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "fft.h"


int main()
{
    int FS_FFT = 50 * SAMFREQ_FFT;
    arduinoFFT FFT = arduinoFFT(50 * 32); /* Create FFT object */

    COMPLEX_T a;
    a.real = 1;
    a.image = 1;
    std::cout << "--------------Test begin--------------------\n";
    std::cout << (F_AMP(a)) << "\n";
    std::cout << (F_DEG(a)) << "\n";
    std::cout << "--------------Test End----------------------\n";

    /* 构造2个函数，分别输入实部和虚部 */
    for (int j = 0; j < SAMFREQ_FFT; j++)
    {                               //直流70   400:50Hz                                         500:100Hz                                          600:150Hz
        FFT.g_ComputerFFT.faSam[j].real = 70 + 400 * cos(PI2 * j * 50.0 / FS_FFT + 0 * C_RPD) + 500 * cos(PI2 * j * 100.0 / FS_FFT + 40 * C_RPD) + 600 * cos(PI2 * j * 150.0 / FS_FFT - 50 * C_RPD);
        //直流10   20.5:50Hz                                          30.6:100Hz                                          40.7:150Hz
        FFT.g_ComputerFFT.faSam[j].image = 10 + 20.5 * cos(PI2 * j * 50.0 / FS_FFT + 30 * C_RPD) + 30.6 * cos(PI2 * j * 100.0 / FS_FFT + 40 * C_RPD) + 40.7 * cos(PI2 * j * 150.0 / FS_FFT - 50 * C_RPD);
        printf("%d      %f      %f \n", j, FFT.g_ComputerFFT.faSam[j].real, FFT.g_ComputerFFT.faSam[j].image);


    }
    FFT.FFT_Computer();

    printf("xiebo:i      U      Angle_U      I      Angle_I\n");
    int i = 0;
    for (int j = 0; j < SAMFREQ_FFT / 2; j++)
    {
        printf("%d      %f      %f", j, FFT.g_ComputerFFT.faValue_U[i][j], FFT.g_ComputerFFT.faAngle_U[i][j]);
        printf("      %f      %f\n", FFT.g_ComputerFFT.faValue_I[i][j], FFT.g_ComputerFFT.faAngle_I[i][j]);
    }



    printf("jun fang geng U:%f      %f \n", FFT.g_ComputerFFT.faRms_U[i], FFT.g_ComputerFFT.faRms_I[i]);
    std::cout << "程序执行结束\n";

}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
