// SVDTest.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <cmath>
#include <intrin.h>
#include <iostream>
#include <iomanip> 


float rsqrt(const float f)
{
    float buf[4];
    buf[0] = f;
    __m128 v = _mm_loadu_ps(buf);
    v = _mm_rsqrt_ss(v);
    _mm_storeu_ps(buf, v);
    return buf[0];
}

int main(int argc, char* argv[])
{
    //if (argc != 2) { printf("Must specify an integer random seed as argument\n"); exit(1); }
    //srand(atoi(argv[1]));



    float test_A11, test_A21, test_A31, test_A12, test_A22, test_A32, test_A13, test_A23, test_A33;



    //test_A11 = 2. * (float)rand() / (float)0x7fff - 1.;
    //test_A21 = 2. * (float)rand() / (float)0x7fff - 1.;
    //test_A31 = 2. * (float)rand() / (float)0x7fff - 1.;
    //test_A12 = 2. * (float)rand() / (float)0x7fff - 1.;
    //test_A22 = 2. * (float)rand() / (float)0x7fff - 1.;
    //test_A32 = 2. * (float)rand() / (float)0x7fff - 1.;
    //test_A13 = 2. * (float)rand() / (float)0x7fff - 1.;
    //test_A23 = 2. * (float)rand() / (float)0x7fff - 1.;
    //test_A33 = 2. * (float)rand() / (float)0x7fff - 1.;
    //float norm_inverse = (float)(1. / sqrt((double)test_A11 * (double)test_A11 + (double)test_A21 * (double)test_A21 + (double)test_A31 * (double)test_A31
    //    + (double)test_A12 * (double)test_A12 + (double)test_A22 * (double)test_A22 + (double)test_A32 * (double)test_A32
    //    + (double)test_A13 * (double)test_A13 + (double)test_A23 * (double)test_A23 + (double)test_A33 * (double)test_A33));
    //test_A11 *= norm_inverse;
    //test_A21 *= norm_inverse;
    //test_A31 *= norm_inverse;
    //test_A12 *= norm_inverse;
    //test_A22 *= norm_inverse;
    //test_A32 *= norm_inverse;
    //test_A13 *= norm_inverse;
    //test_A23 *= norm_inverse;
    //test_A33 *= norm_inverse;



    test_A11 = 0; test_A12 = 0; test_A13 = 0;
    test_A21 = 0; test_A22 = 0; test_A23 = 0;
    test_A31 = 0; test_A32 = 0; test_A33 = 0;

    test_A11 = 1; test_A12 = 0; test_A13 = 0;
    test_A21 = 0; test_A22 = 1; test_A23 = 0;
    test_A31 = 0; test_A32 = 0; test_A33 = 1;


    test_A11 =  8.0; test_A12 = -6.0; test_A13 =  2.0;
    test_A21 = -6.0; test_A22 =  7.0; test_A23 = -4.0;
    test_A31 =  2.0; test_A32 = -4.0; test_A33 =  3.0;


    /* answers:

        U:
        0.66666667	-0.66666667	0.33333333
        -0.66666667	-0.33333333	0.66666667
        0.33333333	0.66666667	0.66666667

        Sigma:
        15	0	0
        0	3	0
        0	0	0

        V:
        0.66666667	-0.66666667	0.33333333
        -0.66666667	-0.33333333	0.66666667
        0.33333333	0.66666667	0.66666667

        or

        U:
        -0.66666667	-0.66666667	0.33333333
        0.66666667	-0.33333333	0.66666667
        -0.33333333	0.66666667	0.66666667

        Sigma:
        15	0	0
        0	3	0
        0	0	0

        V:
        -0.66666667	-0.66666667	0.33333333
        0.66666667	-0.33333333	0.66666667
        -0.33333333	0.66666667	0.66666667
    */







//#line 1 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Kernel_Declarations.hpp"


























































//#line 60 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Kernel_Declarations.hpp"

    const float Four_Gamma_Squared = sqrt(8.) + 3.;
    const float Sine_Pi_Over_Eight = .5 * sqrt(2. - sqrt(2.));
    const float Cosine_Pi_Over_Eight = .5 * sqrt(2. + sqrt(2.));

    union { float f; unsigned int ui; } Sfour_gamma_squared;
    union { float f; unsigned int ui; } Ssine_pi_over_eight;
    union { float f; unsigned int ui; } Scosine_pi_over_eight;
    union { float f; unsigned int ui; } Sone_half;
    union { float f; unsigned int ui; } Sone;
    union { float f; unsigned int ui; } Stiny_number;
    union { float f; unsigned int ui; } Ssmall_number;

    Sfour_gamma_squared.f = Four_Gamma_Squared;
    Ssine_pi_over_eight.f = Sine_Pi_Over_Eight;
    Scosine_pi_over_eight.f = Cosine_Pi_Over_Eight;
    Sone_half.f = .5;
    Sone.f = 1.;
    Stiny_number.f = 1.e-20;
    Ssmall_number.f = 1.e-12;

    union { float f; unsigned int ui; } Sa11;
    union { float f; unsigned int ui; } Sa21;
    union { float f; unsigned int ui; } Sa31;
    union { float f; unsigned int ui; } Sa12;
    union { float f; unsigned int ui; } Sa22;
    union { float f; unsigned int ui; } Sa32;
    union { float f; unsigned int ui; } Sa13;
    union { float f; unsigned int ui; } Sa23;
    union { float f; unsigned int ui; } Sa33;


    union { float f; unsigned int ui; } Sv11;
    union { float f; unsigned int ui; } Sv21;
    union { float f; unsigned int ui; } Sv31;
    union { float f; unsigned int ui; } Sv12;
    union { float f; unsigned int ui; } Sv22;
    union { float f; unsigned int ui; } Sv32;
    union { float f; unsigned int ui; } Sv13;
    union { float f; unsigned int ui; } Sv23;
    union { float f; unsigned int ui; } Sv33;
//#line 102 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Kernel_Declarations.hpp"


    union { float f; unsigned int ui; } Sqvs;
    union { float f; unsigned int ui; } Sqvvx;
    union { float f; unsigned int ui; } Sqvvy;
    union { float f; unsigned int ui; } Sqvvz;
//#line 109 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Kernel_Declarations.hpp"


    union { float f; unsigned int ui; } Su11;
    union { float f; unsigned int ui; } Su21;
    union { float f; unsigned int ui; } Su31;
    union { float f; unsigned int ui; } Su12;
    union { float f; unsigned int ui; } Su22;
    union { float f; unsigned int ui; } Su32;
    union { float f; unsigned int ui; } Su13;
    union { float f; unsigned int ui; } Su23;
    union { float f; unsigned int ui; } Su33;
//#line 121 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Kernel_Declarations.hpp"


    union { float f; unsigned int ui; } Squs;
    union { float f; unsigned int ui; } Squvx;
    union { float f; unsigned int ui; } Squvy;
    union { float f; unsigned int ui; } Squvz;
//#line 128 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Kernel_Declarations.hpp"

    union { float f; unsigned int ui; } Sc;
    union { float f; unsigned int ui; } Ss;
    union { float f; unsigned int ui; } Sch;
    union { float f; unsigned int ui; } Ssh;
    union { float f; unsigned int ui; } Stmp1;
    union { float f; unsigned int ui; } Stmp2;
    union { float f; unsigned int ui; } Stmp3;
    union { float f; unsigned int ui; } Stmp4;
    union { float f; unsigned int ui; } Stmp5;
//#line 63 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Unit_Test.cpp"

    Sa11.f = test_A11;
    Sa21.f = test_A21;
    Sa31.f = test_A31;
    Sa12.f = test_A12;
    Sa22.f = test_A22;
    Sa32.f = test_A32;
    Sa13.f = test_A13;
    Sa23.f = test_A23;
    Sa33.f = test_A33;

//#line 1 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"























    {








        {

            union { float f; unsigned int ui; } Ss11;
            union { float f; unsigned int ui; } Ss21;
            union { float f; unsigned int ui; } Ss31;
            union { float f; unsigned int ui; } Ss22;
            union { float f; unsigned int ui; } Ss32;
            union { float f; unsigned int ui; } Ss33;

            Sqvs.f = 1.;
            Sqvvx.f = 0.;
            Sqvvy.f = 0.;
            Sqvvz.f = 0.;





            Ss11.f = Sa11.f * Sa11.f;
            Stmp1.f = Sa21.f * Sa21.f;
            Ss11.f = Stmp1.f + Ss11.f;
            Stmp1.f = Sa31.f * Sa31.f;
            Ss11.f = Stmp1.f + Ss11.f;

            Ss21.f = Sa12.f * Sa11.f;
            Stmp1.f = Sa22.f * Sa21.f;
            Ss21.f = Stmp1.f + Ss21.f;
            Stmp1.f = Sa32.f * Sa31.f;
            Ss21.f = Stmp1.f + Ss21.f;

            Ss31.f = Sa13.f * Sa11.f;
            Stmp1.f = Sa23.f * Sa21.f;
            Ss31.f = Stmp1.f + Ss31.f;
            Stmp1.f = Sa33.f * Sa31.f;
            Ss31.f = Stmp1.f + Ss31.f;

            Ss22.f = Sa12.f * Sa12.f;
            Stmp1.f = Sa22.f * Sa22.f;
            Ss22.f = Stmp1.f + Ss22.f;
            Stmp1.f = Sa32.f * Sa32.f;
            Ss22.f = Stmp1.f + Ss22.f;

            Ss32.f = Sa13.f * Sa12.f;
            Stmp1.f = Sa23.f * Sa22.f;
            Ss32.f = Stmp1.f + Ss32.f;
            Stmp1.f = Sa33.f * Sa32.f;
            Ss32.f = Stmp1.f + Ss32.f;

            Ss33.f = Sa13.f * Sa13.f;
            Stmp1.f = Sa23.f * Sa23.f;
            Ss33.f = Stmp1.f + Ss33.f;
            Stmp1.f = Sa33.f * Sa33.f;
            Ss33.f = Stmp1.f + Ss33.f;





            for (int sweep = 1; sweep <= 4; sweep++) {





























//#line 1 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Jacobi_Conjugation_Kernel.hpp"




















                Ssh.f = Ss21.f * Sone_half.f;
                Stmp5.f = Ss11.f - Ss22.f;

                Stmp2.f = Ssh.f * Ssh.f;
                Stmp1.ui = (Stmp2.f >= Stiny_number.f) ? 0xffffffff : 0;
                Ssh.ui = Stmp1.ui & Ssh.ui;
                Sch.ui = Stmp1.ui & Stmp5.ui;
                Stmp2.ui = ~Stmp1.ui & Sone.ui;
                Sch.ui = Sch.ui | Stmp2.ui;

                Stmp1.f = Ssh.f * Ssh.f;
                Stmp2.f = Sch.f * Sch.f;
                Stmp3.f = Stmp1.f + Stmp2.f;
                Stmp4.f = rsqrt(Stmp3.f);










                Ssh.f = Stmp4.f * Ssh.f;
                Sch.f = Stmp4.f * Sch.f;

                Stmp1.f = Sfour_gamma_squared.f * Stmp1.f;
                Stmp1.ui = (Stmp2.f <= Stmp1.f) ? 0xffffffff : 0;

                Stmp2.ui = Ssine_pi_over_eight.ui & Stmp1.ui;
                Ssh.ui = ~Stmp1.ui & Ssh.ui;
                Ssh.ui = Ssh.ui | Stmp2.ui;
                Stmp2.ui = Scosine_pi_over_eight.ui & Stmp1.ui;
                Sch.ui = ~Stmp1.ui & Sch.ui;
                Sch.ui = Sch.ui | Stmp2.ui;

                Stmp1.f = Ssh.f * Ssh.f;
                Stmp2.f = Sch.f * Sch.f;
                Sc.f = Stmp2.f - Stmp1.f;
                Ss.f = Sch.f * Ssh.f;
                Ss.f = Ss.f + Ss.f;






                Stmp3.f = Stmp1.f + Stmp2.f;
                Ss33.f = Ss33.f * Stmp3.f;
                Ss31.f = Ss31.f * Stmp3.f;
                Ss32.f = Ss32.f * Stmp3.f;
                Ss33.f = Ss33.f * Stmp3.f;
//#line 75 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Jacobi_Conjugation_Kernel.hpp"

                Stmp1.f = Ss.f * Ss31.f;
                Stmp2.f = Ss.f * Ss32.f;
                Ss31.f = Sc.f * Ss31.f;
                Ss32.f = Sc.f * Ss32.f;
                Ss31.f = Stmp2.f + Ss31.f;
                Ss32.f = Ss32.f - Stmp1.f;

                Stmp2.f = Ss.f * Ss.f;
                Stmp1.f = Ss22.f * Stmp2.f;
                Stmp3.f = Ss11.f * Stmp2.f;
                Stmp4.f = Sc.f * Sc.f;
                Ss11.f = Ss11.f * Stmp4.f;
                Ss22.f = Ss22.f * Stmp4.f;
                Ss11.f = Ss11.f + Stmp1.f;
                Ss22.f = Ss22.f + Stmp3.f;
                Stmp4.f = Stmp4.f - Stmp2.f;
                Stmp2.f = Ss21.f + Ss21.f;
                Ss21.f = Ss21.f * Stmp4.f;
                Stmp4.f = Sc.f * Ss.f;
                Stmp2.f = Stmp2.f * Stmp4.f;
                Stmp5.f = Stmp5.f * Stmp4.f;
                Ss11.f = Ss11.f + Stmp2.f;
                Ss21.f = Ss21.f - Stmp5.f;
                Ss22.f = Ss22.f - Stmp2.f;





                Stmp1.f = Ssh.f * Sqvvx.f;
                Stmp2.f = Ssh.f * Sqvvy.f;
                Stmp3.f = Ssh.f * Sqvvz.f;
                Ssh.f = Ssh.f * Sqvs.f;

                Sqvs.f = Sch.f * Sqvs.f;
                Sqvvx.f = Sch.f * Sqvvx.f;
                Sqvvy.f = Sch.f * Sqvvy.f;
                Sqvvz.f = Sch.f * Sqvvz.f;

                Sqvvz.f = Sqvvz.f + Ssh.f;
                Sqvs.f = Sqvs.f - Stmp3.f;
                Sqvvx.f = Sqvvx.f + Stmp2.f;
                Sqvvy.f = Sqvvy.f - Stmp1.f;
//#line 122 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"























































//#line 1 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Jacobi_Conjugation_Kernel.hpp"




















                Ssh.f = Ss32.f * Sone_half.f;
                Stmp5.f = Ss22.f - Ss33.f;

                Stmp2.f = Ssh.f * Ssh.f;
                Stmp1.ui = (Stmp2.f >= Stiny_number.f) ? 0xffffffff : 0;
                Ssh.ui = Stmp1.ui & Ssh.ui;
                Sch.ui = Stmp1.ui & Stmp5.ui;
                Stmp2.ui = ~Stmp1.ui & Sone.ui;
                Sch.ui = Sch.ui | Stmp2.ui;

                Stmp1.f = Ssh.f * Ssh.f;
                Stmp2.f = Sch.f * Sch.f;
                Stmp3.f = Stmp1.f + Stmp2.f;
                Stmp4.f = rsqrt(Stmp3.f);










                Ssh.f = Stmp4.f * Ssh.f;
                Sch.f = Stmp4.f * Sch.f;

                Stmp1.f = Sfour_gamma_squared.f * Stmp1.f;
                Stmp1.ui = (Stmp2.f <= Stmp1.f) ? 0xffffffff : 0;

                Stmp2.ui = Ssine_pi_over_eight.ui & Stmp1.ui;
                Ssh.ui = ~Stmp1.ui & Ssh.ui;
                Ssh.ui = Ssh.ui | Stmp2.ui;
                Stmp2.ui = Scosine_pi_over_eight.ui & Stmp1.ui;
                Sch.ui = ~Stmp1.ui & Sch.ui;
                Sch.ui = Sch.ui | Stmp2.ui;

                Stmp1.f = Ssh.f * Ssh.f;
                Stmp2.f = Sch.f * Sch.f;
                Sc.f = Stmp2.f - Stmp1.f;
                Ss.f = Sch.f * Ssh.f;
                Ss.f = Ss.f + Ss.f;






                Stmp3.f = Stmp1.f + Stmp2.f;
                Ss11.f = Ss11.f * Stmp3.f;
                Ss21.f = Ss21.f * Stmp3.f;
                Ss31.f = Ss31.f * Stmp3.f;
                Ss11.f = Ss11.f * Stmp3.f;
//#line 75 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Jacobi_Conjugation_Kernel.hpp"

                Stmp1.f = Ss.f * Ss21.f;
                Stmp2.f = Ss.f * Ss31.f;
                Ss21.f = Sc.f * Ss21.f;
                Ss31.f = Sc.f * Ss31.f;
                Ss21.f = Stmp2.f + Ss21.f;
                Ss31.f = Ss31.f - Stmp1.f;

                Stmp2.f = Ss.f * Ss.f;
                Stmp1.f = Ss33.f * Stmp2.f;
                Stmp3.f = Ss22.f * Stmp2.f;
                Stmp4.f = Sc.f * Sc.f;
                Ss22.f = Ss22.f * Stmp4.f;
                Ss33.f = Ss33.f * Stmp4.f;
                Ss22.f = Ss22.f + Stmp1.f;
                Ss33.f = Ss33.f + Stmp3.f;
                Stmp4.f = Stmp4.f - Stmp2.f;
                Stmp2.f = Ss32.f + Ss32.f;
                Ss32.f = Ss32.f * Stmp4.f;
                Stmp4.f = Sc.f * Ss.f;
                Stmp2.f = Stmp2.f * Stmp4.f;
                Stmp5.f = Stmp5.f * Stmp4.f;
                Ss22.f = Ss22.f + Stmp2.f;
                Ss32.f = Ss32.f - Stmp5.f;
                Ss33.f = Ss33.f - Stmp2.f;





                Stmp1.f = Ssh.f * Sqvvx.f;
                Stmp2.f = Ssh.f * Sqvvy.f;
                Stmp3.f = Ssh.f * Sqvvz.f;
                Ssh.f = Ssh.f * Sqvs.f;

                Sqvs.f = Sch.f * Sqvs.f;
                Sqvvx.f = Sch.f * Sqvvx.f;
                Sqvvy.f = Sch.f * Sqvvy.f;
                Sqvvz.f = Sch.f * Sqvvz.f;

                Sqvvx.f = Sqvvx.f + Ssh.f;
                Sqvs.f = Sqvs.f - Stmp1.f;
                Sqvvy.f = Sqvvy.f + Stmp3.f;
                Sqvvz.f = Sqvvz.f - Stmp2.f;
//#line 178 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"























































//#line 1 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Jacobi_Conjugation_Kernel.hpp"




















                Ssh.f = Ss31.f * Sone_half.f;
                Stmp5.f = Ss33.f - Ss11.f;

                Stmp2.f = Ssh.f * Ssh.f;
                Stmp1.ui = (Stmp2.f >= Stiny_number.f) ? 0xffffffff : 0;
                Ssh.ui = Stmp1.ui & Ssh.ui;
                Sch.ui = Stmp1.ui & Stmp5.ui;
                Stmp2.ui = ~Stmp1.ui & Sone.ui;
                Sch.ui = Sch.ui | Stmp2.ui;

                Stmp1.f = Ssh.f * Ssh.f;
                Stmp2.f = Sch.f * Sch.f;
                Stmp3.f = Stmp1.f + Stmp2.f;
                Stmp4.f = rsqrt(Stmp3.f);










                Ssh.f = Stmp4.f * Ssh.f;
                Sch.f = Stmp4.f * Sch.f;

                Stmp1.f = Sfour_gamma_squared.f * Stmp1.f;
                Stmp1.ui = (Stmp2.f <= Stmp1.f) ? 0xffffffff : 0;

                Stmp2.ui = Ssine_pi_over_eight.ui & Stmp1.ui;
                Ssh.ui = ~Stmp1.ui & Ssh.ui;
                Ssh.ui = Ssh.ui | Stmp2.ui;
                Stmp2.ui = Scosine_pi_over_eight.ui & Stmp1.ui;
                Sch.ui = ~Stmp1.ui & Sch.ui;
                Sch.ui = Sch.ui | Stmp2.ui;

                Stmp1.f = Ssh.f * Ssh.f;
                Stmp2.f = Sch.f * Sch.f;
                Sc.f = Stmp2.f - Stmp1.f;
                Ss.f = Sch.f * Ssh.f;
                Ss.f = Ss.f + Ss.f;






                Stmp3.f = Stmp1.f + Stmp2.f;
                Ss22.f = Ss22.f * Stmp3.f;
                Ss32.f = Ss32.f * Stmp3.f;
                Ss21.f = Ss21.f * Stmp3.f;
                Ss22.f = Ss22.f * Stmp3.f;
//#line 75 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Jacobi_Conjugation_Kernel.hpp"

                Stmp1.f = Ss.f * Ss32.f;
                Stmp2.f = Ss.f * Ss21.f;
                Ss32.f = Sc.f * Ss32.f;
                Ss21.f = Sc.f * Ss21.f;
                Ss32.f = Stmp2.f + Ss32.f;
                Ss21.f = Ss21.f - Stmp1.f;

                Stmp2.f = Ss.f * Ss.f;
                Stmp1.f = Ss11.f * Stmp2.f;
                Stmp3.f = Ss33.f * Stmp2.f;
                Stmp4.f = Sc.f * Sc.f;
                Ss33.f = Ss33.f * Stmp4.f;
                Ss11.f = Ss11.f * Stmp4.f;
                Ss33.f = Ss33.f + Stmp1.f;
                Ss11.f = Ss11.f + Stmp3.f;
                Stmp4.f = Stmp4.f - Stmp2.f;
                Stmp2.f = Ss31.f + Ss31.f;
                Ss31.f = Ss31.f * Stmp4.f;
                Stmp4.f = Sc.f * Ss.f;
                Stmp2.f = Stmp2.f * Stmp4.f;
                Stmp5.f = Stmp5.f * Stmp4.f;
                Ss33.f = Ss33.f + Stmp2.f;
                Ss31.f = Ss31.f - Stmp5.f;
                Ss11.f = Ss11.f - Stmp2.f;





                Stmp1.f = Ssh.f * Sqvvx.f;
                Stmp2.f = Ssh.f * Sqvvy.f;
                Stmp3.f = Ssh.f * Sqvvz.f;
                Ssh.f = Ssh.f * Sqvs.f;

                Sqvs.f = Sch.f * Sqvs.f;
                Sqvvx.f = Sch.f * Sqvvx.f;
                Sqvvy.f = Sch.f * Sqvvy.f;
                Sqvvz.f = Sch.f * Sqvvz.f;

                Sqvvy.f = Sqvvy.f + Ssh.f;
                Sqvs.f = Sqvs.f - Stmp2.f;
                Sqvvz.f = Sqvvz.f + Stmp1.f;
                Sqvvx.f = Sqvvx.f - Stmp3.f;
//#line 234 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"


























            }



            std::cout << "Scalar S =" << std::endl;
            std::cout << std::setw(12) << Ss11.f << std::endl;
            std::cout << std::setw(12) << Ss21.f << "  " << std::setw(12) << Ss22.f << std::endl;
            std::cout << std::setw(12) << Ss31.f << "  " << std::setw(12) << Ss32.f << "  " << std::setw(12) << Ss33.f << std::endl;
//#line 269 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"
























//#line 294 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"

        }







        Stmp2.f = Sqvs.f * Sqvs.f;
        Stmp1.f = Sqvvx.f * Sqvvx.f;
        Stmp2.f = Stmp1.f + Stmp2.f;
        Stmp1.f = Sqvvy.f * Sqvvy.f;
        Stmp2.f = Stmp1.f + Stmp2.f;
        Stmp1.f = Sqvvz.f * Sqvvz.f;
        Stmp2.f = Stmp1.f + Stmp2.f;

        Stmp1.f = rsqrt(Stmp2.f);
        Stmp4.f = Stmp1.f * Sone_half.f;
        Stmp3.f = Stmp1.f * Stmp4.f;
        Stmp3.f = Stmp1.f * Stmp3.f;
        Stmp3.f = Stmp2.f * Stmp3.f;
        Stmp1.f = Stmp1.f + Stmp4.f;
        Stmp1.f = Stmp1.f - Stmp3.f;

        Sqvs.f = Sqvs.f * Stmp1.f;
        Sqvvx.f = Sqvvx.f * Stmp1.f;
        Sqvvy.f = Sqvvy.f * Stmp1.f;
        Sqvvz.f = Sqvvz.f * Stmp1.f;



        std::cout << "Scalar qV =" << std::endl;
        std::cout << std::setw(12) << Sqvs.f << "  " << std::setw(12) << Sqvvx.f << "  " << std::setw(12) << Sqvvy.f << "  " << std::setw(12) << Sqvvz.f << std::endl;
//#line 329 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"
















//#line 346 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"

//#line 348 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"

        {

















            Stmp1.f = Sqvvx.f * Sqvvx.f;
            Stmp2.f = Sqvvy.f * Sqvvy.f;
            Stmp3.f = Sqvvz.f * Sqvvz.f;
            Sv11.f = Sqvs.f * Sqvs.f;
            Sv22.f = Sv11.f - Stmp1.f;
            Sv33.f = Sv22.f - Stmp2.f;
            Sv33.f = Sv33.f + Stmp3.f;
            Sv22.f = Sv22.f + Stmp2.f;
            Sv22.f = Sv22.f - Stmp3.f;
            Sv11.f = Sv11.f + Stmp1.f;
            Sv11.f = Sv11.f - Stmp2.f;
            Sv11.f = Sv11.f - Stmp3.f;
            Stmp1.f = Sqvvx.f + Sqvvx.f;
            Stmp2.f = Sqvvy.f + Sqvvy.f;
            Stmp3.f = Sqvvz.f + Sqvvz.f;
            Sv32.f = Sqvs.f * Stmp1.f;
            Sv13.f = Sqvs.f * Stmp2.f;
            Sv21.f = Sqvs.f * Stmp3.f;
            Stmp1.f = Sqvvy.f * Stmp1.f;
            Stmp2.f = Sqvvz.f * Stmp2.f;
            Stmp3.f = Sqvvx.f * Stmp3.f;
            Sv12.f = Stmp1.f - Sv21.f;
            Sv23.f = Stmp2.f - Sv32.f;
            Sv31.f = Stmp3.f - Sv13.f;
            Sv21.f = Stmp1.f + Sv21.f;
            Sv32.f = Stmp2.f + Sv32.f;
            Sv13.f = Stmp3.f + Sv13.f;




            std::cout << "Scalar V =" << std::endl;
            std::cout << std::setw(12) << Sv11.f << "  " << std::setw(12) << Sv12.f << "  " << std::setw(12) << Sv13.f << std::endl;
            std::cout << std::setw(12) << Sv21.f << "  " << std::setw(12) << Sv22.f << "  " << std::setw(12) << Sv23.f << std::endl;
            std::cout << std::setw(12) << Sv31.f << "  " << std::setw(12) << Sv32.f << "  " << std::setw(12) << Sv33.f << std::endl;
//#line 403 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"






























//#line 434 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"
//#line 435 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"





            Stmp2.f = Sa12.f;
            Stmp3.f = Sa13.f;
            Sa12.f = Sv12.f * Sa11.f;
            Sa13.f = Sv13.f * Sa11.f;
            Sa11.f = Sv11.f * Sa11.f;
            Stmp1.f = Sv21.f * Stmp2.f;
            Sa11.f = Sa11.f + Stmp1.f;
            Stmp1.f = Sv31.f * Stmp3.f;
            Sa11.f = Sa11.f + Stmp1.f;
            Stmp1.f = Sv22.f * Stmp2.f;
            Sa12.f = Sa12.f + Stmp1.f;
            Stmp1.f = Sv32.f * Stmp3.f;
            Sa12.f = Sa12.f + Stmp1.f;
            Stmp1.f = Sv23.f * Stmp2.f;
            Sa13.f = Sa13.f + Stmp1.f;
            Stmp1.f = Sv33.f * Stmp3.f;
            Sa13.f = Sa13.f + Stmp1.f;

            Stmp2.f = Sa22.f;
            Stmp3.f = Sa23.f;
            Sa22.f = Sv12.f * Sa21.f;
            Sa23.f = Sv13.f * Sa21.f;
            Sa21.f = Sv11.f * Sa21.f;
            Stmp1.f = Sv21.f * Stmp2.f;
            Sa21.f = Sa21.f + Stmp1.f;
            Stmp1.f = Sv31.f * Stmp3.f;
            Sa21.f = Sa21.f + Stmp1.f;
            Stmp1.f = Sv22.f * Stmp2.f;
            Sa22.f = Sa22.f + Stmp1.f;
            Stmp1.f = Sv32.f * Stmp3.f;
            Sa22.f = Sa22.f + Stmp1.f;
            Stmp1.f = Sv23.f * Stmp2.f;
            Sa23.f = Sa23.f + Stmp1.f;
            Stmp1.f = Sv33.f * Stmp3.f;
            Sa23.f = Sa23.f + Stmp1.f;

            Stmp2.f = Sa32.f;
            Stmp3.f = Sa33.f;
            Sa32.f = Sv12.f * Sa31.f;
            Sa33.f = Sv13.f * Sa31.f;
            Sa31.f = Sv11.f * Sa31.f;
            Stmp1.f = Sv21.f * Stmp2.f;
            Sa31.f = Sa31.f + Stmp1.f;
            Stmp1.f = Sv31.f * Stmp3.f;
            Sa31.f = Sa31.f + Stmp1.f;
            Stmp1.f = Sv22.f * Stmp2.f;
            Sa32.f = Sa32.f + Stmp1.f;
            Stmp1.f = Sv32.f * Stmp3.f;
            Sa32.f = Sa32.f + Stmp1.f;
            Stmp1.f = Sv23.f * Stmp2.f;
            Sa33.f = Sa33.f + Stmp1.f;
            Stmp1.f = Sv33.f * Stmp3.f;
            Sa33.f = Sa33.f + Stmp1.f;



            std::cout << "Scalar A (after multiplying with V) =" << std::endl;
            std::cout << std::setw(12) << Sa11.f << "  " << std::setw(12) << Sa12.f << "  " << std::setw(12) << Sa13.f << std::endl;
            std::cout << std::setw(12) << Sa21.f << "  " << std::setw(12) << Sa22.f << "  " << std::setw(12) << Sa23.f << std::endl;
            std::cout << std::setw(12) << Sa31.f << "  " << std::setw(12) << Sa32.f << "  " << std::setw(12) << Sa33.f << std::endl;
//#line 501 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"






























//#line 532 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"

        }

    }





    Stmp1.f = Sa11.f * Sa11.f;
    Stmp4.f = Sa21.f * Sa21.f;
    Stmp1.f = Stmp1.f + Stmp4.f;
    Stmp4.f = Sa31.f * Sa31.f;
    Stmp1.f = Stmp1.f + Stmp4.f;

    Stmp2.f = Sa12.f * Sa12.f;
    Stmp4.f = Sa22.f * Sa22.f;
    Stmp2.f = Stmp2.f + Stmp4.f;
    Stmp4.f = Sa32.f * Sa32.f;
    Stmp2.f = Stmp2.f + Stmp4.f;

    Stmp3.f = Sa13.f * Sa13.f;
    Stmp4.f = Sa23.f * Sa23.f;
    Stmp3.f = Stmp3.f + Stmp4.f;
    Stmp4.f = Sa33.f * Sa33.f;
    Stmp3.f = Stmp3.f + Stmp4.f;



    Stmp4.ui = (Stmp1.f < Stmp2.f) ? 0xffffffff : 0;
    Stmp5.ui = Sa11.ui ^ Sa12.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sa11.ui = Sa11.ui ^ Stmp5.ui;
    Sa12.ui = Sa12.ui ^ Stmp5.ui;

    Stmp5.ui = Sa21.ui ^ Sa22.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sa21.ui = Sa21.ui ^ Stmp5.ui;
    Sa22.ui = Sa22.ui ^ Stmp5.ui;

    Stmp5.ui = Sa31.ui ^ Sa32.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sa31.ui = Sa31.ui ^ Stmp5.ui;
    Sa32.ui = Sa32.ui ^ Stmp5.ui;


    Stmp5.ui = Sv11.ui ^ Sv12.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sv11.ui = Sv11.ui ^ Stmp5.ui;
    Sv12.ui = Sv12.ui ^ Stmp5.ui;

    Stmp5.ui = Sv21.ui ^ Sv22.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sv21.ui = Sv21.ui ^ Stmp5.ui;
    Sv22.ui = Sv22.ui ^ Stmp5.ui;

    Stmp5.ui = Sv31.ui ^ Sv32.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sv31.ui = Sv31.ui ^ Stmp5.ui;
    Sv32.ui = Sv32.ui ^ Stmp5.ui;
//#line 593 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"

    Stmp5.ui = Stmp1.ui ^ Stmp2.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Stmp1.ui = Stmp1.ui ^ Stmp5.ui;
    Stmp2.ui = Stmp2.ui ^ Stmp5.ui;



    Stmp5.f = -2.;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Stmp4.f = 1.;
    Stmp4.f = Stmp4.f + Stmp5.f;

    Sa12.f = Sa12.f * Stmp4.f;
    Sa22.f = Sa22.f * Stmp4.f;
    Sa32.f = Sa32.f * Stmp4.f;


    Sv12.f = Sv12.f * Stmp4.f;
    Sv22.f = Sv22.f * Stmp4.f;
    Sv32.f = Sv32.f * Stmp4.f;
//#line 615 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"




    Stmp4.f = Stmp4.f * Sone_half.f;
    Stmp4.f = Stmp4.f - Sone_half.f;

    Stmp5.f = Stmp4.f * Sqvvz.f;
    Stmp5.f = Stmp5.f + Sqvs.f;
    Sqvs.f = Sqvs.f * Stmp4.f;
    Sqvvz.f = Sqvvz.f - Sqvs.f;
    Sqvs.f = Stmp5.f;

    Stmp5.f = Stmp4.f * Sqvvx.f;
    Stmp5.f = Stmp5.f + Sqvvy.f;
    Sqvvy.f = Sqvvy.f * Stmp4.f;
    Sqvvx.f = Sqvvx.f - Sqvvy.f;
    Sqvvy.f = Stmp5.f;
//#line 634 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"



    Stmp4.ui = (Stmp1.f < Stmp3.f) ? 0xffffffff : 0;
    Stmp5.ui = Sa11.ui ^ Sa13.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sa11.ui = Sa11.ui ^ Stmp5.ui;
    Sa13.ui = Sa13.ui ^ Stmp5.ui;

    Stmp5.ui = Sa21.ui ^ Sa23.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sa21.ui = Sa21.ui ^ Stmp5.ui;
    Sa23.ui = Sa23.ui ^ Stmp5.ui;

    Stmp5.ui = Sa31.ui ^ Sa33.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sa31.ui = Sa31.ui ^ Stmp5.ui;
    Sa33.ui = Sa33.ui ^ Stmp5.ui;


    Stmp5.ui = Sv11.ui ^ Sv13.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sv11.ui = Sv11.ui ^ Stmp5.ui;
    Sv13.ui = Sv13.ui ^ Stmp5.ui;

    Stmp5.ui = Sv21.ui ^ Sv23.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sv21.ui = Sv21.ui ^ Stmp5.ui;
    Sv23.ui = Sv23.ui ^ Stmp5.ui;

    Stmp5.ui = Sv31.ui ^ Sv33.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sv31.ui = Sv31.ui ^ Stmp5.ui;
    Sv33.ui = Sv33.ui ^ Stmp5.ui;
//#line 669 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"

    Stmp5.ui = Stmp1.ui ^ Stmp3.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Stmp1.ui = Stmp1.ui ^ Stmp5.ui;
    Stmp3.ui = Stmp3.ui ^ Stmp5.ui;



    Stmp5.f = -2.;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Stmp4.f = 1.;
    Stmp4.f = Stmp4.f + Stmp5.f;

    Sa11.f = Sa11.f * Stmp4.f;
    Sa21.f = Sa21.f * Stmp4.f;
    Sa31.f = Sa31.f * Stmp4.f;


    Sv11.f = Sv11.f * Stmp4.f;
    Sv21.f = Sv21.f * Stmp4.f;
    Sv31.f = Sv31.f * Stmp4.f;
//#line 691 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"




    Stmp4.f = Stmp4.f * Sone_half.f;
    Stmp4.f = Stmp4.f - Sone_half.f;

    Stmp5.f = Stmp4.f * Sqvvy.f;
    Stmp5.f = Stmp5.f + Sqvs.f;
    Sqvs.f = Sqvs.f * Stmp4.f;
    Sqvvy.f = Sqvvy.f - Sqvs.f;
    Sqvs.f = Stmp5.f;

    Stmp5.f = Stmp4.f * Sqvvz.f;
    Stmp5.f = Stmp5.f + Sqvvx.f;
    Sqvvx.f = Sqvvx.f * Stmp4.f;
    Sqvvz.f = Sqvvz.f - Sqvvx.f;
    Sqvvx.f = Stmp5.f;
//#line 710 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"



    Stmp4.ui = (Stmp2.f < Stmp3.f) ? 0xffffffff : 0;
    Stmp5.ui = Sa12.ui ^ Sa13.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sa12.ui = Sa12.ui ^ Stmp5.ui;
    Sa13.ui = Sa13.ui ^ Stmp5.ui;

    Stmp5.ui = Sa22.ui ^ Sa23.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sa22.ui = Sa22.ui ^ Stmp5.ui;
    Sa23.ui = Sa23.ui ^ Stmp5.ui;

    Stmp5.ui = Sa32.ui ^ Sa33.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sa32.ui = Sa32.ui ^ Stmp5.ui;
    Sa33.ui = Sa33.ui ^ Stmp5.ui;


    Stmp5.ui = Sv12.ui ^ Sv13.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sv12.ui = Sv12.ui ^ Stmp5.ui;
    Sv13.ui = Sv13.ui ^ Stmp5.ui;

    Stmp5.ui = Sv22.ui ^ Sv23.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sv22.ui = Sv22.ui ^ Stmp5.ui;
    Sv23.ui = Sv23.ui ^ Stmp5.ui;

    Stmp5.ui = Sv32.ui ^ Sv33.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Sv32.ui = Sv32.ui ^ Stmp5.ui;
    Sv33.ui = Sv33.ui ^ Stmp5.ui;
//#line 745 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"

    Stmp5.ui = Stmp2.ui ^ Stmp3.ui;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Stmp2.ui = Stmp2.ui ^ Stmp5.ui;
    Stmp3.ui = Stmp3.ui ^ Stmp5.ui;



    Stmp5.f = -2.;
    Stmp5.ui = Stmp5.ui & Stmp4.ui;
    Stmp4.f = 1.;
    Stmp4.f = Stmp4.f + Stmp5.f;

    Sa13.f = Sa13.f * Stmp4.f;
    Sa23.f = Sa23.f * Stmp4.f;
    Sa33.f = Sa33.f * Stmp4.f;


    Sv13.f = Sv13.f * Stmp4.f;
    Sv23.f = Sv23.f * Stmp4.f;
    Sv33.f = Sv33.f * Stmp4.f;
//#line 767 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"




    Stmp4.f = Stmp4.f * Sone_half.f;
    Stmp4.f = Stmp4.f - Sone_half.f;

    Stmp5.f = Stmp4.f * Sqvvx.f;
    Stmp5.f = Stmp5.f + Sqvs.f;
    Sqvs.f = Sqvs.f * Stmp4.f;
    Sqvvx.f = Sqvvx.f - Sqvs.f;
    Sqvs.f = Stmp5.f;

    Stmp5.f = Stmp4.f * Sqvvy.f;
    Stmp5.f = Stmp5.f + Sqvvz.f;
    Sqvvz.f = Sqvvz.f * Stmp4.f;
    Sqvvy.f = Sqvvy.f - Sqvvz.f;
    Sqvvz.f = Stmp5.f;
//#line 786 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"




    std::cout << "Scalar V =" << std::endl;
    std::cout << std::setw(12) << Sv11.f << "  " << std::setw(12) << Sv12.f << "  " << std::setw(12) << Sv13.f << std::endl;
    std::cout << std::setw(12) << Sv21.f << "  " << std::setw(12) << Sv22.f << "  " << std::setw(12) << Sv23.f << std::endl;
    std::cout << std::setw(12) << Sv31.f << "  " << std::setw(12) << Sv32.f << "  " << std::setw(12) << Sv33.f << std::endl;
//#line 795 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"






























//#line 826 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"
//#line 827 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"



    std::cout << "Scalar A (after multiplying with V) =" << std::endl;
    std::cout << std::setw(12) << Sa11.f << "  " << std::setw(12) << Sa12.f << "  " << std::setw(12) << Sa13.f << std::endl;
    std::cout << std::setw(12) << Sa21.f << "  " << std::setw(12) << Sa22.f << "  " << std::setw(12) << Sa23.f << std::endl;
    std::cout << std::setw(12) << Sa31.f << "  " << std::setw(12) << Sa32.f << "  " << std::setw(12) << Sa33.f << std::endl;
//#line 835 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"






























//#line 866 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"






    Stmp2.f = Sqvs.f * Sqvs.f;
    Stmp1.f = Sqvvx.f * Sqvvx.f;
    Stmp2.f = Stmp1.f + Stmp2.f;
    Stmp1.f = Sqvvy.f * Sqvvy.f;
    Stmp2.f = Stmp1.f + Stmp2.f;
    Stmp1.f = Sqvvz.f * Sqvvz.f;
    Stmp2.f = Stmp1.f + Stmp2.f;
    Stmp1.f = rsqrt(Stmp2.f);










    Sqvs.f = Sqvs.f * Stmp1.f;
    Sqvvx.f = Sqvvx.f * Stmp1.f;
    Sqvvy.f = Sqvvy.f * Stmp1.f;
    Sqvvz.f = Sqvvz.f * Stmp1.f;



    std::cout << "Scalar qV =" << std::endl;
    std::cout << std::setw(12) << Sqvs.f << "  " << std::setw(12) << Sqvvx.f << "  " << std::setw(12) << Sqvvy.f << "  " << std::setw(12) << Sqvvz.f << std::endl;
//#line 900 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"
















//#line 917 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"
//#line 918 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"






    Su11.f = 1.;
    Su21.f = 0.;
    Su31.f = 0.;
    Su12.f = 0.;
    Su22.f = 1.;
    Su32.f = 0.;
    Su13.f = 0.;
    Su23.f = 0.;
    Su33.f = 1.;
//#line 934 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"


    Squs.f = 1.;
    Squvx.f = 0.;
    Squvy.f = 0.;
    Squvz.f = 0.;
//#line 941 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"

































//#line 1 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Givens_QR_Factorization_Kernel.hpp"




















    Ssh.f = Sa21.f * Sa21.f;
    Ssh.ui = (Ssh.f >= Ssmall_number.f) ? 0xffffffff : 0;
    Ssh.ui = Ssh.ui & Sa21.ui;

    Stmp5.f = 0.;
    Sch.f = Stmp5.f - Sa11.f;
    Sch.f = std::max(Sch.f, Sa11.f);
    Sch.f = std::max(Sch.f, Ssmall_number.f);
    Stmp5.ui = (Sa11.f >= Stmp5.f) ? 0xffffffff : 0;

    Stmp1.f = Sch.f * Sch.f;
    Stmp2.f = Ssh.f * Ssh.f;
    Stmp2.f = Stmp1.f + Stmp2.f;
    Stmp1.f = rsqrt(Stmp2.f);

    Stmp4.f = Stmp1.f * Sone_half.f;
    Stmp3.f = Stmp1.f * Stmp4.f;
    Stmp3.f = Stmp1.f * Stmp3.f;
    Stmp3.f = Stmp2.f * Stmp3.f;
    Stmp1.f = Stmp1.f + Stmp4.f;
    Stmp1.f = Stmp1.f - Stmp3.f;
    Stmp1.f = Stmp1.f * Stmp2.f;

    Sch.f = Sch.f + Stmp1.f;

    Stmp1.ui = ~Stmp5.ui & Ssh.ui;
    Stmp2.ui = ~Stmp5.ui & Sch.ui;
    Sch.ui = Stmp5.ui & Sch.ui;
    Ssh.ui = Stmp5.ui & Ssh.ui;
    Sch.ui = Sch.ui | Stmp1.ui;
    Ssh.ui = Ssh.ui | Stmp2.ui;

    Stmp1.f = Sch.f * Sch.f;
    Stmp2.f = Ssh.f * Ssh.f;
    Stmp2.f = Stmp1.f + Stmp2.f;
    Stmp1.f = rsqrt(Stmp2.f);

    Stmp4.f = Stmp1.f * Sone_half.f;
    Stmp3.f = Stmp1.f * Stmp4.f;
    Stmp3.f = Stmp1.f * Stmp3.f;
    Stmp3.f = Stmp2.f * Stmp3.f;
    Stmp1.f = Stmp1.f + Stmp4.f;
    Stmp1.f = Stmp1.f - Stmp3.f;

    Sch.f = Sch.f * Stmp1.f;
    Ssh.f = Ssh.f * Stmp1.f;

    Sc.f = Sch.f * Sch.f;
    Ss.f = Ssh.f * Ssh.f;
    Sc.f = Sc.f - Ss.f;
    Ss.f = Ssh.f * Sch.f;
    Ss.f = Ss.f + Ss.f;





    Stmp1.f = Ss.f * Sa11.f;
    Stmp2.f = Ss.f * Sa21.f;
    Sa11.f = Sc.f * Sa11.f;
    Sa21.f = Sc.f * Sa21.f;
    Sa11.f = Sa11.f + Stmp2.f;
    Sa21.f = Sa21.f - Stmp1.f;

    Stmp1.f = Ss.f * Sa12.f;
    Stmp2.f = Ss.f * Sa22.f;
    Sa12.f = Sc.f * Sa12.f;
    Sa22.f = Sc.f * Sa22.f;
    Sa12.f = Sa12.f + Stmp2.f;
    Sa22.f = Sa22.f - Stmp1.f;

    Stmp1.f = Ss.f * Sa13.f;
    Stmp2.f = Ss.f * Sa23.f;
    Sa13.f = Sc.f * Sa13.f;
    Sa23.f = Sc.f * Sa23.f;
    Sa13.f = Sa13.f + Stmp2.f;
    Sa23.f = Sa23.f - Stmp1.f;






    Stmp1.f = Ss.f * Su11.f;
    Stmp2.f = Ss.f * Su12.f;
    Su11.f = Sc.f * Su11.f;
    Su12.f = Sc.f * Su12.f;
    Su11.f = Su11.f + Stmp2.f;
    Su12.f = Su12.f - Stmp1.f;

    Stmp1.f = Ss.f * Su21.f;
    Stmp2.f = Ss.f * Su22.f;
    Su21.f = Sc.f * Su21.f;
    Su22.f = Sc.f * Su22.f;
    Su21.f = Su21.f + Stmp2.f;
    Su22.f = Su22.f - Stmp1.f;

    Stmp1.f = Ss.f * Su31.f;
    Stmp2.f = Ss.f * Su32.f;
    Su31.f = Sc.f * Su31.f;
    Su32.f = Sc.f * Su32.f;
    Su31.f = Su31.f + Stmp2.f;
    Su32.f = Su32.f - Stmp1.f;
//#line 125 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Givens_QR_Factorization_Kernel.hpp"
//#line 975 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"


































    Squs.f = Sch.f;
    Squvx.f = 0.;
    Squvy.f = 0.;
    Squvz.f = Ssh.f;
//#line 1014 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"

































//#line 1 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Givens_QR_Factorization_Kernel.hpp"




















    Ssh.f = Sa31.f * Sa31.f;
    Ssh.ui = (Ssh.f >= Ssmall_number.f) ? 0xffffffff : 0;
    Ssh.ui = Ssh.ui & Sa31.ui;

    Stmp5.f = 0.;
    Sch.f = Stmp5.f - Sa11.f;
    Sch.f = std::max(Sch.f, Sa11.f);
    Sch.f = std::max(Sch.f, Ssmall_number.f);
    Stmp5.ui = (Sa11.f >= Stmp5.f) ? 0xffffffff : 0;

    Stmp1.f = Sch.f * Sch.f;
    Stmp2.f = Ssh.f * Ssh.f;
    Stmp2.f = Stmp1.f + Stmp2.f;
    Stmp1.f = rsqrt(Stmp2.f);

    Stmp4.f = Stmp1.f * Sone_half.f;
    Stmp3.f = Stmp1.f * Stmp4.f;
    Stmp3.f = Stmp1.f * Stmp3.f;
    Stmp3.f = Stmp2.f * Stmp3.f;
    Stmp1.f = Stmp1.f + Stmp4.f;
    Stmp1.f = Stmp1.f - Stmp3.f;
    Stmp1.f = Stmp1.f * Stmp2.f;

    Sch.f = Sch.f + Stmp1.f;

    Stmp1.ui = ~Stmp5.ui & Ssh.ui;
    Stmp2.ui = ~Stmp5.ui & Sch.ui;
    Sch.ui = Stmp5.ui & Sch.ui;
    Ssh.ui = Stmp5.ui & Ssh.ui;
    Sch.ui = Sch.ui | Stmp1.ui;
    Ssh.ui = Ssh.ui | Stmp2.ui;

    Stmp1.f = Sch.f * Sch.f;
    Stmp2.f = Ssh.f * Ssh.f;
    Stmp2.f = Stmp1.f + Stmp2.f;
    Stmp1.f = rsqrt(Stmp2.f);

    Stmp4.f = Stmp1.f * Sone_half.f;
    Stmp3.f = Stmp1.f * Stmp4.f;
    Stmp3.f = Stmp1.f * Stmp3.f;
    Stmp3.f = Stmp2.f * Stmp3.f;
    Stmp1.f = Stmp1.f + Stmp4.f;
    Stmp1.f = Stmp1.f - Stmp3.f;

    Sch.f = Sch.f * Stmp1.f;
    Ssh.f = Ssh.f * Stmp1.f;

    Sc.f = Sch.f * Sch.f;
    Ss.f = Ssh.f * Ssh.f;
    Sc.f = Sc.f - Ss.f;
    Ss.f = Ssh.f * Sch.f;
    Ss.f = Ss.f + Ss.f;





    Stmp1.f = Ss.f * Sa11.f;
    Stmp2.f = Ss.f * Sa31.f;
    Sa11.f = Sc.f * Sa11.f;
    Sa31.f = Sc.f * Sa31.f;
    Sa11.f = Sa11.f + Stmp2.f;
    Sa31.f = Sa31.f - Stmp1.f;

    Stmp1.f = Ss.f * Sa12.f;
    Stmp2.f = Ss.f * Sa32.f;
    Sa12.f = Sc.f * Sa12.f;
    Sa32.f = Sc.f * Sa32.f;
    Sa12.f = Sa12.f + Stmp2.f;
    Sa32.f = Sa32.f - Stmp1.f;

    Stmp1.f = Ss.f * Sa13.f;
    Stmp2.f = Ss.f * Sa33.f;
    Sa13.f = Sc.f * Sa13.f;
    Sa33.f = Sc.f * Sa33.f;
    Sa13.f = Sa13.f + Stmp2.f;
    Sa33.f = Sa33.f - Stmp1.f;






    Stmp1.f = Ss.f * Su11.f;
    Stmp2.f = Ss.f * Su13.f;
    Su11.f = Sc.f * Su11.f;
    Su13.f = Sc.f * Su13.f;
    Su11.f = Su11.f + Stmp2.f;
    Su13.f = Su13.f - Stmp1.f;

    Stmp1.f = Ss.f * Su21.f;
    Stmp2.f = Ss.f * Su23.f;
    Su21.f = Sc.f * Su21.f;
    Su23.f = Sc.f * Su23.f;
    Su21.f = Su21.f + Stmp2.f;
    Su23.f = Su23.f - Stmp1.f;

    Stmp1.f = Ss.f * Su31.f;
    Stmp2.f = Ss.f * Su33.f;
    Su31.f = Sc.f * Su31.f;
    Su33.f = Sc.f * Su33.f;
    Su31.f = Su31.f + Stmp2.f;
    Su33.f = Su33.f - Stmp1.f;
//#line 125 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Givens_QR_Factorization_Kernel.hpp"
//#line 1048 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"


































    Squvx.f = Ssh.f * Squvz.f;
    Ssh.f = Ssh.f * Squs.f;
    Squvy.f = Squvy.f - Ssh.f;
    Squs.f = Sch.f * Squs.f;
    Squvz.f = Sch.f * Squvz.f;
//#line 1088 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"

































//#line 1 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Givens_QR_Factorization_Kernel.hpp"




















    Ssh.f = Sa32.f * Sa32.f;
    Ssh.ui = (Ssh.f >= Ssmall_number.f) ? 0xffffffff : 0;
    Ssh.ui = Ssh.ui & Sa32.ui;

    Stmp5.f = 0.;
    Sch.f = Stmp5.f - Sa22.f;
    Sch.f = std::max(Sch.f, Sa22.f);
    Sch.f = std::max(Sch.f, Ssmall_number.f);
    Stmp5.ui = (Sa22.f >= Stmp5.f) ? 0xffffffff : 0;

    Stmp1.f = Sch.f * Sch.f;
    Stmp2.f = Ssh.f * Ssh.f;
    Stmp2.f = Stmp1.f + Stmp2.f;
    Stmp1.f = rsqrt(Stmp2.f);

    Stmp4.f = Stmp1.f * Sone_half.f;
    Stmp3.f = Stmp1.f * Stmp4.f;
    Stmp3.f = Stmp1.f * Stmp3.f;
    Stmp3.f = Stmp2.f * Stmp3.f;
    Stmp1.f = Stmp1.f + Stmp4.f;
    Stmp1.f = Stmp1.f - Stmp3.f;
    Stmp1.f = Stmp1.f * Stmp2.f;

    Sch.f = Sch.f + Stmp1.f;

    Stmp1.ui = ~Stmp5.ui & Ssh.ui;
    Stmp2.ui = ~Stmp5.ui & Sch.ui;
    Sch.ui = Stmp5.ui & Sch.ui;
    Ssh.ui = Stmp5.ui & Ssh.ui;
    Sch.ui = Sch.ui | Stmp1.ui;
    Ssh.ui = Ssh.ui | Stmp2.ui;

    Stmp1.f = Sch.f * Sch.f;
    Stmp2.f = Ssh.f * Ssh.f;
    Stmp2.f = Stmp1.f + Stmp2.f;
    Stmp1.f = rsqrt(Stmp2.f);

    Stmp4.f = Stmp1.f * Sone_half.f;
    Stmp3.f = Stmp1.f * Stmp4.f;
    Stmp3.f = Stmp1.f * Stmp3.f;
    Stmp3.f = Stmp2.f * Stmp3.f;
    Stmp1.f = Stmp1.f + Stmp4.f;
    Stmp1.f = Stmp1.f - Stmp3.f;

    Sch.f = Sch.f * Stmp1.f;
    Ssh.f = Ssh.f * Stmp1.f;

    Sc.f = Sch.f * Sch.f;
    Ss.f = Ssh.f * Ssh.f;
    Sc.f = Sc.f - Ss.f;
    Ss.f = Ssh.f * Sch.f;
    Ss.f = Ss.f + Ss.f;





    Stmp1.f = Ss.f * Sa21.f;
    Stmp2.f = Ss.f * Sa31.f;
    Sa21.f = Sc.f * Sa21.f;
    Sa31.f = Sc.f * Sa31.f;
    Sa21.f = Sa21.f + Stmp2.f;
    Sa31.f = Sa31.f - Stmp1.f;

    Stmp1.f = Ss.f * Sa22.f;
    Stmp2.f = Ss.f * Sa32.f;
    Sa22.f = Sc.f * Sa22.f;
    Sa32.f = Sc.f * Sa32.f;
    Sa22.f = Sa22.f + Stmp2.f;
    Sa32.f = Sa32.f - Stmp1.f;

    Stmp1.f = Ss.f * Sa23.f;
    Stmp2.f = Ss.f * Sa33.f;
    Sa23.f = Sc.f * Sa23.f;
    Sa33.f = Sc.f * Sa33.f;
    Sa23.f = Sa23.f + Stmp2.f;
    Sa33.f = Sa33.f - Stmp1.f;






    Stmp1.f = Ss.f * Su12.f;
    Stmp2.f = Ss.f * Su13.f;
    Su12.f = Sc.f * Su12.f;
    Su13.f = Sc.f * Su13.f;
    Su12.f = Su12.f + Stmp2.f;
    Su13.f = Su13.f - Stmp1.f;

    Stmp1.f = Ss.f * Su22.f;
    Stmp2.f = Ss.f * Su23.f;
    Su22.f = Sc.f * Su22.f;
    Su23.f = Sc.f * Su23.f;
    Su22.f = Su22.f + Stmp2.f;
    Su23.f = Su23.f - Stmp1.f;

    Stmp1.f = Ss.f * Su32.f;
    Stmp2.f = Ss.f * Su33.f;
    Su32.f = Sc.f * Su32.f;
    Su33.f = Sc.f * Su33.f;
    Su32.f = Su32.f + Stmp2.f;
    Su33.f = Su33.f - Stmp1.f;
//#line 125 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Givens_QR_Factorization_Kernel.hpp"
//#line 1122 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"


































    Stmp1.f = Ssh.f * Squvx.f;
    Stmp2.f = Ssh.f * Squvy.f;
    Stmp3.f = Ssh.f * Squvz.f;
    Ssh.f = Ssh.f * Squs.f;
    Squs.f = Sch.f * Squs.f;
    Squvx.f = Sch.f * Squvx.f;
    Squvy.f = Sch.f * Squvy.f;
    Squvz.f = Sch.f * Squvz.f;
    Squvx.f = Squvx.f + Ssh.f;
    Squs.f = Squs.f - Stmp1.f;
    Squvy.f = Squvy.f + Stmp3.f;
    Squvz.f = Squvz.f - Stmp2.f;
//#line 1169 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"




    std::cout << "Scalar U =" << std::endl;
    std::cout << std::setw(12) << Su11.f << "  " << std::setw(12) << Su12.f << "  " << std::setw(12) << Su13.f << std::endl;
    std::cout << std::setw(12) << Su21.f << "  " << std::setw(12) << Su22.f << "  " << std::setw(12) << Su23.f << std::endl;
    std::cout << std::setw(12) << Su31.f << "  " << std::setw(12) << Su32.f << "  " << std::setw(12) << Su33.f << std::endl;
//#line 1178 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"






























//#line 1209 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"
//#line 1210 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"



    std::cout << "Scalar A (after multiplying with U-transpose and V) =" << std::endl;
    std::cout << std::setw(12) << Sa11.f << "  " << std::setw(12) << Sa12.f << "  " << std::setw(12) << Sa13.f << std::endl;
    std::cout << std::setw(12) << Sa21.f << "  " << std::setw(12) << Sa22.f << "  " << std::setw(12) << Sa23.f << std::endl;
    std::cout << std::setw(12) << Sa31.f << "  " << std::setw(12) << Sa32.f << "  " << std::setw(12) << Sa33.f << std::endl;
//#line 1218 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"






























//#line 1249 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"




    std::cout << "Scalar qU =" << std::endl;
    std::cout << std::setw(12) << Squs.f << "  " << std::setw(12) << Squvx.f << "  " << std::setw(12) << Squvy.f << "  " << std::setw(12) << Squvz.f << std::endl;
//#line 1256 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"
















//#line 1273 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"
//#line 1274 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Main_Kernel_Body.hpp"




//#line 75 "C:\\Users\\fatfatzhang\\source\\repos\\SVD\\SVD\\Singular_Value_Decomposition_Unit_Test.cpp"

    return 0;
}

