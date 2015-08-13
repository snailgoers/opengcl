//
//  ImageEdge.cpp
//  Kernel
//
//  Created by JiaYonglei on 15/8/11.
//  Copyright (c) 2015年 JiaYonglei. All rights reserved.
//

#include <stdio.h>
#include "Image.h"



void BwThreshold(const uchar *data, uchar *out, int height, int width, double threhold, bool Binary_Inv)
{
    assert(data != NULL && width >= 0 && height >= 0 && threhold >= 0 && threhold <= 255);
    memset(out, 0, width * height);
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (data[i * width + j] > threhold)
            {
                out[i * width + j] = uchar(255);
            }
        }
    }
    if (Binary_Inv)
    {
        for (int i = 0; i < height; i++)
        {
            for (int j = 0; j < width; j++)
            {
                out[i * width + j] = uchar(255 - out[i * width + j]);
            }
        }
    }
}
///大津法阈值分割
void Otsu(uchar *data, uchar *bw, int width, int height, bool Binary_Inv)
{
    assert(data != NULL && bw != NULL);
    int *hist = new int[256];
    memset(hist, 0, sizeof(int) * 256);
    int total = width * height;
    for (int i = 0; i < total; i++)
    {
        hist[data[i]]++;
    }
    double N0 = 0;
    double N1 = 0;
    double w0, w1;
    double u0, u1;
    double gSum0, gSum1;
    double g;
    double maxG = -DBL_MAX;
    int threshold = 128;
    for (int i = 0; i < 256; i++)
    {
        N0 += hist[i];
        N1 = total - N0;
        w0 = N0 / total;
        w1 = 1 - w0;
        gSum0 = 0;
        gSum1 = 0;
        for (int j = 0; j < i; j++)
        {
            gSum0 += hist[j] * j;
        }
        u0 = gSum0 / N0;
        for (int k = i + 1; k < 256; k++)
        {
            gSum1 += hist[k] * k;
        }
        u1 = gSum1 / N1;
        g = w0 * w1 * (u0 - u1) * (u0 - u1);
        if (g > maxG)
        {
            maxG = g;
            threshold = i;
        }
    }
    BwThreshold(data, bw, height, width, threshold, Binary_Inv);
    return;
}
///自适应阈值分割
void AdaptationThreshold(const uchar *data, uchar *bw, int width, int height, bool Binary_Inv)
{
    assert(data != NULL && width >= 0 && height >= 0);
    int maxpixel = 0;
    int minmixel = 255;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (data[i * width + j] > maxpixel)
            {
                maxpixel = data[i * width + j];
            }
            if (data[i * width + j] < minmixel)
            {
                minmixel = data[i * width + j];
            }
        }
    }
    if (maxpixel == minmixel)
    {
        memcpy(bw, data, width * height);
    }
    else
    {
        double inithreshold = double(maxpixel + minmixel) / 2.0;
        // threshold the gray imag
        int flag = 0;
        while (flag != 1)
        {
            double sumwrite = 0;
            double sumblack = 0;
            int wsum = 0;
            int bsum = 0;
            // threshold the imag
            BwThreshold(data, bw, width, height, inithreshold, Binary_Inv);
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    if (bw[i * width + j] == 255)
                    {
                        sumwrite += data[i * width + j];
                        wsum++;
                    }
                    else
                    {
                        sumblack += data[i * width + j];
                        bsum++;
                    }
                }
            }
            double Tnext = (sumwrite / wsum + sumblack / bsum) / 2;
            //printf("%f \n",abs(inithreshold-Tnext));
            if (abs(inithreshold - Tnext) < 0.5)
            {
                flag = 1;
            }
            inithreshold = Tnext;
        }
    }
}


/// Gauss 边缘检测
void EdgeDetect_Gauss(uchar *data, double *edgeX, double *edgeY, double *edgeXY, double *theta, int width, int height, double sigma)
{
    memset(edgeX, 0, width * height);
    memset(edgeY, 0, width * height);
    assert(data != NULL);
    int len = MinGaussSize;
    // create the filter windows
    double *kernelX = new double[len * len];
    double *kernelY = new double[len * len];
    int halfLen = len / 2;

    for (int i = 0; i < len; i++)
    {
        for (int j = 0; j < len; j++)
        {
            int x = j - halfLen;
            int y = i - halfLen;
            double temp = -x * exp(-(x * x + y * y) / (2 * sigma * sigma)) / (M_PI * sigma * sigma);
            kernelX[i * len + j] = temp;
            kernelY[j * len + i] = temp;
        }
    }
    ImageFilter2(data, edgeX, height, width, kernelX, MinGaussSize);
    ImageFilter2(data, edgeY, height, width, kernelY, MinGaussSize);
    for (int i = 0; i < height * width;  i++) {
        edgeXY[i] = sqrt(edgeX[i] * edgeX[i] + edgeY[i] * edgeY[i]);
        theta[i] = atan2(edgeY[i], edgeX[i]);
    }
    
   
    delete [] kernelX;
    delete [] kernelY;
}
/// 创建高斯核
double* CreateGaussKernel(double sigma, int &size)
{
    size = MinGaussSize;
    int center = size / 2;
    double *kernel = new double[size];
    memset(kernel, 0, size * sizeof(double));
    double sumValue = 0.0;
    for (int i = 0; i < size; i++)
    {
        double dis = double(i - center);
        double dValue = exp(-dis * dis / (2 * sigma * sigma)) / (sqrt(2 * M_PI) * sigma);
        kernel[i] = dValue;
        sumValue += dValue;
    }
    // normalization
    for (int i = 0; i < size; i++)
    {
        kernel[i] /= sumValue;
    }
    return kernel;
}
/// 高斯平滑滤波
void GaussFilter(const uchar *data, uchar *out, int width, int height, double sigma)
{
    // 固定大小9X9
    assert(data != NULL);
    memcpy(out, data, width * height);
    // 创建高斯核
    int size = 9;
    double *kernel = CreateGaussKernel(sigma, size);
    int halfSize = size / 2;
    double * pdTemp = new double[width * height];
    // X filtting
    for (int i = halfSize; i < height - halfSize; i++)
    {
        for (int j = halfSize; j < width - halfSize; j++)
        {
            pdTemp[i * width + j] = data[i * width + j - 4] * kernel[0] + data[i * width + j - 3] * kernel[1] + data[i * width + j - 2] * kernel[2] + data[i * width + j - 1] * kernel[3] + data[i * width + j + 0] * kernel[4] + data[i * width + j + 1] * kernel[5] + data[i * width + j + 2] * kernel[6] + data[i * width + j + 3] * kernel[7] + data[i * width + j + 4] * kernel[8];
        }
    }
    // Y fillting
    for (int i = halfSize; i < width - halfSize; i++)
    {
        for (int j = halfSize; j < height - halfSize; j++)
        {
            double temp = pdTemp[(j - 4) * width + i] * kernel[0] + pdTemp[(j - 3) * width + i] * kernel[1] + pdTemp[(j - 2) * width + i] * kernel[2] + pdTemp[(j - 1) * width + i] * kernel[3] + pdTemp[(j + 0) * width + i] * kernel[4] + pdTemp[(j + 1) * width + i] * kernel[5] + pdTemp[(j + 2) * width + i] * kernel[6] + pdTemp[(j + 3) * width + j] * kernel[7] + pdTemp[(j + 4) * width + i] * kernel[8];
            out[j * width + i] = uchar(temp);
        }
    }
    delete[] kernel;
    delete[] pdTemp;
}
/// 自动确定canny边缘检测高低阈值
void AutoCalHighLowThreshold(const double *data, double *dataNormal, int width, int height,
                             double &highThreshold, double &lowThreshold)
{
    assert(data != NULL);
    double maxData = -DBL_MAX;
    memcpy(dataNormal, data, width * height * sizeof(double));
    for (int i = 0; i < height * width; i++)
    {
        if (data[i] > maxData)
        {
            maxData = data[i];
        }
    }
    // normlization
    for (int i = 0; i < height * width; i++)
    {
        dataNormal[i] /= maxData;
    }
    // 统计直方图
    const int binNum = 64;
    double *bin = new double[binNum];
    memset(bin, 0, binNum * sizeof(double));
    for (int i = 0; i < height * width; i++)
    {
        bin[(int)(dataNormal[i] * (binNum - 1))]++;
    }
    double sumTemp = 0;
    double thresholdScale = 0.7; // 确定高阈值用
    double lowHightScale = 0.4;
    double highDex = 0;
    for (int i = 0; i < binNum; i++)
    {
        sumTemp += bin[i];
        if (sumTemp > thresholdScale * height * width)
        {
            highDex = (double)i;
            break;
        }
    }
    highThreshold = (highDex + 1) / binNum;
    lowThreshold = lowHightScale * highThreshold;
    delete[] bin;
}
/// 非极大抑制
void NonMaximumInhibiion(double *edgeX, double *edgeY, double *edgeXYNormal, uchar *nonMaximum,
                         int width, int height, double lowThreshold, double highThreshold)
{
    assert(edgeX != NULL && edgeY != NULL && edgeXYNormal != NULL);
    memset(nonMaximum, 0, width * height);
    for (int i = 1; i < height - 1; i++)
    {
        for (int j = 1; j < width - 1; j++)
        {
            double dx = edgeX[i * width + j];
            double dy = edgeY[i * width + j];
            double gradMag1 = 0, gradMag2 = 0;
            double gradMag = edgeXYNormal[i * width + j];
            if ((dy < 0 && dx > -dy) || (dy > 0 && dx < -dy))
            {
                // 0 ~ 45°
                double d = abs(dy / dx);
                gradMag1 = edgeXYNormal[i * width + j + 1] * (1 - d) + edgeXYNormal[(i - 1) * width + j + 1] * d;
                gradMag2 = edgeXYNormal[i * width + j - 1] * (1 - d) + edgeXYNormal[(i + 1) * width + j - 1] * d;
            }
            else if ((dx > 0 && -dy >= dx) || (dx < 0 && -dy <= dx))
            {
                // 45 ~ 90°
                double d = abs(dx / dy);
                gradMag1 = edgeXYNormal[(i - 1) * width + j] * (1 - d) + edgeXYNormal[(i - 1) * width + j + 1] * d;
                gradMag2 = edgeXYNormal[(i + 1) * width + j] * (1 - d) + edgeXYNormal[(i + 1) * width + j - 1] * d;
            }
            else if ((dx <= 0 && dx > dy) || (dx >= 0 && dx < dy))
            {
                // 90 ~ 135°
                double d = abs(dx / dy);
                gradMag1 = edgeXYNormal[(i - 1) * width + j] * (1 - d) + edgeXYNormal[(i - 1) * width + j - 1] * d;
                gradMag2 = edgeXYNormal[(i + 1) * width + j] * (1 - d) + edgeXYNormal[(i + 1) * width + j + 1] * d;
            }
            else if ((dy < 0 && dx <= dy) || (dy > 0 && dx >= dy))
            {
                // 135 ~ 180°
                double d = abs(dy / dx);
                gradMag1 = edgeXYNormal[i * width + j - 1] * (1 - d) + edgeXYNormal[(i - 1) * width + j - 1] * d;
                gradMag2 = edgeXYNormal[i * width + j + 1] * (1 - d) + edgeXYNormal[(i + 1) * width + j + 1] * d;
            }
            if (gradMag >= gradMag1 && gradMag >= gradMag2)
            {
                if (gradMag >= highThreshold)
                {
                    nonMaximum[i * width + j] = (uchar)255;
                }
                else if (gradMag >= lowThreshold)
                {
                    nonMaximum[i * width + j] = (uchar)127;
                }
                else
                {
                    nonMaximum[i * width + j] = (uchar)0;
                }
            }
            else
            {
                // 非极大抑制
                nonMaximum[i * width + j] = (uchar)0;
            }
        }
    }
}

void EdgeConnectBwlabel(uchar *nonMaximum, uchar *conEd, int width, int height)
{
    memset(conEd, 0, height * width);
    for (int i = 0; i < height * width; i++) {
        if (nonMaximum[i] > 0) {
            conEd[i] = uchar(255);
        }
    }
    
    int numRegion = 0;
    int *label = new int[height * width];
    memset(label, 0, sizeof(int) * height * width);
    Mor_Bwlabel(conEd, label, height, width, 8, numRegion);
    
    
    bool *flag = new bool[numRegion];
    for (int i = 0; i < numRegion; i++) {
        flag[i] = true;
    }
    for (int i = 0; i < height * width; i++) {
        if (label[i] == 0) {
            continue;
        }
        flag[label[i] - 1] &= nonMaximum[i] == 127;
    }
    for (int i = 0; i < height * width; i++) {
        if (label[i] == 0) {
            continue;
        }
        if (flag[label[i] - 1]) {
            conEd[i] = uchar(0);
        }
    }
    
    
    delete [] flag;
    delete [] label;
}


/// canny 边缘检测
void EdgeDetect_Canny(const uchar *data, uchar *bw, double *edgeX, double *edgeY, int height, int width)
{
    assert(data != NULL);
    // 高斯平滑
    double sigma = 1.0;
    uchar *fdata = new uchar[width * height];
    GaussFilter(data, fdata, width, height, sigma);
    
    // 求梯度
    double *edgeXY = new double[width * height];
    double *theta = new double[width * height];
    EdgeDetect_Gauss(fdata, edgeX, edgeY, edgeXY, theta, width, height, sigma);
    
    // 自适应计算高低阈值 
    double lowThreshold = 0.0;
    double highThreshold = 1.0;
    double * edgeXYNormal = new double[width * height];
    AutoCalHighLowThreshold(edgeXY, edgeXYNormal, width, height, highThreshold, lowThreshold);
    
    // 非极大抑制
    uchar *nonMaximum = new uchar[width * height];
    NonMaximumInhibiion(edgeX, edgeY, edgeXYNormal, nonMaximum, width, height, lowThreshold, highThreshold);
   
    
    
    
    // 连接边缘
    memset(bw, 0, height * width);
    EdgeConnectBwlabel(nonMaximum, bw, width, height);
    
    delete[] nonMaximum;
    delete[] edgeXYNormal;
    delete[] edgeX;
    delete[] edgeY;
    delete[] edgeXY;
    delete[] theta;
    delete[] fdata;
}
void Qsort(double a[], int low, int high)
{
    if (low >= high)
    {
        return;
    }
    int first = low;
    int last = high;
    double key = a[first];/*用字表的第一个记录作为枢轴*/
    while (first < last)
    {
        while (first < last && a[last] >= key)
            last--;
        a[first] = a[last];/*将比第一个小的移到低端*/
        while (first < last && a[first] <= key)
            first++;
        a[last] = a[first];/*将比第一个大的移到高端*/
    }
    a[first] = key;/*枢轴记录到位*/
    Qsort(a, low, first - 1);
    Qsort(a, first + 1, high);
}


void EdgeThinning( WORD *edgeMag, short *edgeX, short *edgeY,
                    uchar *edgeThin, int width, int height, float threshLow, float threshHigh )
{
    assert( edgeMag && edgeX && edgeY && edgeThin );
    assert( width % 4 == 0 );
    assert( width > 0 && height > 0 );
    const float TAN_225 = 0.41421356f;      //==tan(22.5)
    const float TAN_675 = 2.41421356f;      //==tan(67.5)
    int r, c, p;
    WORD *edgeTmp = new WORD[width*height];
    memcpy( edgeTmp, edgeMag, width*height*sizeof(WORD) );
    for( p=0; p<width*height; p++ )
    {
        if( edgeTmp[p] >= threshLow && edgeTmp[p] <= threshHigh )
        {
            edgeThin[p] = 255;
        }
        else
        {
            edgeThin[p] = 0;
        }
    }
    //non-maxima supression along the normal direction
    for( r=1; r<height-1; r++ )
        for( c=1; c<width-1; c++ )
        {
            p = r * width + c;
            if( edgeThin[p] == 0 )
            {
                continue;
            }
            float tmp;
            if( edgeX[p] ) tmp = -((float)edgeY[p]) / ((float)edgeX[p]);
            else tmp = 10000;   //a large value
            if( tmp > -TAN_225 && tmp < TAN_225 )
            {
                if( edgeTmp[p] <= max( edgeTmp[p-1], edgeTmp[p+1] ) )
                {
                    edgeTmp[p] -= 1;
                    edgeThin[p] = 0;
                }
            }
            else if( tmp > TAN_225 && tmp < TAN_675 )
            {
                if( edgeTmp[p] <= max( edgeTmp[p-width+1], edgeTmp[p+width-1] ) )
                {
                    edgeTmp[p] -= 1;
                    edgeThin[p] = 0;
                }
            }
            else if( tmp > TAN_675 || tmp < -TAN_675 )
            {
                if( edgeTmp[p] <= max( edgeTmp[p-width], edgeTmp[p+width] ) )
                {
                    edgeTmp[p] -= 1;
                    edgeThin[p] = 0;
                }
            }
            else
            {
                if( edgeTmp[p] <= max( edgeTmp[p-width-1], edgeTmp[p+width+1] ) )
                {
                    edgeTmp[p] -= 1;
                    edgeThin[p] = 0;
                }
            }
        }
    delete []edgeTmp;
}
void EdgeDetectAEX( const uchar *input, short *edgeX, int width, int height )
{
    assert( input && edgeX );
    assert( width % 4 == 0 );
    assert( width > 0 && height > 0 );
    const short maskX[25] = { 0, -6, 0, 6, 0, -12,
        -20, 0, 20, 12, -22, -20, 0, 20, 22, -12, -20, 0, 20, 12, 0, -6, 0, 6, 0 };
    memset( edgeX, 0, width*height*sizeof(short) );
    int r, c, p;
    for( r=2; r<height-2; r++ )
        for( c=2; c<width-2; c++ ) {
            p = r*width + c;
            edgeX[p] =
            ( input[p-2*width-2] * maskX[0] + input[p-2*width-1] * maskX[1] +
             input[p-2*width] * maskX[2] + input[p-2*width+1] * maskX[3] + input[p-2*width+2] * maskX[4] +
             input[p-width-2] * maskX[5] + input[p-width-1] * maskX[6] +
             input[p-width] * maskX[7] + input[p-width+1] * maskX[8] + input[p-width+2] * maskX[9] +
             input[p-2] * maskX[10] + input[p-1] * maskX[11] +
             input[p] * maskX[12] + input[p+1] * maskX[13] + input[p+2] * maskX[14] +
             input[p+width-2] * maskX[15] + input[p+width-1] * maskX[16] +
             input[p+width] * maskX[17] + input[p+width+1] * maskX[18] + input[p+width+2] * maskX[19] +
             input[p+2*width-2] * maskX[20] + input[p+2*width-1] * maskX[21] +
             input[p+2*width] * maskX[22] + input[p+2*width+1] * maskX[23] + input[p+2*width+2] * maskX[24] ) / 2;
        }
}
void EdgeDetectAEY( const uchar *input, short *edgeY, int width, int height )
{
    assert( input && edgeY );
    assert( width % 4 == 0 );
    assert( width > 0 && height > 0 );
    const short maskY[25] = { 0, -12, -22, -12, 0,
        -6, -20, -20, -20, -6, 0, 0, 0, 0, 0, 6, 20, 20, 20, 6, 0, 12, 22, 12, 0 };
    memset( edgeY, 0, width*height*sizeof(short) );
    int r, c, p;
    for( r=2; r<height-2; r++ )
        for( c=2; c<width-2; c++ ) {
            p = r*width + c;
            edgeY[p] =
            ( input[p-2*width-2] * maskY[0] + input[p-2*width-1] * maskY[1] +
             input[p-2*width] * maskY[2] + input[p-2*width+1] * maskY[3] + input[p-2*width+2] * maskY[4] +
             input[p-width-2] * maskY[5] + input[p-width-1] * maskY[6] +
             input[p-width] * maskY[7] + input[p-width+1] * maskY[8] + input[p-width+2] * maskY[9] +
             input[p-2] * maskY[10] + input[p-1] * maskY[11] +
             input[p] * maskY[12] + input[p+1] * maskY[13] + input[p+2] * maskY[14] +
             input[p+width-2] * maskY[15] + input[p+width-1] * maskY[16] +
             input[p+width] * maskY[17] + input[p+width+1] * maskY[18] + input[p+width+2] * maskY[19] +
             input[p+2*width-2] * maskY[20] + input[p+2*width-1] * maskY[21] +
             input[p+2*width] * maskY[22] + input[p+2*width+1] * maskY[23] + input[p+2*width+2] * maskY[24] ) / 2;
        }
}
float calc_threshMag(WORD *edgeMag, int width, int height, double threshRatio)
{
    const int MIN_MAG = 256; //for AE, AC detector
    int number  = 0;
    double mean = 0;
    int max  = 0;
    //double sum = 0;
    for( int r = 0; r < height; r++ )
    {
        for( int c = 0; c < width; c++ )
        {
            WORD mag = edgeMag[r * width + c];
            if( mag > MIN_MAG )
            {
                //sum += mag;
                number ++;
                if (mag > max)
                {
                    max = mag;
                }
            }
        }
    }
    if( number > 0 )
    {
        if (threshRatio >= 1)
            return (float)max;
        //mean = sum / number;//考虑摄取强度最大的一些点后计算
        // Get Histogram
        const static int numHist = 1000;
        int HistBin[numHist];
        memset(HistBin, 0, numHist * sizeof(int));
        for( int r = 0; r < height; r++ )
        {
            for( int c = 0; c < width; c++ )
            {
                WORD mag = edgeMag[r * width + c];
                if (mag > MIN_MAG)
                {
                    HistBin[(mag - MIN_MAG) * numHist / (max - MIN_MAG + 1)]++;
                }
            }
        }
        // Statistically find the max magnitude threshold: x which satisfy numOf(edgeMag > x) / totalNum >= threshRatio (percentage)
        int sum = 0;
        int minAllowNum = (int)((1 - threshRatio) * number);
        for (int i = numHist - 1; i >= 0; i--)
        {
            sum += HistBin[i];
            if (sum >= minAllowNum)
            {
                mean = i * (max - MIN_MAG + 1) / numHist + MIN_MAG;
                break;
            }
        }
    }
    else
    {
        mean = MIN_MAG;
    }
    //return (float)MIN_MAG;
    return (float)mean;
}

void EdgeDetect(uchar *data, uchar *edge, short *edgeX, short *edgeY, int height, int width)
{
    EdgeDetectAEX(data, edgeX, width, height);
    EdgeDetectAEY(data, edgeY, width, height);
    
    WORD *edgeMag = new WORD[height * width];
    
    for (int i = 0; i < height * width; i++) {
        edgeMag[i] = sqrt(edgeX[i] * edgeX[i] + edgeY[i] * edgeY[i]);
    }
    double thresholdLow = calc_threshMag(edgeMag, width, height, 0);
    double thresholdHigh = calc_threshMag(edgeMag, width, height, 1);
    EdgeThinning(edgeMag, edgeX, edgeY, edge, width, height, thresholdLow, thresholdHigh);
    
    delete [] edgeMag;
    
}


