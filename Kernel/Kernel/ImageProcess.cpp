//
//  ImageProcess.cpp
//  Kernel
//
//  Created by JiaYonglei on 15/8/11.
//  Copyright (c) 2015年 JiaYonglei. All rights reserved.
//

#include "Image.h"


void ImageFilter2(uchar *data, double *out, int height, int width, double *kernel, int len)
{
    double *sdata = new double[width * height];
    memset(out, 0, sizeof(double) * height * width);
   
    
    for (int i = 0; i < height * width ; i ++) {
        sdata[i] = double(data[i]);
    }
    Mat src(height, width, CV_64FC1, sdata);
    
    Mat dst(height, width, CV_64FC1, out);
    Mat Mkernel(len, len, CV_64FC1, kernel);
    filter2D(src, dst, src.depth(), Mkernel);

    memcpy(out, dst.data, height * width);
   
    delete [] sdata;
}