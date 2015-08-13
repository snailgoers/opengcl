//
//  main.cpp
//  Kernel
//
//  Created by JiaYonglei on 15/8/11.
//  Copyright (c) 2015年 JiaYonglei. All rights reserved.
//

#include <iostream>
#include "Image.h"


void DetectMainArea(uchar *src, uchar *bw, int height, int width, gclRect &rect)
{
    int strelSize = 11;
    
    memset(bw, 0, width * height);
    AdaptationThreshold(src, bw, height, width, true);
    
    
    Mor_BwDilate(bw, height, width, strelSize);
    
    
    int numRegion = 0;
    int *label = new int[height * width];
    Mor_Bwlabel(bw, label, height, width, 8, numRegion);
    
    // calc the region area and boundingbox
    int *area = new int[numRegion];
    gclRect *boundingbox = new gclRect[numRegion];
    Mor_BwRegionProps(label, height, width, numRegion, area, boundingbox);
    
    int *edgeLabel = new int[numRegion + 1];
    Mor_BwClearBorder(bw, label, height, width, edgeLabel, numRegion);
    
    // select the max area pattern
    double maxArea = -DBL_MAX;
    int maxareaLabel  = 0;
    for (int i = 0; i < numRegion; i ++) {
        if (edgeLabel[i + 1] == 1) {
            continue;
        }
        if (area[i] > maxArea) {
            maxArea = area[i];
            maxareaLabel = i + 1;
            rect = boundingbox[i];
        }
    }
    for (int i = 0; i < height * width;  i++) {
        if (label[i] != maxareaLabel) {
            bw[i] = uchar(0);
        }
    }
    Mor_BwImfillOpencv(bw, height, width);
    

    delete [] boundingbox;
    delete [] area;
    delete [] edgeLabel;
    delete [] label;
}

void DetectKeyPosition(uchar *src, int height, int width, int zoomInSize, int zoomOutSize)
{
    // detect mainArea
    gclRect rect;
    uchar *bw = new uchar[height * width];
    DetectMainArea(src, bw, height, width, rect);
    
    int row = int(rect.height);
    int col = int(rect.width);
    uchar *bwCut = new uchar[row * col];
    for (int i = 0; i < row; i ++) {
        for (int j = 0; j < col; j++) {
            bwCut[i * col + j] = bw[(i + int(rect.y)) * width + j + int(rect.x)];
        }
    }
    
    
    // calc the templete models that contain the edge only
    uchar *zoomIn = new uchar[row * col];
    uchar *zoomOut = new uchar[row * col];
    uchar *templete = new uchar[row * col];
    
    
    memcpy(zoomIn, bwCut, row * col);
    memcpy(zoomOut, bwCut, row * col);
    memset(templete, 0, row * col);
    

    Mor_BwDilateOpencv(zoomIn, row, col, zoomInSize);
    Mor_BwErodeOpencv(zoomOut, row, col, zoomOutSize);
    
    for (int i = 0; i < row * col; i++) {
        if (zoomIn[i] != zoomOut[i]) {
            templete[i] = uchar(255);
        }
    }
    
    // calc the initial fore line
    uchar *bwEdge = new uchar[row * col];
    double *edgeX = new double[row * col];
    double *edgeY = new double[row * col];
    
    
    
    printf(" Detect key Position...");
    
    delete [] edgeY;
    delete [] edgeX;
    delete [] bwEdge;
    delete [] templete;
    delete [] zoomOut;
    delete [] zoomIn;
    delete [] bwCut;
    delete [] bw;
}






int main(void) {
    // insert code here...
    
//    string name = "/Users/jiayonglei/Desktop/imag/allied/8.bmp";
    string name = "/Users/jiayonglei/Desktop/imag/a1.bmp";
    Mat src = imread(name, 0);
    
    uchar *data = new uchar[src.rows * src.cols];
    memcpy(data, src.data, src.rows * src.cols);
    
    
    //
    int zoomInSize = 11;    //膨胀图像，以消除边缘外部区域无效像素
    int zoomOutSize = 171;   // 腐蚀图像，以消除边缘内部区域无效像素
   
//    DetectKeyPosition(data, src.rows, src.cols, zoomInSize, zoomOutSize);
    
    
    uchar *bwedge = new uchar[src.cols * src.rows];
    double *ex = new double[src.cols * src.rows];
    double *ey = new double[src.cols * src.rows];
    clock_t start = clock();
    EdgeDetect_Canny(data, bwedge, ex, ey, src.rows, src.cols);
    
    
    
    
    clock_t end = clock();
    printf("System run time is %f\n", double(end - start) / CLOCKS_PER_SEC);
    Mat out(src.rows, src.cols, CV_8UC1, bwedge);
    imwrite("/Users/jiayonglei/Desktop/out.bmp", out);
    
    delete [] data;
    
    return 0;
}
