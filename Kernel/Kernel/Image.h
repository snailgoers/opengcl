//
//  Image.h
//  Kernel
//
//  Created by JiaYonglei on 15/8/11.
//  Copyright (c) 2015年 JiaYonglei. All rights reserved.
//

#ifndef Kernel_Image_h
#define Kernel_Image_h

#include <opencv2/opencv.hpp>
#include <string>
#include <vector>
#include <assert.h>
#include <time.h>

using namespace std;
using namespace cv;




#define EPS 10e-16
#define MinGaussSize 9


typedef unsigned short WORD;


/// 图算法相关
typedef struct EDGE{
    int start, end;
    EDGE()
    {
        start = 0;
        end = 0;
    }
    EDGE(int a, int b)
    {
        start = a;
        end = b;
    }
}edge;
typedef struct GraphicNode{
    int vertex; //顶点数据信息
    GraphicNode *nextNode;
}graphNode;

typedef struct GCLRECT{
    double x, y, width, height;
}gclRect;
typedef struct GCLPOINT{
    double x, y;
}gclPoint;



// image process
void AdaptationThreshold(const uchar *data, uchar *bw, int width, int height, bool Binary_Inv);

void DfnComponent(vector<edge> pt, vector<vector<int> > &outdata);



// image morphology
bool Mor_BwlabelSlow(const uchar *data, int *label, int height, int width, int mode, int &numRegions);
void Mor_BwClearBorder(uchar *bw, int *label, int height, int width, int *edgeLabel, int num);
void Mor_BwDilate(uchar *data, int width, int height, int winsize);
void Mor_BwErode(uchar *data, int height, int width, int winsize);
void Mor_BwRegionProps(int *label, int height, int width, int num, int *area, gclRect *boundingBox);

bool Mor_Bwlabel(const uchar *data, int *label, int height, int width, int mode, int &numRegions);




// image edge detection
void EdgeDetect_Canny(const uchar *data, uchar *bw, double *edgeX, double *edgeY, int height, int width);

void EdgeDetect(uchar *data, uchar *edge, short *edgeX, short *edgeY, int height, int width);



// opencv Function
void ImageFilter2(uchar *data, double *out, int height, int width, double *kernel, int len);
void Mor_BwImfillOpencv(uchar *data, int height, int width);
void Mor_BwDilateOpencv(uchar *data, int height, int width, int size);
void Mor_BwErodeOpencv(uchar *data, int height, int width, int size);










#endif
