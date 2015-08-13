//
//  ImageMor.cpp
//  Kernel
//
//  Created by JiaYonglei on 15/8/11.
//  Copyright (c) 2015年 JiaYonglei. All rights reserved.
//
#include "Image.h"

int RunlabelDeleteError(int *runlabel, int *outlabel, int numRuns, vector<vector<int> > sameLabel)
{
    for (int i = 0; i < sameLabel.size(); i++)
    {
        
        if (sameLabel[i].size() < 2)
        {
            continue;
        }
        for (int j = 0; j < sameLabel[i].size(); j++)
        {
            for (int k = 0; k < numRuns; k++)
            {
                if (runlabel[k] == sameLabel[i][j])
                {
                    runlabel[k] = sameLabel[i][0];
                }
            }
        }
    }
    // label 去重，排序
    vector<int > temp;
    temp.push_back(runlabel[0]);
    for (int i = 1; i < numRuns; i++)
    {
        bool flag = false;
        for (int j = 0; j < temp.size(); j++)
        {
            flag |= runlabel[i] == temp[j];
        }
        if (!flag)
        {
            temp.push_back(runlabel[i]);
        }
    }
    int *data = new int[temp.size()];
    int *index = new int[temp.size()];
    for (int i = 0; i < temp.size(); i++)
    {
        data[i] = temp[i];
        index[i] = i;
    }
    for (int i = 0; i < temp.size() - 1; i++)
    {
        for (int j = i + 1; j < temp.size(); j++)
        {
            if (data[i] > data[j])
            {
                int temp =data[i];
                data[i] = data[j];
                data[j] = temp;
                temp = index[i];
                index[i] = index[j];
                index[j] = temp;
            }
        }
    }
    for (int i = 0; i < numRuns; i++)
    {
        int index = 0;
        for (int j = 0; j < temp.size(); j++)
        {
            if (runlabel[i] == data[j])
            {
                index = j;
                break;
            }
        }
        outlabel[i] = index + 1;
    }
    
    delete[] data;
    delete[] index;
    
    return int(temp.size());
}


void CalNumberOfRuns(const uchar *data, int height, int width, int &numRuns)
{
    assert(data != NULL);
    int count = 0;
    for (int j = 0; j < width; j++)
    {
        if (data[j] != 0)
        {
            count++;
        }
        for (int i = 1; i < height; i++)
        {
            if (data[i * width + j] != 0 && data[(i - 1) * width + j] == 0)
            {
                count++;
            }
        }
    }
    numRuns = count;
}
// 记录每一个runs的开始行、结束行和runs所在的列
void FillRunVectors(const uchar* data, int height, int width, int *firstRowDex, int *lastRowDex, int *colDex, int numRuns)
{
    assert(data != NULL  && firstRowDex != NULL && lastRowDex != NULL && colDex != NULL);
    int k = -1;
    for (int j = 0; j < width; j++)
    {
        if (data[j] != 0)
        {
            k++;
            firstRowDex[k] = 0;
            colDex[k] = j;
        }
        for (int i = 1; i < height; i++)
        {
            int index = i * width + j;
            int _index = index - width;
            if (data[index] == 0 && data[_index] != 0)
            {
                lastRowDex[k] = i - 1;
            }
            else
            {
                if (data[index] != 0)
                {
                    if (data[_index] == 0)
                    {
                        k++;
                        firstRowDex[k] = i;
                        colDex[k] = j;
                    }
                    
                    if (i == height - 1)
                    {
                        lastRowDex[k] = i;
                    }
                }
            }
        }
    }
    assert(k == numRuns - 1);
}
/// 标记二值连通区域，返回连通区域个数
// 速度太慢
bool Mor_BwlabelSlow(const uchar *data, int *label, int height, int width, int mode, int &numRegions)
{
    assert(data != NULL && label != NULL);
//    int totalPixels = height * width;
    // calculate the number of runs in data
    int numRuns = 0;
    CalNumberOfRuns(data, height, width, numRuns);
    // calculate run vectors, including firstRowDex, lastRowDex, colDex of each run
    int *firstRowDex = new int [numRuns];
    int *lastRowDex = new int [numRuns];
    int *colDex = new int [numRuns];
    FillRunVectors(data, height, width, firstRowDex, lastRowDex, colDex, numRuns);
    
    // Initiate runLabels (the label of each run)
    int *runLabels = new int[numRuns];
    for (int i = 0; i < numRuns; i++)
    {
        runLabels[i] = 0;
    }
    
    // Initiate params
    int currentColumn = -1;
    int nextLabel = 1;
    int firstRunOnPreviousColumn = -1;
    int lastRunOnPreviousColumn = -1;
    int firstRunOnThisColumn = -1;
    int offset = 0;
    if (mode == 8)
    {
        offset = 1;
    }
    
    //////////////////////////////////////////////////////////////////////////
    // 利用图深度优先搜索合并错误标签
    vector<edge> sameLabel;
    int p;
    
    for (int k = 0; k < numRuns; k++)
    {
        if (colDex[k] == currentColumn + 1)
        {
            // 第k个runs跟第k-1个runs在相邻列上
            firstRunOnPreviousColumn = firstRunOnThisColumn;
            firstRunOnThisColumn = k;
            lastRunOnPreviousColumn = k-1;
            currentColumn = colDex[k];
        }
        else if (colDex[k] > currentColumn + 1)
        {
            firstRunOnPreviousColumn = -1;
            firstRunOnThisColumn = k;
            lastRunOnPreviousColumn = -1;
            currentColumn = colDex[k];
        }
        
        if (firstRunOnPreviousColumn >= 0)
        {
            p = firstRunOnPreviousColumn;
            while (p <= lastRunOnPreviousColumn && lastRowDex[k] >= firstRowDex[p] - offset)
            {
                if (firstRowDex[k] <= lastRowDex[p] + offset)
                {
                    if (runLabels[k] == 0)
                    {
                        runLabels[k] = runLabels[p];
                    }
                    else if (runLabels[k] != runLabels[p])
                    {
                        edge etemp;
                        etemp.start = runLabels[k];
                        etemp.end = runLabels[p];
                        sameLabel.push_back(etemp);
                    }
                }
                p++;
            }
        }
        
        if (runLabels[k] == 0)
        {
            runLabels[k] = nextLabel;
            nextLabel++;
        }
    }
    
    vector<vector<int> > tempLabel;
    DfnComponent(sameLabel, tempLabel);
    
    int *newlabel = new int[numRuns];
    numRegions = RunlabelDeleteError(runLabels, newlabel, numRuns, tempLabel);
    
    // 赋值label
    memset(label, 0, sizeof(int) * width * height);
    for (int i = 0; i < numRuns; i++)
    {
        for (int k = firstRowDex[i]; k <= lastRowDex[i]; k++)
        {
            label[k * width + colDex[i]] = newlabel[i];
        }
    }
    
    delete[] newlabel;
    delete[] firstRowDex;
    delete[] lastRowDex;
    delete[] colDex;
    delete[] runLabels;
    return true;
}
void Mor_BwClearBorder(uchar *bw, int *label, int height, int width, int *edgeLabel, int num)
{
    memset(edgeLabel, 0, (num + 1) * sizeof(int));
    for (int i = 0; i < height; i++) {
        if (label[i * width + 0] != 0) {
            edgeLabel[label[i * width + 0]] = 1;
        }
        if (label[i * width + width - 1] != 0) {
            edgeLabel[label[i * width + width - 1]] = 1;
        }
    }
    for (int i = 0; i < width; i++) {
        if (label[i] != 0) {
            edgeLabel[label[i]] = 1;
        }
        if (label[(height - 1) * width + i] != 0) {
            edgeLabel[label[(height - 1) * width + i]] = 1;
        }
    }
    for (int i = 0; i < height; i ++) {
        for (int j = 0; j < width; j++) {
            int LL = label[i * width + j];
            if (LL == 0) {
                continue;
            }
            if (edgeLabel[LL] == 1) {
                bw[i * width + j] = 0;
            }
        }
    }
}
void Mor_BwDilate(uchar *data, int height, int width, int winsize)
{
    assert(winsize % 2 != 0);
    uchar *out = new uchar[height * width];
    memcpy(out, data, height * width);
  
    
    int halfLen = winsize / 2;
    for (int i = halfLen; i < height - halfLen; i++) {
        for (int j = halfLen; j < width - halfLen; j++) {
            if (data[i * width + j] == 255) {
                continue;
            }
            bool flag = true;
            for (int m = i - halfLen; m <= i + halfLen; m++) {
                for (int n = j - halfLen; n <= j + halfLen; n++) {
                    flag &= data[m * width + n] == 0;
                }
            }
            if (flag) {
                out[i * width + j] = uchar(0);
            }
            else{
                out[i * width + j] = uchar(255);
            }
        }
    }
    memcpy(data, out, height * width);
    delete [] out;
}
void Mor_BwErode(uchar *data, int height, int width, int winsize)
{
    assert(winsize % 2 != 0);
    uchar *out = new uchar[height * width];
    memcpy(out, data, height * width);
    
    int halfLen = winsize / 2;
    for (int i = halfLen; i < height - halfLen; i ++) {
        for (int j = halfLen; j < width - halfLen; j++) {
            if (data[i * width + j] == 0) {
                continue;
            }
            bool flag = true;
            for (int m = i - halfLen; m <= i + halfLen ; i++) {
                for (int n = j - halfLen; n <= j + halfLen; j++) {
                    flag |= data[m * width + n] == 255;
                }
            }
            if (flag) {
                out[i * width + j] = uchar(255);
            }
            else{
                out[i * width + j] = uchar(0);
            }
        }
    }
    memcpy(data, out, height * width);
    delete [] out;
   
}
void Mor_BwRegionProps(int *label, int height, int width, int num, int *area, gclRect *boundingBox)
{
    memset(area, 0, sizeof(int) * num);
    vector< vector<gclPoint> > pt(num, vector<gclPoint> (0));
    for (int i = 0; i < height; i ++) {
        for (int j = 0; j < width; j++) {
            int index = label[i * width + j];
            if (index == 0) {
                continue;
            }
            area[index - 1] ++;
            gclPoint ptemp;
            ptemp.x = j;
            ptemp.y = i;
            pt[index - 1].push_back(ptemp);
        }
    }
    // bounding box
    for (int i = 0;  i < num; i++) {
        double minX = DBL_MAX;
        double maxX = -DBL_MAX;
        double minY = DBL_MAX;
        double maxY = -DBL_MAX;
        for (int j = 0; j < pt[i].size(); j++) {
            if (pt[i][j].x > maxX) {
                maxX = pt[i][j].x;
            }
            if (pt[i][j].x < minX) {
                minX = pt[i][j].x;
            }
            if (pt[i][j].y > maxY) {
                maxY = pt[i][j].y;
            }
            if (pt[i][j].y < minY) {
                minY = pt[i][j].y;
            }
        }
        boundingBox[i].x = minX;
        boundingBox[i].y = minY;
        boundingBox[i].height = maxY - minY + 1;
        boundingBox[i].width = maxX - minX + 1;
    }
}

// 暂时由opencv实现imfill函数
void Mor_BwImfillOpencv(uchar *data, int height, int width)
{
    //step 1: make a border
    Mat src(height, width, CV_8UC1, data);
    std::vector< std::vector< cv::Point>> contours;
    std::vector< Vec4i > hierarchy;
    findContours(src, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_NONE);
    
    
    if( !contours.empty() && !hierarchy.empty() )
    {
        for (int idx=0;idx < contours.size();idx++)
        {
            drawContours(src,contours,idx,Scalar::all(255),CV_FILLED,8);
        }
    }
    
    
    memcpy(data, src.data, height * width);
}

void Mor_BwDilateOpencv(uchar *data, int height, int width, int size)
{
    
    Mat element = getStructuringElement(MORPH_RECT, Size(size, size));
    Mat src(height, width, CV_8UC1, data);
    Mat dst;
    dilate(src, dst, element);
    memcpy(data, dst.data, height * width);
}

void Mor_BwErodeOpencv(uchar *data, int height, int width, int size)
{
    Mat element = getStructuringElement(MORPH_RECT, Size(size, size));
    Mat src(height, width, CV_8UC1, data);
    Mat dst;
    erode(src, dst, element);
    memcpy(data, dst.data, height * width);
}
void CalNumberOfRuns1(const uchar *data, int height, int width, int &numRuns)
{
    assert(data != NULL);
    int count = 0;
    for (int j = 0; j < width; j++)
    {
        if (data[j] != 0)
        {
            count++;
        }
        for (int i = 1; i < height; i++)
        {
            if (data[i * width + j] != 0 && data[(i - 1) * width + j] == 0)
            {
                count++;
            }
        }
    }
    numRuns = count;
}

void FillRunVectors1(const uchar* data, int *runs, int height, int width, int *firstRowDex, int *lastRowDex, int *colDex, int numRuns)

{
    assert(data != NULL && runs != NULL && firstRowDex != NULL && lastRowDex != NULL && colDex != NULL);
    int k = -1;
    memset(runs, 0, height * width * sizeof(int));
    for (int j = 0; j < width; j++)
    {
        if (data[j] != 0)
        {
            k++;
            firstRowDex[k] = 0;
            colDex[k] = j;
            runs[j] = k + 1;
        }
        for (int i = 1; i < height; i++)
        {
            int index = i * width + j;
            int _index = index - width;
            if (data[index] == 0 && data[_index] != 0)
            {
                lastRowDex[k] = i - 1;
            }
            else
            {
                if (data[index] != 0)
                {
                    if (data[_index] == 0)
                    {
                        k++;
                        firstRowDex[k] = i;
                        colDex[k] = j;
                    }
                    runs[index] = k + 1;
                    if (i == height - 1)
                    {
                        lastRowDex[k] = i;
                    }
                }
            }
        }
    }
    assert(k == numRuns - 1);
}

// 标记二值连通区域，返回连通区域个数

bool Mor_Bwlabel(const uchar *data, int *label, int height, int width, int mode, int &numRegions)
{
    int minarea = 1;
    assert(data != NULL && label != NULL);
    int totalPixels = height * width;
    // calculate the number of runs in data
    int numRuns = 0;
    CalNumberOfRuns1(data, height, width, numRuns);
    // calculate run vectors, including firstRowDex, lastRowDex, colDex of each run
    int *firstRowDex = new int [numRuns];
    int *lastRowDex = new int [numRuns];
    int *colDex = new int [numRuns];
    int *runs = new int[totalPixels];
    FillRunVectors1(data, runs, height, width, firstRowDex, lastRowDex, colDex, numRuns);
    // Initiate runLabels (the label of each run)
    int *runLabels = new int[numRuns];
    for (int i = 0; i < numRuns; i++)
    {
        runLabels[i] = 0;
    }
    // 初始化等价集合表
    int *equalSet = new int[numRuns];
    for (int i = 0; i < numRuns; i++)
    {
        equalSet[i] = i;
    }
    // Initiate params
    int currentColumn = -1;
    int nextLabel = 1;
    int firstRunOnPreviousColumn = -1;
    int lastRunOnPreviousColumn = -1;
    int firstRunOnThisColumn = -1;
    int offset = 0;
    if (mode == 8)
    {
        offset = 1;
    }
    for (int k = 0; k < numRuns; k++)
    {
        if (colDex[k] == currentColumn + 1)
        {
            firstRunOnPreviousColumn = firstRunOnThisColumn;
            firstRunOnThisColumn = k;
            lastRunOnPreviousColumn = k-1;
            currentColumn = colDex[k];
        }
        else if (colDex[k] > currentColumn + 1)
        {
            firstRunOnPreviousColumn = -1;
            firstRunOnThisColumn = k;
            lastRunOnPreviousColumn = -1;
            currentColumn = colDex[k];
        }
        if (firstRunOnPreviousColumn >= 0)
        {
            int p = firstRunOnPreviousColumn;
            while (p <= lastRunOnPreviousColumn && lastRowDex[k] >= firstRowDex[p] - offset)
            {
                if (firstRowDex[k] <= lastRowDex[p] + offset)
                {
                    if (runLabels[k] == 0)
                    {
                        runLabels[k] = runLabels[p];
                    }
                    else if (runLabels[k] != runLabels[p])
                    {
                        int temp1 = equalSet[runLabels[p]];
                        int temp2 = equalSet[runLabels[k]];
                        if (temp2 < temp1)
                        {
                            for (int i=0; i<nextLabel; i++)
                            {
                                if (equalSet[i] == temp1)
                                {
                                    equalSet[i] = temp2;
                                }
                            }
                        }
                        else
                        {
                            for (int i=0; i<nextLabel; i++)
                            {
                                if (equalSet[i] == temp2)
                                {
                                    equalSet[i] = temp1;
                                }
                            }
                        }
                    }
                }
                p++;
            }
        }
        if (runLabels[k] == 0)
        {
            runLabels[k] = nextLabel;
            nextLabel++;
            if (nextLabel > 10e15)
            {
                return false;
            }
        }
    }
    for (int i = 0; i < numRuns; i++)
    {
        runLabels[i] = equalSet[runLabels[i]];
    }
    // 标记连通区域标号
    for (int i = 0; i < totalPixels; i++)
    {
        if (runs[i] != 0)
        {
            label[i] = runLabels[runs[i]-1];
        }
        else
        {
            label[i] = 0;
        }
    }
    int *area = new int[nextLabel];
    memset(area, 0, sizeof(int) * nextLabel);
    for (int i = 0; i < totalPixels; i++)
    {
        area[label[i]]++;
    }
    numRegions = 0;
    for (int i = 1; i < nextLabel; i++)
    {
        if (area[i] < minarea)
        {
            equalSet[i] = 0;
        }
        else
        {
            numRegions++;
            equalSet[i] = numRegions;
        }
    }
    for (int i = 0; i < totalPixels; i++)
    {
        label[i] = equalSet[label[i]];
    }
    delete[] firstRowDex;
    delete[] lastRowDex;
    delete[] colDex;
    delete[] runs;
    delete[] runLabels;
    delete[] equalSet;
    return true;
}






