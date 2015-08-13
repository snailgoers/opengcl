//
//  Graphic.cpp
//  Kernel
//
//  Created by JiaYonglei on 15/8/12.
//  Copyright (c) 2015年 JiaYonglei. All rights reserved.
//
#include "Image.h"

void CreateGpaphic(vector<edge> pt, GraphicNode *G)
{
    for (int i = 0; i < pt.size(); i++)
    {
        int start = pt[i].start;
        int end = pt[i].end;
        // create a new node
        GraphicNode *node, *ptr;
        node = new GraphicNode();
        node->vertex = end;
        node->nextNode = NULL;
        
        ptr = &(G[start]);
        while (ptr->nextNode != NULL)
        {
            ptr = ptr->nextNode;
        }
        ptr->nextNode = node;
    }
}

// 图深度优先搜索算法
void dfn(GraphicNode *G, bool *visited, int current, vector<int> &index)
{
    //printf("vertex = %d\n", current);
    index.push_back(current);
    visited[current] = true;
    GraphicNode *p = &G[current];
    while (p != NULL)
    {
        if (!visited[p->vertex])
        {
            dfn(G, visited, p->vertex, index);
        }
        p = p->nextNode;
    }
}

// 图连通域分量
void DfnComponent(vector<edge> pt, vector<vector<int> > &outdata)
{
    // 统计节点个数, 以空间消耗替代时间消耗
    int nodeNum = max(pt[0].start, pt[0].end);
    int startIndex = min(pt[0].start, pt[0].end);
    for (int i = 1; i < pt.size(); i++)
    {
        if (pt[i].start > nodeNum)
        {
            nodeNum = pt[i].start;
        }
        if (pt[i].end > nodeNum)
        {
            nodeNum = pt[i].end;
        }
        if (pt[i].start < startIndex)
        {
            startIndex = pt[i].start;
        }
        if (pt[i].end < startIndex)
        {
            startIndex = pt[i].end;
        }
    }
    
    // 构造无向图
    vector<edge> pt2; // 无向图应该为双向的，输入pt为单向的，复制反向
    for (int i = 0; i < pt.size(); i++)
    {
        edge ptemp(pt[i].end, pt[i].start);
        pt2.push_back(pt[i]);
        pt2.push_back(ptemp);
    }
    nodeNum = nodeNum + 1; // 开始索引为0
    
    GraphicNode *G = new GraphicNode[nodeNum];
    bool *visited = new bool[nodeNum];
    memset(visited, 0, nodeNum);
    // initial the G
    for (int i = 0; i < nodeNum; i++)
    {
        G[i].vertex = i;
        G[i].nextNode = NULL;
    }
    CreateGpaphic(pt2, G);
    int count = 0;
    for (int i = startIndex; i < nodeNum; i++)
    {
        if (!visited[i] && G[i].nextNode != NULL)
        {
            vector<int> index;
            dfn(G, visited, i, index);
            count++;
            outdata.push_back(index);
        }
    }
    
    
    delete[] G;
    delete[] visited;
}