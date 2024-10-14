#include <thread>
#include <queue>
#include <functional>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <pcl/common/common.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/kdtree/kdtree.h>

#include <DGtal/math/Statistic.h>
#include <DGtal/math/BasicMathFunctions.h>

#include "CenterlineMultiCore.h"
#include "BSplines.h"
#include "Regression.h"

#include "../Common/Rosin.h"


using namespace DGtal;


std::vector<Z3i::RealPoint>
CenterlineMultiCore::compute(){
    Z3i::RealPoint upPt = points[0];
    Z3i::RealPoint lowPt = points[0];
    for(Z3i::RealPoint p: points){
        lowPt = lowPt.inf(p);
        upPt = upPt.sup(p);
    }

    //get main centerline
    std::vector<Z3i::RealPoint> mainLine = getRawLine(lowPt, upPt);

    
    //IOHelper::export2Text(mainLine, "raw.xyz");
    //std::cout<<"found " << mainLine.size()<<" points "<<std::endl;

    std::vector<Z3i::RealPoint> smoothLine;
    if(mainLine.size() < nbControlPoint + 2){
        std::cout<<"too few points for bsplines: " << mainLine.size()<<" points "<<std::endl;
        return smoothLine;
    }

    PointOrderByZ pointOrderByZ;
    //Spline with 8 break
    std::sort(mainLine.begin(), mainLine.end(), pointOrderByZ);

    int nbKnots = nbControlPoint + 2;

    //IOHelper::export2Text(mainLine, "rawsort.xyz");

    smoothLine = BSplines::bsplines(mainLine, nbKnots);
    std::sort(smoothLine.begin(), smoothLine.end(), pointOrderByZ);
    
    Z3i::RealPoint firstVect = (mainLine[0] - mainLine[10]);
    firstVect=firstVect/firstVect.norm();
    
    Z3i::RealPoint firstPoint = smoothLine[0];
    
    Z3i::RealPoint lastVect = (mainLine[mainLine.size()-1] - mainLine[mainLine.size()-10]);
    lastVect=lastVect/lastVect.norm();
    Z3i::RealPoint lastPoint = smoothLine[smoothLine.size()-1];
   
    std::vector<Z3i::RealPoint> frontVect;

    Z3i::RealPoint nextFront = firstPoint + (firstVect*0.5);

    //extend to the 2 ends of trunk
    //mm -> m
    double zLow = lowPt[2] - 0.1;
    double zUp = upPt[2] + 0.1;


    while(nextFront[2] > zLow && nextFront[2] < zUp){
        frontVect.push_back(nextFront);
        nextFront = nextFront + (firstVect*0.5);
    }

    std::vector<Z3i::RealPoint> backVect;
    Z3i::RealPoint nextBack = lastPoint + (lastVect*0.5);
    while(nextBack[2] > zLow && nextBack[2] < zUp){
        backVect.push_back(nextBack);
        nextBack = nextBack + (lastVect*0.5);
    }

    std::reverse(frontVect.begin(), frontVect.end());

    frontVect.insert(frontVect.end(), smoothLine.begin(), smoothLine.end());
    frontVect.insert(frontVect.end(), backVect.begin(), backVect.end());
    smoothLine = frontVect;
    //IOHelper::export2Text(smoothLine, "smoothLinecompleted.xyz");

    return smoothLine;

}

std::vector<Z3i::RealPoint>
CenterlineMultiCore::getRawLine(const Z3i::RealPoint &lowPt, const Z3i::RealPoint &upPt){

    //if voxel size is 3 -> 900m
    //double segLengthInitial = 300.0;
    //separe into multiple parts
    //zLength =

    double zLength = upPt[2] -lowPt[2];
    //int nbSegment = zLength / segLengthInitial;

    double zMin = lowPt[2];


    //adjusted
    double segLength = zLength / nbSegment;
    double segLengthWithBuff = segLength*1.1;
    double segOverlapLength = segLength*0.1;


    std::vector<double> beginSegments(nbSegment);
    std::vector<double> endSegments(nbSegment);

    for(int i = 0; i < nbSegment; i++){
        beginSegments[i] = zMin + segLength * i;
        endSegments[i] = zMin + segLength * (i + 1) + segOverlapLength;
        //endSegments[i] = zMin + segLengthWithBuff * (i + 1);
    }

    std::vector< std::vector<Z3i::RealPoint> > segPoints(nbSegment);
    //for presentation
    std::vector< std::vector<int> > segPids(nbSegment);
    std::vector< std::vector<Z3i::RealPoint> > segNormals(nbSegment);
    Z3i::RealPoint zero(0,0,0);
    float dMax = std::numeric_limits<float>::max();
    std::vector<Z3i::RealPoint> segLowPts(nbSegment, Z3i::RealPoint(dMax, dMax, dMax));
    std::vector<Z3i::RealPoint> segUpPts(nbSegment, Z3i::RealPoint(-dMax, -dMax, -dMax));
    std::vector<Z3i::Domain> segDomains(nbSegment);

    for(int pId = 0; pId < points.size(); pId++){
        Z3i::RealPoint p = points[pId];
        Z3i::RealPoint n = normals[pId];
        for(int segId = 0; segId < nbSegment; segId++){
            if( p[2] > beginSegments[segId] && p[2] < endSegments[segId]) {
                //for presentation
                segPids[segId].push_back(pId);
                segPoints[segId].push_back(p);
                segNormals[segId].push_back(n);
                segLowPts[segId] = segLowPts[segId].inf(p);
                segUpPts[segId] = segUpPts[segId].sup(p);
            }
        }
    }

    std::vector<Z3i::RealPoint>  rawLine;
    

    std::vector<std::thread> ts(nbSegment);
    //std::vector<std::future<std::vector<Z3i::RealPoint> > > futures(nbSegment);
    std::vector<std::vector<Z3i::RealPoint> > raws(nbSegment);
    for(int segId = 0; segId < nbSegment; segId++){
        if(segPoints[segId].size() == 0){
            continue;
        }
        //trace.info()<<"segLowPts i : " << segLowPts[segId]- DGtal::Z3i::RealPoint::diagonal(1)<< "segUpPts i : "<< segUpPts[segId]+ DGtal::Z3i::RealPoint::diagonal(1) <<std::endl;

        Z3i::Domain sDomain = Z3i::Domain(DGtal::PointVector<3, int>(segLowPts[segId] - DGtal::Z3i::RealPoint::diagonal(1)), DGtal::PointVector<3, int>(segUpPts[segId] + DGtal::Z3i::RealPoint::diagonal(1)));
        //trace.info()<<accRadius << "domaine : "<< sDomain << " "<< segNormals[segId].size()<<"  "<<segPoints[segId].size()<< "  "<< nbControlPoint<<std::endl;
        //threshold is deprecated

        Centerline cen(accRadius, sDomain, segNormals[segId], segPoints[segId], threshold, nbControlPoint);

        //ts[segId] = std::thread(&Centerline::computeRawLine, &cen, std::ref(raws[segId]));
        //std::vector<Z3i::RealPoint> raw;
        //ts[segId] = std::thread(&CenterlineMultiCore::computeRaw, this, accRadius, sDomain, segNormals[segId], segPoints[segId], nbControlPoint, raws[segId]);
        //ts[segId] = std::thread([&cen] { computeRawLine(&raws[segId]); });
        cen.computeRawLine(raws[segId]);

        //for presentation
        //std::string segPointName = "segPoints" + std::to_string(segId) + ".xyz";
        //IOHelper::export2Text(segPoints[segId], segPointName);

        //std::string segPiDName = "segPid" + std::to_string(segId) + ".id";
        //IOHelper::export2Text(segPids[segId], segPiDName);

        //std::string rawLineName = "rawline" + std::to_string(segId) + ".xyz";
        //IOHelper::export2Text(raws[segId], rawLineName);
    }

    //for(int i = 0; i < nbSegment; i++){
    //    ts[i].join();
    //}

    //std::cout<<"nbSegment " << nbSegment<<std::endl;
    for(int i = 0; i < nbSegment; i++){

        //std::cout<<"found " << raws[i].size()<<std::endl;
        rawLine.insert(rawLine.end(), raws[i].begin(), raws[i].end());
    }
    //std::cout<<"found " << rawLine.size()<<" points "<<std::endl;
    return rawLine;
}
