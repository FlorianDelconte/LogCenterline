#include <iostream>
#include <fstream>
#include <utility>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/PointListReader.h"



#include "DGtal/io/writers/MeshWriter.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/shapes/Mesh.h"

#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"


#include "Centerline/Centerline.h"
#include "Centerline/CenterlineHelper.h"

#include "Common/CylindricalCoordinateSystem.h"


#include "CLI11.hpp"


using namespace DGtal;

typedef typename Mesh<Z3i::RealPoint>::MeshFace Face;

void scaleCloud(std::vector<Z3i::RealPoint> &pts, double scale){
  for(unsigned int i = 0; i <  pts.size(); i++){
    Z3i::RealPoint &pi =  pts[i];
    pi *= scale;
  }
}

void export2Text(const std::vector<DGtal::Z3i::RealPoint> &v, const std::string &filename){
  std::ofstream outStream;
  outStream.open(filename.c_str(), std::ofstream::out);
  for(unsigned int i = 0; i < v.size();i++){
      Z3i::RealPoint p = v[i];
      outStream << p[0] << " "<< p[1] << " "<< p[2]<<std::endl;
  }
  outStream.close();
}

void getHomogenityCloud(const std::vector<Z3i::RealPoint> &pointCloud, const double &resolution,
            const Z3i::RealPoint &ptLow, const Z3i::RealPoint &ptUp, std::map<unsigned int, unsigned int> &indexMap){

  double minx = ptLow[0];
  double miny = ptLow[1];
  double minz = ptLow[2];
  double maxx = ptUp[0];
  double maxy = ptUp[1];
  double maxz = ptUp[2];

  double dimx = ceil((maxx - minx)/resolution);
  double dimy = ceil((maxy - miny)/resolution);
  double dimz = ceil((maxz - minz)/resolution);

  for (unsigned int i = 0; i < pointCloud.size(); i++){
    Z3i::RealPoint point = pointCloud.at(i);
    unsigned int colx = (unsigned int) floor((point[0] - minx) / resolution);
    unsigned int liny = (unsigned int) floor((point[1] - miny) / resolution);
    unsigned int levz = (unsigned int) floor((point[2] - minz) / resolution);
    unsigned int grdIndex = levz * dimx * dimy + liny * dimx + colx;
    if( indexMap.find(grdIndex) != indexMap.end()){
      indexMap[grdIndex] = i;
    }else{
      unsigned int previousPointIndex  = indexMap[grdIndex];

      double gridx = minx + colx*resolution + resolution/2;
      double gridy = miny + liny*resolution + resolution/2;
      double gridz = minz + levz*resolution + resolution/2;

      double distance = pow(point[0] - gridx, 2) + pow(point[1] - gridy, 2) + pow(point[2] - gridz, 2);

      Z3i::RealPoint previousPoint = pointCloud.at(previousPointIndex);

      double previousDistance = pow(previousPoint[0] - gridx, 2) + pow(previousPoint[1] - gridy, 2) + pow(previousPoint[2] - gridz, 2);

      if (distance < previousDistance){
        indexMap[grdIndex] = i;
      }

    }
  }
}


void simpleSubSample(const std::vector<Z3i::RealPoint> &pointCloud, const double &resolution, std::vector<Z3i::RealPoint> &subsampleCloud){
  assert(pointCloud.size() > 0);
  Z3i::RealPoint p0 = pointCloud.at(0);
  Z3i::RealPoint ptLow = p0;
  Z3i::RealPoint ptUp = p0;
        

  for(unsigned int i = 0; i < pointCloud.size(); i++){
    Z3i::RealPoint pi = pointCloud.at(i);
    ptLow = ptLow.inf(pi);
    ptUp = ptUp.sup(pi);
  }

  double minZ = ptLow[2];
  double maxZ = ptUp[2];

  std::map<unsigned int, unsigned int> selectedFaces;
  getHomogenityCloud(pointCloud, resolution, ptLow, ptUp, selectedFaces);

  for (std::map<unsigned int, unsigned int>::iterator it=selectedFaces.begin(); it!=selectedFaces.end(); ++it){
    unsigned int index = it->second;
    subsampleCloud.push_back(pointCloud[index]);
  }
        
}



int
main(int argc,char **argv)
{
    std::string inputFileName;
    std::string outputPrefix;
    double searchRadius {30.0};
    double accRadius {350.0};
    int voxelSize {30};
    int nbControlPoint {5};
    int nbSegment {10};
    double confThreshold {0.0};
    double scale=1000.;

    CLI::App app;
    app.description("Allowed options are: ");
    app.add_option("-i,--input,1", inputFileName , "input pcl.")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("--voxelSize,-v", voxelSize, "Voxel size", true);
    app.add_option("--scale,-s", scale, "scale parameter for pcl", true);
    app.add_option("--output,-o,2", outputPrefix, "output prefix: output-defect.off, output-def-faces-ids, ...")
      ->required();
    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    trace.info()<<inputFileName<<std::endl;
    /*Step1 : Launch XYZ file into vector of RealPoint of DGtal*/
    trace.info()<<"Step1 : Launch XYZ file into vector of RealPoint of DGtal..."<<std::endl;
    std::vector<Z3i::RealPoint> pointCloud;
    std::ifstream input{inputFileName};
    if (!input.is_open()) {
      std::cout << "Couldn't read file: " << inputFileName << "\n"<<std::endl;
      return 1;
    }
    std::string line;
    std::getline(input, line);
    while(std::getline(input, line)){
      std::istringstream s(line);
      std::string field;
      int i = 0;
      std::vector<double> tempXYZPoints;
      int tempLuminance;
      while (getline(s, field,' ')) {
        if( i==0 || i==1 || i==2 ){
          tempXYZPoints.push_back(atof(field.c_str()));
        }
        i++;
      }
      pointCloud.push_back(Z3i::RealPoint(tempXYZPoints[0],tempXYZPoints[1],tempXYZPoints[2]));
    }
    scaleCloud(pointCloud,scale);
    trace.info()<<"Step2 : Subsample point cloud..."<<std::endl;
    std::vector<Z3i::RealPoint> sampledCloud;
    simpleSubSample(pointCloud, voxelSize, sampledCloud);
    trace.info()<<"taille du nuage de point sampled : "<<sampledCloud.size()<<std::endl;
    /**step 3 : compute centerline with subsampled pointcloud**/
    trace.info()<<"step 3 : compute centerline with subsampled pointcloud"<<std::endl;
    std::vector<Z3i::RealPoint> centerline;
    centerline = CenterlineHelper::computeCenterline(sampledCloud, voxelSize, searchRadius, accRadius, nbControlPoint, nbSegment, confThreshold);
    trace.info()<<"step 4 : write centerline..."<<std::endl;
    export2Text(centerline,outputPrefix+"_centerline.xyz");
    export2Text(sampledCloud,outputPrefix+"_subsampledPCL.xyz");
    return 0;




}
