#ifndef CYLINDRICAL_COORDINATE_SYSTEM_H
#define CYLINDRICAL_COORDINATE_SYSTEM_H

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

#include "CylindricalPoint.h"


///////////////////////////////////////////////////////////////////////////////
// class CylindricalCoordinateSystem
/**
 * Description of class 'CylindricalCoordinateSystem' <p>
 * 
 * @brief Helper class to compute the cylindrical coordinate system from a centerline (multisegment line)and a mark point,
 * starting point to compute the theta (theta = 0)
 *
 *
 */


using namespace DGtal;


class CylindricalCoordinateSystem{

// ----------------------- Standard methods ------------------------------
public:
    
    CylindricalCoordinateSystem(const std::vector<Z3i::RealPoint> &cen, const Z3i::RealPoint &seedPoint) : centerline(cen), markPoint(seedPoint){
        init();
    }

    CylindricalPoint xyz2Cylindrical(const Z3i::RealPoint &xyz);

    Z3i::RealPoint getDirectionVector(const unsigned int &segmentId);

    unsigned int getSegment(const Z3i::RealPoint &aPoint);

    //Z3i::RealPoint getRadialVector(const CylindricalPoint &cp);


//protected functions
protected:
    void init();

    void computeThetaAxis();

    void computeOxyNormals();

    Z3i::RealPoint getRadialVector(const Z3i::RealPoint &aPoint, const Z3i::RealPoint &aDirection, const Z3i::RealPoint &p0);

    //////////////////////////////////////////
    //z coordinate
    std::vector<Z3i::RealPoint> centerline;
    //used to determine the theta
    Z3i::RealPoint markPoint;

    std::vector<Z3i::RealPoint> thetaAxis;
    std::vector<Z3i::RealPoint> oxyNormals;
};

#endif //end CENTERLINE_H
