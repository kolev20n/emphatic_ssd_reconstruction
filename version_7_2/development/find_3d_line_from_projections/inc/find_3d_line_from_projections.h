#ifndef find_3d_line_from_projections_h
#define find_3d_line_from_projections_h

#include <TMath.h>
#include <Math/Vector3D.h>

class line_3d
{
public:
  ROOT::Math::XYZVector point;
  ROOT::Math::XYZVector direction;
  ROOT::Math::XYZVector error_point;
  ROOT::Math::XYZVector error_direction;
};

void find_3d_line_from_projections(double x[], double z_x[], double y[], double z_y[], ROOT::Math::XYZVector& point, ROOT::Math::XYZVector& direction);

bool find_closest_point(line_3d line_1, line_3d line_2, ROOT::Math::XYZVector& closest_point);


#endif
