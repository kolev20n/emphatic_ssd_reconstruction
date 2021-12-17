#include <iostream>

#include <TMath.h>
#include <Math/Vector3D.h>

#include "find_3d_line_from_projections.h"

using namespace std;
using namespace ROOT::Math;

void find_3d_line_from_projections(double x[], double z_x[], double y[], double z_y[], XYZVector& point, XYZVector& direction)
{
  XYZVector t_x, t_y;
  
  t_x.SetCoordinates(x[1] - x[0], 0., z_x[1] - z_x[0]);
  t_x = t_x / (t_x.Z());
  
  t_y.SetCoordinates(0., y[1] - y[0], z_y[1] - z_y[0]);
  t_y = t_y / (t_y.Z());
  
  direction.SetCoordinates(t_x.X(), t_y.Y(), t_x.Z());
    
  point.SetCoordinates(x[0] - z_x[0] * (x[1] - x[0]) / (z_x[1] - z_x[0]), y[0] - z_y[0] * (y[1] - y[0]) / (z_y[1] - z_y[0]), 0.);
}

bool find_closest_point(line_3d line_1, line_3d line_2, XYZVector& closest_point)
{
  XYZVector d;
  double t1, t2;
  
  line_1.direction = line_1.direction.Unit();
  line_2.direction = line_2.direction.Unit();

  d = line_1.point - line_2.point;
  
  t1 = (line_1.direction.Cross(line_2.direction)).Dot(d.Cross(line_2.direction)) / sqrt((line_1.direction.Cross(line_2.direction)).Mag2());
  
  t2 = (line_1.direction.Cross(line_2.direction)).Dot(d.Cross(line_1.direction)) / sqrt((line_1.direction.Cross(line_2.direction)).Mag2());
  
  closest_point = 0.5 * (line_1.point + t1 * line_1.direction + line_2.point + t2 * line_2.direction);
  
  return true;
}
