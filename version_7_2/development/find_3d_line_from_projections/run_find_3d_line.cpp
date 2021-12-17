#include <iostream>

#include <TMath.h>
#include <Math/Vector3D.h>

#include "find_3d_line_from_projections.h"

using namespace std;
using namespace ROOT::Math;

int main()
{
  //double z_x[2] = {-10., -5.};
  //double z_y[2] = {-8., -4.};
  //double x[2] = {2.0, 1.0};
  //double y[2] = {-3.2, -1.6};
  
  line_3d line_upstream, line_midstream;
  
  XYZVector vertex;
  
  double z_x[2] = {-31.4, -1.8};
  double z_y[2] = {-32., -2.4};
  double x[2] = {0.237, 0.291};
  double y[2] = {-0.543, -0.579};
  
  //double z_wanted = 0.169686;
  //double z_wanted = -7.5;

  XYZVector point;
  XYZVector direction;
  
  find_3d_line_from_projections(x, z_x, y, z_y, point, direction);
  
  //cout << point << " " << direction << endl;
  
  //XYZVector vertex;
  
  //vertex.SetCoordinates(point.X() + (z_wanted - point.Z()) * direction.X() / direction.Z(), point.Y() + (z_wanted - point.Z()) * direction.Y() / direction.Z(), z_wanted);
  
  //cout << vertex << endl;
  
  line_upstream.point = point;
  line_upstream.direction = direction;

  z_x[0] = 1.8;
  z_x[1] = 6.6;
  z_y[0] = 1.2;
  z_y[1] = 6.0;
  x[0] = 0.323573;
  x[1] = 0.416571;
  y[0] = -0.585;
  y[1] = -0.627;
  
  find_3d_line_from_projections(x, z_x, y, z_y, point, direction);
  
  line_midstream.point = point;
  line_midstream.direction = direction;
  
  find_closest_point(line_upstream, line_midstream, vertex);
  
  cout << vertex << endl;
  
  line_upstream.point.SetCoordinates(0., 0.001, 0.002);
  line_midstream.point.SetCoordinates(1., 0., 0.);

  line_upstream.direction.SetCoordinates(0., 0., 1.);
  line_midstream.direction.SetCoordinates(0., 0.0174, 1.);
  
  find_closest_point(line_upstream, line_midstream, vertex);
  
  cout << vertex << endl;

  return 0;
}
