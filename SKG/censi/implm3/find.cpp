#include "find.h"


/*******************************************************************************
*/
double Find::area(const std::vector< std::pair<double,double> >& polygon)
{
  double a = 0.0;

  std::vector< std::pair<double,double> > polygon_closed = polygon;
  polygon_closed.push_back(polygon[0]);

  for (int i = 0; i < polygon_closed.size()-1; i++)
  {
    a += (polygon_closed[i+1].first + polygon_closed[i].first)
      * (polygon_closed[i+1].second - polygon_closed[i].second);
  }

  return a/2;
}


/*******************************************************************************
*/
void Find::boundingEllipse(
  const std::vector< std::pair<double,double> >& points,
  std::vector<double>* coefficients)
{
  coefficients->clear();

  // Construct polygons from points
  Polygon_2 poly;
  for (int p = 0; p < points.size(); p++)
    poly.push_back(Point_2(points[p].first, points[p].second));

  // Construct the convex hull
  std::vector<Point_2> convex_hull;
  convex_hull_2(poly.vertices_begin(), poly.vertices_end(),
    std::back_inserter(convex_hull));

  // Construct bounding ellipse
  Min_ellipse m_ellipse(poly.vertices_begin(), poly.vertices_end(), true);

  // Get coefficients for  scan ellipse
  double a0, a1, a2, a3, a4, a5;
  m_ellipse.ellipse().double_coefficients(a0,a1,a2,a3,a4,a5);

  coefficients->push_back(a0);
  coefficients->push_back(a1);
  coefficients->push_back(a2);
  coefficients->push_back(a3);
  coefficients->push_back(a4);
  coefficients->push_back(a5);
}


/*******************************************************************************
*/
std::pair<double,double> Find::centroid(
  const std::vector< std::pair<double,double> >& polygon)
{
  double a = area(polygon);

  std::vector< std::pair<double,double> > polygon_closed = polygon;
  polygon_closed.push_back(polygon[0]);

  double x = 0.0;
  double y = 0.0;

  for (int i = 0; i < polygon_closed.size()-1; i++)
  {
    x += (polygon_closed[i+1].first + polygon_closed[i].first)
      * (polygon_closed[i].first * polygon_closed[i+1].second
        -polygon_closed[i+1].first * polygon_closed[i].second);

    y += (polygon_closed[i+1].second + polygon_closed[i].second)
      * (polygon_closed[i].first * polygon_closed[i+1].second
        -polygon_closed[i+1].first * polygon_closed[i].second);
  }

  return std::make_pair(x / (6*a), y / (6*a));
}


/*******************************************************************************
*/
std::pair<double, double> Find::closestPoint(
  const std::pair<double, double>& point,
  const std::vector< std::pair<double,double> >& points)
{
  int min_idx = -1;
  double min_dist = 10000000000000000000.0;
  for (int i = 0; i < points.size(); i++)
  {
    double dx = point.first - points[i].first;
    double dy = point.second - points[i].second;

    double dist = dx*dx + dy*dy;

    if (dist < min_dist)
    {
      min_dist = dist;
      min_idx = i;
    }
  }

  assert(min_idx >= 0 && min_idx < points.size());

  return points[min_idx];
}


/*******************************************************************************
*/
std::pair<double, double> Find::furthestPoint(
  const std::pair<double, double>& point,
  const std::vector< std::pair<double,double> >& points)
{
  int min_idx = -1;
  double min_dist = -1.0;
  for (int i = 0; i < points.size(); i++)
  {
    double dx = point.first - points[i].first;
    double dy = point.second - points[i].second;

    double dist = dx*dx + dy*dy;

    if (dist >= min_dist)
    {
      min_dist = dist;
      min_idx = i;
    }
  }

  assert(min_idx >= 0 && min_idx < points.size());

  return points[min_idx];
}


/*******************************************************************************
*/
std::pair<double, double> Find::ellipseAxesPoints(
  const std::vector<double>& coefficients)
{
  assert(coefficients.size() == 6);

  double a_p = coefficients[0]; // a_p x^2 +
  double b_p = coefficients[2]; // b_p xy  +
  double c_p = coefficients[1]; // c_p y^2 +
  double d_p = coefficients[3]; // d_p x   +
  double e_p = coefficients[4]; // e_p y   +
  double f_p = coefficients[5]; // f_p     = 0

  double a = a_p;
  double b = b_p / 2;
  double c = c_p;
  double d = d_p / 2;
  double f = e_p / 2;
  double g = f_p;

  double long semi_a =
    sqrt( 2*(a*f*f + c*d*d + g*b*b - 2*b*d*f - a*c*g) /
      ((b*b-a*c)*(sqrt( (a-c)*(a-c) +4*b*b ) -a-c)));

  double long semi_b =
    sqrt( 2*(a*f*f + c*d*d + g*b*b - 2*b*d*f - a*c*g) /
      ((b*b-a*c)*(-sqrt( (a-c)*(a-c) +4*b*b ) -a-c)));

  std::pair<double, double> ret_pair;
  ret_pair.first = semi_a;
  ret_pair.second = semi_b;

  return ret_pair;
}


/*******************************************************************************
*/
std::pair<double, double> Find::ellipseCenter(
  const std::vector<double>& coefficients)
{
  /* The commented section is equivalent to the uncommented one below
     assert(coefficients.size() == 6);

     double a = coefficients[0];
     double b = coefficients[2];
     double c = coefficients[1];
     double d = coefficients[3];
     double e = coefficients[4];
     double f = coefficients[5];

     double t = 2*M_PI;

     if (b == 0.0 && a < c)
     t = 0.0;
     else if (a == 0.0 && a > c)
     t = M_PI/2;
     else if (b != 0.0 && a < c)
     t = 0.5 * atan(b/(a-c));
     else if (a != 0.0 && a > c)
     t = M_PI/2 + 0.5 * atan(b/(a-c));

     double a_p = a * cos(t)*cos(t) + b * cos(t) * sin(t) + c * sin(t)*sin(t);
     double b_p = 0.0;
     double c_p = a * sin(t)*sin(t) - b * cos(t) * sin(t) + c * cos(t)*cos(t);
     double d_p = d * cos(t) + e * sin(t);
     double e_p = -d * sin(t) + e * cos(t);
     double f_p = f;

     double x0_p = -d_p / (2*a_p);
     double y0_p = -e_p / (2*c_p);

     double x0 = x0_p * cos(t) - y0_p * sin(t);
     double y0 = x0_p * sin(t) + y0_p * cos(t);
     */


  double a_p = coefficients[0]; // a_p x^2 +
  double b_p = coefficients[2]; // b_p xy  +
  double c_p = coefficients[1]; // c_p y^2 +
  double d_p = coefficients[3]; // d_p x   +
  double e_p = coefficients[4]; // e_p y   +
  double f_p = coefficients[5]; // f_p     = 0

  double a = a_p;
  double b = b_p / 2;
  double c = c_p;
  double d = d_p / 2;
  double f = e_p / 2;
  double g = f_p;

  double x0 = (c*d - b*f) / (b*b - a*c);
  double y0 = (a*f - b*d) / (b*b - a*c);

  std::pair<double, double> ret_pair;
  ret_pair.first = x0;
  ret_pair.second = y0;

  return ret_pair;
}


/*******************************************************************************
*/
double Find::ellipseAngle(
  const std::vector<double>& coefficients)
{
  assert(coefficients.size() == 6);

  double a = coefficients[0];
  double b = coefficients[2];
  double c = coefficients[1];
  double d = coefficients[3];
  double e = coefficients[4];
  double f = coefficients[5];

  double t;

  if (fabs(a-c) > 0.00001)
    t = 0.5 * atan(b/(a-c));
  else
    t = 0.0;

  return t;
}


/*******************************************************************************
*/
  std::vector< std::pair<double,double> > Find::points2convexHullPoints
(const std::vector< std::pair<double,double> >& points)
{
  // Construct CGAL polygon from points
  Polygon_2 poly;
  for (int p = 0; p < points.size(); p++)
    poly.push_back(Point_2(points[p].first, points[p].second));

  // Construct the CGAL convex hull
  std::vector<Point_2> convex_hull_cgal;
  convex_hull_2(poly.vertices_begin(), poly.vertices_end(),
    std::back_inserter(convex_hull_cgal));

  // CGAL convex hull to points
  std::vector< std::pair<double,double> > convex_hull;
  for (unsigned int i = 0; i < convex_hull_cgal.size(); i++)
  {
    convex_hull.push_back(
      std::make_pair(convex_hull_cgal[i].x(), convex_hull_cgal[i].y()));
  }

  return convex_hull;
}


/*******************************************************************************
*/
void Find::scansFromConvexHull(
  const std::vector<double>& real_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  std::vector<double>* real_scan2,
  std::vector<double>* virtual_scan2)
{
  std::tuple<double,double,double> zero_pose;
  std::get<0>(zero_pose) = 0.0;
  std::get<1>(zero_pose) = 0.0;
  std::get<2>(zero_pose) = 0.0;

  // Real scan points
  std::vector< std::pair<double,double> > real_scan_points;
  Utils::scan2points(real_scan, zero_pose, &real_scan_points);

  // Virtual scan points
  std::vector< std::pair<double,double> > virtual_scan_points =
    X::find(virtual_pose, map, real_scan.size());

  // Convex hull of real scan points
  std::vector< std::pair<double,double> > real_convex_hull =
    points2convexHullPoints(real_scan_points);

  // Convex hull of virtual scan points
  std::vector< std::pair<double,double> > virtual_convex_hull =
    points2convexHullPoints(virtual_scan_points);

  // The intersections of rays and convex hulls
  std::vector< std::pair<double,double> > real_scan_points2 =
    X::find(zero_pose, real_convex_hull, real_scan.size());
  std::vector< std::pair<double,double> > virtual_scan_points2 =
    X::find(virtual_pose, virtual_convex_hull, real_scan.size());

  // Their corresponding scans
  Utils::points2scan(real_scan_points2, zero_pose, real_scan2);
  Utils::points2scan(virtual_scan_points2, virtual_pose, virtual_scan2);
}
