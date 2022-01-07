#include "utils.h"

/*******************************************************************************
*/
std::pair<double,double> Utils::computeDeltaXY(
  const std::vector< std::pair<double,double> >& real_scan_points,
  const std::vector< std::pair<double,double> >& virtual_scan_points)
{
  assert(real_scan_points.size() == virtual_scan_points.size());

  unsigned int N = real_scan_points.size();

  double delta_x = 0.0;
  double delta_y = 0.0;
  for (int i = 0; i < N; i++)
  {
    delta_x += real_scan_points[i].first - virtual_scan_points[i].first;
    delta_y += real_scan_points[i].second - virtual_scan_points[i].second;
  }

  delta_x /= N;
  delta_y /= N;

  return std::make_pair(delta_x, delta_y);
}


/*******************************************************************************
*/
std::pair<double,double> Utils::computeDeltaXY(
  const std::vector<double>& real_scan,
  const std::tuple<double,double,double>& real_pose,
  const std::vector<double>& virtual_scan,
  const std::tuple<double,double,double>& virtual_pose)
{
  assert(real_scan.size() == virtual_scan.size());

  double rx0 = std::get<0>(real_pose);
  double ry0 = std::get<1>(real_pose);
  double rt0 = std::get<2>(real_pose);

  double vx0 = std::get<0>(virtual_pose);
  double vy0 = std::get<1>(virtual_pose);
  double vt0 = std::get<2>(virtual_pose);

  unsigned int N = real_scan.size();

  double delta_x = 0.0;
  double delta_y = 0.0;
  for (int i = 0; i < real_scan.size(); i++)
  {
    double x_r = rx0 + real_scan[i]*cos(-M_PI + i*2*M_PI/N + rt0);
    double y_r = ry0 + real_scan[i]*sin(-M_PI + i*2*M_PI/N + rt0);

    double x_v = vx0 + virtual_scan[i]*cos(-M_PI + i*2*M_PI/N + vt0);
    double y_v = vy0 + virtual_scan[i]*sin(-M_PI + i*2*M_PI/N + vt0);

    delta_x += x_r - x_v;
    delta_y += y_r - y_v;
  }

  //delta_x /= N;
  //delta_y /= N;

  return std::make_pair(delta_x, delta_y);
}


/*******************************************************************************
*/
std::vector< std::pair<double,double> > Utils::conjugate(
  const std::vector< std::pair<double,double> >& vec)
{
#ifdef TIMES
  std::chrono::high_resolution_clock::time_point start =
    std::chrono::high_resolution_clock::now();
#endif

  std::vector< std::pair<double,double> > ret_vector;
  for (int i = 0; i < vec.size(); i++)
    ret_vector.push_back(std::make_pair(vec[i].first, -vec[i].second));

#ifdef TIMES
  std::chrono::high_resolution_clock::time_point end =
    std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(end-start);

  printf("%f [conjugate]\n", elapsed.count());
#endif

  return ret_vector;
}


/*******************************************************************************
*/
void Utils::diffScansPerRay(
  const std::vector<double>& scan1, const std::vector<double>& scan2,
  const double& inclusion_bound, std::vector<double>* diff,
  std::vector<double>* diff_true)
{
  assert (scan1.size() == scan2.size());

  diff->clear();
  diff_true->clear();

  double eps = 0.000001;
  if (inclusion_bound < 0.0001)
    eps = 1.0;

#ifdef DEBUG
  printf("inclusion_bound = %f\n", inclusion_bound + eps);
#endif

  double d = 0.0;
  for (unsigned int i = 0; i < scan1.size(); i++)
  {
    d = scan1[i] - scan2[i];

    if (fabs(d) <= inclusion_bound + eps)
      diff->push_back(d);
    else
      diff->push_back(0.0);

    diff_true->push_back(d);
  }
}


/*******************************************************************************
*/
void Utils::generatePose(
  const std::tuple<double,double,double>& real_pose,
  const double& dxy, const double& dt,
  std::tuple<double,double,double>* virtual_pose)
{
  assert(dxy >= 0);
  assert(dt >= 0);

  std::random_device rand_dev;
  std::mt19937 generator_x(rand_dev());
  std::mt19937 generator_y(rand_dev());
  std::mt19937 generator_t(rand_dev());
  std::mt19937 generator_sign(rand_dev());

  std::uniform_real_distribution<double> distribution_x(-dxy, dxy);
  std::uniform_real_distribution<double> distribution_y(-dxy, dxy);
  std::uniform_real_distribution<double> distribution_t(-dt, dt);

  double rx = distribution_x(generator_x);
  double ry = distribution_y(generator_y);
  double rt = distribution_t(generator_t);

  std::get<0>(*virtual_pose) = std::get<0>(real_pose) + rx;
  std::get<1>(*virtual_pose) = std::get<1>(real_pose) + ry;
  std::get<2>(*virtual_pose) = std::get<2>(real_pose) + rt;

  Utils::wrapAngle(&std::get<2>(*virtual_pose));
}


/*******************************************************************************
*/
bool Utils::generatePose(
  const std::tuple<double,double,double>& base_pose,
  const std::vector< std::pair<double,double> >& map,
  const double& dxy, const double& dt, const double& dist_threshold,
  std::tuple<double,double,double>* real_pose)
{
  assert(dxy >= 0.0);
  assert(dt >= 0.0);

  std::random_device rand_dev_x;
  std::random_device rand_dev_y;
  std::random_device rand_dev_t;
  std::mt19937 generator_x(rand_dev_x());
  std::mt19937 generator_y(rand_dev_y());
  std::mt19937 generator_t(rand_dev_t());

  std::uniform_real_distribution<double> distribution_x(-dxy, dxy);
  std::uniform_real_distribution<double> distribution_y(-dxy, dxy);
  std::uniform_real_distribution<double> distribution_t(-dt, dt);

  // A temp real pose
  std::tuple<double,double,double> real_pose_ass;

  // Fill in the orientation regardless
  double rt = distribution_t(generator_t);
  std::get<2>(real_pose_ass) = std::get<2>(base_pose) + rt;
  double t = std::get<2>(real_pose_ass);
  Utils::wrapAngle(&t);
  std::get<2>(real_pose_ass) = t;

  // We assume that the lidar sensor is distanced from the closest obstacle
  // by a certain amount (e.g. the radius of a circular base)
  bool pose_found = false;
  while (!pose_found)
  {
    pose_found = true;
    double rx = distribution_x(generator_x);
    double ry = distribution_y(generator_y);

    std::get<0>(real_pose_ass) = std::get<0>(base_pose) + rx;
    std::get<1>(real_pose_ass) = std::get<1>(base_pose) + ry;

    if (isPositionInMap(real_pose_ass, map))
    {
      for (unsigned int i = 0; i < map.size(); i++)
      {
        double dx = std::get<0>(real_pose_ass) - map[i].first;
        double dy = std::get<1>(real_pose_ass) - map[i].second;

        if (dx*dx + dy*dy < dist_threshold*dist_threshold)
        {
          pose_found = false;
          break;
        }
      }
    }
    else pose_found = false;
  }

  *real_pose = real_pose_ass;

  // Verify distance threshold
  std::vector< std::pair<double,double> > intersections =
    X::find(real_pose_ass, map, map.size());
  std::vector<double> real_scan;
  points2scan(intersections, real_pose_ass, &real_scan);

  unsigned int min_dist_idx =
    std::min_element(real_scan.begin(), real_scan.end()) - real_scan.begin();

  return real_scan[min_dist_idx] > dist_threshold;
}


/*******************************************************************************
*/
bool Utils::generatePoseWithinMap(
  const std::vector< std::pair<double,double> >& map,
  const double& dist_threshold,
  std::tuple<double,double,double>* pose)
{
  // A temp real pose
  std::tuple<double,double,double> real_pose_ass;

  // Generate orientation
  std::random_device rand_dev_t;
  std::mt19937 generator_t(rand_dev_t());

  std::uniform_real_distribution<double> distribution_t(-M_PI, M_PI);

  // Fill in the orientation regardless
  std::get<2>(real_pose_ass) = distribution_t(generator_t);

  // Find the bounding box of the map
  double max_x = -1000.0;
  double min_x = +1000.0;
  double max_y = -1000.0;
  double min_y = +1000.0;

  for (unsigned int i = 0; i < map.size(); i++)
  {
    if (map[i].first > max_x)
      max_x = map[i].first;

    if (map[i].first < min_x)
      min_x = map[i].first;

    if (map[i].second > max_y)
      max_y = map[i].second;

    if (map[i].second < min_y)
      min_y = map[i].second;
  }

  std::random_device rand_dev_x;
  std::random_device rand_dev_y;
  std::mt19937 generator_x(rand_dev_x());
  std::mt19937 generator_y(rand_dev_y());

  std::uniform_real_distribution<double> distribution_x(min_x, max_x);
  std::uniform_real_distribution<double> distribution_y(min_y, max_y);

  // We assume that the lidar sensor is distanced from the closest obstacle
  // by a certain amount (e.g. the radius of a circular base)
  bool pose_found = false;
  while (!pose_found)
  {
    pose_found = true;
    double rx = distribution_x(generator_x);
    double ry = distribution_y(generator_y);

    std::get<0>(real_pose_ass) = rx;
    std::get<1>(real_pose_ass) = ry;

    if (isPositionInMap(real_pose_ass, map))
    {
      for (unsigned int i = 0; i < map.size(); i++)
      {
        double dx = std::get<0>(real_pose_ass) - map[i].first;
        double dy = std::get<1>(real_pose_ass) - map[i].second;

        if (dx*dx + dy*dy < dist_threshold*dist_threshold)
        {
          pose_found = false;
          break;
        }
      }
    }
    else pose_found = false;
  }

  *pose = real_pose_ass;

  // Verify distance threshold
  std::vector< std::pair<double,double> > intersections =
    X::find(real_pose_ass, map, map.size());
  std::vector<double> real_scan;
  points2scan(intersections, real_pose_ass, &real_scan);

  unsigned int min_dist_idx =
    std::min_element(real_scan.begin(), real_scan.end()) - real_scan.begin();

  return real_scan[min_dist_idx] > dist_threshold;
}


/*******************************************************************************
*/
std::vector<double> Utils::innerProduct(
  const std::vector<double>& vec1,
  const std::vector<double>& vec2)
{
  assert(vec1.size() == vec2.size());

  std::vector<double> ret_vector;

  for (int i = 0; i < vec1.size(); i++)
  {
    ret_vector.push_back(vec1[i] * vec2[i]);
  }

  return ret_vector;
}


/*******************************************************************************
*/
std::vector< std::pair<double, double> > Utils::innerProductComplex(
  const std::vector< std::pair<double, double> >& vec1,
  const std::vector< std::pair<double, double> >& vec2)
{
#ifdef TIMES
  std::chrono::high_resolution_clock::time_point start =
    std::chrono::high_resolution_clock::now();
#endif

  assert(vec1.size() == vec2.size());

  std::vector< std::pair<double, double> > ret_vector;

  for (int i = 0; i < vec1.size(); i++)
  {
    double re =
      vec1[i].first * vec2[i].first - vec1[i].second * vec2[i].second;
    double im =
      vec1[i].first * vec2[i].second + vec1[i].second * vec2[i].first;

    ret_vector.push_back(std::make_pair(re,im));
  }

#ifdef TIMES
  std::chrono::high_resolution_clock::time_point end =
    std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(end-start);

  printf("%f [innerProductComplex]\n", elapsed.count());
#endif

  return ret_vector;
}


/*******************************************************************************
*/
bool Utils::isPositionInMap(
  const std::tuple<double, double, double>& pose,
  const std::vector< std::pair<double,double> >& map)
{
  Point_2 point(std::get<0>(pose), std::get<1>(pose));

  // Construct polygon from map
  Polygon_2 poly;
  for (int p = 0; p < map.size(); p++)
    poly.push_back(Point_2(map[p].first, map[p].second));

  poly.push_back(Point_2(map[map.size()-1].first, map[map.size()-1].second));

  bool inside = false;
  if(CGAL::bounded_side_2(poly.vertices_begin(),
      poly.vertices_end(),
      point, Kernel()) == CGAL::ON_BOUNDED_SIDE)
  {
    inside = true;
  }

  return inside;
}


/*******************************************************************************
*/
bool Utils::isPositionFartherThan(
  const std::tuple<double, double, double>& pose,
  const std::vector< std::pair<double,double> >& map,
  const double& dist)
{
  for (unsigned int i = 0; i < map.size(); i++)
  {
    double dx = std::get<0>(pose) - map[i].first;
    double dy = std::get<1>(pose) - map[i].second;
    double d = sqrt(dx*dx + dy*dy);

    if (d < dist)
      return false;
  }

  return true;
}


/*******************************************************************************
*/
std::pair<double,double> Utils::multiplyWithRotationMatrix(
  const std::pair<double,double>& point, const double& angle)
{
  double R11 = cos(angle);
  double R12 = -sin(angle);
  double R21 = -R12;
  double R22 = R11;

  double x = R11 * point.first + R12 * point.second;
  double y = R21 * point.first + R22 * point.second;

  return std::make_pair(x,y);

}


/*******************************************************************************
*/
std::vector< std::pair<double,double> > Utils::multiplyWithRotationMatrix(
  const std::vector< std::pair<double,double> >& points,
  const double& angle)
{
  std::vector< std::pair<double,double> > return_vector;

  for (int i = 0; i < points.size(); i++)
    return_vector.push_back(multiplyWithRotationMatrix(points[i], angle));

  return return_vector;
}


/*******************************************************************************
*/
double Utils::norm(const std::pair<double,double>& vec)
{
  return sqrt(vec.first*vec.first + vec.second*vec.second);
}


/*******************************************************************************
*/
std::vector<double> Utils::norm(
  const std::vector< std::pair<double,double> >& vec)
{
  std::vector<double> ret_vector;

  for (int i = 0; i < vec.size(); i++)
    ret_vector.push_back(norm(vec[i]));

  return ret_vector;
}


/*******************************************************************************
*/
double Utils::norm2(
  const std::vector< std::pair<double,double> >& vec)
{
  std::vector<double> ret_vector;

  for (int i = 0; i < vec.size(); i++)
    ret_vector.push_back(norm(vec[i]));

  return accumulate(ret_vector.begin(), ret_vector.end(), 0.0);
}


/*******************************************************************************
*/
std::pair<double,double> Utils::pairDiff(
  const std::pair<double,double>& pair1,
  const std::pair<double,double>& pair2)
{
  std::pair<double,double> ret_pair;
  ret_pair.first = pair2.first - pair1.first;
  ret_pair.second = pair2.second - pair1.second;

  return ret_pair;
}


/*******************************************************************************
*/
void Utils::points2scan(
  const std::vector< std::pair<double,double> >& points,
  const std::tuple<double,double,double>& pose,
  std::vector<double>* scan)
{
#ifdef TIMES
  std::chrono::high_resolution_clock::time_point start =
    std::chrono::high_resolution_clock::now();
#endif

  scan->clear();

  double px = std::get<0>(pose);
  double py = std::get<1>(pose);

  double dx = 0.0;
  double dy = 0.0;
  for (int i = 0; i < points.size(); i++)
  {
    dx = points[i].first - px;
    dy = points[i].second - py;
    scan->push_back(sqrt(dx*dx+dy*dy));
  }

#ifdef TIMES
  std::chrono::high_resolution_clock::time_point end =
    std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(end-start);

  printf("%f [points2scan]\n", elapsed.count());
#endif
}


/*******************************************************************************
*/
void Utils::scan2points(
  const std::vector<double>& scan,
  const std::tuple<double,double,double> pose,
  std::vector< std::pair<double,double> >* points,
  const double& angle_span)
{
#ifdef TIMES
  std::chrono::high_resolution_clock::time_point start =
    std::chrono::high_resolution_clock::now();
#endif

  points->clear();

  double px = std::get<0>(pose);
  double py = std::get<1>(pose);
  double pt = std::get<2>(pose);

  // The angle of the first ray (in the local coordinate system)
  double sa = -angle_span/2;

  for (int i = 0; i < scan.size(); i++)
  {
    double x =
      px + scan[i] * cos(i * angle_span / scan.size() + pt + sa);
    double y =
      py + scan[i] * sin(i * angle_span / scan.size() + pt + sa);

    points->push_back(std::make_pair(x,y));
  }

#ifdef TIMES
  std::chrono::high_resolution_clock::time_point end =
    std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(end-start);

  printf("%f [scan2points]\n", elapsed.count());
#endif
}


/*******************************************************************************
*/
void Utils::scanFromPose(
  const std::tuple<double,double,double>& pose,
  const std::vector< std::pair<double,double> >& points,
  const unsigned int& num_rays,
  std::vector<double>* scan)
{
  scan->clear();

  std::vector< std::pair<double,double> > intersections =
    X::find(pose, points, num_rays);

  points2scan(intersections, pose, scan);
}


/*******************************************************************************
*/
int Utils::sgn(const double& a)
{
  return (a > 0.0) - (a < 0.0);
}


/*******************************************************************************
*/
std::vector< std::pair<double,double> > Utils::vectorDiff(
  const std::vector< std::pair<double,double> >& vec)
{

  std::vector< std::pair<double,double> > ret_vector;

  for (int i = 0; i < vec.size()-1; i++)
    ret_vector.push_back(pairDiff(vec[i], vec[i+1]));

  return ret_vector;
}


/*******************************************************************************
*/
std::pair<double,double> Utils::vectorStatistics(
  const std::vector< double >& v)
{
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  double mean = sum / v.size();

  std::vector<double> diff(v.size());
  std::transform(v.begin(), v.end(), diff.begin(),
    std::bind2nd(std::minus<double>(), mean));
  double sq_sum =
    std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double stdev = std::sqrt(sq_sum / v.size());

  return std::make_pair(mean, stdev);
}


/*******************************************************************************
*/
void Utils::wrapAngle(double* angle)
{
  *angle = fmod(*angle + 5*M_PI, 2*M_PI) - M_PI;
}
