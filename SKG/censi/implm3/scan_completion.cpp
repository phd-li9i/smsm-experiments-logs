#include "scan_completion.h"


/*******************************************************************************
*/
void ScanCompletion::completeScan(std::vector<double>* scan, const int& method)
{
  if (method == 1)
    completeScan1(scan);
  else if (method == 3)
    completeScan3(scan);
  else if (method == 4)
    completeScan4(scan);
  else
    completeScan1(scan);
}


/*******************************************************************************
 * Symmetric to the line passing from the first and last scan points
 */
void ScanCompletion::completeScan1(std::vector<double>* scan)
{
  std::vector<double> scan_copy = *scan;

  for (int i = scan_copy.size()-2; i > 0; i--)
    scan->push_back(scan_copy[i]);

  // Rotate so that it starts from -M_PI rather than -M_PI / 2
  int num_pos = scan->size() / 4;

  std::rotate(scan->begin(),
    scan->begin() + scan->size() - num_pos,
    scan->end());
}


/*******************************************************************************
 * Completes a circular arc around 'pose'
 **/
void ScanCompletion::completeScan2(std::vector<double>* scan,
  const std::tuple<double,double,double>& pose)
{
  std::vector<double> scan_copy = *scan;

  // Locate the first and last points of the scan in the 2D plane
  std::vector< std::pair<double,double> > points;
  Utils::scan2points(scan_copy, pose, &points);
  std::pair<double,double> start_point = points[0];
  std::pair<double,double> end_point = points[points.size()-1];

  double dx = start_point.first - end_point.first;
  double dy = start_point.second - end_point.second;
  double d = sqrt(dx*dx + dy*dy);
  double r = d/2;

  for (int i = scan_copy.size()-2; i > 0; i--)
    scan->push_back(r);

  // Rotate so that it starts from -M_PI rather than -M_PI / 2
  int num_pos = scan->size() / 4;

  std::rotate(scan->begin(),
    scan->begin() + scan->size() - num_pos,
    scan->end());
}


/*******************************************************************************
 * Mirrored
 **/
void ScanCompletion::completeScan3(std::vector<double>* scan)
{
  std::vector<double> scan_copy = *scan;

  for (int i = 1; i < scan_copy.size()-1; i++)
    scan->push_back(scan_copy[i]);

  // Rotate so that it starts from -M_PI rather than -M_PI / 2
  int num_pos = scan->size() / 4;

  std::rotate(scan->begin(),
    scan->begin() + scan->size() - num_pos,
    scan->end());
}


/*******************************************************************************
 * min range arc around true pose
 **/
void ScanCompletion::completeScan4(std::vector<double>* scan)
{
  // Find closest and furthest points in original scan
  double min_range = *std::min_element(scan->begin(), scan->end());
  double max_range = *std::max_element(scan->begin(), scan->end());
  double fill_range = min_range;

  unsigned int scan_size = scan->size();

  for (int i = 1; i < scan_size-1; i++)
    scan->push_back(fill_range);

  // Rotate so that it starts from -M_PI rather than -M_PI / 2
  assert(fmod(scan->size(), 2) == 0);
  int num_pos = scan->size() / 4;

  std::rotate(scan->begin(),
    scan->begin() + scan->size() - num_pos,
    scan->end());
}


/*******************************************************************************
  * straight
 **/
void ScanCompletion::completeScan5(
  const std::tuple<double,double,double>& pose,
  const std::vector<double>& scan_in,
  const unsigned int& num_rays,
  std::vector<double>* scan_out,
  std::vector< std::pair<double,double> >* map,
  std::tuple<double,double,double>* map_origin)
{
  std::vector< std::pair<double,double> > scan_points;
  Utils::scan2points(scan_in, pose, &scan_points, M_PI);

  std::tuple<double,double,double> pose_within_points = pose;

  double farther_than = 0.01;
  bool is_farther_than = false;

  while (!is_farther_than)
  {
    do Utils::generatePose(pose,
      0.05, 0.0, &pose_within_points);
    while(!Utils::isPositionInMap(pose_within_points, scan_points));

    *map = X::find(pose_within_points, scan_points, num_rays);

    is_farther_than =
      Utils::isPositionFartherThan(pose_within_points, *map, farther_than);
  }

  *map_origin = pose_within_points;
  Utils::points2scan(*map, *map_origin, scan_out);
}
