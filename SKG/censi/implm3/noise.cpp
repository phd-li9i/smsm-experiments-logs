#include "noise.h"


/*******************************************************************************
*/
std::vector<double> Noise::smooth(const std::vector<double>& scan,
  const double& window_angle)
{
#ifdef TIMES
  std::chrono::high_resolution_clock::time_point a =
    std::chrono::high_resolution_clock::now();
#endif

  assert(window_angle >= 0);
  assert(window_angle <= 2*M_PI);

  double window_length_d = window_angle / (2*M_PI) * scan.size();
  int window_length = static_cast<int>(window_length_d);

#ifdef TIMES
  std::chrono::high_resolution_clock::time_point b =
    std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(b-a);

  printf("%f [smooth]\n", elapsed.count());
#endif

  return smooth(scan, window_length);
}


/*******************************************************************************
*/
std::vector<double> Noise::smooth(const std::vector<double>& scan,
  const int& window_length)
{
  assert(window_length >= 0);
  assert(window_length < scan.size());

  std::vector<double> smooth_scan;

  for (int i = 0; i < scan.size(); i++)
  {
    double smooth_ray = window(scan, window_length, i);
    smooth_scan.push_back(smooth_ray);
  }

  return smooth_scan;
}


/*******************************************************************************
*/
double Noise::window(const std::vector<double>& vec,
  const int& sz, const int& mid_id)
{
  assert(sz >= 0);
  assert(mid_id >= 0 && mid_id < vec.size());
  assert(sz <= vec.size());
  //        ___________
  // HACK: |vec|vec|vec|
  std::vector<double> vec_x3 = vec;
  vec_x3.insert(vec_x3.end(), vec.begin(), vec.end());
  vec_x3.insert(vec_x3.end(), vec.begin(), vec.end());

  // New mid_point
  int m_id = vec.size() + mid_id;

  return accumulate(
    vec_x3.begin() + m_id -(sz/2-1),
    vec_x3.begin() + m_id + sz/2,
    0.0) / sz;
}
