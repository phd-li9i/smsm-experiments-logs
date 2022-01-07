#include "translation.h"


/*******************************************************************************
*/
double Translation::tff(
  const std::vector< double >& real_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  const int& max_iterations,
  const double& dist_bound,
  const bool& pick_min,
  const fftw_plan& r2rp,
  int* result_iterations,
  std::chrono::duration<double>* intersections_time,
  std::tuple<double,double,double>* result_pose)
{
#ifdef PRINTS
  printf("input pose  (%f,%f,%f) [Translation::tff]\n",
    std::get<0>(virtual_pose),
    std::get<1>(virtual_pose),
    std::get<2>(virtual_pose));
#endif

  std::tuple<double,double,double> current_pose = virtual_pose;

  std::vector<double> errors_xy;

  std::vector<double> deltas;
  std::vector<double> sum_d_vs;
  std::vector<double> x_es;
  std::vector<double> y_es;
  double norm_x1;

  // Start the clock
  std::chrono::high_resolution_clock::time_point start =
    std::chrono::high_resolution_clock::now();

  // Iterate
  unsigned int it = 1;
  double inclusion_bound = 1000.0;
  double err = 1.0 / real_scan.size();
  std::vector<double> d_v;
  double sum_d_v = 1.0 / real_scan.size();

  for (it = 1; it <= max_iterations; it++)
  {
    // Measure the time to find intersections
    std::chrono::high_resolution_clock::time_point int_start =
      std::chrono::high_resolution_clock::now();

    // Find the intersections of the rays from the estimated pose and
    // the map.
    std::vector< std::pair<double,double> > virtual_scan_intersections =
      X::find(current_pose, map, real_scan.size());

    std::chrono::high_resolution_clock::time_point int_end =
      std::chrono::high_resolution_clock::now();
    *intersections_time =
      std::chrono::duration_cast< std::chrono::duration<double> >(int_end-int_start);

    // Find the corresponding ranges
    std::vector<double> virtual_scan_it;
    Utils::points2scan(virtual_scan_intersections, current_pose, &virtual_scan_it);

    assert(virtual_scan_it.size() == real_scan.size());

    inclusion_bound = real_scan.size()/4*err;
    //inclusion_bound = 0.01*sum_d;
    //inclusion_bound = M_PI * (sum_d + err) / real_scan.size();
    //inclusion_bound = 2*M_PI * sum_d_v / real_scan.size();

    // Obtain the correction vector
    std::pair<double,double> errors_xy =
      tffCore(real_scan, virtual_scan_it, std::get<2>(current_pose),
        inclusion_bound, r2rp, &d_v, &norm_x1);



    // These are the corrections
    double x_e = errors_xy.first;
    double y_e = errors_xy.second;


    // The norm of the correction vector
    double err_sq = x_e*x_e + y_e*y_e;
    err = sqrt(err_sq);

    // Correct the position
    std::get<0>(current_pose) += x_e;
    std::get<1>(current_pose) += y_e;

    double dx = std::get<0>(current_pose) - std::get<0>(virtual_pose);
    double dy = std::get<1>(current_pose) - std::get<1>(virtual_pose);

    // Check constraints
    if(!Utils::isPositionInMap(current_pose, map) ||
      //fabs(x_e) > 2*dist_bound || fabs(y_e) > 2*dist_bound ||
      fabs(dx) > dist_bound+0.01 || fabs(dy) > dist_bound+0.01)
    {
#ifdef DEBUG
      printf("OUT OF BOUNDS\n");
#endif

      *result_iterations= it;
      *result_pose = current_pose;
      return -2.0;
    }


    //inclusion_bound =
    //pow(2,2)*(sum_d + err_sq)*(2*it + max_iterations) / max_iterations / real_scan.size(); 1125
    //inclusion_bound = pow(2,2) * (sum_d + err_sq) / real_scan.size(); 1142 3436
    //inclusion_bound = pow(2,2) * (sum_d + err) / real_scan.size(); 1144 3407
    //inclusion_bound = pow(2,2) * sum_d / real_scan.size(); 1155 3454
    //inclusion_bound = 0.01*sum_d; 1168 3487
    //inclusion_bound = 100*err;

    for (unsigned int d = 0; d < d_v.size(); d++)
      d_v[d] = fabs(d_v[d]);

    sum_d_v = std::accumulate(d_v.begin(), d_v.end(), 0.0);

#ifdef DEBUG
    printf("err = %f\n", err);
    printf("norm_x1 = %f\n", norm_x1);
    printf("sum_d_v = %f\n", sum_d_v);
#endif


    if (pick_min)
    {
      x_es.push_back(x_e);
      y_es.push_back(y_e);
      sum_d_vs.push_back(sum_d_v);
    }

    // Break if translation is negligible
    double eps = 0.0000001;
    if (fabs(x_e) < eps && fabs(y_e) < eps)
      break;
  }

  if (pick_min)
  {
    std::vector<double> crit_v = sum_d_vs;
    double min_sum_d_idx =
      std::min_element(crit_v.begin(), crit_v.end()) -crit_v.begin();
    sum_d_v = sum_d_vs[min_sum_d_idx];
    double x_tot = std::accumulate(x_es.begin(), x_es.begin()+min_sum_d_idx, 0.0);
    double y_tot = std::accumulate(y_es.begin(), y_es.begin()+min_sum_d_idx, 0.0);

    std::get<0>(*result_pose) = x_tot + std::get<0>(virtual_pose);
    std::get<1>(*result_pose) = y_tot + std::get<1>(virtual_pose);
  }
  else
    *result_pose = current_pose;

  *result_iterations= it;

  // Stop the clock
  std::chrono::high_resolution_clock::time_point end =
    std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(end-start);

#ifdef PRINTS
  printf("output pose (%f,%f,%f) [Translation::tff]\n",
    std::get<0>(*result_pose),
    std::get<1>(*result_pose),
    std::get<2>(*result_pose));
#endif

  return sum_d_v / real_scan.size();
}


/*******************************************************************************
*/
std::pair<double,double> Translation::tffCore(
  const std::vector< double >& real_scan,
  const std::vector< double >& virtual_scan,
  const double& current_t,
  const double& inclusion_bound,
  const fftw_plan& r2rp,
  std::vector<double>* d_v,
  double* norm_x1)
{
  assert(inclusion_bound >= 0);

  std::vector<double> diff;
  Utils::diffScansPerRay(real_scan, virtual_scan, inclusion_bound, &diff, d_v);

  // X1
  std::vector<double> X1 = DFTUtils::getFirstDFTCoefficient(diff, r2rp);

  *norm_x1 = sqrtf(X1[0]*X1[0] + X1[1]*X1[1]);

  // Find the x-wise and y-wise errors
  double t = M_PI + current_t;
  std::vector<double> errors_xy = turnDFTCoeffsIntoErrors(X1, diff.size(), t);

  double x_e = errors_xy[0];
  double y_e = errors_xy[1];

#ifdef DEBUG
  printf("(x_e,y_e) = (%f,%f)\n", x_e, y_e);
#endif

  return std::make_pair(x_e,y_e);
}


/*******************************************************************************
*/
std::vector<double> Translation::turnDFTCoeffsIntoErrors(
  const std::vector<double>& dft_coeff,
  const int& num_valid_rays,
  const double& starting_angle)
{
  double x_err = 0.0;
  double y_err = 0.0;

  if (num_valid_rays > 0)
  {
    // The error in the x- direction
    x_err = 1.0 / num_valid_rays *
      (-dft_coeff[0] * cos(starting_angle)
       -dft_coeff[1] * sin(starting_angle));

    // The error in the y- direction
    y_err = 1.0 / num_valid_rays *
      (-dft_coeff[0] * sin(starting_angle)
       +dft_coeff[1] * cos(starting_angle));
  }

  std::vector<double> errors;
  errors.push_back(x_err);
  errors.push_back(y_err);

  return errors;
}
