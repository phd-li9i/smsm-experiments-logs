#include "rotation.h"

/*******************************************************************************
*/
double Rotation::angleById(const unsigned int& rotation_id,
  const unsigned int scan_size)
{
  double dt = 2*M_PI*rotation_id / scan_size;

  Utils::wrapAngle(&dt);

  return dt;
}

/*******************************************************************************
*/
void Rotation::coarse(
  const std::vector< double >& real_scan,
  const std::vector<double>& real_scan_points_vectors_norm,
  const double& real_ip,
  const std::tuple<double,double,double>& current_pose,
  const std::vector< std::pair<double,double> >& map,
  const unsigned int& num_rays,
  std::tuple<double,double,double>* result_pose)
{
  *result_pose = current_pose;

#if defined (PRINTS)
  printf("input pose  (%f,%f,%f) [coarse rotation]\n",
    std::get<0>(*result_pose),
    std::get<1>(*result_pose),
    std::get<2>(*result_pose)
    );
#endif

  // For diff purposes
  double prev_error = 1000.0;
  double prev_dt_error = 1000.0;
  std::tuple<double,double,double> prev_pose = current_pose;
  std::get<0>(prev_pose) += 1000.0;
  std::get<1>(prev_pose) += 1000.0;
  std::get<2>(prev_pose) += 1000.0;

  // For projection/abstraction purposes
  std::tuple<double,double,double> zero_pose;
  std::get<0>(zero_pose) = 0.0;
  std::get<1>(zero_pose) = 0.0;
  std::get<2>(zero_pose) = 0.0;



  for (int it = 0; it < 1; it++)
  {
    // Find the intersections of the rays from the estimated pose and
    // the map.
    std::vector< std::pair<double,double> > virtual_scan_intersections =
      X::find(*result_pose, map, num_rays);

    // Find the corresponding ranges

    // The map scan for each iteration
    std::vector<double> virtual_scan_it;
    Utils::points2scan(virtual_scan_intersections, *result_pose, &virtual_scan_it);

    assert(virtual_scan_it.size() == real_scan.size());


    // Find how many shifts need to be made to the virtual scan points *vector*
    // in order for it to be matched against the real scan points vector
    unsigned int rotation_id = findRotationId(real_scan, virtual_scan_it,
      real_scan_points_vectors_norm, real_ip, 0);

    // Calculate the angle corresponding to this shift.
    // max(rotation_angle - actual rotation angle) = k * lidar angle increment
    double rotation_angle = angleById(rotation_id, virtual_scan_it.size());

    // Rotate the virtual pose by that angle
    if (fabs(rotation_angle) < M_PI/4)
      std::get<2>(*result_pose) += rotation_angle;

    double error_dt_it = std::get<2>(*result_pose) - std::get<2>(prev_pose);
    if (fabs(prev_dt_error - error_dt_it) < 0.000001)
      break;

    prev_pose = *result_pose;
    prev_dt_error = std::get<2>(*result_pose) - std::get<2>(prev_pose);
  }


#if defined (PRINTS)
  printf("output pose (%f,%f,%f) [coarse rotation]\n",
    std::get<0>(*result_pose),
    std::get<1>(*result_pose),
    std::get<2>(*result_pose)
    );
#endif
}


/*******************************************************************************
*/
std::vector<double> Rotation::dbh(
  const std::vector< double >& real_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  const unsigned int& magnification_size,
  const std::string& batch_or_sequential,
  const fftw_plan& r2rp, const fftw_plan& c2rp,
  std::vector<double>* rc0, std::vector<double>* rc1,
  std::chrono::duration<double>* intersections_time)
{
  if (batch_or_sequential.compare("batch") == 0)
   return dbh2Batch(real_scan, virtual_pose, map, magnification_size,
     r2rp, c2rp, rc0, rc1, intersections_time);
  else
    if (batch_or_sequential.compare("sequential") == 0)
   return dbh2Sequential(real_scan, virtual_pose, map, magnification_size,
     rc0, rc1, intersections_time);
  else
    printf("[Rotation::dbh] Use 'batch' or 'sequential' instead \n");
}


/*******************************************************************************
*/
std::vector<double> Rotation::dbh2Sequential(
  const std::vector< double >& real_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  const unsigned int& magnification_size,
  std::vector<double>* rc0, std::vector<double>* rc1,
  std::chrono::duration<double>* intersections_time)
{
#if defined (PRINTS)
  printf("input pose  (%f,%f,%f) [Rotation::dbh2Sequential]\n",
    std::get<0>(virtual_pose),
    std::get<1>(virtual_pose),
    std::get<2>(virtual_pose));
#endif

  rc0->clear();
  rc1->clear();

  std::tuple<double,double,double> zero_pose;
  std::get<0>(zero_pose) = 0.0;
  std::get<1>(zero_pose) = 0.0;
  std::get<2>(zero_pose) = 0.0;


  unsigned int num_virtual_scans = pow(2,magnification_size);
  int virtual_scan_size_max = num_virtual_scans * real_scan.size();

  // Measure the time to find intersections
  std::chrono::high_resolution_clock::time_point int_start =
    std::chrono::high_resolution_clock::now();

  std::vector< std::pair<double,double> > virtual_scan_points =
    X::find(virtual_pose, map, virtual_scan_size_max);

  std::chrono::high_resolution_clock::time_point int_end =
    std::chrono::high_resolution_clock::now();
  *intersections_time =
    std::chrono::duration_cast< std::chrono::duration<double> >(int_end-int_start);

  std::vector<double> virtual_scan_fine;
  Utils::points2scan(virtual_scan_points, virtual_pose, &virtual_scan_fine);

  // Downsample from upper limit:
  // construct the upper-most resolution and downsample from there.
  std::vector< std::vector< double> > virtual_scans(num_virtual_scans);

  for (int i = 0; i < virtual_scan_fine.size(); i++)
  {
    unsigned int k = fmod(i,num_virtual_scans);
    virtual_scans[k].push_back(virtual_scan_fine[i]);
  }

  // Make sure that all virtual scans are equal to the real scan in terms of
  // size
  for (unsigned int i = 0; i < virtual_scans.size(); i++)
    assert(virtual_scans[i].size() == real_scan.size());

  // The real scan's (the original) angle increment
  double ang_inc = 2*M_PI / real_scan.size();
  double mul = 1.0 / num_virtual_scans;


  std::vector<double> orientations;
  std::vector<double> snrs;
  std::vector<double> fahms;
  std::vector<double> pds;
  std::vector<double> rot_criteria;

  std::vector< std::pair<double,double> > real_scan_points;
  Utils::scan2points(real_scan, zero_pose, &real_scan_points);

  for (unsigned int a = 0; a < num_virtual_scans; a++)
  {
    std::vector< std::pair<double,double> > virtual_scan_points_a;
    Utils::scan2points(virtual_scans[a], zero_pose, &virtual_scan_points_a);

    double angle = 0.0;
    double snr = 0.0;
    double fahm = 0.0;
    double pd = 0.0;

    dbh1Sequential(
      real_scan_points, virtual_scan_points_a, &angle, &snr, &fahm, &pd);

    double ornt_a = -angle + a*mul*ang_inc;
    Utils::wrapAngle(&ornt_a);

    orientations.push_back(ornt_a);
    snrs.push_back(snr);
    fahms.push_back(fahm);
    pds.push_back(pd);

#if defined (DEBUG)
    printf("a = %u\n", a);
    printf("angle to out = %f\n", std::get<2>(virtual_pose) + ornt_a);
    printf("snr = %.10f\n", snr);
    printf("fahm = %f\n", fahm);
    printf("pd = %.20f\n", pd);
#endif
  }

  // Select some of all the angles based on criteria enforced by rankDBHOutput
  std::vector<unsigned int> optimal_ids =
    rankDBHOutput(snrs, fahms, pds, 3, magnification_size, 0.00001);

  std::vector<double> angles;
  for (unsigned int i = 0; i < optimal_ids.size(); i++)
  {
    double angle = orientations[optimal_ids[i]];
    Utils::wrapAngle(&angle);
    angles.push_back(angle);

    rc0->push_back(pds[optimal_ids[i]]);
    rc1->push_back(snrs[optimal_ids[i]] / fahms[optimal_ids[i]]);
  }

#if defined (PRINTS)
  for (unsigned int i = 0; i < angles.size(); i++)
  {
    printf("cand. poses (%f,%f,%f) [Rotation::dbh2Sequential]\n",
      std::get<0>(virtual_pose),
      std::get<1>(virtual_pose),
      std::get<2>(virtual_pose)+angles[i]);
  }
#endif

  return angles;
}


/*******************************************************************************
*/
void Rotation::dbh1Sequential(
  const std::vector< std::pair<double,double> >& real_scan_points,
  const std::vector< std::pair<double,double> >& virtual_scan_points,
  double* angle, double* snr, double* fahm, double* pd)
{
  std::vector<double> traces;
  unsigned int traces_max_id;
  dbh0Sequential(real_scan_points, virtual_scan_points, &traces, &traces_max_id);

  // Calculate angle -----------------------------------------------------------
  int rot_id = traces_max_id;
  *angle = static_cast<double>(
    (real_scan_points.size()-rot_id))/(real_scan_points.size())*2*M_PI;

  Utils::wrapAngle(angle);

  // Calculate SNR -------------------------------------------------------------
  std::vector<double> traces_background = traces;
  traces_background.erase(traces_background.begin() + traces_max_id);

  std::pair<double,double> traces_mmnts =
    Utils::vectorStatistics(traces_background);

  *snr = fabs((traces[traces_max_id] - traces_mmnts.first)) / traces_mmnts.second;

  // Calculate FAHM ------------------------------------------------------------
  unsigned int count = 0;
  for (unsigned int i = 0; i < traces.size(); i++)
  {
    if (traces[i] >= 0.5 * traces[traces_max_id])
      count++;
  }

  *fahm = static_cast<double>(count) / traces.size();

  // Calculate PD --------------------------------------------------------------
  std::vector<double> traces_ss;
  unsigned int traces_ss_max_id;
  dbh0AutoSequential(real_scan_points, &traces_ss, &traces_ss_max_id);

  std::vector<double> traces_rr;
  unsigned int traces_rr_max_id;
  dbh0AutoSequential(virtual_scan_points, &traces_rr, &traces_rr_max_id);

  *pd = 2*traces[traces_max_id] /
    (traces_ss[traces_ss_max_id] + traces_rr[traces_rr_max_id]);
}


/*******************************************************************************
*/
void Rotation::dbh0Sequential(
  const std::vector< std::pair<double,double> >& real_scan_points,
  const std::vector< std::pair<double,double> >& virtual_scan_points_in,
  std::vector<double>* traces, unsigned int* traces_max_id)
{
  traces->clear();

  std::vector<double> B11;
  std::vector<double> B12;
  std::vector<double> B21;
  std::vector<double> B22;

  std::vector< std::pair<double,double> > virtual_scan_points =
    virtual_scan_points_in;
  std::reverse(virtual_scan_points.begin(), virtual_scan_points.end());

  for (int i = 0; i < real_scan_points.size(); i++)
  {
    B11.push_back(virtual_scan_points[i].first);
    B12.push_back(virtual_scan_points[i].second);

    B21.push_back(real_scan_points[i].first);
    B22.push_back(real_scan_points[i].second);
  }

  std::vector<double> B11_coeffs_ar = DFTUtils::dft(B11);
  std::vector<double> B12_coeffs_ar = DFTUtils::dft(B12);
  std::vector<double> B21_coeffs_ar = DFTUtils::dft(B21);
  std::vector<double> B22_coeffs_ar = DFTUtils::dft(B22);

  std::vector< std::pair<double, double> > B11_coeffs =
    DFTUtils::getDFTCoefficientsPairs(B11_coeffs_ar);
  std::vector< std::pair<double, double> > B12_coeffs =
    DFTUtils::getDFTCoefficientsPairs(B12_coeffs_ar);
  std::vector< std::pair<double, double> > B21_coeffs =
    DFTUtils::getDFTCoefficientsPairs(B21_coeffs_ar);
  std::vector< std::pair<double, double> > B22_coeffs =
    DFTUtils::getDFTCoefficientsPairs(B22_coeffs_ar);

  std::vector<double> A11 =
    DFTUtils::idft(Utils::innerProductComplex(B11_coeffs, B21_coeffs));
  std::vector<double> A12 =
    DFTUtils::idft(Utils::innerProductComplex(B11_coeffs, B22_coeffs));
  std::vector<double> A21 =
    DFTUtils::idft(Utils::innerProductComplex(B12_coeffs, B21_coeffs));
  std::vector<double> A22 =
    DFTUtils::idft(Utils::innerProductComplex(B12_coeffs, B22_coeffs));

  for (int i = 0; i < A11.size(); i++)
  {
    // Only consider angular deviations +/- 45deg TODO SPEEDUP
    //if (i > A11.size()/4 && i < 3*A11.size()/4)
    //continue;

    Eigen::Matrix2d A(2,2);
    A(0,0) = A11[i];
    A(0,1) = A12[i];
    A(1,0) = A21[i];
    A(1,1) = A22[i];

    Eigen::JacobiSVD<Eigen::Matrix2d> svd_of_A(A,
      Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::Matrix2d _S_;
    if (svd_of_A.matrixU().determinant() * svd_of_A.matrixV().determinant() > 0)
    {
      _S_(0,0) = 1;
      _S_(1,0) = 0;
      _S_(0,1) = 0;
      _S_(1,1) = 1;
    }
    else
    {
      _S_(0,0) = 1;
      _S_(1,0) = 0;
      _S_(0,1) = 0;
      _S_(1,1) = -1;
    }

    Eigen::Matrix2d R =
      svd_of_A.matrixU() * _S_ * svd_of_A.matrixV().transpose();


    Eigen::Matrix2d RAt = R * A.transpose();

    traces->push_back(RAt.trace());
  }

  // Reverse the reversion
  std::reverse(traces->begin(), traces->end());

  *traces_max_id =
    std::max_element(traces->begin(), traces->end()) - traces->begin();
}


/*******************************************************************************
*/
void Rotation::dbh0AutoSequential(
  const std::vector< std::pair<double,double> >& real_scan_points,
  std::vector<double>* traces, unsigned int* traces_max_id)
{
  std::vector< std::vector<double> > B_v;

  // Reverse the input virtual scan points
  std::vector< std::pair<double,double> > real_scan_points_inv =
    real_scan_points;
  std::reverse(real_scan_points_inv.begin(), real_scan_points_inv.end());

  std::vector<double> B11;
  std::vector<double> B12;
  std::vector<double> B21;
  std::vector<double> B22;

  for (int i = 0; i < real_scan_points.size(); i++)
  {
    B11.push_back(real_scan_points_inv[i].first);
    B12.push_back(real_scan_points_inv[i].second);

    B21.push_back(real_scan_points[i].first);
    B22.push_back(real_scan_points[i].second);
  }

  B_v.push_back(B11);
  B_v.push_back(B12);
  B_v.push_back(B21);
  B_v.push_back(B22);

  std::vector< std::vector<double> > B_coeffs_ar = DFTUtils::dftBatch(B_v);

  std::vector< std::vector< std::pair<double, double> > > BB_v;

  for (unsigned int b = 0; b < B_coeffs_ar.size(); b=b+4)
  {
    std::vector< std::pair<double, double> > B11_coeffs =
      DFTUtils::getDFTCoefficientsPairs(B_coeffs_ar[b+0]);
    std::vector< std::pair<double, double> > B12_coeffs =
      DFTUtils::getDFTCoefficientsPairs(B_coeffs_ar[b+1]);
    std::vector< std::pair<double, double> > B21_coeffs =
      DFTUtils::getDFTCoefficientsPairs(B_coeffs_ar[b+2]);
    std::vector< std::pair<double, double> > B22_coeffs =
      DFTUtils::getDFTCoefficientsPairs(B_coeffs_ar[b+3]);

    BB_v.push_back(Utils::innerProductComplex(B11_coeffs, B21_coeffs));
    BB_v.push_back(Utils::innerProductComplex(B11_coeffs, B22_coeffs));
    BB_v.push_back(Utils::innerProductComplex(B12_coeffs, B21_coeffs));
    BB_v.push_back(Utils::innerProductComplex(B12_coeffs, B22_coeffs));
  }

  std::vector< std::vector<double> > A_v = DFTUtils::idftBatch(BB_v);

  for (unsigned int a = 0; a < A_v.size(); a=a+4)
  {
    std::vector<double> A11 = A_v[a+0];
    std::vector<double> A12 = A_v[a+1];
    std::vector<double> A21 = A_v[a+2];
    std::vector<double> A22 = A_v[a+3];

    for (int i = 0; i < A11.size(); i++)
    {
      // Only consider angular deviations +/- 45deg TODO SPEEDUP
      //if (i > A11.size()/4 && i < 3*A11.size()/4)
      //continue;

      Eigen::Matrix2d A(2,2);
      A(0,0) = A11[i];
      A(0,1) = A12[i];
      A(1,0) = A21[i];
      A(1,1) = A22[i];

      Eigen::JacobiSVD<Eigen::Matrix2d> svd_of_A(A,
        Eigen::ComputeFullU | Eigen::ComputeFullV);

      Eigen::Matrix2d _S_;
      if (svd_of_A.matrixU().determinant() * svd_of_A.matrixV().determinant() > 0)
      {
        _S_(0,0) = 1;
        _S_(1,0) = 0;
        _S_(0,1) = 0;
        _S_(1,1) = 1;
      }
      else
      {
        _S_(0,0) = 1;
        _S_(1,0) = 0;
        _S_(0,1) = 0;
        _S_(1,1) = -1;
      }

    Eigen::Matrix2d R =
      svd_of_A.matrixU() * _S_ * svd_of_A.matrixV().transpose();


      Eigen::Matrix2d RAt = R * A.transpose();

      traces->push_back(RAt.trace());
    }

    // Revert the reversion
    std::reverse(traces->begin(), traces->end());

    *traces_max_id =
      std::max_element(traces->begin(), traces->end()) - traces->begin();
  }
}


/*******************************************************************************
*/
std::vector<double> Rotation::dbh2Batch(
  const std::vector< double >& real_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  const unsigned int& magnification_size,
  const fftw_plan& r2rp, const fftw_plan& c2rp,
  std::vector<double>* rc0, std::vector<double>* rc1,
  std::chrono::duration<double>* intersections_time)
{
#if defined (PRINTS)
  printf("input pose  (%f,%f,%f) [Rotation::dbh2Batch]\n",
    std::get<0>(virtual_pose),
    std::get<1>(virtual_pose),
    std::get<2>(virtual_pose));
#endif

  rc0->clear();
  rc1->clear();

  std::tuple<double,double,double> zero_pose;
  std::get<0>(zero_pose) = 0.0;
  std::get<1>(zero_pose) = 0.0;
  std::get<2>(zero_pose) = 0.0;


  unsigned int num_virtual_scans = pow(2,magnification_size);
  int virtual_scan_size_max = num_virtual_scans * real_scan.size();

  // Measure the time to find intersections
  std::chrono::high_resolution_clock::time_point int_start =
    std::chrono::high_resolution_clock::now();

  std::vector< std::pair<double,double> > virtual_scan_points =
    X::find(virtual_pose, map, virtual_scan_size_max);

  std::chrono::high_resolution_clock::time_point int_end =
    std::chrono::high_resolution_clock::now();
  *intersections_time =
    std::chrono::duration_cast< std::chrono::duration<double> >(int_end-int_start);

  std::vector<double> virtual_scan_fine;
  Utils::points2scan(virtual_scan_points, virtual_pose, &virtual_scan_fine);

  // Downsample from upper limit:
  // construct the upper-most resolution and downsample from there.
  std::vector< std::vector< double> > virtual_scans(num_virtual_scans);

  for (int i = 0; i < virtual_scan_fine.size(); i++)
  {
    unsigned int k = fmod(i,num_virtual_scans);
    virtual_scans[k].push_back(virtual_scan_fine[i]);
  }

  // Make sure that all virtual scans are equal to the real scan in terms of
  // size
  for (unsigned int i = 0; i < virtual_scans.size(); i++)
    assert(virtual_scans[i].size() == real_scan.size());

  // Turn scans into points
  std::vector< std::pair<double,double> > real_scan_points;
  Utils::scan2points(real_scan, zero_pose, &real_scan_points);

  std::vector< std::vector< std::pair<double,double> > > virtual_scan_points_v;

  for (unsigned int a = 0; a < num_virtual_scans; a++)
  {
    std::vector< std::pair<double,double> > virtual_scan_points_a;
    Utils::scan2points(virtual_scans[a], zero_pose, &virtual_scan_points_a);
    virtual_scan_points_v.push_back(virtual_scan_points_a);
  }

  // The real scan's (the original) angle increment
  double ang_inc = 2*M_PI / real_scan.size();
  double mul = 1.0 / num_virtual_scans;


  // Compute the angles and metrics of matching the real scan against each and
  // all virtual scans
  std::vector<double> un_angles;
  std::vector<double> snrs;
  std::vector<double> fahms;
  std::vector<double> pds;

  dbh1Batch(real_scan_points, virtual_scan_points_v, r2rp, c2rp,
    &un_angles, &snrs, &fahms, &pds);

  // Correct the angles returned to get the proper pose from which each
  // virtual scan was taken (needed due to over-sampling the map)
  std::vector<double> angles;
  for (unsigned int a = 0; a < num_virtual_scans; a++)
  {
    double ornt_a = -un_angles[a] + a*mul*ang_inc;
    Utils::wrapAngle(&ornt_a);

    angles.push_back(ornt_a);
  }

  // Select some of all the angles based on criteria enforced by rankFMTOutput
  std::vector<unsigned int> optimal_ids =
    rankDBHOutput(snrs, fahms, pds, 3, magnification_size, 0.00001);

  std::vector<double> cand_angles;
  for (unsigned int i = 0; i < optimal_ids.size(); i++)
  {
    cand_angles.push_back(angles[optimal_ids[i]]);

    rc0->push_back(pds[optimal_ids[i]]);
    rc1->push_back(snrs[optimal_ids[i]] / fahms[optimal_ids[i]]);
  }

#if defined (TIMES)
  printf("%f [Rotation::dbh2Batch]\n", elapsed.count());
#endif

#if defined (PRINTS)
  for (unsigned int i = 0; i < cand_angles.size(); i++)
  {
    printf("cand. poses (%f,%f,%f) [Rotation::dbh2Batch]\n",
      std::get<0>(virtual_pose),
      std::get<1>(virtual_pose),
      std::get<2>(virtual_pose)+cand_angles[i]);
  }
#endif

  return cand_angles;
}


/*******************************************************************************
*/
void Rotation::dbh1Batch(
  const std::vector< std::pair<double,double> >& real_scan_points,
  const std::vector< std::vector< std::pair<double,double> > >& virtual_scan_points_v,
  const fftw_plan& r2rp, const fftw_plan& c2rp,
  std::vector<double>* angles,
  std::vector<double>* snrs,
  std::vector<double>* fahms,
  std::vector<double>* pds)
{
  std::vector< std::vector<double> > traces_v;
  std::vector< unsigned int > traces_max_id_v;
  dbh0Batch(real_scan_points, virtual_scan_points_v, r2rp, c2rp,
    &traces_v, &traces_max_id_v);

  // Calculate PD --------------------------------------------------------------
  std::vector<double> traces_ss;
  unsigned int traces_ss_max_id;
  dbh0AutoSequential(real_scan_points, &traces_ss, &traces_ss_max_id);

  std::vector< std::vector<double> > traces_rr_v;
  std::vector< unsigned int > traces_rr_max_id_v;

  dbh0AutoBatch(virtual_scan_points_v, r2rp, c2rp,
    &traces_rr_v, &traces_rr_max_id_v);

  for (unsigned int i = 0; i < traces_v.size(); i++)
  {
    // Calculate angle ---------------------------------------------------------
    int rot_id = traces_max_id_v[i];
    double angle = static_cast<double>(
      (real_scan_points.size()-rot_id))/(real_scan_points.size())*2*M_PI;

    Utils::wrapAngle(&angle);
    angles->push_back(angle);

    // Calculate pd ------------------------------------------------------------
    double pd = 2*traces_v[i][traces_max_id_v[i]] /
      (traces_ss[traces_ss_max_id] + traces_rr_v[i][traces_rr_max_id_v[i]]);

    pds->push_back(pd);

    // Calculate SNR -----------------------------------------------------------
    std::vector<double> traces_background = traces_v[i];
    traces_background.erase(traces_background.begin() + traces_max_id_v[i]);

    std::pair<double,double> traces_mmnts =
      Utils::vectorStatistics(traces_background);

    double snr =
      fabs((traces_v[i][traces_max_id_v[i]] - traces_mmnts.first)) / traces_mmnts.second;
    snrs->push_back(snr);

    // Calculate FAHM ----------------------------------------------------------
    unsigned int count = 0;
    for (unsigned int t = 0; t < traces_v[i].size(); t++)
    {
      if (traces_v[i][t] >= 0.5 * traces_v[i][traces_max_id_v[i]])
        count++;
    }

    double fahm = static_cast<double>(count) / traces_v[i].size();
    fahms->push_back(fahm);
  }
}


/*******************************************************************************
*/
void Rotation::dbh0Batch(
  const std::vector< std::pair<double,double> >& real_scan_points,
  const std::vector< std::vector< std::pair<double,double> > >&
  virtual_scan_points_in_v,
  const fftw_plan& r2rp, const fftw_plan& c2rp,
  std::vector< std::vector<double> >* traces_v,
  std::vector< unsigned int>* traces_max_id_v)
{
  std::vector< std::vector<double> > B_v;

  for (unsigned int v = 0; v < virtual_scan_points_in_v.size(); v++)
  {
    // Reverse the input virtual scan points
    std::vector< std::pair<double,double> > virtual_scan_points =
      virtual_scan_points_in_v[v];
    std::reverse(virtual_scan_points.begin(), virtual_scan_points.end());

    std::vector<double> B11;
    std::vector<double> B12;
    std::vector<double> B21;
    std::vector<double> B22;

    for (int i = 0; i < virtual_scan_points.size(); i++)
    {
      B11.push_back(virtual_scan_points[i].first);
      B12.push_back(virtual_scan_points[i].second);

      B21.push_back(real_scan_points[i].first);
      B22.push_back(real_scan_points[i].second);
    }

    B_v.push_back(B11);
    B_v.push_back(B12);
    B_v.push_back(B21);
    B_v.push_back(B22);
  }

  std::vector< std::vector<double> > B_coeffs_ar =
    DFTUtils::dftBatch(B_v, r2rp);


  std::vector< std::vector< std::pair<double, double> > > BB_v;

  for (unsigned int b = 0; b < B_coeffs_ar.size(); b=b+4)
  {
    std::vector< std::pair<double, double> > B11_coeffs =
      DFTUtils::getDFTCoefficientsPairs(B_coeffs_ar[b+0]);
    std::vector< std::pair<double, double> > B12_coeffs =
      DFTUtils::getDFTCoefficientsPairs(B_coeffs_ar[b+1]);
    std::vector< std::pair<double, double> > B21_coeffs =
      DFTUtils::getDFTCoefficientsPairs(B_coeffs_ar[b+2]);
    std::vector< std::pair<double, double> > B22_coeffs =
      DFTUtils::getDFTCoefficientsPairs(B_coeffs_ar[b+3]);

    BB_v.push_back(Utils::innerProductComplex(B11_coeffs, B21_coeffs));
    BB_v.push_back(Utils::innerProductComplex(B11_coeffs, B22_coeffs));
    BB_v.push_back(Utils::innerProductComplex(B12_coeffs, B21_coeffs));
    BB_v.push_back(Utils::innerProductComplex(B12_coeffs, B22_coeffs));
  }

  std::vector< std::vector<double> > A_v = DFTUtils::idftBatch(BB_v, c2rp);

  for (unsigned int a = 0; a < A_v.size(); a=a+4)
  {
    std::vector<double> A11 = A_v[a+0];
    std::vector<double> A12 = A_v[a+1];
    std::vector<double> A21 = A_v[a+2];
    std::vector<double> A22 = A_v[a+3];

    std::vector<double> traces;
    for (int i = 0; i < A11.size(); i++)
    {
      // Only consider angular deviations +/- 45deg TODO SPEEDUP
      //if (i > A11.size()/4 && i < 3*A11.size()/4)
      //continue;

      Eigen::Matrix2d A(2,2);
      A(0,0) = A11[i];
      A(0,1) = A12[i];
      A(1,0) = A21[i];
      A(1,1) = A22[i];

      Eigen::JacobiSVD<Eigen::Matrix2d> svd_of_A(A,
        Eigen::ComputeFullU | Eigen::ComputeFullV);

      Eigen::Matrix2d _S_;
      if (svd_of_A.matrixU().determinant() * svd_of_A.matrixV().determinant() > 0)
      {
        _S_(0,0) = 1;
        _S_(1,0) = 0;
        _S_(0,1) = 0;
        _S_(1,1) = 1;
      }
      else
      {
        _S_(0,0) = 1;
        _S_(1,0) = 0;
        _S_(0,1) = 0;
        _S_(1,1) = -1;
      }

    Eigen::Matrix2d R =
      svd_of_A.matrixU() * _S_ * svd_of_A.matrixV().transpose();


      Eigen::Matrix2d RAt = R * A.transpose();

      traces.push_back(RAt.trace());
    }

    // Revert the reversion
    std::reverse(traces.begin(), traces.end());

    unsigned int traces_max_id =
      std::max_element(traces.begin(), traces.end()) - traces.begin();

    traces_v->push_back(traces);
    traces_max_id_v->push_back(traces_max_id);
  }
}



/*******************************************************************************
*/
void Rotation::dbh0AutoBatch(
  const std::vector< std::vector< std::pair<double,double> > >&
  virtual_scan_points_in_v,
  const fftw_plan& r2rp, const fftw_plan& c2rp,
  std::vector< std::vector<double> >* traces_v,
  std::vector< unsigned int>* traces_max_id_v)
{
  std::vector< std::vector<double> > B_v;

  for (unsigned int v = 0; v < virtual_scan_points_in_v.size(); v++)
  {
    // Reverse the input virtual scan points
    std::vector< std::pair<double,double> > virtual_scan_points =
      virtual_scan_points_in_v[v];
    std::reverse(virtual_scan_points.begin(), virtual_scan_points.end());

    std::vector<double> B11;
    std::vector<double> B12;
    std::vector<double> B21;
    std::vector<double> B22;

    for (int i = 0; i < virtual_scan_points.size(); i++)
    {
      B11.push_back(virtual_scan_points[i].first);
      B12.push_back(virtual_scan_points[i].second);

      B21.push_back(virtual_scan_points_in_v[v][i].first);
      B22.push_back(virtual_scan_points_in_v[v][i].second);
    }

    B_v.push_back(B11);
    B_v.push_back(B12);
    B_v.push_back(B21);
    B_v.push_back(B22);
  }

  std::vector< std::vector<double> > B_coeffs_ar =
    DFTUtils::dftBatch(B_v, r2rp);


  std::vector< std::vector< std::pair<double, double> > > BB_v;

  for (unsigned int b = 0; b < B_coeffs_ar.size(); b=b+4)
  {
    std::vector< std::pair<double, double> > B11_coeffs =
      DFTUtils::getDFTCoefficientsPairs(B_coeffs_ar[b+0]);
    std::vector< std::pair<double, double> > B12_coeffs =
      DFTUtils::getDFTCoefficientsPairs(B_coeffs_ar[b+1]);
    std::vector< std::pair<double, double> > B21_coeffs =
      DFTUtils::getDFTCoefficientsPairs(B_coeffs_ar[b+2]);
    std::vector< std::pair<double, double> > B22_coeffs =
      DFTUtils::getDFTCoefficientsPairs(B_coeffs_ar[b+3]);

    BB_v.push_back(Utils::innerProductComplex(B11_coeffs, B21_coeffs));
    BB_v.push_back(Utils::innerProductComplex(B11_coeffs, B22_coeffs));
    BB_v.push_back(Utils::innerProductComplex(B12_coeffs, B21_coeffs));
    BB_v.push_back(Utils::innerProductComplex(B12_coeffs, B22_coeffs));
  }

  std::vector< std::vector<double> > A_v = DFTUtils::idftBatch(BB_v, c2rp);

  for (unsigned int a = 0; a < A_v.size(); a=a+4)
  {
    std::vector<double> A11 = A_v[a+0];
    std::vector<double> A12 = A_v[a+1];
    std::vector<double> A21 = A_v[a+2];
    std::vector<double> A22 = A_v[a+3];

    std::vector<double> traces;
    for (int i = 0; i < A11.size(); i++)
    {
      // Only consider angular deviations +/- 45deg TODO SPEEDUP
      //if (i > A11.size()/4 && i < 3*A11.size()/4)
      //continue;

      Eigen::Matrix2d A(2,2);
      A(0,0) = A11[i];
      A(0,1) = A12[i];
      A(1,0) = A21[i];
      A(1,1) = A22[i];

      Eigen::JacobiSVD<Eigen::Matrix2d> svd_of_A(A,
        Eigen::ComputeFullU | Eigen::ComputeFullV);

      Eigen::Matrix2d _S_;
      if (svd_of_A.matrixU().determinant() * svd_of_A.matrixV().determinant() > 0)
      {
        _S_(0,0) = 1;
        _S_(1,0) = 0;
        _S_(0,1) = 0;
        _S_(1,1) = 1;
      }
      else
      {
        _S_(0,0) = 1;
        _S_(1,0) = 0;
        _S_(0,1) = 0;
        _S_(1,1) = -1;
      }

    Eigen::Matrix2d R =
      svd_of_A.matrixU() * _S_ * svd_of_A.matrixV().transpose();


      Eigen::Matrix2d RAt = R * A.transpose();

      traces.push_back(RAt.trace());
    }

    // Revert the reversion
    std::reverse(traces.begin(), traces.end());

    unsigned int traces_max_id =
      std::max_element(traces.begin(), traces.end()) - traces.begin();

    traces_v->push_back(traces);
    traces_max_id_v->push_back(traces_max_id);
  }
}


/*******************************************************************************
*/
unsigned int Rotation::findRotationId(
  const std::vector<double>& real_scan,
  const std::vector<double>& virtual_scan_it,
  const std::vector<double>& real_scan_points_vectors_norm,
  const double& real_ip,
  const unsigned int& rotation_fashion)
{
  unsigned int rotation_id = 0;

  std::tuple<double,double,double> zero_pose;
  std::get<0>(zero_pose) = 0.0;
  std::get<1>(zero_pose) = 0.0;
  std::get<2>(zero_pose) = 0.0;

  // cpp
  if (rotation_fashion == 0)
  {
    std::vector< std::pair<double,double> > virtual_scan_points;
    Utils::scan2points(virtual_scan_it, zero_pose, &virtual_scan_points);

    double min_error = 10000.0;
    int min_error_idx = -1;
    for (int s = 0; s < virtual_scan_points.size(); s++)
    {
      rotate(virtual_scan_points.begin(),
        virtual_scan_points.begin()+1,
        virtual_scan_points.end());

      // Do no search for angles over pi/4 or under -pi/4
      if (s > virtual_scan_points.size()/4 && s < 3*virtual_scan_points.size()/4)
        continue;

      // virtual scan points as vectors
      std::vector< std::pair<double,double> >
        virtual_scan_points_vectors = Utils::vectorDiff(virtual_scan_points);

      // normed
      std::vector<double> virtual_scan_points_vectors_norm =
        Utils::norm(virtual_scan_points_vectors);

      // The inner product between them
      std::vector<double> virtual_scan_points_vectors_norm_ip =
        Utils::innerProduct(real_scan_points_vectors_norm,
          virtual_scan_points_vectors_norm);

      // A total sum
      double virtual_ip = accumulate(
        virtual_scan_points_vectors_norm_ip.begin(),
        virtual_scan_points_vectors_norm_ip.end(), 0.0);

      double error = fabs(real_ip - virtual_ip);
      if (error < min_error)
      {
        min_error = error;
        min_error_idx = s+1;
      }
    }

    assert(min_error_idx != -1);

    // Uncomment when translation is absent
    //min_error_idx -= 4;

    //printf("min_error_id = %d\n", min_error_idx);

    rotation_id = min_error_idx;
  }

  // Octave
  /*
     if (rotation_fashion == 1)
     {
     Dump::scan(real_scan, zero_pose, virtual_scan_it, zero_pose,
     base_path_ + "/../matlab/scan_dump.m");

     octave_value_list full_path;
     full_path(0) = "/home/li9i/Desktop/fs2msm/matlab/";
     feval("addpath", full_path);

     octave_value_list input_arguments;

     const octave_value_list result =
     feval("function_rotation_diff", input_arguments);

     rotation_id = result(0).int_value();
     }
     */

  return rotation_id;
}


/*******************************************************************************
*/
std::vector<double> Rotation::fmt(
  const std::vector< double >& real_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  const unsigned int& magnification_size,
  const std::string& batch_or_sequential,
  const fftw_plan& r2rp, const fftw_plan& c2rp,
  std::vector<double>* rc0, std::vector<double>* rc1,
  std::chrono::duration<double>* intersections_time)
{
  if (batch_or_sequential.compare("batch") == 0)
    return fmt2Batch(real_scan, virtual_pose, map, magnification_size,
      r2rp, c2rp, rc0, rc1, intersections_time);
  else if (batch_or_sequential.compare("sequential") == 0)
    return fmt2Sequential(real_scan, virtual_pose, map, magnification_size,
      rc0, rc1, intersections_time);
  else
    printf("[Rotation::fmt] Use 'batch' or 'sequential' instead \n");
}




/*******************************************************************************
*/
std::vector<double> Rotation::fmt2Sequential(
  const std::vector< double >& real_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  const unsigned int& magnification_size,
  std::vector<double>* rc0, std::vector<double>* rc1,
  std::chrono::duration<double>* intersections_time)
{
#if defined (PRINTS)
  printf("input pose  (%f,%f,%f) [Rotation::fmt2]\n",
    std::get<0>(virtual_pose),
    std::get<1>(virtual_pose),
    std::get<2>(virtual_pose));
#endif

  rc0->clear();
  rc1->clear();

  std::tuple<double,double,double> zero_pose;
  std::get<0>(zero_pose) = 0.0;
  std::get<1>(zero_pose) = 0.0;
  std::get<2>(zero_pose) = 0.0;


  unsigned int num_virtual_scans = pow(2,magnification_size);
  int virtual_scan_size_max = num_virtual_scans * real_scan.size();

  // Measure the time to find intersections
  std::chrono::high_resolution_clock::time_point int_start =
    std::chrono::high_resolution_clock::now();

  std::vector< std::pair<double,double> > virtual_scan_points =
    X::find(virtual_pose, map, virtual_scan_size_max);

  std::chrono::high_resolution_clock::time_point int_end =
    std::chrono::high_resolution_clock::now();
  *intersections_time =
    std::chrono::duration_cast< std::chrono::duration<double> >(int_end-int_start);

  std::vector<double> virtual_scan_fine;
  Utils::points2scan(virtual_scan_points, virtual_pose, &virtual_scan_fine);

  // Downsample from upper limit:
  // construct the upper-most resolution and downsample from there.
  std::vector< std::vector< double> > virtual_scans(num_virtual_scans);

  for (int i = 0; i < virtual_scan_fine.size(); i++)
  {
    unsigned int k = fmod(i,num_virtual_scans);
    virtual_scans[k].push_back(virtual_scan_fine[i]);
  }

  // Make sure that all virtual scans are equal to the real scan in terms of
  // size
  for (unsigned int i = 0; i < virtual_scans.size(); i++)
    assert(virtual_scans[i].size() == real_scan.size());

  // The real scan's (the original) angle increment
  double ang_inc = 2*M_PI / real_scan.size();
  double mul = 1.0 / num_virtual_scans;


  std::vector<double> orientations;
  std::vector<double> snrs;
  std::vector<double> fahms;
  std::vector<double> pds;

  for (unsigned int a = 0; a < num_virtual_scans; a++)
  {
    double angle = 0.0;
    double snr = 0.0;
    double fahm = 0.0;
    double pd = 0.0;

    fmt1Sequential(real_scan, virtual_scans[a], &angle, &snr, &fahm, &pd);

    double ornt_a = -angle + a*mul*ang_inc;
    Utils::wrapAngle(&ornt_a);

    orientations.push_back(ornt_a);
    snrs.push_back(snr);
    fahms.push_back(fahm);
    pds.push_back(pd);

#if defined (DEBUG)
    printf("a = %u\n", a);
    printf("angle to out = %f\n", std::get<2>(virtual_pose) + ornt_a);
    printf("snr = %.10f\n", snr);
    printf("fahm = %f\n", fahm);
    printf("pd = %.20f\n", pd);
#endif
  }

  // Select some of all the angles based on criteria enforced by rankFMTOutput
  std::vector<unsigned int> optimal_ids =
    rankFMTOutput(snrs, fahms, pds, 3, magnification_size, 0.00001);

  std::vector<double> angles;
  for (unsigned int i = 0; i < optimal_ids.size(); i++)
  {
    double angle = orientations[optimal_ids[i]];
    Utils::wrapAngle(&angle);
    angles.push_back(angle);

    rc0->push_back(pds[optimal_ids[i]]);
    rc1->push_back(snrs[optimal_ids[i]] / fahms[optimal_ids[i]]);
  }

#if defined (TIMES)
  printf("%f [Rotation::fmt2]\n", elapsed.count());
#endif

#if defined (PRINTS)
  for (unsigned int i = 0; i < angles.size(); i++)
  {
    printf("cand. poses (%f,%f,%f) [Rotation::fmt2]\n",
      std::get<0>(virtual_pose),
      std::get<1>(virtual_pose),
      std::get<2>(virtual_pose)+angles[i]);
  }
#endif

  return angles;
}


/*******************************************************************************
*/
void Rotation::fmt1Sequential(
  const std::vector< double >& real_scan,
  const std::vector< double >& virtual_scan,
  double* angle, double* snr, double* fahm, double* pd)
{
  std::vector<double> q_0;
  unsigned int q_0_max_id;
  fmt0Sequential(real_scan, virtual_scan, &q_0, &q_0_max_id);


  // Calculate angle -----------------------------------------------------------
  *angle = static_cast<double>(
    (real_scan.size()-q_0_max_id))/(real_scan.size())*2*M_PI;
  Utils::wrapAngle(angle);

  // Calculate SNR -------------------------------------------------------------
  std::vector<double> q_0_background = q_0;
  q_0_background.erase(q_0_background.begin() + q_0_max_id);

  std::pair<double,double> q_0_mmnts = Utils::vectorStatistics(q_0_background);

  *snr = fabs((q_0[q_0_max_id] - q_0_mmnts.first)) / q_0_mmnts.second;

  // Calculate FAHM ------------------------------------------------------------
  unsigned int count = 0;
  for (unsigned int i = 0; i < q_0.size(); i++)
  {
    if (q_0[i] >= 0.5 * q_0[q_0_max_id])
      count++;
  }

  *fahm = static_cast<double>(count) / q_0.size();

  // Calculate PD --------------------------------------------------------------
  std::vector<double> q_ss;
  unsigned int q_ss_max_id;
  fmt0Sequential(real_scan, real_scan, &q_ss, &q_ss_max_id);

  std::vector<double> q_rr;
  unsigned int q_rr_max_id;
  fmt0Sequential(virtual_scan, virtual_scan, &q_rr, &q_rr_max_id);

  *pd = 2*q_0[q_0_max_id] / (q_ss[q_ss_max_id] + q_rr[q_rr_max_id]);
}


/*******************************************************************************
*/
void Rotation::fmt0Sequential(
  const std::vector< double >& real_scan,
  const std::vector< double >& virtual_scan,
  std::vector<double>* q_0,
  unsigned int* q_0_max_id)
{
  assert(real_scan.size() == virtual_scan.size());

  // Find fft of real scan
  std::vector<double> fft_real = DFTUtils::dft(real_scan);
  //DFTUtils::fftshift(&fft_real);

  // Find fft of virtual scan
  std::vector<double> fft_virtual = DFTUtils::dft(virtual_scan);
  //DFTUtils::fftshift(&fft_virtual);

  // fft_real is in halfcomplex format; fft_real_coeffs is in normal format
  // (you get the full complex transform)
  std::vector< std::pair<double, double> > fft_real_coeffs =
    DFTUtils::getDFTCoefficientsPairs(fft_real);
  std::vector< std::pair<double, double> > fft_virtual_coeffs =
    DFTUtils::getDFTCoefficientsPairs(fft_virtual);

  // Find conjugates of real coefficients
  std::vector< std::pair<double, double> > fft_real_coeffs_conj =
    Utils::conjugate(fft_real_coeffs);

  // The numerator of Q_0
  std::vector< std::pair<double, double> > numerator =
    Utils::innerProductComplex(fft_real_coeffs_conj, fft_virtual_coeffs);

  // The denominator of Q_0
  double denominator =
    Utils::norm2(fft_real_coeffs) * Utils::norm2(fft_virtual_coeffs);

  /*
     for (int i = 0; i < numerator.size(); i++)
     {
     numerator[i].first /= denominator;
     numerator[i].second /= denominator;
     }
     */

  std::vector< std::pair<double, double> > Q_0 = numerator;

  *q_0 = DFTUtils::idft(Q_0);

  *q_0_max_id = std::max_element(q_0->begin(), q_0->end()) - q_0->begin();
}


/*******************************************************************************
*/
void Rotation::fmt0AutoSequential(
  const std::vector< double >& real_scan,
  std::vector<double>* q_0,
  unsigned int* q_0_max_id)
{
  // Find fft of real scan
  std::vector<double> fft_real = DFTUtils::dft(real_scan);

  std::vector< std::pair<double, double> > fft_real_coeffs =
    DFTUtils::getDFTCoefficientsPairs(fft_real);

  // Find conjugates of real coefficients
  std::vector< std::pair<double, double> > fft_real_coeffs_conj =
    Utils::conjugate(fft_real_coeffs);

  // The numerator of Q_0
  std::vector< std::pair<double, double> > numerator =
    Utils::innerProductComplex(fft_real_coeffs_conj, fft_real_coeffs);

  std::vector< std::pair<double, double> > Q_0 = numerator;

  *q_0 = DFTUtils::idft(Q_0);

  *q_0_max_id = std::max_element(q_0->begin(), q_0->end()) - q_0->begin();
}


/*******************************************************************************
*/
std::vector<double> Rotation::fmt2Batch(
  const std::vector< double >& real_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  const unsigned int& magnification_size,
  const fftw_plan& r2rp, const fftw_plan& c2rp,
  std::vector<double>* rc0, std::vector<double>* rc1,
  std::chrono::duration<double>* intersections_time)
{
#if defined(PRINTS)
  printf("input pose  (%f,%f,%f) [Rotation::fmt2]\n",
    std::get<0>(virtual_pose),
    std::get<1>(virtual_pose),
    std::get<2>(virtual_pose));
#endif

  rc0->clear();
  rc1->clear();

  std::tuple<double,double,double> zero_pose;
  std::get<0>(zero_pose) = 0.0;
  std::get<1>(zero_pose) = 0.0;
  std::get<2>(zero_pose) = 0.0;


  unsigned int num_virtual_scans = pow(2,magnification_size);
  int virtual_scan_size_max = num_virtual_scans * real_scan.size();

  // Measure the time to find intersections
  std::chrono::high_resolution_clock::time_point int_start =
    std::chrono::high_resolution_clock::now();

  std::vector< std::pair<double,double> > virtual_scan_points =
    X::find(virtual_pose, map, virtual_scan_size_max);

  std::chrono::high_resolution_clock::time_point int_end =
    std::chrono::high_resolution_clock::now();
  *intersections_time =
    std::chrono::duration_cast< std::chrono::duration<double> >(int_end-int_start);

  std::vector<double> virtual_scan_fine;
  Utils::points2scan(virtual_scan_points, virtual_pose, &virtual_scan_fine);

  // Downsample from upper limit:
  // construct the upper-most resolution and downsample from there.
  std::vector< std::vector< double> > virtual_scans(num_virtual_scans);

  for (int i = 0; i < virtual_scan_fine.size(); i++)
  {
    unsigned int k = fmod(i,num_virtual_scans);
    virtual_scans[k].push_back(virtual_scan_fine[i]);
  }

  // Make sure that all virtual scans are equal to the real scan in terms of
  // size
  for (unsigned int i = 0; i < virtual_scans.size(); i++)
    assert(virtual_scans[i].size() == real_scan.size());

  // The real scan's (the original) angle increment
  double ang_inc = 2*M_PI / real_scan.size();
  double mul = 1.0 / num_virtual_scans;

  // Compute the angles and metrics of matching the real scan against each and
  // all virtual scans
  std::vector<double> un_angles;
  std::vector<double> snrs;
  std::vector<double> fahms;
  std::vector<double> pds;

  fmt1Batch(real_scan, virtual_scans, r2rp, c2rp,
    &un_angles, &snrs, &fahms, &pds);

  // Correct the angles returned to get the proper pose from which each
  // virtual scan was taken (needed due to over-sampling the map)
  std::vector<double> angles;
  for (unsigned int a = 0; a < num_virtual_scans; a++)
  {
    double angle_a = -un_angles[a] + a*mul*ang_inc;
    Utils::wrapAngle(&angle_a);

    angles.push_back(angle_a);
  }

  // Select some of all the angles based on criteria enforced by rankFMTOutput
  std::vector<unsigned int> optimal_ids =
    rankFMTOutput(snrs, fahms, pds, 3, magnification_size, 0.00001);

  std::vector<double> cand_angles;
  for (unsigned int i = 0; i < optimal_ids.size(); i++)
  {
    cand_angles.push_back(angles[optimal_ids[i]]);

    rc0->push_back(pds[optimal_ids[i]]);
    rc1->push_back(snrs[optimal_ids[i]] / fahms[optimal_ids[i]]);
  }

#if defined (TIMES)
  printf("%f [Rotation::fmt2]\n", elapsed.count());
#endif

#if defined (PRINTS)
  for (unsigned int i = 0; i < cand_angles.size(); i++)
  {
    printf("cand. poses (%f,%f,%f) [Rotation::fmt2]\n",
      std::get<0>(virtual_pose),
      std::get<1>(virtual_pose),
      std::get<2>(virtual_pose)+cand_angles[i]);
  }
#endif

  return cand_angles;
}


/*******************************************************************************
*/
void Rotation::fmt1Batch(
  const std::vector< double >& real_scan,
  const std::vector< std::vector< double > >& virtual_scans,
  const fftw_plan& r2rp, const fftw_plan& c2rp,
  std::vector<double>* angles,
  std::vector<double>* snrs,
  std::vector<double>* fahms,
  std::vector<double>* pds)
{
  std::vector< std::vector<double> > q_0_v;
  std::vector< unsigned int > q_0_max_id_v;
  fmt0Batch(real_scan, virtual_scans, r2rp, c2rp, &q_0_v, &q_0_max_id_v);

  // Calculate PD --------------------------------------------------------------
  std::vector<double> q_ss;
  unsigned int q_ss_max_id;
  fmt0AutoSequential(real_scan, &q_ss, &q_ss_max_id);

  std::vector< std::vector<double> > q_rr_v;
  std::vector< unsigned int > q_rr_max_id_v;
  fmt0AutoBatch(virtual_scans, r2rp, c2rp, &q_rr_v, &q_rr_max_id_v);

  for (unsigned int i = 0; i < virtual_scans.size(); i++)
  {
    // Calculate angle ---------------------------------------------------------
    double angle = static_cast<double>(
      (real_scan.size()-q_0_max_id_v[i]))/(real_scan.size())*2*M_PI;
    Utils::wrapAngle(&angle);

    angles->push_back(angle);

    // Calculate pd ------------------------------------------------------------
    double pd = 2*q_0_v[i][q_0_max_id_v[i]]
      / (q_ss[q_ss_max_id] + q_rr_v[i][q_rr_max_id_v[i]]);
    pds->push_back(pd);

    // Calculate SNR -----------------------------------------------------------
    std::vector<double> q_0_background = q_0_v[i];
    q_0_background.erase(q_0_background.begin() + q_0_max_id_v[i]);

    std::pair<double,double> q_0_mmnts = Utils::vectorStatistics(q_0_background);

    double snr =
      fabs((q_0_v[i][q_0_max_id_v[i]] - q_0_mmnts.first)) / q_0_mmnts.second;
    snrs->push_back(snr);

    // Calculate FAHM ----------------------------------------------------------
    unsigned int count = 0;
    for (unsigned int f = 0; f < q_0_v[i].size(); f++)
    {
      if (q_0_v[i][f] >= 0.5 * q_0_v[i][q_0_max_id_v[i]])
        count++;
    }

    double fahm = static_cast<double>(count) / q_0_v[i].size();
    fahms->push_back(fahm);
  }
}


/*******************************************************************************
*/
void Rotation::fmt0Batch(
  const std::vector<double>& real_scan,
  const std::vector< std::vector<double> > & virtual_scans,
  const fftw_plan& r2rp, const fftw_plan& c2rp,
  std::vector< std::vector<double> >* q_0_v,
  std::vector<unsigned int>* q_0_max_id_v)
{
  assert(virtual_scans.size() > 0);
  assert(real_scan.size() == virtual_scans[0].size());

  // Find fft of real scan
  std::vector<double> fft_real = DFTUtils::dft(real_scan);

  // Find fft of virtual scan
  std::vector< std::vector<double> > fft_virtuals =
    DFTUtils::dftBatch(virtual_scans, r2rp);

  // fft_real is in halfcomplex format; fft_real_coeffs is in normal format
  // (you get the full complex transform)
  std::vector< std::pair<double, double> > fft_real_coeffs =
    DFTUtils::getDFTCoefficientsPairs(fft_real);

  // Find conjugates of real coefficients
  std::vector< std::pair<double, double> > fft_real_coeffs_conj =
    Utils::conjugate(fft_real_coeffs);

  std::vector< std::vector< std::pair<double, double> > > Q_0_v;
  for (unsigned int i = 0; i < virtual_scans.size(); i++)
  {
    std::vector< std::pair<double, double> > fft_virtual_coeffs =
      DFTUtils::getDFTCoefficientsPairs(fft_virtuals[i]);

    // The numerator of Q_0
    std::vector< std::pair<double, double> > numerator =
      Utils::innerProductComplex(fft_real_coeffs_conj, fft_virtual_coeffs);

    /*
    // The denominator of Q_0
    double denominator =
    Utils::norm2(fft_real_coeffs) * Utils::norm2(fft_virtual_coeffs);

    for (int i = 0; i < numerator.size(); i++)
    {
    numerator[i].first /= denominator;
    numerator[i].second /= denominator;
    }
    */

    Q_0_v.push_back(numerator);
  }

  *q_0_v = DFTUtils::idftBatch(Q_0_v, c2rp);

  for (unsigned int i = 0; i < q_0_v->size(); i++)
  {
    unsigned int q_0_max_id =
      std::max_element(q_0_v->at(i).begin(), q_0_v->at(i).end())
      - q_0_v->at(i).begin();

    q_0_max_id_v->push_back(q_0_max_id);
  }
}


/*******************************************************************************
*/
void Rotation::fmt0AutoBatch(
  const std::vector< std::vector<double> > & virtual_scans,
  const fftw_plan& r2rp, const fftw_plan& c2rp,
  std::vector< std::vector<double> >* q_0_v,
  std::vector<unsigned int>* q_0_max_id_v)
{
  assert(virtual_scans.size() > 0);

  // Find fft of virtual scan
  std::vector< std::vector<double> > fft_virtuals =
    DFTUtils::dftBatch(virtual_scans, r2rp);

  std::vector< std::vector< std::pair<double, double> > > Q_0_v;

  for (unsigned int i = 0; i < virtual_scans.size(); i++)
  {
    // Virtual scan dft coefficients
    std::vector< std::pair<double, double> > fft_virtual_coeffs =
      DFTUtils::getDFTCoefficientsPairs(fft_virtuals[i]);

    // Virtual scan dft coefficients conjugates
    std::vector< std::pair<double, double> > fft_virtual_coeffs_conj =
      Utils::conjugate(fft_virtual_coeffs);

    // The numerator of Q_0
    std::vector< std::pair<double, double> > numerator =
      Utils::innerProductComplex(fft_virtual_coeffs_conj, fft_virtual_coeffs);

    Q_0_v.push_back(numerator);
  }

  *q_0_v = DFTUtils::idftBatch(Q_0_v, c2rp);

  for (unsigned int i = 0; i < q_0_v->size(); i++)
  {
    unsigned int q_0_max_id =
      std::max_element(q_0_v->at(i).begin(), q_0_v->at(i).end())
      -q_0_v->at(i).begin();

    q_0_max_id_v->push_back(q_0_max_id);
  }
}


/*******************************************************************************
*/
bool Rotation::fromEllipse2(
  const std::vector< double >& real_scan,
  std::vector<double> real_ellipse_coefficients,
  const std::tuple<double,double,double>& current_pose,
  const std::vector< std::pair<double,double> >& map,
  const unsigned int& num_rays,
  std::tuple<double,double,double>* result_pose)
{
#if defined (TIMES)
  std::chrono::high_resolution_clock::time_point start =
    std::chrono::high_resolution_clock::now();
#endif

  // For projection/abstraction purposes
  std::tuple<double,double,double> zero_pose;
  std::get<0>(zero_pose) = 0.0;
  std::get<1>(zero_pose) = 0.0;
  std::get<2>(zero_pose) = 0.0;

  std::tuple<double,double,double> backup_pose = current_pose;
  *result_pose = current_pose;

#if defined (PRINTS)
  printf("input pose  (%f,%f,%f) [from ellipse2]\n",
    std::get<0>(*result_pose), std::get<1>(*result_pose),
    std::get<2>(*result_pose));
#endif

  std::vector<double> virtual_ellipse_coefficients;
  for (int i = 0; i < 10; i++)
  {
    std::vector< std::pair<double,double> >
      virtual_scan_intersections_after_rot =
      X::find(*result_pose, map, num_rays);

    // We will draw the above points around zero, so first find the virtual scan
    std::vector<double> virtual_scan_it;
    Utils::points2scan(virtual_scan_intersections_after_rot, *result_pose,
      &virtual_scan_it);

    //Dump::scan(real_scan, zero_pose, virtual_scan_it, zero_pose,
    //base_path_ + "/../matlab/scan_dump.m");

    // The are the points with respect to zero
    std::vector< std::pair<double,double> > virtual_scan_points;
    Utils::scan2points(virtual_scan_it, zero_pose, &virtual_scan_points);



    // Find the virtual points bounding ellipse
    Find::boundingEllipse(virtual_scan_points, &virtual_ellipse_coefficients);

    /*
       printf("virtual_coefficients:\n");
       printf("v0 = %f;\n", virtual_ellipse_coefficients[0]);
       printf("v1 = %f;\n", virtual_ellipse_coefficients[1]);
       printf("v2 = %f;\n", virtual_ellipse_coefficients[2]);
       printf("v3 = %f;\n", virtual_ellipse_coefficients[3]);
       printf("v4 = %f;\n", virtual_ellipse_coefficients[4]);
       printf("v5 = %f;\n", virtual_ellipse_coefficients[5]);
       */

    double real_t = Find::ellipseAngle(real_ellipse_coefficients);
    double virtual_t = Find::ellipseAngle(virtual_ellipse_coefficients);

    double t = real_t - virtual_t;
    Utils::wrapAngle(&t);


    while (fabs(t) > M_PI/4)
    {
      t += M_PI/2;
      Utils::wrapAngle(&t);
    }

    std::get<2>(*result_pose) -= t;
    Utils::wrapAngle(&std::get<2>(*result_pose));
  }

  bool ret = true;
  if (fabs(std::get<2>(*result_pose) - std::get<2>(backup_pose)) > M_PI/4)
  {
    *result_pose = backup_pose;
    ret = false;
  }

#if defined (TIMES)
  std::chrono::high_resolution_clock::time_point end =
    std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(end-start);

  printf("%f [findRotationFromEllipse2]\n", elapsed.count());
#endif

#if defined (PRINTS)
  printf("output pose (%f,%f,%f) [from ellipse2]\n",
    std::get<0>(*result_pose), std::get<1>(*result_pose),
    std::get<2>(*result_pose));
#endif

  return ret;
}


/*******************************************************************************
*/
std::vector<double> Rotation::ku2Sequential(
  const std::vector< double >& real_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  const unsigned int& magnification_size,
  std::vector<double>* rc0, std::vector<double>* rc1,
  std::chrono::duration<double>* intersections_time)
{
#if defined (PRINTS)
  printf("input pose  (%f,%f,%f) [Rotation::ku2Sequential]\n",
    std::get<0>(virtual_pose),
    std::get<1>(virtual_pose),
    std::get<2>(virtual_pose));
#endif

  rc0->clear();
  rc1->clear();

  std::tuple<double,double,double> zero_pose;
  std::get<0>(zero_pose) = 0.0;
  std::get<1>(zero_pose) = 0.0;
  std::get<2>(zero_pose) = 0.0;


  unsigned int num_virtual_scans = pow(2,magnification_size);
  int virtual_scan_size_max = num_virtual_scans * real_scan.size();

  // Measure the time to find intersections
  std::chrono::high_resolution_clock::time_point int_start =
    std::chrono::high_resolution_clock::now();

  std::vector< std::pair<double,double> > virtual_scan_points =
    X::find(virtual_pose, map, virtual_scan_size_max);

  std::chrono::high_resolution_clock::time_point int_end =
    std::chrono::high_resolution_clock::now();
  *intersections_time =
    std::chrono::duration_cast< std::chrono::duration<double> >(int_end-int_start);

  std::vector<double> virtual_scan_fine;
  Utils::points2scan(virtual_scan_points, virtual_pose, &virtual_scan_fine);

  // Downsample from upper limit:
  // construct the upper-most resolution and downsample from there.
  std::vector< std::vector< double> > virtual_scans(num_virtual_scans);

  for (int i = 0; i < virtual_scan_fine.size(); i++)
  {
    unsigned int k = fmod(i,num_virtual_scans);
    virtual_scans[k].push_back(virtual_scan_fine[i]);
  }

  // Make sure that all virtual scans are equal to the real scan in terms of
  // size
  for (unsigned int i = 0; i < virtual_scans.size(); i++)
    assert(virtual_scans[i].size() == real_scan.size());

  // The real scan's (the original) angle increment
  double ang_inc = 2*M_PI / real_scan.size();
  double mul = 1.0 / num_virtual_scans;


  std::vector<double> orientations;
  std::vector<double> snrs;
  std::vector<double> fahms;
  std::vector<double> pds;
  std::vector<double> rot_criteria;

  std::vector< std::pair<double,double> > real_scan_points;
  Utils::scan2points(real_scan, zero_pose, &real_scan_points);

  for (unsigned int a = 0; a < num_virtual_scans; a++)
  {
    std::vector< std::pair<double,double> > virtual_scan_points_a;
    Utils::scan2points(virtual_scans[a], zero_pose, &virtual_scan_points_a);

    double angle = 0.0;
    double snr = 0.0;
    double fahm = 0.0;
    double pd = 0.0;

    ku1Sequential(
      real_scan_points, virtual_scan_points_a, &angle, &snr, &fahm, &pd);

    double ornt_a = -angle + a*mul*ang_inc;
    Utils::wrapAngle(&ornt_a);

    orientations.push_back(ornt_a);
    snrs.push_back(snr);
    fahms.push_back(fahm);
    pds.push_back(pd);

#if defined (DEBUG)
    printf("a = %u\n", a);
    printf("angle to out = %f\n", std::get<2>(virtual_pose) + ornt_a);
    printf("snr = %.10f\n", snr);
    printf("fahm = %f\n", fahm);
    printf("pd = %.20f\n", pd);
#endif
  }

  // Select some of all the angles based on criteria enforced by rankKUOutput
  std::vector<unsigned int> optimal_ids =
    rankKUOutput(snrs, fahms, pds, 3, magnification_size, 0.00001);

  std::vector<double> angles;
  for (unsigned int i = 0; i < optimal_ids.size(); i++)
  {
    double angle = orientations[optimal_ids[i]];
    Utils::wrapAngle(&angle);
    angles.push_back(angle);

    rc0->push_back(pds[optimal_ids[i]]);
    rc1->push_back(snrs[optimal_ids[i]] / fahms[optimal_ids[i]]);
  }

#if defined (PRINTS)
  for (unsigned int i = 0; i < angles.size(); i++)
  {
    printf("cand. poses (%f,%f,%f) [Rotation::ku2Sequential]\n",
      std::get<0>(virtual_pose),
      std::get<1>(virtual_pose),
      std::get<2>(virtual_pose)+angles[i]);
  }
#endif

  return angles;
}


/*******************************************************************************
*/
void Rotation::ku1Sequential(
  const std::vector< std::pair<double,double> >& real_scan_points,
  const std::vector< std::pair<double,double> >& virtual_scan_points,
  double* angle, double* snr, double* fahm, double* pd)
{
  std::vector<double> traces;
  unsigned int traces_max_id;
  ku0Sequential(real_scan_points, virtual_scan_points, &traces, &traces_max_id);

  // Calculate angle -----------------------------------------------------------
  int rot_id = traces_max_id;
  *angle = static_cast<double>(
    (real_scan_points.size()-rot_id))/(real_scan_points.size())*2*M_PI;

  Utils::wrapAngle(angle);

  // Calculate SNR -------------------------------------------------------------
  std::vector<double> traces_background = traces;
  traces_background.erase(traces_background.begin() + traces_max_id);

  std::pair<double,double> traces_mmnts =
    Utils::vectorStatistics(traces_background);

  *snr = fabs((traces[traces_max_id] - traces_mmnts.first)) / traces_mmnts.second;

  // Calculate FAHM ------------------------------------------------------------
  unsigned int count = 0;
  for (unsigned int i = 0; i < traces.size(); i++)
  {
    if (traces[i] >= 0.5 * traces[traces_max_id])
      count++;
  }

  *fahm = static_cast<double>(count) / traces.size();

  // Calculate PD --------------------------------------------------------------
  std::vector<double> traces_ss;
  unsigned int traces_ss_max_id;
  ku0AutoSequential(real_scan_points, &traces_ss, &traces_ss_max_id);

  std::vector<double> traces_rr;
  unsigned int traces_rr_max_id;
  ku0AutoSequential(virtual_scan_points, &traces_rr, &traces_rr_max_id);

  *pd = 2*traces[traces_max_id] /
    (traces_ss[traces_ss_max_id] + traces_rr[traces_rr_max_id]);
}


/*******************************************************************************
*/
void Rotation::ku0Sequential(
  const std::vector< std::pair<double,double> >& real_scan_points,
  const std::vector< std::pair<double,double> >& virtual_scan_points_in,
  std::vector<double>* traces, unsigned int* traces_max_id)
{
  traces->clear();

  Eigen::MatrixXd R(2, real_scan_points.size());
  for (int i = 0; i < real_scan_points.size(); i++)
  {
    R(0,i) = real_scan_points[i].first;
    R(1,i) = real_scan_points[i].second;
  }

  std::vector< std::pair<double,double> > virtual_scan_points =
    virtual_scan_points_in;
  //std::reverse(virtual_scan_points.begin(), virtual_scan_points.end());

  std::vector< std::pair<double,double> > virtual_scan_points_0 =
    virtual_scan_points;



  std::vector<double> A11;
  std::vector<double> A12;
  std::vector<double> A21;
  std::vector<double> A22;

  std::vector<Eigen::MatrixXd> As;
  Eigen::MatrixXd V(2, virtual_scan_points.size());

  for (int i = 0; i < virtual_scan_points.size(); i++)
  {
    std::rotate(virtual_scan_points.begin(),
      virtual_scan_points.begin()+i,
      virtual_scan_points.end());

    // TODO implement shift of columns for faster execution
    for (int j = 0; j < virtual_scan_points.size(); j++)
    {
      V(0,j) = virtual_scan_points[j].first;
      V(1,j) = virtual_scan_points[j].second;
    }

    Eigen::Matrix2d A = V * R.transpose();

    A11.push_back(A(0,0));
    A12.push_back(A(0,1));
    A21.push_back(A(1,0));
    A22.push_back(A(1,1));


    virtual_scan_points = virtual_scan_points_0;
  }




  for (int i = 0; i < A11.size(); i++)
  {
    // Only consider angular deviations +/- 45deg TODO SPEEDUP
    //if (i > A11.size()/4 && i < 3*A11.size()/4)
    //continue;

    Eigen::Matrix2d A(2,2);
    A(0,0) = A11[i];
    A(0,1) = A12[i];
    A(1,0) = A21[i];
    A(1,1) = A22[i];

    Eigen::JacobiSVD<Eigen::Matrix2d> svd_of_A(A,
      Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::Matrix2d _S_;
    if (svd_of_A.matrixU().determinant() * svd_of_A.matrixV().determinant() > 0)
    {
      _S_(0,0) = 1;
      _S_(1,0) = 0;
      _S_(0,1) = 0;
      _S_(1,1) = 1;
    }
    else
    {
      _S_(0,0) = 1;
      _S_(1,0) = 0;
      _S_(0,1) = 0;
      _S_(1,1) = -1;
    }

    Eigen::Matrix2d R =
      svd_of_A.matrixU() * _S_ * svd_of_A.matrixV().transpose();

    Eigen::Matrix2d RAt = R * A.transpose();

    traces->push_back(RAt.trace());
  }

  // Reverse the reversion
  //std::reverse(traces->begin(), traces->end());

  *traces_max_id =
    std::max_element(traces->begin(), traces->end()) - traces->begin();
}


/*******************************************************************************
*/
void Rotation::ku0AutoSequential(
  const std::vector< std::pair<double,double> >& real_scan_points,
  std::vector<double>* traces, unsigned int* traces_max_id)
{
  traces->clear();

  Eigen::MatrixXd R(2, real_scan_points.size());
  for (int i = 0; i < real_scan_points.size(); i++)
  {
    R(0,i) = real_scan_points[i].first;
    R(1,i) = real_scan_points[i].second;
  }


  Eigen::Matrix2d A(2,2);
  A = R * R.transpose();


  Eigen::JacobiSVD<Eigen::Matrix2d> svd_of_A(A,
    Eigen::ComputeFullU | Eigen::ComputeFullV);

  Eigen::Matrix2d _S_;
  if (svd_of_A.matrixU().determinant() * svd_of_A.matrixV().determinant() > 0)
  {
    _S_(0,0) = 1;
    _S_(1,0) = 0;
    _S_(0,1) = 0;
    _S_(1,1) = 1;
  }
  else
  {
    _S_(0,0) = 1;
    _S_(1,0) = 0;
    _S_(0,1) = 0;
    _S_(1,1) = -1;
  }

  Eigen::Matrix2d Rot =
    svd_of_A.matrixU() * _S_ * svd_of_A.matrixV().transpose();

  Eigen::Matrix2d RAt = Rot * A.transpose();

  traces->push_back(RAt.trace());

  // Reverse the reversion
  //std::reverse(traces->begin(), traces->end());

  *traces_max_id =
    std::max_element(traces->begin(), traces->end()) - traces->begin();
}


/*******************************************************************************
*/
std::vector<unsigned int> Rotation::rankDBHOutput(
  const std::vector<double>& snr,
  const std::vector<double>& fahm,
  const std::vector<double>& pd,
  const unsigned int& method,
  const unsigned int& magnification_size,
  const double& pd_threshold)
{
  assert (snr.size() == fahm.size());
  assert (fahm.size() == pd.size());
  assert (pd_threshold >= 0);
  assert (method >= 0 && method <= 3);

  // Return the indices of those angles for which criteria are near
  // the maximum criterion
  std::vector<unsigned int> best_ids;

  // Simply the one please
  if (method == 0)
  {
    // What are the criteria for ranking angles?
    std::vector<double> criteria = pd;

    // Identify maximum criterion
    double max_c = *std::max_element(criteria.begin(), criteria.end());

    best_ids.push_back(
      std::max_element(criteria.begin(), criteria.end()) -criteria.begin());

  }

  // The one + those within pd_threshold around it
  if (method == 1)
  {
    // What are the criteria for ranking angles?
    std::vector<double> criteria = pd;

    // Identify maximum criterion
    double max_c = *std::max_element(criteria.begin(), criteria.end());

    for (unsigned int i = 0; i < criteria.size(); i++)
    {
      if (fabs(criteria[i]-max_c) <= pd_threshold)
        best_ids.push_back(i);
    }
  }

  // The one + those within (max critetia - min crtieria)/2
  if (method == 2)
  {
    // What are the criteria for ranking angles?
    std::vector<double> criteria = pd;

    // Identify maximum criterion
    double max_c = *std::max_element(criteria.begin(), criteria.end());
    double min_c = *std::min_element(criteria.begin(), criteria.end());


    for (unsigned int i = 0; i < criteria.size(); i++)
    {
      if (fabs(criteria[i]-max_c) <= (max_c - min_c)/2)
        best_ids.push_back(i);
    }
  }

  // Pick `pick_num_surr` around max criterion every time
  std::set<unsigned int> best_ids_set;
  if (method == 3)
  {
    // What are the criteria for ranking angles?
    std::vector<double> criteria = pd;


    // Identify maximum criterion
    int max_c_idx =
      std::max_element(criteria.begin(), criteria.end()) - criteria.begin();
    double max_c = criteria[max_c_idx];

#if defined (DEBUG)
    printf("best id = %d\n", max_c_idx);
#endif

    int vendalia_method = 1;

    int pick_num_surr = 0;
    if (vendalia_method == 0)
    {
      pick_num_surr = pow(2,magnification_size) / pow(2,3);
      if (pick_num_surr == 0)
        pick_num_surr = 1;
    }
    if (vendalia_method == 1)
      pick_num_surr = pow(2,7) / pow(2,magnification_size);
    if (vendalia_method == 2)
      pick_num_surr = 4;


    for (int i =  -pick_num_surr + max_c_idx;
      i <= +pick_num_surr + max_c_idx; i++)
    {
      int k = i;

      while(k < 0)
        k += criteria.size();

      while (k > criteria.size())
        k -= criteria.size();

      if (k == criteria.size())
        k = 0;

#if defined (DEBUG)
      printf("k = %d\n", k);
#endif
      best_ids_set.insert(k);
    }

    /*
       for (unsigned int i = 0; i < criteria.size(); i++)
       {
       if (fabs(criteria[i]-max_c) <= pd_threshold)
       best_ids_set.insert(i);
       }
       */

    for (std::set<unsigned int>::iterator it = best_ids_set.begin();
      it != best_ids_set.end(); it++) best_ids.push_back(*it);
  }

#if defined (DEBUG)
  printf("BEST IDS = [");
  for (unsigned int i = 0; i < best_ids.size(); i++)
    printf("%u ", best_ids[i]);

  printf("]\n");
#endif

  return best_ids;
}


/*******************************************************************************
*/
std::vector<unsigned int> Rotation::rankFMTOutput(
  const std::vector<double>& snr,
  const std::vector<double>& fahm,
  const std::vector<double>& pd,
  const unsigned int& method,
  const unsigned int& magnification_size,
  const double& pd_threshold)
{
  assert (snr.size() == fahm.size());
  assert (fahm.size() == pd.size());
  assert (pd_threshold >= 0);
  assert (method >= 0 && method <= 3);

  // Return the indices of those angles for which criteria are near
  // the maximum criterion
  std::vector<unsigned int> best_ids;

  // Simply the one please
  if (method == 0)
  {
    // What are the criteria for ranking angles?
    std::vector<double> criteria = pd;

    // Identify maximum criterion
    double max_c = *std::max_element(criteria.begin(), criteria.end());

    best_ids.push_back(
      std::max_element(criteria.begin(), criteria.end()) -criteria.begin());

  }

  // The one + those within pd_threshold around it
  if (method == 1)
  {
    // What are the criteria for ranking angles?
    std::vector<double> criteria = pd;

    // Identify maximum criterion
    double max_c = *std::max_element(criteria.begin(), criteria.end());

    for (unsigned int i = 0; i < criteria.size(); i++)
    {
      if (fabs(criteria[i]-max_c) <= pd_threshold)
        best_ids.push_back(i);
    }
  }

  // The one + those within (max critetia - min crtieria)/2
  if (method == 2)
  {
    // What are the criteria for ranking angles?
    std::vector<double> criteria = pd;

    // Identify maximum criterion
    double max_c = *std::max_element(criteria.begin(), criteria.end());
    double min_c = *std::min_element(criteria.begin(), criteria.end());


    for (unsigned int i = 0; i < criteria.size(); i++)
    {
      if (fabs(criteria[i]-max_c) <= (max_c - min_c)/2)
        best_ids.push_back(i);
    }
  }

  // Pick `pick_num_surr` around max criterion every time
  std::set<unsigned int> best_ids_set;
  if (method == 3)
  {
    // What are the criteria for ranking angles?
    std::vector<double> criteria = pd;


    // Identify maximum criterion
    int max_c_idx =
      std::max_element(criteria.begin(), criteria.end()) - criteria.begin();
    double max_c = criteria[max_c_idx];

#if defined (DEBUG)
    printf("best id = %d\n", max_c_idx);
#endif

    int vendalia_method = 1;

    int pick_num_surr = 0;
    if (vendalia_method == 0)
    {
      pick_num_surr = pow(2,magnification_size) / pow(2,3);
      if (pick_num_surr == 0)
        pick_num_surr = 1;
    }
    if (vendalia_method == 1)
      pick_num_surr = pow(2,7) / pow(2,magnification_size);
    if (vendalia_method == 2)
      pick_num_surr = 4;


    for (int i =  -pick_num_surr + max_c_idx;
      i <= +pick_num_surr + max_c_idx; i++)
    {
      int k = i;

      while(k < 0)
        k += criteria.size();

      while (k > criteria.size())
        k -= criteria.size();

      if (k == criteria.size())
        k = 0;

#if defined (DEBUG)
      printf("k = %d\n", k);
#endif
      best_ids_set.insert(k);
    }

    /*
       for (unsigned int i = 0; i < criteria.size(); i++)
       {
       if (fabs(criteria[i]-max_c) <= pd_threshold)
       best_ids_set.insert(i);
       }
       */

    for (std::set<unsigned int>::iterator it = best_ids_set.begin();
      it != best_ids_set.end(); it++) best_ids.push_back(*it);
  }

#if defined (DEBUG)
  printf("BEST IDS = [");
  for (unsigned int i = 0; i < best_ids.size(); i++)
    printf("%u ", best_ids[i]);

  printf("]\n");
#endif

  return best_ids;
}


/*******************************************************************************
*/
std::vector<unsigned int> Rotation::rankKUOutput(
  const std::vector<double>& snr,
  const std::vector<double>& fahm,
  const std::vector<double>& pd,
  const unsigned int& method,
  const unsigned int& magnification_size,
  const double& pd_threshold)
{
  assert (snr.size() == fahm.size());
  assert (fahm.size() == pd.size());
  assert (pd_threshold >= 0);
  assert (method >= 0 && method <= 3);

  // Return the indices of those angles for which criteria are near
  // the maximum criterion
  std::vector<unsigned int> best_ids;

  // Simply the one please
  if (method == 0)
  {
    // What are the criteria for ranking angles?
    std::vector<double> criteria = pd;

    // Identify maximum criterion
    double max_c = *std::max_element(criteria.begin(), criteria.end());

    best_ids.push_back(
      std::max_element(criteria.begin(), criteria.end()) -criteria.begin());

  }

  // The one + those within pd_threshold around it
  if (method == 1)
  {
    // What are the criteria for ranking angles?
    std::vector<double> criteria = pd;

    // Identify maximum criterion
    double max_c = *std::max_element(criteria.begin(), criteria.end());

    for (unsigned int i = 0; i < criteria.size(); i++)
    {
      if (fabs(criteria[i]-max_c) <= pd_threshold)
        best_ids.push_back(i);
    }
  }

  // The one + those within (max critetia - min crtieria)/2
  if (method == 2)
  {
    // What are the criteria for ranking angles?
    std::vector<double> criteria = pd;

    // Identify maximum criterion
    double max_c = *std::max_element(criteria.begin(), criteria.end());
    double min_c = *std::min_element(criteria.begin(), criteria.end());


    for (unsigned int i = 0; i < criteria.size(); i++)
    {
      if (fabs(criteria[i]-max_c) <= (max_c - min_c)/2)
        best_ids.push_back(i);
    }
  }

  // Pick `pick_num_surr` around max criterion every time
  std::set<unsigned int> best_ids_set;
  if (method == 3)
  {
    // What are the criteria for ranking angles?
    std::vector<double> criteria = pd;


    // Identify maximum criterion
    int max_c_idx =
      std::max_element(criteria.begin(), criteria.end()) - criteria.begin();
    double max_c = criteria[max_c_idx];

#if defined (DEBUG)
    printf("best id = %d\n", max_c_idx);
#endif

    int vendalia_method = 1;

    int pick_num_surr = 0;
    if (vendalia_method == 0)
    {
      pick_num_surr = pow(2,magnification_size) / pow(2,3);
      if (pick_num_surr == 0)
        pick_num_surr = 1;
    }
    if (vendalia_method == 1)
      pick_num_surr = pow(2,7) / pow(2,magnification_size);
    if (vendalia_method == 2)
      pick_num_surr = 4;


    for (int i =  -pick_num_surr + max_c_idx;
      i <= +pick_num_surr + max_c_idx; i++)
    {
      int k = i;

      while(k < 0)
        k += criteria.size();

      while (k > criteria.size())
        k -= criteria.size();

      if (k == criteria.size())
        k = 0;

#if defined (DEBUG)
      printf("k = %d\n", k);
#endif
      best_ids_set.insert(k);
    }

    /*
       for (unsigned int i = 0; i < criteria.size(); i++)
       {
       if (fabs(criteria[i]-max_c) <= pd_threshold)
       best_ids_set.insert(i);
       }
       */

    for (std::set<unsigned int>::iterator it = best_ids_set.begin();
      it != best_ids_set.end(); it++) best_ids.push_back(*it);
  }

#if defined (DEBUG)
  printf("BEST IDS = [");
  for (unsigned int i = 0; i < best_ids.size(); i++)
    printf("%u ", best_ids[i]);

  printf("]\n");
#endif

  return best_ids;
}


/*******************************************************************************
*/
std::vector<double> Rotation::skg(
  const std::vector< double >& real_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  const unsigned int& magnification_size,
  const fftw_plan& r2rp,
  std::vector<double>* rc0, std::vector<double>* rc1,
  std::chrono::duration<double>* intersections_time)
{
#if defined (PRINTS)
  printf("input pose  (%f,%f,%f) [Rotation::skg]\n",
    std::get<0>(virtual_pose),
    std::get<1>(virtual_pose),
    std::get<2>(virtual_pose));
#endif

  rc0->clear();
  rc1->clear();

  std::tuple<double,double,double> zero_pose;
  std::get<0>(zero_pose) = 0.0;
  std::get<1>(zero_pose) = 0.0;
  std::get<2>(zero_pose) = 0.0;


  unsigned int num_virtual_scans = pow(2,magnification_size);
  int virtual_scan_size_max = num_virtual_scans * real_scan.size();

  // Measure the time to find intersections
  std::chrono::high_resolution_clock::time_point int_start =
    std::chrono::high_resolution_clock::now();

  std::vector< std::pair<double,double> > virtual_scan_points =
    X::find(virtual_pose, map, virtual_scan_size_max);

  std::chrono::high_resolution_clock::time_point int_end =
    std::chrono::high_resolution_clock::now();
  *intersections_time =
    std::chrono::duration_cast< std::chrono::duration<double> >(int_end-int_start);

  std::vector<double> virtual_scan_fine;
  Utils::points2scan(virtual_scan_points, virtual_pose, &virtual_scan_fine);

  // Downsample from upper limit:
  // construct the upper-most resolution and downsample from there.
  std::vector< std::vector< double> > virtual_scans(num_virtual_scans);

  for (int i = 0; i < virtual_scan_fine.size(); i++)
  {
    unsigned int k = fmod(i,num_virtual_scans);
    virtual_scans[k].push_back(virtual_scan_fine[i]);
  }

  // Make sure that all virtual scans are equal to the real scan in terms of
  // size
  for (unsigned int i = 0; i < virtual_scans.size(); i++)
    assert(virtual_scans[i].size() == real_scan.size());

  // The real scan's (the original) angle increment
  double ang_inc = 2*M_PI / real_scan.size();
  double mul = 1.0 / num_virtual_scans;


  std::vector<double> orientations;
  std::vector<double> snrs;
  std::vector<double> fahms;
  std::vector<double> pds;

  for (unsigned int a = 0; a < num_virtual_scans; a++)
  {
    double angle = 0.0;
    double snr = 1.0;
    double fahm = 1.0;
    double pd = 1.0;

    skg0(real_scan, virtual_scans[a], r2rp, &angle);

    double ornt_a = angle + a*mul*ang_inc;
    Utils::wrapAngle(&ornt_a);

    orientations.push_back(ornt_a);

    rc0->push_back(1.0);
    rc1->push_back(1.0);

#if defined (DEBUG)
    printf("a = %u\n", a);
    printf("angle to out = %f\n", std::get<2>(virtual_pose) + ornt_a);
#endif
  }

#if defined (TIMES)
  printf("%f [Rotation::skg]\n", elapsed.count());
#endif

#if defined (PRINTS)
  for (unsigned int i = 0; i < orientations.size(); i++)
  {
    printf("cand. poses (%f,%f,%f) [Rotation::skg]\n",
      std::get<0>(virtual_pose),
      std::get<1>(virtual_pose),
      std::get<2>(virtual_pose)+orientations[i]);
  }
#endif

  return orientations;
}


/*******************************************************************************
*/
void Rotation::skg0(
  const std::vector< double >& real_scan,
  const std::vector< double >& virtual_scan,
  const fftw_plan& r2rp,
  double* angle)
{
  // Compute R1, V1
  std::vector<double> R1 = DFTUtils::getFirstDFTCoefficient(real_scan, r2rp);
  std::vector<double> V1 = DFTUtils::getFirstDFTCoefficient(virtual_scan, r2rp);

  // Rotate
  double mov_t = atan2(R1[1],R1[0]) - atan2(V1[1],V1[0]);

  //if (mov_t < -M_PI/2) mov_t += M_PI;
  //if (mov_t > +M_PI/2) mov_t -= M_PI;

  *angle = mov_t;
}
