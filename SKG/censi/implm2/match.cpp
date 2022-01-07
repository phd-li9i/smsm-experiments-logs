#include "match.h"


/*******************************************************************************
*/
bool Match::canGiveNoMore(
  const std::vector<double>& xs,
  const std::vector<double>& ys,
  const std::vector<double>& ts,
  const double& xy_eps,
  const double& t_eps)
{
  assert(xs.size() == ys.size());

  unsigned int sz = xs.size();
  bool xy_converged = false;
  bool t_converged = false;

  if (sz < 2)
    return false;
  else
  {
    for (unsigned int i = 2; i < sz; i++)
    {
      if (fabs(ts[sz-1] - ts[sz-i]) < t_eps)
        t_converged = true;

      if (fabs(xs[sz-1] - xs[sz-i]) < xy_eps &&
          fabs(ys[sz-1] - ys[sz-i]) < xy_eps)
        xy_converged = true;

      if (xy_converged && t_converged)
        return true;
    }

    return false;
  }
}


/*******************************************************************************
*/
void Match::csm(
  const std::vector< double >& real_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  sm_params* input_, sm_result* output_,
  const input_params& ip, output_params* op,
  std::tuple<double,double,double>* result_pose)
{
  // Start the clock
  std::chrono::high_resolution_clock::time_point start =
    std::chrono::high_resolution_clock::now();

  std::tuple<double,double,double> current_pose = virtual_pose;

  // The map scan
  std::vector<double> virtual_scan;

  // The total number of inside iterations
  unsigned int num_iterations = 0;

  // The total correspondence error
  std::vector<double> c_errors;
  std::vector<double> traj_x;
  std::vector<double> traj_y;
  std::vector<double> traj_t;

  LDP real_scan_ldp;
  LDP virtual_scan_ldp;
  // --------------------- ROTATION SEGMENT ----------------------------------
#if (defined TIMES) || (defined LOGS)
  std::chrono::high_resolution_clock::time_point start_rotation_x =
    std::chrono::high_resolution_clock::now();
#endif

  Utils::scanFromPose(current_pose, map, real_scan.size(), &virtual_scan);

  assert(virtual_scan.size() == real_scan.size());

  // Convert scans to LDP
  CS2MSM::convertRealScanToLDP(real_scan, current_pose, real_scan_ldp);
  CS2MSM::convertVirtualScanToLDP(virtual_scan, current_pose, virtual_scan_ldp);

  input_->laser_ref = real_scan_ldp;
  input_->laser_sens = virtual_scan_ldp;

  // DO sm_gpm
  sm_gpm(input_, output_);

  num_iterations += output_->iterations;

  // Update orientation
  double dt = output_->x[2];
  Utils::wrapAngle(&dt);

  std::get<2>(current_pose) -= dt;
  Utils::wrapAngle(&std::get<2>(current_pose));

  ld_free(real_scan_ldp);
  ld_free(virtual_scan_ldp);

#if (defined TIMES) || (defined LOGS)
  std::chrono::high_resolution_clock::time_point end_rotation_x =
    std::chrono::high_resolution_clock::now();

  op->intersections_times += std::chrono::duration_cast<
    std::chrono::duration<double> >(end_rotation_x-start_rotation_x).count();
#endif

  for(int i = 0; i < ip.num_iterations; i++)
  {
    // -------------------- PLICP SEGMENT --------------------------------------
    virtual_scan.clear();

#if (defined TIMES) || (defined LOGS)
    std::chrono::high_resolution_clock::time_point start_translation_x =
      std::chrono::high_resolution_clock::now();
#endif

    Utils::scanFromPose(current_pose, map, real_scan.size(), &virtual_scan);

#if (defined TIMES) || (defined LOGS)
    std::chrono::high_resolution_clock::time_point end_translation_x =
      std::chrono::high_resolution_clock::now();

    op->intersections_times +=std::chrono::duration_cast<
      std::chrono::duration<double> >(end_translation_x-start_translation_x).count();
#endif

    // Convert real scan to LDP
    CS2MSM::convertRealScanToLDP(real_scan, current_pose, real_scan_ldp);

    // Convert virtual scan LDP
    CS2MSM::convertVirtualScanToLDP(virtual_scan, current_pose, virtual_scan_ldp);

    // DO sm_icp
    input_->laser_ref = real_scan_ldp;
    input_->laser_sens = virtual_scan_ldp;
    sm_icp(input_, output_);

    num_iterations += output_->iterations;

    // Update location
    std::get<0>(current_pose) -= output_->x[0];
    std::get<1>(current_pose) -= output_->x[1];

    // Update orientation as well (works better than leaving it to gpm alone)
    dt = output_->x[2];
    Utils::wrapAngle(&dt);
    std::get<2>(current_pose) -= dt;
    Utils::wrapAngle(&std::get<2>(current_pose));

    ld_free(real_scan_ldp);
    ld_free(virtual_scan_ldp);

    c_errors.push_back(output_->error);

    traj_x.push_back(std::get<0>(current_pose));
    traj_y.push_back(std::get<1>(current_pose));
    traj_t.push_back(std::get<2>(current_pose));

    if(!Utils::isPositionInMap(current_pose, map) || output_->valid != 1)
      l2recovery(virtual_pose, map, ip.xy_bound, ip.t_bound, &current_pose);

  }

  // Pick the pose with the least matching error
  int min_idx =
    std::min_element(c_errors.begin(), c_errors.end()) - c_errors.begin();
  std::get<0>(*result_pose) = traj_x[min_idx];
  std::get<1>(*result_pose) = traj_y[min_idx];
  std::get<2>(*result_pose) = traj_t[min_idx];

  // Stop the clock
  std::chrono::high_resolution_clock::time_point end =
    std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(end-start);

  op->exec_time = elapsed.count();

  //*result_pose = current_pose;

  op->rotation_iterations = num_iterations;
  op->translation_iterations = num_iterations;
}


/*******************************************************************************
  The rationale of this function is the same for both matching methods, FMT and
  DBH: Overall there are two discrete stages: one is (A) rotation and the other
  is (B) translation.
  (A): You take some scans from the map, feed those and the real scan to the
  core matching method (fmt2 or dbh2), and in return you get a number of
  best orientation estimates (best in the sense of some metric which
  quantifies the confidence in accurate matching), along with these metrics.
  Although you could have taken the best candidate of these scans, due to the
  position displacement it is actually not the best candidate despite the fact
  that the metric says otherwise. So how do you sift through them? You update
  the current pose estimate with each candidate orientation estimate and then
  perform one step of translation for each of them. In theory the metric of
  summing the pair-wise ray differences of the real scan and the virtual scan
  that corresponds to the newly-moved candidate pose estimate should be lower
  the best alignment was performed.
  (B): So the pose candidate with the lowest pair-wise difference sum wins and
  gets to be translated all the way.

  In the meantime we are storing the best orientation estimate and feeding it
  also as a candidate angle at each round.

*/
void Match::fmtdbh(
  const std::vector< double >& real_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  const std::string& match_method,
  const fftw_plan& r2rp, const fftw_plan& c2rp,
  const input_params& ip, output_params* op,
  std::tuple<double,double,double>* result_pose)
{
  std::chrono::high_resolution_clock::time_point start =
    std::chrono::high_resolution_clock::now();

  *result_pose = virtual_pose;

  // Maximum counter value means a new recovery attempt
  int min_counter = 0;
  int max_counter = 10;
  int counter = min_counter;

  // By a factor of what do you need to over-sample angularly?
  unsigned int min_magnification_size = 2;
  unsigned int max_magnification_size = 5;
  unsigned int current_magnification_size = min_magnification_size;

  // How many times do I attempt recovery?
  unsigned int num_recoveries = 0;
  unsigned int max_recoveries = 1000;

  // These three vectors hold the trajectory for each iteration
  std::vector<double> xs;
  std::vector<double> ys;
  std::vector<double> ts;

  // Two rotation criteria
  std::vector<double> rc0_v;
  std::vector<double> rc1_v;

  // One translation criterion
  std::vector<double> tc_v;

  std::vector<double> dxys;
  std::chrono::duration<double> intersections_time;

  // The best candidate angle found at each iterations is stored and made a
  // candidate each time. Its criterion is its translation criterion after
  // ni-1 translations
  double best_cand_angle = 0.0;
  double best_min_tc = 100000.0;

  // A lock for going overdrive when the rotation criterion is near-excellent
  bool up_lock = false;
  int total_iterations = 0;
  int num_iterations = 0;


  // ROTATION ONLY TEST; (same location) ---------------------------------------
#if defined (TEST_ROTATION_ONLY_DISC) || defined (TEST_ROTATION_ONLY_CONT)
  while (current_magnification_size <= max_magnification_size)
  {
    printf("current_magnification_size = %d ---\n", current_magnification_size);
    printf("counter                    = %d ---\n", counter);

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // ------------------ Rotation correction phase ----------------------------
    std::vector<double> rc0;
    std::vector<double> rc1;
    std::vector<double> dts;

    if (match_method.compare("FMT") == 0)
      dts = Rotation::fmt(real_scan, *result_pose, map,
        current_magnification_size, "batch", r2rp, c2rp,
        &rc0, &rc1, &intersections_time);

    if (match_method.compare("DBH") == 0)
      dts = Rotation::dbh(real_scan, *result_pose, map,
        current_magnification_size, "batch", r2rp, c2rp,
        &rc0, &rc1, &intersections_time);

    unsigned int max_rc0_idx = std::max_element(rc0.begin(), rc0.end())
      - rc0.begin();
    std::get<2>(*result_pose) += dts[max_rc0_idx];
    Utils::wrapAngle(&std::get<2>(*result_pose));

    current_magnification_size++;
  }
  return;
#endif

  // TRANSLATION ONLY TEST; (same orientation) ---------------------------------
#if defined (TEST_TRANSLATION_ONLY)
  int tr_iterations = -1;
  double trans_criterion = 0.0;
  do
  {
    current_magnification_size = max_magnification_size;
    double int_time_trans = 0.0;

    trans_criterion = Translation::tff(real_scan, *result_pose, map, 60, false,
      ip.xy_bound, r2rp, &tr_iterations, &int_time_trans, result_pose);

    if (trans_criterion != -2.0)
      current_magnification_size++;
    else
      while(!Utils::generatePose(virtual_pose, map,
          ip.xy_bound, 0.0, 0.0, result_pose));

  } while (trans_criterion == -2.0);

  return;
#endif


  while (current_magnification_size <= max_magnification_size)
  {
#if defined (DEBUG)
    printf("current_magnification_size = %d ---\n", current_magnification_size);
    printf("counter                    = %d ---\n", counter);
#endif

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // ------------------ Rotation correction phase ----------------------------
    std::vector<double> rc0;
    std::vector<double> rc1;
    std::vector<double> cand_angles;

#if (defined TIMES) || (defined LOGS)
    std::chrono::high_resolution_clock::time_point start_rotation =
      std::chrono::high_resolution_clock::now();
#endif

    if (match_method.compare("FMT") == 0)
      cand_angles = Rotation::fmt(real_scan, *result_pose, map,
        current_magnification_size, "batch", r2rp, c2rp,
        &rc0, &rc1, &intersections_time);

    if (match_method.compare("DBH") == 0)
      cand_angles = Rotation::dbh(real_scan, *result_pose, map,
        current_magnification_size, "batch", r2rp, c2rp,
        &rc0, &rc1, &intersections_time);

    if (match_method.compare("KU") == 0)
      cand_angles = Rotation::ku2Sequential(real_scan, *result_pose, map,
        current_magnification_size,
        &rc0, &rc1, &intersections_time);

#if (defined TIMES) || (defined LOGS)
    std::chrono::high_resolution_clock::time_point end_rotation =
      std::chrono::high_resolution_clock::now();

    op->rotation_times += std::chrono::duration_cast<
      std::chrono::duration<double> >(end_rotation-start_rotation).count();

    op->intersections_times += intersections_time.count();
#endif

    bool ca_exists = false;
    for (unsigned int i = 0; i < cand_angles.size(); i++)
    {
      if (cand_angles[i] == best_cand_angle)
      {
        ca_exists = true;
        break;
      }
    }
    if (!ca_exists)
      cand_angles.push_back(best_cand_angle);

    // ------------------ Candidate angles sifting -----------------------------
    unsigned int min_tc_idx = 0;
    if (cand_angles.size() > 1)
    {
      std::vector<double> tcs_sift;
      for (unsigned int ca = 0; ca < cand_angles.size(); ca++)
      {
        // How many test iterations?
        unsigned int ni = 2;
        int tr_i = 0;

        std::tuple<double,double,double> cand_pose = *result_pose;
        std::get<2>(cand_pose) += cand_angles[ca];

#if (defined TIMES) || (defined LOGS)
        std::chrono::high_resolution_clock::time_point start_translation =
          std::chrono::high_resolution_clock::now();
#endif

        double tc = Translation::tff(real_scan, cand_pose, map,
          ni, ip.xy_bound, false, r2rp, &tr_i, &intersections_time, &cand_pose);

#if (defined TIMES) || (defined LOGS)
        std::chrono::high_resolution_clock::time_point end_translation =
          std::chrono::high_resolution_clock::now();

        op->translation_times += std::chrono::duration_cast<
          std::chrono::duration<double> >(end_translation-start_translation).count();

        op->intersections_times += intersections_time.count();
#endif

#if (defined LOGS)
        op->translation_iterations += tr_i;
#endif

        if (tc == -2.0)
          tcs_sift.push_back(1000000.0);
        else
          tcs_sift.push_back(tc);
      }

      // The index of the angle with the least translation criterion
      min_tc_idx =
        std::min_element(tcs_sift.begin(), tcs_sift.end()) - tcs_sift.begin();


      // Check if the newly-found angle is the angle with the least
      // translation criterion so far
      if (tcs_sift[min_tc_idx] < best_min_tc)
      {
        best_min_tc = tcs_sift[min_tc_idx];
        best_cand_angle = cand_angles[min_tc_idx];
      }
    }

    rc0_v.push_back(rc0[min_tc_idx]);
    rc1_v.push_back(rc1[min_tc_idx]);

    // Update the current orientation estimate with the angle that sports the
    // least translation criterion overall
    std::get<2>(*result_pose) += cand_angles[min_tc_idx];
    Utils::wrapAngle(&std::get<2>(*result_pose));

    // ... and store it
    ts.push_back(std::get<2>(*result_pose));

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // ---------------- Translation correction phase ---------------------------
    num_iterations =
      (current_magnification_size+1)*ip.num_iterations;

    // Place a very heavy burden on the first few iterations:
    // if a shitty angle has been calculated then the position estimate will
    // either fall out of the map, or positional bounds will be exceeded;
    // on the contrary, an angle estimate with low error will not fall for
    // these traps.
    /*
       if (current_magnification_size == 0)
       num_iterations = 0;
       if (current_magnification_size == 1)
       num_iterations = 2;
       if (current_magnification_size == 2)
       num_iterations = 3;
       if (current_magnification_size == 3)
       num_iterations = 3;
       if (current_magnification_size == 4)
       num_iterations = 3;
       if (current_magnification_size == 5)
       num_iterations = 4;
       if (current_magnification_size == 6)
       num_iterations = 20;
       */

    /*
    if (rc0_v.size() > 1)
    {
      if (rc0_v[rc0_v.size()-1] >= rc0_v[rc0_v.size()-2] && rc0_v.back() > 0.9)
        num_iterations++;

      if (rc0_v[rc0_v.size()-1] < rc0_v[rc0_v.size()-2])
        num_iterations--;

      if (num_iterations < 0)
        num_iterations = 0;
    }
    */



    int tr_iterations = -1;
    double int_time_trans = 0.0;

#if (defined TIMES) || (defined LOGS)
    std::chrono::high_resolution_clock::time_point start_translation =
      std::chrono::high_resolution_clock::now();
#endif

    double trans_criterion = Translation::tff(real_scan,
      *result_pose, map, num_iterations, ip.xy_bound, true, r2rp,
      &tr_iterations, &intersections_time, result_pose);

#if (defined TIMES) || (defined LOGS)
    std::chrono::high_resolution_clock::time_point end_translation =
      std::chrono::high_resolution_clock::now();

    op->translation_times += std::chrono::duration_cast<
      std::chrono::duration<double> >(end_translation-start_translation).count();

    op->intersections_times += intersections_time.count();
#endif

#if (defined LOGS)
    op->translation_iterations += tr_iterations;
#endif

    tc_v.push_back(trans_criterion);

#if defined (DEBUG)
    printf("rc0 = %f\n", rc0_v.back());
    printf("rc1 = %f\n", rc1_v.back());
    printf("tc  = %f\n", tc_v.back());
#endif

    xs.push_back(std::get<0>(*result_pose));
    ys.push_back(std::get<1>(*result_pose));

#if (defined LOGS)
    std::tuple<double,double,double> traj_i;
    std::get<0>(traj_i) = xs.back();
    std::get<1>(traj_i) = ys.back();
    std::get<2>(traj_i) = ts.back();
    op->trajectory.push_back(traj_i);
#endif


    // ----------------------- Recovery modes ----------------------------------
    bool l2_recovery = false;

    // Perilous pose at exterior of map's bounds detected
    if (tc_v.back() == -2.0)
    {
#if defined (DEBUG)
      printf("Will trigger recovery due to condition 0\n");
#endif
      l2_recovery = true;
    }

/*
    // Detect non-smooth entry into maximum magnification size
    if (current_magnification_size >= max_magnification_size-1)
    {
      if(rc0_v.back() < 0.99)
      {
#if defined (DEBUG)
        printf("Will trigger recovery due to condition 2\n");
#endif
        l2_recovery = true;
      }
    }
*/

    // Impose strict measures when on the final straight
    if (current_magnification_size >= max_magnification_size)
    {
      /*
      if(rc0_v.back() < 0.995)
      {
#if defined (DEBUG)
        printf("Will trigger recovery due to condition 2.2\n");
#endif
        l2_recovery = true;
      }
      */

      // Detect when stuck at awkward pose
      // trans_criterion is a measure of the deviation between rays from the
      // same pose; wherefore this should be proportionate to the
      // square root of the sum of variance estimates of the laser's rays
      // and the rays of the virtual scan
      // (assuming they are distributed normally)
      if (tc_v.back() > 4*sqrtf(ip.sigma_noise_real*ip.sigma_noise_real+
                                ip.sigma_noise_map*ip.sigma_noise_map)
                          + 0.001)
      {
#if defined (DEBUG)
        printf("Will trigger recovery due to condition 3\n");
#endif
        l2_recovery = true;
      }
    }

    // Do not allow more than `max_counter` iterations per resolution
    if (counter > max_counter)
    {
#if defined (DEBUG)
      printf("Will trigger recovery due to condition 4\n");
#endif
      //l2_recovery = true;

      counter = 0;
      current_magnification_size++;
    }


    // Recover if need be
    if (l2_recovery)
    {
      if (num_recoveries > max_recoveries)
      {
#if defined (DEBUG)
        printf("ERROR: MAXIMUM RECOVERIES\n");
#endif
        break;
      }

      num_recoveries++;
      l2recovery(virtual_pose, map, ip.xy_bound, ip.t_bound, result_pose);

      counter = min_counter;
      current_magnification_size = min_magnification_size;
    }
    else
    {
      counter++;

      // -------------------------- Level-up -------------------------------------
      double xy_eps = 10.1;
      double t_eps = 0.00001; // 0.0001
      if (canGiveNoMore(xs,ys,ts, xy_eps, t_eps) && counter > min_counter)
      {
        current_magnification_size += 1;
        counter = 0;

        if (rc0_v.back() > 0.99999 && up_lock == false)
        {
          current_magnification_size = max_magnification_size;
          up_lock = true;
        }
      }
    }

    total_iterations++;
  }

  std::chrono::high_resolution_clock::time_point end =
    std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(end-start);


#if defined (TIMES)
  printf("%f [Match::fmt]\n", elapsed.count());
#endif

  op->exec_time = elapsed.count();
  op->rc = rc0_v.back();
  op->tc = tc_v.back();
#if defined (LOGS)
  op->rotation_iterations = total_iterations;
  op->num_recoveries = num_recoveries;
#endif
}


/*******************************************************************************
*/
void Match::ndt(
  const std::vector< double >& real_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  const input_params& ip, output_params* op,
  std::tuple<double,double,double>* result_pose)
{
  // Start the clock
  std::chrono::high_resolution_clock::time_point start =
    std::chrono::high_resolution_clock::now();

  std::tuple<double,double,double> current_pose = virtual_pose;


  // The total correspondence error
  std::vector<double> c_errors;
  std::vector<double> traj_x;
  std::vector<double> traj_y;
  std::vector<double> traj_t;

  pcl::PointCloud<pcl::PointXYZ>::Ptr
    real_scan_pcl(new pcl::PointCloud<pcl::PointXYZ>);

  double x0 = 0.0;
  double y0 = 0.0;
  double t0 = 0.0;
  double z0 = 0;
  for (int i = 0; i < real_scan.size(); i++)
  {
    double px = x0 + real_scan[i]*cos(-M_PI + i*2*M_PI/real_scan.size() + t0);
    double py = y0 + real_scan[i]*sin(-M_PI + i*2*M_PI/real_scan.size() + t0);

    pcl::PointXYZ p;
    p.x = px;
    p.y = py;
    p.z = z0;

    real_scan_pcl->push_back(p);
  }





  for(int k = 0; k < ip.num_iterations; k++)
  {
    std::vector<double> virtual_scan;
    Utils::scanFromPose(current_pose, map, real_scan.size(), &virtual_scan);

    // https://pointclouds.org/documentation/tutorials/normal_distributions_transform.html
    pcl::PointCloud<pcl::PointXYZ>::Ptr
      virtual_scan_pcl(new pcl::PointCloud<pcl::PointXYZ>);

    for (int i = 0; i < virtual_scan.size(); i++)
    {
      double px =
        x0 + virtual_scan[i]*cos(-M_PI + i*2*M_PI/virtual_scan.size() + t0);
      double py =
        y0 + virtual_scan[i]*sin(-M_PI + i*2*M_PI/virtual_scan.size() + t0);

      pcl::PointXYZ p;
      p.x = px;
      p.y = py;
      p.z = z0;

      virtual_scan_pcl->push_back(p);
    }


    //Filter the input scan to about 10% of the original size to improve the matching speed.
    pcl::PointCloud<pcl::PointXYZ>::Ptr filtered_cloud(new pcl::PointCloud<pcl::PointXYZ>);

    pcl::ApproximateVoxelGrid<pcl::PointXYZ> approximate_voxel_filter;

    approximate_voxel_filter.setLeafSize(0.01, 0.01, 0.01); // 0.01, fixed
    approximate_voxel_filter.setInputCloud(virtual_scan_pcl);
    approximate_voxel_filter.filter(*filtered_cloud);

    /*
       cout << "Filtered cloud contains " << filtered_cloud->size()
       << " data points " << endl;
       */

    //Initialize normal distribution transformation (NDT)
    pcl::NormalDistributionsTransform<pcl::PointXYZ, pcl::PointXYZ> ndt;

    //Set the NDT parameter dependent on the scale
    //Set the minimum conversion difference for the termination condition
    ndt.setTransformationEpsilon(0.0001);

    //Set the maximum step size for More-Thuente line search
    ndt.setStepSize(0.5);

    //Set the resolution of the NDT grid structure (VoxelGridCovariance)
    ndt.setResolution(0.5);

    //Set the maximum number of matching iterations
    int max_ndt_iterations = 35;
    ndt.setMaximumIterations(max_ndt_iterations);

    // Set the point cloud to be registered
    ndt.setInputSource(filtered_cloud);

    //Set point cloud registration target
    ndt.setInputTarget(real_scan_pcl);

    /*
    //Set the initial alignment estimation result obtained by using the robot ranging method
    Eigen::AngleAxisf init_rotation(0.6931, Eigen::Vector3f::UnitZ());
    Eigen::Translation3f init_translation(1.79387, 0.720047, 0);
    Eigen::Matrix4f init_guess = (init_translation * init_rotation).matrix();*/

    //Calculate the required rigid body transformation to match the input point cloud to the target point cloud
    pcl::PointCloud<pcl::PointXYZ>::Ptr output_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    ndt.align(*output_cloud);
    /*
       cout << "NDT has converged:" << ndt.hasConverged()
       << " score: " << ndt.getFitnessScore() << endl;
       cout << "Transformation matrix:\n" << ndt.getFinalTransformation() << endl;
       */

    //std::cout << "Transformation matrix:\n" << ndt.getFinalTransformation() << std::endl;

    /*
       if (fmod(i,2) == 1)
       {
       std::get<0>(current_pose) -= ndt.getFinalTransformation()(0,3);
       std::get<1>(current_pose) -= ndt.getFinalTransformation()(1,3);
       }
       else
       {
       double t = asin(ndt.getFinalTransformation()(0,1));
       Utils::wrapAngle(&t);
       std::get<2>(current_pose) += t;
       Utils::wrapAngle(&std::get<2>(current_pose));
       }
       */
    double t = asin(ndt.getFinalTransformation()(0,1));
    Utils::wrapAngle(&t);
    std::get<2>(current_pose) += t;
    Utils::wrapAngle(&std::get<2>(current_pose));

    double x_pre = ndt.getFinalTransformation()(0,3);
    double y_pre = ndt.getFinalTransformation()(1,3);
    double t_est = std::get<2>(current_pose);
    double dx = cos(t_est) * x_pre - sin(t_est) * y_pre;
    double dy = sin(t_est) * x_pre + cos(t_est) * y_pre;

    std::get<0>(current_pose) -= dx;
    std::get<1>(current_pose) -= dy;



    if(!Utils::isPositionInMap(current_pose, map))
      l2recovery(virtual_pose, map, ip.xy_bound, ip.t_bound, &current_pose);

  }


  // Stop the clock
  std::chrono::high_resolution_clock::time_point end =
    std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(end-start);

  op->exec_time = elapsed.count();

  *result_pose = current_pose;

  op->rotation_iterations = ip.num_iterations;
  op->translation_iterations = ip.num_iterations;
}


/*******************************************************************************
*/
void Match::l2recovery(
  const std::tuple<double,double,double>& input_pose,
  const std::vector< std::pair<double,double> >& map,
  const double& xy_bound, const double& t_bound,
  std::tuple<double,double,double>* output_pose)
{
#if defined (PRINTS)
  printf("*********************************\n");
  printf("************CAUTION**************\n");
  printf("Level 2 recovery mode activated\n");
  printf("*********************************\n");
#endif

  /*
     do Utils::generatePose(input_pose,
     xy_bound, t_bound, output_pose);
     while(!Utils::isPositionInMap(*output_pose, map));
     */

  while(!Utils::generatePose(input_pose, map,
      1*xy_bound, t_bound, 0.0, output_pose));
}


/*******************************************************************************
*/
void Match::skg(
  const std::vector< double >& real_scan,
  const std::tuple<double,double,double>& real_pose,
  const std::tuple<double,double,double>& virtual_pose,
  const std::vector< std::pair<double,double> >& map,
  const fftw_plan& r2rp,
  const input_params& ip, output_params* op,
  std::tuple<double,double,double>* result_pose)
{
  std::chrono::high_resolution_clock::time_point start =
    std::chrono::high_resolution_clock::now();

  *result_pose = virtual_pose;

  // Maximum counter value means a new recovery attempt
  int min_counter = 0;
  int max_counter = 20;
  int counter = min_counter;

  // By a factor of what do you need to over-sample angularly?
  unsigned int min_magnification_size = 2;
  unsigned int max_magnification_size = 4;
  unsigned int current_magnification_size = min_magnification_size;

  // How many times do I attempt recovery?
  unsigned int num_recoveries = 0;
  unsigned int max_recoveries = 1000;

  // These three vectors hold the trajectory for each iteration
  std::vector<double> xs;
  std::vector<double> ys;
  std::vector<double> ts;

  // Two rotation criteria
  std::vector<double> rc0_v;
  std::vector<double> rc1_v;

  // One translation criterion
  std::vector<double> tc_v;

  std::vector<double> dxys;
  std::chrono::duration<double> intersections_time;

  // The best candidate angle found at each iterations is stored and made a
  // candidate each time. Its criterion is its translation criterion after
  // ni-1 translations
  double best_min_tc = 100000.0;
  std::tuple<double,double,double> best_cand_pose = *result_pose;

  // A lock for going overdrive when the rotation criterion is near-excellent
  int total_iterations = 0;
  int num_iterations = 0;


  while (current_magnification_size <= max_magnification_size)
  {

#if defined (DEBUG)
    printf("current_magnification_size = %d ---\n", current_magnification_size);
    printf("counter                    = %d ---\n", counter);
    printf("real pose (%f,%f,%f) [skg]\n",
      std::get<0>(real_pose),
      std::get<1>(real_pose),
      std::get<2>(real_pose));
    printf("     pose (%f,%f,%f) [skg]\n",
      std::get<0>(*result_pose),
      std::get<1>(*result_pose),
      std::get<2>(*result_pose));
#endif

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // ------------------ Rotation correction phase ----------------------------
    std::vector<double> rc0;
    std::vector<double> rc1;
    std::vector<double> cand_angles;
    std::vector< std::tuple<double,double,double> > cand_poses;

    cand_angles = Rotation::skg(real_scan, *result_pose, map,
      current_magnification_size, r2rp, &rc0, &rc1, &intersections_time);

    for (unsigned int a = 0; a < cand_angles.size(); a++)
    {
      std::tuple<double,double,double> cand_pose_a = *result_pose;
      std::get<2>(cand_pose_a) += cand_angles[a];
      Utils::wrapAngle(&std::get<2>(cand_pose_a));
      cand_poses.push_back(cand_pose_a);
    }
    cand_poses.push_back(best_cand_pose);

    // ------------------ Candidate angles sifting -----------------------------
    unsigned int min_tc_idx = 0;
    if (cand_angles.size() > 1)
    {
      std::vector<double> tcs_sift;
      for (unsigned int ca = 0; ca < cand_poses.size(); ca++)
      {
        // How many test iterations?
        unsigned int ni = 2;
        int tr_i = 0;

        std::tuple<double,double,double> cand_pose = cand_poses[ca];

        double tc = Translation::tff(real_scan, cand_pose, map,
          ni, ip.xy_bound, false, r2rp, &tr_i, &intersections_time, &cand_pose);
        cand_poses.at(ca) = cand_pose;

        if (tc == -2.0)
          tcs_sift.push_back(1000000.0);
        else
          tcs_sift.push_back(tc);
      }

      // The index of the angle with the least translation criterion
      min_tc_idx =
        std::min_element(tcs_sift.begin(), tcs_sift.end()) - tcs_sift.begin();


      // Check if the newly-found angle is the angle with the least
      // translation criterion so far
      if (tcs_sift[min_tc_idx] < best_min_tc)
      {
        best_min_tc = tcs_sift[min_tc_idx];
        best_cand_pose = cand_poses[min_tc_idx];
      }
    }

    // Update the current estimate with the one that sports the
    // least translation criterion overall
    *result_pose = cand_poses[min_tc_idx]; // results in loops; avoid
    //std::get<2>(*result_pose) = std::get<2>(cand_poses[min_tc_idx]);
    //Utils::wrapAngle(&std::get<2>(*result_pose));

    // ... and store it
    ts.push_back(std::get<2>(*result_pose));

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // ---------------- Translation correction phase ---------------------------
    num_iterations = ip.num_iterations;
      //(current_magnification_size)*ip.num_iterations;


    int tr_iterations = -1;
    double int_time_trans = 0.0;


    double trans_criterion = Translation::tff(real_scan,
      *result_pose, map, num_iterations, ip.xy_bound, false, r2rp,
      &tr_iterations, &intersections_time, result_pose);


    tc_v.push_back(trans_criterion);


    xs.push_back(std::get<0>(*result_pose));
    ys.push_back(std::get<1>(*result_pose));


    // ----------------------- Recovery modes ----------------------------------
    bool l2_recovery = false;

    // Perilous pose at exterior of map's bounds detected
    if (tc_v.back() == -2.0)
      l2_recovery = true;


    // Impose strict measures when on the final straight
    if (current_magnification_size >= max_magnification_size)
    {
      // Detect when stuck at awkward pose
      // trans_criterion is a measure of the deviation between rays from the
      // same pose; wherefore this should be proportionate to the
      // square root of the sum of variance estimates of the laser's rays
      // and the rays of the virtual scan
      // (assuming they are distributed normally)
      if (tc_v.back() > 4*sqrtf(ip.sigma_noise_real*ip.sigma_noise_real+
          ip.sigma_noise_map*ip.sigma_noise_map)
        + 0.001)
      {
        l2_recovery = true;
      }
    }

    // Do not allow more than `max_counter` iterations per resolution
    if (counter > max_counter)
    {
      counter = 0;
      current_magnification_size++;
    }

    double dx = fabs(std::get<0>(*result_pose) - std::get<0>(virtual_pose));
    double dy = fabs(std::get<1>(*result_pose) - std::get<1>(virtual_pose));
    double dt = fabs(std::get<2>(*result_pose) - std::get<2>(virtual_pose));
    Utils::wrapAngle(&dt);
    if (dx > 2*ip.xy_bound || dy > 2*ip.xy_bound || dt > 2*ip.t_bound)
      l2_recovery = true;


    // Recover if need be
    if (l2_recovery)
    {
      if (num_recoveries > max_recoveries)
        break;

      num_recoveries++;
      l2recovery(virtual_pose, map, ip.xy_bound, ip.t_bound, result_pose);

      counter = min_counter;
      current_magnification_size = min_magnification_size;

      best_cand_pose = *result_pose;
    }
    else
    {
      counter++;

      // -------------------------- Level-up -------------------------------------
      double xy_eps = 10.01;
      double t_eps = 0.00001; // 0.00001
      if (canGiveNoMore(xs,ys,ts, xy_eps, t_eps) && counter > min_counter+1)
      {
        current_magnification_size += 1;
        counter = 0;
      }
    }

    total_iterations++;
  }

  std::chrono::high_resolution_clock::time_point end =
    std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(end-start);

  op->exec_time = elapsed.count();
}
