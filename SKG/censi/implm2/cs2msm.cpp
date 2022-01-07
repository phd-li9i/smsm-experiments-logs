#include <cs2msm.h>


/*******************************************************************************
*/
void CS2MSM::convertRealScanToLDP(const std::vector<double>& scan,
  const std::tuple<double,double,double>&pose,
  LDP& ldp)
{
  int n = scan.size();
  ldp = ld_alloc_new(n);

  float angle_range_deg = 360;

  float min_range = 0.001;
  float max_range = 100.0;
  float angle_min = -angle_range_deg / 2 * M_PI / 180 + std::get<2>(pose);
  float angle_inc = angle_range_deg * M_PI / 180 / n;

  for (int i = 0; i < n; i++)
  {
    ldp->readings[i] = scan[i];
    ldp->valid[i] = 1;
    ldp->cluster[i] = -1;
    ldp->theta[i] = angle_min + i * angle_inc;
  }

  ldp->min_theta = ldp->theta[0];
  ldp->max_theta = ldp->theta[n-1];

  ldp->odometry[0] = 0.0;
  ldp->odometry[1] = 0.0;
  ldp->odometry[2] = 0.0;

  ldp->true_pose[0] = 0.0;
  ldp->true_pose[1] = 0.0;
  ldp->true_pose[2] = 0.0;
}


/*******************************************************************************
*/
void CS2MSM::convertVirtualScanToLDP(const std::vector<double>& scan,
  const std::tuple<double,double,double>& virtual_pose,
  LDP& ldp)
{
  int n = scan.size();
  ldp = ld_alloc_new(n);

  float angle_range_deg = 360;

  float min_range = 0.001;
  float max_range = 100.0;
  float angle_min = -angle_range_deg / 2 * M_PI / 180 + std::get<2>(virtual_pose);
  float angle_inc = angle_range_deg * M_PI / 180 / n;

  for (int i = 0; i < n; i++)
  {
    ldp->readings[i] = scan[i];
    ldp->valid[i] = 1;
    ldp->cluster[i] = -1;
    ldp->theta[i] = angle_min + i * angle_inc;
  }

  ldp->min_theta = ldp->theta[0];
  ldp->max_theta = ldp->theta[n-1];

  ldp->odometry[0] = 0.0;
  ldp->odometry[1] = 0.0;
  ldp->odometry[2] = 0.0;

  ldp->true_pose[0] = 0.0;
  ldp->true_pose[1] = 0.0;
  ldp->true_pose[2] = 0.0;
}
