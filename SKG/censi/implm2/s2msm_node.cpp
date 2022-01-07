#include <s2msm.h>

int main(int argc, char** argv)
{
  if (argc < 11)
    S2MSM *obj = new S2MSM();
  else
  {
    unsigned int max_iterations = std::stod(argv[1]);
    unsigned int num_iterations = std::stod(argv[2]);
    unsigned int start_sample = std::stod(argv[3]);
    unsigned int end_sample = std::stod(argv[4]);
    double xy_uniform_displacement = std::stod(argv[5]);
    double t_uniform_displacement = std::stod(argv[6]);
    double sigma_noise_real = std::stod(argv[7]);
    double sigma_noise_map = std::stod(argv[8]);
    double invalid_rays_randomly_percent = std::stod(argv[9]);
    double invalid_rays_sequentially_percent = std::stod(argv[10]);
    double size_real_scan = std::stod(argv[11]);
    double size_map = std::stod(argv[12]);
    std::string method = argv[13];

    assert(max_iterations > 0);
    assert(num_iterations > 0);
    assert(start_sample >= 0);
    assert(end_sample > 0);
    assert(xy_uniform_displacement >= 0.0);
    assert(t_uniform_displacement >= 0.0);
    assert(sigma_noise_real >= 0.0);
    assert(sigma_noise_map >= 0.0);
    assert(invalid_rays_randomly_percent >= 0.0);
    assert(invalid_rays_sequentially_percent >= 0.0);
    assert(size_real_scan > 0.0);
    assert(size_map > 0.0);

    S2MSM *obj = new S2MSM(
      max_iterations,
      num_iterations,
      start_sample,
      end_sample,
      xy_uniform_displacement,
      t_uniform_displacement,
      sigma_noise_real,
      sigma_noise_map,
      invalid_rays_randomly_percent,
      invalid_rays_sequentially_percent,
      size_real_scan,
      size_map,
      method);
  }

  return 0;
}
