#include "dft_utils.h"


/*******************************************************************************
*/
void DFTUtils::fftshift(std::vector<double>* vec)
{
#ifdef TIMES
    std::chrono::high_resolution_clock::time_point a =
      std::chrono::high_resolution_clock::now();
#endif

  std::rotate(
    vec->begin(),
    vec->begin() + static_cast<unsigned int>(vec->size()/2),
    vec->end());

#ifdef TIMES
    std::chrono::high_resolution_clock::time_point b =
      std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed =
      std::chrono::duration_cast< std::chrono::duration<double> >(b-a);

    printf("%f [fftshift]\n", elapsed.count());
#endif
}

/*******************************************************************************
 * @brief Calculates the X1 coefficient of the rays_diff input vector.
 * @param[in] rays_diff [const std::vector<double>&] The difference in range
 * between a world and a map scan.
 * @param[in] num_valid_rays [const int&] The number of valid rays (rays
 * whose difference in range is lower than a set threshold) between the
 * world and map scans.
 * @return [std::vector<double>] A vector of size two, of which the first
 * position holds the real part of the first DFT coefficient, and the
 * second the imaginary part of it.
 */
std::vector<double> DFTUtils::getFirstDFTCoefficient(
  const std::vector<double>& rays_diff)
{
#ifdef TIMES
    std::chrono::high_resolution_clock::time_point a =
      std::chrono::high_resolution_clock::now();
#endif

  // A vector holding the coefficients of the DFT
  std::vector<double> dft_coeff_vector;

  // Do the DFT thing
  std::vector<double> dft_coeffs = dft(rays_diff);

  // The real and imaginary part of the first coefficient are
  // out[1] and out[N-1] respectively

  // The real part of the first coefficient
  double x1_r = dft_coeffs[1];

  // The imaginary part of the first coefficient
  double x1_i = dft_coeffs[rays_diff.size()-1];

  // Is x1_r finite?
  if (std::isfinite(x1_r))
    dft_coeff_vector.push_back(x1_r);
  else
    dft_coeff_vector.push_back(0.0);

  // Is x1_i finite?
  if (std::isfinite(x1_i))
    dft_coeff_vector.push_back(x1_i);
  else
    dft_coeff_vector.push_back(0.0);

#ifdef TIMES
    std::chrono::high_resolution_clock::time_point b =
      std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed =
      std::chrono::duration_cast< std::chrono::duration<double> >(b-a);

    printf("%f [getFirstDFTCoefficient]\n", elapsed.count());
#endif

  return dft_coeff_vector;
}


/*******************************************************************************
 */
std::vector<double> DFTUtils::getFirstDFTCoefficient(
  const std::vector<double>& rays_diff,
  const fftw_plan& r2rp)
{
#ifdef TIMES
    std::chrono::high_resolution_clock::time_point a =
      std::chrono::high_resolution_clock::now();
#endif

  // A vector holding the coefficients of the DFT
  std::vector<double> dft_coeff_vector;

  // Do the DFT thing
  std::vector<double> dft_coeffs = dft(rays_diff, r2rp);

  // The real and imaginary part of the first coefficient are
  // out[1] and out[N-1] respectively

  // The real part of the first coefficient
  double x1_r = dft_coeffs[1];

  // The imaginary part of the first coefficient
  double x1_i = dft_coeffs[rays_diff.size()-1];

  // Is x1_r finite?
  if (std::isfinite(x1_r))
    dft_coeff_vector.push_back(x1_r);
  else
    dft_coeff_vector.push_back(0.0);

  // Is x1_i finite?
  if (std::isfinite(x1_i))
    dft_coeff_vector.push_back(x1_i);
  else
    dft_coeff_vector.push_back(0.0);

#ifdef TIMES
    std::chrono::high_resolution_clock::time_point b =
      std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed =
      std::chrono::duration_cast< std::chrono::duration<double> >(b-a);

    printf("%f [getFirstDFTCoefficient]\n", elapsed.count());
#endif

  return dft_coeff_vector;
}


/*******************************************************************************
*/
  std::vector< std::pair<double, double> >
DFTUtils::getDFTCoefficientsPairs(const std::vector<double>& coeffs)
{
#ifdef TIMES
    std::chrono::high_resolution_clock::time_point a =
      std::chrono::high_resolution_clock::now();
#endif

  std::vector< std::pair<double, double> > fft_coeff_pairs;
  for (int i = 0; i <= coeffs.size()/2; i++)
  {
    if (i == 0 || i == coeffs.size()/2)
      fft_coeff_pairs.push_back(std::make_pair(coeffs[i], 0.0));
    else
    {
      fft_coeff_pairs.push_back(
        std::make_pair(coeffs[i], coeffs[coeffs.size()-i]));
    }
  }

  std::vector< std::pair<double, double> > fft_coeff_pairs_bak =
    fft_coeff_pairs;
  for (int i = fft_coeff_pairs_bak.size()-2; i > 0; i--)
  {
    fft_coeff_pairs.push_back(
      std::make_pair(fft_coeff_pairs_bak[i].first, -fft_coeff_pairs_bak[i].second));
  }

#ifdef TIMES
    std::chrono::high_resolution_clock::time_point b =
      std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed =
      std::chrono::duration_cast< std::chrono::duration<double> >(b-a);

    printf("%f [getDFTCoefficientsPairs]\n", elapsed.count());
#endif

  return fft_coeff_pairs;
}


/*******************************************************************************
 * @brief Performs DFT in a vector of doubles via fftw. Returns the DFT
 * coefficients vector in the order described in
 * http://www.fftw.org/fftw3_doc/Real_002dto_002dReal-Transform-Kinds.html#Real_002dto_002dReal-Transform-Kinds.
 * @param[in] rays_diff [const std::vector<double>&] The vector of differences
 * in range between a world scan and a map scan.
 * @return [std::vector<double>] The vector's DFT coefficients.
 */
std::vector<double> DFTUtils::dft(const std::vector<double>& rays_diff)
{
#ifdef TIMES
    std::chrono::high_resolution_clock::time_point a =
      std::chrono::high_resolution_clock::now();
#endif

  double* in;
  double* out;

  const size_t num_rays = rays_diff.size();

  in = (double*) fftw_malloc(num_rays * sizeof(double));
  out = (double*) fftw_malloc(num_rays * sizeof(double));

  // Create plan
  fftw_plan p = fftw_plan_r2r_1d(num_rays, in, out, FFTW_R2HC, FFTW_MEASURE);

  // Transfer the input vector to a structure preferred by fftw
  for (unsigned int i = 0; i < num_rays; i++)
    in[i] = rays_diff[i];

  // Execute plan
  fftw_execute(p);

  // Store all DFT coefficients
  std::vector<double> dft_coeff_vector;
  for (unsigned int i = 0; i < num_rays; i++)
    dft_coeff_vector.push_back(out[i]);

  // Free memory
  fftw_destroy_plan(p);
  fftw_free(out);
  fftw_free(in);

#ifdef TIMES
    std::chrono::high_resolution_clock::time_point b =
      std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed =
      std::chrono::duration_cast< std::chrono::duration<double> >(b-a);

    printf("%f [dft]\n", elapsed.count());
#endif

  return dft_coeff_vector;
}


/*******************************************************************************
 */
std::vector<double> DFTUtils::dft(const std::vector<double>& rays_diff,
  const fftw_plan& r2rp)
{
#ifdef TIMES
    std::chrono::high_resolution_clock::time_point a =
      std::chrono::high_resolution_clock::now();
#endif

  double* in;
  double* out;

  const size_t num_rays = rays_diff.size();

  in = (double*) fftw_malloc(num_rays * sizeof(double));
  out = (double*) fftw_malloc(num_rays * sizeof(double));

  // Create plan
  //fftw_plan p = fftw_plan_r2r_1d(num_rays, in, out, FFTW_R2HC, FFTW_MEASURE);

  // Transfer the input vector to a structure preferred by fftw
  for (unsigned int i = 0; i < num_rays; i++)
    in[i] = rays_diff[i];


  // Execute plan
  fftw_execute_r2r(r2rp, in, out);

  // Store all DFT coefficients
  std::vector<double> dft_coeff_vector;
  for (unsigned int i = 0; i < num_rays; i++)
    dft_coeff_vector.push_back(out[i]);

  // Free memory
  //fftw_destroy_plan(p);
  fftw_free(out);
  fftw_free(in);

#ifdef TIMES
    std::chrono::high_resolution_clock::time_point b =
      std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed =
      std::chrono::duration_cast< std::chrono::duration<double> >(b-a);

    printf("%f [dft]\n", elapsed.count());
#endif

  return dft_coeff_vector;
}


/*******************************************************************************
 */
std::vector< std::vector<double> > DFTUtils::dftBatch(
  const std::vector< std::vector<double> >& scans)
{
#ifdef TIMES
  std::chrono::high_resolution_clock::time_point a =
    std::chrono::high_resolution_clock::now();
#endif

  assert(scans.size() > 0);

  // What will be returned
  std::vector< std::vector<double> > coeff_vector_v;

  // Input/output arrays for fftw
  double* in;
  double* out;

  const size_t num_rays = scans[0].size();

  in = (double*) fftw_malloc(num_rays * sizeof(double));
  out = (double*) fftw_malloc(num_rays * sizeof(double));

  // Create plan once
  fftw_plan p = fftw_plan_r2r_1d(num_rays, in, out, FFTW_R2HC, FFTW_MEASURE);

  for (unsigned int v = 0; v < scans.size(); v++)
  {
    // Transfer the input vector to a structure preferred by fftw
    for (unsigned int i = 0; i < num_rays; i++)
      in[i] = scans[v][i];

    // Execute plan with new input/output arrays
    fftw_execute_r2r(p, in, out);

    // Store all DFT coefficients for the v-th scan
    std::vector<double> dft_coeffs;
    for (unsigned int i = 0; i < num_rays; i++)
      dft_coeffs.push_back(out[i]);

    coeff_vector_v.push_back(dft_coeffs);
  }

  // Free memory
  fftw_destroy_plan(p);
  fftw_free(out);
  fftw_free(in);


#ifdef TIMES
  std::chrono::high_resolution_clock::time_point b =
    std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(b-a);

  printf("%f [dftBatch]\n", elapsed.count());
#endif

  return coeff_vector_v;
}


/*******************************************************************************
 */
std::vector< std::vector<double> > DFTUtils::dftBatch(
  const std::vector< std::vector<double> >& scans,
  const fftw_plan& r2rp)
{
#ifdef TIMES
  std::chrono::high_resolution_clock::time_point a =
    std::chrono::high_resolution_clock::now();
#endif

  assert(scans.size() > 0);

  // What will be returned
  std::vector< std::vector<double> > coeff_vector_v;

  // Input/output arrays for fftw
  double* in;
  double* out;

  const size_t num_rays = scans[0].size();

  in = (double*) fftw_malloc(num_rays * sizeof(double));
  out = (double*) fftw_malloc(num_rays * sizeof(double));

  // Create plan once
  //fftw_plan p = fftw_plan_r2r_1d(num_rays, in, out, FFTW_R2HC, FFTW_MEASURE);

  for (unsigned int v = 0; v < scans.size(); v++)
  {
    // Transfer the input vector to a structure preferred by fftw
    for (unsigned int i = 0; i < num_rays; i++)
      in[i] = scans[v][i];

    // Execute plan with new input/output arrays
    fftw_execute_r2r(r2rp, in, out);

    // Store all DFT coefficients for the v-th scan
    std::vector<double> dft_coeffs;
    for (unsigned int i = 0; i < num_rays; i++)
      dft_coeffs.push_back(out[i]);

    coeff_vector_v.push_back(dft_coeffs);
  }

  // Free memory
  //fftw_destroy_plan(p);
  fftw_free(out);
  fftw_free(in);


#ifdef TIMES
  std::chrono::high_resolution_clock::time_point b =
    std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(b-a);

  printf("%f [dftBatch]\n", elapsed.count());
#endif

  return coeff_vector_v;
}


/*******************************************************************************
*/
std::vector<double> DFTUtils::idft(
  const std::vector<std::pair<double, double> >& rays_diff)
{
#ifdef TIMES
  std::chrono::high_resolution_clock::time_point a =
    std::chrono::high_resolution_clock::now();
#endif

  fftw_complex* in;
  double* out;

  const size_t num_rays = rays_diff.size();

  in = (fftw_complex*) fftw_malloc(num_rays * sizeof(fftw_complex));
  out = (double*) fftw_malloc(num_rays * sizeof(double));

  // Create plan
  fftw_plan p = fftw_plan_dft_c2r_1d(num_rays, in, out, FFTW_MEASURE);

  // Transfer the input vector to a structure preferred by fftw
  for (unsigned int i = 0; i < num_rays; i++)
  {
    in[i][0] = rays_diff[i].first;
    in[i][1] = rays_diff[i].second;
  }

  // Execute plan
  fftw_execute(p);

  // Store all DFT coefficients
  std::vector<double> dft_coeff_vector;
  for (unsigned int i = 0; i < num_rays; i++)
    dft_coeff_vector.push_back(out[i]/num_rays);

  // Free memory
  fftw_destroy_plan(p);
  fftw_free(out);
  fftw_free(in);

#ifdef TIMES
  std::chrono::high_resolution_clock::time_point b =
    std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(b-a);

  printf("%f [idft]\n", elapsed.count());
#endif

  return dft_coeff_vector;
}


/*******************************************************************************
*/
std::vector< std::vector<double> > DFTUtils::idftBatch(
  const std::vector< std::vector<std::pair<double, double> > >& scans)
{
#ifdef TIMES
  std::chrono::high_resolution_clock::time_point a =
    std::chrono::high_resolution_clock::now();
#endif

  assert(scans.size() > 0);

  // What will be returned
  std::vector< std::vector<double> > dft_coeffs_v;

  fftw_complex* in;
  double* out;

  const size_t num_rays = scans[0].size();

  in = (fftw_complex*) fftw_malloc(num_rays * sizeof(fftw_complex));
  out = (double*) fftw_malloc(num_rays * sizeof(double));

  // Create plan once
  fftw_plan p = fftw_plan_dft_c2r_1d(num_rays, in, out, FFTW_MEASURE);


  for (unsigned int v = 0; v < scans.size(); v++)
  {
    // Transfer the input vector to a structure preferred by fftw
    for (unsigned int i = 0; i < num_rays; i++)
    {
      in[i][0] = scans[v][i].first;
      in[i][1] = scans[v][i].second;
    }

    // Execute plan
    fftw_execute_dft_c2r(p, in, out);

    // Store all DFT coefficients
    std::vector<double> dft_coeffs;
    for (unsigned int i = 0; i < num_rays; i++)
      dft_coeffs.push_back(out[i]/num_rays);

    dft_coeffs_v.push_back(dft_coeffs);
  }

  // Free memory
  fftw_destroy_plan(p);
  fftw_free(out);
  fftw_free(in);

#ifdef TIMES
  std::chrono::high_resolution_clock::time_point b =
    std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(b-a);

  printf("%f [idftBatch]\n", elapsed.count());
#endif

  return dft_coeffs_v;
}


/*******************************************************************************
*/
std::vector< std::vector<double> > DFTUtils::idftBatch(
  const std::vector< std::vector<std::pair<double, double> > >& scans,
  const fftw_plan& c2rp)
{
#ifdef TIMES
  std::chrono::high_resolution_clock::time_point a =
    std::chrono::high_resolution_clock::now();
#endif

  assert(scans.size() > 0);

  // What will be returned
  std::vector< std::vector<double> > dft_coeffs_v;

  fftw_complex* in;
  double* out;

  const size_t num_rays = scans[0].size();

  in = (fftw_complex*) fftw_malloc(num_rays * sizeof(fftw_complex));
  out = (double*) fftw_malloc(num_rays * sizeof(double));

  // Create plan once
  //fftw_plan p = fftw_plan_dft_c2r_1d(num_rays, in, out, FFTW_MEASURE);


  for (unsigned int v = 0; v < scans.size(); v++)
  {
    // Transfer the input vector to a structure preferred by fftw
    for (unsigned int i = 0; i < num_rays; i++)
    {
      in[i][0] = scans[v][i].first;
      in[i][1] = scans[v][i].second;
    }

    // Execute plan
    fftw_execute_dft_c2r(c2rp, in, out);

    // Store all DFT coefficients
    std::vector<double> dft_coeffs;
    for (unsigned int i = 0; i < num_rays; i++)
      dft_coeffs.push_back(out[i]/num_rays);

    dft_coeffs_v.push_back(dft_coeffs);
  }

  // Free memory
  //fftw_destroy_plan(p);
  fftw_free(out);
  fftw_free(in);

#ifdef TIMES
  std::chrono::high_resolution_clock::time_point b =
    std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed =
    std::chrono::duration_cast< std::chrono::duration<double> >(b-a);

  printf("%f [idftBatch]\n", elapsed.count());
#endif

  return dft_coeffs_v;
}
