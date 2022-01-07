#include "dataset_utils.h"

/*******************************************************************************
*/
  std::vector< std::vector< std::pair<double,double> > >
DatasetUtils::dataset2points(const char* dataset_filepath)
{
  std::vector< std::vector<double> > ranges;
  std::vector< std::tuple<double,double,double> > poses;

  readDataset(dataset_filepath, &ranges, &poses);

  int num_scans = ranges.size();
  int num_rays = ranges[0].size();
  double angle_span = M_PI;

  std::vector< std::vector< std::pair<double,double> > > polygons;
  for (int s = 0; s < ranges.size(); s++)
  {
    double px = std::get<0>(poses[s]);
    double py = std::get<1>(poses[s]);
    double pt = std::get<2>(poses[s]);

    std::vector< std::pair<double,double> > polygon;

    for (int r = 0; r < ranges[s].size(); r++)
    {
      double x =
        px + ranges[s][r] * cos(r*angle_span/(num_rays-1) + pt -angle_span/2);
      double y =
        py + ranges[s][r] * sin(r*angle_span/(num_rays-1) + pt -angle_span/2);

      polygon.push_back(std::make_pair(x,y));
    }

    polygons.push_back(polygon);
  }

  return polygons;
}


/*******************************************************************************
*/
void DatasetUtils::dataset2rangesAndPose(
  const char* dataset_filepath,
  std::vector<double>* ranges,
  std::tuple<double,double,double>* pose)
{
  readDataset(dataset_filepath, ranges, pose);
}


/*******************************************************************************
*/
void DatasetUtils::readDataset(
  const char* filepath,
  std::vector<double>* ranges,
  std::tuple<double,double,double>* pose)
{
  // First read the first two number: they show
  // (1) the number of scans and
  // (2) the number of rays per scan.
  FILE* fp = fopen(filepath, "r");
  if (fp == NULL)
    exit(EXIT_FAILURE);

  char* line = NULL;
  size_t len = 0;

  unsigned int line_number = 0;
  long int num_scans = 0;
  long int num_rays = 0;
  while ((getline(&line, &len, fp)) != -1 && line_number < 1)
  {
    line_number++;

    char * pEnd;
    num_scans = strtol (line, &pEnd, 10);
    num_rays = strtol (pEnd, &pEnd, 10);
  }

  fclose(fp);

  if (line)
    free(line);


  // Begin for all scans
  fp = fopen(filepath, "r");
  line = NULL;
  len = 0;

  // The line number read at each iteration
  line_number = 0;

  // loop
  while ((getline(&line, &len, fp)) != -1)
  {
    line_number++;

    // We don't have to care about the first line now
    if (line_number == 1)
      continue;

    // These lines host the poses from which the scans were taken
    if ((line_number-1) % (num_rays+1) == 0)
    {
      // The pose from which the scan_number-th scan was taken
      std::string pose_d(line); // convert from char to string
      std::string::size_type sz; // alias of size_t

      double px = std::stod(pose_d,&sz);
      pose_d = pose_d.substr(sz);
      double py = std::stod(pose_d,&sz);
      double pt = std::stod(pose_d.substr(sz));
      Utils::wrapAngle(&pt);
      *pose = std::make_tuple(px,py,pt);

      continue;
    }

    // At this point we are in a line holding a range measurement; fo sho
    double range_d;
    assert(sscanf(line, "%lf", &range_d) == 1);
    ranges->push_back(range_d);
  }

  fclose(fp);

  if (line)
    free(line);
}


/*******************************************************************************
*/
void DatasetUtils::readDataset(
  const char* filepath,
  std::vector< std::vector<double> >* ranges,
  std::vector< std::tuple<double,double,double> >* poses)
{
  // First read the first two number: they show
  // (1) the number of scans and
  // (2) the number of rays per scan.
  FILE* fp = fopen(filepath, "r");
  if (fp == NULL)
    exit(EXIT_FAILURE);

  char* line = NULL;
  size_t len = 0;

  unsigned int line_number = 0;
  long int num_scans = 0;
  long int num_rays = 0;
  while ((getline(&line, &len, fp)) != -1 && line_number < 1)
  {
    line_number++;

    char * pEnd;
    num_scans = strtol (line, &pEnd, 10);
    num_rays = strtol (pEnd, &pEnd, 10);
  }

  fclose(fp);

  if (line)
    free(line);


  // Begin for all scans
  fp = fopen(filepath, "r");
  line = NULL;
  len = 0;

  // The line number read at each iteration
  line_number = 0;

  // A vector holding scan ranges for one scan
  std::vector<double> ranges_one_scan;

  // loop
  while ((getline(&line, &len, fp)) != -1)
  {
    line_number++;

    // We don't have to care about the first line now
    if (line_number == 1)
      continue;

    // These lines host the poses from which the scans were taken
    if ((line_number-1) % (num_rays+1) == 0)
    {
      // Finished with this scan
      ranges->push_back(ranges_one_scan);

      // Clear the vector so we can begin all over
      ranges_one_scan.clear();

      // The pose from which the scan_number-th scan was taken
      std::string pose(line); // convert from char to string
      std::string::size_type sz; // alias of size_t

      double px = std::stod(pose,&sz);
      pose = pose.substr(sz);
      double py = std::stod(pose,&sz);
      double pt = std::stod(pose.substr(sz));
      Utils::wrapAngle(&pt);
      poses->push_back(std::make_tuple(px,py,pt));

      continue;
    }

    // At this point we are in a line holding a range measurement; fo sho
    double range;
    assert(sscanf(line, "%lf", &range) == 1);
    ranges_one_scan.push_back(range);
  }

  fclose(fp);

  if (line)
    free(line);
}


/*******************************************************************************
*/
void DatasetUtils::printDataset(const char* dataset_filepath)
{
  std::vector< std::vector<double> > ranges;
  std::vector< std::tuple<double,double,double> > poses;

  readDataset(dataset_filepath, &ranges, &poses);

  for (int s = 0; s < ranges.size(); s++)
  {
    printf("NEW SCAN\n");
    for (int r = 0; r < ranges[s].size(); r++)
    {
      printf("r[%d] = %f\n", r, ranges[s][r]);
    }

    printf("FROM POSE (%f,%f,%f)\n",
      std::get<0>(poses[s]), std::get<1>(poses[s]), std::get<2>(poses[s]));
  }
}


/*******************************************************************************
void DatasetUtils::splitDataset()
{
  // First read the first two numbers: they show
  // (1) the number of scans and
  // (2) the number of rays per scan.
  FILE* fp = fopen(dataset_filepath_, "r");
  if (fp == NULL)
    exit(EXIT_FAILURE);

  char* line = NULL;
  size_t len = 0;

  unsigned int line_number = 0;
  long int num_scans = 0;
  long int num_rays = 0;
  while ((getline(&line, &len, fp)) != -1 && line_number < 1)
  {
    line_number++;

    char * pEnd;
    num_scans = strtol (line, &pEnd, 10);
    num_rays = strtol (pEnd, &pEnd, 10);
  }

  fclose(fp);

  if (line)
    free(line);

  // Begin for all scans
  fp = fopen(dataset_filepath_, "r");
  line = NULL;
  len = 0;

  // The line number read at each iteration
  line_number = 0;

  // loop
  int scan_id = 0;
  bool new_scan = true;

  while ((getline(&line, &len, fp)) != -1)
  {
    // Open output file
    std::string dataset_scan = base_path_ + "/../dataset/dataset_" +
      std::to_string(scan_id) + ".txt";

    std::ofstream file(dataset_scan.c_str(), std::ios::app);


    line_number++;

    // We don't have to care about the first line now
    if (line_number == 1)
      continue;

    if ((line_number-2) % 362 == 0 || line_number == 2)
      new_scan = true;
    else
      new_scan = false;

    if (new_scan)
      if (file.is_open())
      {
        file << "1 361" << std::endl;
      }


    // These lines host the poses from which the scans were taken
    if ((line_number-1) % (num_rays+1) == 0)
    {
      // Finished with this scan
      scan_id++;

      // The pose from which the scan_number-th scan was taken
      std::string pose(line); // convert from char to string
      std::string::size_type sz; // alias of size_t

      double px = std::stod(pose,&sz);
      pose = pose.substr(sz);
      double py = std::stod(pose,&sz);
      double pt = std::stod(pose.substr(sz));

      if (file.is_open())
      {
        file << px << " " << py << " " << pt << std::endl;
        file.close();
      }

      continue;
    }

    // At this point we are in a line holding a range measurement; fo sho
    double range;
    assert(sscanf(line, "%lf", &range) == 1);

    if (file.is_open())
    {
      file << range << std::endl;
    }
    else
      printf("[S2MSM] Could not split dataset\n");
  }
}
*/
