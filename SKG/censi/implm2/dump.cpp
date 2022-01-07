#include "dump.h"


/*******************************************************************************
*/
void Dump::map(const std::vector< std::pair<double,double> >& map,
  const std::string& dump_filepath)
{
  std::ofstream file(dump_filepath.c_str(), std::ios::trunc);

  if (file.is_open())
  {
    file << "mx = [];" << std::endl;
    file << "my = [];" << std::endl;
    for (int i = 0; i < map.size(); i++)
    {
      file << "mx = [mx " << map[i].first << "];" << std::endl;
      file << "my = [my " << map[i].second << "];" << std::endl;
    }

    file.close();
  }
  else
    printf("Could not log scans\n");
}

/*******************************************************************************
*/
void Dump::rangeScan(
  const std::vector<double>& real_scan,
  const std::vector<double>& virtual_scan,
  const std::string& dump_filepath)
{
  std::ofstream file(dump_filepath.c_str(), std::ios::trunc);

  if (file.is_open())
  {
    file << "rr = [];" << std::endl;
    for (int i = 0; i < real_scan.size(); i++)
      file << "rr = [rr " << real_scan[i] << "];" << std::endl;

    file << "rt = [];" << std::endl;
    for (int i = 0; i < real_scan.size(); i++)
      file << "rt = [rt " << i * 2 * M_PI / real_scan.size() << "];" << std::endl;

    file << "vr = [];" << std::endl;
    for (int i = 0; i < virtual_scan.size(); i++)
      file << "vr = [vr " << virtual_scan[i] << "];" << std::endl;

    file << "vt = [];" << std::endl;
    for (int i = 0; i < virtual_scan.size(); i++)
      file << "vt = [vt " << i * 2 * M_PI / virtual_scan.size() << "];" << std::endl;

    file.close();
  }
  else
    printf("Could not log range scans\n");
}


/*******************************************************************************
*/
void Dump::scan(
  const std::vector<double>& real_scan,
  const std::tuple<double,double,double>& real_pose,
  const std::vector<double>& virtual_scan,
  const std::tuple<double,double,double>& virtual_pose,
  const std::string& dump_filepath)
{
  std::vector< std::pair<double,double> > real_scan_points;
  Utils::scan2points(real_scan, real_pose, &real_scan_points);

  std::vector< std::pair<double,double> > virtual_scan_points;
  Utils::scan2points(virtual_scan, virtual_pose, &virtual_scan_points);

  std::ofstream file(dump_filepath.c_str(), std::ios::trunc);

  if (file.is_open())
  {
    file << "rx = [];" << std::endl;
    file << "ry = [];" << std::endl;

    for (int i = 0; i < real_scan.size(); i++)
    {
      file << "rx = [rx " << real_scan_points[i].first << "];" << std::endl;
      file << "ry = [ry " << real_scan_points[i].second << "];" << std::endl;
    }

    file << "vx = [];" << std::endl;
    file << "vy = [];" << std::endl;
    for (int i = 0; i < virtual_scan.size(); i++)
    {
      file << "vx = [vx " << virtual_scan_points[i].first << "];" << std::endl;
      file << "vy = [vy " << virtual_scan_points[i].second << "];" << std::endl;
    }

    file << "r00 = [" << std::get<0>(real_pose) <<
      ", " << std::get<1>(real_pose) << "];" << std::endl;
    file << "v00 = [" << std::get<0>(virtual_pose) <<
      ", " << std::get<1>(virtual_pose) << "];" << std::endl;

    file.close();
  }
  else
    printf("Could not log scans\n");
}


/*******************************************************************************
*/
void Dump::points(const std::vector< std::pair<double,double> >& real_points,
  const std::vector< std::pair<double,double> >& virtual_points,
  const unsigned int& id,
  const std::string& dump_filepath)
{
  std::ofstream file(dump_filepath.c_str(), std::ios::trunc);

  if (file.is_open())
  {
    file << "rx = [];" << std::endl;
    file << "ry = [];" << std::endl;
    for (int i = 0; i < real_points.size(); i++)
    {
      file << "rx = [rx " << real_points[i].first << "];" << std::endl;
      file << "ry = [ry " << real_points[i].second << "];" << std::endl;
    }

    file << "vx = [];" << std::endl;
    file << "vy = [];" << std::endl;
    for (int i = 0; i < virtual_points.size(); i++)
    {
      file << "vx = [vx " << virtual_points[i].first << "];" << std::endl;
      file << "vy = [vy " << virtual_points[i].second << "];" << std::endl;
    }

    file.close();
  }
  else
    printf("Could not log points\n");
}


/*******************************************************************************
*/
void Dump::polygon(const Polygon_2& poly, const std::string& dump_filepath)
{
  std::ofstream file(dump_filepath.c_str(), std::ios::trunc);

  if (file.is_open())
  {
    file << "px = [];" << std::endl;
    file << "py = [];" << std::endl;

    for (VertexIterator vi = poly.vertices_begin();
      vi != poly.vertices_end(); vi++)
    {
      file << "px = [px " << vi->x() << "];" << std::endl;
      file << "py = [py " << vi->y() << "];" << std::endl;
    }

    file.close();
  }
  else
    printf("Could not log polygon\n");
}


/*******************************************************************************
*/
void Dump::polygons(const Polygon_2& real_poly, const Polygon_2& virtual_poly,
  const std::string& dump_filepath)
{
  std::ofstream file(dump_filepath.c_str(), std::ios::trunc);

  if (file.is_open())
  {
    file << "p_rx = [];" << std::endl;
    file << "p_ry = [];" << std::endl;

    for (VertexIterator vi = real_poly.vertices_begin();
      vi != real_poly.vertices_end(); vi++)
    {
      file << "p_rx = [p_rx " << vi->x() << "];" << std::endl;
      file << "p_ry = [p_ry " << vi->y() << "];" << std::endl;
    }

    file << "p_vx = [];" << std::endl;
    file << "p_vy = [];" << std::endl;

    for (VertexIterator vi = virtual_poly.vertices_begin();
      vi != virtual_poly.vertices_end(); vi++)
    {
      file << "p_vx = [p_vx " << vi->x() << "];" << std::endl;
      file << "p_vy = [p_vy " << vi->y() << "];" << std::endl;
    }

    file.close();
  }
  else
    printf("Could not log polygons \n");
}


/*******************************************************************************
*/
void Dump::convexHulls(const std::vector<Point_2>& real_hull,
  const std::vector<Point_2>& virtual_hull,
  const std::string& dump_filepath)
{
  std::ofstream file(dump_filepath.c_str(), std::ios::trunc);

  if (file.is_open())
  {
    file << "h_rx = [];" << std::endl;
    file << "h_ry = [];" << std::endl;

    for (int i = 0; i < real_hull.size(); i++)
    {
      file << "h_rx = [h_rx " << real_hull[i].x() << "];" << std::endl;
      file << "h_ry = [h_ry " << real_hull[i].y() << "];" << std::endl;
    }

    file << "h_vx = [];" << std::endl;
    file << "h_vy = [];" << std::endl;

    for (int i = 0; i < virtual_hull.size(); i++)
    {
      file << "h_vx = [h_vx " << virtual_hull[i].x() << "];" << std::endl;
      file << "h_vy = [h_vy " << virtual_hull[i].y() << "];" << std::endl;
    }

    file.close();
  }
  else
    printf("Could not log hulls \n");
}
