#ifndef __POINTWISE_ASSEMBLY_DATA__
#define __POINTWISE_ASSEMBLY_DATA__

#include <vector>

class PointwiseAssemblyData
{
public:
  PointwiseAssemblyData();
  PointwiseAssemblyData(const PointwiseAssemblyData& other);

  void SetCurrentCell(int cell_id, int n_points, int n_data_per_point);
  void SetCurrentPoint(int point_id);

  void Advance();

  static const double* GetOldData();
  static double* GetCurrentData();

private:
  std::vector<std::vector<double>> data_arrays_per_cell;

  int n_points;
  int n_data_per_point;
  int cell_id;

  bool swap_arrays;

  static double* old_data;
  static double* current_data;
};

#endif