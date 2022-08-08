#include <PointwiseAssemblyData.h>
#include <MooNMD_Io.h>

PointwiseAssemblyData::PointwiseAssemblyData()
: data_arrays_per_cell(),
  n_points(-1), n_data_per_point(-1), cell_id(-1),
  swap_arrays(false)
{
}

PointwiseAssemblyData::PointwiseAssemblyData(const PointwiseAssemblyData& other)
: data_arrays_per_cell(other.data_arrays_per_cell.size()),
  n_points(other.n_points), n_data_per_point(other.n_data_per_point), cell_id(other.cell_id),
  swap_arrays(other.swap_arrays)
{
  for (size_t i = 0; i < data_arrays_per_cell.size(); i++)
  {
    data_arrays_per_cell[i] = other.data_arrays_per_cell[i];
  }
}

void PointwiseAssemblyData::SetCurrentCell(int cell_id, int n_points, int n_data_per_point)
{
  if (this->n_points >= 0 && this->n_points != n_points)
  {
    ErrThrow("Cannot change quadrature formula when using pointwise persistent data!");
  }

  if (this->n_data_per_point >= 0 && this->n_data_per_point != n_data_per_point)
  {
    ErrThrow("Mismatched amount of persistent data per quadrature point!");
  }

  this->n_points = n_points;
  this->n_data_per_point = n_data_per_point;
  this->cell_id = cell_id;

  if (n_data_per_point > 0)
  {
    if (cell_id >= (int)data_arrays_per_cell.size())
    {
      data_arrays_per_cell.resize(cell_id + 1);
    }

    data_arrays_per_cell[cell_id].resize(2 * n_points * n_data_per_point, 0.0);
  }
}

void PointwiseAssemblyData::SetCurrentPoint(int point_id)
{
  if (cell_id < 0 || n_points < 0 || n_data_per_point < 0)
  {
    ErrThrow("Tried to use pointwise persistent data without setup!");
  }

  if (cell_id >= (int)data_arrays_per_cell.size())
  {
    old_data = nullptr;
    current_data = nullptr;
    return;
  }

  auto& cell_data = data_arrays_per_cell[cell_id];

  if (swap_arrays)
  {
    old_data = cell_data.data() + 2 * n_data_per_point * point_id;
    current_data = old_data + n_data_per_point;
  }
  else
  {
    current_data = cell_data.data() + 2 * n_data_per_point * point_id;
    old_data = current_data + n_data_per_point;
  }
}

void PointwiseAssemblyData::Advance()
{
  swap_arrays = !swap_arrays;
}

const double* PointwiseAssemblyData::GetOldData()
{
  return old_data;
}

double* PointwiseAssemblyData::GetCurrentData()
{
  return current_data;
}

double* PointwiseAssemblyData::old_data = nullptr;
double* PointwiseAssemblyData::current_data = nullptr;