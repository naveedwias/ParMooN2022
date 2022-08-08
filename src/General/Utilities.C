#include <Utilities.h>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <unistd.h>
#include <MooNMD_Io.h>

std::string utilities::get_host_name()
{
  char buf[80];
  gethostname(buf,80);
  return std::string(buf);
}


std::string utilities::get_date_and_time()
{
  auto time = std::chrono::system_clock::now();
  std::time_t t = std::chrono::system_clock::to_time_t(time);
  return std::ctime(&t);
}

double utilities::minmod(const std::vector< double >& input_vec)
{
  if (input_vec.empty())
  {
    ErrThrow("The input vector input_vec is empty. Verify that it is filled ",
        "with at least one entry.");
  }
  else
  {
    double return_val = (input_vec[0] != 0) ? input_vec[0] : 0;
    if (input_vec[0] != 0 && input_vec.size() > 1)
    {
      // if one entry is zero minmod has to return 0 directly.
      double sign = input_vec[0] / std::abs(input_vec[0]);
      for (auto elem : input_vec)
      {
        if (elem == input_vec[0])
        {
          // the first entry or entries that have been checked before do not
          // have to be tested again
          continue;
        }
        if (elem == 0)
        {
          // if one entry is zero minmod has to return 0 directly.
          return_val = 0;
          break;
        }
        else if ( (sign > 0 && elem < 0) || (sign < 0 && elem > 0) )
        {
          // signs differ
          return_val = 0;
          break;
        }
        else if (std::abs(elem) < sign * return_val)
        {
          // new "minimal" value
          return_val = elem;
        }
      }
    }
    return return_val;
  }
}

bool utilities::are_equal(double a, double b, double tol)
{
  bool are_equal = (std::abs(a) < 1) ? std::abs(a-b) < tol :
    std::abs(a-b) < tol * std::abs(a);
  return are_equal;
}
