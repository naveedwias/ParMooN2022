#include <cmath> // std::isnan
#include <AnalyticalFunction.h>
#include <MooNMD_Io.h>

using namespace parmoon;

/* ************************************************************************* */
AnalyticalFunction::AnalyticalFunction()
  : depends_on_time(false), function(), constant(true)
{
}

/* ************************************************************************* */
AnalyticalFunction::AnalyticalFunction(double a)
  : depends_on_time(false), 
    function(std::function<void(const Point&, double time,
                                FunctionEvaluation&)>(
        [a](const Point&, double, FunctionEvaluation& v)
        {
          v.reset();
          v.set(a, 0,
                0); // value of this function, all derivatives remain zero.
        })),
    constant(true)
{
}

/* ************************************************************************* */
AnalyticalFunction::AnalyticalFunction(
    std::function<void(const Point&, FunctionEvaluation&)> f)
  : depends_on_time(false), 
    function([f](const Point& p, double, FunctionEvaluation& v)
             { return f(p,v); }),
    constant(false)
{
}
/* ************************************************************************* */
AnalyticalFunction::AnalyticalFunction(
    std::function<void(const Point&, double, FunctionEvaluation&)> f)
  : function(f), constant(false)
{
}

/* ************************************************************************* */
AnalyticalFunction::AnalyticalFunction(const AnalyticalFunction& other,
                                       double time)
 : depends_on_time(false), function(), 
   constant(other.is_constant())
{
  if(other.exists())
  {
    function = [other, time](const Point& point, double, FunctionEvaluation& v)
               { other.function(point, time, v); };
  }
}


/* ************************************************************************* */
void AnalyticalFunction::get(const Point& point, FunctionEvaluation& v) const
{
  if(v.space_dimension() != point.dimension())
  {
    ErrThrow("FunctionEvaluation object not suitable to evaluate at this "
             "point, dimension mismatch");
  }
  if(this->depends_on_time)
  {
    ErrThrow("unable to evaluate an AnalyticalFunction without a point in "
             "time");
  }
  if(this->function) // callable
  {
    this->function(point, 0., v); // call at some arbitrary time t = 0
  }
  else
  {
    v.reset();
  }
}

/* ************************************************************************* */
void AnalyticalFunction::get(const Point& point, double t,
                             FunctionEvaluation& v) const
{
  if(v.space_dimension() != point.dimension())
  {
    ErrThrow("FunctionEvaluation object not suitable");
  }
  if(this->function) // callable
  {
    this->function(point, t, v);
  }
  else
  {
    v.reset();
  }
}
