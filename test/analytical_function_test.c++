#include <MooNMD_Io.h>
#include <cmath>
#include <AnalyticalFunction.h>

int main(int, char **)
{
  // empty object to indicate no such function exists.
  AnalyticalFunction f;
  
  if(f.exists())
  {
    ErrThrow("empty f claims it is a real AnalyticalFunction");
  }
  // some arbitrary point in (2D) space
  Point p(1.0, 2.0);
  FunctionEvaluation f_eval(2, 1);
  f.get(p, f_eval);
  auto error_message = "empty AnalyticalFunction returns nonzero value ";
  if(f_eval(0, MultiIndex::D00) != 0.) ErrThrow(error_message, MultiIndex::D00);
  if(f_eval(0, MultiIndex::D10) != 0.) ErrThrow(error_message, MultiIndex::D10);
  if(f_eval(0, MultiIndex::D01) != 0.) ErrThrow(error_message, MultiIndex::D01);
  if(f_eval(0, MultiIndex::D20) != 0.) ErrThrow(error_message, MultiIndex::D20);
  if(f_eval(0, MultiIndex::D11) != 0.) ErrThrow(error_message, MultiIndex::D11);
  if(f_eval(0, MultiIndex::D02) != 0.) ErrThrow(error_message, MultiIndex::D02);
  
  f_eval = FunctionEvaluation(3, 1);
  p = Point(1.0, 2.0, 3.0);
  f.get(p, 1.23, f_eval); // at some point in time
  if(f_eval(0, MultiIndex::D000)!=0.) ErrThrow(error_message, MultiIndex::D000);
  if(f_eval(0, MultiIndex::D100)!=0.) ErrThrow(error_message, MultiIndex::D100);
  if(f_eval(0, MultiIndex::D010)!=0.) ErrThrow(error_message, MultiIndex::D010);
  if(f_eval(0, MultiIndex::D001)!=0.) ErrThrow(error_message, MultiIndex::D001);
  if(f_eval(0, MultiIndex::D200)!=0.) ErrThrow(error_message, MultiIndex::D200);
  if(f_eval(0, MultiIndex::D110)!=0.) ErrThrow(error_message, MultiIndex::D110);
  if(f_eval(0, MultiIndex::D101)!=0.) ErrThrow(error_message, MultiIndex::D101);
  if(f_eval(0, MultiIndex::D020)!=0.) ErrThrow(error_message, MultiIndex::D020);
  if(f_eval(0, MultiIndex::D011)!=0.) ErrThrow(error_message, MultiIndex::D011);
  if(f_eval(0, MultiIndex::D002)!=0.) ErrThrow(error_message, MultiIndex::D002);
  
  // try to evaluate at a 2D point
  try
  {
    f.get(Point(1.0,2.0), f_eval);
    Output::warn<>("AnalyticalFunction::get", "I was able to evaluate a 3D "
                   "AnalyticalFunction at a 2D point without an exception");
    return 1;
  }
  catch(std::runtime_error err)
  {
    // this was expected
  }
  catch(...)
  {
    Output::warn<>("caught wrong exception");
    return 1;
  }
  
  auto f_generator = [](const Point& point, FunctionEvaluation& v)
    {
      const double x = point.x;
      const double y = point.y;
      const double z = point.z;
      v.set(x*y*z, 0, MultiIndex::D000);
      v.set(  y*z, 0, MultiIndex::D100);
      v.set(x*  z, 0, MultiIndex::D010);
      v.set(x*y  , 0, MultiIndex::D001);
      v.set(0.   , 0, MultiIndex::D200);
      v.set(0.   , 0, MultiIndex::D110);
      v.set(0.   , 0, MultiIndex::D101);
      v.set(0.   , 0, MultiIndex::D020);
      v.set(0.   , 0, MultiIndex::D011);
      v.set(0.   , 0, MultiIndex::D002);
    };
  f = AnalyticalFunction(f_generator);
  p = Point(1.0, 2.0, 3.0);
  f.get(p, 1.23, f_eval); // at some point in time
  if(f_eval(0, MultiIndex::D000)!=6.) ErrThrow(error_message, MultiIndex::D000);
  if(f_eval(0, MultiIndex::D100)!=6.) ErrThrow(error_message, MultiIndex::D100);
  if(f_eval(0, MultiIndex::D010)!=3.) ErrThrow(error_message, MultiIndex::D010);
  if(f_eval(0, MultiIndex::D001)!=2.) ErrThrow(error_message, MultiIndex::D001);
  if(f_eval(0, MultiIndex::D200)!=0.) ErrThrow(error_message, MultiIndex::D200);
  if(f_eval(0, MultiIndex::D110)!=0.) ErrThrow(error_message, MultiIndex::D110);
  if(f_eval(0, MultiIndex::D101)!=0.) ErrThrow(error_message, MultiIndex::D101);
  if(f_eval(0, MultiIndex::D020)!=0.) ErrThrow(error_message, MultiIndex::D020);
  if(f_eval(0, MultiIndex::D011)!=0.) ErrThrow(error_message, MultiIndex::D011);
  if(f_eval(0, MultiIndex::D002)!=0.) ErrThrow(error_message, MultiIndex::D002);
  
  
  
  auto f_generator2 = [](const Point& point, double t, FunctionEvaluation& v)
    {
      const double x = point.x;
      const double y = point.y;
      const double z = point.z;
      v.set(t*x*y*z, 0, MultiIndex::D000);
      v.set(t*  y*z, 0, MultiIndex::D100);
      v.set(t*x*  z, 0, MultiIndex::D010);
      v.set(t*x*y  , 0, MultiIndex::D001);
      v.set(0.     , 0, MultiIndex::D200);
      v.set(0.     , 0, MultiIndex::D110);
      v.set(0.     , 0, MultiIndex::D101);
      v.set(0.     , 0, MultiIndex::D020);
      v.set(0.     , 0, MultiIndex::D011);
      v.set(0.     , 0, MultiIndex::D002);
    };
  f = AnalyticalFunction(f_generator2);
  p = Point(1.0, 2.0, 3.0);
  double t = 12.34; // some point in time
  f.get(p, t, f_eval);
  if(f_eval(0, MultiIndex::D000)!=t*6.)ErrThrow(error_message,MultiIndex::D000);
  if(f_eval(0, MultiIndex::D100)!=t*6.)ErrThrow(error_message,MultiIndex::D100);
  if(f_eval(0, MultiIndex::D010)!=t*3.)ErrThrow(error_message,MultiIndex::D010);
  if(f_eval(0, MultiIndex::D001)!=t*2.)ErrThrow(error_message,MultiIndex::D001);
  if(f_eval(0, MultiIndex::D200)!=0.)  ErrThrow(error_message,MultiIndex::D200);
  if(f_eval(0, MultiIndex::D110)!=0.)  ErrThrow(error_message,MultiIndex::D110);
  if(f_eval(0, MultiIndex::D101)!=0.)  ErrThrow(error_message,MultiIndex::D101);
  if(f_eval(0, MultiIndex::D020)!=0.)  ErrThrow(error_message,MultiIndex::D020);
  if(f_eval(0, MultiIndex::D011)!=0.)  ErrThrow(error_message,MultiIndex::D011);
  if(f_eval(0, MultiIndex::D002)!=0.) ErrThrow(error_message, MultiIndex::D002);
  
  try
  {
    f.get(p, f_eval);
    Output::warn<>("AnalyticalFunction::get", "I could evaluate a time "
                   "dependent AnalyticalFunction without providing a time");
    return 1;
  }
  catch(std::runtime_error)
  {
    // this was expected
  }
  catch(...)
  {
    Output::warn<>("AnalyticalFunction::get", "wrong exception thrown");
  }
}
