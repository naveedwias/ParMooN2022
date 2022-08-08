#include <MooNMD_Io.h>
#include <cmath>
#include <FunctionEvaluation.h>

int main(int, char**)
{
  unsigned int space_dim = 2;
  unsigned int n_functions = 3;
  FunctionEvaluation f_eval(space_dim, n_functions, 1);
  // sanity checks
  if(f_eval.n_multi_indices() != 6)
    ErrThrow("FunctionEvaluation::n_multi_indices() seems wrong ",
             f_eval.n_multi_indices(), " != 6");
  if(f_eval.space_dimension() != space_dim)
    ErrThrow("FunctionEvaluation::space_dimension() seems wrong ",
             f_eval.space_dimension(), " != ", space_dim);
  if(f_eval.n_functions() != n_functions)
    ErrThrow("FunctionEvaluation::n_functions() seems wrong ",
             f_eval.n_functions(), " != ", n_functions);
  
  // fill with some numbers
  for(unsigned int i = 0; i < n_functions; ++i)
  {
    f_eval.set(1. + i*10., i, MultiIndex::D00);
    f_eval.set(2. + i*10., i, MultiIndex::D10);
    f_eval.set(3. + i*10., i, MultiIndex::D01);
    f_eval.set(4. + i*10., i, MultiIndex::D20);
    f_eval.set(5. + i*10., i, MultiIndex::D11);
    f_eval.set(6. + i*10., i, MultiIndex::D02);
  }
  f_eval.reset();
  // same with the other set method
  for(unsigned int i = 0; i < n_functions; ++i)
  {
    f_eval.set(1. + i*10., i, 0);
    f_eval.set(2. + i*10., i, 1);
    f_eval.set(3. + i*10., i, 2);
    f_eval.set(4. + i*10., i, 3);
    f_eval.set(5. + i*10., i, 4);
    f_eval.set(6. + i*10., i, 5);
  }
  
  // try out wrong MultiIndex
  try
  {
    f_eval.set(0., 0, MultiIndex::D000);
    Output::warn<>("FunctionEvaluation::set with wrong MultiIndex",
                   "no exception thrown ");
    return 1;
  }
  catch(std::runtime_error err)
  {
    // this was expected
  }
  catch(...)
  {
    Output::warn<>("testing FunctionEvaluation::set", "wrong exception thrown");
    return 1;
  }
  
  // test the getters:
  auto error_message = "wrong value in FunctionEvaluation::get ";
  for(unsigned int i = 0; i < n_functions; ++i)
  {
    if( f_eval.get(i, MultiIndex::D00) != 1.+i*10.) ErrThrow(error_message, i);
    if( f_eval.get(i, MultiIndex::D10) != 2.+i*10.) ErrThrow(error_message, i);
    if( f_eval.get(i, MultiIndex::D01) != 3.+i*10.) ErrThrow(error_message, i);
    if( f_eval.get(i, MultiIndex::D20) != 4.+i*10.) ErrThrow(error_message, i);
    if( f_eval.get(i, MultiIndex::D11) != 5.+i*10.) ErrThrow(error_message, i);
    if( f_eval.get(i, MultiIndex::D02) != 6.+i*10.) ErrThrow(error_message, i);
    
    if( f_eval.get(i, 0) != 1.+i*10.) ErrThrow(error_message, i);
    if( f_eval.get(i, 1) != 2.+i*10.) ErrThrow(error_message, i);
    if( f_eval.get(i, 2) != 3.+i*10.) ErrThrow(error_message, i);
    if( f_eval.get(i, 3) != 4.+i*10.) ErrThrow(error_message, i);
    if( f_eval.get(i, 4) != 5.+i*10.) ErrThrow(error_message, i);
    if( f_eval.get(i, 5) != 6.+i*10.) ErrThrow(error_message, i);
    
    if( f_eval(i, MultiIndex::D00) != 1.+i*10.) ErrThrow(error_message, i);
    if( f_eval(i, MultiIndex::D10) != 2.+i*10.) ErrThrow(error_message, i);
    if( f_eval(i, MultiIndex::D01) != 3.+i*10.) ErrThrow(error_message, i);
    if( f_eval(i, MultiIndex::D20) != 4.+i*10.) ErrThrow(error_message, i);
    if( f_eval(i, MultiIndex::D11) != 5.+i*10.) ErrThrow(error_message, i);
    if( f_eval(i, MultiIndex::D02) != 6.+i*10.) ErrThrow(error_message, i);
  }
  // try out wrong MultiIndex
  try
  {
    f_eval.get(0, MultiIndex::D000);
    Output::warn<>("FunctionEvaluation::get with wrong MultiIndex",
                   "no exception thrown ");
    return 1;
  }
  catch(std::runtime_error err)
  {
    // this was expected
  }
  catch(...)
  {
    Output::warn<>("testing FunctionEvaluation::get", "wrong exception thrown");
    return 1;
  }
  
  
  // testing the iterators
  auto begin_D00 = f_eval.begin(MultiIndex::D00);
  auto end_D00 = f_eval.end(MultiIndex::D00);
  Output::print("values of the functions are ");
  for(auto it = begin_D00; it != end_D00; ++it)
    Output::print(*it);
  // try out wrong MultiIndex
  try
  {
    f_eval.get(0, MultiIndex::D000);
    Output::warn<>("FunctionEvaluation::get with wrong MultiIndex",
                   "no exception thrown ");
    return 1;
  }
  catch(std::runtime_error err)
  {
    // this was expected
  }
  catch(...)
  {
    Output::warn<>("testing FunctionEvaluation::get", "wrong exception thrown");
    return 1;
  }
  
  // copy construct
  FunctionEvaluation f_eval2(f_eval);
  
  // call scale only to make sure it doesn't crash
  f_eval.scale(0, -1.);
  
  f_eval2.add(f_eval, -1.);
  for(unsigned int i = 1; i < n_functions; ++i)
  {
    if( f_eval2.get(i, MultiIndex::D00) != 0.) ErrThrow(error_message, i);
    if( f_eval2.get(i, MultiIndex::D10) != 0.) ErrThrow(error_message, i);
    if( f_eval2.get(i, MultiIndex::D01) != 0.) ErrThrow(error_message, i);
    if( f_eval2.get(i, MultiIndex::D20) != 0.) ErrThrow(error_message, i);
    if( f_eval2.get(i, MultiIndex::D11) != 0.) ErrThrow(error_message, i);
    if( f_eval2.get(i, MultiIndex::D02) != 0.) ErrThrow(error_message, i);
  }
  
  // undo the scale from before
  f_eval.scale(0, -1.);
  f_eval.merge();
  if(f_eval.n_functions() != 1)
  {
    ErrThrow("FunctionEvaluation::n_functions() seems wrong after merge",
             f_eval.n_functions(), " != 1");
  }
  double x = 0.5*n_functions*(n_functions-1);
  if(f_eval.get(0, MultiIndex::D00) != n_functions+10*x)
  {
    Output::print("wrong value after merge ", f_eval.get(0, MultiIndex::D00),
                  "  ", x , " ", n_functions+10*x);
    return 1;
  }
  if(f_eval.get(0, MultiIndex::D10) != 2*n_functions+10*x)
  {
    Output::print("wrong D10 after merge ", f_eval.get(0, MultiIndex::D10),
                  "  ", x , " ", 2*n_functions+10*x);
    return 1;
  }
  if(f_eval.get(0, MultiIndex::D01) != 3*n_functions+10*x)
  {
    Output::print("wrong D01 after merge ", f_eval.get(0, MultiIndex::D01),
                  "  ", x , " ", 3*n_functions+10*x);
    return 1;
  }
  if(f_eval.get(0, MultiIndex::D20) != 4*n_functions+10*x)
  {
    Output::print("wrong D20 after merge ", f_eval.get(0, MultiIndex::D20),
                  "  ", x , " ", 4*n_functions+10*x);
    return 1;
  }
  if(f_eval.get(0, MultiIndex::D11) != 5*n_functions+10*x)
  {
    Output::print("wrong D11 after merge ", f_eval.get(0, MultiIndex::D11),
                  "  ", x , " ", 5*n_functions+10*x);
    return 1;
  }
  if(f_eval.get(0, MultiIndex::D02) != 6*n_functions+10*x)
  {
    Output::print("wrong D02 after merge ", f_eval.get(0, MultiIndex::D02),
                  "  ", x , " ", 6*n_functions+10*x);
    return 1;
  }
  
  
  
  
  // a few things in 3D:
  f_eval2 = FunctionEvaluation(3, n_functions);
  try
  {
    f_eval.add(f_eval2, -5.);
    Output::warn<>("FunctionEvaluation::add with wrong FunctionEvaluation",
                   "no exception thrown ");
    return 1;
  }
  catch(std::runtime_error err)
  {
    // this was expected
  }
  catch(...)
  {
    Output::warn<>("testing FunctionEvaluation::add", "wrong exception thrown");
    return 1;
  }
  for(unsigned int i = 0; i < n_functions; ++i)
  {
    f_eval2.set(1. + i*10., i, MultiIndex::D000);
    f_eval2.set(2. + i*10., i, MultiIndex::D100);
    f_eval2.set(3. + i*10., i, MultiIndex::D010);
    f_eval2.set(4. + i*10., i, MultiIndex::D001);
    f_eval2.set(5. + i*10., i, MultiIndex::D200);
    f_eval2.set(6. + i*10., i, MultiIndex::D110);
    f_eval2.set(7. + i*10., i, MultiIndex::D101);
    f_eval2.set(8. + i*10., i, MultiIndex::D020);
    f_eval2.set(9. + i*10., i, MultiIndex::D011);
    f_eval2.set(10.+ i*10., i, MultiIndex::D002);
  }
  
  
}
