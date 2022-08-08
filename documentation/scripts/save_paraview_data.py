//! [block0] 
## Python script to save data from Paraview
import os
import sys
#### import the simple module from the paraview
from paraview.simple import *
//! [block0] 

//! [block5] 
def make_directory(dir_):
# check if directories exist and if not create it.
  if not os.path.exists(dir_):
      os.makedirs(dir_)
  else:
      print("A directory " + dir_ + " already exists!")
      exit(0)
//! [block5] 
      
      
//! [block8]       
def clear_pipeline() :
    srcs = GetSources()
    for key, val in srcs.items() :
        #print( "key = " + str(key) + ", value = " + str(val) )
        #print( "   (" + str(type(key)) + "),  (" + str(type(val)) + ")" )
        Delete(val)
        del val
    return

def ResetSession():
    # RESET A PARAVIEW SESSION
    clear_pipeline()
    pxm = servermanager.ProxyManager()
    pxm.UnRegisterProxies()
    del pxm
    Disconnect()
    Connect()
//! [block8] 

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------
file_counter = 0
## LOOP OVER SIMULATIONS
for p in range(0,len(parameters)):
    //! [block9] 
    ResetSession()
    //! [block9] 
    
    //! [block1] 
    File = []
    output_directory = "/home/myname/my/vtk/output/dir"
    parameters = ["m1", "m2", "m3"]
    ntimesteps = [60, 100, 50]
    p = 0
    for i in range(0,ntimesteps[p]):
        index = '%05d' % i  
        readin_filepath = output_directory + "/simulation_" + str(parameters[p]) + "." + index + ".pvtu" 
        File.append(readin_filepath)
    //! [block1] 

    //! [block2]  
    # create a new 'XML Partitioned Unstructured Grid Reader'
    results = XMLPartitionedUnstructuredGridReader(FileName = File)
    results.CellArrayStatus = ['SubDomain', 'RegionID']
    results.PointArrayStatus = ['u', 'u0', 'u1', 'u2', 'p']
    SetActiveSource(results)
    //! [block2]  
    
    //! [block3] 
    # create a new 'Plot Over Line'
    tsteps = []
    tsteps = results.TimestepValues
    plotOverLine1 = PlotOverLine(Input=results,Source='High Resolution Line Source')
    plotOverLine1.Tolerance = 2.22044604925031e-16
    plotOverLine1.Source.Point1 = [-0.105, -0.105, 0.0]
    plotOverLine1.Source.Point2 = [-0.105, -0.105, 0.649999976158142]
    source = plotOverLine1
    //! [block3] 

    //! [block4] 
    results_directory = "./my/extracted/data/dir/simulation_"+ str(parameters[p]) + "/"
    make_directory(results_directory)
    //! [block4] 

    //! [block6]
        writer = CreateWriter(results_directory + "velocity_timestep.csv", source)
        writer.FieldAssociation = "Points"
        writer.WriteAllTimeSteps = 1
        writer.UpdatePipeline()
        del writer
    //! [block6]
    
    //! [block7]
    Delete(plotOverLine1)
    Delete(results)
    del tsteps
    del plotOverLine1
    del results
    //! [block7]

    file_counter = file_counter + 1
