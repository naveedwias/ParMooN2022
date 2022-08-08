/*
 * Read_write_metis.C
 *
 *  Created on: Mar 13, 2017
 *      Author: rominger
 */
#ifdef _MPI

#include <Domain.h>
#include <MeshPartitionInOut.h>
#include <MooNMD_Io.h>

#include "mpi.h"



void MeshPartitionInOut::read_file(const TDomain& Domain,int size,
                                   std::vector<idx_t>& Cell_Rank)
{
    int N_Cells = Cell_Rank.size();
    std::string file_name = Domain.get_database()["read_metis_file"];

    Output::info("PARTITIONING", "Reading mesh partition information "
         "from file ",file_name);

	  std::ifstream ifs; //Deklarieren von ifs zum Oeffnen der Textdatei
	  ifs.open(file_name); //Textdatei wird geoeffnet
	  if(ifs.is_open())
	      {
		      std::string line;
		      int l =0; //a line counter
		      while (std::getline(ifs, line)) //Reading line by line
		      {
		    	  if(l==0)//we expect this line to hold "ParMooN Mesh Partition File"
		    	  {
		    		  if (line!="ParMooN Mesh Partition File")
		    		  {
		    		    throw std::runtime_error("The first line of the "
		    		        "written-file does not match.");
		    		  }
		    		  ++l;
		    		  continue;
		    	  }
		    	  else if(l==1)
		    		  // Reads the total Number of Cells and the Number of Partitions
		    		  //  of the written file and compares them to current configuration
		    		  // throws error if they don't match
		    	  {
		    		  int position = line.find_first_of(','); //finds the ',' in file
		    		  std::string number_partition = line.substr(position+1); //writes 'N_Processors: ...' into string
		    		  std::string number_cells = line.substr(0, position); //writes 'N_Cells: ...' into string
		    		  number_partition.erase (0,14); // Erases the first 14 characters, so the string is only the total number of Processors
		    		  number_cells.erase (0,9); //Erases the first 9 characters, so the string is only the total number of Cells
		    		  int total_partition_number;
		    		  try
		    		  {
		    			  total_partition_number=std::stoi(number_partition); //Turns string of partition number into integer
		    			  if (total_partition_number != size)
		    			  {
		    				  throw std::runtime_error("The total number of partitions of "
		    				      "the written file does not match the total number of "
		    				      "partitions that is used");
		    			  }
		    		  }
		    		  catch (std::exception& e)
		    		  {
		    			  Output::print(e.what());
		    		  }
		    		  int total_number_of_cells;
		    		  try
		    		  {
		    			  total_number_of_cells=std::stoi(number_cells);//Turns string of cell number into integer
		    			  if (total_number_of_cells!=N_Cells)
		    			  {
		    				  throw std::runtime_error("The total number of cells of the "
		    				      "written file does not match the total number of "
		    				      "cells that are used.");
		    			  }
		    		  }
		    		  catch (std::exception& e)
		    		  {
		    			  Output::print(e.what());
		    		  }
		    		  ++l;
		    		  continue;
		    	  }
		    	  else
		    	  {
		    	           //Reads in the partitioning from the read file and checks if the file has the right size and possible values for the partitioning
		    		  int Comment = line.find_first_of('#'); //Looks for lines that start with #
		    		  if (Comment == 0){cout<<line<<endl;continue;} // If a Line starts with # it prints the comment
		    		  else
		    			  //Here we write out the Cellnumber and the matching Cell_Rank.
		    			  //If the partition number or number of cells is too big, than an error is thrown
		    		  {
		    			  int pos = line.find_first_of(',');//finds the ',' in file
		    		  	  std::string Rank = line.substr(pos+1);//writes 'Rank: ...' into string
		    		  	  std::string Cell = line.substr(0, pos);//writes 'Cells: ...' into string
		    		  	  Rank.erase (0,5);//Erases the first 5 characters, so the string is only the partition number
		    		  	  Cell.erase (0,5);//Erases the first 5 characters, so the string is only the cell number
		    		  	  int Rank_Number;
		    		  	  try
		    		  	  {
		    		  		  Rank_Number=std::stoi(Rank);//Turns string of rank into integer
		    		  	  }
		    		  	  catch (std::exception& e)
		    		  	  {
		    		  		  ErrThrow("std::stoi had a problem converting 'Rank': ",Rank,
		    		  				  ". Here is the error message:", e.what());
		    		  	  }
		    		  	  if (Rank_Number >= size)
		    		  	  {
		    		  		  ErrThrow("The input file contains a rank greater than mpi_size.");
		    		  	  }


		    		  	  int Cellnumber = 0;
		    		  	  try
		    		  	  {
		    		  		  Cellnumber=std::stoi(Cell);//Turns cell number of rank into integer
		    		  	  }
		    		  	  catch (std::exception& e)
		    		  	  {
		    		  		  ErrThrow("std::stoi had a problem converting 'Cell': ",Cell,
		    		  				  ". Here is the error message:", e.what());
		    		  	  }

	    		  		  if (Cellnumber>=N_Cells)
	    		  		  {
	    		  			ErrThrow("The input file contains a cell number greater than the maximal cell number.");
	    		  		  }
		    		  	  Cell_Rank[l-2]=Rank_Number;

		    			  if (l > N_Cells + 1){ErrThrow("The file has more lines than it should.");}
		    			  ++l;
		    		  }
		    	  }
	      }
		      if (l-2!=N_Cells) //After the while loop l is the total number of read-in lines.
		      {
		    	Output::print(l, " ", N_Cells);
		    	ErrThrow("The file has less lines than it should.");
		      }
	      }
	  else
	  {
	     ErrThrow("The file was not read. ");
	  }

}
void MeshPartitionInOut::write_file(const TDomain& Domain, int size,
                                    const std::vector<idx_t>& Cell_Rank)
{
    int N_Cells = Cell_Rank.size();
	  std::string file_name = Domain.get_database()["write_metis_file"];

    Output::info("MeshPartitionInOut", "Writing mesh partition information "
        "to file ",file_name);

	  std::ofstream ofs (file_name); //Create write-txt-file
	  //first two lines are headlines and give the total number of Cells and Processors
	  ofs << "ParMooN Mesh Partition File\n" << "N_Cells: "<< N_Cells << " ,N_Processors: " << size <<"\n";

	  for(int i =0; i < N_Cells ; i++) // Loop cell entries of cell rank
	  {
	    //Write entries of each Cell of Cell_Rank in seperate line in txt-file
	    ofs <<"Cell "<< i<< " ,Rank "<< Cell_Rank[i] <<endl;
	  }
	  ofs.close(); //closes write in file
	}

#endif
