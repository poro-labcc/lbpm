/*
 * Pre-processor to generate signed distance function from segmented data
 * segmented data should be stored in a raw binary file as 1-byte integer (type char)
 * will output distance functions for phases
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "common/Array.h"
#include "common/Domain.h"
#include "analysis/distance.h"
#include "analysis/morphology.h"

//*************************************************************************
// Morpohologica pre-processor
//   Initialize phase distribution using morphological approach
//   Signed distance function is used to determine fluid configuration
//*************************************************************************

int main(int argc, char **argv)
{
	// Initialize MPI
    Utilities::startup( argc, argv );
	Utilities::MPI comm( MPI_COMM_WORLD );
    int rank = comm.getRank();
	{
		//.......................................................................
		// Reading the domain information file
		//.......................................................................
		char LocalRankFilename[40];

		string filename;

		if (argc > 1){
			filename=argv[1];
		}
		else ERROR("No input database provided\n");

		auto db = std::make_shared<Database>( filename );
		auto domain_db = db->getDatabase( "Domain" );

		// Read domain parameters
		auto READFILE = domain_db->getScalar<std::string>( "Filename" );
		auto size = domain_db->getVector<int>( "n" );
		auto nproc = domain_db->getVector<int>( "nproc" );
		auto ReadValues = domain_db->getVector<int>( "ReadValues" );
		auto WriteValues = domain_db->getVector<int>( "WriteValues" );

		int nx = size[0];
		int ny = size[1];
		int nz = size[2];

		size_t N = (nx+2)*(ny+2)*(nz+2);

		std::shared_ptr<Domain> Dm (new Domain(domain_db,comm));
		std::shared_ptr<Domain> Mask (new Domain(domain_db,comm));

		for (size_t n=0; n<N; n++) Dm->id[n]=1;
		Dm->CommInit();

		signed char *id;
		id = new signed char [N];
		Mask->Decomp(READFILE);
		Mask->CommInit();
		
		nx+=2; ny+=2; nz+=2;

		Array<char> id_solid(nx,ny,nz);
		IntArray SignDist(nx,ny,nz);

		// Solve for the position of the solid phase
		for (int k=0;k<nz;k++){
			for (int j=0;j<ny;j++){
				for (int i=0;i<nx;i++){
					int n = k*nx*ny+j*nx+i;
					id[n] = Mask->id[n];
					// Initialize the solid phase
					if (Mask->id[n] > 0){
						id_solid(i,j,k) = 1;
					}
					else	    
						id_solid(i,j,k) = 0;
				}
			}
		}

		if (rank==0) printf("Initialized solid phase -- Converting to Signed Distance function \n");
		CalcClassicEDT(SignDist,id_solid,*Dm);

		comm.barrier();

		FILE *ID = fopen("distance.raw","wb");
		fwrite(SignDist.data(),sizeof(int),nx*ny*nz,ID);
		fclose(ID);
		
		// write the geometry to a single file
		for (int k=0;k<nz;k++){
			for (int j=0;j<ny;j++){ 
				for (int i=0;i<nx;i++){
					int n = k*nx*ny+j*nx+i;
					Mask->id[n] = id[n];
				}
			}
		}
		comm.barrier();
	}

    Utilities::shutdown();
}
