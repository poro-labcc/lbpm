#include "IO/HDF5_IO.h"
#include "IO/IOHelpers.h"
#include "IO/MeshDatabase.h"
#include "IO/Writer.h"
#include "IO/silo.h"
#include "IO/xmlvtk.h"
#include "common/MPI.h"
#include "common/Utilities.h"

#include <algorithm>
#include <memory>
#include <set>
#include <sys/stat.h>
#include <vector>

// std::vector<IO::MeshDatabase> writeMeshesVTI( const std::vector<IO::MeshDataStruct> &meshData,
//     const std::string &path, int rank )
// {
//     const IO::DomainMesh &mesh = dynamic_cast<IO::DomainMesh &>( *meshData.mesh );
//     // RankInfoStruct info( mesh.rank, mesh.nprocx, mesh.nprocy, mesh.nprocz );
        
//     char filename[100], fullpath[200];
//     sprintf( filename, "%05i.vti", rank );
//     sprintf( fullpath, "%s/%s", path.c_str(), filename );
    // VTIWriter vti( std::string(fullpath) );

    // vti.setWholeExtent( info.ix * info.nx, info.jy * info.ny , info.kz * info.nz, 
    //                   ( info.ix + 1 ) * info.nx, ( info.jy + 1 ) * info.ny, ( info.kz + 1 ) * info.nz);
    // vti.setSpacing( 1.0 , 1.0 , 1.0 );
    // vti.setOrigin(0,0,0);
    // vti.setCompress();

    // for ( size_t i = 0; i < meshData.vars.size(); i++ ) 
    // {
    //     const auto &var = *meshData.vars[i];
    //     if ( var.precision == IO::DataType::Double ) {
    //         vti.addCellData( var.name , "Float64" , "binary" , var.dim , (unsigned char*) var.data() );
    //     } else if ( var.precision == IO::DataType::Float ) {
    //         Array<float> data2( var.data.size() );
    //         data2.copy( var.data );
    //         vti.addCellData( var.name , "Float32" , "binary" , var.dim , (unsigned char*) var.data() );
    //     } else if ( var.precision == IO::DataType::Int ) {
    //         Array<int> data2( var.data.size() );
    //         data2.copy( var.data );
    //         vti.addCellData( var.name , "Int32" , "binary" , var.dim , (unsigned char*) var.data() );
    //     } else {
    //         ERROR( "Unsupported format" );
    //     }
    // }

    // vti.write();
//     return std::vector<IO::MeshDatabase>();
// }


// Write a DomainMesh mesh (and variables) to a file
static void writeVtiDomainMesh(
    const std::string &fullpath, const IO::MeshDataStruct &meshData, IO::MeshDatabase database )
{
    const IO::DomainMesh &mesh = dynamic_cast<IO::DomainMesh &>( *meshData.mesh );
    RankInfoStruct info( mesh.rank, mesh.nprocx, mesh.nprocy, mesh.nprocz );

    VTIWriter vti = VTIWriter(std::string(fullpath));
    vti.setWholeExtent( info.ix * mesh.nx, info.jy * mesh.ny , info.kz * mesh.nz, 
                      ( info.ix + 1 ) * mesh.nx, ( info.jy + 1 ) * mesh.ny, ( info.kz + 1 ) * mesh.nz);    
    vti.setSpacing( 1.0 , 1.0 , 1.0 );
    vti.setOrigin(0,0,0);
    vti.setCompress();                      

    for ( size_t i = 0; i < meshData.vars.size(); i++ ) 
    {
        const auto &var = *meshData.vars[i];
        if ( var.precision == IO::DataType::Double ) {
            vti.addCellData( var.name , "Float64" , "binary" , var.dim , (unsigned char*) var.data.begin() );
        } else if ( var.precision == IO::DataType::Float ) {
            Array<float> data2( var.data.size() );
            data2.copy( var.data );
            vti.addCellData( var.name , "Float32" , "binary" , var.dim , (unsigned char*)  data2.begin());
        } else if ( var.precision == IO::DataType::Int ) {
            Array<int> data2( var.data.size() );
            data2.copy( var.data );
            vti.addCellData( var.name , "Int32" , "binary" , var.dim , (unsigned char*) var.data.begin() );
        } else {
            ERROR( "Unsupported format" );
        }
    }

    vti.write();
}
// Write a mesh (and variables) to a file
static IO::MeshDatabase write_domain_vti( const std::string &filename,
    const IO::MeshDataStruct &mesh, IO::FileFormat format, int rank )
{
    // Create the MeshDatabase
    auto database = getDatabase( filename, mesh, format, rank );
    if ( database.meshClass == "DomainMesh" ) {
        writeVtiDomainMesh( filename, mesh, database );
    } else {
        ERROR( "Unknown mesh class or not implemented for vtk/vti output" );
    }
    return database;
}

std::vector<IO::MeshDatabase> writeMeshesVti( const std::vector<IO::MeshDataStruct> &meshData,
    const std::string &path, IO::FileFormat format, int rank )
{
    std::vector<IO::MeshDatabase> meshes_written;
    char filename[100], fullpath[200];
    sprintf( filename, "%05i.vti", rank );
    sprintf( fullpath, "%s/%s", path.c_str(), filename );
    for ( size_t i = 0; i < meshData.size(); i++ ) {
        auto mesh = meshData[i].mesh;
        meshes_written.push_back( write_domain_vti( fullpath , meshData[i], format, rank ) );
    }
    return meshes_written;
}

// void writeVtiSummary(
//     const std::vector<IO::MeshDatabase> &meshes_written, const std::string &filename )
//     {
//         PVTIWriter pvti = PVTIWriter( filename );

//         for ( const auto &data : meshes_written ) 
//         {
//             for ( const auto &tmp : data.domains )
//             {
//                 VTIWriter vti = VTIWriter(  tmp.file  );
//                 pvti.addVTIWriter(vti);
//             }
//         }
//             pvti.write();
            
//     }


void writeVtiSummary(
    const std::vector<IO::MeshDataStruct> &meshData, const std::string &filename )
    {
        PVTIWriter pvti = PVTIWriter( filename );

        pvti.write();
            
    }

