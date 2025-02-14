#include "xmlvtk.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <zlib.h>

using namespace std;

bool Element::compress = false;
headerType cacheSize = 1024*1024;
static char b64table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

unsigned char *spc_base64_encode(unsigned char *input, size_t len)
{
	unsigned char *output, *p;
	size_t        i = 0, mod = len % 3;

	p = output = (unsigned char *) malloc(((len / 3) + (mod ? 1 : 0)) * 4 + 1);
	if (!p) return 0;
	while (i < len - mod)
	{
		*p++ = b64table[input[i++] >> 2];
		*p++ = b64table[((input[i - 1] << 4) | (input[i] >> 4)) & 0x3f];
		*p++ = b64table[((input[i] << 2) | (input[i + 1] >> 6)) & 0x3f];
		*p++ = b64table[input[i + 1] & 0x3f];
		i += 2;
	}
	if (!mod)
	{
		*p = 0;
		return output;
	}
	else
	{
		*p++ = b64table[input[i++] >> 2];
		*p++ = b64table[((input[i - 1] << 4) | (input[i] >> 4)) & 0x3f];
		if (mod == 1)
		{
		*p++ = '=';
		*p++ = '=';
		*p = 0;
		return output;
		}
		else
		{
			*p++ = b64table[(input[i] << 2) & 0x3f];
			*p++ = '=';
			*p = 0;
			return output;
		}
	}
}

std::string CellData::header()
{
    return "<CellData>\n";
}

std::string CellData::footer()
{
    return "</CellData>\n";
}

std::string PointData::header()
{
    return "<PointData>\n";
}

std::string PointData::footer()
{
    return "</PointData>\n";
}

std::string AppendData::header()
{
    return "<AppendedData encoding=\"raw\">\n";
}

std::string AppendData::footer()
{
    return "</AppendedData>\n";
}

AppendData::AppendData()
{
    this -> totalSize = 0;
}

unsigned int AppendData::addData(unsigned char* pointer, unsigned int size)
{
    int offset = this->totalSize;
    this->sizeList.push_back(size);
    this->pointerList.push_back(pointer);
    this->totalSize += (size+4);
    return offset;
}

DataArray::DataArray(const std::string&  name,const std::string&  type, const std::string&  format, int components, unsigned char* pointer, uint64_t points)
{
    this -> components = components;
    this -> name = name;
    this -> type = type;
    this -> format = format;
    this -> pointer = pointer;
    this -> points = points;

    if ( (this -> type == "Int8") || (this -> type == "UInt8") ) this -> typeSize = 1;
    else if ( (this -> type == "Int16") || (this -> type == "UInt16") ) this -> typeSize = 2;
    else if ( (this -> type == "Int32") || (this -> type == "UInt32") ||  (this -> type == "Float32") ) this -> typeSize = 4;
    else if ( (this -> type == "Int64") || (this -> type == "UInt64") ||  (this -> type == "Float64") ) this -> typeSize = 8;
    this -> dataSize = components * points;
}

std::string DataArray::header()
{
    std::ostringstream stringStream;
    stringStream << "<DataArray Name=\"" <<  this-> name << "\"" << " type=\"";
    stringStream << this -> type << "\"";

    stringStream << " NumberOfComponents=\"" << this->components << "\"";

    stringStream << " format=\"";
    stringStream << this -> format << "\" ";

    if (this->format == "appended")
    {
        stringStream  << " offset=\"" << this-> offset << "\" ";
        stringStream << "/>";
    }
    else
    {
        stringStream << ">";
    }

    return stringStream.str();
}

std::string DataArray::footer()
{
    return "\n</DataArray>\n";
}

void VTIWriter::setCompress()
{
    this->compress = true;
}

VTIWriter::VTIWriter(const std::string& filename)
{
    this -> filename = filename;
    this -> vtkVersion = 0.1;
    this -> originX = 0;
    this -> originY = 0;
    this -> originZ = 0;
    this -> compress = false;
    if (sizeof(headerType) == 8) this -> headerTypeName = "UInt64";
    else if (sizeof(headerType) == 4) this -> headerTypeName = "UInt32";

    this -> byteOrder = "LittleEndian";
}

void VTIWriter::addPointData(const std::string&  name,const std::string&  type, const std::string&  format, int components, unsigned char* pointer)
{
    DataArray* data = new DataArray( name , type , format , components , (unsigned char*) pointer ,  pd.points  );
    data -> compress = this -> compress;

    if (format == "appended")
    {
        data -> offset = ad.addData( data->pointer , data->dataSize * data->typeSize );
        ad.compress = this->compress;
    }

    pd.addChild( data );
}

void VTIWriter::addCellData(const std::string&  name,const std::string&  type, const std::string&  format, int components, unsigned char* pointer)
{
    DataArray* data = new DataArray( name , type , format , components , (unsigned char*) pointer ,  cd.cells  );
    data -> compress = this -> compress;

    if (format == "appended")
    {
        data -> offset = ad.addData( data->pointer , data->dataSize * data->typeSize );
        ad.compress = this->compress;
    }

    cd.addChild( data );
}




std::string VTIWriter::header()
{
    std::ostringstream stringStream;
    stringStream << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl;
    stringStream << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"" << this -> byteOrder << "\"";
    stringStream << " header_type=\"" << this -> headerTypeName << "\"";
    if ( this -> compress ) stringStream  << " compressor=\"vtkZLibDataCompressor\"";
    stringStream << ">" << endl;

    stringStream << "<ImageData WholeExtent=\"";
    stringStream << " " << wholeMinX << " " << wholeMaxX;
    stringStream << " " << wholeMinY << " " << wholeMaxY;
    stringStream << " " << wholeMinZ << " " << wholeMaxZ;
    stringStream << "\"";
    stringStream << " Origin=\""  << this -> originX << " " << this -> originY << " " << this-> originZ << "\"";
    stringStream << " Spacing=\"" << this-> spaceX << " " << this->spaceY << " " << this->spaceZ << "\">" << endl;

    stringStream << "<Piece Extent= \"";
    stringStream << " " << pieceMinX << " " << pieceMaxX;
    stringStream << " " << pieceMinY << " " << pieceMaxY;
    stringStream << " " << pieceMinZ << " " << pieceMaxZ;
    stringStream << "\"";
    stringStream << ">" << endl;
    return stringStream.str();
}

std::string VTIWriter::footer()
{
    std::ostringstream stringStream;
	stringStream << "</VTKFile>" << endl;
    return stringStream.str();
}

void VTIWriter::setWholeExtent(int64_t minX_,int64_t minY_, int64_t minZ_,int64_t maxX_,int64_t maxY_, int64_t maxZ_)
{
    this -> wholeMinX = minX_;
    this -> wholeMinY = minY_;
    this -> wholeMinZ = minZ_;
    this -> wholeMaxX = maxX_;
    this -> wholeMaxY = maxY_;
    this -> wholeMaxZ = maxZ_;
    setPiece(minX_,minY_,minZ_,maxX_,maxY_, maxZ_);
}

std::ostream& operator<<(std::ostream& os, AppendData& obj)
{
    if (obj.totalSize > 0)
    {
        os << obj.header();
        os << "_";
        for (std::vector<unsigned char*>::size_type k = 0; k < obj.pointerList.size(); k++)
        {
            if (obj.compress == false)
            {
                os.write( (char*) &obj.sizeList[k], sizeof(headerType) );
                os.write( (char*) obj.pointerList[k], obj.sizeList[k]  );
            }
            else
            {
                streampos beginPos = os.tellp();

                char* compressedData = (char*) malloc(cacheSize);
                headerType totalByteSize = obj.sizeList[k];
                headerType numberOfBlocks =  totalByteSize / cacheSize  +  (totalByteSize % cacheSize != 0 );
                headerType* compressedInfo =  (headerType*) malloc( sizeof(headerType) * (3 + numberOfBlocks) );

                compressedInfo[0] = numberOfBlocks;
                compressedInfo[1] = (numberOfBlocks > 1) ?  cacheSize : totalByteSize;
                compressedInfo[2] = (totalByteSize % cacheSize == 0) ? cacheSize : totalByteSize % cacheSize ;

		os.write( (char*) compressedInfo, (3 + numberOfBlocks) *sizeof(headerType) );

                for (u_int64_t n = 1; n <= numberOfBlocks;  n++)
                {
                    uLongf compressedLength = cacheSize;
                    int numberOfBytesInBlock = (n < numberOfBlocks) ? compressedInfo[1] : compressedInfo[2];
                    int ret = compress( (Bytef *) (compressedData), &compressedLength, obj.pointerList[k] + (n-1) * cacheSize, numberOfBytesInBlock );
                    if (ret != Z_OK)
                    {
                        cout << "Erro No: " << ret;
                        exit(3);
                    }
                    compressedInfo[n+2] = compressedLength;

                    os.write( (char*) compressedData, compressedLength );
                }

                streampos endPos = os.tellp();

                os.seekp( beginPos );

                os.write( (char*) compressedInfo, (3 + numberOfBlocks) *sizeof(headerType) );
                
                os.seekp( endPos );

                free( compressedData );
                free( compressedInfo );
            }
        }        
        os << endl;
        os << obj.footer();

    }
    return os;
}

std::ostream& operator<<(std::ostream& os, DataArray& obj)
{
    os << obj.header() << endl;

    headerType totalByteSize = static_cast<headerType>( (obj.dataSize) * (obj.typeSize) );
    headerType count = min<headerType >( 12 - sizeof(headerType) , totalByteSize );
    unsigned char* base64Data;

    if  (obj.format == "binary")
    {
        if (obj.compress == false)
        {
            os << spc_base64_encode( (unsigned char *) &totalByteSize, sizeof(headerType) );
            count = 0;

            int leftOverSize = 0;
            unsigned char leftOverBuffer[3];

            for (; count < totalByteSize ; count += cacheSize)
            {
                unsigned char* pointer = obj.pointer + count;
                size_t size = min<size_t>( cacheSize, totalByteSize - count);
                if (leftOverSize > 0)
                {
                        for (; leftOverSize < 3; leftOverSize++)
                        {
                            if (size > 0)
                            {
                                leftOverBuffer[leftOverSize] = *(pointer++);
                                size--;
                            }
                            else break;
                        }

                        base64Data = spc_base64_encode( leftOverBuffer , leftOverSize );
                        os << base64Data;
                        free(base64Data);
                        leftOverSize = 0;
                }

                leftOverSize = (size % 3);
                size -= leftOverSize;
                for (int i = 0; i < leftOverSize; i++)
                {
                     leftOverBuffer[i] = pointer[size + i];
                }

                base64Data = spc_base64_encode( pointer , size );
                os << base64Data;
                free(base64Data);
            }

            if (leftOverSize > 0)
            {
                base64Data = spc_base64_encode( leftOverBuffer , leftOverSize );           
                os << base64Data;
                free(base64Data);
                leftOverSize = 0;
            }
        }
        else
        {
            streampos beginPos = os.tellp();

            char* compressedData = (char*) malloc(cacheSize + 3);
            int leftOver = 0;

            headerType numberOfBlocks =  totalByteSize / cacheSize  +  (totalByteSize % cacheSize != 0 );
            headerType* compressedInfo =  (headerType*) malloc( sizeof(headerType) * (3 + numberOfBlocks) );

            compressedInfo[0] = numberOfBlocks;
            compressedInfo[1] = (numberOfBlocks > 1) ?  cacheSize : totalByteSize;
            compressedInfo[2] = (totalByteSize % cacheSize == 0) ? cacheSize : totalByteSize % cacheSize ;

            base64Data =  spc_base64_encode( (unsigned char*) compressedInfo , (3 + numberOfBlocks) *sizeof(headerType) );
            os << base64Data;
            free(base64Data);

            for (u_int64_t n = 1; n <= numberOfBlocks;  n++)
            {
                uLongf compressedLength = cacheSize;
                int numberOfBytesInBlock = (n < numberOfBlocks) ? compressedInfo[1] : compressedInfo[2];
                compress( (Bytef *) (compressedData+3), &compressedLength, obj.pointer + (n-1) * cacheSize, numberOfBytesInBlock );
                compressedInfo[n+2] = compressedLength;
                
                int encodeSize = 3* ( (compressedLength + leftOver)/3 );
                base64Data = spc_base64_encode( (unsigned char*) (compressedData + 3 - leftOver) , encodeSize );
                os << base64Data;
                free(base64Data);
                
                int newLeftOver = compressedLength + leftOver - encodeSize ;
                if (newLeftOver > 0) 
                {
                    memcpy( compressedData + 3 - newLeftOver, compressedData + 3 - leftOver + encodeSize, newLeftOver);
                }
                leftOver = newLeftOver;
            }

            if (leftOver > 0)
            {
                base64Data = spc_base64_encode( (unsigned char*) (compressedData + 3 - leftOver) , leftOver );
                os << base64Data;
                free(base64Data);
            }

            free(compressedData);              

            streampos endPos = os.tellp();
            os.seekp( beginPos );

            base64Data =  spc_base64_encode( (unsigned char*) compressedInfo , (3 + numberOfBlocks) *sizeof(headerType) );                      
            os << base64Data;
            free(base64Data);
            free(compressedInfo);
            os.seekp( endPos );
            
        }

        os << obj.footer();
    }
    
    return os;
}

void VTIWriter::setPiece(int64_t minX_,int64_t minY_, int64_t minZ_,int64_t maxX_,int64_t maxY_, int64_t maxZ_)
{
    this -> pieceMinX = minX_;
    this -> pieceMinY = minY_;
    this -> pieceMinZ = minZ_;
    this -> pieceMaxX = maxX_;
    this -> pieceMaxY = maxY_;
    this -> pieceMaxZ = maxZ_;
    this -> sizeX = maxX_ - minX_ + 1;
    this -> sizeY = maxY_ - minY_ + 1;
    this -> sizeZ = maxZ_ - minZ_ + 1;
    this -> pd.points = sizeX * sizeY * sizeZ;
    this -> cd.cells  = (sizeX-1) * (sizeY-1) * (sizeZ-1);
}

void VTIWriter::setOrigin(double x_,double y_, double z_)
{
    this -> originX = x_;
    this -> originY = y_;
    this -> originZ = z_;
}

void VTIWriter::setSpacing(double sx_,double sy_, double sz_)
{
    this -> spaceX = sx_;
    this -> spaceY = sy_;
    this -> spaceZ = sz_;
}

void VTIWriter::write()
{
    file.open(this -> filename);
    file << header();

    if (pd.sizeChild() > 0)
    {
        file << pd.header();
        for (int n = 0; n < pd.sizeChild() ; n++)
        {
            DataArray* array = (DataArray*) pd.getChild(n);
            file << *array;
        }
        file << pd.footer();
    }

    if (cd.sizeChild() > 0)
    {
        file << cd.header();
        for (int n = 0; n < cd.sizeChild() ; n++)
        {
            DataArray* array = (DataArray*) cd.getChild(n);
            file << *array;
        }
        file << cd.footer();
    }

    file << "</Piece>" << endl;
    file << "</ImageData>" << endl;

    file << ad;
    file << footer();
    file.close();
}
