/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <PCU.h>
#include <lionBase64.h>
#include <lionCompress.h>
#include "apfMesh.h"
#include "apfNumbering.h"
#include "apfNumberingClass.h"
#include "apfShape.h"
#include "apfFieldData.h"
#include "apfConvert.h"
#include "../mds/apfMDS.h"
#include <sstream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <stdint.h>
#include <vector>

// === includes for safe_mkdir ===
#include <reel.h>
#include <sys/types.h> /*required for mode_t for mkdir on some systems*/
#include <sys/stat.h> /*using POSIX mkdir call for SMB "foo/" path*/
#include <errno.h> /* for checking the error from mkdir */
// ===============================

namespace apf {

class HasAll : public FieldOp
{
  public:
    virtual bool inEntity(MeshEntity* e)
    {
      if (!f->getData()->hasEntity(e))
        ok = false;
      return false;
    }
    bool run(FieldBase* f_)
    {
      f = f_;
      ok = true;
      this->apply(f);
      return ok;
    }
  private:
    bool ok;
    FieldBase* f;
};

static bool shouldPrint(
    FieldBase* f,
    std::vector<std::string> writeFields)
{
  bool print = false;
  std::string fieldName = f->getName();
  for (size_t i=0; i < writeFields.size(); ++i)
    if (fieldName == writeFields[i])
      print = true;
  return print;
}

static bool isPrintable(FieldBase* f)
{
  HasAll op;
  return PCU_And(op.run(f));
}

static bool isNodal(FieldBase* f)
{
  return f->getShape() == f->getMesh()->getShape();
}

static bool isIP(FieldBase* f, int cellDim)
{
  FieldShape* s = f->getShape();
  for (int i=0; i < cellDim; ++i)
    if (s->hasNodesIn(i))
      return false;
  return s->hasNodesIn(cellDim);
}

static void describeArray(
    std::ostream& file,
    const char* name,
    int type,
    int size,
    bool isWritingBinary = false)
{
  file << "type=\"";
  const char* typeNames[3] = {"Float64","Int32","Int64"};
  file << typeNames[type];
  file << "\" Name=\"" << name;
  file << "\" NumberOfComponents=\"" << size;
  if (isWritingBinary)
  {
    file << "\" format=\"binary\"";
  }
  else
  {
    file << "\" format=\"ascii\"";
  }
}

static void writePDataArray(
    std::ostream& file,
    const char* name,
    int type,
    int size,
    bool isWritingBinary = false)
{
  file << "<PDataArray ";
  describeArray(file,name,type,size,isWritingBinary);
  file << "/>\n";
}

static void writePDataArray(
    std::ostream& file,
    FieldBase* f,
    bool isWritingBinary = false)
{
  file << "<PDataArray ";
  describeArray(file,
      f->getName(),
      f->getScalarType(),
      f->countComponents(),
      isWritingBinary);
  file << "/>\n";
}

static void writePPoints(std::ostream& file,
    Field* f,
    bool isWritingBinary = false)
{
  file << "<PPoints>\n";
  writePDataArray(file,f,isWritingBinary);
  file << "</PPoints>\n";
}

static void writePPointData(std::ostream& file,
    Mesh* m,
    std::vector<std::string> writeFields,
    bool isWritingBinary = false)
{
  file << "<PPointData>\n";
  for (int i=0; i < m->countFields(); ++i)
  {
    Field* f = m->getField(i);
    if (isNodal(f) && shouldPrint(f,writeFields))
    {
      writePDataArray(file,f,isWritingBinary);
    }
  }
  for (int i=0; i < m->countNumberings(); ++i)
  {
    Numbering* n = m->getNumbering(i);
    if (isNodal(n) && shouldPrint(n,writeFields))
    {
      writePDataArray(file,n,isWritingBinary);
    }
  }
  for (int i=0; i < m->countGlobalNumberings(); ++i)
  {
    GlobalNumbering* n = m->getGlobalNumbering(i);
    if (isNodal(n) && shouldPrint(n,writeFields))
    {
      writePDataArray(file,n,isWritingBinary);
    }
  }
  file << "</PPointData>\n";
}

static int countIPs(FieldBase* f, int cellDim)
{
  /* assuming a non-mixed mesh here. for the already strained capabilities
     of VTK to accept IP fields, this is the best we can do */
  Mesh* m = f->getMesh();
  MeshIterator* it = m->begin(cellDim);
  MeshEntity* e = m->iterate(it);
  m->end(it);
  return f->countNodesOn(e);
}

static std::string getIPName(FieldBase* f, int point)
{
  std::stringstream ss;
  /* People looking at these files get scared of 0-based indexing
                                     V                        */
  ss << f->getName() << '_' << (point+1);
  return ss.str();
}

static void writeIP_PCellData(std::ostream& file,
    FieldBase* f,
    bool isWritingBinary,
    int cellDim)
{
  int n = countIPs(f, cellDim);
  for (int p=0; p < n; ++p)
  {
    std::string s= getIPName(f,p);
    writePDataArray(file,
        s.c_str(),f->getScalarType(),
        f->countComponents(),
        isWritingBinary);
  }
}

static void writePCellParts(std::ostream& file, bool isWritingBinary = false)
{
  writePDataArray(file, "apf_part", apf::Mesh::INT, 1, isWritingBinary);
}

static void writePCellData(std::ostream& file,
    Mesh* m,
    std::vector<std::string> writeFields,
    bool isWritingBinary,
    int cellDim)
{
  file << "<PCellData>\n";
  for (int i=0; i < m->countFields(); ++i)
  {
    Field* f = m->getField(i);
    if (isIP(f, cellDim) && shouldPrint(f,writeFields))
    {
      writeIP_PCellData(file, f, isWritingBinary, cellDim);
    }
  }
  for (int i=0; i < m->countNumberings(); ++i)
  {
    Numbering* n = m->getNumbering(i);
    if (isIP(n, cellDim) && shouldPrint(n,writeFields))
    {
      writeIP_PCellData(file, n, isWritingBinary, cellDim);
    }
  }
  for (int i=0; i < m->countGlobalNumberings(); ++i)
  {
    GlobalNumbering* n = m->getGlobalNumbering(i);
    if (isIP(n, cellDim) && shouldPrint(n,writeFields))
    {
      writeIP_PCellData(file, n, isWritingBinary, cellDim);
    }
  }
  writePCellParts(file,isWritingBinary);
  file << "</PCellData>\n";
}

static std::string getPieceFileName(int id)
{
  std::stringstream ss;
  ss << id << ".vtu";
  return ss.str();
}

static std::string stripPath(std::string const& s)
{
  size_t i = s.rfind('/');
  if (i == std::string::npos)
  {
    return s;
  }
  return s.substr(i + 1, std::string::npos);
}

static std::string getRelativePathPSource(int id)
{
  std::stringstream ss;
  int dirNum = id/1024;
  ss << dirNum << '/';
  return ss.str();
}

static void writePSources(std::ostream& file)
{
  for (int i=0; i < PCU_Comm_Peers(); ++i)
  {
    std::string fileName = stripPath(getPieceFileName(i));
    std::string fileNameAndPath = getRelativePathPSource(i) + fileName;
    file << "<Piece Source=\"" << fileNameAndPath << "\"/>\n";
  }
}

static void writePvtuFile(const char* prefix,
    Mesh* m,
    std::vector<std::string> writeFields,
    bool isWritingBinary,
    int cellDim)
{
  std::string fileName = stripPath(prefix);
  fileName += ".pvtu";
  std::stringstream ss;
  ss << prefix << '/' << fileName;
  std::string fileNameAndPath = ss.str();
  std::ofstream file(fileNameAndPath.c_str());
  assert(file.is_open());
  file << "<VTKFile type=\"PUnstructuredGrid\">\n";
  file << "<PUnstructuredGrid GhostLevel=\"0\">\n";
  writePPoints(file,m->getCoordinateField(),isWritingBinary);
  writePPointData(file,m,writeFields,isWritingBinary);
  writePCellData(file, m, writeFields, isWritingBinary, cellDim);
  writePSources(file);
  file << "</PUnstructuredGrid>\n";
  file << "</VTKFile>\n";
}

static void readPvtuFile(const char* prefix, bool& isWritingBinary)
{
  std::string fileName = stripPath(prefix);
  fileName += ".pvtu";
  std::stringstream ss;
  ss << prefix << '/' << fileName;
  std::string fileNameAndPath = ss.str();
  std::ifstream file(fileNameAndPath.c_str());
  assert(file.is_open());

  std::string line;
  while (std::getline(file, line)) {
      if (line.find("PDataArray") != std::string::npos) {
          if (line.find("binary") != std::string::npos)
              isWritingBinary = true;
          break;
      }
  }
  file.close();
}

static void writeDataHeader(std::ostream& file,
    const char* name,
    int type,
    int size,
    bool isWritingBinary = false)
{
  file << "<DataArray ";
  describeArray(file,name,type,size,isWritingBinary);
  file << ">\n";
}

static void writeEncodedArray(std::ostream& file,
    unsigned int dataLenBytes,
    char* dataToEncode)
{
  if ( lion::can_compress )
  {
    //build the data header and compress dataToEncode
    long lensToEncode[4];
    lensToEncode[0] = 1; //data is compressed in one block
    lensToEncode[1] = dataLenBytes; //size of each block before compression
    lensToEncode[2] = dataLenBytes; //size of the final block before compression
    unsigned long dataCompressedLen = lion::compressBound(dataLenBytes);
        //gets the size to allocate for the compressed data array
    char* dataCompressed = new char[dataCompressedLen]();
    lion::compress( (void*)dataCompressed,
        dataCompressedLen,
        (void*)dataToEncode,
        dataLenBytes);
        //compresses and sets dataCompressedLen to the size of the correct length
    lensToEncode[3] = dataCompressedLen; //size of the compressed block
    //encode and output the compressed data
    file << lion::base64Encode( (char*)lensToEncode, sizeof(lensToEncode) );
    file << lion::base64Encode( dataCompressed, dataCompressedLen ) << '\n';
    delete[] dataCompressed;
  }
  else
  {
    //not compressing, encode and output
    file << lion::base64Encode( (char*)&dataLenBytes, sizeof(dataLenBytes) );
    file << lion::base64Encode( dataToEncode, dataLenBytes ) << '\n';
  }
}

/* Paraview/VTK has trouble with sub-normal double precision floating point
 * ASCII values.
 *
 * http://www.paraview.org/Bug/view.php?id=15925
 *
 * This function exists to cast "double" to "float"
 * before writing it to file, the others are to maintain the
 * templated design of writeNodalField and others */

static float workaround(double v)
{
  return static_cast<float>(v);
}
static int workaround(int v)
{
  return v;
}
static long workaround(long v)
{
  return v;
}

template <class T>
static void writeNodalField(std::ostream& file,
    FieldBase* f,
    DynamicArray<Node>& nodes,
    bool isWritingBinary = false)
{
  int nc = f->countComponents();
  writeDataHeader(file,f->getName(),f->getScalarType(),nc,isWritingBinary);
  NewArray<T> nodalData(nc);
  FieldDataOf<T>* data = static_cast<FieldDataOf<T>*>(f->getData());
  if (isWritingBinary)
  {
    unsigned int dataLen = nc * nodes.getSize();
    unsigned int dataLenBytes = dataLen*sizeof(T);
    T* dataToEncode = new T[dataLen]();
    unsigned int dataIndex = 0;
    for (size_t i = 0; i < nodes.getSize(); ++i)
    {
      data->getNodeComponents(nodes[i].entity,nodes[i].node,&(nodalData[0]));
      for (int j = 0; j < nc; ++j)
      {
        dataToEncode[dataIndex] = nodalData[j];
        dataIndex++;
      }
    }
    writeEncodedArray(file, dataLenBytes, (char*)dataToEncode);
    delete [] dataToEncode;
  }
  else
  {
    for (size_t i = 0; i < nodes.getSize(); ++i)
    {
      data->getNodeComponents(nodes[i].entity,nodes[i].node,&(nodalData[0]));
      for (int j = 0; j < nc; ++j)
      {
        file << workaround(nodalData[j]) << ' ';
      }
      file << '\n';
    }
  }
  file << "</DataArray>\n";
}

static void writePoints(std::ostream& file,
    Mesh* m,
    DynamicArray<Node>& nodes,
    bool isWritingBinary = false)
{
  file << "<Points>\n";
  writeNodalField<double>(file,m->getCoordinateField(),nodes,isWritingBinary);
  file << "</Points>\n";
}

static int countElementNodes(Numbering* n, MeshEntity* e)
{
  return n->getShape()->getEntityShape(n->getMesh()->getType(e))->countNodes();
}

static void writeConnectivity(std::ostream& file,
    Numbering* n,
    bool isWritingBinary,
    int cellDim)
{
  file << "<DataArray type=\"Int32\" Name=\"connectivity\"";
  if (isWritingBinary)
  {
    file << " format=\"binary\"";
  }
  else
  {
    file << " format=\"ascii\"";
  }
  file << ">\n";
  Mesh* m = n->getMesh();
  MeshEntity* e;
  if (isWritingBinary)
  {
    MeshIterator* elements = m->begin(cellDim);
    unsigned int dataLen = 0;
    while ((e = m->iterate(elements)))
    {
      dataLen += countElementNodes(n,e);
    }
    m->end(elements);
    unsigned int dataLenBytes = dataLen*sizeof(int);
    int* dataToEncode = new int[dataLen]();
    elements = m->begin(cellDim);
    unsigned int dataIndex = 0;
    while ((e = m->iterate(elements)))
    {
      int nen = countElementNodes(n,e);
      NewArray<int> numbers(nen);
      getElementNumbers(n,e,numbers);
      for (int i=0; i < nen; ++i)
      {
        dataToEncode[dataIndex] = numbers[i];
        dataIndex++;
      }
    }
    m->end(elements);
    writeEncodedArray(file, dataLenBytes, (char*)dataToEncode);
    delete [] dataToEncode;
  }
  else
  {
    MeshIterator* elements = m->begin(cellDim);
    while ((e = m->iterate(elements)))
    {
      int nen = countElementNodes(n,e);
      NewArray<int> numbers(nen);
      getElementNumbers(n,e,numbers);
      for (int i=0; i < nen; ++i)
      {
        file << numbers[i] << ' ';
      }
      file << '\n';
    }
    m->end(elements);
  }
  file << "</DataArray>\n";
}

static void writeOffsets(std::ostream& file,
    Numbering* n,
    bool isWritingBinary,
    int cellDim)
{
  file << "<DataArray type=\"Int32\" Name=\"offsets\"";
  if (isWritingBinary)
  {
    file << " format=\"binary\"";
  }
  else
  {
    file << " format=\"ascii\"";
  }
  file << ">\n";
  Mesh* m = n->getMesh();
  MeshEntity* e;
  if (isWritingBinary)
  {
    MeshIterator* elements = m->begin(cellDim);
    unsigned int dataLen = 0;
    while ((e = m->iterate(elements)))
    {
      dataLen++;
    }
    m->end(elements);
    unsigned int dataLenBytes = dataLen*sizeof(int);
    int* dataToEncode = new int[dataLen]();
    elements = m->begin(cellDim);
    unsigned int dataIndex = 0;
    int offset = 0;
    while ((e = m->iterate(elements)))
    {
      offset += countElementNodes(n,e);
      dataToEncode[dataIndex] = offset;
      dataIndex++;
    }
    m->end(elements);
    writeEncodedArray(file, dataLenBytes, (char*)dataToEncode);
    delete [] dataToEncode;
  }
  else
  {
    MeshIterator* elements = m->begin(cellDim);
    int offset = 0;
    while ((e = m->iterate(elements)))
    {
      offset += countElementNodes(n,e);
      file << offset << '\n';
    }
    m->end(elements);
  }
  file << "</DataArray>\n";
}

static void writeTypes(std::ostream& file,
    Mesh* m,
    bool isWritingBinary,
    int cellDim)
{
  file << "<DataArray type=\"UInt8\" Name=\"types\"";
  if (isWritingBinary)
  {
    file << " format=\"binary\"";
  }
  else
  {
    file << " format=\"ascii\"";
  }
  file << ">\n";
  MeshEntity* e;
  int order = m->getShape()->getOrder();
  static int vtkTypes[Mesh::TYPES][2] =
    /* order
       linear,quadratic
       V  V */
    {{ 1,-1}//vertex
    ,{ 3,21}//edge
    ,{ 5,22}//triangle
    ,{ 9,23}//quad
    ,{10,24}//tet
    ,{12,25}//hex
    ,{13,-1}//prism
    ,{14,-1}//pyramid
  };
  if (isWritingBinary)
  {
    unsigned int dataLen = 0;
    MeshIterator* elements = m->begin(cellDim);
    while ((e = m->iterate(elements)))
    {
      dataLen++;
    }
    m->end(elements);
    unsigned int dataLenBytes = dataLen*sizeof(uint8_t);
    uint8_t* dataToEncode = new uint8_t[dataLen]();
    elements = m->begin(cellDim);
    unsigned int dataIndex = 0;
    while ((e = m->iterate(elements)))
    {
      dataToEncode[dataIndex] = vtkTypes[m->getType(e)][order-1];
      dataIndex++;
    }
    m->end(elements);
    writeEncodedArray(file, dataLenBytes, (char*)dataToEncode);
    delete [] dataToEncode;
  }
  else
  {
    MeshIterator* elements = m->begin(cellDim);
    while ((e = m->iterate(elements)))
    {
      file << vtkTypes[m->getType(e)][order-1] << '\n';
    }
    m->end(elements);
  }
  file << "</DataArray>\n";
}

static void writeCells(std::ostream& file,
    Numbering* n,
    bool isWritingBinary,
    int cellDim)
{
  file << "<Cells>\n";
  writeConnectivity(file, n, isWritingBinary, cellDim);
  writeOffsets(file, n, isWritingBinary, cellDim);
  writeTypes(file, n->getMesh(), isWritingBinary, cellDim);
  file << "</Cells>\n";
}

static void writePointData(std::ostream& file,
    Mesh* m,
    DynamicArray<Node>& nodes,
    std::vector<std::string> writeFields,
    bool isWritingBinary = false)
{
  file << "<PointData>\n";
  for (int i=0; i < m->countFields(); ++i)
  {
    Field* f = m->getField(i);
    if (isNodal(f) && shouldPrint(f,writeFields))
    {
      writeNodalField<double>(file,f,nodes,isWritingBinary);
    }
  }
  for (int i=0; i < m->countNumberings(); ++i)
  {
    Numbering* n = m->getNumbering(i);
    if (isNodal(n) && shouldPrint(n,writeFields))
    {
      writeNodalField<int>(file,n,nodes,isWritingBinary);
    }
  }
  for (int i=0; i < m->countGlobalNumberings(); ++i)
  {
    GlobalNumbering* n = m->getGlobalNumbering(i);
    if (isNodal(n) && shouldPrint(n,writeFields))
    {
      writeNodalField<long>(file,n,nodes,isWritingBinary);
    }
  }
  file << "</PointData>\n";
}

template <class T>
class WriteIPField : public FieldOp
{
  public:
    int point;
    int components;
    NewArray<T> ipData;
    FieldDataOf<T>* data;
    MeshEntity* entity;
    std::ostream* fp;
    bool isWritingBinary;
    int cellDim;

    T* dataToEncode;
    int dataIndex;

    virtual bool inEntity(MeshEntity* e)
    {
      entity = e;
      return true;
    }
    virtual void atNode(int node)
    {
      if (node != point)
        return;
      data->getNodeComponents(entity,node,&(ipData[0]));
      for (int i=0; i < components; ++i)
      {
        if (isWritingBinary)
        {
          // if we are writing base64 then populate the array to be encoded
          dataToEncode[dataIndex] = ipData[i];
          dataIndex++;
        }
        else
        {
          //otherwise simply write to the file
          (*fp) << workaround(ipData[i]) << ' ';
        }
      }
      if (!isWritingBinary)
      {
        (*fp) << '\n'; //newline for each node when not writing base64
      }
    }
    void runOnce(FieldBase* f)
    {
      std::string s = getIPName(f,point);
      components = f->countComponents();
      writeDataHeader(*fp,
        s.c_str(),
        f->getScalarType(),
        f->countComponents(),
        isWritingBinary);
      ipData.allocate(components);
      data = static_cast<FieldDataOf<T>*>(f->getData());

      if (isWritingBinary) //extra steps if we are writing base64
      {
        //calculate size of array
        Mesh* m = f->getMesh();
        int arraySize = m->count(cellDim)*components;

        dataToEncode = new T[arraySize](); //allocate space for array
        apply(f); //populate array

        //encode and write to file
        int dataLenBytes = arraySize * sizeof(T);
        writeEncodedArray( (*fp), dataLenBytes, (char*)dataToEncode);

        //free array
        delete [] dataToEncode;
        dataIndex = 0;
      }
      else
      {
        apply(f); //same function call if writing ASCII
      }
      (*fp) << "</DataArray>\n";
    }
    void run(std::ostream& file,
      FieldBase* f,
      bool isWritingBinaryArg,
      int cellDimArg)
    {
      isWritingBinary = isWritingBinaryArg;
      cellDim = cellDimArg;
      fp = &file;
      dataIndex = 0;
      int n = countIPs(f, cellDim);
      for (point=0; point < n; ++point)
      {
        runOnce(f);
      }
    }
};

static void writeCellParts(std::ostream& file,
    Mesh* m,
    bool isWritingBinary,
    int cellDim)
{
  writeDataHeader(file, "apf_part", apf::Mesh::INT, 1, isWritingBinary);
  size_t n = m->count(cellDim);
  int id = m->getId();
  if (isWritingBinary)
  {
    unsigned int dataLen = n;
    unsigned int dataLenBytes = dataLen*sizeof(int);
    int* dataToEncode = new int[dataLen]();
    for (size_t i = 0; i < n; ++i )
    {
      dataToEncode[i] = id;
    }
    writeEncodedArray(file, dataLenBytes, (char*)dataToEncode);
    file << "</DataArray>\n";
    delete [] dataToEncode;
  }
  else
  {
    for (size_t i = 0; i < n; ++i)
    {
      file << id << '\n';
    }
    file << "</DataArray>\n";
  }
}

static void writeCellData(std::ostream& file,
    Mesh* m,
    std::vector<std::string> writeFields,
    bool isWritingBinary,
    int cellDim)
{
  file << "<CellData>\n";
  WriteIPField<double> wd;
  for (int i=0; i < m->countFields(); ++i)
  {
    Field* f = m->getField(i);
    if (isIP(f, cellDim) && shouldPrint(f,writeFields))
    {
      wd.run(file, f, isWritingBinary, cellDim);
    }
  }
  WriteIPField<int> wi;
  for (int i=0; i < m->countNumberings(); ++i)
  {
    Numbering* n = m->getNumbering(i);
    if (isIP(n, cellDim) && shouldPrint(n,writeFields))
    {
      wi.run(file, n, isWritingBinary, cellDim);
    }
  }
  WriteIPField<long> wl;
  for (int i=0; i < m->countGlobalNumberings(); ++i)
  {
    GlobalNumbering* n = m->getGlobalNumbering(i);
    if (isIP(n, cellDim) && shouldPrint(n,writeFields))
    {
      wl.run(file, n, isWritingBinary, cellDim);
    }
  }
  writeCellParts(file, m, isWritingBinary, cellDim);
  file << "</CellData>\n";
}

//function checks if the machine is a big or little endian machine
bool isBigEndian()
{
  union {
      unsigned int i;
      char c[4];
  } bint = {0x01020304};
  return bint.c[0] == 1;
}

static std::string getFileNameAndPathVtu(const char* prefix,
    std::string fileName,
    int id)
{
  int dirNum = id/1024;
  std::stringstream ss;
  ss << prefix << '/' << dirNum << '/' << fileName;
  return ss.str();
}

static void writeVtuFile(const char* prefix,
    Numbering* n,
    std::vector<std::string> writeFields,
    bool isWritingBinary,
    int cellDim)
{
  double t0 = PCU_Time();
  std::string fileName = getPieceFileName(PCU_Comm_Self());
  std::string fileNameAndPath = getFileNameAndPathVtu(prefix, fileName, PCU_Comm_Self());
  std::stringstream buf;
  Mesh* m = n->getMesh();
  DynamicArray<Node> nodes;
  getNodes(n,nodes);
  buf << "<VTKFile type=\"UnstructuredGrid\"";
  if (isWritingBinary)
  {
    buf << " byte_order=";
    if (isBigEndian())
    {
      buf << "\"BigEndian\"";
    }
    else
    {
      buf << "\"LittleEndian\"";
    }
    if (lion::can_compress )
    {
      //TODO determine what the header_type should be definitively
      buf << " header_type=\"UInt64\"";
      buf << " compressor=\"vtkZLibDataCompressor\"";
    }
    else
    {
      buf << " header_type=\"UInt32\"";
    }
  }
  buf<< ">\n";
  buf << "<UnstructuredGrid>\n";
  buf << "<Piece NumberOfPoints=\"" << nodes.getSize();
  buf << "\" NumberOfCells=\"" << m->count(cellDim);
  buf << "\">\n";
  writePoints(buf,m,nodes,isWritingBinary);
  writeCells(buf, n, isWritingBinary, cellDim);
  writePointData(buf,m,nodes,writeFields,isWritingBinary);
  writeCellData(buf, m, writeFields, isWritingBinary, cellDim);
  buf << "</Piece>\n";
  buf << "</UnstructuredGrid>\n";
  buf << "</VTKFile>\n";
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
  {
    printf("writeVtuFile into buffers: %f seconds\n", t1 - t0);
  }
  { //block forces std::ofstream destructor call
    std::ofstream file(fileNameAndPath.c_str());
    assert(file.is_open());
    file << buf.rdbuf();
  }
  double t2 = PCU_Time();
  if (!PCU_Comm_Self())
  {
    printf("writeVtuFile buffers to disk: %f seconds\n", t2 - t1);
  }
}

static void readVtuFile(const char* prefix, Mesh2* m)
{
    std::string fileName = getPieceFileName(PCU_Comm_Self());
    std::string fileNameAndPath = getFileNameAndPathVtu(prefix, fileName, PCU_Comm_Self());
    std::ifstream file(fileNameAndPath.c_str());
    assert(file.is_open());

    int nPoints, nCells;
    bool lendian = false, rbinary = false;
    int dim = 1;
    std::string line;

    /// Line with byte order (big endian by default)
    /// TODO header type (compression)
    std::getline(file, line);
    if (line.find("LittleEndian") != std::string::npos)
        lendian = true;
    if(line.find("compressor") != std::string::npos){
        std::cerr << "VTK read does not support compression yet, aborting..." << std::endl;
        abort();
    }
    /// Skip line
    std::getline(file, line);
    /// Line with NumberOfPoints and NumberOfCells
    std::getline(file, line);
    sscanf(line.c_str(),"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">", &nPoints, &nCells);
    /// Skip line
    std::getline(file, line);
    /// Line for DataArray with type, name, space dim and format
    /// TODO type (Float64 by default)
    std::getline(file, line);
    if (line.find("binary") != std::string::npos)
        rbinary = true;
    if (line.find("NumberOfComponents=\"2\"") != std::string::npos)
        dim = 2;
    else if (line.find("NumberOfComponents=\"3\"") != std::string::npos)
        dim = 3;

    if (dim < 3) {
        std::cerr << "VTK read supports 3D meshes for now, aborting..." << std::endl;
        abort();
    }

    if (!PCU_Comm_Self())
    {
        std::cout << "Reading VTK mesh. Format: ";
        if (rbinary) {
            std::cout << "binary. Endian: ";
            if (lendian)
                std::cout << "little";
            else
                std::cout << "big";
        }
        else
            std::cout << "ascii.";

        std::cout << ". Mesh dim: " << dim << ". Number of vertices: " << nPoints;
        std::cout << ". Number of regions: " << nCells << std::endl;
    }

    /// Read the vertex coordinates
    double* coords;
    int ndoubles = nPoints*3;
    coords = new double[ndoubles];
    if (rbinary) {
        std::getline(file, line);

        /// Decoding the line
        std::string dline = lion::base64Decode(line);

        char* cline = new char[8];

        /// Getting the header with the decoded buffer size
        for (int i = 0; i < 4; ++i)
            cline[i] = dline[i];
        int nbytes = *(int*)(cline);
        assert(nbytes == ndoubles*8);

        /// Shift to the data buffer
        dline = dline.substr(4);
        for (int i = 0; i < ndoubles; ++i) {
            for (int j = 0 ; j < 8; ++j)
                cline[j] = dline[i*8+j];
            coords[i] = *(double*)cline;
        }
    }
    else {
        for (int i = 0; i < nPoints; ++i) {
            std::getline(file, line);
            sscanf(line.c_str(), "%lf %lf %lf\n", coords+i*3, coords+i*3+1, coords+i*3+2);
        }
    }

    /// Skip 3 lines
    for(int i = 0; i < 3; ++i)
        std::getline(file, line);
    /// TODO difference betwee type Int32 and Int64
    std::getline(file, line);

    /// Read the connectivity
    int* conn;
    int nelem;
    if (dim == 2)
        nelem = 3;
    else
        nelem = 4;
    conn = new int[nCells*nelem];

    /// Suppot 3D meshes for the moment
    if (rbinary) {
        std::getline(file, line);

        /// Decoding the line
        std::string dline = lion::base64Decode(line);

        char* cline = new char[8];

        /// Getting the header with the decoded buffer size
        for (int i = 0; i < 4; ++i)
            cline[i] = dline[i];
        int nbytes = *(int*)(cline);
        assert(nbytes == nCells*nelem*4);

        /// Shift to the data buffer
        dline = dline.substr(4);
        for (int i = 0; i < nCells*nelem; ++i) {
            for (int j = 0 ; j < 4; ++j)
                cline[j] = dline[i*4+j];
            conn[i] = *(int*)cline;
        }
    }
    else {
        for (int i = 0; i < nCells; ++i) {
            std::getline(file, line);
            sscanf(line.c_str(), "%d %d %d %d\n", conn+i*nelem, conn+i*nelem+1, conn+i*nelem+2, conn+i*nelem+3);
        }
    }

    /// TODO read fields
    file.close();

    /// Constructing the mesh
    GlobalToVert outMap;
    construct(m, conn, nCells, Mesh::TET, outMap);
    delete [] conn;
    alignMdsRemotes(m);
    deriveMdsModel(m);
    setCoords(m, coords, nPoints, outMap);
    delete [] coords;
    outMap.clear();
}


static void safe_mkdir(const char* path)
{
  mode_t const mode = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;
  int err;
  errno = 0;
  err = mkdir(path, mode);
  if (err != 0 && errno != EEXIST)
  {
    reel_fail("MDS: could not create directory \"%s\"\n", path);
  }
}

static void makeVtuSubdirectories(const char* prefix, int numParts)
{
  std::stringstream ss1;
  ss1 << prefix;
  std::string prefixStr = ss1.str();
  int numDirectories = numParts/1024;
  if (numParts % 1024 != 0)
  {
    numDirectories++;
  }
  for (int i = 0; i < numDirectories; i++)
  {
    std::stringstream ss2;
    ss2 << prefix << '/' << i;
    safe_mkdir(ss2.str().c_str());
  }
}

void writeVtkFilesRunner(const char* prefix,
    Mesh* m,
    std::vector<std::string> writeFields,
    bool isWritingBinary,
    int cellDim)
{
  if (cellDim == -1) cellDim = m->getDimension();
  double t0 = PCU_Time();
  if (!PCU_Comm_Self())
  {
    safe_mkdir(prefix);
    makeVtuSubdirectories(prefix, PCU_Comm_Peers());
    writePvtuFile(prefix, m, writeFields, isWritingBinary, cellDim);
  }
  PCU_Barrier();
  Numbering* n = numberOverlapNodes(m,"apf_vtk_number");
  m->removeNumbering(n);
  writeVtuFile(prefix, n, writeFields, isWritingBinary, cellDim);
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
  {
    printf("vtk files %s written in %f seconds\n", prefix, t1 - t0);
  }
  delete n;
}

std::vector<std::string> populateWriteFields(Mesh* m)
{
  std::vector<std::string> writeFields;
  for (int i=0; i < m->countFields(); ++i)
  {
    Field* f = m->getField(i);
    if (isPrintable(f))
      writeFields.push_back(f->getName());
  }
  for (int i=0; i < m->countNumberings(); ++i)
  {
    Numbering* n = m->getNumbering(i);
    if (isPrintable(n))
      writeFields.push_back(n->getName());
  }
  for (int i=0; i < m->countGlobalNumberings(); ++i)
  {
    GlobalNumbering* n = m->getGlobalNumbering(i);
    if (isPrintable(n))
      writeFields.push_back(n->getName());
  }
  return writeFields;
}

void writeOneVtkFile(const char* prefix, Mesh* m)
{
  /* creating a non-collective numbering is
     a tad bit risky, but we should be fine
     given the current state of the code */

  // bool isWritingBinary = true;
  Numbering* n = numberOverlapNodes(m,"apf_vtk_number");
  m->removeNumbering(n);
  std::vector<std::string> writeFields = populateWriteFields(m);
  writeVtuFile(prefix, n, writeFields, false, m->getDimension());
  delete n;
}

void writeVtkFiles(
    const char* prefix,
    Mesh* m,
    std::vector<std::string> writeFields,
    int cellDim)
{
  writeVtkFilesRunner(prefix, m, writeFields, false, cellDim);
}

void writeVtkFiles(const char* prefix, Mesh* m, int cellDim)
{
  std::vector<std::string> writeFields = populateWriteFields(m);
  writeVtkFiles(prefix, m, writeFields, cellDim);
}

void readVtkFiles(const char* prefix, Mesh2* m)
{
    bool isWritingBinary = false;
    double t0 = PCU_Time();
    if (!PCU_Comm_Self())
    {
      readPvtuFile(prefix, isWritingBinary);
    }
    PCU_Barrier();

    readVtuFile(prefix, m);
    double t1 = PCU_Time();
    if (!PCU_Comm_Self())
    {
      printf("vtk files %s read in %f seconds\n", prefix, t1 - t0);
    }
}

void writeASCIIVtkFiles(
    const char* prefix,
    Mesh* m,
    std::vector<std::string> writeFields)
{
  writeVtkFilesRunner(prefix, m, writeFields, false, -1);
}

void writeASCIIVtkFiles(const char* prefix, Mesh* m)
{
  std::vector<std::string> writeFields = populateWriteFields(m);
  writeASCIIVtkFiles(prefix, m, writeFields);
}

}
