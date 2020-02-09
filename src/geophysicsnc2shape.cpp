/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cstdio>
#include <netcdf>
#include <vector>
#include <limits>

using namespace netCDF;
using namespace netCDF::exceptions;

#ifdef USEGLOBALSTACKTRACE
#include "stacktrace.h"
class cStackTrace gtrace;
#endif

#define _PROGRAM_ "geophysicsnc2shape"
#define _VERSION_ "1.0"

#include <mpi.h>

#include "general_utils.h"
#include "file_utils.h"
#include "logger.h"
#include "ogr_utils.h"
#include "geophysics_netcdf.h"

#ifdef HAVE_GDAL
	#include "crs.h"
#endif

class cLogger glog; //The instance of the global log file manager

class cNcToShapefileConverter {
	int MPISize;
	int MPIRank;	
	std::string NCPath;
	std::string ShapePath;	
	std::string LogFile;	

public:

	cNcToShapefileConverter(const std::string& ncfilepath, const std::string& shapefilepath, 		const int mpisize, const int mpirank) {
		_GSTITEM_ 
		MPISize = mpisize;
		MPIRank = mpirank;
				
		NCPath    = fixseparator(ncfilepath);
		ShapePath = fixseparator(shapefilepath);
		LogFile   = ShapePath + ".log";
		
		glog.open(LogFile);
		glog.log("Program %s \n", _PROGRAM_);
		glog.log("Version %s Compiled at %s on %s\n", _VERSION_, __TIME__, __DATE__);
		glog.log("Working directory %s\n", getcurrentdirectory().c_str());		
		bool status = process();	
		if (status == false) {
			glog.logmsg(MPIRank,"Error 0: creating shapefile %s from %s\n",ShapePath.c_str(),NCPath.c_str());
		}
		glog.close();		
	};

	~cNcToShapefileConverter() {};

	bool process() {		
		cGeoDataset D = cGeoDataset::create_shapefile(ShapePath);
		cLayer L = D.create_layer("flight_lines", OGRwkbGeometryType::wkbLineString);
		std::vector<cAttribute> atts;
		atts.push_back(cAttribute("linenumber", (int)0));
		L.add_fields(atts);

		cGeophysicsNcFile N(NCPath);
		std::vector<unsigned int> ln;
		bool status = N.getLineNumbers(ln);
		//std::vector<unsigned int> ln = N.getLineNumbers();
		for (size_t i = 0; i < ln.size(); i++) {
			
			//std::cout << i << std::endl;
			//if (i == 267){
			//	int dummy = 0;
			//}

			std::vector<double> x;
			std::vector<double> y;
			bool status = false;
			status = N.getDataByLineIndex("longitude", i, x);
			status = N.getDataByLineIndex("latitude", i, y);
			std::vector<double> xout;
			std::vector<double> yout;

			double null = defaultmissingvalue(ncDouble);
			size_t ns = x.size();
			int k = 0;
			int kstart, kend;
			while (x[k] == null || y[k] == null) {				
				k++;
				if (k == ns)break;
			}
			kstart = k;
			
			k = ns - 1;
			while (x[k] == null || y[k] == null) {				
				k--;
				if (k == -1)break;
			}
			kend = k;
			
			if(kend >= kstart){
				int minpoints = 5;
				int ss = (kend - kstart) / (minpoints - 2);
				if (ss < 1) ss = 1;
				for (k = kstart; k < kend; k += ss) {
					xout.push_back(x[k]);
					yout.push_back(y[k]);
				}
				xout.push_back(x[kend]);
				yout.push_back(y[kend]);			
				atts[0].value = (int)ln[i];
				L.add_linestring_feature(atts, xout, yout);
			}
		}		
		return true;
	}
};

int main(int argc, char** argv)
{
	_GSTITEM_

	GDALAllRegister();

	int mpisize;
	int mpirank;	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

	glog.logmsg(0,"Program %s \n", _PROGRAM_);
	glog.logmsg(0,"Version %s Compiled at %s on %s\n", _VERSION_, __TIME__, __DATE__);
	glog.logmsg(0,"Working directory %s\n", getcurrentdirectory().c_str());

	try
	{
		std::string ncdir     = argv[1];
		std::string shapedir  = argv[2];
		std::string listfile  = argv[3];
		std::ifstream file(listfile);
		addtrailingseparator(ncdir);
		addtrailingseparator(shapedir);
		int k = 0;
		while (file.eof() == false) {
			std::string s;	
			file >> s;
			s = trim(s);
			if (s.size() > 0 && s[0] != '#') {
				sFilePathParts fpp = getfilepathparts(s);
				std::string NCPath    = ncdir    + fpp.directory + fpp.prefix + ".nc";
				std::string ShapePath = shapedir + fpp.directory + fpp.prefix + ".shp";
				
				if (k % mpisize == mpirank) {
					if (exists(ShapePath) == false) {
						std::cout << "[" << mpirank << "] " << NCPath << " " << ShapePath << std::endl << std::flush;
						cNcToShapefileConverter C(NCPath, ShapePath, mpisize, mpirank);
					}
				}
				k++;
			}
		}		
	}
	catch (NcException& e)
	{
		_GSTPRINT_
		std::cout << e.what() << std::endl;
		MPI_Finalize();
		return 1;
	}
	catch (std::exception& e)
	{
		_GSTPRINT_
		std::cout << e.what() << std::endl;
		MPI_Finalize();
		return 1;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	glog.logmsg(0,"Finished\n");
	MPI_Finalize();	
	return 0;
}



