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
#include "gdal_utils.h"
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

	std::string look_for_var(cGeophysicsNcFile& N,  const std::vector<std::string>& candidates)
	{
		std::vector<NcVar> v = N.getAllVars();		
		for(size_t ci=0; ci<candidates.size(); ci++){
			for (size_t vi = 0; vi < v.size(); vi++) {
				if(strcasecmp(v[vi].getName(), candidates[ci]) == 0){
					return v[vi].getName();
				}
			}
		}
		return std::string();
	}

	std::vector<double> get_linetype(cGeophysicsNcFile& N)
	{		
		std::vector<std::string> candidates;
		candidates.push_back("linetype");
		candidates.push_back("line_type");
		candidates.push_back("ltype");
		std::string vname = look_for_var(N, candidates);
		std::vector<double> ltype;

		if (vname.size() > 0) {
			NcVar var = N.getVar(vname);
			if (N.isLineVar(var)) {
				cLineVar v = N.getLineVar(vname);
				v.getAll(ltype);
			}
			else if (N.isSampleVar(var)) {
				ltype.resize(N.nlines());
				cSampleVar v = N.getSampleVar(vname);
				for (size_t li = 0; li < N.nlines(); li++) {
					v.getSample(li, 0, 0, ltype[li]);
				}
			}
		}
		return ltype;
	}

	bool process() {		
		cGeoDataset D = cGeoDataset::create_shapefile(ShapePath);
		cLayer L = D.create_layer("flight_lines", OGRwkbGeometryType::wkbLineString);
		std::vector<cAttribute> atts;
		atts.push_back(cAttribute("linenumber", (int)0));
		atts.push_back(cAttribute("linetype", (int)0));
		L.add_fields(atts);

		cGeophysicsNcFile N(NCPath);
		std::vector<unsigned int> ln;
		N.getLineNumbers(ln);
		const size_t nl = N.nlines();
				
		std::vector<double> ltype = get_linetype(N);				
		std::vector<double> lkm0(nl);
		std::vector<double> lkm2(nl);
		std::vector<double> lkm4(nl);
		for (size_t li = 0; li < nl; li++) {	
			lkm0[li] = 0.0; 
			lkm2[li] = 0.0;
			lkm4[li] = 0.0;

			std::vector<double> x;
			std::vector<double> y;			
			N.getDataByLineIndex("longitude", li, x);
			N.getDataByLineIndex("latitude", li, y);
			std::vector<double> xout;
			std::vector<double> yout;

			double null = defaultmissingvalue(ncDouble);

			int ns = (int) x.size();
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
				int minpoints = 20;
				int ss = (kend - kstart) / (minpoints - 2);
				if (ss < 1) ss = 1;
				for (k = kstart; k < kend; k += ss) {
					if(x[k] != null && y[k] != null) {
						xout.push_back(x[k]);
						yout.push_back(y[k]);						
					}					
				}
				xout.push_back(x[kend]);
				yout.push_back(y[kend]);			

				atts[0].value = (int)ln[li];
				atts[1].value = (int)0;
				if (ltype.size() == nl) {
					atts[1].value = (int)ltype[li];
				}
				L.add_linestring_feature(atts, xout, yout);

				continue;
				
				double meanlong = mean(xout);
				int zone = std::ceil((meanlong + 180.0) / 6.0);
				
				//std::string zonestr = strprint("MGA%02d", zone);
				//int epsgcode_geodetic = getepsgcode("GDA94", "GEODETIC");
				//int epsgcode_utm      = getepsgcode("GDA94", zonestr);

				std::string zonestr   = strprint("UTM%02dS", zone);
				int epsgcode_geodetic = getepsgcode("WGS84", "GEODETIC");
				//int epsgcode_utm      = getepsgcode("WGS84", zonestr);
				int epsgcode_utm = 32752;

			    printf("zone=%d\n",zone);
			    printf("zonestr=%s\n",zonestr.c_str());
				printf("epsgcode_geodetic=%d\n",epsgcode_geodetic);
			    printf("epsgcode_utm     =%d\n",epsgcode_utm);
			    
				for (size_t j = 1; j < xout.size(); j++) {
				    printf("%lf,%lf\n",xout[j],yout[j]);
				}

				std::vector<double> e;
				std::vector<double> n;
				transform(epsgcode_geodetic, xout, yout, epsgcode_utm, e, n);
				
				for (size_t j = 0; j < xout.size(); j++) {
				    printf("%lf,%lf\n",xout[j],yout[j]);
				}

				double d = 0.0;
				for (size_t j = 1; j < e.size(); j++) {
					d += std::hypot(e[j] - e[j - 1], n[j] - n[j - 1]);
				}

				lkm0[li] = d / 1000.0;				
				if (ltype.size() == nl) {
					if (ltype[li] == 2) lkm2[li] = lkm0[li];
					if (ltype[li] == 4) lkm4[li] = lkm0[li];
				}		
				
			}
		}		
		
		//std::string a = strprint("%s,%.1lf,%.1lf,%.1lf\n",NCPath.c_str(),sum(lkm0), sum(lkm2), sum(lkm4));
		//glog.logmsg(a);
		//std::printf(a.c_str());
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
					//if (exists(ShapePath) == false) {
						std::cout << "[" << mpirank << "] " << NCPath << " " << ShapePath << std::endl << std::flush;
						cNcToShapefileConverter C(NCPath, ShapePath, mpisize, mpirank);
					//}
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



