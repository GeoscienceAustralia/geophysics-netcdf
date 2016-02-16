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
using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

#define USEGLOBALSTACKTRACE
#ifdef USEGLOBALSTACKTRACE
#include "stacktrace.h"
cStackTrace globalstacktrace;
#endif

#include "general_utils.h"
#include "file_utils.h"
#include "blocklanguage.h"

class cGeophysicsNcFile : public NcFile{
	
private:

	std::vector<int> line_index_start;
	std::vector<int> line_index_count;
	
	bool setLineIndex(){

		bool status;
		size_t ns = getDim("sample").getSize();
		status = getVar("index_lines", line_index_start);
		line_index_start.pop_back();

		size_t nl = nlines();
		line_index_count.resize(nl);
		for (size_t i = 0; i < nl-1; i++){
			line_index_count[i] = line_index_start[i + 1] - line_index_start[i];
		}	
		line_index_count[nl-1] = ns - line_index_start[nl-1];
		return true;

	}

public:

	cGeophysicsNcFile(const std::string& ncpath, const FileMode& filemode)
		: NcFile(ncpath, filemode)
	{
		setLineIndex();
	};

	size_t nlines(){ return line_index_start.size(); }

	template<typename T>
	bool getLineNumbers(std::vector<T>& vals){
		return getVar("LINE", vals);
	}

	template<typename T>
	bool getFlightNumbers(std::vector<T>& vals){
		return getVar("flight", vals);
	}

	template<typename T>
	bool getVar(const std::string& varname, std::vector<T>& vals){
		NcVar var = NcFile::getVar(varname);
		std::vector<NcDim> dims = var.getDims();
		size_t ne = 1;
		for (size_t i = 0; i < dims.size(); i++) ne *= dims[i].getSize();		
		vals.resize(ne);
		var.getVar(vals.data());
		return true;
	}

	template<typename T>
	bool getVarByLineIndex(const std::string& varname, std::vector<T>& vals, const size_t& lineindex){
		NcVar var = NcFile::getVar(varname);
		std::vector<NcDim> dims = var.getDims();

		std::vector<size_t> start(1);
		std::vector<size_t> count(1);
		start[0] = line_index_start[lineindex];
		count[0] = line_index_count[lineindex];

		size_t ne = count[0];
		for (size_t i = 1; i < dims.size(); i++){
			ne *= dims[i].getSize();
		}
		vals.resize(ne);				
		var.getVar(start, count, vals.data());
		return true;
	}

	

	
	
};

int main(int argc, char** argv)
{
	_GSTITEM_
	//if (argc != 2){
	//	printf("Usage: %s nc_file_path\n", argv[0]);
	//	return 1;
	//}

	try
	{
		//std::string ncpath = argv[1];
		//std::string ncpath = "http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/magrad_tests_indexed_v2/GSQP1029MAG.nc";
		std::string ncpath   = "http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/magrad_tests_indexed_v2/P583MAG.nc";
		//std::string ncpath = "Z:\\projects\\intrepid2netcdf\\ncfiles\\P583MAG.nc";

		cGeophysicsNcFile ncfile(ncpath, NcFile::read);
		
		//int ndims = ncfile.getDimCount();
		//std::multimap<string, NcDim> dims = ncfile.getDims();
		//NcDim dim_line = dims.find("line")->second;
		//NcDim dim_line_index = dims.find("line_index")->second;		
		//auto it = dims.find("sample");
		//if (it != dims.end()){
		//	NcDim d = dims.find("sample")->second;
		//	size_t dimname = d.getSize();
		//}

		std::vector<int> linenumber;
		ncfile.getLineNumbers(linenumber);

		std::vector<int> flightnumber;
		ncfile.getFlightNumbers(flightnumber);

		std::vector<double> datavector;
		ncfile.getVar("mag_microLevelled",datavector);
		
		std::vector<double> v1,v2,v3;
		for (size_t i = 0; i < ncfile.nlines(); i=i+10){			
			double t1 = gettime();			
			ncfile.getVarByLineIndex("mag_microLevelled", v1, i);
			ncfile.getVarByLineIndex("latitude_GDA94", v2, i);
			ncfile.getVarByLineIndex("longitude_GDA94", v3, i);
			double t2 = gettime();
			printf("%lu nsamples=%lu time=%lf\n", i, v1.size(), t2 - t1);
		}
		prompttocontinue();
		
	}

	catch (NcException& e)
	{
		_GSTPRINT_
		std::cout << e.what() << endl;
		return 1;
	}
	catch (std::exception& e)
	{
		_GSTPRINT_
		std::cout << e.what() << endl;
		return 1;
	}
	std::printf("Success\n");
	return 0;  // successfully terminated
}


