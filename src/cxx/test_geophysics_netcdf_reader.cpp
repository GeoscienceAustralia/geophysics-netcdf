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

#define USEGLOBALSTACKTRACE
#ifdef USEGLOBALSTACKTRACE
#include "stacktrace.h"
cStackTrace globalstacktrace;
#endif

using namespace netCDF;
using namespace netCDF::exceptions;

#include "general_utils.h"
#include "file_utils.h"
#include "vector_utils.h"
#include "geophysics_netcdf.h"


bool test_read1(){

	//std::string ncpath = argv[1];
	//std::string ncpath = "http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/magrad_tests_indexed_v2/GSQP1029MAG.nc";
	std::string ncpath   = "http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/magrad_tests_indexed_v2/P583MAG.nc";

	cGeophysicsNcFile ncfile(ncpath, NcFile::read);

	std::vector<int> linenumber;
	ncfile.getLineNumbers(linenumber);
	
	size_t lnum = linenumber[20];
	size_t index = ncfile.getLineIndex(lnum);
	std::vector<double> v0;
	bool status = ncfile.getVarByLineNumber("mag_microLevelled", 2440, v0);

	std::vector<std::vector<double>> v00;
	status = ncfile.getVarByLineNumber("mag_microLevelled", 2440, v00);

	std::vector<size_t> flightnumber;
	ncfile.getFlightNumbers(flightnumber);

	std::vector<double> v1, v2, v3;
	for (size_t i = 0; i < ncfile.nlines(); i++){
		double t1 = gettime();
		ncfile.getVarByLineIndex("mag_microLevelled", i, v1);
		ncfile.getVarByLineIndex("latitude_GDA94", i, v2);
		ncfile.getVarByLineIndex("longitude_GDA94", i, v3);
		double t2 = gettime();
		printf("%lu nsamples=%lu time=%lf\n", i, v1.size(), t2 - t1);
	}	
	return true;
}

bool test_read2(){

	std::string ncpath = "test.nc";

	cGeophysicsNcFile ncfile(ncpath, NcFile::read);

	std::vector<int> linenumber;
	ncfile.getLineNumbers(linenumber);

	size_t lnum = linenumber[3];
	size_t index = ncfile.getLineIndex(lnum);
	std::vector<double> v0;
	bool status;

	status = ncfile.getVarByLineNumber("easting", 300, v0);
	status = ncfile.getVarByLineNumber("northing", 300, v0);

	std::vector<size_t> flightnumber;
	ncfile.getFlightNumbers(flightnumber);

	std::vector<double> v1, v2, v3;
	for (size_t i = 0; i < ncfile.nlines(); i++){		
		ncfile.getVarByLineIndex("fiducial", i, v1);						
	}
	return true;
}

bool test_create(){
	std::string ncpath = "test.nc";
	deletefile(ncpath);
	std::vector<size_t> linenumbers  = { 100, 200, 300, 400 };
	std::vector<size_t> linensamples = { 10,   20,  30,  40 };
	std::vector<size_t> flightnumbers = {11, 22, 33, 44 };

	cGeophysicsNcFile   nc(ncpath, linenumbers, linensamples);
	size_t nwindows = 45;
	size_t nlayers  = 30;
	size_t nrxcomponents = 3;

	size_t ntotalsamples = nc.ntotalsamples();
	std::vector<int> fid     = increment(ntotalsamples,0,1);
	std::vector<int> layers  = increment(nlayers, 0, 1);
	std::vector<int> windows = increment(nwindows, 0, 1);
	std::vector<int> rxcomponents = increment(nrxcomponents, 0, 1);	
	
	NcDim dim_rxcomponent = nc.addDimVar("rxcomponents", rxcomponents);
	NcDim dim_window = nc.addDimVar("windows", windows);	
	NcDim dim_layer  = nc.addDimVar("layers", layers);
	
	cSampleVar vfid = nc.addSampleVar("fiducial", ncInt);
	vfid.add_standard_name("fiducial");
	vfid.add_units("1");
	vfid.add_fillvalue(34);


	cSampleVar vx   = nc.addSampleVar("easting",  ncDouble);
	vx.add_standard_name("X");
	vx.add_units("m");
	

	cSampleVar vy   = nc.addSampleVar("northing", ncDouble);	
	cSampleVar vconductivity = nc.addSampleVar("conductivity", ncDouble, dim_layer);
	cSampleVar vthickness    = nc.addSampleVar("thickness", ncDouble, dim_layer);
		
	cLineVar vflight = nc.addLineVar("flight", ncInt);
	vflight.putAll(flightnumbers);

	std::vector<NcDim> emdims = { dim_rxcomponent, dim_window };
	cSampleVar vem = nc.addSampleVar("em", ncDouble, emdims);
		
	const size_t n = ntotalsamples*nrxcomponents*nwindows;
	std::vector<double> em = increment(n,0.0,1.0);		
	vfid.putAll(fid);
	vem.putAll(em);

	for (size_t li = 0; li < nc.nlines(); li++){
		size_t nls = nc.nlinesamples(li);
		std::vector<double> x = increment(nls,500000.0,10.0);
		std::vector<double> y = increment(nls,6500000.0,10.0);
		
		vx.putLine(x, li);
		vy.putLine(y, li);
		
		for (size_t bi = 0; bi < nlayers; bi++){
			std::vector<double> c(nls, li*10+bi);
			vconductivity.putLine(c, li, bi);
		}

		for (size_t bi = 0; bi < nlayers; bi++){
			std::vector<double> t(nls, li * 100 + bi);
			vthickness.putLine(t, li, bi);
		}
		
	}		
	return true;
};

int main(int argc, char** argv)
{
	_GSTITEM_
	
	try
	{		
		test_create();
		test_read2();
	}

	catch (NcException& e)
	{
		_GSTPRINT_
		std::cout << e.what() << std::endl;
		return 1;
	}
	catch (std::exception& e)
	{
		_GSTPRINT_
		std::cout << e.what() << std::endl;
		return 1;
	}
	std::printf("Success\n");
	return 0;  // successfully terminated
}


