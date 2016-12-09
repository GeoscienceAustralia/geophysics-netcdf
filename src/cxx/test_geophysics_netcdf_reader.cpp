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

FILE* global_log_file = NULL;

bool test_magnetics(){

	bool status;		
	std::string indir,ncpath;
	//indir = "z:\\projects\\geophysics_netcdf\\conversion_scripts\\ncfiles\\";
	indir = "http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/mag_database_reformat_2016_adjusted/netcdf/";
	ncpath = indir + "GSSA_P1255MAG_Marree.nc";

	//Open the file and initialise the indexes
	cGeophysicsNcFile ncfile(ncpath, NcFile::read);

	//Get the line numbers using the standard_name attributes
	std::vector<int> linenumber = ncfile.getLineNumbers();
	std::vector<int> flightnumber = ncfile.getFlightNumbers();

	std::string stdname_x         = "longitude";
	std::string stdname_y         = "latitude";
	std::string stdname_tielev    = "total_magnetic_intensity_anomaly_tie_line_levelled";
	std::string stdname_mlev      = "total_magnetic_intensity_anomaly_micro_levelled";
	std::string stdname_awagslev  = "total_magnetic_intensity_anomaly_datum_levelled";
	
	std::string xvarname = ncfile.getVarNameByStandardName(stdname_x);
	std::string yvarname = ncfile.getVarNameByStandardName(stdname_y);
	std::string mvarname = ncfile.getVarNameByStandardName(stdname_awagslev);

	//Get 21th line number and index
	size_t lnum = linenumber[20];
	size_t lind = ncfile.getLineIndex(lnum);
	
	//Get mag data for 21th line by number and index
	std::vector<double> v1darray;
	status = ncfile.getDataByLineNumber(mvarname, lnum, v1darray);
	status = ncfile.getDataByLineIndex(mvarname, lind, v1darray);

	//Loop over all lines getting the x, y, and mag values
	std::vector<double> x, y;
	std::vector<float> mag;
	int date;
	cSampleVar xvar = ncfile.getSampleVar(xvarname);
	cSampleVar yvar = ncfile.getSampleVar(yvarname);
	cSampleVar mvar = ncfile.getSampleVar(mvarname);
	
	double t, ta, tb;
	ta = gettime(); status = xvar.getAll(x); tb = gettime();
	logmsg("Get all of double x nsample=%lu time=%lf\n", x.size(), tb - ta);
	ta = gettime(); status = yvar.getAll(x); tb = gettime();
	logmsg("Get all of double y nsample=%lu time=%lf\n", x.size(), tb - ta);
	ta = gettime(); status = mvar.getAll(mag); tb = gettime();
	logmsg("Get all of float  mag nsample=%lu time=%lf\n", x.size(), tb - ta);

	t = 0;
	//for (size_t li = 0; li < ncfile.nlines(); li++){
	for (size_t li = 0; li < 100; li++){		
		ta = gettime();	
		status = xvar.getLine(li, x);
		status = yvar.getLine(li, y);
		status = mvar.getLine(li, mag);		
		tb = gettime();
		logmsg("%lu nsamples=%lu time=%lf\n", li, x.size(), tb - ta);
		t += (tb - ta);
	}		
	logmsg("Total time=%lf\n", t);
	return true;
}

bool test_aem_conductivity(){

	bool status;
	std::string indir  = "Z:\\projects\\geophysics_netcdf\\aem\\temp\\";			
	std::string ncpath = indir + "AUS_10009_Ord_Bonaparte_LCI.nc";

	//Open the file
	cGeophysicsNcFile ncfile(ncpath, NcFile::write);			
		
	//Get the line numbers
	std::vector<int> linenumber = ncfile.getLineNumbers();

	//Get conductivity variable name by its standard name attribute
	std::string stdname_conductivity = "layer_conductivity";
	std::string varname = ncfile.getVarNameByStandardName(stdname_conductivity);
	
	//Get conductivity all at once in 1d array
	std::vector<float> c1darray;
	cSampleVar vc = ncfile.getSampleVar(varname);
	status = vc.getAll(c1darray);

	//Get conductivity line by line and band by band in 1d array
	for (size_t li = 0; li < ncfile.nlines(); li++){
		for (size_t bi = 0; bi < vc.nbands(); bi++){			
			status = vc.getLine(li, bi, c1darray);
		}
	}	

	//Get conductivity line by line in 2d array
	std::vector<std::vector<float>> c2darray;
	for (size_t li = 0; li < ncfile.nlines(); li++){
		status = ncfile.getDataByLineIndex(varname, li, c2darray);
	}

	return true;	
}

bool test_create(){
	std::string ncpath = "test.nc";
	deletefile(ncpath);
	std::vector<size_t> linenumbers  = { 100, 200, 300, 400 };
	std::vector<size_t> linensamples = { 10,   20,  30,  40 };
	std::vector<size_t> flightnumbers = {11, 22, 33, 44 };

	cGeophysicsNcFile   nc(ncpath,NcFile::FileMode::replace);
	nc.InitialiseNew(linenumbers, linensamples);
	size_t nwindows = 45;
	size_t nlayers  = 30;
	size_t nrxcomponents = 3;

	size_t ntotalsamples = nc.ntotalsamples();
	std::vector<int> fid = increment((int)ntotalsamples);
	std::vector<int> layers = increment((int)nlayers);
	std::vector<int> windows = increment((int)nwindows);
	std::vector<int> rxcomponents = increment((int)nrxcomponents);
	
	NcDim dim_rxcomponent = nc.addDimVar("rxcomponents", rxcomponents);
	NcDim dim_window = nc.addDimVar("windows", windows);	
	NcDim dim_layer  = nc.addDimVar("layers", layers);
	
	cSampleVar vfid = nc.addSampleVar("fiducial", ncInt);
	vfid.add_standard_name("fiducial");
	vfid.add_units("1");
	vfid.add_missing_value(34);


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
	std::vector<double> em = increment((double)n,0.0,1.0);		
	vfid.putAll(fid);
	vem.putAll(em);

	for (size_t li = 0; li < nc.nlines(); li++){
		size_t nls = nc.nlinesamples(li);
		std::vector<double> x = increment((double)nls,500000.0,10.0);
		std::vector<double> y = increment((double)nls,6500000.0,10.0);
		
		vx.putLine(li, x);
		vy.putLine(li, y);
		
		for (size_t bi = 0; bi < nlayers; bi++){
			std::vector<double> c(nls, li*10+bi);
			vconductivity.putLine(li, bi, c);
		}

		for (size_t bi = 0; bi < nlayers; bi++){
			std::vector<double> t(nls, li * 100 + bi);
			vthickness.putLine(li, bi, t);
		}
		
	}		
	return true;
};

int main(int argc, char** argv)
{
	_GSTITEM_	
	try
	{	
		logmsg("Opening log file\n");
		global_log_file = fopen("test.log", "w");		
		logmsg("Log file opened\n");

		test_magnetics();
		//test_aem_conductivity();
		
		logmsg("Closing log file\n");
		fclose(global_log_file);
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


